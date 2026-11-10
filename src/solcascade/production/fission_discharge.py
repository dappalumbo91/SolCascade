"""Name the heavy-metal discharge as a share of the ranked U-235 density.

The notes figure stays 45 GWd per tonne. The ranked density stays the filed
ledger. The page enrichment, the fuel-year companion, and the built-plant
efficiency stay outside this score. The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay_dress import _base_pieces, _dressings
from solcascade.production.fission_decay_seed import NamedGap, _five_pieces, _nearest, _size
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.fission_process import (
    ENRICHED_PERCENT,
    FABRICATION_CAP_PERCENT,
    build_process,
)
from solcascade.production.fission_safety import BUILT_EFFICIENCY
from solcascade.production.refine import project_root

_ANCHOR = "e^-3-e^-6"
_STEPS = (
    "(e^-3-e^-6)*(1+alpha*phi^-4)",
    "((e^-3-e^-6)*(1+alpha*phi^-4))*(1-alpha*pi^-5)",
    "(((e^-3-e^-6)*(1+alpha*phi^-4))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
)
_COUNT = 116
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4


@dataclass(frozen=True)
class DischargeStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class DischargeReport:
    anchor: NamedGap
    steps: tuple[DischargeStep, ...]
    target: float
    discharge_gwd: float
    enriched_percent: float
    fabrication_cap_percent: float
    companion_ratio: float
    built_efficiency: float
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_discharge.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _gap(derived: float, witness: float) -> float:
    if witness == 0.0:
        raise RuntimeError("the witness is zero")
    return abs(derived - witness) / abs(witness)


def _second_gap(rows, target: float, best: float) -> float:
    gaps = sorted(_gap(value, target) for _name, value in rows)
    if len(gaps) < 2 or abs(gaps[0] - best) > 1.0e-12:
        raise RuntimeError("the nearest dressing is missing from its own list")
    return gaps[1]


def _walk(target: float, anchor: NamedGap, alpha: float, helium: float, base) -> tuple[DischargeStep, ...]:
    name = anchor.name
    value = anchor.value
    steps: list[DischargeStep] = []
    for step_name in _STEPS:
        rows = _dressings(name, value, alpha, base)
        hit = _nearest(rows, target)
        if hit.name != step_name:
            raise RuntimeError("the discharge dressing changed")
        if hit.relative_gap >= _gap(value, target):
            raise RuntimeError("the discharge dressing did not tighten")
        size = _size(rows)
        if size.expressions != _COUNT or size.inside != _COUNT:
            raise RuntimeError("the discharge dressing count changed")
        if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
            raise RuntimeError("the discharge spacing is not the base-piece spacing")
        inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
        if inside_helium > 1:
            raise RuntimeError("several discharge dressings meet the helium standard")
        steps.append(
            DischargeStep(
                dressing=hit,
                inside_review=sum(_gap(number, target) <= _REVIEW_GAP for _row, number in rows),
                inside_tenth=sum(_gap(number, target) <= _TENTH for _row, number in rows),
                inside_helium=inside_helium,
                second_gap=_second_gap(rows, target, hit.relative_gap),
                expressions=size.expressions,
                median_spacing=size.median_spacing,
            )
        )
        name = hit.name
        value = hit.value
    closing = steps[-1]
    if closing.inside_helium != 1 or closing.dressing.relative_gap > helium:
        raise RuntimeError("the closing discharge figure misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError("a second discharge dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier discharge figure already meets helium")
    return tuple(steps)


def _leave(target: float, alpha: float, helium: float, five, base) -> str:
    anchor = _nearest(five, target)
    name = anchor.name
    value = anchor.value
    for _step in range(3):
        rows = _dressings(name, value, alpha, base)
        size = _size(rows)
        if size.expressions != _COUNT or _num(size.median_spacing) != _SPACING:
            return "count"
        hit = _nearest(rows, target)
        inside = sum(_gap(number, target) <= helium for _row, number in rows)
        if inside > 1:
            return "several"
        if hit.relative_gap >= _gap(value, target):
            return "loose"
        if inside == 1 and hit.relative_gap <= helium and _second_gap(rows, target, hit.relative_gap) > helium:
            return "unique"
        name = hit.name
        value = hit.value
    return "open"


def build_discharge() -> DischargeReport:
    process = build_process()
    module, digest = load_engine()
    if digest[:16] != process.pin_prefix:
        raise RuntimeError("the discharge pin does not match the plant pin")
    if process.discharge_gwd != 45.0:
        raise RuntimeError("the discharge is not the 45 GWd material balance")
    if FABRICATION_CAP_PERCENT != 5.0:
        raise RuntimeError("the fabrication cap left 5 percent")
    if process.discharge_seed.relative_gap <= _REVIEW_GAP:
        raise RuntimeError("the process discharge seed entered the review bar")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    five = _five_pieces(module)
    base = _base_pieces(module)
    target = process.discharge_over_ranked
    anchor = _nearest(five, target)
    if anchor.name != _ANCHOR:
        raise RuntimeError("the discharge anchor changed")
    steps = _walk(target, anchor, alpha, helium, base)
    companion = process.fuel_side_j / process.nameplate_j
    for label, share in (
        ("enrichment", ENRICHED_PERCENT / 100.0),
        ("companion", companion),
        ("built efficiency", BUILT_EFFICIENCY),
    ):
        if _leave(share, alpha, helium, five, base) == "unique":
            raise RuntimeError(f"the {label} share became a unique helium hit")
    return DischargeReport(
        anchor=anchor,
        steps=steps,
        target=target,
        discharge_gwd=process.discharge_gwd,
        enriched_percent=ENRICHED_PERCENT,
        fabrication_cap_percent=FABRICATION_CAP_PERCENT,
        companion_ratio=companion,
        built_efficiency=BUILT_EFFICIENCY,
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=process.nuclear_scalar,
        thermodynamics_scalar=process.thermodynamics_scalar,
        gamma_fsot=process.gamma_fsot,
        valve_beta=process.valve_beta,
        pin_prefix=process.pin_prefix,
        digest_prefix=digest[:16],
    )


def _step_lines(step: DischargeStep, closing: bool) -> list[str]:
    lines = [
        f"The dressing is `{step.dressing.name}` at `{_num(step.dressing.value)}`.",
        f"The relative gap is `{_num(step.dressing.relative_gap)}`.",
    ]
    if closing:
        lines.append("The relative gap meets the helium quotient match.")
        lines.append("One dressing of that anchor meets the helium quotient match.")
        lines.append(f"The next dressing sits at relative gap `{_num(step.second_gap)}`.")
    else:
        lines.append("The dressing stays wider than the helium quotient match.")
    lines.append(
        f"The anchor has `{step.expressions}` dressings, median spacing `{_num(step.median_spacing)}`."
    )
    return lines


def render() -> str:
    report = build_discharge()
    lines = [
        "# Fission discharge",
        "",
        "The heavy-metal discharge is the notes figure of 45 GWd per tonne.",
        "Per kilogram it is a share of the ranked pure U-235 density.",
        "The ranked density stays the filed ledger.",
        "The chain stays at k_eff of 1.",
        "",
        "## Share",
        "",
        f"The filed share is `{_num(report.target)}`.",
        f"The five-atom nearest is `{report.anchor.name}` at `{_num(report.anchor.value)}`.",
        f"The anchor gap is `{_num(report.anchor.relative_gap)}`.",
    ]
    for step in report.steps[:-1]:
        lines.extend(_step_lines(step, False))
    lines.extend(_step_lines(report.steps[-1], True))
    lines.extend(
        [
            f"The closing share is `{_num(report.steps[-1].dressing.value)}`.",
            "The 45 GWd/t discharge stays the material balance.",
            "The fabrication cap stays 5 percent.",
            "The 60 GWd/t stack stays unrun.",
            f"The page enrichment stays `{_num(report.enriched_percent)}` percent.",
            "The fuel-year companion stays the page pair.",
            "The plants-being-built efficiency stays the page figure.",
            "The temperature Carnot stays the page temperatures.",
            "This run does not install a seed as a new burn-up or a new plant row.",
            "",
            "## FSOT record",
            "",
            f"Pin `{report.pin_prefix}` matches the loaded engine prefix `{report.digest_prefix}`.",
            "The scalars below are dimensionless.",
            "A domain scalar is not a plant efficiency.",
            "The nuclear scalar stays out of the fission joules.",
            f"Nuclear_Physics: `{_num(report.nuclear_scalar)}`.",
            f"Thermodynamics: `{_num(report.thermodynamics_scalar)}`.",
            f"gamma_FSOT: `{_num(report.gamma_fsot)}`.",
            "gamma_FSOT does not multiply the plant.",
            "The engine valve stays out of the fission joules.",
            f"Engine valve: `{_num(report.valve_beta)}`.",
            f"alpha in the dressing is the fine-structure leaf `{_num(report.alpha)}`.",
            f"The helium quotient match is `{_num(report.helium_mismatch)}`.",
            "Every closing gap meets the helium quotient match.",
            "The process score's foundational seed stays on that score.",
            "",
            "The checked-in schedule is not rewritten.",
            "",
        ]
    )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    del argv
    text = render()
    output = default_output_path()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text, encoding="utf-8")
    print(text, end="" if text.endswith("\n") else "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
