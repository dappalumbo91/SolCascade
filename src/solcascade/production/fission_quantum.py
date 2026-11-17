"""Name the filed 200 MeV ledger the pinned waves do not carry.

The eight wave rows keep their own targets. The positive dressing names the
file quantum. The expression stays off the wave list, and the channel search
stays closed. The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.constants import FISSION_RECOVERABLE_ENERGY_MEV
from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay_dress import _base_pieces
from solcascade.production.fission_decay_heat import (
    _dressings,
    _nearest_row,
    _positive_nearest,
    _positive_size,
)
from solcascade.production.fission_decay_seed import NamedGap, _five_pieces
from solcascade.production.fission_observables import (
    NUCLEAR_NAMES,
    _adopted_alpha,
    _helium_mismatch,
    _readings,
)
from solcascade.production.fission_paths import ADDED_MEV, FRAGMENT_SHARE, LEDGER_MEV
from solcascade.production.fission_process import build_process
from solcascade.production.plants import FAST_ID, THERMAL_ID, build_plants
from solcascade.production.refine import catalog_record, project_root

_ANCHOR = "e^2*gamma^-6"
_STEPS = (
    "(e^2*gamma^-6)*(1+alpha*phi^-4)",
    "((e^2*gamma^-6)*(1+alpha*phi^-4))*(1+alpha*e^-6)",
)
_COUNT = 116
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4


@dataclass(frozen=True)
class QuantumStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class QuantumReport:
    anchor: NamedGap
    steps: tuple[QuantumStep, ...]
    target: float
    wave_targets: tuple[float, ...]
    ranked_density_j: float
    fragment_share: float
    added_mev: float
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_quantum.md"


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


def _walk(target: float, alpha: float, helium: float, five, base) -> tuple[NamedGap, tuple[QuantumStep, ...]]:
    name, value, anchor_gap = _positive_nearest(five, target)
    if name != _ANCHOR:
        raise RuntimeError("the quantum anchor changed")
    anchor = NamedGap(name, value, anchor_gap, anchor_gap <= _REVIEW_GAP)
    steps: list[QuantumStep] = []
    for step_name in _STEPS:
        rows = _dressings(name, value, alpha, base)
        hit_name, hit_value, hit_gap = _nearest_row(rows, target)
        if hit_name != step_name:
            raise RuntimeError("the quantum dressing changed")
        if hit_gap >= _gap(value, target):
            raise RuntimeError("the quantum dressing did not tighten")
        size = _positive_size(rows)
        if size.expressions != _COUNT or size.inside != _COUNT:
            raise RuntimeError("the quantum dressing count changed")
        if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
            raise RuntimeError("the quantum spacing is not the base-piece spacing")
        inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
        if inside_helium > 1:
            raise RuntimeError("several quantum dressings meet the helium standard")
        steps.append(
            QuantumStep(
                dressing=NamedGap(hit_name, hit_value, hit_gap, hit_gap <= _REVIEW_GAP),
                inside_review=sum(_gap(number, target) <= _REVIEW_GAP for _row, number in rows),
                inside_tenth=sum(_gap(number, target) <= _TENTH for _row, number in rows),
                inside_helium=inside_helium,
                second_gap=_second_gap(rows, target, hit_gap),
                expressions=size.expressions,
                median_spacing=size.median_spacing,
            )
        )
        name = hit_name
        value = hit_value
    closing = steps[-1]
    if closing.inside_helium != 1 or closing.dressing.relative_gap > helium:
        raise RuntimeError("the closing quantum misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError("a second quantum dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier quantum dressing already meets helium")
    return anchor, tuple(steps)


def build_quantum() -> QuantumReport:
    if FISSION_RECOVERABLE_ENERGY_MEV != 200.0 or LEDGER_MEV != 200.0:
        raise RuntimeError("the file ledger left 200 MeV")
    if ADDED_MEV != 30.0 or FRAGMENT_SHARE != (LEDGER_MEV - ADDED_MEV) / LEDGER_MEV:
        raise RuntimeError("the fragment share left the round split of the ledger")
    plants = build_plants()
    process = build_process()
    module, digest = load_engine()
    if digest[:16] != process.pin_prefix or plants.pin_prefix != process.pin_prefix:
        raise RuntimeError("the quantum pin does not match the plant pin")
    record = catalog_record(THERMAL_ID)
    if "k_eff is held at 1" not in str(record["governing_equation"]):
        raise RuntimeError("the governing equation does not hold k_eff at 1")
    if "Sodium" not in str(catalog_record(FAST_ID).get("engineering_barrier") or ""):
        raise RuntimeError("the fast barrier does not name sodium")
    if process.ranked_density_j != plants.thermal.energy_density_j:
        raise RuntimeError("the ranked density left the filed ledger")
    readings = _readings(module)
    if tuple(row.name for row in readings) != NUCLEAR_NAMES:
        raise RuntimeError("the wave rows left the pinned names")
    if any("fission" in row.name.lower() for row in readings):
        raise RuntimeError("a wave name carries the ledger")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    for row in readings:
        if _gap(row.computed, LEDGER_MEV) <= helium or _gap(row.wave_target, LEDGER_MEV) <= helium:
            raise RuntimeError(f"{row.name} carries the ledger quantum")
    five = _five_pieces(module)
    base = _base_pieces(module)
    anchor, steps = _walk(LEDGER_MEV, alpha, helium, five, base)
    if anchor.review:
        raise RuntimeError("the quantum anchor entered the review bar")
    return QuantumReport(
        anchor=anchor,
        steps=steps,
        target=LEDGER_MEV,
        wave_targets=tuple(row.wave_target for row in readings),
        ranked_density_j=process.ranked_density_j,
        fragment_share=FRAGMENT_SHARE,
        added_mev=ADDED_MEV,
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=process.nuclear_scalar,
        thermodynamics_scalar=process.thermodynamics_scalar,
        gamma_fsot=process.gamma_fsot,
        valve_beta=process.valve_beta,
        pin_prefix=process.pin_prefix,
        digest_prefix=digest[:16],
        )


def _step_lines(step: QuantumStep, closing: bool) -> list[str]:
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
    report = build_quantum()
    closing = report.steps[-1]
    lines = [
        "# Fission quantum",
        "",
        "The pinned waves keep their own targets.",
        "The file ledger is the plant energy.",
        "The chain stays at k_eff of 1.",
        "",
        "## Wave rows",
        "",
        f"The file value is `{_num(report.target)}` MeV.",
        "No wave result carries that quantum.",
        f"The wave targets are `{_num(report.wave_targets[0])}`, `{_num(report.wave_targets[1])}`, `{_num(report.wave_targets[2])}`, `{_num(report.wave_targets[3])}`, `{_num(report.wave_targets[4])}`, `{_num(report.wave_targets[5])}`, `{_num(report.wave_targets[6])}`, and `{_num(report.wave_targets[7])}`.",
        "The largest binding wave target is the helium-4 row.",
        "The channel search stays closed.",
        "",
        "## Ledger",
        "",
        f"The positive five-atom nearest is `{report.anchor.name}` at `{_num(report.anchor.value)}`.",
        f"The anchor gap is `{_num(report.anchor.relative_gap)}`.",
    ]
    for step in report.steps[:-1]:
        lines.extend(_step_lines(step, False))
    lines.extend(_step_lines(closing, True))
    if _num(closing.dressing.value) == _num(report.target):
        lines.append(f"The closing value agrees with `{_num(report.target)}` at six figures.")
    else:
        lines.append(f"The closing print is `{_num(closing.dressing.value)}`.")
    lines.extend(
        [
            "The closing value names the file quantum.",
            "This expression is not a wave row.",
            f"The ranked density stays `{_num(report.ranked_density_j)}` J per kilogram.",
            "The catalog density stays on the 200 MeV ledger.",
            f"The about `{_num(report.added_mev)}` MeV addendum stays the paths score.",
            "The fragment share stays the round split of this ledger.",
            "The demonstrated plant stays at its filed efficiency.",
            "The fabrication cap stays 5 percent.",
            "The 60 GWd/t stack stays unrun.",
            "The catalog ceiling stays unrun.",
            "The fast row's filed barrier is the sodium system.",
            "This score adds no joules.",
            "The published stamp records this close.",
            "",
            "## FSOT record",
            "",
            f"Pin `{report.pin_prefix}` matches the loaded engine prefix `{report.digest_prefix}`.",
            "The scalars below are dimensionless.",
            "A domain scalar is not a megaelectronvolt.",
            "The nuclear scalar stays out of the fission joules.",
            f"Nuclear_Physics: `{_num(report.nuclear_scalar)}`.",
            f"Thermodynamics: `{_num(report.thermodynamics_scalar)}`.",
            f"gamma_FSOT: `{_num(report.gamma_fsot)}`.",
            "gamma_FSOT does not multiply the plant.",
            "The engine valve stays out of the fission joules.",
            f"Engine valve: `{_num(report.valve_beta)}`.",
            f"alpha in the dressing is the fine-structure leaf `{_num(report.alpha)}`.",
            f"The helium quotient match is `{_num(report.helium_mismatch)}`.",
            "This run does not install a seed as a new megaelectronvolt or a new plant row.",
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
