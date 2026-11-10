"""Derive the decay operating heat from the FSOT heat ratio.

The ratio is operating heat per watt of filed nameplate. The bare value is
the positive five-atom nearest. Two base-piece factors close it. The
nameplate stays the published electric boundary. The 1952 fit stays the
calorimetry witness. The checked-in schedule is not rewritten.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay import build_decay
from solcascade.production.fission_decay_dress import _base_pieces
from solcascade.production.fission_decay_load import build_decay_load
from solcascade.production.fission_decay_seed import (
    GrammarSize,
    NamedGap,
    _REVIEW_GAP,
    _ROUND,
    _five_pieces,
)
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.fission_process import NAMEPLATE_MWE
from solcascade.production.fission_safety import SMALL_CLASS_MWE
from solcascade.production.refine import project_root

_ANCHOR = "pi-gamma^4"
_STEPS = (
    "(pi-gamma^4)*(1-alpha*pi^-4)",
    "((pi-gamma^4)*(1-alpha*pi^-4))*(1-alpha*e^-6)",
)
_COUNT = 116
_SPACING = "0.000511364"
_TENTH = 1.0e-4


@dataclass(frozen=True)
class HeatStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    size: GrammarSize


@dataclass(frozen=True)
class HeatLoad:
    label: str
    expression: str
    watts: float
    fit_watts: float
    small_watts: float
    fit_small_watts: float
    relative_gap: float


@dataclass(frozen=True)
class HeatReport:
    anchor_name: str
    anchor_value: float
    anchor_gap: float
    steps: tuple[HeatStep, ...]
    ratio: float
    filed_ratio: float
    efficiency: float
    filed_efficiency: float
    thermal_w: float
    filed_thermal_w: float
    small_thermal_w: float
    filed_small_thermal_w: float
    loads: tuple[HeatLoad, ...]
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_decay_heat.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _gap(derived: float, witness: float) -> float:
    if witness == 0.0:
        raise RuntimeError("the witness is zero")
    return abs(derived - witness) / abs(witness)


def _positive_nearest(pieces, target: float) -> tuple[str, float, float]:
    best: tuple[float, str, float] | None = None
    for name, value in pieces:
        number = float(value)
        if not math.isfinite(number) or number <= 0.0:
            continue
        gap = abs(number - target) / target
        if best is None or (gap, name) < (best[0], best[1]):
            best = (gap, name, number)
    if best is None:
        raise RuntimeError("the grammar has no positive value")
    gap, name, number = best
    return name, number, gap


def _dressings(anchor_name: str, anchor: float, alpha: float, base) -> tuple[tuple[str, float], ...]:
    rows = []
    for name, value in base:
        for sign in (1.0, -1.0):
            number = anchor * (1.0 + sign * alpha * float(value))
            if math.isfinite(number) and number > 0.0:
                mark = "-" if sign < 0 else "+"
                rows.append((f"({anchor_name})*(1{mark}alpha*{name})", number))
    if not rows:
        raise RuntimeError("the heat dressing is empty")
    return tuple(rows)


def _positive_size(rows) -> GrammarSize:
    vals = sorted({round(float(value), _ROUND) for _name, value in rows if float(value) > 0.0})
    spacings = [(right - left) / right for left, right in zip(vals, vals[1:]) if right > left]
    if not spacings:
        raise RuntimeError("the heat grammar has no spacing")
    spacings.sort()
    return GrammarSize(len(rows), len(vals), spacings[len(spacings) // 2])


def _nearest_row(rows, target: float) -> tuple[str, float, float]:
    best: tuple[float, str, float] | None = None
    for name, value in rows:
        gap = abs(value - target) / target
        if best is None or (gap, name) < (best[0], best[1]):
            best = (gap, name, value)
    if best is None:
        raise RuntimeError("the heat dressing is empty")
    gap, name, value = best
    return name, value, gap


def build_decay_heat() -> HeatReport:
    decay = build_decay()
    loads = build_decay_load()
    module, _digest = load_engine()
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    nameplate_w = NAMEPLATE_MWE * 1.0e6
    small_w = SMALL_CLASS_MWE * 1.0e6
    if not 0.0 < small_w < nameplate_w:
        raise RuntimeError("the small class is not under the nameplate")
    filed_ratio = decay.thermal_w / nameplate_w
    if abs(filed_ratio - decay.small_thermal_w / small_w) > 1.0e-12:
        raise RuntimeError("the two filed heats do not share one ratio")
    anchor_name, anchor_value, anchor_gap = _positive_nearest(_five_pieces(module), filed_ratio)
    if anchor_name != _ANCHOR:
        raise RuntimeError("the heat-ratio anchor changed")
    base = _base_pieces(module)
    name = anchor_name
    value = anchor_value
    steps: list[HeatStep] = []
    for expected in _STEPS:
        rows = _dressings(name, value, alpha, base)
        hit_name, hit_value, hit_gap = _nearest_row(rows, filed_ratio)
        if hit_name != expected:
            raise RuntimeError("the heat-ratio dressing changed")
        if hit_gap >= _gap(value, filed_ratio):
            raise RuntimeError("the heat-ratio dressing did not tighten")
        size = _positive_size(rows)
        if size.expressions != _COUNT or size.inside != _COUNT:
            raise RuntimeError("the heat-ratio dressing count changed")
        if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
            raise RuntimeError("the heat-ratio spacing is not the base-piece spacing")
        inside_review = sum(_gap(number, filed_ratio) <= _REVIEW_GAP for _row, number in rows)
        inside_tenth = sum(_gap(number, filed_ratio) <= _TENTH for _row, number in rows)
        inside_helium = sum(_gap(number, filed_ratio) <= helium for _row, number in rows)
        if inside_helium > 1:
            raise RuntimeError("several heat ratios meet the helium standard")
        steps.append(
            HeatStep(
                dressing=NamedGap(hit_name, hit_value, hit_gap, hit_gap <= _REVIEW_GAP),
                inside_review=inside_review,
                inside_tenth=inside_tenth,
                inside_helium=inside_helium,
                size=size,
            )
        )
        name = hit_name
        value = hit_value
    closing = steps[-1]
    if closing.inside_helium != 1 or closing.dressing.relative_gap > helium:
        raise RuntimeError("the heat ratio does not meet the helium standard")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier heat ratio already meets helium")
    thermal_w = nameplate_w * closing.dressing.value
    small_thermal_w = small_w * closing.dressing.value
    if _gap(thermal_w, decay.thermal_w) > helium or _gap(small_thermal_w, decay.small_thermal_w) > helium:
        raise RuntimeError("the derived operating heat misses the helium match")
    efficiency = 1.0 / closing.dressing.value
    filed_efficiency = nameplate_w / decay.thermal_w
    if _gap(efficiency, filed_efficiency) > helium:
        raise RuntimeError("the implied efficiency misses the helium match")
    derived_loads: list[HeatLoad] = []
    for load in loads.loads:
        watts = load.fraction.value * thermal_w
        small = load.fraction.value * small_thermal_w
        gap = _gap(watts, load.fit_watts)
        if gap > helium or _gap(small, load.fit_small_watts) > helium:
            raise RuntimeError(f"the {load.label} misses the helium match on the derived heat")
        derived_loads.append(
            HeatLoad(
                label=load.label,
                expression=load.fraction.name,
                watts=watts,
                fit_watts=load.fit_watts,
                small_watts=small,
                fit_small_watts=load.fit_small_watts,
                relative_gap=gap,
            )
        )
    return HeatReport(
        anchor_name=anchor_name,
        anchor_value=anchor_value,
        anchor_gap=anchor_gap,
        steps=tuple(steps),
        ratio=closing.dressing.value,
        filed_ratio=filed_ratio,
        efficiency=efficiency,
        filed_efficiency=filed_efficiency,
        thermal_w=thermal_w,
        filed_thermal_w=decay.thermal_w,
        small_thermal_w=small_thermal_w,
        filed_small_thermal_w=decay.small_thermal_w,
        loads=tuple(derived_loads),
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=decay.nuclear_scalar,
        thermodynamics_scalar=decay.thermodynamics_scalar,
        gamma_fsot=decay.gamma_fsot,
        valve_beta=decay.valve_beta,
        pin_prefix=decay.pin_prefix,
        digest_prefix=decay.digest_prefix,
    )


def _heat_lines(report: HeatReport) -> list[str]:
    added, closing = report.steps
    return [
        f"The filed heat ratio is `{_num(report.filed_ratio)}`.",
        f"Its five-atom nearest is `{report.anchor_name}` at `{_num(report.anchor_value)}`.",
        f"The anchor gap is `{_num(report.anchor_gap)}`.",
        (
            f"The added dressing is `{added.dressing.name}` at `{_num(added.dressing.value)}`, "
            f"relative gap `{_num(added.dressing.relative_gap)}`."
        ),
        "The added dressing stays wider than the helium quotient match.",
        (
            f"The closing ratio is `{closing.dressing.name}` at `{_num(closing.dressing.value)}`, "
            f"relative gap `{_num(closing.dressing.relative_gap)}`."
        ),
        "The closing ratio meets the helium quotient match.",
        "One dressing of the closing anchor meets the helium quotient match.",
        f"{closing.inside_review} dressings of the closing anchor meet the review bar.",
        f"{closing.inside_tenth} dressings of the closing anchor are inside one part in ten thousand.",
        (
            f"The closing anchor has `{closing.size.expressions}` dressings, "
            f"`{closing.size.inside}` distinct, median spacing `{_num(closing.size.median_spacing)}`."
        ),
        f"The large nameplate is `{_num(NAMEPLATE_MWE)}` MWe.",
        f"The derived operating heat is `{_num(report.thermal_w)}` W.",
        f"The filed operating heat is `{_num(report.filed_thermal_w)}` W.",
        f"The relative gap is `{_num(_gap(report.thermal_w, report.filed_thermal_w))}`.",
        f"The small-class boundary is `{_num(SMALL_CLASS_MWE)}` MWe.",
        f"The derived small-class operating heat is `{_num(report.small_thermal_w)}` W.",
        f"The filed small-class operating heat is `{_num(report.filed_small_thermal_w)}` W.",
        f"The implied efficiency is `{_num(report.efficiency)}`.",
        f"The filed efficiency is `{_num(report.filed_efficiency)}`.",
        f"The efficiency gap is `{_num(_gap(report.efficiency, report.filed_efficiency))}`.",
        "The implied efficiency meets the helium quotient match.",
    ]


def _load_lines(load: HeatLoad) -> list[str]:
    return [
        f"The {load.label} is `{load.expression}`.",
        f"The derived cooling load is `{_num(load.watts)}` W.",
        f"The 1952 fit load is `{_num(load.fit_watts)}` W.",
        f"The relative gap is `{_num(load.relative_gap)}`.",
        f"The derived small-class load is `{_num(load.small_watts)}` W.",
        f"The 1952 small-class load is `{_num(load.fit_small_watts)}` W.",
    ]


def render() -> str:
    report = build_decay_heat()
    lines = [
        "# Fission decay heat",
        "",
        "The operating heat is the filed nameplate times the FSOT heat ratio.",
        "The heat ratio is operating heat per watt of nameplate.",
        "The bare value is pi minus gamma to the fourth.",
        "Two base-piece factors close that ratio.",
        "A base piece is an atom, its square root, or a power from 2 through 6 in either sign.",
        f"alpha in the dressing is the fine-structure leaf `{_num(report.alpha)}`.",
        f"The helium quotient match is `{_num(report.helium_mismatch)}`.",
        "The nameplate stays the published electric boundary.",
        "The filed efficiency stays the plant witness.",
        "Every scored gap meets the helium quotient match.",
        "",
        "## Heat",
        "",
    ]
    lines.extend(_heat_lines(report))
    lines.extend(["", "## Loads", ""])
    lines.append("Each cooling load is its FSOT fraction times this operating heat.")
    lines.append("")
    for load in report.loads:
        lines.extend(_load_lines(load))
        lines.append("")
    lines.extend(
        [
            "## Ledger",
            "",
            "The derived operating heat is the heat in the cooling load.",
            "Decay heat is a cooling load.",
            "It is not a second electricity product.",
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
