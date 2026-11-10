"""Name the filed bounds that still had no helium-tight expression.

The catalog ceiling stays above the temperature Carnot and is not run.
The fast row keeps its sodium barrier. The 1952 fit remains the witness
for its stated bands. The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay import _LEAD, stated_band
from solcascade.production.fission_decay_dress import _base_pieces
from solcascade.production.fission_decay_dress import _dressings as unit_dressings
from solcascade.production.fission_decay_heat import _dressings as positive_dressings
from solcascade.production.fission_decay_heat import _nearest_row, _positive_nearest, _positive_size
from solcascade.production.fission_decay_seed import NamedGap, _five_pieces, _nearest, _size
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.fission_process import (
    FABRICATION_CAP_PERCENT,
    NAMEPLATE_MWE,
    build_process,
)
from solcascade.production.fission_safety import SMALL_CLASS_MWE
from solcascade.production.plants import FAST_ID, THERMAL_ID, build_plants
from solcascade.production.refine import catalog_record, project_root

_COUNT = 116
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4
_CEILING_ANCHOR = "sqrt(gamma)*G^6"
_CEILING_STEPS = (
    "(sqrt(gamma)*G^6)*(1+alpha*phi^-2)",
    "((sqrt(gamma)*G^6)*(1+alpha*phi^-2))*(1+alpha*e^-4)",
)
_HEAD_ANCHOR = "sqrt(phi)+phi^-5"
_HEAD_STEPS = ("(sqrt(phi)+phi^-5)*(1+alpha*phi^-4)",)
_FAST_EFF_ANCHOR = "phi^-3*G^-6"
_FAST_EFF_STEPS = (
    "(phi^-3*G^-6)*(1+alpha*phi^-5)",
    "((phi^-3*G^-6)*(1+alpha*phi^-5))*(1+alpha*pi^-5)",
    "(((phi^-3*G^-6)*(1+alpha*phi^-5))*(1+alpha*pi^-5))*(1+alpha*pi^-6)",
)
_FAST_CEIL_ANCHOR = "sqrt(pi)-sqrt(phi)"
_FAST_CEIL_STEPS = (
    "(sqrt(pi)-sqrt(phi))*(1-alpha*gamma^4)",
    "((sqrt(pi)-sqrt(phi))*(1-alpha*gamma^4))*(1-alpha*e^-5)",
    "(((sqrt(pi)-sqrt(phi))*(1-alpha*gamma^4))*(1-alpha*e^-5))*(1-alpha*pi^-6)",
)
_FAST_HEAD_ANCHOR = "gamma^2+G"
_FAST_HEAD_STEPS = (
    "(gamma^2+G)*(1+alpha*phi^-5)",
    "((gamma^2+G)*(1+alpha*phi^-5))*(1+alpha*pi^-5)",
)
_SMALL_ANCHOR = "phi^-3+gamma^5"
_SMALL_STEPS = (
    "(phi^-3+gamma^5)*(1-alpha*gamma^5)",
    "((phi^-3+gamma^5)*(1-alpha*gamma^5))*(1-alpha*pi^-6)",
)
_BAND_ANCHOR = "e-phi^2"
_BAND_STEPS = (
    "(e-phi^2)*(1-alpha*gamma^2)",
    "((e-phi^2)*(1-alpha*gamma^2))*(1-alpha*e^-5)",
    "(((e-phi^2)*(1-alpha*gamma^2))*(1-alpha*e^-5))*(1+alpha*pi^-6)",
)
_THERMAL_RETURN_ANCHOR = "pi^6/sqrt(phi)"
_THERMAL_RETURN_STEPS = (
    "(pi^6/sqrt(phi))*(1-alpha*gamma)",
    "((pi^6/sqrt(phi))*(1-alpha*gamma))*(1-alpha*e^-5)",
    "(((pi^6/sqrt(phi))*(1-alpha*gamma))*(1-alpha*e^-5))*(1+alpha*pi^-6)",
)
_FAST_RETURN_ANCHOR = "pi^6-e^4"
_FAST_RETURN_STEPS = (
    "(pi^6-e^4)*(1-alpha*G^-5)",
    "((pi^6-e^4)*(1-alpha*G^-5))*(1+alpha*gamma^5)",
    "(((pi^6-e^4)*(1-alpha*G^-5))*(1+alpha*gamma^5))*(1-alpha*pi^-5)",
)
_STEP_ANCHOR = "e^-4+G^-2"
_STEP_STEPS = (
    "(e^-4+G^-2)*(1+alpha*phi^-3)",
    "((e^-4+G^-2)*(1+alpha*phi^-3))*(1-alpha*e^-4)",
    "(((e^-4+G^-2)*(1+alpha*phi^-3))*(1-alpha*e^-4))*(1-alpha*e^-6)",
)


@dataclass(frozen=True)
class BoundStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class BoundFigure:
    anchor: NamedGap
    steps: tuple[BoundStep, ...]
    target: float


@dataclass(frozen=True)
class BoundsReport:
    ceiling: BoundFigure
    headroom: BoundFigure
    fast_efficiency: BoundFigure
    fast_ceiling: BoundFigure
    fast_headroom: BoundFigure
    small_class: BoundFigure
    late_band: BoundFigure
    thermal_return: BoundFigure
    fast_return: BoundFigure
    efficiency_step: BoundFigure
    early_band: float
    middle_band: float
    carnot: float
    demonstrated_efficiency: float
    fabrication_cap_percent: float
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_bounds.md"


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


def _walk(target, anchor_name, expected, alpha, helium, five, base, positive: bool) -> BoundFigure:
    if positive:
        name, value, anchor_gap = _positive_nearest(five, target)
        anchor = NamedGap(name, value, anchor_gap, anchor_gap <= _REVIEW_GAP)
    else:
        anchor = _nearest(five, target)
        name, value = anchor.name, anchor.value
    if anchor.name != anchor_name:
        raise RuntimeError("the bounds anchor changed")
    steps: list[BoundStep] = []
    for step_name in expected:
        rows = (positive_dressings if positive else unit_dressings)(name, value, alpha, base)
        if positive:
            hit_name, hit_value, hit_gap = _nearest_row(rows, target)
            hit = NamedGap(hit_name, hit_value, hit_gap, hit_gap <= _REVIEW_GAP)
        else:
            hit = _nearest(rows, target)
        if hit.name != step_name:
            raise RuntimeError("the bounds dressing changed")
        if hit.relative_gap >= _gap(value, target):
            raise RuntimeError("the bounds dressing did not tighten")
        size = (_positive_size if positive else _size)(rows)
        if size.expressions != _COUNT or size.inside != _COUNT:
            raise RuntimeError("the bounds dressing count changed")
        if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
            raise RuntimeError("the bounds spacing is not the base-piece spacing")
        inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
        if inside_helium > 1:
            raise RuntimeError("several bounds dressings meet the helium standard")
        steps.append(
            BoundStep(
                dressing=hit,
                inside_review=sum(_gap(number, target) <= _REVIEW_GAP for _row, number in rows),
                inside_tenth=sum(_gap(number, target) <= _TENTH for _row, number in rows),
                inside_helium=inside_helium,
                second_gap=_second_gap(rows, target, hit.relative_gap),
                expressions=size.expressions,
                median_spacing=size.median_spacing,
            )
        )
        name, value = hit.name, hit.value
    closing = steps[-1]
    if closing.inside_helium != 1 or closing.dressing.relative_gap > helium:
        raise RuntimeError("the closing bounds figure misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError("a second bounds dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier bounds figure already meets helium")
    return BoundFigure(anchor, tuple(steps), target)


def build_bounds() -> BoundsReport:
    plants = build_plants()
    process = build_process()
    module, digest = load_engine()
    if digest[:16] != plants.pin_prefix or plants.pin_prefix != process.pin_prefix:
        raise RuntimeError("the bounds pin does not match the plant pin")
    thermal = catalog_record(THERMAL_ID)
    fast = catalog_record(FAST_ID)
    if "k_eff is held at 1" not in str(thermal["governing_equation"]):
        raise RuntimeError("the governing equation does not hold k_eff at 1")
    if "Sodium" not in str(fast.get("engineering_barrier") or ""):
        raise RuntimeError("the fast barrier does not name sodium")
    if FABRICATION_CAP_PERCENT != 5.0:
        raise RuntimeError("the fabrication cap left 5 percent")
    if plants.burnup_demonstrated != 1.0 or plants.burnup_cap != 1.0:
        raise RuntimeError("the phase-two burnup reading left one")
    thermal_plant = plants.thermal
    fast_plant = plants.fast
    if not thermal_plant.efficiency < process.carnot < thermal_plant.ceiling:
        raise RuntimeError("the catalog ceiling is not above the temperature Carnot")
    if not fast_plant.efficiency < fast_plant.ceiling:
        raise RuntimeError("the fast ceiling is not above the fast efficiency")
    small = SMALL_CLASS_MWE / NAMEPLATE_MWE
    early = stated_band(1.0)
    middle = stated_band(100.0)
    late = stated_band(1.0e4)
    if early != fast_plant.ceiling or middle != small or late != 0.10:
        raise RuntimeError("a stated band left its filed number")
    if _LEAD != late:
        raise RuntimeError("the equation (2) lead left the ten percent band")
    if _gap(thermal_plant.headroom, thermal_plant.ceiling / thermal_plant.efficiency) > 1.0e-12:
        raise RuntimeError("the thermal headroom left the ceiling quotient")
    if _gap(fast_plant.headroom, fast_plant.ceiling / fast_plant.efficiency) > 1.0e-12:
        raise RuntimeError("the fast headroom left the ceiling quotient")
    if _gap(plants.efficiency_step, fast_plant.efficiency / thermal_plant.efficiency) > 1.0e-12:
        raise RuntimeError("the efficiency step left the fast-over-thermal quotient")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    five = _five_pieces(module)
    base = _base_pieces(module)
    ceiling = _walk(thermal_plant.ceiling, _CEILING_ANCHOR, _CEILING_STEPS, alpha, helium, five, base, True)
    headroom = _walk(thermal_plant.headroom, _HEAD_ANCHOR, _HEAD_STEPS, alpha, helium, five, base, True)
    fast_efficiency = _walk(fast_plant.efficiency, _FAST_EFF_ANCHOR, _FAST_EFF_STEPS, alpha, helium, five, base, True)
    fast_ceiling = _walk(fast_plant.ceiling, _FAST_CEIL_ANCHOR, _FAST_CEIL_STEPS, alpha, helium, five, base, True)
    fast_headroom = _walk(fast_plant.headroom, _FAST_HEAD_ANCHOR, _FAST_HEAD_STEPS, alpha, helium, five, base, True)
    small_class = _walk(small, _SMALL_ANCHOR, _SMALL_STEPS, alpha, helium, five, base, True)
    late_band = _walk(late, _BAND_ANCHOR, _BAND_STEPS, alpha, helium, five, base, False)
    thermal_return = _walk(
        thermal_plant.output_over_input, _THERMAL_RETURN_ANCHOR, _THERMAL_RETURN_STEPS, alpha, helium, five, base, True
    )
    fast_return = _walk(
        fast_plant.output_over_input, _FAST_RETURN_ANCHOR, _FAST_RETURN_STEPS, alpha, helium, five, base, True
    )
    efficiency_step = _walk(plants.efficiency_step, _STEP_ANCHOR, _STEP_STEPS, alpha, helium, five, base, True)
    if not process.carnot < ceiling.steps[-1].dressing.value:
        raise RuntimeError("the named catalog ceiling fell to the temperature Carnot")
    if not fast_efficiency.steps[-1].dressing.value < fast_ceiling.steps[-1].dressing.value:
        raise RuntimeError("the named fast efficiency met the fast ceiling")
    if _gap(headroom.steps[-1].dressing.value * thermal_plant.efficiency, thermal_plant.ceiling) > helium:
        raise RuntimeError("the named headroom left the ceiling quotient")
    if _gap(small_class.steps[-1].dressing.value * NAMEPLATE_MWE, SMALL_CLASS_MWE) > helium:
        raise RuntimeError("the named small class left the 300 MWe boundary")
    return BoundsReport(
        ceiling=ceiling,
        headroom=headroom,
        fast_efficiency=fast_efficiency,
        fast_ceiling=fast_ceiling,
        fast_headroom=fast_headroom,
        small_class=small_class,
        late_band=late_band,
        thermal_return=thermal_return,
        fast_return=fast_return,
        efficiency_step=efficiency_step,
        early_band=early,
        middle_band=middle,
        carnot=process.carnot,
        demonstrated_efficiency=thermal_plant.efficiency,
        fabrication_cap_percent=FABRICATION_CAP_PERCENT,
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=plants.nuclear_scalar,
        thermodynamics_scalar=plants.thermodynamics_scalar,
        gamma_fsot=plants.gamma_fsot,
        valve_beta=plants.valve_beta,
        pin_prefix=plants.pin_prefix,
        digest_prefix=digest[:16],
    )


def _step_lines(step: BoundStep, closing: bool) -> list[str]:
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


def _figure_lines(title: str, intro: list[str], figure: BoundFigure, nearest: str) -> list[str]:
    lines = ["", f"## {title}", "", *intro]
    lines.append(f"The {nearest} five-atom nearest is `{figure.anchor.name}` at `{_num(figure.anchor.value)}`.")
    lines.append(f"The anchor gap is `{_num(figure.anchor.relative_gap)}`.")
    for step in figure.steps[:-1]:
        lines.extend(_step_lines(step, False))
    lines.extend(_step_lines(figure.steps[-1], True))
    return lines


def render() -> str:
    report = build_bounds()
    lines = [
        "# Fission bounds",
        "",
        "The catalog ceiling stays above the temperature Carnot and is not run.",
        "The fast row keeps its sodium barrier.",
        "The demonstrated thermal efficiency stays the heat-ratio score.",
        "Ratios that leave the unit interval use the positive dressing from the heat-ratio score.",
        "The chain stays at k_eff of 1.",
    ]
    lines.extend(
        _figure_lines(
            "Catalog ceiling",
            [f"The filed ceiling is `{_num(report.ceiling.target)}`."],
            report.ceiling,
            "positive",
        )
    )
    lines.extend(
        _figure_lines(
            "Thermal headroom",
            [
                f"The ceiling divided by the demonstrated efficiency is `{_num(report.headroom.target)}`.",
            ],
            report.headroom,
            "positive",
        )
    )
    lines.extend(
        _figure_lines(
            "Fast efficiency",
            [f"The filed fast efficiency is `{_num(report.fast_efficiency.target)}`."],
            report.fast_efficiency,
            "positive",
        )
    )
    lines.extend(
        _figure_lines(
            "Fast ceiling",
            [
                f"The filed fast ceiling is `{_num(report.fast_ceiling.target)}`.",
                "The 1952 band under 100 seconds is that same number.",
                "The fast ceiling source states no temperature pair, so this score does not tighten that row.",
            ],
            report.fast_ceiling,
            "positive",
        )
    )
    lines.extend(
        _figure_lines(
            "Fast headroom",
            [f"The fast ceiling divided by the fast efficiency is `{_num(report.fast_headroom.target)}`."],
            report.fast_headroom,
            "positive",
        )
    )
    lines.extend(
        _figure_lines(
            "Small class",
            [
                f"The 300 MWe boundary over the 1000 MWe nameplate is `{_num(report.small_class.target)}`.",
                "The 1952 band from 100 seconds to 10,000 seconds is that same number.",
                "The band stays an uncertainty on the fit.",
                "The small class stays the published electric boundary.",
            ],
            report.small_class,
            "positive",
        )
    )
    lines.extend(
        _figure_lines(
            "Ten percent band",
            [
                f"The 1952 band from 10,000 seconds to 1,000,000 seconds is `{_num(report.late_band.target)}`.",
                "Four hours and twenty hours sit in this band.",
                "The equation (2) lead coefficient is that same number.",
                "The 1952 fit remains the witness.",
            ],
            report.late_band,
            "unit-interval",
        )
    )
    lines.extend(
        _figure_lines(
            "Thermal output over input",
            [
                f"The ranked thermal output over the filed enrichment input is `{_num(report.thermal_return.target)}`.",
                "This ratio compares those two filed energies.",
            ],
            report.thermal_return,
            "positive",
        )
    )
    lines.extend(
        _figure_lines(
            "Fast output over input",
            [
                f"The ranked fast output over the filed enrichment input is `{_num(report.fast_return.target)}`.",
                "This ratio compares those two filed energies.",
            ],
            report.fast_return,
            "positive",
        )
    )
    lines.extend(
        _figure_lines(
            "Efficiency step",
            [
                f"The fast efficiency divided by the thermal efficiency is `{_num(report.efficiency_step.target)}`.",
            ],
            report.efficiency_step,
            "positive",
        )
    )
    lines.extend(
        [
            "The phase-two burnup readings stay at one.",
            "The fabrication cap stays 5 percent.",
            "This run does not install a seed as a new efficiency or a new plant row.",
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
