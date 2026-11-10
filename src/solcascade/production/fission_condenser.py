"""Name the condenser shares, the service multiples, and the linear scale.

The shares are the scored cooling loads over the running condenser. The
multiples are the four-hour and twenty-hour loads over the scaled
design-point duty. The linear scale is the cube root of the small-class
share. The unit-interval dressing of that scale drops below 116 values.
The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay_dress import _base_pieces
from solcascade.production.fission_decay_dress import _dressings as unit_dressings
from solcascade.production.fission_decay_heat import _dressings as positive_dressings
from solcascade.production.fission_decay_heat import _nearest_row, _positive_nearest, _positive_size
from solcascade.production.fission_decay_seed import NamedGap, _five_pieces, _nearest, _size
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.fission_removal import build_removal
from solcascade.production.fission_safety import build_safety
from solcascade.production.fission_service import build_service
from solcascade.production.fission_stress import build_stress
from solcascade.production.refine import project_root

_COUNT = 116
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4
_SECOND_ANCHOR = "e^-3+gamma^6"
_SECOND_STEPS = (
    "(e^-3+gamma^6)*(1+alpha*gamma^3)",
    "((e^-3+gamma^6)*(1+alpha*gamma^3))*(1+alpha*pi^-5)",
)
_DAY_ANCHOR = "pi^-4*sqrt(G)"
_DAY_STEPS = (
    "(pi^-4*sqrt(G))*(1+alpha*G^4)",
    "((pi^-4*sqrt(G))*(1+alpha*G^4))*(1-alpha*e^-4)",
    "(((pi^-4*sqrt(G))*(1+alpha*G^4))*(1-alpha*e^-4))*(1+alpha*e^-6)",
)
_TWO_ANCHOR = "pi^-4*G^3"
_TWO_STEPS = (
    "(pi^-4*G^3)*(1+alpha*e^-4)",
    "((pi^-4*G^3)*(1+alpha*e^-4))*(1+alpha*e^-6)",
)
_FOUR_ANCHOR = "e+phi^-3"
_FOUR_STEPS = ("(e+phi^-3)*(1+alpha*e^-3)",)
_TWENTY_ANCHOR = "G+sqrt(G)"
_TWENTY_STEPS = (
    "(G+sqrt(G))*(1-alpha*phi^-5)",
    "((G+sqrt(G))*(1-alpha*phi^-5))*(1+alpha*pi^-6)",
)
_SCALE_ANCHOR = "phi^-5+gamma"
_SCALE_STEPS = (
    "(phi^-5+gamma)*(1+alpha*phi^-2)",
    "((phi^-5+gamma)*(1+alpha*phi^-2))*(1+alpha*gamma^6)",
    "(((phi^-5+gamma)*(1+alpha*phi^-2))*(1+alpha*gamma^6))*(1+alpha*pi^-6)",
)


@dataclass(frozen=True)
class CondenserStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class CondenserFigure:
    anchor: NamedGap
    steps: tuple[CondenserStep, ...]
    target: float


@dataclass(frozen=True)
class CondenserReport:
    one_second: CondenserFigure
    one_day: CondenserFigure
    two_day: CondenserFigure
    four_hour: CondenserFigure
    twenty_hour: CondenserFigure
    linear_scale: CondenserFigure
    condenser_w: float
    small_condenser_w: float
    plant_duty_w: float
    size_ratio: float
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_condenser.md"


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


def _share(stress, label: str) -> float:
    large = next(case.share for case in stress.shutdown if case.name == f"large {label}")
    small = next(case.share for case in stress.shutdown if case.name == f"small {label}")
    if _gap(large, small) > 1.0e-12:
        raise RuntimeError("the small-class condenser share left the large share")
    return large


def _walk(target, anchor_name, expected, alpha, helium, five, base, positive: bool) -> CondenserFigure:
    if positive:
        name, value, anchor_gap = _positive_nearest(five, target)
        anchor = NamedGap(name, value, anchor_gap, anchor_gap <= _REVIEW_GAP)
    else:
        anchor = _nearest(five, target)
        name, value = anchor.name, anchor.value
    if anchor.name != anchor_name:
        raise RuntimeError(f"the condenser anchor changed to {anchor.name}")
    if anchor.relative_gap <= helium:
        raise RuntimeError("the condenser anchor already meets helium")
    if positive and not target > 1.0:
        dropped = len(unit_dressings(anchor.name, anchor.value, alpha, base))
        if dropped >= _COUNT:
            raise RuntimeError("the unit grammar still has 116 dressings")
    steps: list[CondenserStep] = []
    for step_name in expected:
        rows = (positive_dressings if positive else unit_dressings)(name, value, alpha, base)
        if positive:
            hit_name, hit_value, hit_gap = _nearest_row(rows, target)
            hit = NamedGap(hit_name, hit_value, hit_gap, hit_gap <= _REVIEW_GAP)
        else:
            hit = _nearest(rows, target)
        if hit.name != step_name:
            raise RuntimeError(f"the condenser dressing changed to {hit.name}")
        if hit.relative_gap >= _gap(value, target):
            raise RuntimeError("the condenser dressing did not tighten")
        size = (_positive_size if positive else _size)(rows)
        if size.expressions != _COUNT or size.inside != _COUNT:
            raise RuntimeError("the condenser dressing count changed")
        if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
            raise RuntimeError("the condenser spacing is not the base-piece spacing")
        inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
        if inside_helium > 1:
            raise RuntimeError("several condenser dressings meet the helium standard")
        steps.append(
            CondenserStep(
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
        raise RuntimeError("the closing condenser figure misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError("a second condenser dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier condenser figure already meets helium")
    return CondenserFigure(anchor, tuple(steps), target)


def build_condenser() -> CondenserReport:
    stress = build_stress()
    safety = build_safety()
    service = build_service()
    removal = build_removal()
    module, digest = load_engine()
    if digest[:16] != stress.pin_prefix or stress.digest_prefix != stress.pin_prefix:
        raise RuntimeError("the condenser pin does not match the plant pin")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    if abs(helium - stress.helium_mismatch) > 0.0:
        raise RuntimeError("the helium standard changed")
    five = _five_pieces(module)
    base = _base_pieces(module)
    one_second = _share(stress, "one-second fraction")
    one_day = _share(stress, "one-day fraction")
    two_day = _share(stress, "two-day fraction")
    if not 0.0 < two_day < one_day < one_second < 1.0:
        raise RuntimeError("the condenser shares left their order")
    four_hour = service.service_w / removal.plant_duty_w
    four_small = service.service_small_w / removal.small_duty_w
    twenty_hour = service.cooldown_w / removal.plant_duty_w
    twenty_small = service.cooldown_small_w / removal.small_duty_w
    if _gap(four_hour, four_small) > 1.0e-12 or _gap(twenty_hour, twenty_small) > 1.0e-12:
        raise RuntimeError("the small-class service multiple left the large multiple")
    if not four_hour > twenty_hour > 1.0:
        raise RuntimeError("a service load sits inside the design-point duty")
    scale = safety.linear_scale
    if safety.size_ratio != 0.3 or scale != safety.size_ratio ** (1.0 / 3.0):
        raise RuntimeError("the linear scale left the cube root of the small-class share")
    if not 0.0 < scale < 1.0:
        raise RuntimeError("the linear scale left the unit interval")
    return CondenserReport(
        one_second=_walk(one_second, _SECOND_ANCHOR, _SECOND_STEPS, alpha, helium, five, base, False),
        one_day=_walk(one_day, _DAY_ANCHOR, _DAY_STEPS, alpha, helium, five, base, False),
        two_day=_walk(two_day, _TWO_ANCHOR, _TWO_STEPS, alpha, helium, five, base, False),
        four_hour=_walk(four_hour, _FOUR_ANCHOR, _FOUR_STEPS, alpha, helium, five, base, True),
        twenty_hour=_walk(twenty_hour, _TWENTY_ANCHOR, _TWENTY_STEPS, alpha, helium, five, base, True),
        linear_scale=_walk(scale, _SCALE_ANCHOR, _SCALE_STEPS, alpha, helium, five, base, True),
        condenser_w=stress.power[0].waste_w,
        small_condenser_w=stress.power[1].waste_w,
        plant_duty_w=removal.plant_duty_w,
        size_ratio=safety.size_ratio,
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=stress.nuclear_scalar,
        thermodynamics_scalar=stress.thermodynamics_scalar,
        gamma_fsot=stress.gamma_fsot,
        valve_beta=stress.valve_beta,
        pin_prefix=stress.pin_prefix,
        digest_prefix=stress.digest_prefix,
    )


def _step_lines(step: CondenserStep, closing: bool) -> list[str]:
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


def _figure_lines(figure: CondenserFigure, lead: str, role: str) -> list[str]:
    lines = [
        lead,
        f"The five-atom nearest is `{figure.anchor.name}` at `{_num(figure.anchor.value)}`.",
        f"The anchor gap is `{_num(figure.anchor.relative_gap)}`.",
    ]
    if len(figure.steps) == 1:
        lines.append("The first factor closes that quotient.")
    for step in figure.steps[:-1]:
        lines.extend(_step_lines(step, False))
    lines.extend(_step_lines(figure.steps[-1], True))
    closing = figure.steps[-1].dressing.value
    if _num(closing) == _num(figure.target):
        lines.append(f"The closing value agrees with `{_num(figure.target)}` at six figures.")
    lines.append("The closing value names that quotient.")
    lines.append(role)
    return lines


def render() -> str:
    report = build_condenser()
    lines = [
        "# Fission condenser",
        "",
        "The stress score prints three cooling loads as shares of the running condenser.",
        "The service score prints the four-hour and twenty-hour loads.",
        "The safety score prints the linear scale of a core that keeps its shape.",
        "FSOT names those quotients.",
        "The chain stays at k_eff of 1.",
        "",
        "## Condenser",
        "",
        f"The running condenser duty is `{_num(report.condenser_w)}` W.",
        f"The small-class condenser duty is `{_num(report.small_condenser_w)}` W.",
        "The small-class shares match the large-class shares.",
    ]
    lines.extend(
        _figure_lines(
            report.one_second,
            f"The one-second load is `{_num(report.one_second.target)}` of that condenser duty.",
            "The operating share stays the one-second load over the running condenser.",
        )
    )
    lines.extend(
        _figure_lines(
            report.one_day,
            f"The one-day load is `{_num(report.one_day.target)}` of that condenser duty.",
            "The operating share stays the one-day load over the running condenser.",
        )
    )
    lines.extend(
        _figure_lines(
            report.two_day,
            f"The two-day load is `{_num(report.two_day.target)}` of that condenser duty.",
            "The operating share stays the two-day load over the running condenser.",
        )
    )
    lines.extend(
        [
            "Decay heat is a cooling load.",
            "It is not a second electricity product.",
            "",
            "## Service",
            "",
            f"The scaled design-point duty is `{_num(report.plant_duty_w)}` W.",
            "The small-class multiples match the large-class multiples.",
            "Both service loads stay above that duty.",
        ]
    )
    lines.extend(
        _figure_lines(
            report.four_hour,
            f"The four-hour load is `{_num(report.four_hour.target)}` times that duty.",
            "The operating multiple stays the four-hour load over the scaled design-point duty.",
        )
    )
    lines.extend(
        _figure_lines(
            report.twenty_hour,
            f"The twenty-hour load is `{_num(report.twenty_hour.target)}` times that duty.",
            "The operating multiple stays the twenty-hour load over the scaled design-point duty.",
        )
    )
    lines.extend(
        [
            "The UFSAR does not print a watt rate at the 350°F entry.",
            "This score assigns that hotter rate no figure.",
            "A loss of condenser circulation has no second watt capacity on file.",
            "This score assigns that loss no figure.",
            "",
            "## Scale",
            "",
            f"The small-class share of the nameplate stays `{_num(report.size_ratio)}`.",
            "The linear scale is the cube root of that share.",
            "The unit-interval dressing of that scale drops below 116 values.",
            "The positive dressing names it.",
        ]
    )
    lines.extend(
        _figure_lines(
            report.linear_scale,
            f"The linear scale is `{_num(report.linear_scale.target)}`.",
            "The operating scale stays the cube root of the small-class share.",
        )
    )
    lines.extend(
        [
            "The 60 GWd/t stack stays unrun.",
            "The fabrication cap stays 5 percent.",
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
            "This run does not install a seed as a new enrichment or a new plant row.",
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
