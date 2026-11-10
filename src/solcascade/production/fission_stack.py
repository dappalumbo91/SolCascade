"""Name the higher-burnup page ratios and leave the dense page figures.

The 60 GWd/t stack is not a filed plant. The fabrication cap stays 5 percent.
The 4.5 percent enrichment and the U-233 count stay page figures because more
than one dressing meets helium. The enrichment fraction meets the helium
quotient match on two neighbors. Both neighbors stay uninstalled. The operating
fraction stays the printed page figure. The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_balance import _DELAYED_PERCENT
from solcascade.production.fission_decay_dress import _base_pieces
from solcascade.production.fission_decay_dress import _dressings as unit_dressings
from solcascade.production.fission_decay_heat import _dressings as positive_dressings
from solcascade.production.fission_decay_heat import _nearest_row, _positive_nearest, _positive_size
from solcascade.production.fission_decay_seed import NamedGap, _five_pieces, _nearest, _size
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.fission_process import (
    ATOMIC_PERCENT,
    ENRICHED_PERCENT,
    FABRICATION_CAP_PERCENT,
    HIGHER_BURNUP_GWD,
    build_process,
)
from solcascade.production.refine import project_root

_COUNT = 116
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4
_ATOMIC_ANCHOR = "pi^-6+gamma^5"
_ATOMIC_STEPS = (
    "(pi^-6+gamma^5)*(1-alpha*phi^-3)",
    "((pi^-6+gamma^5)*(1-alpha*phi^-3))*(1-alpha*e^-5)",
)
_BURN_ANCHOR = "pi^-3+G^-3"
_BURN_STEPS = ("(pi^-3+G^-3)*(1-alpha*e^-4)",)
_STACK_ANCHOR = "gamma^4+G^-6"
_STACK_STEPS = (
    "(gamma^4+G^-6)*(1-alpha*gamma^6)",
    "((gamma^4+G^-6)*(1-alpha*gamma^6))*(1-alpha*e^-6)",
)
_HIGHER_ANCHOR = "e^-3*sqrt(phi)"
_HIGHER_STEPS = (
    "(e^-3*sqrt(phi))*(1-alpha*phi^-2)",
    "((e^-3*sqrt(phi))*(1-alpha*phi^-2))*(1-alpha*e^-4)",
    "(((e^-3*sqrt(phi))*(1-alpha*phi^-2))*(1-alpha*e^-4))*(1-alpha*e^-5)",
    "((((e^-3*sqrt(phi))*(1-alpha*phi^-2))*(1-alpha*e^-4))*(1-alpha*e^-5))*(1+alpha*pi^-6)",
)
_ELECTRIC_ANCHOR = "gamma^4*G^3"
_ELECTRIC_STEPS = (
    "(gamma^4*G^3)*(1+alpha*gamma^3)",
    "((gamma^4*G^3)*(1+alpha*gamma^3))*(1-alpha*pi^-4)",
    "(((gamma^4*G^3)*(1+alpha*gamma^3))*(1-alpha*pi^-4))*(1-alpha*pi^-6)",
    "((((gamma^4*G^3)*(1+alpha*gamma^3))*(1-alpha*pi^-4))*(1-alpha*pi^-6))*(1-alpha*pi^-6)",
)
_ENRICHED_ANCHOR = "gamma^-3-G^4"
_ENRICHED_STEPS = (
    "(gamma^-3-G^4)*(1+alpha*e^-2)",
    "((gamma^-3-G^4)*(1+alpha*e^-2))*(1-alpha*pi^-4)",
)


@dataclass(frozen=True)
class StackStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class StackFigure:
    anchor: NamedGap
    steps: tuple[StackStep, ...]
    target: float


@dataclass(frozen=True)
class Collision:
    parent_name: str
    parent_value: float
    parent_gap: float
    e_name: str
    e_value: float
    e_gap: float
    pi_name: str
    pi_value: float
    pi_gap: float
    required: float
    required_e_gap: float
    required_pi_gap: float
    nearest_piece: str
    nearest_piece_gap: float
    separation: float
    piece_separation: float
    scale: float
    scale_to_low: float
    scale_to_high: float
    target: float
    printed_place: float
    half_step: float
    nested_spacing: float


@dataclass(frozen=True)
class StackReport:
    atomic: StackFigure
    burnup: StackFigure
    stack: StackFigure
    higher: StackFigure
    electric: StackFigure
    enrichment: StackFigure
    fraction_pair: Collision
    enriched_percent: float
    uranium_233_percent: float
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
    return project_root() / "runs" / "production" / "fission_stack.md"


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


def _walk(target, anchor_name, expected, alpha, helium, five, base, positive: bool) -> StackFigure:
    if positive:
        name, value, anchor_gap = _positive_nearest(five, target)
        anchor = NamedGap(name, value, anchor_gap, anchor_gap <= _REVIEW_GAP)
    else:
        anchor = _nearest(five, target)
        name, value = anchor.name, anchor.value
    if anchor.name != anchor_name:
        raise RuntimeError("the stack anchor changed")
    steps: list[StackStep] = []
    for step_name in expected:
        rows = (positive_dressings if positive else unit_dressings)(name, value, alpha, base)
        if positive:
            hit_name, hit_value, hit_gap = _nearest_row(rows, target)
            hit = NamedGap(hit_name, hit_value, hit_gap, hit_gap <= _REVIEW_GAP)
        else:
            hit = _nearest(rows, target)
        if hit.name != step_name:
            raise RuntimeError("the stack dressing changed")
        if hit.relative_gap >= _gap(value, target):
            raise RuntimeError("the stack dressing did not tighten")
        size = (_positive_size if positive else _size)(rows)
        if size.expressions != _COUNT or size.inside != _COUNT:
            raise RuntimeError("the stack dressing count changed")
        if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
            raise RuntimeError("the stack spacing is not the base-piece spacing")
        inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
        if inside_helium > 1:
            raise RuntimeError("several stack dressings meet the helium standard")
        steps.append(
            StackStep(
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
        raise RuntimeError("the closing stack figure misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError("a second stack dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier stack figure already meets helium")
    return StackFigure(anchor, tuple(steps), target)


def _leave(target: float, alpha: float, helium: float, five, base) -> str:
    anchor = _nearest(five, target)
    name, value = anchor.name, anchor.value
    for _step in range(3):
        rows = unit_dressings(name, value, alpha, base)
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
        name, value = hit.name, hit.value
    return "open"


def _fourth_loosens(target: float, alpha: float, helium: float, five, base) -> None:
    anchor = _nearest(five, target)
    name, value = anchor.name, anchor.value
    for _step in range(2):
        hit = _nearest(unit_dressings(name, value, alpha, base), target)
        name, value = hit.name, hit.value
    rows = unit_dressings(name, value, alpha, base)
    hit = _nearest(rows, target)
    inside = sum(_gap(number, target) <= helium for _row, number in rows)
    if inside != 2:
        raise RuntimeError("the double helium hit left two dressings")
    child = _nearest(unit_dressings(hit.name, hit.value, alpha, base), target)
    if child.relative_gap <= hit.relative_gap:
        raise RuntimeError("another factor on the double hit tightened")


def _printed_place(percent: float) -> float:
    place = 1.0
    for _digit in range(8):
        quot = percent / place
        nearest = round(quot)
        if abs(quot - nearest) <= 1.0e-8 and f"{nearest * place:.6g}" == f"{percent:.6g}":
            return place
        place /= 10.0
    raise RuntimeError("the printed place left eight digits")


def _collision(target: float, alpha: float, helium: float, five, base, scale: float, percent: float) -> Collision:
    anchor = _nearest(five, target)
    name = anchor.name
    value = anchor.value
    for _step in range(2):
        rows = unit_dressings(name, value, alpha, base)
        inside = sum(_gap(number, target) <= helium for _row, number in rows)
        if inside != 0:
            raise RuntimeError("an earlier pair step meets helium")
        hit = _nearest(rows, target)
        if hit.relative_gap >= _gap(value, target):
            raise RuntimeError("the pair walk did not tighten")
        name = hit.name
        value = hit.value
    parent_gap = _gap(value, target)
    if parent_gap <= helium:
        raise RuntimeError("the pair parent meets the helium standard")
    rows = unit_dressings(name, value, alpha, base)
    size = _size(rows)
    if size.expressions != _COUNT or size.inside != _COUNT:
        raise RuntimeError("the pair dressing count changed")
    if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
        raise RuntimeError("the pair spacing changed")
    hits = [
        (row, number, _gap(number, target))
        for row, number in rows
        if _gap(number, target) <= helium
    ]
    if len(hits) != 2:
        raise RuntimeError("the pair left two helium dressings")
    sign = -1.0 if value > target else 1.0
    mark = "-" if sign < 0.0 else "+"
    expected = {
        f"({name})*(1{mark}alpha*e^-6)": "e^-6",
        f"({name})*(1{mark}alpha*pi^-5)": "pi^-5",
    }
    if {row for row, _number, _hit_gap in hits} != set(expected):
        raise RuntimeError("the pair left e^-6 and pi^-5")
    found = {expected[row]: (row, number, hit_gap) for row, number, hit_gap in hits}
    e_name, e_value, e_gap = found["e^-6"]
    pi_name, pi_value, pi_gap = found["pi^-5"]
    if not min(e_value, pi_value) < target < max(e_value, pi_value):
        raise RuntimeError("the page figure left the interval between the two dressings")
    required = (target / value - 1.0) / (sign * alpha)
    if required <= 0.0:
        raise RuntimeError("the required factor left the positive pieces")
    pieces = {piece: float(number) for piece, number in base}
    required_e_gap = _gap(required, pieces["e^-6"])
    required_pi_gap = _gap(required, pieces["pi^-5"])
    if required_e_gap <= helium or required_pi_gap <= helium:
        raise RuntimeError("a neighboring piece meets the required factor")
    nearest_piece, nearest_piece_gap = min(
        ((piece, _gap(required, number)) for piece, number in pieces.items() if number > 0.0),
        key=lambda item: (item[1], item[0]),
    )
    if nearest_piece_gap <= helium:
        raise RuntimeError("a base piece meets the required factor")
    separation = abs(e_value - pi_value) / target
    piece_separation = alpha * abs(pieces["e^-6"] - pieces["pi^-5"])
    if _gap(separation, piece_separation * value / target) > 1.0e-9:
        raise RuntimeError("the separation left alpha times the piece difference")
    if not helium < separation < 2.0 * helium:
        raise RuntimeError("the two dressings left the one-window overlap")
    if _gap(scale, target) > helium:
        raise RuntimeError("the named percent left the page figure")
    if scale == e_value or scale == pi_value:
        raise RuntimeError("the named percent selected a neighbor")
    low = min(e_value, pi_value)
    high = max(e_value, pi_value)
    place = _printed_place(percent)
    if _gap(percent / 100.0, target) > 1.0e-12:
        raise RuntimeError("the printed percent left the page figure")
    half_step = (place / 2.0) / percent
    if half_step <= helium:
        raise RuntimeError("the printed half-step met the helium match")
    if e_gap > half_step or pi_gap > half_step or separation > half_step:
        raise RuntimeError("a neighbor left the printed half-step")
    nested = alpha * size.median_spacing
    if not nested < helium:
        raise RuntimeError("the dressed-factor spacing left the fine side of helium")
    return Collision(
        parent_name=name,
        parent_value=value,
        parent_gap=parent_gap,
        e_name=e_name,
        e_value=e_value,
        e_gap=e_gap,
        pi_name=pi_name,
        pi_value=pi_value,
        pi_gap=pi_gap,
        required=required,
        required_e_gap=required_e_gap,
        required_pi_gap=required_pi_gap,
        nearest_piece=nearest_piece,
        nearest_piece_gap=nearest_piece_gap,
        separation=separation,
        piece_separation=piece_separation,
        scale=scale,
        scale_to_low=abs(scale - low) / target,
        scale_to_high=abs(scale - high) / target,
        target=target,
        printed_place=place,
        half_step=half_step,
        nested_spacing=nested,
    )


def _pair_lines(pair: Collision, helium: float, subject: str) -> list[str]:
    operating = subject.removeprefix("page ")
    lines = [
        f"The parent of the {subject} is `{pair.parent_name}` at `{_num(pair.parent_value)}`.",
        f"The parent gap is `{_num(pair.parent_gap)}`.",
        "The parent stays wider than the helium quotient match.",
        f"The `e^-6` neighbor is `{pair.e_name}` at `{_num(pair.e_value)}`.",
        f"The `e^-6` neighbor gap is `{_num(pair.e_gap)}`.",
        f"The `pi^-5` neighbor is `{pair.pi_name}` at `{_num(pair.pi_value)}`.",
        f"The `pi^-5` neighbor gap is `{_num(pair.pi_gap)}`.",
        f"The {subject} sits between those two values.",
        f"Their separation is `{_num(pair.separation)}`.",
        f"alpha times the difference of `e^-6` and `pi^-5` is `{_num(pair.piece_separation)}`.",
        "The separation equals that product times the parent over the page figure.",
        "The separation is wider than one helium quotient match.",
        "The separation is narrower than two helium quotient matches.",
        f"The factor the {subject} requires is `{_num(pair.required)}`.",
        f"Its relative gap to `e^-6` is `{_num(pair.required_e_gap)}`.",
        f"Its relative gap to `pi^-5` is `{_num(pair.required_pi_gap)}`.",
        "Both factor gaps stay wider than the helium quotient match.",
        "No base piece meets the required factor.",
        f"The page prints the percent to the place `{_num(pair.printed_place)}`.",
        f"Half of that place over the printed percent is `{_num(pair.half_step)}`.",
        "Both neighbors sit inside that half-step.",
        "Their separation sits inside that half-step.",
        f"Dressing the factor by another base piece scales the spacing by alpha to `{_num(pair.nested_spacing)}`.",
        "That spacing is finer than the helium quotient match.",
        "The dressed factor stays unused.",
        f"The named percent over 100 sits `{_num(pair.scale_to_low)}` from the lower neighbor.",
        f"It sits `{_num(pair.scale_to_high)}` from the upper neighbor.",
    ]
    if pair.scale_to_low <= helium and pair.scale_to_high <= helium:
        lines.append("The named percent over 100 meets both neighbors inside the helium quotient match.")
    elif pair.scale_to_low <= helium:
        lines.append("The named percent over 100 meets the lower neighbor inside the helium quotient match.")
        lines.append("It stays wider than the helium quotient match from the upper neighbor.")
    elif pair.scale_to_high <= helium:
        lines.append("The named percent over 100 meets the upper neighbor inside the helium quotient match.")
        lines.append("It stays wider than the helium quotient match from the lower neighbor.")
    else:
        lines.append("The named percent over 100 stays wider than the helium quotient match from both neighbors.")
    lines.append("Both neighbors stay uninstalled.")
    lines.append(f"The operating {operating} stays the printed page figure.")
    lines.append("Both routes stay for a later experiment or a later published figure sharper than their separation.")
    return lines


def build_stack() -> StackReport:
    process = build_process()
    module, digest = load_engine()
    if digest[:16] != process.pin_prefix:
        raise RuntimeError("the stack pin does not match the plant pin")
    if process.discharge_gwd != 45.0:
        raise RuntimeError("the discharge is not the 45 GWd material balance")
    if FABRICATION_CAP_PERCENT != 5.0 or not ATOMIC_PERCENT > FABRICATION_CAP_PERCENT:
        raise RuntimeError("the atomic figure is not outside the fabrication cap")
    if _gap(process.burnup_step, HIGHER_BURNUP_GWD / process.discharge_gwd) > 1.0e-12:
        raise RuntimeError("the burnup step left 60 over 45")
    if _gap(process.stack_factor, process.burnup_step * process.carnot_step) > 1.0e-12:
        raise RuntimeError("the stack factor left the burnup times the Carnot step")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    five = _five_pieces(module)
    base = _base_pieces(module)
    if not ENRICHED_PERCENT > 1.0:
        raise RuntimeError("the page enrichment percent stayed inside the unit interval")
    if _leave(ENRICHED_PERCENT / 100.0, alpha, helium, five, base) != "several":
        raise RuntimeError("the 4.5 percent enrichment left its recorded dressing")
    _fourth_loosens(ENRICHED_PERCENT / 100.0, alpha, helium, five, base)
    uranium_233 = next(percent for nucleus, percent in _DELAYED_PERCENT if nucleus == "U-233")
    if _leave(uranium_233 / 100.0, alpha, helium, five, base) != "several":
        raise RuntimeError("the U-233 proportion left its recorded dressing")
    atomic = _walk(ATOMIC_PERCENT / 100.0, _ATOMIC_ANCHOR, _ATOMIC_STEPS, alpha, helium, five, base, False)
    burnup = _walk(process.burnup_step, _BURN_ANCHOR, _BURN_STEPS, alpha, helium, five, base, True)
    stack = _walk(process.stack_factor, _STACK_ANCHOR, _STACK_STEPS, alpha, helium, five, base, True)
    higher = _walk(
        process.higher_over_ranked, _HIGHER_ANCHOR, _HIGHER_STEPS, alpha, helium, five, base, False
    )
    electric = _walk(
        process.stack_over_ranked_electric, _ELECTRIC_ANCHOR, _ELECTRIC_STEPS, alpha, helium, five, base, False
    )
    enrichment = _walk(ENRICHED_PERCENT, _ENRICHED_ANCHOR, _ENRICHED_STEPS, alpha, helium, five, base, True)
    if _gap(enrichment.steps[-1].dressing.value / 100.0, ENRICHED_PERCENT / 100.0) > helium:
        raise RuntimeError("the named enrichment left the plant fraction")
    if not enrichment.steps[-1].dressing.value < FABRICATION_CAP_PERCENT:
        raise RuntimeError("the named enrichment met the fabrication cap")
    if not atomic.steps[-1].dressing.value > FABRICATION_CAP_PERCENT / 100.0:
        raise RuntimeError("the named atomic figure fell to the fabrication cap")
    fraction_pair = _collision(
        ENRICHED_PERCENT / 100.0,
        alpha,
        helium,
        five,
        base,
        enrichment.steps[-1].dressing.value / 100.0,
        ENRICHED_PERCENT,
    )
    if enrichment.steps[-1].dressing.name in (fraction_pair.e_name, fraction_pair.pi_name):
        raise RuntimeError("the named percent selected a neighbor")
    if (
        fraction_pair.e_value * 100.0 >= FABRICATION_CAP_PERCENT
        or fraction_pair.pi_value * 100.0 >= FABRICATION_CAP_PERCENT
    ):
        raise RuntimeError("a neighbor met the fabrication cap")
    return StackReport(
        atomic=atomic,
        burnup=burnup,
        stack=stack,
        higher=higher,
        electric=electric,
        enrichment=enrichment,
        fraction_pair=fraction_pair,
        enriched_percent=ENRICHED_PERCENT,
        uranium_233_percent=uranium_233,
        fabrication_cap_percent=FABRICATION_CAP_PERCENT,
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=process.nuclear_scalar,
        thermodynamics_scalar=process.thermodynamics_scalar,
        gamma_fsot=process.gamma_fsot,
        valve_beta=process.valve_beta,
        pin_prefix=process.pin_prefix,
        digest_prefix=digest[:16],
    )


def _step_lines(step: StackStep, closing: bool) -> list[str]:
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


def _figure_lines(title: str, intro: list[str], figure: StackFigure, nearest: str) -> list[str]:
    lines = ["", f"## {title}", "", *intro]
    lines.append(f"The {nearest} five-atom nearest is `{figure.anchor.name}` at `{_num(figure.anchor.value)}`.")
    lines.append(f"The anchor gap is `{_num(figure.anchor.relative_gap)}`.")
    for step in figure.steps[:-1]:
        lines.extend(_step_lines(step, False))
    lines.extend(_step_lines(figure.steps[-1], True))
    return lines


def render() -> str:
    report = build_stack()
    lines = [
        "# Fission stack",
        "",
        "The overview anchors 60 GWd/t at 6.5 percent of the uranium atoms.",
        "That figure stays outside the fabrication cap.",
        "The 60 GWd/t stack stays unrun.",
        "The chain stays at k_eff of 1.",
    ]
    lines.extend(
        _figure_lines(
            "Atomic anchor",
            [f"The page figure is `{_num(report.atomic.target * 100.0)}` percent."],
            report.atomic,
            "unit-interval",
        )
    )
    lines.extend(
        _figure_lines(
            "Burnup step",
            [f"Sixty over 45 is `{_num(report.burnup.target)}`."],
            report.burnup,
            "positive",
        )
    )
    lines.extend(
        _figure_lines(
            "Stack factor",
            [
                f"The burnup step times the Carnot step is `{_num(report.stack.target)}`.",
                "That product is not a filed plant.",
            ],
            report.stack,
            "positive",
        )
    )
    lines.extend(
        _figure_lines(
            "Ranked share",
            [
                f"The 60 GWd/t discharge over the ranked density is `{_num(report.higher.target)}`.",
                "This share compares those two filed energies.",
            ],
            report.higher,
            "unit-interval",
        )
    )
    lines.extend(
        _figure_lines(
            "Stacked electricity",
            [
                f"The 60 GWd/t discharge at the temperature Carnot, over the ranked electricity, is `{_num(report.electric.target)}`.",
                "That quotient is not a filed plant.",
            ],
            report.electric,
            "unit-interval",
        )
    )
    lines.extend(
        _figure_lines(
            "Page enrichment",
            [
                f"The page percent is `{_num(report.enrichment.target)}`.",
                "The plant fraction is that percent over 100.",
                "The fraction's own dressing meets the helium quotient match more than once.",
                "Another factor on that fraction moves farther from the page figure.",
            ],
            report.enrichment,
            "positive",
        )
    )
    lines.extend(_pair_lines(report.fraction_pair, report.helium_mismatch, "page fraction"))
    lines.append("Both neighbors stay under the fabrication cap.")
    lines.extend(
        [
            f"The U-233 proportion stays `{_num(report.uranium_233_percent)}` percent.",
            "Its proportion dressing meets the helium quotient match more than once.",
            "The page percent is named with the delayed counts.",
            "The fabrication cap stays 5 percent.",
            "This run does not install a seed as a new enrichment or a new plant row.",
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
