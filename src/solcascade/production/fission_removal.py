"""Name the fabrication cap, the published duty, and the cooling-load multiples.

The cap stays the fabrication page's 5 percent handling limit. The duty is
the Surry UFSAR cold-end design point, one residual-heat exchanger. The
one-second multiple closes on three factors. The two-day multiple closes
on a fourth factor of the same form. The one-day multiple keeps both
routes. FSOT names those printed figures. It does not raise the cap, and
it does not assign a watt rate the safety report does not print. The
checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay import decay_fraction, stated_band
from solcascade.production.fission_decay_dress import _base_pieces
from solcascade.production.fission_decay_dress import _dressings as unit_dressings
from solcascade.production.fission_decay_heat import _dressings as positive_dressings
from solcascade.production.fission_decay_heat import (
    _nearest_row,
    _positive_nearest,
    _positive_size,
    build_decay_heat,
)
from solcascade.production.fission_decay_seed import NamedGap, _five_pieces, _nearest, _size
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.fission_process import (
    FABRICATION_CAP_PERCENT,
    FABRICATION_SOURCE,
    FULL_POWER_HOURS,
    NAMEPLATE_MWE,
)
from solcascade.production.fission_safety import SMALL_CLASS_MWE
from solcascade.production.refine import project_root

HALEU_SOURCE = (
    "https://world-nuclear.org/information-library/nuclear-fuel-cycle/"
    "conversion-enrichment-and-fabrication/high-assay-low-enriched-uranium-haleu"
)
RHR_SOURCE = "https://www.nrc.gov/docs/ML2121/ML21210A126.pdf"
CORE_SOURCE = "https://www.nrc.gov/docs/ML2120/ML21208A424.pdf"
# International Table Btu. The UFSAR prints the duty in Btu/hr.
_BTU_IT_J = 1055.05585262
_DUTY_BTU_H = 33.0e6
_CORE_W = 2587.0e6
_SERVICE_S = 4.0 * 3600.0
_COOLDOWN_S = 16.0 * 3600.0
_CAP_ANCHOR = "pi^-3*G^-5"
_CAP_STEPS = (
    "(pi^-3*G^-5)*(1-alpha*phi^-6)",
    "((pi^-3*G^-5)*(1-alpha*phi^-6))*(1-alpha*e^-6)",
)
_DUTY_ANCHOR = "pi^-2*gamma^6"
_DUTY_STEPS = (
    "(pi^-2*gamma^6)*(1-alpha*gamma^2)",
    "((pi^-2*gamma^6)*(1-alpha*gamma^2))*(1+alpha*e^-5)",
    "(((pi^-2*gamma^6)*(1-alpha*gamma^2))*(1+alpha*e^-5))*(1-alpha*pi^-6)",
)
_SECOND_ANCHOR = "gamma^-5-gamma^6"
_SECOND_STEPS = (
    "(gamma^-5-gamma^6)*(1+alpha*pi^-3)",
    "((gamma^-5-gamma^6)*(1+alpha*pi^-3))*(1+alpha*pi^-6)",
    "(((gamma^-5-gamma^6)*(1+alpha*pi^-3))*(1+alpha*pi^-6))*(1+alpha*pi^-6)",
)
_TWO_DAY_ANCHOR = "G^3+G^5"
_TWO_DAY_STEPS = (
    "(G^3+G^5)*(1+alpha*phi^-5)",
    "((G^3+G^5)*(1+alpha*phi^-5))*(1-alpha*e^-5)",
    "(((G^3+G^5)*(1+alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6)",
    "((((G^3+G^5)*(1+alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6))*(1+alpha*pi^-6)",
)
_DAY_ANCHOR = "sqrt(pi)-e^-6"
_DAY_PARENT = "(sqrt(pi)-e^-6)*(1-alpha*e^-4)"
_DAY_E_CHILD = "(((sqrt(pi)-e^-6)*(1-alpha*e^-4))*(1-alpha*e^-6))*(1-alpha*pi^-6)"
_DAY_PI_CHILD = "(((sqrt(pi)-e^-6)*(1-alpha*e^-4))*(1-alpha*pi^-5))*(1+alpha*pi^-6)"
_COUNT = 116
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4


@dataclass(frozen=True)
class RemovalStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class DutyPair:
    anchor: NamedGap
    parent: RemovalStep
    e_name: str
    e_value: float
    e_gap: float
    pi_name: str
    pi_value: float
    pi_gap: float
    e_child_name: str
    e_child_gap: float
    pi_child_name: str
    pi_child_gap: float
    required: float
    required_e_gap: float
    required_pi_gap: float
    nearest_piece: str
    nearest_piece_gap: float
    separation: float
    piece_separation: float
    expressions: int
    median_spacing: float
    nested_spacing: float


@dataclass(frozen=True)
class RemovalReport:
    cap_anchor: NamedGap
    cap_steps: tuple[RemovalStep, ...]
    cap_value: float
    duty_anchor: NamedGap
    duty_steps: tuple[RemovalStep, ...]
    duty_w: float
    core_w: float
    duty_fraction: float
    derived_fraction: float
    surry_duty_w: float
    plant_duty_w: float
    small_duty_w: float
    service_fraction: float
    cooldown_fraction: float
    service_band: float
    cross_hours: float
    one_second_over_duty: float
    one_day_over_duty: float
    two_day_over_duty: float
    one_second_anchor: NamedGap
    one_second_steps: tuple[RemovalStep, ...]
    two_day_anchor: NamedGap
    two_day_steps: tuple[RemovalStep, ...]
    one_day_pair: DutyPair
    condenser_w: float
    one_second_w: float
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_removal.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _gap(derived: float, witness: float) -> float:
    if witness == 0.0:
        raise RuntimeError("the witness is zero")
    return abs(derived - witness) / abs(witness)


def _second_gap(rows, target: float, best: float) -> float:
    gaps = sorted(_gap(value, target) for _name, value in rows)
    if len(gaps) < 2 or gaps[0] != best and abs(gaps[0] - best) > 1.0e-15:
        raise RuntimeError("the nearest dressing is missing from its own list")
    return gaps[1]


def _walk(target: float, anchor: NamedGap, expected: tuple[str, ...], alpha: float, helium: float, base) -> tuple[RemovalStep, ...]:
    name = anchor.name
    value = anchor.value
    steps: list[RemovalStep] = []
    for step_name in expected:
        rows = unit_dressings(name, value, alpha, base)
        hit = _nearest(rows, target)
        if hit.name != step_name:
            raise RuntimeError("the removal dressing changed")
        if hit.relative_gap >= _gap(value, target):
            raise RuntimeError("the removal dressing did not tighten")
        size = _size(rows)
        if size.expressions != _COUNT or size.inside != _COUNT:
            raise RuntimeError("the removal dressing count changed")
        if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
            raise RuntimeError("the removal spacing is not the base-piece spacing")
        inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
        if inside_helium > 1:
            raise RuntimeError("several removal dressings meet the helium standard")
        steps.append(
            RemovalStep(
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
        raise RuntimeError("the closing removal figure misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError("a second removal dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier removal figure already meets helium")
    return tuple(steps)


def _positive_record(rows, target: float, hit: NamedGap, helium: float) -> RemovalStep:
    size = _positive_size(rows)
    if size.expressions != _COUNT or size.inside != _COUNT:
        raise RuntimeError("the removal dressing count changed")
    if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
        raise RuntimeError("the removal spacing is not the base-piece spacing")
    inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
    if inside_helium > 1:
        raise RuntimeError("several removal dressings meet the helium standard")
    return RemovalStep(
        dressing=hit,
        inside_review=sum(_gap(number, target) <= _REVIEW_GAP for _row, number in rows),
        inside_tenth=sum(_gap(number, target) <= _TENTH for _row, number in rows),
        inside_helium=inside_helium,
        second_gap=_second_gap(rows, target, hit.relative_gap),
        expressions=size.expressions,
        median_spacing=size.median_spacing,
    )


def _walk_positive(target: float, anchor_name: str, expected: tuple[str, ...], alpha: float, helium: float, five, base) -> tuple[NamedGap, tuple[RemovalStep, ...]]:
    if not target > 1.0:
        raise RuntimeError("the positive removal walk left the ratios above one")
    name, value, anchor_gap = _positive_nearest(five, target)
    if name != anchor_name:
        raise RuntimeError(f"the removal anchor changed to {name}")
    anchor = NamedGap(name, value, anchor_gap, anchor_gap <= _REVIEW_GAP)
    steps: list[RemovalStep] = []
    for step_name in expected:
        rows = positive_dressings(name, value, alpha, base)
        hit_name, hit_value, hit_gap = _nearest_row(rows, target)
        if hit_name != step_name:
            raise RuntimeError(f"the removal dressing changed to {hit_name}")
        if hit_gap >= _gap(value, target):
            raise RuntimeError("the removal dressing did not tighten")
        hit = NamedGap(hit_name, hit_value, hit_gap, hit_gap <= _REVIEW_GAP)
        steps.append(_positive_record(rows, target, hit, helium))
        name = hit_name
        value = hit_value
    closing = steps[-1]
    if closing.inside_helium != 1 or closing.dressing.relative_gap > helium:
        raise RuntimeError("the closing removal figure misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError("a second removal dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier removal figure already meets helium")
    return anchor, tuple(steps)


def _duty_child(parent_name: str, parent_value: float, parent_gap: float, expected: str, target: float, alpha: float, helium: float, base) -> tuple[str, float]:
    rows = positive_dressings(parent_name, parent_value, alpha, base)
    size = _positive_size(rows)
    if size.expressions != _COUNT or size.inside != _COUNT:
        raise RuntimeError("the one-day child count changed")
    if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
        raise RuntimeError("the one-day child spacing changed")
    name, _value, gap = _nearest_row(rows, target)
    if name != expected:
        raise RuntimeError(f"the one-day child changed to {name}")
    if gap <= parent_gap or gap <= helium:
        raise RuntimeError("another factor on the one-day pair tightened")
    return name, gap


def _duty_pair(target: float, anchor_name: str, parent_name: str, alpha: float, helium: float, five, base) -> DutyPair:
    if not target > 1.0:
        raise RuntimeError("the one-day multiple left the ratios above one")
    name, value, anchor_gap = _positive_nearest(five, target)
    if name != anchor_name:
        raise RuntimeError(f"the one-day anchor changed to {name}")
    anchor = NamedGap(name, value, anchor_gap, anchor_gap <= _REVIEW_GAP)
    rows = positive_dressings(name, value, alpha, base)
    hit_name, hit_value, hit_gap = _nearest_row(rows, target)
    if hit_name != parent_name:
        raise RuntimeError(f"the one-day parent changed to {hit_name}")
    if hit_gap >= anchor_gap:
        raise RuntimeError("the one-day parent did not tighten")
    parent = _positive_record(
        rows, target, NamedGap(hit_name, hit_value, hit_gap, hit_gap <= _REVIEW_GAP), helium
    )
    if parent.inside_helium != 0 or parent.dressing.relative_gap <= helium:
        raise RuntimeError("the one-day parent meets the helium standard")
    rows = positive_dressings(hit_name, hit_value, alpha, base)
    size = _positive_size(rows)
    if size.expressions != _COUNT or size.inside != _COUNT:
        raise RuntimeError("the one-day dressing count changed")
    if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
        raise RuntimeError("the one-day spacing changed")
    hits = [
        (row, number, _gap(number, target))
        for row, number in rows
        if _gap(number, target) <= helium
    ]
    if len(hits) != 2:
        raise RuntimeError("the one-day pair left two helium dressings")
    sign = -1.0 if hit_value > target else 1.0
    mark = "-" if sign < 0.0 else "+"
    expected = {
        f"({hit_name})*(1{mark}alpha*e^-6)": "e^-6",
        f"({hit_name})*(1{mark}alpha*pi^-5)": "pi^-5",
    }
    if {row for row, _number, _hit_gap in hits} != set(expected):
        raise RuntimeError("the one-day pair left e^-6 and pi^-5")
    found = {expected[row]: (row, number, hit_gap) for row, number, hit_gap in hits}
    e_name, e_value, e_gap = found["e^-6"]
    pi_name, pi_value, pi_gap = found["pi^-5"]
    if not pi_value < target < e_value:
        raise RuntimeError("the one-day quotient left the interval between the two dressings")
    if _num(pi_value) != _num(target) or pi_value == target:
        raise RuntimeError("the pi neighbor left the six-digit print of the quotient")
    if _num(e_value) == _num(target):
        raise RuntimeError("the e neighbor printed the quotient")
    required = (target / hit_value - 1.0) / (sign * alpha)
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
    if nearest_piece != "pi^-5" or nearest_piece_gap <= helium:
        raise RuntimeError("the one-day required factor left pi^-5")
    separation = abs(e_value - pi_value) / target
    piece_separation = alpha * abs(pieces["e^-6"] - pieces["pi^-5"])
    if _gap(separation, piece_separation * hit_value / target) > 1.0e-9:
        raise RuntimeError("the separation left alpha times the piece difference")
    if not helium < separation < 2.0 * helium:
        raise RuntimeError("the two dressings left the one-window overlap")
    nested = alpha * size.median_spacing
    if not nested < helium:
        raise RuntimeError("the dressed-factor spacing left the fine side of helium")
    e_child_name, e_child_gap = _duty_child(
        e_name, e_value, e_gap, _DAY_E_CHILD, target, alpha, helium, base
    )
    pi_child_name, pi_child_gap = _duty_child(
        pi_name, pi_value, pi_gap, _DAY_PI_CHILD, target, alpha, helium, base
    )
    return DutyPair(
        anchor=anchor,
        parent=parent,
        e_name=e_name,
        e_value=e_value,
        e_gap=e_gap,
        pi_name=pi_name,
        pi_value=pi_value,
        pi_gap=pi_gap,
        e_child_name=e_child_name,
        e_child_gap=e_child_gap,
        pi_child_name=pi_child_name,
        pi_child_gap=pi_child_gap,
        required=required,
        required_e_gap=required_e_gap,
        required_pi_gap=required_pi_gap,
        nearest_piece=nearest_piece,
        nearest_piece_gap=nearest_piece_gap,
        separation=separation,
        piece_separation=piece_separation,
        expressions=size.expressions,
        median_spacing=size.median_spacing,
        nested_spacing=nested,
    )


def _cross_hours(target: float, operation_s: float) -> float:
    lo, hi = 1.0, 30.0 * 86400.0
    if not decay_fraction(lo, operation_s) > target > decay_fraction(hi, operation_s):
        raise RuntimeError("the decay fraction does not cross the design-point duty")
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if decay_fraction(mid, operation_s) > target:
            lo = mid
        else:
            hi = mid
    return lo / 3600.0


def build_removal() -> RemovalReport:
    if FABRICATION_CAP_PERCENT != 5.0:
        raise RuntimeError("the fabrication cap is not 5 percent")
    heat = build_decay_heat()
    module, digest = load_engine()
    if digest[:16] != heat.pin_prefix or heat.digest_prefix != heat.pin_prefix:
        raise RuntimeError("the removal pin does not match the plant pin")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    if abs(helium - heat.helium_mismatch) > 0.0:
        raise RuntimeError("the helium standard changed")
    five = _five_pieces(module)
    base = _base_pieces(module)
    cap_target = FABRICATION_CAP_PERCENT / 100.0
    cap_anchor = _nearest(five, cap_target)
    if cap_anchor.name != _CAP_ANCHOR:
        raise RuntimeError("the fabrication-cap anchor changed")
    duty_w = _DUTY_BTU_H * _BTU_IT_J / 3600.0
    duty_fraction = duty_w / _CORE_W
    duty_anchor = _nearest(five, duty_fraction)
    if duty_anchor.name != _DUTY_ANCHOR:
        raise RuntimeError("the removal-duty anchor changed")
    cap_steps = _walk(cap_target, cap_anchor, _CAP_STEPS, alpha, helium, base)
    duty_steps = _walk(duty_fraction, duty_anchor, _DUTY_STEPS, alpha, helium, base)
    derived_fraction = duty_steps[-1].dressing.value
    surry_duty_w = derived_fraction * _CORE_W
    plant_duty_w = derived_fraction * heat.thermal_w
    small_duty_w = derived_fraction * heat.small_thermal_w
    if _gap(surry_duty_w, duty_w) > helium or _gap(derived_fraction, duty_fraction) > helium:
        raise RuntimeError("the derived duty misses the published design point")
    operation_s = FULL_POWER_HOURS * 3600.0
    service_fraction = decay_fraction(_SERVICE_S, operation_s)
    cooldown_fraction = decay_fraction(_SERVICE_S + _COOLDOWN_S, operation_s)
    service_band = stated_band(_SERVICE_S)
    cooldown_band = stated_band(_SERVICE_S + _COOLDOWN_S)
    if service_band != 0.10 or cooldown_band != 0.10:
        raise RuntimeError("the UFSAR times are outside the 10 percent band")
    if cooldown_fraction * (1.0 - cooldown_band) <= duty_fraction:
        raise RuntimeError("the cooldown band reaches the design-point duty")
    loads = heat.loads
    if tuple(load.label for load in loads) != (
        "one-second fraction",
        "one-day fraction",
        "two-day fraction",
    ):
        raise RuntimeError("the cooling-load order changed")
    one_second, one_day, two_day = loads
    if not one_second.watts > one_day.watts > two_day.watts > plant_duty_w:
        raise RuntimeError("a scored cooling load sits inside the design-point duty")
    condenser_w = heat.thermal_w - NAMEPLATE_MWE * 1.0e6
    if not 0.0 < one_second.watts < condenser_w:
        raise RuntimeError("the one-second load is outside the running condenser")
    cross_hours = _cross_hours(duty_fraction, operation_s)
    if not cross_hours > 48.0:
        raise RuntimeError("the design-point duty covers the two-day load")
    one_second_over_duty = one_second.watts / plant_duty_w
    one_day_over_duty = one_day.watts / plant_duty_w
    two_day_over_duty = two_day.watts / plant_duty_w
    one_second_anchor, one_second_steps = _walk_positive(
        one_second_over_duty, _SECOND_ANCHOR, _SECOND_STEPS, alpha, helium, five, base
    )
    two_day_anchor, two_day_steps = _walk_positive(
        two_day_over_duty, _TWO_DAY_ANCHOR, _TWO_DAY_STEPS, alpha, helium, five, base
    )
    one_day_pair = _duty_pair(one_day_over_duty, _DAY_ANCHOR, _DAY_PARENT, alpha, helium, five, base)
    if len(one_second_steps) != 3:
        raise RuntimeError("the one-second multiple left three factors")
    if len(two_day_steps) != 4 or two_day_steps[-2].inside_helium != 0:
        raise RuntimeError("the two-day multiple met helium before the fourth factor")
    return RemovalReport(
        cap_anchor=cap_anchor,
        cap_steps=cap_steps,
        cap_value=cap_steps[-1].dressing.value,
        duty_anchor=duty_anchor,
        duty_steps=duty_steps,
        duty_w=duty_w,
        core_w=_CORE_W,
        duty_fraction=duty_fraction,
        derived_fraction=derived_fraction,
        surry_duty_w=surry_duty_w,
        plant_duty_w=plant_duty_w,
        small_duty_w=small_duty_w,
        service_fraction=service_fraction,
        cooldown_fraction=cooldown_fraction,
        service_band=service_band,
        cross_hours=cross_hours,
        one_second_over_duty=one_second_over_duty,
        one_day_over_duty=one_day_over_duty,
        two_day_over_duty=two_day_over_duty,
        one_second_anchor=one_second_anchor,
        one_second_steps=one_second_steps,
        two_day_anchor=two_day_anchor,
        two_day_steps=two_day_steps,
        one_day_pair=one_day_pair,
        condenser_w=condenser_w,
        one_second_w=one_second.watts,
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=heat.nuclear_scalar,
        thermodynamics_scalar=heat.thermodynamics_scalar,
        gamma_fsot=heat.gamma_fsot,
        valve_beta=heat.valve_beta,
        pin_prefix=heat.pin_prefix,
        digest_prefix=heat.digest_prefix,
    )


def _step_lines(step: RemovalStep, closing: bool) -> list[str]:
    lines = [
        f"The dressing is `{step.dressing.name}` at `{_num(step.dressing.value)}`.",
        f"The relative gap is `{_num(step.dressing.relative_gap)}`.",
    ]
    if closing:
        lines.append("The relative gap meets the helium quotient match.")
        lines.append("One dressing of that anchor meets the helium quotient match.")
        lines.append(
            f"The next dressing sits at relative gap `{_num(step.second_gap)}`."
        )
    else:
        lines.append("The dressing stays wider than the helium quotient match.")
    lines.append(
        f"The anchor has `{step.expressions}` dressings, median spacing `{_num(step.median_spacing)}`."
    )
    return lines


def _multiple_lines(anchor: NamedGap, steps: tuple[RemovalStep, ...], notes: tuple[str, ...], quotient_note: str) -> list[str]:
    lines = [
        f"The five-atom nearest is `{anchor.name}` at `{_num(anchor.value)}`.",
        f"The anchor gap is `{_num(anchor.relative_gap)}`.",
    ]
    for step in steps[:-1]:
        lines.extend(_step_lines(step, False))
    lines.extend(notes)
    lines.extend(_step_lines(steps[-1], True))
    lines.append("The closing value names that multiple.")
    lines.append(quotient_note)
    return lines


def _day_lines(pair: DutyPair) -> list[str]:
    lines = [
        f"The five-atom nearest is `{pair.anchor.name}` at `{_num(pair.anchor.value)}`.",
        f"The anchor gap is `{_num(pair.anchor.relative_gap)}`.",
    ]
    lines.extend(_step_lines(pair.parent, False))
    lines.extend(
        [
            "That dressing is the parent of the one-day multiple.",
            f"The `e^-6` neighbor is `{pair.e_name}` at `{_num(pair.e_value)}`.",
            f"The `e^-6` neighbor gap is `{_num(pair.e_gap)}`.",
            f"The `pi^-5` neighbor is `{pair.pi_name}` at `{_num(pair.pi_value)}`.",
            f"The `pi^-5` neighbor gap is `{_num(pair.pi_gap)}`.",
            "The one-day multiple sits between those two values.",
            "Two dressings of that parent meet the helium quotient match.",
            f"Their separation is `{_num(pair.separation)}`.",
            f"alpha times the difference of `e^-6` and `pi^-5` is `{_num(pair.piece_separation)}`.",
            "The separation equals that product times the parent over the quotient.",
            "The separation is wider than one helium quotient match.",
            "The separation is narrower than two helium quotient matches.",
            f"The factor the one-day multiple requires is `{_num(pair.required)}`.",
            f"Its relative gap to `e^-6` is `{_num(pair.required_e_gap)}`.",
            f"Its relative gap to `pi^-5` is `{_num(pair.required_pi_gap)}`.",
            "Both factor gaps stay wider than the helium quotient match.",
            f"The nearest base piece is `{pair.nearest_piece}`.",
            "No base piece meets the required factor.",
            f"Another factor on the `e^-6` neighbor is `{pair.e_child_name}`.",
            f"Its relative gap is `{_num(pair.e_child_gap)}`.",
            f"Another factor on the `pi^-5` neighbor is `{pair.pi_child_name}`.",
            f"Its relative gap is `{_num(pair.pi_child_gap)}`.",
            "Another factor on either neighbor moves farther from the quotient.",
            f"Dressing the factor by another base piece scales the spacing by alpha to `{_num(pair.nested_spacing)}`.",
            "That spacing is finer than the helium quotient match.",
            "The dressed factor stays unused.",
            "The `pi^-5` neighbor prints the same six digits as the quotient.",
            "The two values stay distinct.",
            "Both neighbors stay uninstalled.",
            "The operating multiple stays the quotient of the named one-day load and the named duty.",
            "Both routes stay for a later experiment or a later published figure sharper than their separation.",
            "This score makes no prediction that chooses between them.",
        ]
    )
    return lines


def render() -> str:
    report = build_removal()
    cap_added, cap_closing = report.cap_steps
    duty_added, duty_wider, duty_closing = report.duty_steps
    lines = [
        "# Fission removal",
        "",
        "The fabrication cap and the second heat-removal duty are the two bounds the stress run left open.",
        "Both figures are read from public safety and fuel-cycle reports.",
        "FSOT names the printed figures.",
        "The comparison names the one-second multiple and the two-day multiple of the design-point duty.",
        "The one-day multiple keeps both routes.",
        "The chain stays at k_eff of 1.",
        "",
        "## Cap",
        "",
        f"The fabrication page is {FABRICATION_SOURCE}.",
        "That page keeps a normal handling limit of 5 percent U-235 for conventional light-water fuel.",
        "Light-water fuel on the page is enriched up to about 4.8 percent.",
        "Burn-up on the page runs from about 40 GWd/t to more than 60 GWd/t as enrichment runs toward that 5 percent limit.",
        f"The fabrication cap remains `{_num(FABRICATION_CAP_PERCENT)}` percent U-235.",
        f"The five-atom nearest is `{report.cap_anchor.name}` at `{_num(report.cap_anchor.value)}`.",
        f"The anchor gap is `{_num(report.cap_anchor.relative_gap)}`.",
    ]
    lines.extend(_step_lines(cap_added, False))
    lines.extend(_step_lines(cap_closing, True))
    lines.extend(
        [
            f"The closing value agrees with `{_num(FABRICATION_CAP_PERCENT / 100.0)}` at six figures.",
            "The 60 GWd/t stack stays unrun.",
            f"The HALEU page is {HALEU_SOURCE}.",
            "That page defines HALEU as enrichment above 5 percent and below 20 percent.",
            "It records test assemblies up to 6 percent U-235 at Vogtle 2.",
            "Those assemblies sit outside this plant's fabrication cap.",
            "This score does not install them.",
            "",
            "## Removal",
            "",
            f"The residual-heat section is {RHR_SOURCE}.",
            "Surry UFSAR Revision 52 prints two residual-heat exchangers.",
            "The design heat transfer rate of each is 33 times 10 to the sixth Btu per hour.",
            "The table's tube-side temperatures are 140°F in and 124°F out.",
            "That pair is the cold end of the cooldown.",
            "One pump and one exchanger are the train the section says is enough for the decay-heat function.",
            "The section places that train in service about four hours after shutdown, near 350°F.",
            "The coolant then falls from 350°F to 140°F in 16 hours.",
            f"The core rating in the same revision is {CORE_SOURCE}.",
            f"The reactor core average thermal power there is `{_num(report.core_w)}` W.",
            "The watt figure uses the International Table Btu, 1055.05585262 J.",
            f"One exchanger at the design point is `{_num(report.duty_w)}` W.",
            f"That duty is `{_num(report.duty_fraction)}` of the core rating.",
            f"The five-atom nearest is `{report.duty_anchor.name}` at `{_num(report.duty_anchor.value)}`.",
            f"The anchor gap is `{_num(report.duty_anchor.relative_gap)}`.",
        ]
    )
    lines.extend(_step_lines(duty_added, False))
    lines.extend(_step_lines(duty_wider, False))
    lines.extend(_step_lines(duty_closing, True))
    lines.extend(
        [
            f"The closing fraction is `{_num(report.derived_fraction)}`.",
            f"The derived Surry duty is `{_num(report.surry_duty_w)}` W.",
            f"The same fraction times the FSOT operating heat is `{_num(report.plant_duty_w)}` W.",
            f"The small-class design-point duty is `{_num(report.small_duty_w)}` W.",
            "This is the cold-end design point scaled by core power.",
            "The UFSAR does not print a watt rate at the 350°F entry.",
            "This score assigns that hotter rate no figure.",
            "",
            "## Comparison",
            "",
            f"At four hours the 1952 fraction is `{_num(report.service_fraction)}`.",
            f"At twenty hours, the end of the published cooldown, it is `{_num(report.cooldown_fraction)}`.",
            "Both times sit in the report's 10 percent band.",
            "The low edge of that band at twenty hours stays above the design-point fraction.",
            f"On this plant the one-second load is `{_num(report.one_second_over_duty)}` times the scaled design-point duty.",
        ]
    )
    lines.extend(
        _multiple_lines(
            report.one_second_anchor,
            report.one_second_steps,
            ("The third factor closes the one-second multiple.",),
            "The operating multiple stays the quotient of the named one-second load and the named duty.",
        )
    )
    lines.extend(
        [
            f"The one-day load is `{_num(report.one_day_over_duty)}` times that duty.",
        ]
    )
    lines.extend(_day_lines(report.one_day_pair))
    lines.extend(
        [
            f"The two-day load is `{_num(report.two_day_over_duty)}` times that duty.",
        ]
    )
    lines.extend(
        _multiple_lines(
            report.two_day_anchor,
            report.two_day_steps,
            (
                "Three factors of the two-day multiple stay wider than the helium quotient match.",
                "The fourth factor is one more base-piece factor of the same form.",
            ),
            "The operating multiple stays the quotient of the named two-day load and the named duty.",
        )
    )
    lines.extend(
        [
            f"The one-year fraction falls through the design-point fraction at `{_num(report.cross_hours)}` hours.",
            f"The running condenser duty is `{_num(report.condenser_w)}` W.",
            f"The one-second load, `{_num(report.one_second_w)}` W, sits inside that running condenser.",
            "The condenser remains the early removal path.",
            "A loss of condenser circulation has no second watt capacity on file.",
            "This score assigns that loss no figure.",
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
            f"alpha in the dressing is the fine-structure leaf `{_num(report.alpha)}`.",
            f"The helium quotient match is `{_num(report.helium_mismatch)}`.",
            "Every closing gap meets the helium quotient match.",
            "The page efficiencies and the discharge ratio stay on their earlier scores.",
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
