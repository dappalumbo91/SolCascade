"""Name the filed quotients the ceiling, the year, and the residual-heat table left open.

The round-file margin is the catalog ceiling minus the temperature Carnot.
The rejected fraction keeps both routes. The separative rate, the year
shares, the offset ratio, and the two temperature ratios close on one
dressing each. The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_ceiling import build_ceiling
from solcascade.production.fission_decay import _EXPONENT, _LONG_OFFSET_S, _SHORT_OFFSET_S, build_decay
from solcascade.production.fission_decay_dress import _base_pieces
from solcascade.production.fission_decay_dress import _dressings as unit_dressings
from solcascade.production.fission_decay_heat import _dressings as positive_dressings
from solcascade.production.fission_decay_heat import _nearest_row, _positive_nearest, _positive_size
from solcascade.production.fission_decay_seed import NamedGap, _five_pieces, _nearest, _size
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.fission_process import build_process
from solcascade.production.plants import FAST_ID
from solcascade.production.refine import catalog_record, project_root

_COUNT = 116
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4
_MARGIN_ANCHOR = "phi^-6*gamma^5"
_MARGIN_STEPS = (
    "(phi^-6*gamma^5)*(1+alpha*e^-4)",
    "((phi^-6*gamma^5)*(1+alpha*e^-4))*(1+alpha*e^-5)",
    "(((phi^-6*gamma^5)*(1+alpha*e^-4))*(1+alpha*e^-5))*(1-alpha*pi^-6)",
)
_RATE_ANCHOR = "e^4*G"
_RATE_STEPS = (
    "(e^4*G)*(1-alpha*pi^-3)",
    "((e^4*G)*(1-alpha*pi^-3))*(1+alpha*pi^-5)",
    "(((e^4*G)*(1-alpha*pi^-3))*(1+alpha*pi^-5))*(1+alpha*pi^-6)",
)
_DAY_ANCHOR = "pi^-5*G^2"
_DAY_STEPS = (
    "(pi^-5*G^2)*(1-alpha*phi^-5)",
    "((pi^-5*G^2)*(1-alpha*phi^-5))*(1-alpha*e^-5)",
    "(((pi^-5*G^2)*(1-alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6)",
    "((((pi^-5*G^2)*(1-alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6))*(1+alpha*pi^-6)",
)
_COOL_ANCHOR = "e^-6*G"
_COOL_STEPS = (
    "(e^-6*G)*(1+alpha*sqrt(gamma))",
    "((e^-6*G)*(1+alpha*sqrt(gamma)))*(1+alpha*pi^-5)",
    "(((e^-6*G)*(1+alpha*sqrt(gamma)))*(1+alpha*pi^-5))*(1+alpha*pi^-6)",
)
_OFFSET_ANCHOR = "pi^6*pi^6"
_OFFSET_STEPS = (
    "(pi^6*pi^6)*(1+alpha*e^5)",
    "((pi^6*pi^6)*(1+alpha*e^5))*(1+alpha*gamma^-3)",
    "(((pi^6*pi^6)*(1+alpha*e^5))*(1+alpha*gamma^-3))*(1+alpha*gamma^4)",
    "((((pi^6*pi^6)*(1+alpha*e^5))*(1+alpha*gamma^-3))*(1+alpha*gamma^4))*(1+alpha*pi^-5)",
)
_TUBE_ANCHOR = "phi-G^5"
_TUBE_STEPS = (
    "(phi-G^5)*(1+alpha*e^-5)",
    "((phi-G^5)*(1+alpha*e^-5))*(1-alpha*pi^-6)",
)
_ENDS_ANCHOR = "gamma^6+G^4"
_ENDS_STEPS = (
    "(gamma^6+G^4)*(1-alpha*e^-3)",
    "((gamma^6+G^4)*(1-alpha*e^-3))*(1+alpha*e^-6)",
    "(((gamma^6+G^4)*(1-alpha*e^-3))*(1+alpha*e^-6))*(1-alpha*pi^-6)",
)
_REJECT_ANCHOR = "gamma*sqrt(G)"
_REJECT_FIRST = "(gamma*sqrt(G))*(1+alpha*phi^-3)"
_REJECT_PARENT = "((gamma*sqrt(G))*(1+alpha*phi^-3))*(1+alpha*e^-3)"
_REJECT_E_CHILD = "((((gamma*sqrt(G))*(1+alpha*phi^-3))*(1+alpha*e^-3))*(1-alpha*e^-6))*(1-alpha*pi^-6)"
_REJECT_PI_CHILD = "((((gamma*sqrt(G))*(1+alpha*phi^-3))*(1+alpha*e^-3))*(1-alpha*pi^-5))*(1+alpha*pi^-6)"
_TUBE_IN_F = 140.0
_TUBE_OUT_F = 124.0
_COOLDOWN_HOT_F = 350.0
_COOLDOWN_COLD_F = 140.0


@dataclass(frozen=True)
class FiledStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class FiledFigure:
    anchor: NamedGap
    steps: tuple[FiledStep, ...]
    target: float


@dataclass(frozen=True)
class RejectedPair:
    anchor: NamedGap
    first: FiledStep
    parent: FiledStep
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
    nested_spacing: float
    target: float


@dataclass(frozen=True)
class FiledReport:
    margin: FiledFigure
    rejected: RejectedPair
    rate: FiledFigure
    day_share: FiledFigure
    cooldown_share: FiledFigure
    offset_ratio: FiledFigure
    tube_ratio: FiledFigure
    cooldown_ratio: FiledFigure
    hot_k: float
    cold_k: float
    margin_j: float
    rate_kwh: float
    short_offset_s: float
    long_offset_s: float
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_filed.md"


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


def _kelvin(fahrenheit: float) -> float:
    return (fahrenheit - 32.0) * 5.0 / 9.0 + 273.15


def _record(rows, target: float, hit: NamedGap, helium: float, positive: bool) -> FiledStep:
    size = (_positive_size if positive else _size)(rows)
    if size.expressions != _COUNT or size.inside != _COUNT:
        raise RuntimeError("the filed dressing count changed")
    if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
        raise RuntimeError("the filed spacing is not the base-piece spacing")
    inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
    if inside_helium > 1:
        raise RuntimeError("several filed dressings meet the helium standard")
    return FiledStep(
        dressing=hit,
        inside_review=sum(_gap(number, target) <= _REVIEW_GAP for _row, number in rows),
        inside_tenth=sum(_gap(number, target) <= _TENTH for _row, number in rows),
        inside_helium=inside_helium,
        second_gap=_second_gap(rows, target, hit.relative_gap),
        expressions=size.expressions,
        median_spacing=size.median_spacing,
    )


def _walk(target, anchor_name, expected, alpha, helium, five, base, positive: bool) -> FiledFigure:
    if positive:
        name, value, anchor_gap = _positive_nearest(five, target)
        anchor = NamedGap(name, value, anchor_gap, anchor_gap <= _REVIEW_GAP)
    else:
        anchor = _nearest(five, target)
        name, value = anchor.name, anchor.value
    if anchor.name != anchor_name:
        raise RuntimeError(f"the filed anchor changed to {anchor.name}")
    if anchor.relative_gap <= helium:
        raise RuntimeError("the filed anchor already meets helium")
    if positive and not target > 1.0:
        dropped = len(unit_dressings(anchor.name, anchor.value, alpha, base))
        if dropped >= _COUNT:
            raise RuntimeError("the unit grammar still has 116 dressings")
    steps: list[FiledStep] = []
    for step_name in expected:
        rows = (positive_dressings if positive else unit_dressings)(name, value, alpha, base)
        if positive:
            hit_name, hit_value, hit_gap = _nearest_row(rows, target)
            hit = NamedGap(hit_name, hit_value, hit_gap, hit_gap <= _REVIEW_GAP)
        else:
            hit = _nearest(rows, target)
        if hit.name != step_name:
            raise RuntimeError(f"the filed dressing changed to {hit.name}")
        if hit.relative_gap >= _gap(value, target):
            raise RuntimeError("the filed dressing did not tighten")
        recorded = _record(rows, target, hit, helium, positive)
        if recorded.inside_helium > 1:
            raise RuntimeError("several filed dressings meet the helium standard")
        steps.append(recorded)
        name, value = hit.name, hit.value
    closing = steps[-1]
    if closing.inside_helium != 1 or closing.dressing.relative_gap > helium:
        raise RuntimeError("the closing filed figure misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError("a second filed dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier filed figure already meets helium")
    return FiledFigure(anchor, tuple(steps), target)


def _rejected(target: float, alpha: float, helium: float, five, base) -> RejectedPair:
    if not 0.0 < target < 1.0:
        raise RuntimeError("the rejected fraction left the unit interval")
    name, value, anchor_gap = _positive_nearest(five, target)
    if name != _REJECT_ANCHOR:
        raise RuntimeError(f"the rejected anchor changed to {name}")
    if len(unit_dressings(name, value, alpha, base)) >= _COUNT:
        raise RuntimeError("the unit grammar still has 116 dressings")
    anchor = NamedGap(name, value, anchor_gap, anchor_gap <= _REVIEW_GAP)
    if anchor.relative_gap <= helium:
        raise RuntimeError("the rejected anchor already meets helium")
    rows = positive_dressings(name, value, alpha, base)
    hit_name, hit_value, hit_gap = _nearest_row(rows, target)
    if hit_name != _REJECT_FIRST:
        raise RuntimeError(f"the rejected first factor changed to {hit_name}")
    if hit_gap >= anchor_gap:
        raise RuntimeError("the rejected first factor did not tighten")
    first = _record(rows, target, NamedGap(hit_name, hit_value, hit_gap, hit_gap <= _REVIEW_GAP), helium, True)
    if first.inside_helium != 0:
        raise RuntimeError("the rejected first factor meets helium")
    rows = positive_dressings(hit_name, hit_value, alpha, base)
    parent_name, parent_value, parent_gap = _nearest_row(rows, target)
    if parent_name != _REJECT_PARENT:
        raise RuntimeError(f"the rejected parent changed to {parent_name}")
    if parent_gap >= hit_gap:
        raise RuntimeError("the rejected parent did not tighten")
    parent = _record(
        rows, target, NamedGap(parent_name, parent_value, parent_gap, parent_gap <= _REVIEW_GAP), helium, True
    )
    if parent.inside_helium != 0 or parent.dressing.relative_gap <= helium:
        raise RuntimeError("the rejected parent meets helium")
    rows = positive_dressings(parent_name, parent_value, alpha, base)
    size = _positive_size(rows)
    if size.expressions != _COUNT or size.inside != _COUNT:
        raise RuntimeError("the rejected dressing count changed")
    if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
        raise RuntimeError("the rejected spacing changed")
    hits = [
        (row, number, _gap(number, target))
        for row, number in rows
        if _gap(number, target) <= helium
    ]
    if len(hits) != 2:
        raise RuntimeError("the rejected fraction left two helium dressings")
    if not parent_value > target:
        raise RuntimeError("the rejected parent left the high side of the fraction")
    expected = {
        f"({parent_name})*(1-alpha*e^-6)": "e^-6",
        f"({parent_name})*(1-alpha*pi^-5)": "pi^-5",
    }
    if {row for row, _number, _hit_gap in hits} != set(expected):
        raise RuntimeError("the rejected pair left e^-6 and pi^-5")
    found = {expected[row]: (row, number, hit_gap) for row, number, hit_gap in hits}
    e_name, e_value, e_gap = found["e^-6"]
    pi_name, pi_value, pi_gap = found["pi^-5"]
    if not pi_value < target < e_value:
        raise RuntimeError("the rejected fraction left the interval between the two dressings")
    if _num(pi_value) == _num(target) or _num(e_value) == _num(target):
        raise RuntimeError("a rejected neighbor printed the fraction")
    required = (target / parent_value - 1.0) / (-1.0 * alpha)
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
        raise RuntimeError("the rejected required factor left pi^-5")
    separation = abs(e_value - pi_value) / target
    piece_separation = alpha * abs(pieces["e^-6"] - pieces["pi^-5"])
    if _gap(separation, piece_separation * parent_value / target) > 1.0e-9:
        raise RuntimeError("the separation left alpha times the piece difference")
    if not helium < separation < 2.0 * helium:
        raise RuntimeError("the two dressings left the one-window overlap")
    nested = alpha * size.median_spacing
    if not nested < helium or _num(nested) != "3.7316e-06":
        raise RuntimeError("the dressed-factor spacing left the fine side of helium")
    e_child_name, _e_child_value, e_child_gap = _nearest_row(
        positive_dressings(e_name, e_value, alpha, base), target
    )
    pi_child_name, _pi_child_value, pi_child_gap = _nearest_row(
        positive_dressings(pi_name, pi_value, alpha, base), target
    )
    if e_child_name != _REJECT_E_CHILD or pi_child_name != _REJECT_PI_CHILD:
        raise RuntimeError("a rejected child changed")
    if e_child_gap <= e_gap or pi_child_gap <= pi_gap:
        raise RuntimeError("another factor on the rejected pair tightened")
    if e_child_gap <= helium or pi_child_gap <= helium:
        raise RuntimeError("another factor on the rejected pair meets helium")
    return RejectedPair(
        anchor=anchor,
        first=first,
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
        nested_spacing=nested,
        target=target,
    )


def build_filed() -> FiledReport:
    ceiling = build_ceiling()
    process = build_process()
    decay = build_decay()
    module, digest = load_engine()
    if digest[:16] != ceiling.pin_prefix or ceiling.digest_prefix != ceiling.pin_prefix:
        raise RuntimeError("the filed pin does not match the plant pin")
    if "Sodium" not in str(catalog_record(FAST_ID).get("engineering_barrier") or ""):
        raise RuntimeError("the fast barrier does not name sodium")
    hot = ceiling.hot_k
    cold = ceiling.cold_k
    if hot != 560.0 or cold != 310.0:
        raise RuntimeError("the ceiling temperatures left 560 K and 310 K")
    rejected_target = cold / hot
    if _gap(ceiling.carnot, 1.0 - rejected_target) > 0.0:
        raise RuntimeError("the temperature Carnot left one minus the rejected fraction")
    margin = ceiling.catalog_ceiling - ceiling.carnot
    if margin <= 0.0:
        raise RuntimeError("the catalog ceiling is not above the temperature Carnot")
    margin_j = ceiling.density_j * margin
    if _gap(margin_j, ceiling.catalog_margin_j) > 1.0e-12:
        raise RuntimeError("the round-file joules left the margin times the ranked density")
    if process.kwh_per_swu != 50.0:
        raise RuntimeError("the separative rate left 50 kWh per SWU")
    if decay.operation_s != 365.0 * 86400.0:
        raise RuntimeError("the filed year left 365 civil days")
    if _EXPONENT != -0.2:
        raise RuntimeError("the exponent magnitude left 0.2")
    if abs(_SHORT_OFFSET_S - 10.0) > 0.0 or abs(_LONG_OFFSET_S - 2.0e7) > 0.0:
        raise RuntimeError("the equation (2) offsets changed")
    day = 86400.0 / decay.operation_s
    cooldown = 20.0 * 3600.0 / decay.operation_s
    if _gap(2.0 * day, 2.0 * 86400.0 / decay.operation_s) > 0.0:
        raise RuntimeError("two days left twice the one-day share")
    if _gap((4.0 * 3600.0) / (20.0 * 3600.0), 0.2) > 0.0:
        raise RuntimeError("four hours over twenty hours left the exponent magnitude")
    tube = _kelvin(_TUBE_OUT_F) / _kelvin(_TUBE_IN_F)
    ends = _kelvin(_COOLDOWN_COLD_F) / _kelvin(_COOLDOWN_HOT_F)
    if not 0.0 < tube < 1.0 or not 0.0 < ends < 1.0:
        raise RuntimeError("a temperature ratio left the unit interval")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    five = _five_pieces(module)
    base = _base_pieces(module)
    day_share = _walk(day, _DAY_ANCHOR, _DAY_STEPS, alpha, helium, five, base, False)
    offset_ratio = _walk(
        _LONG_OFFSET_S / _SHORT_OFFSET_S, _OFFSET_ANCHOR, _OFFSET_STEPS, alpha, helium, five, base, True
    )
    if len(day_share.steps) != 4 or day_share.steps[-2].inside_helium != 0:
        raise RuntimeError("the one-day share met helium before the fourth factor")
    if len(offset_ratio.steps) != 4 or offset_ratio.steps[-2].inside_helium != 0:
        raise RuntimeError("the offset ratio met helium before the fourth factor")
    return FiledReport(
        margin=_walk(margin, _MARGIN_ANCHOR, _MARGIN_STEPS, alpha, helium, five, base, False),
        rejected=_rejected(rejected_target, alpha, helium, five, base),
        rate=_walk(process.kwh_per_swu, _RATE_ANCHOR, _RATE_STEPS, alpha, helium, five, base, True),
        day_share=day_share,
        cooldown_share=_walk(cooldown, _COOL_ANCHOR, _COOL_STEPS, alpha, helium, five, base, False),
        offset_ratio=offset_ratio,
        tube_ratio=_walk(tube, _TUBE_ANCHOR, _TUBE_STEPS, alpha, helium, five, base, True),
        cooldown_ratio=_walk(ends, _ENDS_ANCHOR, _ENDS_STEPS, alpha, helium, five, base, True),
        hot_k=hot,
        cold_k=cold,
        margin_j=margin_j,
        rate_kwh=process.kwh_per_swu,
        short_offset_s=_SHORT_OFFSET_S,
        long_offset_s=_LONG_OFFSET_S,
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=ceiling.nuclear_scalar,
        thermodynamics_scalar=ceiling.thermodynamics_scalar,
        gamma_fsot=ceiling.gamma_fsot,
        valve_beta=ceiling.valve_beta,
        pin_prefix=ceiling.pin_prefix,
        digest_prefix=ceiling.digest_prefix,
    )


def _step_lines(step: FiledStep, closing: bool) -> list[str]:
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


def _figure_lines(figure: FiledFigure, lead: str, role: str, fourth: bool) -> list[str]:
    lines = [
        lead,
        f"The five-atom nearest is `{figure.anchor.name}` at `{_num(figure.anchor.value)}`.",
        f"The anchor gap is `{_num(figure.anchor.relative_gap)}`.",
    ]
    if len(figure.steps) == 1:
        lines.append("The first factor closes that quotient.")
    for step in figure.steps[:-1]:
        lines.extend(_step_lines(step, False))
    if fourth:
        lines.append("Three factors stay wider than the helium quotient match.")
        lines.append("The fourth factor is one more base-piece factor of the same form.")
    lines.extend(_step_lines(figure.steps[-1], True))
    if _num(figure.steps[-1].dressing.value) == _num(figure.target):
        lines.append(f"The closing value agrees with `{_num(figure.target)}` at six figures.")
    lines.append("The closing value names that quotient.")
    lines.append(role)
    return lines


def _rejected_lines(pair: RejectedPair) -> list[str]:
    return [
        f"The cold temperature over the hot temperature is `{_num(pair.target)}`.",
        f"The five-atom nearest is `{pair.anchor.name}` at `{_num(pair.anchor.value)}`.",
        f"The anchor gap is `{_num(pair.anchor.relative_gap)}`.",
        *_step_lines(pair.first, False),
        *_step_lines(pair.parent, False),
        "That dressing is the parent of the rejected fraction.",
        f"The `e^-6` neighbor is `{pair.e_name}` at `{_num(pair.e_value)}`.",
        f"The `e^-6` neighbor gap is `{_num(pair.e_gap)}`.",
        f"The `pi^-5` neighbor is `{pair.pi_name}` at `{_num(pair.pi_value)}`.",
        f"The `pi^-5` neighbor gap is `{_num(pair.pi_gap)}`.",
        "The rejected fraction sits between those two values.",
        "Two dressings of that parent meet the helium quotient match.",
        f"Their separation is `{_num(pair.separation)}`.",
        f"alpha times the difference of `e^-6` and `pi^-5` is `{_num(pair.piece_separation)}`.",
        "The separation equals that product times the parent over the fraction.",
        "The separation is wider than one helium quotient match.",
        "The separation is narrower than two helium quotient matches.",
        f"The factor the rejected fraction requires is `{_num(pair.required)}`.",
        f"Its relative gap to `e^-6` is `{_num(pair.required_e_gap)}`.",
        f"Its relative gap to `pi^-5` is `{_num(pair.required_pi_gap)}`.",
        "Both factor gaps stay wider than the helium quotient match.",
        f"The nearest base piece is `{pair.nearest_piece}`.",
        "No base piece meets the required factor.",
        f"Another factor on the `e^-6` neighbor is `{pair.e_child_name}`.",
        f"Its relative gap is `{_num(pair.e_child_gap)}`.",
        f"Another factor on the `pi^-5` neighbor is `{pair.pi_child_name}`.",
        f"Its relative gap is `{_num(pair.pi_child_gap)}`.",
        "Another factor on either neighbor moves farther from the fraction.",
        f"Dressing the factor by another base piece scales the spacing by alpha to `{_num(pair.nested_spacing)}`.",
        "That spacing is finer than the helium quotient match.",
        "The dressed factor stays unused.",
        "Both neighbors stay uninstalled.",
        "The operating fraction stays the cold temperature over the hot temperature.",
        "The temperature Carnot is one minus the rejected fraction.",
        "Both routes stay for a later experiment or a later published figure sharper than their separation.",
        "This score makes no prediction that chooses between them.",
    ]


def render() -> str:
    report = build_filed()
    lines = [
        "# Fission filed",
        "",
        "The ceiling file leaves the round-file margin and the rejected temperature fraction.",
        "The material balance leaves the separative rate.",
        "The filed year leaves the one-day share and the twenty-hour share.",
        "Equation (2) leaves the ratio of its two offsets.",
        "The residual-heat table leaves two absolute temperature ratios.",
        "FSOT names each quotient that meets the helium standard on one dressing.",
        "The rejected fraction keeps both routes.",
        "The chain stays at k_eff of 1.",
        "",
        "## Margin",
        "",
        f"The steam cycle is `{_num(report.hot_k)}` K rejected near `{_num(report.cold_k)}` K.",
        "The catalog ceiling stays above the temperature Carnot.",
        "The catalog ceiling stays unrun.",
    ]
    lines.extend(
        _figure_lines(
            report.margin,
            f"The round-file margin is `{_num(report.margin.target)}`.",
            "The operating margin stays the catalog ceiling minus the temperature Carnot.",
            False,
        )
    )
    lines.extend(
        [
            f"The extra joules are `{_num(report.margin_j)}` J.",
            "The extra joules are that margin times the ranked density.",
            "",
            "## Rejection",
            "",
            "The unit-interval dressing of the rejected fraction drops below 116 values.",
            "The positive dressing compares it.",
        ]
    )
    lines.extend(_rejected_lines(report.rejected))
    lines.extend(
        [
            "",
            "## Separative rate",
            "",
            "Ratios that leave the unit interval use the positive dressing from the heat-ratio score.",
            "The separative-work estimate stays on the enrichment basis.",
        ]
    )
    lines.extend(
        _figure_lines(
            report.rate,
            f"The material-balance rate is `{_num(report.rate.target)}` kWh per SWU.",
            "The operating rate stays 50 kWh per SWU.",
            False,
        )
    )
    lines.extend(
        [
            "",
            "## Year",
            "",
        ]
    )
    lines.extend(
        _figure_lines(
            report.day_share,
            f"One day over the filed year is `{_num(report.day_share.target)}`.",
            "The operating share stays one day over the filed year.",
            True,
        )
    )
    lines.extend(
        [
            "Two days over the filed year is twice the one-day share.",
            "",
            "## Cooldown clock",
            "",
            "Four hours over twenty hours stays the exponent magnitude.",
        ]
    )
    lines.extend(
        _figure_lines(
            report.cooldown_share,
            f"Twenty hours over the filed year is `{_num(report.cooldown_share.target)}`.",
            "The operating share stays twenty hours over the filed year.",
            False,
        )
    )
    lines.extend(
        [
            "",
            "## Offsets",
            "",
            f"The short offset stays `{_num(report.short_offset_s)}` seconds.",
            f"The long offset stays `{_num(report.long_offset_s)}` seconds.",
            "Equation (2) stays the witness.",
            "This score does not replace equation (2).",
        ]
    )
    lines.extend(
        _figure_lines(
            report.offset_ratio,
            f"The long offset over the short offset is `{_num(report.offset_ratio.target)}`.",
            "The operating ratio stays the long offset over the short offset.",
            True,
        )
    )
    lines.extend(
        [
            "",
            "## Tubes",
            "",
            "The table's tube-side temperatures are 140°F in and 124°F out.",
            "The unit-interval dressing of that ratio drops below 116 values.",
            "The positive dressing names it.",
        ]
    )
    lines.extend(
        _figure_lines(
            report.tube_ratio,
            f"The outlet over the inlet, on the absolute scale, is `{_num(report.tube_ratio.target)}`.",
            "The operating ratio stays the outlet temperature over the inlet temperature.",
            False,
        )
    )
    lines.extend(
        [
            "",
            "## Cooldown temperatures",
            "",
            "The coolant falls from 350°F to 140°F.",
            "The UFSAR does not print a watt rate at the 350°F entry.",
            "This score assigns that hotter rate no figure.",
            "The unit-interval dressing of that ratio drops below 116 values.",
            "The positive dressing names it.",
        ]
    )
    lines.extend(
        _figure_lines(
            report.cooldown_ratio,
            f"The cold end over the 350°F entry, on the absolute scale, is `{_num(report.cooldown_ratio.target)}`.",
            "The operating ratio stays the cold-end temperature over the cooldown entry.",
            False,
        )
    )
    lines.extend(
        [
            "The fast ceiling source states no temperature pair.",
            "The fast row's filed barrier is the sodium system.",
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
