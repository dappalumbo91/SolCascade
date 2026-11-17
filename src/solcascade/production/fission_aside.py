"""Name the overview figures the operating column left aside.

The comparison column sits on the fabrication cap and stays unrun. The
Red Book assays stay off that column. The CANDU row stays a different
plant. A shared number keeps the expression it already has. The
checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_balance import delayed_proportions
from solcascade.production.fission_decay_dress import _base_pieces
from solcascade.production.fission_decay_dress import _dressings as unit_dressings
from solcascade.production.fission_decay_heat import _dressings as positive_dressings
from solcascade.production.fission_decay_heat import _nearest_row, _positive_nearest, _positive_size
from solcascade.production.fission_decay_seed import NamedGap, _five_pieces, _nearest, _size
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.fission_process import (
    ENRICHED_PERCENT,
    FABRICATION_CAP_PERCENT,
    FULL_POWER_HOURS,
    HIGHER_BURNUP_GWD,
    MATERIAL_BALANCE_SOURCE,
    NAMEPLATE_MWE,
    RELOAD_TONNES,
)
from solcascade.production.fission_tails import TAILS_PERCENT
from solcascade.production.fission_safety import BUILT_EFFICIENCY, SMALL_CLASS_MWE
from solcascade.production.plants import FAST_ID, THERMAL_ID, build_plants
from solcascade.production.refine import catalog_record, project_root

_COUNT = 116
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4
_RATIONAL_GAP = 1.0e-14
_NATURAL_T = 211.0
_COMPARE_NATURAL_T = 163.0
_COMPARE_PRODUCT_T = 16.85
_COMPARE_CONCENTRATE_T = 192.0
_COMPARE_HEX_T = 241.0
_COMPARE_ENRICHED_HEX_T = 25.0
_COMPARE_OXIDE_T = 19.1
_COMPARE_TRANSURANIC_T = 0.200
_COMPARE_USED_OXIDE_T = 18.3
_COMPARE_BURNUP_GWD = 65.0
_OPERATING_GWD = 45.0
_PAGE_PER_TWH = 18.6
_RED_HIGH_PERCENT = 0.30
_RED_LOW_PERCENT = 0.25
_RED_HIGH_TU = 175.0
_RED_LOW_TU = 160.0
_RED_HIGH_OXIDE = 206.0
_RED_LOW_OXIDE = 189.0
_CANDU_U238_PERCENT = 98.8
_CANDU_U235_PERCENT = 0.23
_CANDU_PU_PERCENT = 0.38
_CANDU_PRODUCT_PERCENT = 0.6
_CANDU_EFFICIENCY = 0.31
_CANDU_BURNUP_GWD = 7.5
_CANDU_DEMAND_T = 157.0
_CANDU_CAPACITY_T = 142.0
_CANDU_PAGE_PER_TWH = 17.9
_PAGE_SHARE = 0.80
_LOAD_FACTOR = 0.90
_LOAD_EXAMPLE_T = 190.0
_LOAD_MODERN_T = 147.0
_LOAD_TWH = 7.9


@dataclass(frozen=True)
class _Spec:
    label: str
    rational: Fraction
    positive: bool
    unit_rows: int | None
    anchor: str
    steps: tuple[str, ...]


_SPECS = (
    _Spec(
        "comparison concentrate",
        Fraction(163, 192),
        True,
        108,
        "pi^-4+G^2",
        (
            "(pi^-4+G^2)*(1-alpha*e^-3)",
            "((pi^-4+G^2)*(1-alpha*e^-3))*(1+alpha*pi^-6)",
        ),
    ),
    _Spec(
        "comparison hex",
        Fraction(163, 241),
        True,
        111,
        "pi^-3+G^5",
        (
            "(pi^-3+G^5)*(1-alpha*e^-2)",
            "((pi^-3+G^5)*(1-alpha*e^-2))*(1+alpha*e^-6)",
            "(((pi^-3+G^5)*(1-alpha*e^-2))*(1+alpha*e^-6))*(1-alpha*pi^-6)",
        ),
    ),
    _Spec(
        "comparison enriched hex",
        Fraction(337, 500),
        True,
        111,
        "sqrt(G)*G^4",
        ("(sqrt(G)*G^4)*(1+alpha*gamma^5)",),
    ),
    _Spec(
        "comparison oxide",
        Fraction(337, 382),
        True,
        107,
        "phi^-3+G^5",
        (
            "(phi^-3+G^5)*(1+alpha*gamma^3)",
            "((phi^-3+G^5)*(1+alpha*gamma^3))*(1+alpha*e^-4)",
            "(((phi^-3+G^5)*(1+alpha*gamma^3))*(1+alpha*e^-4))*(1+alpha*pi^-5)",
        ),
    ),
    _Spec(
        "comparison per terawatt-hour",
        Fraction(4075, 219),
        True,
        None,
        "gamma^-2+gamma^-5",
        (
            "(gamma^-2+gamma^-5)*(1-alpha*e^-5)",
            "((gamma^-2+gamma^-5)*(1-alpha*e^-5))*(1+alpha*pi^-6)",
        ),
    ),
    _Spec(
        "comparison transuranic",
        Fraction(2, 191),
        False,
        _COUNT,
        "e^-5*G^-5",
        (
            "(e^-5*G^-5)*(1+alpha*phi^-3)",
            "((e^-5*G^-5)*(1+alpha*phi^-3))*(1+alpha*gamma^6)",
        ),
    ),
    _Spec(
        "comparison used oxide",
        Fraction(183, 191),
        True,
        100,
        "pi^-6+sqrt(G)",
        ("(pi^-6+sqrt(G))*(1+alpha*e^-6)",),
    ),
    _Spec(
        "red book tails percent",
        Fraction(1, 4),
        True,
        115,
        "gamma^3*G^-3",
        (
            "(gamma^3*G^-3)*(1-alpha*e^-2)",
            "((gamma^3*G^-3)*(1-alpha*e^-2))*(1-alpha*e^-6)",
        ),
    ),
    _Spec(
        "red book demand",
        Fraction(175),
        True,
        None,
        "e^5+gamma^-6",
        (
            "(e^5+gamma^-6)*(1-alpha*gamma^2)",
            "((e^5+gamma^-6)*(1-alpha*gamma^2))*(1-alpha*e^-4)",
            "(((e^5+gamma^-6)*(1-alpha*gamma^2))*(1-alpha*e^-4))*(1-alpha*pi^-6)",
        ),
    ),
    _Spec(
        "red book lower demand",
        Fraction(160),
        True,
        None,
        "e^5+phi^5",
        (
            "(e^5+phi^5)*(1+alpha*phi^-2)",
            "((e^5+phi^5)*(1+alpha*phi^-2))*(1+alpha*e^-3)",
            "(((e^5+phi^5)*(1+alpha*phi^-2))*(1+alpha*e^-3))*(1-alpha*e^-5)",
            "((((e^5+phi^5)*(1+alpha*phi^-2))*(1+alpha*e^-3))*(1-alpha*e^-5))*(1+alpha*pi^-6)",
        ),
    ),
    _Spec(
        "red book oxide",
        Fraction(206),
        True,
        None,
        "e^4+e^5",
        (
            "(e^4+e^5)*(1+alpha*sqrt(pi))",
            "((e^4+e^5)*(1+alpha*sqrt(pi)))*(1+alpha*phi^-3)",
            "(((e^4+e^5)*(1+alpha*sqrt(pi)))*(1+alpha*phi^-3))*(1+alpha*e^-5)",
            "((((e^4+e^5)*(1+alpha*sqrt(pi)))*(1+alpha*phi^-3))*(1+alpha*e^-5))*(1-alpha*pi^-6)",
        ),
    ),
    _Spec(
        "red book lower oxide",
        Fraction(189),
        True,
        None,
        "pi^5/phi",
        (
            "(pi^5/phi)*(1-alpha*phi^-5)",
            "((pi^5/phi)*(1-alpha*phi^-5))*(1-alpha*pi^-5)",
            "(((pi^5/phi)*(1-alpha*phi^-5))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
        ),
    ),
    _Spec(
        "candu uranium 238",
        Fraction(494, 5),
        True,
        None,
        "pi^4+G^-4",
        (
            "(pi^4+G^-4)*(1-alpha*gamma^6)",
            "((pi^4+G^-4)*(1-alpha*gamma^6))*(1-alpha*pi^-5)",
            "(((pi^4+G^-4)*(1-alpha*gamma^6))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
        ),
    ),
    _Spec(
        "candu plutonium",
        Fraction(19, 50),
        True,
        113,
        "G^5*G^6",
        (
            "(G^5*G^6)*(1-alpha*phi^-3)",
            "((G^5*G^6)*(1-alpha*phi^-3))*(1-alpha*e^-3)",
            "(((G^5*G^6)*(1-alpha*phi^-3))*(1-alpha*e^-3))*(1+alpha*e^-5)",
            "((((G^5*G^6)*(1-alpha*phi^-3))*(1-alpha*e^-3))*(1+alpha*e^-5))*(1-alpha*pi^-6)",
        ),
    ),
    _Spec(
        "candu efficiency",
        Fraction(31, 100),
        True,
        113,
        "phi^-3/sqrt(gamma)",
        (
            "(phi^-3/sqrt(gamma))*(1-alpha*gamma^2)",
            "((phi^-3/sqrt(gamma))*(1-alpha*gamma^2))*(1+alpha*e^-4)",
            "(((phi^-3/sqrt(gamma))*(1-alpha*gamma^2))*(1+alpha*e^-4))*(1-alpha*e^-6)",
        ),
    ),
    _Spec(
        "candu demand",
        Fraction(157),
        True,
        None,
        "e^5+gamma^-4",
        (
            "(e^5+gamma^-4)*(1-alpha*phi^-2)",
            "((e^5+gamma^-4)*(1-alpha*phi^-2))*(1+alpha*e^-4)",
            "(((e^5+gamma^-4)*(1-alpha*phi^-2))*(1+alpha*e^-4))*(1-alpha*pi^-5)",
        ),
    ),
    _Spec(
        "candu capacity demand",
        Fraction(142),
        True,
        None,
        "e^5*sqrt(G)",
        (
            "(e^5*sqrt(G))*(1-alpha*gamma^6)",
            "((e^5*sqrt(G))*(1-alpha*gamma^6))*(1-alpha*e^-6)",
        ),
    ),
)


@dataclass(frozen=True)
class AsideStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class AsideFigure:
    label: str
    anchor: NamedGap
    steps: tuple[AsideStep, ...]
    target: float
    rational: Fraction
    positive: bool
    unit_rows: int | None


@dataclass(frozen=True)
class AsideReport:
    figures: tuple[AsideFigure, ...]
    compare_burnup_step: float
    candu_burnup_step: float
    red_high_share: float
    red_low_share: float
    candu_per_twh: float
    page_per_twh_gap: float
    candu_page_gap: float
    load_example_gap: float
    helium_mismatch: float
    alpha: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str
    source: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_aside.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _gap(derived: float, witness: float) -> float:
    if witness == 0.0:
        raise RuntimeError("the witness is zero")
    return abs(derived - witness) / abs(witness)


def _require(value: float, rational: Fraction, label: str) -> None:
    if _gap(value, float(rational)) > _RATIONAL_GAP:
        raise RuntimeError(f"the {label} left {rational.numerator}/{rational.denominator}")


def _second_gap(rows, target: float, best: float) -> float:
    gaps = sorted(_gap(value, target) for _name, value in rows)
    if len(gaps) < 2 or abs(gaps[0] - best) > 1.0e-12:
        raise RuntimeError("the nearest dressing is missing from its own list")
    return gaps[1]


def _unit_count(name: str, value: float, alpha: float, base) -> int | None:
    try:
        return len(unit_dressings(name, value, alpha, base))
    except RuntimeError:
        return None


def _rows_for(spec: _Spec, name: str, value: float, alpha: float, base):
    counted = _unit_count(name, value, alpha, base)
    if counted != spec.unit_rows:
        raise RuntimeError(f"the {spec.label} unit dressing left {spec.unit_rows}")
    if spec.positive:
        rows = positive_dressings(name, value, alpha, base)
        hit_name, hit_value, hit_gap = _nearest_row(rows, float(spec.rational))
        hit = NamedGap(hit_name, hit_value, hit_gap, hit_gap <= _REVIEW_GAP)
        size = _positive_size(rows)
    else:
        rows = unit_dressings(name, value, alpha, base)
        hit = _nearest(rows, float(spec.rational))
        size = _size(rows)
    return rows, hit, size


def _walk(spec: _Spec, alpha: float, helium: float, five, base) -> AsideFigure:
    target = float(spec.rational)
    if spec.positive:
        name, value, anchor_gap = _positive_nearest(five, target)
    else:
        anchor = _nearest(five, target)
        name, value, anchor_gap = anchor.name, anchor.value, anchor.relative_gap
    if name != spec.anchor:
        raise RuntimeError(f"the {spec.label} anchor changed to {name}")
    if anchor_gap <= helium:
        raise RuntimeError(f"the {spec.label} anchor already meets helium")
    anchor_named = NamedGap(name, value, anchor_gap, anchor_gap <= _REVIEW_GAP)
    steps: list[AsideStep] = []
    for step_name in spec.steps:
        rows, hit, size = _rows_for(spec, name, value, alpha, base)
        if size.expressions != _COUNT or size.inside != _COUNT:
            raise RuntimeError(f"the {spec.label} dressing count changed")
        if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
            raise RuntimeError(f"the {spec.label} spacing is not the base-piece spacing")
        if hit.name != step_name:
            raise RuntimeError(f"the {spec.label} dressing changed to {hit.name}")
        if hit.relative_gap >= _gap(value, target):
            raise RuntimeError(f"the {spec.label} dressing did not tighten")
        inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
        if inside_helium > 1:
            raise RuntimeError(f"several {spec.label} dressings meet the helium standard")
        steps.append(
            AsideStep(
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
        raise RuntimeError(f"the closing {spec.label} misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError(f"a second {spec.label} dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError(f"an earlier {spec.label} figure already meets helium")
    if len(steps) > 4:
        raise RuntimeError(f"the {spec.label} added a fifth factor")
    return AsideFigure(
        spec.label,
        anchor_named,
        tuple(steps),
        target,
        spec.rational,
        spec.positive,
        spec.unit_rows,
    )


def _by_label(report: AsideReport, label: str) -> AsideFigure:
    for figure in report.figures:
        if figure.label == label:
            return figure
    raise RuntimeError(f"the {label} left the aside report")


def build_aside() -> AsideReport:
    plants = build_plants()
    module, digest = load_engine()
    if digest[:16] != plants.pin_prefix:
        raise RuntimeError("the aside pin does not match the plant pin")
    record = catalog_record(THERMAL_ID)
    if "k_eff is held at 1" not in str(record["governing_equation"]):
        raise RuntimeError("the governing equation does not hold k_eff at 1")
    if "Sodium" not in str(catalog_record(FAST_ID).get("engineering_barrier") or ""):
        raise RuntimeError("the fast barrier does not name sodium")
    if RELOAD_TONNES != 24.3 or ENRICHED_PERCENT != 4.5 or TAILS_PERCENT != 0.22:
        raise RuntimeError("the operating column left 24.3 t, 4.5 percent, and 0.22 percent tails")
    if FABRICATION_CAP_PERCENT != 5.0 or HIGHER_BURNUP_GWD != 60.0:
        raise RuntimeError("the cap left 5 percent or the stack left 60 GWd/t")
    if BUILT_EFFICIENCY != 0.36:
        raise RuntimeError("the plants-being-built efficiency left 0.36")
    if SMALL_CLASS_MWE / NAMEPLATE_MWE != _RED_HIGH_PERCENT:
        raise RuntimeError("the 0.30 percent assay left the small-class share")
    delayed = {row.nucleus: row for row in delayed_proportions()}
    if delayed["Pu-239"].percent != _CANDU_U235_PERCENT:
        raise RuntimeError("the 0.23 percent figure left the delayed plutonium count")
    if delayed["Pu-239"].proportion != _CANDU_U235_PERCENT / 100.0:
        raise RuntimeError("the 0.23 percent fraction left the delayed plutonium proportion")
    observables = (
        project_root() / "src" / "solcascade" / "production" / "fission_observables.py"
    ).read_text(encoding="utf-8")
    if "The channel search stays closed." not in observables:
        raise RuntimeError("the channel sentence changed")
    compare_step = _COMPARE_BURNUP_GWD / _OPERATING_GWD
    candu_step = _CANDU_BURNUP_GWD / _OPERATING_GWD
    _require(compare_step, Fraction(13, 9), "comparison burn-up step")
    _require(candu_step, Fraction(1, 6), "candu burn-up step")
    per_twh = _COMPARE_NATURAL_T * 1000.0 / FULL_POWER_HOURS
    _require(per_twh, Fraction(4075, 219), "comparison uranium per terawatt-hour")
    if Fraction(_COMPARE_NATURAL_T) * 1000 / Fraction(8760) != Fraction(4075, 219):
        raise RuntimeError("the comparison year left 4075/219")
    candu_per_twh = _CANDU_DEMAND_T * 1000.0 / FULL_POWER_HOURS
    _require(candu_per_twh, Fraction(3925, 219), "candu uranium per terawatt-hour")
    red_high_share = _RED_HIGH_TU / _RED_HIGH_OXIDE
    red_low_share = _RED_LOW_TU / _RED_LOW_OXIDE
    _require(red_high_share, Fraction(175, 206), "red book high oxide share")
    _require(red_low_share, Fraction(160, 189), "red book low oxide share")
    page_gap = _gap(_PAGE_PER_TWH, per_twh)
    candu_gap = _gap(_CANDU_PAGE_PER_TWH, candu_per_twh)
    if page_gap > _REVIEW_GAP:
        raise RuntimeError("the page's 18.6 left the review bar")
    if candu_gap <= _REVIEW_GAP:
        raise RuntimeError("the page's 17.9 met the exact quotient inside the review bar")
    if _gap(_PAGE_SHARE, _CANDU_DEMAND_T / _NATURAL_T) <= _REVIEW_GAP:
        raise RuntimeError("the page's 80 percent met 157 over 211 inside the review bar")
    load_gap = _gap(_LOAD_EXAMPLE_T, _LOAD_FACTOR * _NATURAL_T)
    if load_gap > _REVIEW_GAP:
        raise RuntimeError("the example 190 tU left the review bar")
    if _gap(_LOAD_MODERN_T, _LOAD_FACTOR * _COMPARE_NATURAL_T) <= _REVIEW_GAP:
        raise RuntimeError("the example 147 tU met 90 percent of 163 inside the review bar")
    if _gap(_LOAD_TWH, _LOAD_FACTOR * (NAMEPLATE_MWE * FULL_POWER_HOURS / 1_000_000.0)) <= _REVIEW_GAP:
        raise RuntimeError("the example 7.9 TWh met 90 percent of the nameplate year")
    product_gap = _gap(100.0 - _CANDU_U238_PERCENT - _CANDU_U235_PERCENT - _CANDU_PU_PERCENT, _CANDU_PRODUCT_PERCENT)
    if product_gap <= _REVIEW_GAP:
        raise RuntimeError("the 0.6 percent fission products closed the footnote balance")
    if _LOAD_FACTOR * _CANDU_DEMAND_T == _CANDU_CAPACITY_T:
        raise RuntimeError("the 142 t figure became 90 percent of 157 t")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    five = _five_pieces(module)
    base = _base_pieces(module)
    figures = tuple(_walk(spec, alpha, helium, five, base) for spec in _SPECS)
    if any(figure.target in (0.30, 0.0023, 7.5, 18.6, 17.9, 0.80, 0.90) for figure in figures):
        raise RuntimeError("a kept or round figure was installed")
    return AsideReport(
        figures=figures,
        compare_burnup_step=compare_step,
        candu_burnup_step=candu_step,
        red_high_share=red_high_share,
        red_low_share=red_low_share,
        candu_per_twh=candu_per_twh,
        page_per_twh_gap=page_gap,
        candu_page_gap=candu_gap,
        load_example_gap=load_gap,
        helium_mismatch=helium,
        alpha=alpha,
        nuclear_scalar=plants.nuclear_scalar,
        thermodynamics_scalar=plants.thermodynamics_scalar,
        gamma_fsot=plants.gamma_fsot,
        valve_beta=plants.valve_beta,
        pin_prefix=plants.pin_prefix,
        digest_prefix=plants.digest_prefix,
        source=MATERIAL_BALANCE_SOURCE,
    )


def _step_lines(step: AsideStep, closing: bool) -> list[str]:
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


def _figure_lines(figure: AsideFigure, lead: str, rational: str, role: str) -> list[str]:
    lines = [lead]
    if rational:
        lines.append(f"The same quotient is `{rational}`.")
    if figure.unit_rows is None:
        lines.append("The unit-interval dressing is empty.")
        lines.append("The positive dressing names it.")
    elif figure.positive:
        lines.append(f"The unit-interval dressing has `{figure.unit_rows}` values.")
        lines.append("The positive dressing names it.")
    else:
        lines.append("The unit dressing names it.")
    lines.extend(
        [
            f"The five-atom nearest is `{figure.anchor.name}` at `{_num(figure.anchor.value)}`.",
            f"The anchor gap is `{_num(figure.anchor.relative_gap)}`.",
        ]
    )
    if figure.anchor.relative_gap <= _REVIEW_GAP:
        lines.append("The five-atom nearest stays inside the review bar.")
        lines.append("The five-atom nearest is not the installed close.")
    for step in figure.steps[:-1]:
        lines.extend(_step_lines(step, False))
    lines.extend(_step_lines(figure.steps[-1], True))
    if _num(figure.steps[-1].dressing.value) == _num(figure.target):
        lines.append(f"The closing value agrees with `{_num(figure.target)}` at six figures.")
    else:
        lines.append(f"The closing print is `{_num(figure.steps[-1].dressing.value)}`.")
    lines.append(role)
    return lines


def render() -> str:
    report = build_aside()
    concentrate = _by_label(report, "comparison concentrate")
    hex_share = _by_label(report, "comparison hex")
    enriched_hex = _by_label(report, "comparison enriched hex")
    oxide = _by_label(report, "comparison oxide")
    per_twh = _by_label(report, "comparison per terawatt-hour")
    transuranic = _by_label(report, "comparison transuranic")
    used_oxide = _by_label(report, "comparison used oxide")
    tails = _by_label(report, "red book tails percent")
    demand = _by_label(report, "red book demand")
    lower_demand = _by_label(report, "red book lower demand")
    red_oxide = _by_label(report, "red book oxide")
    red_lower = _by_label(report, "red book lower oxide")
    uranium_238 = _by_label(report, "candu uranium 238")
    plutonium = _by_label(report, "candu plutonium")
    efficiency = _by_label(report, "candu efficiency")
    candu_demand = _by_label(report, "candu demand")
    capacity = _by_label(report, "candu capacity demand")
    lines = [
        "# Fission aside",
        "",
        f"The material balance is {report.source}.",
        "The operating column stays 211 t of natural uranium for 24.3 t at 4.5 percent, 0.22 percent tails, 45 GWd/t, and 33 percent.",
        "The chain stays at k_eff of 1.",
        "",
        "## Comparison column",
        "",
        "The comparison column is the 5 percent enrichment at 65 GWd/t.",
        "That enrichment is the fabrication cap.",
        "This comparison stays unrun.",
        f"The concentrate line is `{_num(_COMPARE_CONCENTRATE_T)}` t containing `{_num(_COMPARE_NATURAL_T)}` t of uranium.",
    ]
    lines.extend(
        _figure_lines(
            concentrate,
            f"The uranium share of that concentrate is `{_num(concentrate.target)}`.",
            "163/192",
            "That share keeps the printed comparison concentrate.",
        )
    )
    lines.extend(
        [
            f"The conversion line is `{_num(_COMPARE_HEX_T)}` t of uranium hexafluoride.",
        ]
    )
    lines.extend(
        _figure_lines(
            hex_share,
            f"The uranium share of that hexafluoride is `{_num(hex_share.target)}`.",
            "163/241",
            "That share keeps the printed comparison conversion.",
        )
    )
    lines.extend([f"The enrichment line is `{_num(_COMPARE_ENRICHED_HEX_T)}` t of enriched hexafluoride."])
    lines.extend(
        _figure_lines(
            enriched_hex,
            f"The uranium share of that enriched hexafluoride is `{_num(enriched_hex.target)}`.",
            "337/500",
            "That share keeps the printed 16.85 t and 25 t.",
        )
    )
    lines.extend([f"The fabrication line is `{_num(_COMPARE_OXIDE_T)}` t of uranium dioxide."])
    lines.extend(
        _figure_lines(
            oxide,
            f"The uranium share of that dioxide is `{_num(oxide.target)}`.",
            "337/382",
            "That share keeps the printed 16.85 t and 19.1 t.",
        )
    )
    lines.extend(
        _figure_lines(
            per_twh,
            f"Natural uranium over the nameplate year is `{_num(per_twh.target)}` t per terawatt-hour.",
            "4075/219",
            "That quotient keeps the printed 163 t and the nameplate year.",
        )
    )
    lines.extend(
        [
            "The page's 18.6 tonnes per terawatt-hour sits inside the review bar.",
            "The page's 18.6 is not a second target.",
            f"The used-fuel line is `{_num(_COMPARE_OXIDE_T)}` t containing `{_num(_COMPARE_TRANSURANIC_T * 1000.0)}` kg of transuranics.",
        ]
    )
    lines.extend(
        _figure_lines(
            transuranic,
            f"The transuranic share of that used fuel is `{_num(transuranic.target)}`.",
            "2/191",
            "That share keeps the printed 200 kg and 19.1 t.",
        )
    )
    lines.extend(
        _figure_lines(
            used_oxide,
            f"The uranium-oxide share of that used fuel is `{_num(used_oxide.target)}`.",
            "183/191",
            "That share keeps the printed 18.3 t and 19.1 t.",
        )
    )
    lines.extend(
        [
            f"The comparison burn-up over the operating discharge is `{_num(report.compare_burnup_step)}`.",
            "The same quotient is `13/9`.",
            "That burn-up stays unrun.",
            "The comparison tails assay stays the named 0.22 percent.",
            "The enrichment note's about 140000 SWU stays approximate.",
            "The mining line stays a range.",
            "The feed multiple stays the named rational `3260/337`.",
            "",
            "## Red Book",
            "",
            "The Red Book tails assays are 0.30 percent and 0.25 percent.",
            "Those assays stay off this column.",
            "The 0.30 percent assay keeps the named small-class share.",
            "The 0.30 percent fraction is `3/1000`.",
            "That fraction keeps the small-class share.",
        ]
    )
    lines.extend(
        _figure_lines(
            tails,
            f"The 0.25 percent assay is `{_num(tails.target)}`.",
            "1/4",
            "That assay stays off the operating column.",
        )
    )
    lines.extend(
        [
            "The 0.25 percent fraction is `1/400`.",
            "That fraction keeps this percent.",
            f"The demand at the 0.30 percent assay is `{_num(demand.target)}` t of uranium per gigawatt-year.",
        ]
    )
    lines.extend(
        _figure_lines(
            demand,
            f"The Red Book demand is `{_num(demand.target)}` tU.",
            "",
            "That demand stays off the operating column.",
        )
    )
    lines.extend(
        _figure_lines(
            lower_demand,
            f"The demand at the 0.25 percent assay is `{_num(lower_demand.target)}` tU.",
            "",
            "That lower demand stays off the operating column.",
        )
    )
    lines.extend(
        _figure_lines(
            red_oxide,
            f"The corresponding oxide at the 0.30 percent assay is `{_num(red_oxide.target)}` t.",
            "",
            "That oxide stays off the operating column.",
        )
    )
    lines.extend(
        _figure_lines(
            red_lower,
            f"The corresponding oxide at the 0.25 percent assay is `{_num(red_lower.target)}` t.",
            "",
            "That lower oxide stays off the operating column.",
        )
    )
    lines.extend(
        [
            f"The uranium share of the `{_num(_RED_HIGH_OXIDE)}` t oxide is `{_num(report.red_high_share)}`.",
            "The same quotient is `175/206`.",
            "That share keeps the printed Red Book pair.",
            f"The uranium share of the `{_num(_RED_LOW_OXIDE)}` t oxide is `{_num(report.red_low_share)}`.",
            "The same quotient is `160/189`.",
            "That share keeps the printed lower pair.",
            "",
            "## CANDU",
            "",
            "The CANDU row stays a different plant.",
            "The footnote composition is after 7.5 GWd/t.",
        ]
    )
    lines.extend(
        _figure_lines(
            uranium_238,
            f"The uranium-238 percent is `{_num(uranium_238.target)}`.",
            "494/5",
            "That percent keeps the footnote.",
        )
    )
    lines.extend(
        [
            "The 0.23 percent figure is the delayed plutonium count.",
            "Its fraction keeps the named plutonium margin.",
            "This score does not add a second expression for that count.",
        ]
    )
    lines.extend(
        _figure_lines(
            plutonium,
            f"The plutonium percent is `{_num(plutonium.target)}`.",
            "19/50",
            "That percent keeps the footnote.",
        )
    )
    lines.extend(
        [
            "The 0.38 percent fraction is `19/5000`.",
            "That fraction keeps this percent.",
            "The 0.6 percent fission products stay one digit.",
            "The footnote balance does not close on that digit.",
        ]
    )
    lines.extend(
        _figure_lines(
            efficiency,
            f"The stated thermal efficiency is `{_num(efficiency.target)}`.",
            "31/100",
            "That efficiency stays on the CANDU row.",
        )
    )
    lines.extend(
        [
            f"The 7.5 GWd/t burn-up over the operating discharge is `{_num(report.candu_burnup_step)}`.",
            "The same quotient is `1/6`.",
            "The overview's about 7.5 stays that figure.",
            "The almost 50 GWd/t equivalent stays prose.",
        ]
    )
    lines.extend(
        _figure_lines(
            candu_demand,
            f"The stated demand is `{_num(candu_demand.target)}` tU.",
            "",
            "That demand stays on the CANDU row.",
        )
    )
    lines.extend(
        [
            f"That demand over the nameplate year is `{_num(report.candu_per_twh)}` t per terawatt-hour.",
            "The same quotient is `3925/219`.",
            "The page's 17.9 tonnes per terawatt-hour stays the page's round.",
        ]
    )
    lines.extend(
        _figure_lines(
            capacity,
            f"The 90 percent capacity figure is `{_num(capacity.target)}` tU.",
            "",
            "That figure is not 90 percent of 157 t.",
        )
    )
    lines.extend(
        [
            "The page's 80 percent stays that round.",
            "The 90 percent load factor stays the page's round.",
            "The example 190 tU sits inside the review bar of 90 percent of the operating natural uranium.",
            "The example 190 tU is not a second target.",
            "The example 147 tU stays the page's example.",
            "The example 7.9 TWh stays the page's round.",
            "",
            "## Unrun rows",
            "",
            "The catalog ceiling stays unrun.",
            "The catalog ceiling stays above the temperature Carnot.",
            "The 6.5 percent page stays unrun.",
            "The 60 GWd/t stack stays unrun.",
            "The 55 GWd/t step stays unrun.",
            "The 70 GWd/t step stays unrun.",
            "The comparison column stays unrun.",
            "The fabrication cap stays 5 percent.",
            "The fast row's filed barrier is the sodium system.",
            "The demonstrated plant stays at its filed efficiency.",
            "This score adds no joules.",
            "The published stamp records these closes.",
            "The checked-in schedule is not rewritten.",
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
            "The channel search stays closed.",
            "This run does not install a seed as a new operating plant.",
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
