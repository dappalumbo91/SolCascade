"""Name the operating material-balance quotients on the 1000 MWe example.

The overview's operating column prints 211 t of natural uranium for 24.3 t
of enriched uranium at 4.5 percent and 0.22 percent tails. The comparison
column at the fabrication cap and 65 GWd/t stays unrun. The checked-in
pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay_dress import _base_pieces
from solcascade.production.fission_decay_dress import _dressings as unit_dressings
from solcascade.production.fission_decay_heat import _dressings as positive_dressings
from solcascade.production.fission_decay_heat import _nearest_row, _positive_nearest, _positive_size
from solcascade.production.fission_decay_seed import NamedGap, _five_pieces, _nearest, _size
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.fission_process import (
    ANNUAL_SWU,
    ENRICHED_PERCENT,
    FABRICATION_CAP_PERCENT,
    FULL_POWER_HOURS,
    HIGHER_BURNUP_GWD,
    MATERIAL_BALANCE_SOURCE,
    NAMEPLATE_MWE,
    RELOAD_TONNES,
    build_process,
)
from solcascade.production.fission_safety import BUILT_EFFICIENCY
from solcascade.production.fission_tails import TAILS_PERCENT
from solcascade.production.plants import FAST_ID, THERMAL_ID, build_plants
from solcascade.production.refine import catalog_record, project_root

_COUNT = 116
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4
_RATIONAL_GAP = 1.0e-14
_NATURAL_T = 211.0
_CONCENTRATE_T = 249.0
_FEED_HEX_T = 312.0
_ENRICHED_HEX_T = 35.9
_FABRICATED_T = 27.6
_CORE_T = 73.0
_EIGHTEEN_MONTH_T = 36.5
_TRANSURANIC_T = 0.280
_COMPARE_NATURAL_T = 163.0
_COMPARE_PRODUCT_T = 16.85
_COMPARE_BURNUP_GWD = 65.0
_FEED_MULTIPLE = Fraction(2110, 243)
_PER_TWH = Fraction(5275, 219)
_CONCENTRATE_SHARE = Fraction(211, 249)
_HEX_SHARE = Fraction(211, 312)
_ENRICHED_HEX_SHARE = Fraction(243, 359)
_OXIDE_SHARE = Fraction(81, 92)
_CORE_FRACTION = Fraction(243, 730)
_TRANSURANIC_SHARE = Fraction(7, 690)
_COMPARE_MULTIPLE = Fraction(3260, 337)
_FAST_CEILING = Fraction(1, 2)


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
        "feed multiple",
        _FEED_MULTIPLE,
        True,
        None,
        "pi^2-G^-2",
        (
            "(pi^2-G^-2)*(1+alpha*phi^-5)",
            "((pi^2-G^-2)*(1+alpha*phi^-5))*(1-alpha*pi^-5)",
            "(((pi^2-G^-2)*(1+alpha*phi^-5))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
        ),
    ),
    _Spec(
        "natural uranium per terawatt-hour",
        _PER_TWH,
        True,
        None,
        "pi^3-phi^4",
        (
            "(pi^3-phi^4)*(1-alpha*phi^-2)",
            "((pi^3-phi^4)*(1-alpha*phi^-2))*(1+alpha*pi^-4)",
        ),
    ),
    _Spec(
        "concentrate share",
        _CONCENTRATE_SHARE,
        True,
        108,
        "e^-5+G^2",
        (
            "(e^-5+G^2)*(1+alpha*phi^-3)",
            "((e^-5+G^2)*(1+alpha*phi^-3))*(1+alpha*pi^-3)",
        ),
    ),
    _Spec(
        "feed hexafluoride share",
        _HEX_SHARE,
        True,
        111,
        "pi^-3+G^5",
        (
            "(pi^-3+G^5)*(1-alpha*phi^-4)",
            "((pi^-3+G^5)*(1-alpha*phi^-4))*(1-alpha*pi^-6)",
        ),
    ),
    _Spec(
        "enriched hexafluoride share",
        _ENRICHED_HEX_SHARE,
        True,
        111,
        "pi^-3+G^5",
        (
            "(pi^-3+G^5)*(1-alpha*pi^-3)",
            "((pi^-3+G^5)*(1-alpha*pi^-3))*(1+alpha*e^-5)",
        ),
    ),
    _Spec(
        "oxide share",
        _OXIDE_SHARE,
        True,
        107,
        "sqrt(e)-G^3",
        (
            "(sqrt(e)-G^3)*(1+alpha*pi^-3)",
            "((sqrt(e)-G^3)*(1+alpha*pi^-3))*(1-alpha*pi^-6)",
        ),
    ),
    _Spec(
        "core fraction",
        _CORE_FRACTION,
        True,
        113,
        "gamma*gamma",
        (
            "(gamma*gamma)*(1-alpha*e^-2)",
            "((gamma*gamma)*(1-alpha*e^-2))*(1+alpha*pi^-4)",
            "(((gamma*gamma)*(1-alpha*e^-2))*(1+alpha*pi^-4))*(1+alpha*pi^-6)",
        ),
    ),
    _Spec(
        "transuranic share",
        _TRANSURANIC_SHARE,
        False,
        _COUNT,
        "pi*pi^-5",
        (
            "(pi*pi^-5)*(1-alpha*phi)",
            "((pi*pi^-5)*(1-alpha*phi))*(1+alpha*e^-6)",
        ),
    ),
)


@dataclass(frozen=True)
class FeedStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class FeedFigure:
    label: str
    anchor: NamedGap
    steps: tuple[FeedStep, ...]
    target: float
    rational: Fraction
    positive: bool
    unit_rows: int | None


@dataclass(frozen=True)
class FeedReport:
    figures: tuple[FeedFigure, ...]
    natural_t: float
    reload_tonnes: float
    concentrate_t: float
    feed_hex_t: float
    enriched_hex_t: float
    fabricated_t: float
    core_t: float
    eighteen_month_t: float
    transuranic_t: float
    nameplate_twh: float
    compare_natural_t: float
    compare_product_t: float
    compare_burnup_gwd: float
    compare_multiple: float
    fast_ceiling: float
    enriched_percent: float
    tails_percent: float
    discharge_gwd: float
    higher_burnup_gwd: float
    fabrication_cap_percent: float
    built_efficiency: float
    annual_swu: float
    source: str
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_feed.md"


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


def _walk(spec: _Spec, alpha: float, helium: float, five, base) -> FeedFigure:
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
    steps: list[FeedStep] = []
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
            FeedStep(
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
    return FeedFigure(
        spec.label,
        anchor_named,
        tuple(steps),
        target,
        spec.rational,
        spec.positive,
        spec.unit_rows,
    )


def build_feed() -> FeedReport:
    plants = build_plants()
    process = build_process()
    module, digest = load_engine()
    if digest[:16] != plants.pin_prefix or process.pin_prefix != plants.pin_prefix:
        raise RuntimeError("the feed pin does not match the plant pin")
    record = catalog_record(THERMAL_ID)
    if "k_eff is held at 1" not in str(record["governing_equation"]):
        raise RuntimeError("the governing equation does not hold k_eff at 1")
    if "Sodium" not in str(catalog_record(FAST_ID).get("engineering_barrier") or ""):
        raise RuntimeError("the fast barrier does not name sodium")
    if RELOAD_TONNES != 24.3 or ENRICHED_PERCENT != 4.5 or TAILS_PERCENT != 0.22:
        raise RuntimeError("the operating column left 24.3 t, 4.5 percent, and 0.22 percent tails")
    if ANNUAL_SWU != 182_000.0 or process.discharge_gwd != 45.0:
        raise RuntimeError("the operating column left 182000 SWU and 45 GWd/t")
    if FABRICATION_CAP_PERCENT != 5.0 or HIGHER_BURNUP_GWD != 60.0:
        raise RuntimeError("the cap left 5 percent or the stack left 60 GWd/t")
    if not _COMPARE_BURNUP_GWD > HIGHER_BURNUP_GWD > process.discharge_gwd:
        raise RuntimeError("the comparison burn-up left its place above the stack")
    if BUILT_EFFICIENCY != 0.36 or _gap(process.efficiency, 0.33) > 1.0e-12:
        raise RuntimeError("the efficiencies left 0.36 and 0.33")
    if plants.fast.ceiling != 0.5:
        raise RuntimeError("the fast ceiling left 0.5")
    nameplate_twh = NAMEPLATE_MWE * FULL_POWER_HOURS / 1_000_000.0
    if Fraction(_NATURAL_T) / Fraction(RELOAD_TONNES).limit_denominator(10) != _FEED_MULTIPLE:
        raise RuntimeError("the feed multiple left 2110/243")
    if Fraction(211) * 1_000_000 / (Fraction(int(NAMEPLATE_MWE)) * Fraction(int(FULL_POWER_HOURS))) != _PER_TWH:
        raise RuntimeError("the natural uranium per terawatt-hour left 5275/219")
    if Fraction(36.5) != Fraction(73, 2):
        raise RuntimeError("the 18-month reload left half the core")
    if _EIGHTEEN_MONTH_T / _CORE_T != float(_FAST_CEILING):
        raise RuntimeError("the 18-month reload left the fast ceiling")
    if _COMPARE_NATURAL_T / _COMPARE_PRODUCT_T != float(_COMPARE_MULTIPLE) and _gap(
        _COMPARE_NATURAL_T / _COMPARE_PRODUCT_T, float(_COMPARE_MULTIPLE)
    ) > _RATIONAL_GAP:
        raise RuntimeError("the comparison multiple left 3260/337")
    feed_value = _NATURAL_T / RELOAD_TONNES
    per_twh = _NATURAL_T / nameplate_twh
    concentrate = _NATURAL_T / _CONCENTRATE_T
    hex_share = _NATURAL_T / _FEED_HEX_T
    enriched_hex = RELOAD_TONNES / _ENRICHED_HEX_T
    oxide = RELOAD_TONNES / _FABRICATED_T
    core_fraction = RELOAD_TONNES / _CORE_T
    transuranic = _TRANSURANIC_T / _FABRICATED_T
    compare_multiple = _COMPARE_NATURAL_T / _COMPARE_PRODUCT_T
    for label, value, rational in (
        ("feed multiple", feed_value, _FEED_MULTIPLE),
        ("natural uranium per terawatt-hour", per_twh, _PER_TWH),
        ("concentrate share", concentrate, _CONCENTRATE_SHARE),
        ("feed hexafluoride share", hex_share, _HEX_SHARE),
        ("enriched hexafluoride share", enriched_hex, _ENRICHED_HEX_SHARE),
        ("oxide share", oxide, _OXIDE_SHARE),
        ("core fraction", core_fraction, _CORE_FRACTION),
        ("transuranic share", transuranic, _TRANSURANIC_SHARE),
        ("comparison multiple", compare_multiple, _COMPARE_MULTIPLE),
    ):
        _require(value, rational, label)
    if not compare_multiple > feed_value:
        raise RuntimeError("the comparison multiple left its place above the operating feed")
    if not _COMPARE_NATURAL_T < _NATURAL_T or not _COMPARE_PRODUCT_T < RELOAD_TONNES:
        raise RuntimeError("the comparison masses left their place under the operating column")
    if not core_fraction < float(Fraction(1, 3)):
        raise RuntimeError("the annual reload left its place under one third of the core")
    if abs((_FABRICATED_T - 1.0 - _TRANSURANIC_T) - 26.0) < 0.05:
        raise RuntimeError("the printed 26 t uranium oxide closed the used-fuel row")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    five = _five_pieces(module)
    base = _base_pieces(module)
    figures = tuple(_walk(spec, alpha, helium, five, base) for spec in _SPECS)
    if _FEED_MULTIPLE * _PER_TWH != Fraction(211) * Fraction(211) * 1_000_000 / (
        Fraction(243, 10) * Fraction(1000) * Fraction(8760)
    ):
        raise RuntimeError("the feed product left the printed masses")
    return FeedReport(
        figures=figures,
        natural_t=_NATURAL_T,
        reload_tonnes=RELOAD_TONNES,
        concentrate_t=_CONCENTRATE_T,
        feed_hex_t=_FEED_HEX_T,
        enriched_hex_t=_ENRICHED_HEX_T,
        fabricated_t=_FABRICATED_T,
        core_t=_CORE_T,
        eighteen_month_t=_EIGHTEEN_MONTH_T,
        transuranic_t=_TRANSURANIC_T,
        nameplate_twh=nameplate_twh,
        compare_natural_t=_COMPARE_NATURAL_T,
        compare_product_t=_COMPARE_PRODUCT_T,
        compare_burnup_gwd=_COMPARE_BURNUP_GWD,
        compare_multiple=compare_multiple,
        fast_ceiling=plants.fast.ceiling,
        enriched_percent=ENRICHED_PERCENT,
        tails_percent=TAILS_PERCENT,
        discharge_gwd=process.discharge_gwd,
        higher_burnup_gwd=HIGHER_BURNUP_GWD,
        fabrication_cap_percent=FABRICATION_CAP_PERCENT,
        built_efficiency=BUILT_EFFICIENCY,
        annual_swu=ANNUAL_SWU,
        source=MATERIAL_BALANCE_SOURCE,
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=plants.nuclear_scalar,
        thermodynamics_scalar=plants.thermodynamics_scalar,
        gamma_fsot=plants.gamma_fsot,
        valve_beta=plants.valve_beta,
        pin_prefix=plants.pin_prefix,
        digest_prefix=plants.digest_prefix,
    )


def _step_lines(step: FeedStep, closing: bool) -> list[str]:
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


def _figure_lines(figure: FeedFigure, lead: str, role: str) -> list[str]:
    lines = [lead]
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
    lines.append("The closing value names that quotient.")
    lines.append(role)
    return lines


def _by_label(report: FeedReport, label: str) -> FeedFigure:
    for figure in report.figures:
        if figure.label == label:
            return figure
    raise RuntimeError(f"the {label} left the feed report")


def render() -> str:
    report = build_feed()
    feed = _by_label(report, "feed multiple")
    per_twh = _by_label(report, "natural uranium per terawatt-hour")
    concentrate = _by_label(report, "concentrate share")
    hex_share = _by_label(report, "feed hexafluoride share")
    enriched_hex = _by_label(report, "enriched hexafluoride share")
    oxide = _by_label(report, "oxide share")
    core = _by_label(report, "core fraction")
    transuranic = _by_label(report, "transuranic share")
    lines = [
        "# Fission feed",
        "",
        f"The material balance is {report.source}.",
        "The operating column is 211 t of natural uranium for 24.3 t of enriched uranium at 4.5 percent, 0.22 percent tails, 45 GWd/t, and 33 percent.",
        "FSOT names the operating quotients on that column.",
        "The chain stays at k_eff of 1.",
        "",
        "## Natural uranium",
        "",
    ]
    lines.extend(
        _figure_lines(
            feed,
            f"The feed multiple is `{_num(feed.target)}`.",
            "The operating feed multiple keeps the printed 211 t and the filed reload.",
        )
    )
    lines.extend(
        [
            "The same quotient is `2110/243`.",
            f"The filed reload stays `{_num(report.reload_tonnes)}` t.",
            f"The natural uranium stays `{_num(report.natural_t)}` t.",
            "The plants-being-built efficiency leaves these masses in place.",
            "",
        ]
    )
    lines.extend(
        _figure_lines(
            per_twh,
            f"Natural uranium over the nameplate year is `{_num(per_twh.target)}` t per terawatt-hour.",
            "That quotient keeps the printed 211 t and the nameplate year.",
        )
    )
    lines.extend(
        [
            "The same quotient is `5275/219`.",
            "The page rounds that quotient to 24.",
            "The rounded 24 stays the page's hence.",
            "",
            "## Element shares",
            "",
            f"The concentrate line is `{_num(report.concentrate_t)}` t of uranium oxide containing `{_num(report.natural_t)}` t of uranium.",
        ]
    )
    lines.extend(
        _figure_lines(
            concentrate,
            f"The uranium share of that concentrate is `{_num(concentrate.target)}`.",
            "That share keeps the printed concentrate masses.",
        )
    )
    lines.extend(
        [
            "The same quotient is `211/249`.",
            f"The conversion line is `{_num(report.feed_hex_t)}` t of uranium hexafluoride containing the same `{_num(report.natural_t)}` t of uranium.",
        ]
    )
    lines.extend(
        _figure_lines(
            hex_share,
            f"The uranium share of that hexafluoride is `{_num(hex_share.target)}`.",
            "That share keeps the printed conversion masses.",
        )
    )
    lines.extend(
        [
            "The same quotient is `211/312`.",
            f"The enrichment line is `{_num(report.enriched_hex_t)}` t of enriched hexafluoride containing `{_num(report.reload_tonnes)}` t of uranium.",
        ]
    )
    lines.extend(
        _figure_lines(
            enriched_hex,
            f"The uranium share of that enriched hexafluoride is `{_num(enriched_hex.target)}`.",
            "That share keeps the printed enrichment masses.",
        )
    )
    lines.extend(
        [
            "The same quotient is `243/359`.",
            "The feed hexafluoride and the enriched hexafluoride stay different masses.",
            f"The fabrication line is `{_num(report.fabricated_t)}` t of uranium dioxide containing `{_num(report.reload_tonnes)}` t of uranium.",
        ]
    )
    lines.extend(
        _figure_lines(
            oxide,
            f"The uranium share of that dioxide is `{_num(oxide.target)}`.",
            "That share keeps the printed fabrication masses.",
        )
    )
    lines.extend(
        [
            "The same quotient is `81/92`.",
            "",
            "## Core load",
            "",
            f"The note prints a core load of `{_num(report.core_t)}` tU.",
        ]
    )
    lines.extend(
        _figure_lines(
            core,
            f"The annual reload over that core is `{_num(core.target)}`.",
            "That fraction keeps the filed reload and the printed core load.",
        )
    )
    lines.extend(
        [
            "The same quotient is `243/730`.",
            "The annual reload over the core stays under the burnup gain.",
            "The prose one-third stays approximate.",
            f"The 18-month reload is `{_num(report.eighteen_month_t)}` tU.",
            "The 18-month reload is half the core load.",
            "That half keeps the fast ceiling.",
            "",
            "## Used fuel",
            "",
            f"The used-fuel line is `{_num(report.fabricated_t)}` t containing `{_num(report.transuranic_t * 1000.0)}` kg of transuranics.",
        ]
    )
    lines.extend(
        _figure_lines(
            transuranic,
            f"The transuranic share of that used fuel is `{_num(transuranic.target)}`.",
            "That share keeps the printed used-fuel masses.",
        )
    )
    lines.extend(
        [
            "The same quotient is `7/690`.",
            "The used-fuel uranium oxide is printed as 26 t and does not close the fabricated mass.",
            "The used-fuel fission-product mass stays the printed 1 tonne.",
            "The fast row's filed barrier is the sodium system.",
            "",
            "## Comparison column",
            "",
            f"The comparison column prints `{_num(report.compare_natural_t)}` t of natural uranium for `{_num(report.compare_product_t)}` t of enriched uranium.",
            f"Its feed multiple is `{_num(report.compare_multiple)}`.",
            "The same quotient is `3260/337`.",
            "That multiple stays above the operating feed multiple.",
            "The comparison product stays under the filed reload.",
            "The comparison natural uranium stays under the operating natural uranium.",
            "That column sits on the fabrication cap.",
            f"The `{_num(report.compare_burnup_gwd)}` GWd/t burn-up stays above the operating discharge.",
            "The 60 GWd/t stack stays unrun.",
            "This comparison stays unrun.",
            "The enrichment note's about 140000 SWU stays approximate.",
            "The mining line stays a range.",
            f"The enriched percent stays `{_num(report.enriched_percent)}`.",
            f"The tails percent stays `{_num(report.tails_percent)}`.",
            f"The annual separative work stays `{_num(report.annual_swu)}` SWU.",
            "The demonstrated plant stays at its filed efficiency.",
            f"The plants-being-built efficiency stays `{_num(report.built_efficiency)}`.",
            "The fabrication cap stays 5 percent.",
            "The catalog ceiling stays unrun.",
            "This score adds no joules.",
            "The published stamp records these closes.",
            "The published stamp records the comparison multiple as a rational of the printed masses.",
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
