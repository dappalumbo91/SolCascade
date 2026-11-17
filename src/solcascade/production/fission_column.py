"""Record the operating-column quotients that keep figures already named.

The separative work over the nameplate year is the named separative work
per kilogram times the enriched kilograms in that terawatt-hour. The
product share is the reciprocal of the feed multiple. The page's other
exact burn-ups stay unrun. The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from solcascade.production.fission_ceiling import build_ceiling
from solcascade.production.fission_process import (
    ANNUAL_SWU,
    ATOMIC_PERCENT,
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

_RATIONAL_GAP = 1.0e-14
_JOULE_PER_TWH = 1_000_000_000_000 * 3600
_NATURAL_T = 211.0
_CORE_T = 73.0
_EIGHTEEN_MONTH_T = 36.5
_COMPARE_NATURAL_T = 163.0
_COMPARE_PRODUCT_T = 16.85
_COMPARE_BURNUP_GWD = 65.0
_CAP_BURNUP_GWD = 55.0
_SIGHTED_BURNUP_GWD = 70.0
_SIGHTED_PERCENT = 6.0
_OPERATING_GWD = 45.0
_SWU_PER_TWH = Fraction(4550000, 219)
_PRODUCT_SHARE = Fraction(243, 2110)
_TAILS_STREAM = Fraction(1867, 2110)
_SWU_PER_NATURAL_T = Fraction(182000, 211)
_ENRICHED_PER_TWH = Fraction(405, 146)
_FISSILE_PART = Fraction(26001, 5275000)
_COMPARE_PRODUCT = Fraction(337, 3260)
_COMPARE_MULTIPLE = Fraction(3260, 337)
_CAP_STEP = Fraction(11, 9)
_SIGHTED_STEP = Fraction(14, 9)
_BUILT_STEP = Fraction(12, 11)
_FAST_CEILING = Fraction(1, 2)
_FEED_MULTIPLE = Fraction(2110, 243)
_PER_TWH = Fraction(5275, 219)


@dataclass(frozen=True)
class ColumnReport:
    swu_per_twh: float
    product_share: float
    tails_stream: float
    swu_per_natural_t: float
    enriched_per_twh: float
    fissile_part: float
    cap_step: float
    sighted_step: float
    built_step: float
    compare_product: float
    compare_multiple: float
    nameplate_twh: float
    annual_swu: float
    reload_tonnes: float
    natural_t: float
    discharge_gwd: float
    cap_burnup_gwd: float
    higher_burnup_gwd: float
    compare_burnup_gwd: float
    sighted_burnup_gwd: float
    enriched_percent: float
    tails_percent: float
    fabrication_cap_percent: float
    atomic_percent: float
    sighted_percent: float
    fast_ceiling: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str
    source: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_column.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _gap(derived: float, witness: float) -> float:
    if witness == 0.0:
        raise RuntimeError("the witness is zero")
    return abs(derived - witness) / abs(witness)


def _require(value: float, rational: Fraction, label: str) -> None:
    if _gap(value, float(rational)) > _RATIONAL_GAP:
        raise RuntimeError(f"the {label} left {rational.numerator}/{rational.denominator}")


def build_column() -> ColumnReport:
    plants = build_plants()
    process = build_process()
    ceiling = build_ceiling()
    if process.pin_prefix != plants.pin_prefix or ceiling.pin_prefix != plants.pin_prefix:
        raise RuntimeError("the column pin does not match the plant pin")
    record = catalog_record(THERMAL_ID)
    if "k_eff is held at 1" not in str(record["governing_equation"]):
        raise RuntimeError("the governing equation does not hold k_eff at 1")
    if "Sodium" not in str(catalog_record(FAST_ID).get("engineering_barrier") or ""):
        raise RuntimeError("the fast barrier does not name sodium")
    if ANNUAL_SWU != 182_000.0 or RELOAD_TONNES != 24.3:
        raise RuntimeError("the operating column left 182000 SWU and 24.3 t")
    if ENRICHED_PERCENT != 4.5 or TAILS_PERCENT != 0.22:
        raise RuntimeError("the operating column left 4.5 percent and 0.22 percent tails")
    if process.discharge_gwd != _OPERATING_GWD or HIGHER_BURNUP_GWD != 60.0:
        raise RuntimeError("the burnup pair left 45 GWd/t and 60 GWd/t")
    if FABRICATION_CAP_PERCENT != 5.0 or not ATOMIC_PERCENT > FABRICATION_CAP_PERCENT:
        raise RuntimeError("the atomic percent left its place above the fabrication cap")
    if not _SIGHTED_PERCENT > FABRICATION_CAP_PERCENT:
        raise RuntimeError("the sighted enrichment left its place above the fabrication cap")
    if BUILT_EFFICIENCY != 0.36 or _gap(process.efficiency, 0.33) > 1.0e-12:
        raise RuntimeError("the demonstrated plant left 0.33 and the built plant left 0.36")
    if plants.fast.ceiling != 0.5:
        raise RuntimeError("the fast ceiling left 0.5")
    if Fraction(36.5) != Fraction(73, 2):
        raise RuntimeError("the 18-month reload left half the core")
    if _EIGHTEEN_MONTH_T / _CORE_T != float(_FAST_CEILING):
        raise RuntimeError("the 18-month reload left the fast ceiling")
    if not (
        _OPERATING_GWD
        < _CAP_BURNUP_GWD
        < HIGHER_BURNUP_GWD
        < _COMPARE_BURNUP_GWD
        < _SIGHTED_BURNUP_GWD
    ):
        raise RuntimeError("the page burn-ups left 45, 55, 60, 65, 70")
    nameplate_twh = process.nameplate_j / _JOULE_PER_TWH
    if FULL_POWER_HOURS != 8760.0:
        raise RuntimeError("the nameplate hours left 8760")
    _require(nameplate_twh, Fraction(8760, 1000), "nameplate year")
    if NAMEPLATE_MWE != 1000.0:
        raise RuntimeError("the nameplate left 1000 MWe")
    feed_multiple = _NATURAL_T / RELOAD_TONNES
    per_twh = _NATURAL_T / nameplate_twh
    _require(feed_multiple, _FEED_MULTIPLE, "feed multiple")
    _require(per_twh, _PER_TWH, "natural uranium per terawatt-hour")
    if _PER_TWH / _FEED_MULTIPLE != _ENRICHED_PER_TWH:
        raise RuntimeError("the enriched rate left the per-terawatt-hour quotient over the feed multiple")
    enriched_per_twh = per_twh / feed_multiple
    _require(enriched_per_twh, _ENRICHED_PER_TWH, "enriched uranium per terawatt-hour")
    swu_per_kg = process.swu_per_kg_reload
    _require(swu_per_kg, Fraction(1820, 243), "separative work per kilogram")
    swu_per_twh = swu_per_kg * RELOAD_TONNES * 1000.0 / nameplate_twh
    if _gap(swu_per_twh, ANNUAL_SWU / nameplate_twh) > 0.0:
        raise RuntimeError("the separative year left the annual separative work over the nameplate year")
    _require(swu_per_twh, _SWU_PER_TWH, "separative work per terawatt-hour")
    if Fraction(ANNUAL_SWU) / Fraction(8760, 1000) != _SWU_PER_TWH:
        raise RuntimeError("the separative year left 4550000/219")
    product_share = RELOAD_TONNES / _NATURAL_T
    _require(product_share, _PRODUCT_SHARE, "product share")
    if _PRODUCT_SHARE * _FEED_MULTIPLE != 1:
        raise RuntimeError("the product share left the reciprocal of the feed multiple")
    tails_stream = 1.0 - product_share
    _require(tails_stream, _TAILS_STREAM, "tails stream")
    if _TAILS_STREAM != 1 - _PRODUCT_SHARE:
        raise RuntimeError("the tails stream left one minus the product share")
    swu_per_natural_t = ANNUAL_SWU / _NATURAL_T
    _require(swu_per_natural_t, _SWU_PER_NATURAL_T, "separative work per natural tonne")
    if Fraction(1820, 243) * 1000 / _FEED_MULTIPLE != _SWU_PER_NATURAL_T:
        raise RuntimeError("the natural-tonne separative work left the reload score and the feed multiple")
    fissile_part = (ENRICHED_PERCENT - TAILS_PERCENT) / 100.0 / feed_multiple
    _require(fissile_part, _FISSILE_PART, "fissile part")
    enriched_fraction = Fraction(45, 1000)
    tails_fraction = Fraction(22, 10000)
    if (enriched_fraction - tails_fraction) / _FEED_MULTIPLE != _FISSILE_PART:
        raise RuntimeError("the fissile part left the enriched fraction, the tails fraction, and the feed multiple")
    if _gap(0.0049, fissile_part) <= 1.0e-3:
        raise RuntimeError("the prose 0.49 percent met the fissile part inside the review bar")
    cap_step = _CAP_BURNUP_GWD / process.discharge_gwd
    sighted_step = _SIGHTED_BURNUP_GWD / process.discharge_gwd
    _require(cap_step, _CAP_STEP, "55 GWd/t step")
    _require(sighted_step, _SIGHTED_STEP, "70 GWd/t step")
    if Fraction(HIGHER_BURNUP_GWD) / Fraction(_CAP_BURNUP_GWD) != _BUILT_STEP:
        raise RuntimeError("sixty over 55 left the named built step")
    built_step = BUILT_EFFICIENCY / process.efficiency
    _require(built_step, _BUILT_STEP, "built step")
    compare_multiple = _COMPARE_NATURAL_T / _COMPARE_PRODUCT_T
    compare_product = _COMPARE_PRODUCT_T / _COMPARE_NATURAL_T
    _require(compare_multiple, _COMPARE_MULTIPLE, "comparison multiple")
    _require(compare_product, _COMPARE_PRODUCT, "comparison product share")
    if _COMPARE_PRODUCT * _COMPARE_MULTIPLE != 1:
        raise RuntimeError("the comparison product share left the reciprocal of the unrun multiple")
    return ColumnReport(
        swu_per_twh=swu_per_twh,
        product_share=product_share,
        tails_stream=tails_stream,
        swu_per_natural_t=swu_per_natural_t,
        enriched_per_twh=enriched_per_twh,
        fissile_part=fissile_part,
        cap_step=cap_step,
        sighted_step=sighted_step,
        built_step=built_step,
        compare_product=compare_product,
        compare_multiple=compare_multiple,
        nameplate_twh=nameplate_twh,
        annual_swu=ANNUAL_SWU,
        reload_tonnes=RELOAD_TONNES,
        natural_t=_NATURAL_T,
        discharge_gwd=process.discharge_gwd,
        cap_burnup_gwd=_CAP_BURNUP_GWD,
        higher_burnup_gwd=HIGHER_BURNUP_GWD,
        compare_burnup_gwd=_COMPARE_BURNUP_GWD,
        sighted_burnup_gwd=_SIGHTED_BURNUP_GWD,
        enriched_percent=ENRICHED_PERCENT,
        tails_percent=TAILS_PERCENT,
        fabrication_cap_percent=FABRICATION_CAP_PERCENT,
        atomic_percent=ATOMIC_PERCENT,
        sighted_percent=_SIGHTED_PERCENT,
        fast_ceiling=plants.fast.ceiling,
        nuclear_scalar=process.nuclear_scalar,
        thermodynamics_scalar=process.thermodynamics_scalar,
        gamma_fsot=process.gamma_fsot,
        valve_beta=process.valve_beta,
        pin_prefix=process.pin_prefix,
        digest_prefix=process.digest_prefix,
        source=MATERIAL_BALANCE_SOURCE,
    )


def render() -> str:
    report = build_column()
    lines = [
        "# Fission column",
        "",
        f"The material balance is {report.source}.",
        "The operating column leaves no further free quotient.",
        "The chain stays at k_eff of 1.",
        "",
        "## Separative year",
        "",
        f"The nameplate year is `{_num(report.nameplate_twh)}` TWh.",
        f"The annual separative work stays `{_num(report.annual_swu)}` SWU.",
        f"The separative work over the nameplate year is `{_num(report.swu_per_twh)}` SWU per terawatt-hour.",
        "The same quotient is `4550000/219`.",
        "That quotient is the named separative work per kilogram times the enriched kilograms in that terawatt-hour.",
        "That quotient keeps the reload score and the feed score.",
        f"The enriched uranium in one nameplate terawatt-hour is `{_num(report.enriched_per_twh)}` t.",
        "The same quotient is `405/146`.",
        "That rate keeps the per-terawatt-hour quotient over the feed multiple.",
        f"The filed reload stays `{_num(report.reload_tonnes)}` t.",
        "",
        "## Product share",
        "",
        f"The natural uranium stays `{_num(report.natural_t)}` t.",
        f"The product share of the natural uranium is `{_num(report.product_share)}`.",
        "The same quotient is `243/2110`.",
        "The product share is the reciprocal of the feed multiple.",
        "The product share keeps the feed score.",
        f"The tails stream is `{_num(report.tails_stream)}` of the natural uranium.",
        "The same quotient is `1867/2110`.",
        "The tails stream is one minus the product share.",
        "The tails stream keeps the feed score.",
        "The page prints no tails mass.",
        f"The separative work per tonne of natural uranium is `{_num(report.swu_per_natural_t)}` SWU.",
        "The same quotient is `182000/211`.",
        "That quotient keeps the named separative work per kilogram and the feed multiple.",
        f"The fissile part of the natural uranium is `{_num(report.fissile_part)}`.",
        "The same quotient is `26001/5275000`.",
        "That quotient keeps the enriched fraction, the tails fraction, and the feed multiple.",
        "The prose 0.49 percent stays the page's round.",
        "The prose 0.394 percent stays outside the operating column.",
        "The prose 0.6 percent utilization stays approximate.",
        "",
        "## Page burn-ups",
        "",
        f"The operating discharge stays `{_num(report.discharge_gwd)}` GWd/t.",
        f"The page's 55 GWd/t step over that discharge is `{_num(report.cap_step)}`.",
        "The same quotient is `11/9`.",
        "That burn-up sits on the fabrication cap and stays unrun.",
        f"Sixty over that 55 GWd/t stays the named built step `{_num(report.built_step)}`.",
        "The same quotient is `12/11`.",
        "That quotient keeps the built score.",
        "The 60 GWd/t stack stays unrun.",
        f"The comparison burn-up stays `{_num(report.compare_burnup_gwd)}` GWd/t.",
        f"The comparison product share is `{_num(report.compare_product)}`.",
        "The same quotient is `337/3260`.",
        "The comparison product share is the reciprocal of the unrun multiple.",
        "This comparison stays unrun.",
        f"The page's 70 GWd/t step over the operating discharge is `{_num(report.sighted_step)}`.",
        "The same quotient is `14/9`.",
        f"The page states `{_num(report.sighted_percent)}` percent enrichment for that burn-up.",
        "That enrichment stays above the fabrication cap.",
        "That burn-up stays unrun.",
        f"The atomic percent stays `{_num(report.atomic_percent)}`.",
        "The atomic percent stays above the fabrication cap.",
        "The 18-month reload is half the core load.",
        "That half keeps the fast ceiling.",
        "The enrichment note's about 140000 SWU stays approximate.",
        "The mining line stays a range.",
        "The 90 percent load factor stays the page's round.",
        "The Red Book tails assays are 0.30 percent and 0.25 percent.",
        "Those assays stay off this column.",
        "The CANDU row stays a different plant.",
        f"The enriched percent stays `{_num(report.enriched_percent)}`.",
        f"The tails percent stays `{_num(report.tails_percent)}`.",
        "The fabrication cap stays 5 percent.",
        "The demonstrated plant stays at its filed efficiency.",
        "The plants-being-built efficiency leaves the operating masses in place.",
        "The catalog ceiling stays unrun.",
        "The fast row's filed barrier is the sodium system.",
        "This score adds no joules.",
        "The published stamp records these quotients as rationals of the named scores.",
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
        "This run does not install a seed as a new enrichment or a new plant row.",
        "",
        "The checked-in schedule is not rewritten.",
        "",
    ]
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
