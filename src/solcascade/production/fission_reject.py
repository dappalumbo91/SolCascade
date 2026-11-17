"""Record the cooling-page heat gap between the two named efficiencies.

A plant at 33 percent dumps more heat per kilowatt-hour than a plant at
36 percent. The extra share keeps the named built step. The page's about
14 percent is that share's round. The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from solcascade.production.fission_ceiling import build_ceiling
from solcascade.production.fission_process import build_process
from solcascade.production.fission_safety import BUILT_EFFICIENCY, COOLING_SOURCE, heat_dump_mj_per_kwh
from solcascade.production.plants import FAST_ID, THERMAL_ID, build_plants
from solcascade.production.refine import catalog_record, project_root

_RATIONAL_GAP = 1.0e-14
_REVIEW_GAP = 1.0e-3
_DUMP_RATIO = Fraction(201, 176)
_EXTRA_HEAT = Fraction(25, 176)
_EXTRA_OVER_GAIN = Fraction(25, 16)
_PUBLISHED_GAIN = Fraction(1, 11)
_BUILT_STEP = Fraction(12, 11)
_DEMONSTRATED = Fraction(33, 100)
_BUILT = Fraction(9, 25)
_PAGE_ROUND = 0.14


@dataclass(frozen=True)
class RejectReport:
    demonstrated: float
    built: float
    carnot: float
    catalog: float
    built_step: float
    published_gain: float
    demonstrated_dump: float
    built_dump: float
    dump_ratio: float
    extra_heat: float
    extra_over_gain: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str
    source: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_reject.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _gap(derived: float, witness: float) -> float:
    if witness == 0.0:
        raise RuntimeError("the witness is zero")
    return abs(derived - witness) / abs(witness)


def _require(value: float, rational: Fraction, label: str) -> None:
    if _gap(value, float(rational)) > _RATIONAL_GAP:
        raise RuntimeError(f"the {label} left {rational.numerator}/{rational.denominator}")


def build_reject() -> RejectReport:
    plants = build_plants()
    process = build_process()
    ceiling = build_ceiling()
    if process.pin_prefix != plants.pin_prefix or ceiling.pin_prefix != plants.pin_prefix:
        raise RuntimeError("the reject pin does not match the plant pin")
    record = catalog_record(THERMAL_ID)
    if "k_eff is held at 1" not in str(record["governing_equation"]):
        raise RuntimeError("the governing equation does not hold k_eff at 1")
    if "Sodium" not in str(catalog_record(FAST_ID).get("engineering_barrier") or ""):
        raise RuntimeError("the fast barrier does not name sodium")
    if ceiling.hot_k != 560.0 or ceiling.cold_k != 310.0:
        raise RuntimeError("the ceiling temperatures left 560 K and 310 K")
    if BUILT_EFFICIENCY != 0.36:
        raise RuntimeError("the plants-being-built figure left 0.36")
    if _gap(process.efficiency, 0.33) > 1.0e-12:
        raise RuntimeError("the demonstrated efficiency left 0.33")
    if plants.thermal.ceiling != 0.45:
        raise RuntimeError("the catalog ceiling left 0.45")
    demonstrated = process.efficiency
    built = BUILT_EFFICIENCY
    carnot = ceiling.carnot
    if not demonstrated < built < carnot < plants.thermal.ceiling:
        raise RuntimeError("the plants-being-built figure left the open interval under the Carnot")
    built_step = built / demonstrated
    published_gain = built_step - 1.0
    _require(built_step, _BUILT_STEP, "built step")
    _require(published_gain, _PUBLISHED_GAIN, "published gain")
    if _BUILT_STEP * (1 - _DEMONSTRATED) / (1 - _BUILT) != _DUMP_RATIO:
        raise RuntimeError("the heat ratio left the named built step")
    if _DUMP_RATIO - 1 != _EXTRA_HEAT:
        raise RuntimeError("the extra heat left one less than the heat ratio")
    if _EXTRA_HEAT / _PUBLISHED_GAIN != _EXTRA_OVER_GAIN:
        raise RuntimeError("the extra heat left 25/16 of the published gain")
    demonstrated_dump = heat_dump_mj_per_kwh(demonstrated)
    built_dump = heat_dump_mj_per_kwh(built)
    if not built_dump < demonstrated_dump:
        raise RuntimeError("the heat dump did not fall as the efficiency rose")
    dump_ratio = demonstrated_dump / built_dump
    extra_heat = dump_ratio - 1.0
    extra_over_gain = extra_heat / published_gain
    _require(dump_ratio, _DUMP_RATIO, "heat ratio")
    _require(extra_heat, _EXTRA_HEAT, "extra heat")
    _require(extra_over_gain, _EXTRA_OVER_GAIN, "extra heat over the published gain")
    if _gap(_PAGE_ROUND, extra_heat) <= _REVIEW_GAP:
        raise RuntimeError("the prose 14 percent met the extra heat inside the review bar")
    return RejectReport(
        demonstrated=demonstrated,
        built=built,
        carnot=carnot,
        catalog=plants.thermal.ceiling,
        built_step=built_step,
        published_gain=published_gain,
        demonstrated_dump=demonstrated_dump,
        built_dump=built_dump,
        dump_ratio=dump_ratio,
        extra_heat=extra_heat,
        extra_over_gain=extra_over_gain,
        nuclear_scalar=process.nuclear_scalar,
        thermodynamics_scalar=process.thermodynamics_scalar,
        gamma_fsot=process.gamma_fsot,
        valve_beta=process.valve_beta,
        pin_prefix=process.pin_prefix,
        digest_prefix=process.digest_prefix,
        source=COOLING_SOURCE,
    )


def render() -> str:
    report = build_reject()
    lines = [
        "# Fission reject",
        "",
        f"The cooling page is {report.source}.",
        "The cooling a steam plant requires is set by its thermal efficiency.",
        "The chain stays at k_eff of 1.",
        "",
        "## Extra heat",
        "",
        f"The demonstrated efficiency is `{_num(report.demonstrated)}`.",
        f"The plants-being-built efficiency is `{_num(report.built)}`.",
        f"The demonstrated dump stays `{_num(report.demonstrated_dump)}` MJ per kilowatt-hour.",
        f"The built dump stays `{_num(report.built_dump)}` MJ per kilowatt-hour.",
        f"The demonstrated plant dumps `{_num(report.dump_ratio)}` times that built dump.",
        "The same quotient is `201/176`.",
        "That ratio keeps the named built step.",
        f"The extra heat is `{_num(report.extra_heat)}`.",
        "The same quotient is `25/176`.",
        "The extra heat is one less than the heat ratio.",
        "The extra heat keeps the named built step.",
        "The page's about 14 percent stays that round.",
        f"The named built step stays `{_num(report.built_step)}`.",
        f"The published gain stays `{_num(report.published_gain)}`.",
        f"The extra heat is `{_num(report.extra_over_gain)}` times the published gain.",
        "The same quotient is `25/16`.",
        "That multiple keeps the named built step and the published gain.",
        "The extra heat stays above the published gain.",
        f"The temperature Carnot is `{_num(report.carnot)}`.",
        "The temperature Carnot stays the bound.",
        f"The catalog ceiling is `{_num(report.catalog)}`.",
        "The catalog ceiling stays unrun.",
        "The prose about two-thirds stays approximate.",
        "The tower penalty stays the page's 2 to 5 percent range.",
        "The fan share stays the page's about 1 percent.",
        "The page's evaporative figure stays an upper bound.",
        "The 60 GWd/t stack stays unrun.",
        "The fabrication cap stays 5 percent.",
        "The demonstrated plant stays at its filed efficiency.",
        "The fast row's filed barrier is the sodium system.",
        "This score adds no joules.",
        "The published stamp records these quotients as rationals of the named steps.",
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
        "This run does not install a seed as a new efficiency or a new plant row.",
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
