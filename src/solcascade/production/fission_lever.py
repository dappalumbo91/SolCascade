"""Rank the electricity gains the steam share leaves on the filed figures.

The published gain is the plants-being-built step over the demonstrated plant.
The remaining gain runs from that page figure up to the temperature Carnot.
The burnup gain is the unrun 60 GWd/t step. Each ratio keeps the named steps.
The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from solcascade.production.fission_ceiling import build_ceiling
from solcascade.production.fission_process import (
    ATOMIC_PERCENT,
    FABRICATION_CAP_PERCENT,
    HIGHER_BURNUP_GWD,
    RELOAD_TONNES,
    build_process,
)
from solcascade.production.fission_safety import BUILT_EFFICIENCY
from solcascade.production.plants import FAST_ID, THERMAL_ID, build_plants
from solcascade.production.refine import catalog_record, project_root

_RATIONAL_GAP = 1.0e-14
_PUBLISHED_GAIN = Fraction(1, 11)
_REMAINING_GAIN = Fraction(121, 504)
_GAIN_RATIO = Fraction(1331, 504)
_BURNUP_STEP = Fraction(4, 3)
_BURNUP_GAIN = Fraction(1, 3)
_BURNUP_RATIO = Fraction(11, 3)


@dataclass(frozen=True)
class LeverReport:
    demonstrated: float
    built: float
    carnot: float
    catalog: float
    built_step: float
    built_carnot_step: float
    published_gain: float
    remaining_gain: float
    gain_ratio: float
    burnup_step: float
    burnup_gain: float
    burnup_ratio: float
    discharge_gwd: float
    higher_burnup_gwd: float
    atomic_percent: float
    fabrication_cap_percent: float
    companion: float
    built_companion: float
    reload_tonnes: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_lever.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _gap(derived: float, witness: float) -> float:
    if witness == 0.0:
        raise RuntimeError("the witness is zero")
    return abs(derived - witness) / abs(witness)


def _require(value: float, rational: Fraction, label: str) -> None:
    if _gap(value, float(rational)) > _RATIONAL_GAP:
        raise RuntimeError(f"the {label} left {rational.numerator}/{rational.denominator}")


def build_lever() -> LeverReport:
    plants = build_plants()
    ceiling = build_ceiling()
    process = build_process()
    if process.pin_prefix != plants.pin_prefix or ceiling.pin_prefix != plants.pin_prefix:
        raise RuntimeError("the lever pin does not match the plant pin")
    record = catalog_record(THERMAL_ID)
    if "k_eff is held at 1" not in str(record["governing_equation"]):
        raise RuntimeError("the governing equation does not hold k_eff at 1")
    if "Sodium" not in str(catalog_record(FAST_ID).get("engineering_barrier") or ""):
        raise RuntimeError("the fast barrier does not name sodium")
    if ceiling.hot_k != 560.0 or ceiling.cold_k != 310.0:
        raise RuntimeError("the ceiling temperatures left 560 K and 310 K")
    if _gap(ceiling.carnot, 1.0 - ceiling.cold_k / ceiling.hot_k) > 0.0:
        raise RuntimeError("the temperature Carnot left one minus the rejected fraction")
    if _gap(process.carnot, ceiling.carnot) > _RATIONAL_GAP:
        raise RuntimeError("the process Carnot left the ceiling")
    if BUILT_EFFICIENCY != 0.36:
        raise RuntimeError("the plants-being-built figure left 0.36")
    if _gap(process.efficiency, 0.33) > 1.0e-12:
        raise RuntimeError("the demonstrated efficiency left 0.33")
    if process.discharge_gwd != 45.0 or HIGHER_BURNUP_GWD != 60.0:
        raise RuntimeError("the burnup pair left 45 GWd/t and 60 GWd/t")
    if FABRICATION_CAP_PERCENT != 5.0 or not ATOMIC_PERCENT > FABRICATION_CAP_PERCENT:
        raise RuntimeError("the atomic percent left its place above the fabrication cap")
    if plants.thermal.ceiling != 0.45:
        raise RuntimeError("the catalog ceiling left 0.45")
    demonstrated = process.efficiency
    built = BUILT_EFFICIENCY
    carnot = ceiling.carnot
    if not demonstrated < built < carnot < plants.thermal.ceiling:
        raise RuntimeError("the plants-being-built figure left the open interval under the Carnot")
    built_step = built / demonstrated
    built_carnot_step = carnot / built
    published_gain = built_step - 1.0
    remaining_gain = built_carnot_step - 1.0
    gain_ratio = remaining_gain / published_gain
    burnup_step = HIGHER_BURNUP_GWD / process.discharge_gwd
    burnup_gain = burnup_step - 1.0
    burnup_ratio = burnup_gain / published_gain
    if not published_gain < remaining_gain < burnup_gain:
        raise RuntimeError("the electricity gains left the scored order")
    if (Fraction(625, 504) - 1) / (Fraction(12, 11) - 1) != _GAIN_RATIO:
        raise RuntimeError("the gain ratio left 1331/504")
    if Fraction(625, 504) - 1 != _REMAINING_GAIN or Fraction(12, 11) - 1 != _PUBLISHED_GAIN:
        raise RuntimeError("the electricity gains left 121/504 and 1/11")
    if (Fraction(4, 3) - 1) / (Fraction(12, 11) - 1) != _BURNUP_RATIO:
        raise RuntimeError("the burnup ratio left 11/3")
    if Fraction(4, 3) - 1 != _BURNUP_GAIN:
        raise RuntimeError("the burnup gain left 1/3")
    _require(published_gain, _PUBLISHED_GAIN, "published gain")
    _require(remaining_gain, _REMAINING_GAIN, "remaining gain")
    _require(gain_ratio, _GAIN_RATIO, "gain ratio")
    _require(burnup_step, _BURNUP_STEP, "burnup step")
    _require(burnup_gain, _BURNUP_GAIN, "burnup gain")
    _require(burnup_ratio, _BURNUP_RATIO, "burnup ratio")
    companion = process.fuel_side_j / process.nameplate_j
    built_companion = companion * built_step
    if not 0.0 < companion < 1.0 < built_companion:
        raise RuntimeError("the built companion left its place across the nameplate year")
    if _gap(built_companion, companion * (built / demonstrated)) > 0.0:
        raise RuntimeError("the built companion left the filed companion times the built step")
    if RELOAD_TONNES != 24.3:
        raise RuntimeError("the filed reload left 24.3 t")
    return LeverReport(
        demonstrated=demonstrated,
        built=built,
        carnot=carnot,
        catalog=plants.thermal.ceiling,
        built_step=built_step,
        built_carnot_step=built_carnot_step,
        published_gain=published_gain,
        remaining_gain=remaining_gain,
        gain_ratio=gain_ratio,
        burnup_step=burnup_step,
        burnup_gain=burnup_gain,
        burnup_ratio=burnup_ratio,
        discharge_gwd=process.discharge_gwd,
        higher_burnup_gwd=HIGHER_BURNUP_GWD,
        atomic_percent=ATOMIC_PERCENT,
        fabrication_cap_percent=FABRICATION_CAP_PERCENT,
        companion=companion,
        built_companion=built_companion,
        reload_tonnes=RELOAD_TONNES,
        nuclear_scalar=process.nuclear_scalar,
        thermodynamics_scalar=process.thermodynamics_scalar,
        gamma_fsot=process.gamma_fsot,
        valve_beta=process.valve_beta,
        pin_prefix=process.pin_prefix,
        digest_prefix=process.digest_prefix,
    )


def render() -> str:
    report = build_lever()
    lines = [
        "# Fission lever",
        "",
        "The steam share leaves a larger electricity gain above the plants-being-built efficiency.",
        "The chain stays at k_eff of 1.",
        "",
        "## Electricity gains",
        "",
        f"The demonstrated efficiency is `{_num(report.demonstrated)}`.",
        f"The plants-being-built efficiency is `{_num(report.built)}`.",
        f"The temperature Carnot is `{_num(report.carnot)}`.",
        f"The catalog ceiling is `{_num(report.catalog)}`.",
        f"The named built step stays `{_num(report.built_step)}`.",
        f"The published gain is `{_num(report.published_gain)}`.",
        "The same quotient is `1/11`.",
        "The published gain is one less than the named built step.",
        "The published gain keeps the built score.",
        f"The named Carnot step on the built score stays `{_num(report.built_carnot_step)}`.",
        f"The remaining gain is `{_num(report.remaining_gain)}`.",
        "The same quotient is `121/504`.",
        "The remaining gain is one less than the named Carnot step on the built score.",
        "The remaining gain keeps the built score.",
        f"The remaining gain over the published gain is `{_num(report.gain_ratio)}`.",
        "The same quotient is `1331/504`.",
        "That ratio keeps the two named steps.",
        "The published gain stays under the remaining gain.",
        "The temperature Carnot stays the bound.",
        "",
        "## Burnup gain",
        "",
        f"The 60 GWd/t step over 45 GWd/t is `{_num(report.burnup_step)}`.",
        "The same quotient is `4/3`.",
        "The burnup step stays the stack score.",
        f"The burnup gain is `{_num(report.burnup_gain)}`.",
        "The same quotient is `1/3`.",
        "The burnup gain is one less than the burnup step.",
        "The burnup gain keeps the stack score.",
        f"The burnup gain over the published gain is `{_num(report.burnup_ratio)}`.",
        "The same quotient is `11/3`.",
        "That ratio keeps the burnup step and the named built step.",
        "The remaining gain stays under the burnup gain.",
        f"The atomic percent stays `{_num(report.atomic_percent)}`.",
        "The atomic percent stays above the fabrication cap.",
        "The 60 GWd/t stack stays unrun.",
        "The fabrication cap stays 5 percent.",
        "",
        "## Filed reload",
        "",
        f"The filed companion is `{_num(report.companion)}`.",
        "The filed companion stays the reload score.",
        "The fuel-side year at the demonstrated efficiency stays under the nameplate year.",
        f"The same reload at the plants-being-built efficiency stands at `{_num(report.built_companion)}` of the nameplate year.",
        "That quotient is the filed companion times the named built step.",
        "That quotient keeps the reload score and the built score.",
        "The fuel-side year at the plants-being-built efficiency stands above the nameplate year.",
        f"The filed reload stays `{_num(report.reload_tonnes)}` t.",
        "The demonstrated plant stays at its filed efficiency.",
        "The catalog ceiling stays unrun.",
        "The fast row's filed barrier is the sodium system.",
        "This score adds no joules.",
        "The published stamp records these gains as rationals of the named steps.",
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
