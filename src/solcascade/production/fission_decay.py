"""Evaluate the 1952 decay-heat fit on the filed full-power year.

The cooling page states no fraction. Equation (2) of Untermyer and Weills,
AECD-3454, is the calorimetry fit for irradiated uranium. It includes U-239
and Np-239 decay and leaves delayed-neutron fission out. The run uses the
filed 8760-hour year and the filed nameplate heat. A foundational seed is
not multiplied into the watts. The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_balance import delayed_proportions
from solcascade.production.fission_paths import SeedGap, _nearest_seeds
from solcascade.production.fission_process import (
    FULL_POWER_HOURS,
    NAMEPLATE_MWE,
    build_process,
)
from solcascade.production.fission_safety import SMALL_CLASS_MWE
from solcascade.production.refine import project_root

DECAY_SOURCE = "https://www.osti.gov/servlets/purl/4376281"
# Equation (2) of AECD-3454, page 4. Times are seconds.
_LEAD = 0.1
_LONG_WEIGHT = 0.87
_SHORT_OFFSET_S = 10.0
_LONG_OFFSET_S = 2.0e7
_EXPONENT = -0.2
_SECOND_S = 1.0
_DAY_S = 86_400.0
_TWO_DAY_S = 172_800.0
_REVIEW_GAP = 1.0e-3


@dataclass(frozen=True)
class DecayPoint:
    label: str
    cooling_s: float
    band: float
    year_fraction: float
    limit_fraction: float
    year_watts: float
    limit_watts: float
    small_watts: float
    year_seed: SeedGap
    limit_seed: SeedGap


@dataclass(frozen=True)
class DecayReport:
    operation_s: float
    thermal_w: float
    small_thermal_w: float
    points: tuple[DecayPoint, ...]
    day_over_second: float
    two_over_day: float
    day_seed: SeedGap
    two_seed: SeedGap
    delayed_gap: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_decay.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _brace(time_s: float) -> float:
    """One brace of equation (2). The second brace is this function at Ts+To."""
    return (time_s + _SHORT_OFFSET_S) ** _EXPONENT - _LONG_WEIGHT * (
        time_s + _LONG_OFFSET_S
    ) ** _EXPONENT


def decay_fraction(cooling_s: float, operation_s: float) -> float:
    """Ps/Po from equation (2). Po is the operating heat rate."""
    if cooling_s < _SECOND_S:
        raise RuntimeError("the 1952 report assigns a large error under one second")
    if operation_s < 0.0:
        raise RuntimeError("the operation time is negative")
    return _LEAD * (_brace(cooling_s) - _brace(cooling_s + operation_s))


def limit_fraction(cooling_s: float) -> float:
    """Equation (2) as the operation time grows and the second brace vanishes."""
    if cooling_s < _SECOND_S:
        raise RuntimeError("the 1952 report assigns a large error under one second")
    return _LEAD * _brace(cooling_s)


def stated_band(cooling_s: float) -> float:
    """Fractional band printed by the 1952 report for that cooling time."""
    if cooling_s < _SECOND_S or cooling_s >= 1.0e8:
        raise RuntimeError("the cooling time is outside the report's stated bands")
    if cooling_s < 100.0:
        return 0.50
    if cooling_s < 1.0e4:
        return 0.30
    if cooling_s < 1.0e6:
        return 0.10
    return 0.50


def _watch(module, share: float, label: str) -> SeedGap:
    seed = _nearest_seeds(module, share)[0]
    if seed.relative_gap <= _REVIEW_GAP:
        raise RuntimeError(f"{seed.name} meets the {label} closely enough to review")
    return seed


def build_decay() -> DecayReport:
    process = build_process()
    operation_s = FULL_POWER_HOURS * 3600.0
    if abs(operation_s - 365.0 * _DAY_S) > 1.0:
        raise RuntimeError("the filed year is not 365 civil days")
    if abs(process.nameplate_j - NAMEPLATE_MWE * 1.0e6 * operation_s) > 1.0:
        raise RuntimeError("the nameplate year does not match the process score")
    if not 0.0 < SMALL_CLASS_MWE < NAMEPLATE_MWE:
        raise RuntimeError("the small class is not under the nameplate")
    if not 0.0 < process.efficiency < 1.0:
        raise RuntimeError("the efficiency is outside (0, 1)")
    thermal_w = NAMEPLATE_MWE * 1.0e6 / process.efficiency
    small_thermal_w = SMALL_CLASS_MWE * 1.0e6 / process.efficiency
    module, _digest = load_engine()
    schedule = (
        ("one second", _SECOND_S),
        ("one day", _DAY_S),
        ("two days", _TWO_DAY_S),
    )
    points: list[DecayPoint] = []
    previous = 1.0
    for label, cooling_s in schedule:
        band = stated_band(cooling_s)
        year = decay_fraction(cooling_s, operation_s)
        limit = limit_fraction(cooling_s)
        if not 0.0 < year < limit < 1.0:
            raise RuntimeError(f"the {label} fraction is outside the operating heat")
        if year >= previous:
            raise RuntimeError("the decay load does not fall as the cooling time grows")
        previous = year
        points.append(
            DecayPoint(
                label=label,
                cooling_s=cooling_s,
                band=band,
                year_fraction=year,
                limit_fraction=limit,
                year_watts=year * thermal_w,
                limit_watts=limit * thermal_w,
                small_watts=year * small_thermal_w,
                year_seed=_watch(module, year, f"{label} fraction"),
                limit_seed=_watch(module, limit, f"{label} limit"),
            )
        )
    second, day, two = points
    if second.band != 0.50 or day.band != 0.10 or two.band != 0.10:
        raise RuntimeError("the scored times are not in the report's stated bands")
    day_over_second = day.year_fraction / second.year_fraction
    two_over_day = two.year_fraction / day.year_fraction
    if not 0.0 < day_over_second < two_over_day < 1.0:
        raise RuntimeError("the early drop is not steeper than the later drop")
    u235 = next(row for row in delayed_proportions() if row.nucleus == "U-235")
    delayed_gap = abs(day.year_fraction - u235.proportion) / u235.proportion
    if delayed_gap <= _REVIEW_GAP:
        raise RuntimeError(
            "the one-day fraction meets the delayed-neutron proportion closely enough to review"
        )
    return DecayReport(
        operation_s=operation_s,
        thermal_w=thermal_w,
        small_thermal_w=small_thermal_w,
        points=tuple(points),
        day_over_second=day_over_second,
        two_over_day=two_over_day,
        day_seed=_watch(module, day_over_second, "one-day ratio"),
        two_seed=_watch(module, two_over_day, "two-day ratio"),
        delayed_gap=delayed_gap,
        nuclear_scalar=process.nuclear_scalar,
        thermodynamics_scalar=process.thermodynamics_scalar,
        gamma_fsot=process.gamma_fsot,
        valve_beta=process.valve_beta,
        pin_prefix=process.pin_prefix,
        digest_prefix=process.digest_prefix,
    )


def _point_lines(point: DecayPoint) -> list[str]:
    percent = point.band * 100.0
    return [
        f"At {point.label} after the filed full-power year, the fit gives `{_num(point.year_fraction)}` of operating heat.",
        f"The report's stated band at that time is ±{_num(percent)} percent.",
        f"That cooling load is `{_num(point.year_watts)}` W on the large nameplate.",
        f"On the small-class boundary the same fraction is `{_num(point.small_watts)}` W.",
        f"The formula's limit as operation time grows is `{_num(point.limit_fraction)}` of operating heat, `{_num(point.limit_watts)}` W.",
        f"The nearest seed to the one-year fraction is `{point.year_seed.name}` at `{_num(point.year_seed.value)}`.",
        f"Its relative gap is `{_num(point.year_seed.relative_gap)}`.",
        f"The nearest seed to the formula's limit is `{point.limit_seed.name}` at `{_num(point.limit_seed.value)}`.",
        f"Its relative gap is `{_num(point.limit_seed.relative_gap)}`.",
    ]


def render() -> str:
    report = build_decay()
    lines = [
        "# Fission decay",
        "",
        f"The decay-heat fit is equation (2) of Untermyer and Weills, AECD-3454, 25 February 1952, {DECAY_SOURCE}.",
        "The cooling page states no decay-heat fraction.",
        "The fraction below is that 1952 fit, evaluated on the filed plant.",
        "It is a cooling load.",
        "It is not a second electricity product.",
        "The fit includes the decay heat of U-239 and Np-239.",
        "The 1952 report leaves delayed-neutron fission out of the fit.",
        "The delayed-neutron margin stays in the core.",
        "The report states that equation (2) sits about 22 percent above its fission-product sum.",
        "The fit is scored as published, with U-239 and Np-239 left in.",
        f"The filed operation time is `{_num(report.operation_s)}` s, one full-power year.",
        f"Operating heat on the large nameplate is `{_num(report.thermal_w)}` W.",
        f"Operating heat on the small-class boundary is `{_num(report.small_thermal_w)}` W.",
        "",
        "## Cooling times",
        "",
    ]
    for point in report.points:
        lines.extend(_point_lines(point))
        lines.append("")
    lines.extend(
        [
            "## How the load falls",
            "",
            f"From one second to one day the one-year load falls to `{_num(report.day_over_second)}` of the one-second load.",
            f"The nearest seed to that ratio is `{report.day_seed.name}` at `{_num(report.day_seed.value)}`.",
            f"Its relative gap is `{_num(report.day_seed.relative_gap)}`.",
            f"From one day to two days the one-year load falls to `{_num(report.two_over_day)}` of the one-day load.",
            f"The nearest seed to that ratio is `{report.two_seed.name}` at `{_num(report.two_seed.value)}`.",
            f"Its relative gap is `{_num(report.two_seed.relative_gap)}`.",
            "The closest of these seeds is still wider than one part in a thousand.",
            "No seed is installed as a decay-heat fraction.",
            "The one-day figure counts heat.",
            "The U-235 figure of 0.66 percent counts delayed neutrons.",
            f"Their relative gap is `{_num(report.delayed_gap)}`.",
            "That gap is wider than one part in a thousand.",
            "The one-day heat ratio stays off that neutron count.",
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
