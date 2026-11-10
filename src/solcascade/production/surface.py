"""Replace the distant windings with a spherical surface current.

A surface current K = (3/2) (B / μ0) sinθ φ-hat makes a uniform field inside
a sphere and a dipole outside. The stored energy is 3/2 of the internal
magnetic energy. The quadrature in this module checks the dipole half.

Three arrangements use that integral. A current on the plasma surface holds
the compressed field in the plasma volume. A current on the shell inner wall
holds the vacuum field in the bore, and the checked-in compression still
raises the field in the plasma. The thermal transit is the clock with the
field removed. The checked-in YAML is not rewritten.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from solcascade.constants import (
    JOULE_PER_EV,
    MASS_EXCESS_DEUTERIUM_KEV,
    MASS_EXCESS_TRITIUM_KEV,
    VACUUM_PERMEABILITY_H_PER_M,
    ion_mass_kg,
)
from solcascade.plasma.metrics import gyro_walk_time_s
from solcascade.production.cost import CostReport, build_cost
from solcascade.production.refine import Problem, load_checked_in, project_root

# Exterior dipole energy over internal energy for this surface current.
DIPOLE_EXTERIOR_OVER_INTERIOR = 0.5
SURFACE_OVER_INTERIOR = 1.0 + DIPOLE_EXTERIOR_OVER_INTERIOR


@dataclass(frozen=True)
class SurfaceDesign:
    label: str
    radius_m: float
    field_t: float
    magnetic_j: float
    pulsed_supply_j: float
    equator_current_a_per_m: float


@dataclass(frozen=True)
class SurfaceReport:
    helmholtz_stored_j: float
    helmholtz_practical_j: float
    recovery: float
    beta1_complete_electric_j: float
    beta1_schedule_electric_j: float
    beta1_walk_electric_j: float
    beta1_walk_time_s: float
    transit_time_s: float
    transit_electric_j: float
    bore_schedule_electric_j: float
    bore_walk_electric_j: float
    bore_complete_electric_j: float
    bore_electric_power_w: float
    plasma: SurfaceDesign
    wall: SurfaceDesign
    bore_plasma: SurfaceDesign
    quadrature_exterior_over_interior: float


def default_output_path():
    return project_root() / "runs" / "production" / "surface_field.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _side(left: float, right: float) -> str:
    if left > right:
        return "above"
    if left < right:
        return "below"
    return "level with"


def thermal_speed_m_s(temperature_kev: float) -> float:
    """Mean deuterium-tritium ion speed, v = sqrt(kT / m)."""
    mass = 0.5 * (
        ion_mass_kg(2, MASS_EXCESS_DEUTERIUM_KEV, 1)
        + ion_mass_kg(3, MASS_EXCESS_TRITIUM_KEV, 1)
    )
    thermal = float(temperature_kev) * JOULE_PER_EV * 1.0e3
    return math.sqrt(thermal / mass)


def internal_energy_j(radius_m: float, field_t: float) -> float:
    volume = (4.0 / 3.0) * math.pi * radius_m ** 3
    return field_t ** 2 / (2.0 * VACUUM_PERMEABILITY_H_PER_M) * volume


def equator_current_a_per_m(field_t: float) -> float:
    """Peak of K = (3/2) (B / μ0) sinθ, from the tangential H jump."""
    return 1.5 * field_t / VACUUM_PERMEABILITY_H_PER_M


def exterior_over_interior_quadrature(n_radial: int = 800, n_polar: int = 800) -> float:
    """∫ B² dV outside the sphere, divided by B² times the interior volume.

    The exterior field is the dipole matched to a uniform interior field.
    The radial grid is logarithmic. The result is dimensionless.
    """
    if n_radial < 10 or n_polar < 10:
        raise ValueError("quadrature needs a real grid")
    # r = exp(u). The integrand r^{-4} dr becomes exp(-3u) du, which the
    # trapezoid on a uniform u grid resolves. A trapezoid in r itself
    # overestimates this convex integrand.
    u_max = math.log(1.0e6)
    u = np.linspace(0.0, u_max, n_radial)
    theta = np.linspace(0.0, math.pi, n_polar)
    angular = np.cos(theta) ** 2 + 0.25 * np.sin(theta) ** 2
    angular_integral = np.trapezoid(angular * np.sin(theta), theta)
    radial_integral = np.trapezoid(np.exp(-3.0 * u), u)
    # ∫ (a/r)^6 (cos²θ + sin²θ/4) r² sinθ dθ dφ dr with a = 1.
    integral = radial_integral * angular_integral * 2.0 * math.pi
    volume = (4.0 / 3.0) * math.pi
    return float(integral / volume)


def _design(label: str, radius_m: float, field_t: float, recovery: float) -> SurfaceDesign:
    magnetic = SURFACE_OVER_INTERIOR * internal_energy_j(radius_m, field_t)
    return SurfaceDesign(
        label=label,
        radius_m=float(radius_m),
        field_t=float(field_t),
        magnetic_j=float(magnetic),
        pulsed_supply_j=float((1.0 - recovery) * magnetic),
        equator_current_a_per_m=float(equator_current_a_per_m(field_t)),
    )


def build_surface(problem: Problem | None = None) -> SurfaceReport:
    problem = load_checked_in() if problem is None else problem
    cost: CostReport = build_cost(problem)
    recovery = float(problem.slots["conversion.coil_energy_recovery_fraction"])
    bore = cost.bore
    compression = (bore.initial_radius_m / bore.peak_radius_m) ** 2
    vacuum_field = bore.peak_field_t / compression
    small_radius = bore.peak_radius_m / bore.radius_factor
    small_power = bore.peak_power_w / bore.radius_factor ** 3
    walk = gyro_walk_time_s(small_radius, bore.peak_temperature_kev, bore.peak_field_t)
    speed = thermal_speed_m_s(bore.peak_temperature_kev)
    transit = small_radius / speed
    eta = bore.electric_power_w / bore.peak_power_w
    quadrature = exterior_over_interior_quadrature()
    if abs(quadrature - DIPOLE_EXTERIOR_OVER_INTERIOR) > 1.0e-3:
        raise RuntimeError(f"dipole exterior fraction quadrature is {quadrature}")
    return SurfaceReport(
        helmholtz_stored_j=cost.helmholtz_stored_j,
        helmholtz_practical_j=cost.practical_j,
        recovery=recovery,
        beta1_complete_electric_j=cost.beta1_complete_electric_j,
        beta1_schedule_electric_j=cost.beta1_schedule_electric_j,
        beta1_walk_electric_j=cost.beta1_walk_electric_j,
        beta1_walk_time_s=float(walk),
        transit_time_s=float(transit),
        transit_electric_j=float(eta * small_power * transit),
        bore_schedule_electric_j=bore.electric_j,
        bore_walk_electric_j=bore.walk_electric_j,
        bore_complete_electric_j=bore.complete_electric_j,
        bore_electric_power_w=bore.electric_power_w,
        plasma=_design("plasma surface", small_radius, bore.peak_field_t, recovery),
        wall=_design("shell inner wall", cost.shell_inner_m, vacuum_field, recovery),
        bore_plasma=_design("bore plasma surface", bore.peak_radius_m, bore.peak_field_t, recovery),
        quadrature_exterior_over_interior=quadrature,
    )


def render(problem: Problem | None = None) -> str:
    report = build_surface(problem)
    plasma = report.plasma
    wall = report.wall
    bore = report.bore_plasma
    lines = [
        "# Surface field",
        "",
        "The distant windings are the large store. A spherical surface current",
        "replaces them. Inside the sphere the field is uniform. Outside it is",
        "the dipole that matches the boundary. The stored energy is three halves",
        "of the internal magnetic energy. The exterior half is the dipole",
        f"integral, and the quadrature gives `{_num(report.quadrature_exterior_over_interior)}`.",
        "",
        "The current density on the sphere is `(3/2) (B / μ0) sinθ` in the",
        "azimuthal direction. A persistent current is already established, so",
        "the pulse supply does not recharge it. A pulsed current at the",
        f"configured recovery `{_num(report.recovery)}` leaves the unrecovered",
        "fraction of that magnetic energy as the supply.",
        "",
        "## Arrangements",
        "",
        "| Arrangement | Radius (m) | Field (T) | Magnetic energy (J) | Pulsed supply (J) | Equator current (A/m) |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
        f"| Helmholtz windings |  |  | {report.helmholtz_stored_j:.6g} | {(1.0 - report.recovery) * report.helmholtz_stored_j:.6g} |  |",
        f"| {wall.label} | {wall.radius_m:.6g} | {wall.field_t:.6g} | {wall.magnetic_j:.6g} | {wall.pulsed_supply_j:.6g} | {wall.equator_current_a_per_m:.6g} |",
        f"| {plasma.label} | {plasma.radius_m:.6g} | {plasma.field_t:.6g} | {plasma.magnetic_j:.6g} | {plasma.pulsed_supply_j:.6g} | {plasma.equator_current_a_per_m:.6g} |",
        f"| {bore.label} | {bore.radius_m:.6g} | {bore.field_t:.6g} | {bore.magnetic_j:.6g} | {bore.pulsed_supply_j:.6g} | {bore.equator_current_a_per_m:.6g} |",
        "",
        "The shell-wall current fills the bore at the vacuum field. The",
        "checked-in compression still raises the field inside the plasma.",
        "The plasma-surface current holds the compressed field in the plasma",
        "volume and has no external winding. The bore row is that plasma-surface",
        "current on the larger sphere that fits in the shell.",
        "",
        "The pulsed column is the unrecovered share of the magnetic energy.",
        "The Helmholtz practical supply, ohmic ledger included, is "
        f"`{_num(report.helmholtz_practical_j)}` J.",
        "",
        "Shell-wall magnetic energy versus the Helmholtz store: "
        f"{_side(wall.magnetic_j, report.helmholtz_stored_j)}.",
        "Plasma-surface magnetic energy versus the shell-wall energy: "
        f"{_side(plasma.magnetic_j, wall.magnetic_j)}.",
        "",
        "## Ledger",
        "",
        f"Beta = 1 schedule electricity is `{_num(report.beta1_schedule_electric_j)}` J.",
        f"The gyro-walk electricity is `{_num(report.beta1_walk_electric_j)}` J",
        f"in `{_num(report.beta1_walk_time_s)}` s. Complete-burn electricity of that",
        f"load is `{_num(report.beta1_complete_electric_j)}` J.",
        "",
        "Complete burn versus the pulsed plasma-surface supply: "
        f"{_side(report.beta1_complete_electric_j, plasma.pulsed_supply_j)}.",
        "Schedule electricity versus the pulsed plasma-surface supply: "
        f"{_side(report.beta1_schedule_electric_j, plasma.pulsed_supply_j)}.",
        "Complete burn versus the pulsed shell-wall supply: "
        f"{_side(report.beta1_complete_electric_j, wall.pulsed_supply_j)}.",
        "Schedule electricity versus the pulsed shell-wall supply: "
        f"{_side(report.beta1_schedule_electric_j, wall.pulsed_supply_j)}.",
        "",
        "A persistent surface current draws no magnetic energy from the supply",
        "on the pulse. The interior field is uniform, so the axis is open and",
        "the parallel clock is the thermal transit. The schedule integral and",
        "the gyro-walk remain the cross-field figures.",
        f"On the bore-limited sphere those figures are `{_num(report.bore_schedule_electric_j)}` J",
        f"and `{_num(report.bore_walk_electric_j)}` J. Complete-burn electricity is",
        f"`{_num(report.bore_complete_electric_j)}` J. Peak electric power is",
        f"`{_num(report.bore_electric_power_w)}` W.",
        "Bore complete burn versus the pulsed bore-surface supply: "
        f"{_side(report.bore_complete_electric_j, bore.pulsed_supply_j)}.",
        "Bore schedule electricity versus the pulsed bore-surface supply: "
        f"{_side(report.bore_schedule_electric_j, bore.pulsed_supply_j)}.",
        "",
        "## No field",
        "",
        "With the field removed, the ions cross the sphere at the thermal speed",
        f"`v = sqrt(kT / m)`. The transit time is `{_num(report.transit_time_s)}` s",
        f"and the electricity in that transit is `{_num(report.transit_electric_j)}` J.",
        "The thermal transit is shorter than the gyro-walk.",
        "The magnetic supply on that arrangement is zero, and the clock is the transit.",
        "",
        "The checked-in schedule is not rewritten. The distant windings remain",
        "the ledger until a surface current is built.",
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
