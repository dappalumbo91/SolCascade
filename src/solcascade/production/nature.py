"""Score the axial throat of the fields this architecture already builds.

The gyro-walk is the cross-field clock installed in ``evolve``. A particle
leaves along the field on the thermal transit unless a throat is stronger
than the field it occupies. The uniform surface current, the Helmholtz pair,
and the checked-in four-coil set are scored with that comparison. The
checked-in YAML is not rewritten.

FSOT scalars are recorded beside the ratios. A domain scalar is not a
mirror ratio, and the engine valve is not inverted into joules.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from solcascade.assemble import build_shell_from_slots
from solcascade.constants import (
    ELEMENTARY_CHARGE_C,
    JOULE_PER_EV,
    MASS_EXCESS_DEUTERIUM_KEV,
    MASS_EXCESS_TRITIUM_KEV,
    ion_mass_kg,
)
from solcascade.energy_baseline.phase2.map_fsot import (
    CAPTURED_DOMAINS,
    EXTRA_OUTPUTS,
    PIN_PREFIX,
    capture_scalars,
    load_engine,
)
from solcascade.plasma.metrics import gyro_walk_time_s
from solcascade.production.cost import _helmholtz
from solcascade.production.refine import Problem, load_checked_in, project_root
from solcascade.production.surface import SurfaceReport, build_surface, thermal_speed_m_s

# Points on the axis from the midplane to the spherical shell.
AXIAL_SAMPLES = 2001
WITNESS_DOMAINS = (
    "Electromagnetism",
    "Geophysics",
    "Fluid_Dynamics",
    "Nuclear_Physics",
)


@dataclass(frozen=True)
class AxialField:
    """One installed field, scored from the midplane out to the wall."""

    label: str
    vacuum_rise: float
    throat_over_internal: float
    vacuum_loss_fraction: float
    internal_loss_fraction: float
    axial_electric_j: float
    trapped_electric_j: float


@dataclass(frozen=True)
class NatureReport:
    gyroradii: float
    bore_gyroradii: float
    walk_time_s: float
    walk_electric_j: float
    transit_time_s: float
    transit_electric_j: float
    bore_transit_time_s: float
    bore_transit_electric_j: float
    bore_walk_electric_j: float
    four_coil: AxialField
    helmholtz: AxialField
    surface_current: AxialField
    pin_prefix: str
    digest_prefix: str
    valve_beta: float
    scalars: dict[str, float]


def default_output_path():
    return project_root() / "runs" / "production" / "nature.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _side(left: float, right: float) -> str:
    if left > right:
        return "above"
    if left < right:
        return "below"
    return "level with"


def mean_ion_mass_kg() -> float:
    """Mean of the deuterium and tritium ion masses used by the thermal speed."""
    return 0.5 * (
        ion_mass_kg(2, MASS_EXCESS_DEUTERIUM_KEV, 1)
        + ion_mass_kg(3, MASS_EXCESS_TRITIUM_KEV, 1)
    )


def isotropic_loss_fraction(mirror_ratio: float) -> float:
    """Share of an isotropic distribution inside the loss cone.

    A particle at the midplane mirrors when sin²θ >= 1/R and R is the
    strongest field on its path divided by the field it occupies. The
    isotropic share with a smaller pitch angle is 1 - sqrt(1 - 1/R).
    A path whose strongest field is the occupied field has no trapped set.
    """
    if mirror_ratio < 1.0:
        raise ValueError("mirror ratio is the strongest field over the occupied field")
    if mirror_ratio == 1.0:
        return 1.0
    trapped = math.sqrt(1.0 - 1.0 / mirror_ratio)
    return 1.0 - trapped


def dipole_axis_ratio(radius_over_sphere: float) -> float:
    """On-axis field of the surface-current dipole, in units of the interior field.

    Inside the sphere the field is uniform. Outside, the axis of the matching
    dipole falls as (a/r)³.
    """
    if radius_over_sphere <= 0.0:
        raise ValueError("radius must be positive")
    if radius_over_sphere <= 1.0:
        return 1.0
    return radius_over_sphere ** (-3)


def gyroradii_across(radius_m: float, temperature_kev: float, field_t: float) -> float:
    """Sphere radius over one thermal gyroradius, ρ = sqrt(kT m) / (e B)."""
    if radius_m <= 0.0 or temperature_kev <= 0.0 or field_t <= 0.0:
        raise ValueError("gyroradius needs a positive radius, temperature, and field")
    speed = thermal_speed_m_s(temperature_kev)
    omega = ELEMENTARY_CHARGE_C * float(field_t) / mean_ion_mass_kg()
    return float(radius_m) * omega / speed


def axial_peak(
    coils, currents: np.ndarray, z_wall_m: float, samples: int = AXIAL_SAMPLES
) -> tuple[float, float]:
    """Strongest |B| on the axis from the midplane to the shell, and that rise.

    The spherical shell meets the axis at ``z_wall_m``. The rise is the
    strongest field divided by the midplane field. A rise of 1 means the
    midplane is the strongest field on that path.
    """
    if samples < 11:
        raise ValueError("axial mirror needs a real grid")
    if z_wall_m <= 0.0:
        raise ValueError("the wall meets the axis at a positive radius")
    if len(currents) != len(coils.coils):
        raise ValueError("current vector length does not match the coil set")
    z = np.linspace(0.0, float(z_wall_m), samples)
    field = np.empty(samples, dtype=float)
    for index, position in enumerate(z):
        field[index] = abs(float(coils.on_axis_field(float(position), currents)))
    mid = float(field[0])
    if mid == 0.0:
        raise RuntimeError("midplane field is zero")
    peak = float(np.max(field))
    rise = peak / mid
    if rise <= 1.0 + 1.0e-9:
        rise = 1.0
        peak = mid
    return rise, peak


def _axial(
    label: str,
    vacuum_rise: float,
    throat_t: float,
    internal_t: float,
    transit_electric_j: float,
    walk_electric_j: float,
) -> AxialField:
    if internal_t <= 0.0 or throat_t < 0.0:
        raise ValueError("the throat comparison needs a positive internal field")
    ratio = throat_t / internal_t
    # The occupied field is the internal field. A weaker throat does not mirror.
    occupied = ratio if ratio > 1.0 else 1.0
    loss = isotropic_loss_fraction(occupied)
    trapped = 1.0 - loss
    vacuum_loss = isotropic_loss_fraction(vacuum_rise if vacuum_rise > 1.0 else 1.0)
    stay = max(0.0, walk_electric_j)
    return AxialField(
        label=label,
        vacuum_rise=float(vacuum_rise),
        throat_over_internal=float(ratio),
        vacuum_loss_fraction=float(vacuum_loss),
        internal_loss_fraction=float(loss),
        axial_electric_j=float(transit_electric_j),
        trapped_electric_j=float(trapped ** 2 * stay),
    )


def _witness() -> tuple[str, float, dict[str, float]]:
    module, digest = load_engine()
    if not str(digest).startswith(PIN_PREFIX):
        raise RuntimeError(f"engine pin is {str(digest)[:16]}")
    scalars = capture_scalars(module)
    valve = float(module.BETA)
    if valve <= 0.0 or valve >= 1.0:
        raise RuntimeError("engine valve is outside the recorded range")
    return str(digest)[:16], valve, scalars


def build_nature(problem: Problem | None = None) -> NatureReport:
    problem = load_checked_in() if problem is None else problem
    surface: SurfaceReport = build_surface(problem)
    shell = build_shell_from_slots(problem.config, problem.slots)
    coils = problem.coils
    base = coils.base_currents_a()
    first = base * problem.schedule.current_scales(base.size)[0]
    vacuum = float(coils.on_axis_field(0.0, first))
    internal = surface.plasma.field_t
    if vacuum <= 0.0:
        raise RuntimeError("checked-in midplane field is zero")
    if abs(vacuum - surface.wall.field_t) / vacuum > 1.0e-9:
        raise RuntimeError("surface vacuum field disagrees with the coil axis")
    # The walk identity kT = 2 π a² e B / τ uses the installed clock.
    radius = surface.plasma.radius_m
    walk = surface.beta1_walk_time_s
    thermal = 2.0 * math.pi * ELEMENTARY_CHARGE_C * internal * radius ** 2 / walk
    temperature_kev = thermal / (JOULE_PER_EV * 1.0e3)
    speed = thermal_speed_m_s(temperature_kev)
    transit_speed = radius / surface.transit_time_s
    if abs(speed - transit_speed) / transit_speed > 1.0e-6:
        raise RuntimeError("walk temperature and the transit speed disagree")
    installed = gyro_walk_time_s(radius, temperature_kev, internal)
    if abs(installed - walk) / walk > 1.0e-9:
        raise RuntimeError("reported walk disagrees with gyro_walk_time_s")
    four_rise, four_peak = axial_peak(coils, first, shell.inner_radius_m)
    pair, pair_currents = _helmholtz(shell.outer_radius_m, vacuum)
    helm_rise, helm_peak = axial_peak(pair, pair_currents, shell.inner_radius_m)
    digest, valve, scalars = _witness()
    bore_speed_time = surface.bore_plasma.radius_m / speed
    return NatureReport(
        gyroradii=gyroradii_across(radius, temperature_kev, internal),
        bore_gyroradii=gyroradii_across(surface.bore_plasma.radius_m, temperature_kev, internal),
        walk_time_s=walk,
        walk_electric_j=surface.beta1_walk_electric_j,
        transit_time_s=surface.transit_time_s,
        transit_electric_j=surface.transit_electric_j,
        bore_transit_time_s=float(bore_speed_time),
        bore_transit_electric_j=float(surface.bore_electric_power_w * bore_speed_time),
        bore_walk_electric_j=surface.bore_walk_electric_j,
        four_coil=_axial(
            "four-coil",
            four_rise,
            four_peak,
            internal,
            surface.transit_electric_j,
            surface.beta1_walk_electric_j,
        ),
        helmholtz=_axial(
            "Helmholtz",
            helm_rise,
            helm_peak,
            internal,
            surface.transit_electric_j,
            surface.beta1_walk_electric_j,
        ),
        surface_current=_axial(
            "surface current",
            dipole_axis_ratio(1.0),
            internal * dipole_axis_ratio(1.0),
            internal,
            surface.transit_electric_j,
            surface.beta1_walk_electric_j,
        ),
        pin_prefix=PIN_PREFIX,
        digest_prefix=digest,
        valve_beta=valve,
        scalars=scalars,
    )


def render(problem: Problem | None = None) -> str:
    report = build_nature(problem)
    four = report.four_coil
    helm = report.helmholtz
    surface = report.surface_current
    lines = [
        "# Nature",
        "",
        "The gyro-walk remains the cross-field clock. Along the field, a",
        "particle mirrors only when some throat is stronger than the field",
        "it occupies. The isotropic share that misses that throat is",
        "`1 - sqrt(1 - 1/R)`. A path whose strongest field is the occupied",
        "field has no trapped set, and those particles leave on the thermal transit.",
        "",
        f"The beta = 1 sphere is `{_num(report.gyroradii)}` thermal gyroradii across.",
        f"The bore-limited sphere is `{_num(report.bore_gyroradii)}` thermal gyroradii across.",
        f"The installed walk is `{_num(report.walk_time_s)}` s and its electricity is",
        f"`{_num(report.walk_electric_j)}` J. The thermal transit is",
        f"`{_num(report.transit_time_s)}` s and its electricity is",
        f"`{_num(report.transit_electric_j)}` J.",
        "",
        "## Throats",
        "",
        "Vacuum rise is the strongest on-axis field between the midplane and",
        "the shell, divided by the midplane field. Throat over internal is",
        "that strongest field divided by the compressed field in the plasma.",
        "The surface-current row is the uniform interior: the matching dipole",
        "falls as `(a/r)³` on the axis outside the sphere.",
        "",
        "| Field | Vacuum rise | Throat / internal | Vacuum loss | Internal loss | Axial electricity (J) |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
        (
            f"| {four.label} | {four.vacuum_rise:.6g} | {four.throat_over_internal:.6g} | "
            f"{four.vacuum_loss_fraction:.6g} | {four.internal_loss_fraction:.6g} | "
            f"{four.axial_electric_j:.6g} |"
        ),
        (
            f"| {helm.label} | {helm.vacuum_rise:.6g} | {helm.throat_over_internal:.6g} | "
            f"{helm.vacuum_loss_fraction:.6g} | {helm.internal_loss_fraction:.6g} | "
            f"{helm.axial_electric_j:.6g} |"
        ),
        (
            f"| {surface.label} | {surface.vacuum_rise:.6g} | {surface.throat_over_internal:.6g} | "
            f"{surface.vacuum_loss_fraction:.6g} | {surface.internal_loss_fraction:.6g} | "
            f"{surface.axial_electric_j:.6g} |"
        ),
        "",
        "The four-coil vacuum axis rises above the midplane. "
        "The compressed plasma sits above the throat. The Helmholtz axis is strongest at the",
        "midplane. The surface current is uniform inside the sphere and weaker",
        "on the dipole axis outside it.",
        "",
        "An isotropic fuel leaves on the thermal transit.",
        "Axial electricity versus the gyro-walk electricity: "
        f"{_side(four.axial_electric_j, report.walk_electric_j)}.",
        "Trapped electricity on the four-coil internal field is "
        f"`{_num(four.trapped_electric_j)}` J.",
        f"On the bore-limited sphere the transit is `{_num(report.bore_transit_time_s)}` s",
        f"and its electricity is `{_num(report.bore_transit_electric_j)}` J.",
        "Bore axial electricity versus the bore gyro-walk electricity: "
        f"{_side(report.bore_transit_electric_j, report.bore_walk_electric_j)}.",
        "",
        "## FSOT record",
        "",
        f"Pin `{report.pin_prefix}` matches the loaded engine prefix",
        f"`{report.digest_prefix}`. The scalars below are dimensionless.",
        "A domain scalar is not a mirror ratio.",
        "The engine valve stays out of the mirror ratio.",
        f"Engine valve: `{_num(report.valve_beta)}`.",
        "",
        "| Scalar | Value |",
        "| --- | ---: |",
    ]
    for name in (*CAPTURED_DOMAINS, *EXTRA_OUTPUTS):
        lines.append(f"| {name} | {report.scalars[name]:.6g} |")
    lines.extend(
        [
            "",
            "Electromagnetism is the field witness. The throat ratios are the",
            "coil field and the dipole axis. Geophysics is the bulk witness and",
            "Fluid_Dynamics is the tank witness. Neither scalar is given a",
            "speed. Nuclear_Physics is the orifice witness. It is not a",
            "confinement time. gamma_FSOT does not multiply the walk.",
            "",
            "The checked-in schedule is not rewritten. The gyro-walk remains",
            "the cross-field clock installed on the ledger.",
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
