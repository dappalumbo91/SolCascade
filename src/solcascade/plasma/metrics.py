"""Zero-dimensional compression metrics for a spherical DT plasma.

Assumptions, stated so they can be replaced later:

- The plasma is a uniform sphere. Density scales with the inverse volume.
- Ions and electrons are quasineutral and singly charged, so n_e = n_i.
- The configured adiabatic law sets T V^(gamma-1) constant. For a monatomic
  gas gamma comes from three translational degrees of freedom.
- ``frozen_axial_flux`` holds B * R^2 fixed at the value of the vacuum coil
  field on axis at the first schedule time. ``vacuum_only`` uses the coil
  field at each time and does not amplify B by compression.
- Fusion power uses the Bosch-Hale thermal reactivity and a 50/50-capable
  deuterium fraction. There is no beam-target term.
- The energy confinement time is the smaller of the gyro-period random
  walk and the schedule duration. One thermal gyroradius per cyclotron
  period gives D = kT / (4 π e B). Bohm's 16 is the comparison, not the law.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from solcascade.compression.schedule import PulseSchedule
from solcascade.constants import (
    BOHM_DENOMINATOR,
    DT_TRIPLE_PRODUCT_BENCHMARK_KEV_S_M3,
    ELEMENTARY_CHARGE_C,
    JOULE_PER_EV,
    MONATOMIC_ADIABATIC_INDEX,
    VACUUM_PERMEABILITY_H_PER_M,
    dt_reaction,
)
from solcascade.geometry.coils import CoilSet
from solcascade.plasma.reactivity import dt_reactivity_m3_s


@dataclass(frozen=True)
class PlasmaPoint:
    time_s: float
    radius_m: float
    ion_density_m3: float
    ion_temperature_kev: float
    electron_temperature_kev: float
    vacuum_field_t: float
    internal_field_t: float
    beta: float
    reactivity_m3_s: float
    fusion_power_w: float
    alpha_power_w: float
    neutron_power_w: float
    neutron_rate_per_s: float
    tau_e_s: float
    triple_product_kev_s_m3: float
    triple_product_ratio: float
    pressure_pa: float


@dataclass(frozen=True)
class PlasmaHistory:
    points: tuple[PlasmaPoint, ...]
    flux_conservation: str
    adiabatic_index: float

    def column(self, name: str) -> np.ndarray:
        return np.array([getattr(point, name) for point in self.points], dtype=float)

    def integrate(self, name: str) -> float:
        values = self.column(name)
        times = self.column("time_s")
        return float(np.trapezoid(values, times))


def adiabatic_temperature_exponent(adiabatic_index: float) -> float:
    """T scales as R to this power for a sphere, from T V^(gamma-1) = constant."""
    return -3.0 * (adiabatic_index - 1.0)


def evolve(
    coils: CoilSet,
    schedule: PulseSchedule,
    ion_density_m3: float,
    ion_temperature_kev: float,
    deuterium_fraction: float,
    flux_conservation: str,
    electron_temperature_equals_ion: bool = True,
) -> PlasmaHistory:
    if not 0.0 < deuterium_fraction < 1.0:
        raise ValueError("deuterium fraction must be between 0 and 1")
    if ion_density_m3 <= 0.0 or ion_temperature_kev <= 0.0:
        raise ValueError("initial density and temperature must be positive")
    if flux_conservation not in {"frozen_axial_flux", "vacuum_only"}:
        raise ValueError(f"unknown flux conservation law {flux_conservation!r}")
    if not electron_temperature_equals_ion:
        raise NotImplementedError(
            "Separate electron temperature is not a slot in this version. "
            "Set plasma.electron_temperature_equals_ion true."
        )
    reaction = dt_reaction()
    gamma = MONATOMIC_ADIABATIC_INDEX
    # exponent is the power of radius. For a monatomic sphere it is -2,
    # so T / T0 = (R0 / R) ** (-exponent) = (R0 / R) ** 2.
    exponent = adiabatic_temperature_exponent(gamma)
    base = coils.base_currents_a()
    radius0 = float(schedule.plasma_radius_m[0])
    scales = schedule.current_scales(base.size)
    current0 = base * scales[0]
    vacuum0 = coils.on_axis_field(0.0, current0)
    points: list[PlasmaPoint] = []
    for time_s, radius, scale in zip(
        schedule.time_s, schedule.plasma_radius_m, scales, strict=True
    ):
        currents = base * scale
        vacuum = coils.on_axis_field(0.0, currents)
        if flux_conservation == "frozen_axial_flux":
            internal = vacuum0 * (radius0 / float(radius)) ** 2
        else:
            internal = vacuum
        density = ion_density_m3 * (radius0 / float(radius)) ** 3
        temperature = ion_temperature_kev * (radius0 / float(radius)) ** (-exponent)
        electron = temperature
        # kT in joules is temperature_keV * (1000 eV/keV) * (joule/eV).
        # Quasineutral Z=1 fuel: p = n_i k T_i + n_e k T_e with n_e = n_i.
        thermal_joule = _joule_per_kev() * (temperature + electron)
        pressure = density * thermal_joule
        beta = _beta(pressure, internal)
        reactivity = dt_reactivity_m3_s(temperature)
        n_d = deuterium_fraction * density
        n_t = (1.0 - deuterium_fraction) * density
        volume = (4.0 / 3.0) * math.pi * float(radius) ** 3
        power_density = n_d * n_t * reactivity * reaction.q_joule
        fusion = power_density * volume
        alpha = fusion * reaction.alpha_fraction
        neutron = fusion * reaction.neutron_fraction
        rate = neutron / reaction.neutron_energy_joule if neutron > 0.0 else 0.0
        tau = _confinement_time(float(radius), temperature, internal, schedule.duration_s)
        triple = density * temperature * tau
        points.append(
            PlasmaPoint(
                time_s=float(time_s),
                radius_m=float(radius),
                ion_density_m3=float(density),
                ion_temperature_kev=float(temperature),
                electron_temperature_kev=float(electron),
                vacuum_field_t=float(vacuum),
                internal_field_t=float(internal),
                beta=float(beta),
                reactivity_m3_s=float(reactivity),
                fusion_power_w=float(fusion),
                alpha_power_w=float(alpha),
                neutron_power_w=float(neutron),
                neutron_rate_per_s=float(rate),
                tau_e_s=float(tau),
                triple_product_kev_s_m3=float(triple),
                triple_product_ratio=float(triple / DT_TRIPLE_PRODUCT_BENCHMARK_KEV_S_M3),
                pressure_pa=float(pressure),
            )
        )
    return PlasmaHistory(
        points=tuple(points),
        flux_conservation=flux_conservation,
        adiabatic_index=gamma,
    )


def _joule_per_kev() -> float:
    return JOULE_PER_EV * 1.0e3


def _beta(pressure_pa: float, field_t: float) -> float:
    if field_t == 0.0:
        return math.inf
    magnetic = field_t ** 2 / (2.0 * VACUUM_PERMEABILITY_H_PER_M)
    return pressure_pa / magnetic


def bohm_time_s(radius_m: float, temperature_kev: float, field_t: float) -> float:
    """Bohm comparison time a^2 / (2 D_B), with the historical 16. Not the ledger time."""
    if radius_m <= 0.0 or temperature_kev <= 0.0 or field_t <= 0.0:
        return math.inf
    thermal = temperature_kev * _joule_per_kev()
    diffusion = thermal / (BOHM_DENOMINATOR * ELEMENTARY_CHARGE_C * field_t)
    if diffusion <= 0.0:
        return math.inf
    return radius_m ** 2 / (2.0 * diffusion)


def gyro_walk_time_s(radius_m: float, temperature_kev: float, field_t: float) -> float:
    """Confinement time for one thermal gyroradius per cyclotron period.

    The step is ρ = v/ω_c and the clock is τ_c = 2π/ω_c, with v² = kT/m.
    One-dimensional diffusion is D = ρ² / (2 τ_c), so

        D = kT / (4 π e B),    τ = a² / (2 D) = 2 π a² e B / kT.

    The 4π is that period. The ion mass cancels. This is π/4 of the Bohm
    time at every radius, temperature, and field. It is not capped by the pulse.
    """
    if radius_m <= 0.0 or temperature_kev <= 0.0 or field_t <= 0.0:
        return math.inf
    thermal = temperature_kev * _joule_per_kev()
    diffusion = thermal / (4.0 * math.pi * ELEMENTARY_CHARGE_C * field_t)
    if diffusion <= 0.0:
        return math.inf
    return radius_m ** 2 / (2.0 * diffusion)


def _confinement_time(
    radius_m: float, temperature_kev: float, field_t: float, pulse_s: float
) -> float:
    """Smaller of the gyro-period walk and the schedule duration."""
    if pulse_s <= 0.0:
        return 0.0
    return min(gyro_walk_time_s(radius_m, temperature_kev, field_t), pulse_s)
