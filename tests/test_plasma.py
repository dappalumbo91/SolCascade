"""Compression scalings, reactivity range, and the supply ledger."""

import math

import numpy as np
import pytest

from solcascade.constants import MONATOMIC_ADIABATIC_INDEX
from solcascade.conversion import carnot_efficiency, convert, supply_energy
from solcascade.geometry.coils import Coil, CoilSet
from solcascade.plasma.metrics import (
    adiabatic_temperature_exponent,
    bohm_time_s,
    evolve,
    gyro_walk_time_s,
)
from solcascade.plasma.reactivity import dt_reactivity_m3_s
from solcascade.compression.schedule import PulseSchedule, build_schedule


def test_gyro_walk_is_pi_over_four_of_bohm_at_every_state():
    samples = ((0.018, 3.858, 1.67), (0.05, 0.5, 0.216), (0.02, 10.0, 5.0))
    for radius, temperature, field in samples:
        walk = gyro_walk_time_s(radius, temperature, field)
        bohm = bohm_time_s(radius, temperature, field)
        assert walk == pytest.approx(bohm * math.pi / 4.0)


def test_monatomic_sphere_heats_as_radius_inverse_square():
    exponent = adiabatic_temperature_exponent(MONATOMIC_ADIABATIC_INDEX)
    assert exponent == pytest.approx(-2.0)


def test_reactivity_at_10_kev_is_in_the_published_band():
    reactivity_cm3 = dt_reactivity_m3_s(10.0) / 1.0e-6
    assert 8.0e-17 < reactivity_cm3 < 1.5e-16
    with pytest.raises(ValueError):
        dt_reactivity_m3_s(0.01)


def test_frozen_flux_and_adiabatic_compression_follow_the_radius():
    coil = Coil("c", 0.5, 0.0, 1.0, 1.0e5, 0.01, "copper")
    schedule = PulseSchedule(
        time_s=np.array([0.0, 1.0e-3]),
        plasma_radius_m=np.array([0.04, 0.02]),
        coil_current_scale=np.array([1.0, 1.0]),
    )
    history = evolve(
        CoilSet((coil,)),
        schedule,
        ion_density_m3=1.0e20,
        ion_temperature_kev=1.0,
        deuterium_fraction=0.5,
        flux_conservation="frozen_axial_flux",
    )
    start, stop = history.points
    assert stop.ion_density_m3 == pytest.approx(start.ion_density_m3 * 8.0)
    assert stop.ion_temperature_kev == pytest.approx(start.ion_temperature_kev * 4.0)
    assert stop.internal_field_t == pytest.approx(start.internal_field_t * 4.0)
    assert stop.fusion_power_w > 0.0
    assert math.isfinite(stop.beta)


def test_closed_ramp_returns_no_magnetic_energy_when_recovery_is_perfect():
    time = np.array([0.0, 1.0, 2.0])
    magnetic = np.array([10.0, 40.0, 10.0])
    ohmic = np.zeros(3)
    perfect = supply_energy(time, magnetic, ohmic, recovery_fraction=1.0)
    dumped = supply_energy(time, magnetic, ohmic, recovery_fraction=0.0)
    assert perfect.drawn_j == pytest.approx(0.0)
    assert dumped.drawn_j == pytest.approx(30.0)
    assert dumped.magnetic_energy_change_j == pytest.approx(0.0)


def test_carnot_channel_cannot_exceed_the_carnot_limit():
    assert carnot_efficiency(600.0, 300.0) == pytest.approx(0.5)
    assert carnot_efficiency(300.0, 300.0) == 0.0
    supply = supply_energy(np.array([0.0, 1.0]), np.array([0.0, 0.0]), np.array([0.0, 0.0]), 1.0)
    ledger = convert(
        fusion_j=100.0,
        alpha_fraction=0.2,
        neutron_deposited_fraction=1.0,
        alpha_escape_fraction=0.0,
        alpha_wall_fraction=0.0,
        alpha_direct_efficiency=0.5,
        shell_channel="carnot_fraction",
        carnot_fraction=1.0,
        hot_k=600.0,
        cold_k=300.0,
        supply=supply,
    )
    assert ledger.electric_shell_j == pytest.approx(0.5 * ledger.shell_heat_j)
    assert ledger.alpha_escaped_j + ledger.alpha_wall_j + ledger.alpha_retained_j == pytest.approx(
        ledger.alpha_j
    )


def test_per_coil_scale_freezes_the_first_knot_and_leaves_later_rows_free():
    coils = CoilSet(
        (
            Coil("inner", 0.4, 0.0, 1.0, 1.0e4, 0.01, "copper"),
            Coil("outer", 0.6, 0.2, 1.0, 1.0e4, 0.01, "copper"),
        )
    )
    times = np.array([0.0, 1.0e-3])
    radii = np.array([0.04, 0.02])
    both = PulseSchedule(times, radii, np.array([1.0, 1.0]))
    inner_only = PulseSchedule(
        times,
        radii,
        np.array([1.0, 1.0]),
        per_coil_scale=np.array([[1.0, 0.0], [1.0, 0.0]]),
    )
    later = PulseSchedule(
        times,
        radii,
        np.array([1.0, 1.0]),
        per_coil_scale=np.array([[1.0, 0.0], [1.0, 1.0]]),
    )
    kwargs = {
        "ion_density_m3": 1.0e20,
        "ion_temperature_kev": 1.0,
        "deuterium_fraction": 0.5,
        "flux_conservation": "frozen_axial_flux",
    }
    full = evolve(coils, both, **kwargs)
    staged = evolve(coils, inner_only, **kwargs)
    moved = evolve(coils, later, **kwargs)
    inner_field = coils.on_axis_field(0.0, np.array([1.0e4, 0.0]))
    assert staged.points[0].vacuum_field_t == pytest.approx(inner_field)
    assert staged.points[0].vacuum_field_t < full.points[0].vacuum_field_t
    assert staged.points[1].internal_field_t == pytest.approx(staged.points[0].internal_field_t * 4.0)
    assert moved.points[1].internal_field_t == pytest.approx(staged.points[1].internal_field_t)
    assert moved.points[1].vacuum_field_t > moved.points[0].vacuum_field_t


def test_injected_shared_scale_fills_coils_that_have_no_own_scale():
    slots = {
        "compression.time_s.0": 0.0,
        "compression.time_s.1": 1.0,
        "compression.plasma_radius_m.0": 0.05,
        "compression.plasma_radius_m.1": 0.02,
        "compression.coil_current_scale.0": 0.2,
        "compression.coil_current_scale.1": 0.2,
        "compression.coils.a.scale.0": 1.0,
        "compression.coils.a.scale.1": 1.0,
        "compression.coils.b.scale.0": 0.0,
        "compression.coils.b.scale.1": 0.5,
    }
    provenance = {key: "config" for key in slots}
    provenance["compression.coil_current_scale.0"] = "mapping"
    provenance["compression.coil_current_scale.1"] = "mapping"
    provenance["compression.coils.b.scale.0"] = "mapping"
    provenance["compression.coils.b.scale.1"] = "mapping"
    mixed = build_schedule(slots, coil_names=["a", "b"], provenance=provenance)
    row0, row1 = mixed.current_scales(2)
    assert row0[0] == pytest.approx(0.2)
    assert row0[1] == pytest.approx(0.0)
    assert row1[1] == pytest.approx(0.5)
