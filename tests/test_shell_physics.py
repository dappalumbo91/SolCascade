"""Transport balance, heat, and stress limits."""

import math

import numpy as np
import pytest

from solcascade.geometry.shell import build_shell
from solcascade.materials.library import Material
from solcascade.mechanics.heat import solve_transient, steady_temperature
from solcascade.mechanics.stress import lame_stress, thermal_stress_free_sphere
from solcascade.transport.analytical import first_collision


def _material(name: str, conductivity: float, density: float = 1000.0, heat: float = 500.0) -> Material:
    return Material(
        name=name,
        density_kg_m3=density,
        thermal_conductivity_W_mK=conductivity,
        specific_heat_J_kgK=heat,
        youngs_modulus_Pa=1.0e11,
        poissons_ratio=0.3,
        thermal_expansion_per_K=1.0e-5,
        yield_strength_Pa=1.0e8,
        resistivity_ohm_m=None,
        resistivity_temperature_coeff_per_K=None,
        resistivity_reference_K=None,
        one_group_macroscopic_total_per_m=10.0,
        one_group_kerma_fraction=0.5,
        composition=(("X", 1.0),),
        note="",
    )


def test_first_collision_fractions_account_for_every_neutron():
    shell = build_shell(
        0.2,
        [
            {"name": "inner", "material": "mat", "thickness_m": 0.05},
            {"name": "outer", "material": "mat", "thickness_m": 0.05},
        ],
        radial_cells=8,
    )
    deposition = first_collision(shell, {"mat": _material("mat", 20.0)})
    total = float(deposition.energy_fraction.sum() + deposition.unresolved_fraction)
    assert total == pytest.approx(1.0)
    assert float(deposition.collision_fraction.sum()) + deposition.leakage_uncollided_fraction == pytest.approx(
        1.0
    )
    # Optical depth is 10 / m * 0.10 m = 1, so uncollided leakage is e^-1.
    assert deposition.leakage_uncollided_fraction == pytest.approx(math.exp(-1.0))


def test_short_pulse_with_no_conduction_is_the_adiabatic_rise():
    shell = build_shell(0.2, [{"name": "wall", "material": "mat", "thickness_m": 0.02}], 4)
    library = {"mat": _material("mat", conductivity=0.0, density=2000.0, heat=400.0)}
    source = np.full(len(shell.cells), 50.0)
    duration = 2.0
    result = solve_transient(
        shell,
        library,
        source,
        duration_s=duration,
        steps=5,
        initial_temperature_k=300.0,
        sink_temperature_k=300.0,
        heat_transfer_coefficient=0.0,
    )
    expected = 300.0 + source / np.array([cell.volume_m3 for cell in shell.cells]) * duration / (
        2000.0 * 400.0
    )
    assert result.temperature_k == pytest.approx(expected)


def test_zero_source_stays_at_the_sink_temperature():
    shell = build_shell(0.3, [{"name": "wall", "material": "mat", "thickness_m": 0.05}], 6)
    library = {"mat": _material("mat", conductivity=25.0)}
    source = np.zeros(len(shell.cells))
    steady = steady_temperature(shell, library, source, sink_temperature_k=450.0, heat_transfer_coefficient=800.0)
    assert steady == pytest.approx(np.full(len(shell.cells), 450.0))


def test_steady_shell_tracks_the_closed_form():
    shell = build_shell(0.2, [{"name": "wall", "material": "mat", "thickness_m": 0.1}], 40)
    library = {"mat": _material("mat", conductivity=20.0)}
    volume = np.array([cell.volume_m3 for cell in shell.cells])
    source_density = 1.0e6
    steady = steady_temperature(
        shell,
        library,
        source_density * volume,
        sink_temperature_k=300.0,
        heat_transfer_coefficient=2000.0,
    )
    analytical = _uniform_source_temperature(
        shell.face_radii_m(),
        source_density,
        conductivity=20.0,
        coefficient=2000.0,
        sink=300.0,
    )
    centers = np.array([cell.center_radius_m for cell in shell.cells])
    assert steady == pytest.approx(analytical(centers), rel=0.02, abs=1.0)
    assert float(np.min(steady)) > 300.0


def test_uniform_temperature_has_no_thermal_stress_and_lame_balances_pressure():
    radius = np.linspace(0.22, 0.28, 5)
    radial, hoop = thermal_stress_free_sphere(
        radius,
        np.full(radius.shape, 80.0),
        inner_m=0.2,
        outer_m=0.3,
        youngs_pa=2.0e11,
        poisson=0.3,
        expansion_per_k=1.0e-5,
    )
    assert radial == pytest.approx(0.0, abs=1.0)
    assert hoop == pytest.approx(0.0, abs=1.0)
    pressure = 2.0e6
    radial_p, _hoop_p = lame_stress(np.array([0.2, 0.25, 0.3]), 0.2, 0.3, pressure)
    assert radial_p[0] == pytest.approx(-pressure)
    assert radial_p[-1] == pytest.approx(0.0)


def _uniform_source_temperature(faces, source_density, conductivity, coefficient, sink):
    inner = float(faces[0])
    outer = float(faces[-1])
    constant = source_density * inner ** 3 / 3.0
    outward = source_density / 3.0 * (outer - inner ** 3 / outer ** 2)
    surface = sink + outward / coefficient
    offset = surface + source_density * outer ** 2 / (6.0 * conductivity) + constant / (
        conductivity * outer
    )

    def temperature(radius):
        radius = np.asarray(radius, dtype=float)
        return (
            -source_density * radius ** 2 / (6.0 * conductivity)
            - constant / (conductivity * radius)
            + offset
        )

    return temperature
