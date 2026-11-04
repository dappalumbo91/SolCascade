"""Axisymmetric coil field, inductance, and the segmented cross-check."""

import math

import numpy as np

from solcascade.constants import VACUUM_PERMEABILITY_H_PER_M
from solcascade.geometry.coils import Coil, CoilSet, loop_field, loop_self_inductance_h


def _coil(**overrides) -> Coil:
    data = dict(
        name="a",
        radius_m=0.5,
        z_m=0.2,
        turns=1.0,
        current_a=1.0e5,
        wire_radius_m=0.01,
        material="copper",
    )
    data.update(overrides)
    return Coil(**data)


def test_on_axis_field_matches_the_closed_form():
    radius = 0.4
    current = 2.0e4
    axial = 0.15
    _radial, field = loop_field(0.0, axial, radius, current)
    expected = (
        VACUUM_PERMEABILITY_H_PER_M
        * current
        * radius ** 2
        / (2.0 * (radius ** 2 + axial ** 2) ** 1.5)
    )
    # radius*radius and radius**2 differ by one unit in the last place.
    assert math.isclose(field, expected, rel_tol=1.0e-12, abs_tol=0.0)
    near, near_field = loop_field(1.0e-8, axial, radius, current)
    assert abs(near) < 1.0e-8
    assert near_field == expected or math.isclose(near_field, expected, rel_tol=1.0e-6)


def test_segmented_biot_savart_matches_the_elliptic_field():
    coil = _coil(z_m=0.0, turns=4.0, current_a=3.0e4)
    coils = CoilSet((coil,))
    currents = coils.base_currents_a()
    radial = np.array([0.2])
    axial = np.array([0.15])
    elliptic_r, elliptic_z = coils.field(0.2, 0.15, currents)
    summed_r, summed_z = coils.field_grid(
        radial, axial, currents, method="biot_savart", segments=180
    )
    assert math.isclose(float(summed_z[0, 0]), elliptic_z, rel_tol=5.0e-3)
    assert math.isclose(float(summed_r[0, 0]), elliptic_r, rel_tol=5.0e-3, abs_tol=1.0e-6)


def test_inductance_matrix_is_symmetric_and_energy_is_positive():
    coils = CoilSet(
        (
            _coil(name="minus", z_m=-0.3),
            _coil(name="plus", z_m=0.3, turns=2.0),
        )
    )
    matrix = coils.inductance_matrix_h()
    assert np.allclose(matrix, matrix.T)
    assert matrix[0, 0] > 0.0
    assert matrix[0, 1] > 0.0
    single = loop_self_inductance_h(0.5, 0.01)
    assert math.isclose(matrix[0, 0], single, rel_tol=1.0e-12)
    assert math.isclose(matrix[1, 1], 4.0 * single, rel_tol=1.0e-12)
    energy = coils.magnetic_energy_j(coils.base_currents_a(), matrix)
    assert energy > 0.0
