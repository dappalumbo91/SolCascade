"""Axial throat of the installed fields."""

import math

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.nature import (
    axial_peak,
    build_nature,
    default_output_path,
    dipole_axis_ratio,
    isotropic_loss_fraction,
    main,
    render,
)
from solcascade.production.refine import load_checked_in, project_root


def test_loss_cone_and_dipole_axis_are_geometric():
    assert isotropic_loss_fraction(1.0) == 1.0
    assert isotropic_loss_fraction(4.0) == pytest.approx(1.0 - math.sqrt(0.75))
    assert dipole_axis_ratio(1.0) == 1.0
    assert dipole_axis_ratio(2.0) == pytest.approx(0.125)
    with pytest.raises(ValueError):
        isotropic_loss_fraction(0.5)


@pytest.fixture(scope="module")
def report():
    return build_nature()


def test_compressed_plasma_sits_above_every_throat(report):
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX
    assert report.valve_beta < 1.0e-15
    assert 2.0 < report.gyroradii < 5.0
    assert report.bore_gyroradii > report.gyroradii
    assert report.four_coil.vacuum_rise > 1.0
    assert report.four_coil.vacuum_loss_fraction < 1.0
    assert report.four_coil.throat_over_internal < 1.0
    assert report.four_coil.internal_loss_fraction == 1.0
    assert report.four_coil.trapped_electric_j == 0.0
    assert report.helmholtz.vacuum_rise == 1.0
    assert report.helmholtz.throat_over_internal < 1.0
    assert report.helmholtz.internal_loss_fraction == 1.0
    assert report.surface_current.vacuum_rise == 1.0
    assert report.surface_current.throat_over_internal == 1.0
    assert report.surface_current.internal_loss_fraction == 1.0
    assert report.four_coil.axial_electric_j == pytest.approx(report.transit_electric_j)
    assert report.transit_electric_j < report.walk_electric_j
    assert report.bore_transit_electric_j < report.bore_walk_electric_j
    assert "Electromagnetism" in report.scalars
    assert "Geophysics" in report.scalars
    assert "gamma_FSOT" in report.scalars
    problem = load_checked_in()
    coarse = axial_peak(problem.coils, problem.coils.base_currents_a(), 0.2, samples=101)
    fine = axial_peak(
        problem.coils,
        problem.coils.base_currents_a() * problem.schedule.current_scales(4)[0],
        0.2,
        samples=2001,
    )
    assert fine[0] == pytest.approx(report.four_coil.vacuum_rise, rel=1.0e-6)
    assert coarse[0] == pytest.approx(fine[0], rel=1.0e-4)


def test_report_names_the_open_axis_and_writes_under_runs(report):
    text = render()
    assert text == render()
    assert "The gyro-walk remains the cross-field clock." in text
    assert "The compressed plasma sits above the throat." in text
    assert "An isotropic fuel leaves on the thermal transit." in text
    assert "Axial electricity versus the gyro-walk electricity: below." in text
    assert "Bore axial electricity versus the bore gyro-walk electricity: below." in text
    assert "A domain scalar is not a mirror ratio." in text
    assert "The engine valve stays out of the mirror ratio." in text
    assert "gamma_FSOT does not multiply the walk." in text
    assert "The checked-in schedule is not rewritten." in text
    reciprocal = 1.0 / report.valve_beta
    assert f"{reciprocal:.6e}" not in text
    assert f"{reciprocal:.6g}" not in text
    assert "3.328" not in text
    assert "3.703" not in text
    assert "warp_shell_energy_j" not in text
    output = default_output_path()
    assert output.parent.parent == project_root() / "runs"
    assert main() == 0
    assert output.read_text(encoding="utf-8") == text
    readme = (project_root() / "README.md").read_text(encoding="utf-8")
    surface = readme.index("## Surface field")
    nature = readme.index("## Nature")
    throat = readme.index("## Throat")
    plants = readme.index("## Working plants")
    observables = readme.index("## Fission observables")
    tests = readme.index("## Tests")
    assert surface < nature < throat < plants < observables < tests
    assert "py -3.14 -m solcascade.production.nature" in readme
    source = (project_root() / "src" / "solcascade" / "production" / "nature.py").read_text(
        encoding="utf-8"
    )
    assert "warp_shell_energy_j" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9"):
        assert needle not in source
