"""Surface current in place of the distant windings."""

import pytest

from solcascade.production.refine import project_root
from solcascade.production.surface import (
    DIPOLE_EXTERIOR_OVER_INTERIOR,
    SURFACE_OVER_INTERIOR,
    build_surface,
    default_output_path,
    exterior_over_interior_quadrature,
    internal_energy_j,
    main,
    render,
)


@pytest.fixture(scope="module")
def report():
    return build_surface()


def test_dipole_exterior_is_half_the_interior_energy():
    ratio = exterior_over_interior_quadrature()
    assert ratio == pytest.approx(DIPOLE_EXTERIOR_OVER_INTERIOR, abs=1.0e-3)


def test_surface_current_stores_less_than_the_windings(report):
    plasma = report.plasma
    wall = report.wall
    bore = report.bore_plasma
    assert plasma.magnetic_j == pytest.approx(
        SURFACE_OVER_INTERIOR * internal_energy_j(plasma.radius_m, plasma.field_t)
    )
    assert wall.magnetic_j < report.helmholtz_stored_j
    assert plasma.magnetic_j < wall.magnetic_j
    assert bore.magnetic_j > plasma.magnetic_j
    assert plasma.pulsed_supply_j == pytest.approx((1.0 - report.recovery) * plasma.magnetic_j)
    assert report.beta1_complete_electric_j > plasma.pulsed_supply_j
    assert report.beta1_complete_electric_j > wall.pulsed_supply_j
    assert report.beta1_schedule_electric_j < plasma.pulsed_supply_j
    assert report.beta1_schedule_electric_j < wall.pulsed_supply_j
    assert report.bore_complete_electric_j > bore.pulsed_supply_j
    assert report.bore_schedule_electric_j < bore.pulsed_supply_j
    assert report.transit_time_s < report.beta1_walk_time_s
    assert report.transit_electric_j < report.beta1_walk_electric_j
    assert report.quadrature_exterior_over_interior == pytest.approx(0.5, abs=1.0e-3)


def test_report_names_the_surface_current_and_writes_under_runs():
    text = render()
    assert text == render()
    assert "A surface current on the plasma stores one and a half" not in text
    assert "The distant windings are the large store." in text
    assert "A persistent surface current draws no magnetic energy from the supply" in text
    assert "The thermal transit is shorter than the gyro-walk." in text
    assert "Complete burn versus the pulsed plasma-surface supply: above." in text
    assert "Schedule electricity versus the pulsed plasma-surface supply: below." in text
    assert "Shell-wall magnetic energy versus the Helmholtz store: below." in text
    assert "3.328" not in text
    assert "3.703" not in text
    assert "warp_shell_energy_j" not in text
    output = default_output_path()
    assert output.parent.parent == project_root() / "runs"
    assert main() == 0
    assert output.read_text(encoding="utf-8") == text
    readme = (project_root() / "README.md").read_text(encoding="utf-8")
    cost = readme.index("## Dominant cost")
    surface = readme.index("## Surface field")
    nature = readme.index("## Nature")
    throat = readme.index("## Throat")
    plants = readme.index("## Working plants")
    observables = readme.index("## Fission observables")
    tests = readme.index("## Tests")
    assert cost < surface < nature < throat < plants < observables < tests
    assert "py -3.14 -m solcascade.production.surface" in readme
    source = (project_root() / "src" / "solcascade" / "production" / "surface.py").read_text(encoding="utf-8")
    assert "warp_shell_energy_j" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9"):
        assert needle not in source
