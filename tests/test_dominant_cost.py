"""Dominant cost: the external coil store, and a larger beta = 1 sphere in the bore."""

import pytest

from solcascade.production.cost import BORE_FRACTION, build_cost, default_output_path, main, render
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_cost()


def test_unrecovered_store_is_the_helmholtz_cost_above_the_beta1_ceiling(report):
    assert report.baseline_peak_j > report.baseline_supply_j
    assert report.baseline_peak_j > report.fuel_ceiling_j
    assert report.fuel_ceiling_j > report.helmholtz_stored_j
    assert report.unrecovered_store_j > report.ohmic_ledger_j
    assert report.unrecovered_store_j > report.beta1_complete_electric_j
    assert report.beta1_complete_electric_j > report.ohmic_ledger_j
    assert report.store_over_plasma > 100.0
    assert report.plasma_magnetic_j < report.beta1_complete_electric_j
    assert report.beta1_walk_electric_j < report.beta1_schedule_electric_j < report.beta1_complete_electric_j


def test_larger_spheres_fit_the_bore_and_the_schedule_stays_under_the_store(report):
    match = report.store_match
    bore = report.bore
    assert match.initial_radius_m < report.shell_inner_m
    assert bore.initial_radius_m == pytest.approx(BORE_FRACTION * report.shell_inner_m)
    assert bore.initial_radius_m < report.shell_inner_m
    assert match.radius_factor > 1.0
    assert bore.radius_factor > match.radius_factor
    assert match.peak_beta == pytest.approx(1.0, rel=1.0e-6)
    assert bore.peak_beta == pytest.approx(1.0, rel=1.0e-6)
    assert match.field_edge_over_axis > 0.99
    assert bore.field_edge_over_axis > 0.95
    assert match.complete_electric_j == pytest.approx(report.practical_j, rel=1.0e-9)
    assert match.electric_j < report.practical_j
    assert match.walk_electric_j < match.electric_j
    assert bore.complete_electric_j > report.practical_j
    assert bore.electric_j < report.practical_j
    assert bore.walk_electric_j < bore.electric_j
    assert bore.walk_time_s < bore.efold_s
    assert bore.complete_fusion_j < report.jet_j
    assert bore.electric_power_w < bore.peak_power_w


def test_report_names_the_store_and_writes_under_runs():
    text = render()
    assert text == render()
    assert "The baseline coil store is the largest energy in this ranking." in text
    assert "The unrecovered Helmholtz store is larger than the" in text
    assert "ohmic ledger and larger than the beta = 1 complete-burn electricity." in text
    assert "The store-matching sphere fits inside the shell." in text
    assert "The bore-limited complete burn sits above the practical supply." in text
    assert "The bore-limited schedule stays under the practical supply." in text
    assert "3.328" not in text
    assert "3.703" not in text
    assert "warp_shell_energy_j" not in text
    output = default_output_path()
    assert output.parent.parent == project_root() / "runs"
    assert main() == 0
    assert output.read_text(encoding="utf-8") == text
    readme = (project_root() / "README.md").read_text(encoding="utf-8")
    failure = readme.index("## Failure trace")
    cost = readme.index("## Dominant cost")
    surface = readme.index("## Surface field")
    nature = readme.index("## Nature")
    throat = readme.index("## Throat")
    plants = readme.index("## Working plants")
    observables = readme.index("## Fission observables")
    tests = readme.index("## Tests")
    assert failure < cost < surface < nature < throat < plants < observables < tests
    assert "py -3.14 -m solcascade.production.cost" in readme
    source = (project_root() / "src" / "solcascade" / "production" / "cost.py").read_text(encoding="utf-8")
    assert "warp_shell_energy_j" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9"):
        assert needle not in source
