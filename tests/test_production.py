"""Production ledger: four topologies, the Bosch-Hale rescan, inventory, neutrons."""

import math

import pytest

from solcascade.constants import BOSCH_HALE_DT_T_MAX_KEV, BOSCH_HALE_DT_T_MIN_KEV, dt_reaction
from solcascade.production.refine import (
    SUPPLIES,
    catalog_record,
    default_output_path,
    electric_per_fusion_j,
    fuel_of,
    inventory_mass_kg,
    load_checked_in,
    main,
    neutron_account,
    pb11_burn_fraction,
    project_root,
    render,
    schedule_yield,
    score_topologies,
    temperature_scan,
)


@pytest.fixture(scope="module")
def problem():
    return load_checked_in()


@pytest.fixture(scope="module")
def fuel(problem):
    return fuel_of(problem)


@pytest.fixture(scope="module")
def scan(problem, fuel):
    return temperature_scan(problem, fuel)


def test_checked_in_schedule_matches_the_locked_fusion_joules(problem, fuel, scan):
    checked = scan[0]
    assert checked.temperature0_kev == pytest.approx(0.5)
    assert checked.fusion_j == pytest.approx(0.4977033387900156, rel=0.0, abs=1.0e-12)
    assert checked.burn_fraction == pytest.approx(checked.fusion_j / fuel.ceiling_j)
    assert fuel.ceiling_j == pytest.approx(73778.1142668, rel=1.0e-9)
    assert checked.walk_yield_j < checked.fusion_j


def test_baseline_store_stays_closed_and_same_b_draws_sit_under_the_ceiling(problem, fuel, scan):
    scores = score_topologies(scan[0].fusion_j, fuel, electric_per_fusion_j(problem))
    by_name = {row.name: row for row in scores}
    assert [name for name, _supply in SUPPLIES] == [row.name for row in scores]
    assert by_name["baseline"].closed
    assert by_name["baseline"].ceiling_over_supply < 1.0
    for name in ("same_b_flat", "same_b_optimal", "same_b_helmholtz"):
        assert not by_name[name].closed
        assert by_name[name].ceiling_over_supply > 1.0
        assert by_name[name].net_j < 0.0
    assert by_name["baseline"].net_j == pytest.approx(-by_name["baseline"].supply_j, rel=1.0e-6)
    reaction = dt_reaction()
    expected = (
        reaction.alpha_fraction
        * problem.config["conversion"]["alpha_escape_fraction"]
        * problem.config["conversion"]["alpha_direct_efficiency"]
    )
    assert electric_per_fusion_j(problem) == pytest.approx(expected)
    assert by_name["same_b_helmholtz"].supply_j / fuel.ceiling_j == pytest.approx(0.01121607631232924, rel=1.0e-9)


def test_hotter_schedules_stay_inside_bosch_hale_and_miss_the_helmholtz_draw(fuel, scan):
    assert scan[-1].peak_temperature_kev == pytest.approx(100.0)
    for row in scan:
        assert BOSCH_HALE_DT_T_MIN_KEV <= row.temperature0_kev
        assert row.peak_temperature_kev <= BOSCH_HALE_DT_T_MAX_KEV
        assert 0.0 < row.burn_fraction < 1.0
        assert row.walk_burn_fraction < row.burn_fraction
    assert scan[-1].fusion_j > scan[0].fusion_j
    helmholtz = dict(SUPPLIES)["same_b_helmholtz"]
    assert scan[-1].fusion_j < helmholtz
    assert scan[-1].burn_fraction < helmholtz / fuel.ceiling_j
    walk_20 = next(row for row in scan if math.isclose(row.peak_temperature_kev, 20.0, rel_tol=1.0e-6))
    assert scan[-1].walk_yield_j < walk_20.walk_yield_j
    assert scan[0].peak_beta > 1.0
    assert scan[-1].peak_beta > scan[0].peak_beta


def test_inventory_and_tritium_account_scale_with_reactions(problem, fuel, scan):
    hottest = max(scan, key=lambda row: row.fusion_j)
    jet_j = float(catalog_record("fusion-dt-magnetic")["energy_density_j"])
    fission_j = float(catalog_record("fission-fast")["energy_density_j"])
    jet_mass = inventory_mass_kg(jet_j, hottest, fuel)
    fission_mass = inventory_mass_kg(fission_j, hottest, fuel)
    assert jet_j == pytest.approx(5.9e7)
    assert jet_mass > fuel.mass_kg
    assert fission_mass > jet_mass
    checked_account = neutron_account(
        "checked-in", scan[0].fusion_j, fuel.mass_kg, scan[0].burn_fraction, fuel
    )
    reaction = dt_reaction()
    assert checked_account.neutrons == pytest.approx(scan[0].fusion_j / reaction.q_joule)
    assert checked_account.tritium_consumed_kg < checked_account.tritium_loaded_kg
    assert checked_account.neutron_energy_kev == pytest.approx(reaction.neutron_energy_kev)
    assert checked_account.neutron_energy_kev / 1000.0 == pytest.approx(14.048, rel=1.0e-3)
    hot_jet = neutron_account("jet", jet_j, jet_mass, hottest.burn_fraction, fuel)
    assert hot_jet.tritium_consumed_kg == pytest.approx(
        hot_jet.tritium_loaded_kg * hottest.burn_fraction, rel=1.0e-9
    )
    assert hot_jet.neutrons == pytest.approx(jet_j / reaction.q_joule)
    pb11 = pb11_burn_fraction(scan[0], problem)
    assert 0.0 < pb11 < 1.0e-15


def test_out_of_window_temperature_is_rejected(problem, fuel):
    with pytest.raises(ValueError):
        schedule_yield(problem, 0.1, fuel)
    radius0 = float(problem.schedule.plasma_radius_m[0])
    radius_peak = float(min(problem.schedule.plasma_radius_m))
    too_hot = (BOSCH_HALE_DT_T_MAX_KEV + 1.0) / (radius0 / radius_peak) ** 2
    with pytest.raises(ValueError):
        schedule_yield(problem, too_hot, fuel)


def test_report_states_the_closed_store_and_writes_under_runs(scan):
    text = render()
    assert text == render()
    assert "The baseline store stays closed on this fuel load." in text
    assert "The confinement time does not close the" in text
    assert "Helmholtz gap." in text
    assert "p-B11 burn fraction" in text
    assert "3.328" not in text
    assert "3.703" not in text
    assert "warp_shell_energy_j" not in text
    output = default_output_path()
    assert output.parent.parent == project_root() / "runs"
    assert main() == 0
    assert output.read_text(encoding="utf-8") == text
    readme = (project_root() / "README.md").read_text(encoding="utf-8")
    phase3 = readme.index("## Phase Three — Derived Values")
    standing = readme.index("## Production standing")
    refinement = readme.index("## Production refinement")
    failure = readme.index("## Failure trace")
    cost = readme.index("## Dominant cost")
    surface = readme.index("## Surface field")
    nature = readme.index("## Nature")
    throat = readme.index("## Throat")
    plants = readme.index("## Working plants")
    observables = readme.index("## Fission observables")
    tests = readme.index("## Tests")
    assert phase3 < standing < refinement < failure < cost < surface < nature < throat < plants < observables < tests
    assert "py -3.14 -m solcascade.production.refine" in readme
    assert "py -3.14 -m solcascade.production.trace" in readme
    assert scan[0].fusion_j == pytest.approx(0.4977033387900156)
