"""Operating fission plants, with FSOT recorded beside the catalog joules."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX, READINGS
from solcascade.production.plants import (
    build_plants,
    default_output_path,
    discharge_j_per_kg,
    main,
    render,
)
from solcascade.production.refine import catalog_record, project_root


def test_discharge_converts_the_cited_gigawatt_days():
    notes = str(catalog_record("fission-thermal")["notes"])
    gigawatt_days, per_kg = discharge_j_per_kg(notes)
    assert gigawatt_days == 45.0
    assert per_kg == pytest.approx(45.0 * 1.0e9 * 86400.0 / 1000.0)
    with pytest.raises(RuntimeError):
        discharge_j_per_kg("no discharge here")


@pytest.fixture(scope="module")
def report():
    return build_plants()


def test_working_plants_keep_the_catalog_joules(report):
    thermal = catalog_record("fission-thermal")
    fast = catalog_record("fission-fast")
    assert thermal["fsot_scalar"] is None
    assert fast["fsot_scalar"] is None
    assert report.thermal.trl == 9
    assert report.fast.trl == 8
    assert report.thermal.self_sustaining is True
    assert report.fast.continuous_external_drive is False
    assert report.thermal.delivered_j == pytest.approx(
        report.thermal.energy_density_j * report.thermal.efficiency
    )
    assert report.thermal.ceiling_delivered_j == pytest.approx(
        report.thermal.energy_density_j * report.thermal.ceiling
    )
    assert report.thermal.headroom == pytest.approx(0.45 / 0.33)
    assert report.fast.headroom == pytest.approx(0.5 / 0.4)
    assert report.efficiency_step == pytest.approx(0.4 / 0.33)
    assert report.discharge_j_per_kg == pytest.approx(3.888e12)
    assert report.discharge_over_ranked == pytest.approx(
        report.discharge_j_per_kg / report.thermal.energy_density_j
    )
    assert report.thermal.output_over_input == pytest.approx(
        report.thermal.delivered_j / report.thermal.input_j
    )
    assert report.burnup_demonstrated == 1.0
    assert report.burnup_cap == 1.0
    assert READINGS["fission-thermal"].physical_max == 1.0
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX
    assert report.nuclear_scalar == pytest.approx(report.thermodynamics_scalar)
    assert report.nuclear_scalar < 1.0
    assert report.nuclear_scalar != pytest.approx(report.thermal.efficiency)
    assert report.valve_beta < 1.0e-15
    assert report.thermal.delivered_j > report.thermal.input_j
    assert report.fast.delivered_j > report.fast.input_j


def test_report_leaves_the_scalar_out_of_the_joules(report):
    text = render()
    assert text == render()
    assert "The chain sustains once critical." in text
    assert "A domain scalar is not a plant efficiency." in text
    assert "The nuclear scalar stays out of the fission joules." in text
    assert "gamma_FSOT does not multiply the plant." in text
    assert "The engine valve stays out of the fission joules." in text
    assert "The tenfold window is empty." in text
    assert "A heavy-metal kilogram and a pure U-235 kilogram are different loads." in text
    assert "The checked-in schedule is not rewritten." in text
    scaled = report.nuclear_scalar * report.thermal.delivered_j
    assert f"{scaled:.6g}" not in text
    assert f"{scaled:.6e}" not in text
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
    throat = readme.index("## Throat")
    plants = readme.index("## Working plants")
    observables = readme.index("## Fission observables")
    tests = readme.index("## Tests")
    assert throat < plants < observables < tests
    assert "py -3.14 -m solcascade.production.plants" in readme
    source = (project_root() / "src" / "solcascade" / "production" / "plants.py").read_text(
        encoding="utf-8"
    )
    assert "warp_shell_energy_j" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9"):
        assert needle not in source
    package = (project_root() / "src" / "solcascade" / "production" / "__init__.py").read_text(
        encoding="utf-8"
    )
    assert "plants" not in package
