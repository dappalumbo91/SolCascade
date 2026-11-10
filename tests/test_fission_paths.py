"""Filed fission converters on one U-235 kilogram, with the page share beside them."""

import pytest

from solcascade.constants import FISSION_RECOVERABLE_ENERGY_MEV
from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_paths import (
    ADDED_MEV,
    FAST_ID,
    FRAGMENT_SHARE,
    THERMAL_ID,
    THERMIONIC_ID,
    build_paths,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import catalog_record, project_root


def test_page_addendum_agrees_with_the_ledger_share():
    assert FISSION_RECOVERABLE_ENERGY_MEV == 200.0
    assert ADDED_MEV == 30.0
    assert FRAGMENT_SHARE == (200.0 - 30.0) / 200.0
    assert FRAGMENT_SHARE == pytest.approx(0.85)


@pytest.fixture(scope="module")
def report():
    return build_paths()


def test_same_kilogram_keeps_the_catalog_joules(report):
    thermal = catalog_record(THERMAL_ID)
    thermionic = catalog_record(THERMIONIC_ID)
    fast = catalog_record(FAST_ID)
    assert thermal["fsot_scalar"] is None
    assert thermionic["fsot_scalar"] is None
    assert report.same_ids == (THERMAL_ID, THERMIONIC_ID)
    assert report.thermal.energy_density_j == report.thermionic.energy_density_j
    assert report.thermal.energy_density_j == pytest.approx(float(thermal["energy_density_j"]))
    assert float(fast["energy_density_j"]) != report.thermal.energy_density_j
    assert "Pu-239" in report.fast_reference_unit
    assert "physics-of-nuclear-energy" in report.governing_source
    assert report.thermal.trl == 9
    assert report.thermionic.trl == 7
    assert report.thermal.efficiency == pytest.approx(0.33)
    assert report.thermal.ceiling == pytest.approx(0.45)
    assert report.thermionic.efficiency == pytest.approx(0.05)
    assert report.thermionic.ceiling == pytest.approx(0.3)
    assert report.thermal.delivered_j == pytest.approx(
        report.thermal.energy_density_j * report.thermal.efficiency
    )
    assert report.thermionic.delivered_j == pytest.approx(
        report.thermionic.energy_density_j * report.thermionic.efficiency
    )
    assert report.thermionic.delivered_j < report.thermal.delivered_j
    assert report.thermionic.ceiling_delivered_j < report.thermal.delivered_j
    assert report.thermionic.ceiling < report.thermal.efficiency
    assert report.thermal.input_j == report.thermionic.input_j
    assert report.kinetic_inventory_j == pytest.approx(
        report.fragment_share * report.thermal.energy_density_j
    )
    assert report.kinetic_inventory_j < report.thermal.energy_density_j
    assert report.ledger_mev == 200.0
    assert report.added_mev == 30.0
    assert report.nearest.name == "G^2"
    assert report.nearest.value == pytest.approx(0.8389929697164259)
    assert report.nearest.relative_gap == pytest.approx(0.01294945, rel=1.0e-4)
    assert report.next_seed.name == "e/pi"
    assert report.next_seed.value == pytest.approx(0.8652559794322651)
    assert report.next_seed.relative_gap > report.nearest.relative_gap
    assert report.seeds_named == 45
    assert report.seeds_ranked == 19
    assert report.scalar_gap > report.nearest.relative_gap
    assert report.nuclear_scalar == pytest.approx(report.thermodynamics_scalar)
    assert report.nuclear_scalar != pytest.approx(report.fragment_share)
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX
    assert report.valve_beta < 1.0e-15


def test_report_keeps_the_share_out_of_the_joules(report):
    text = render()
    assert text == render()
    assert "The thermionic row is the other filed converter on this same U-235 kilogram." in text
    assert "The thermionic row delivers less electricity than the thermal plant." in text
    assert "Even at its file ceiling the thermionic row delivers less than the thermal plant." in text
    assert "The fast row stays on that unit." in text
    assert "The refinement inside the operating plants is the file efficiency ceiling." in text
    assert "That page states that about 85 percent of the released energy starts as fission-fragment kinetic energy." in text
    assert "In solid fuel the fragments stop in a microscopic distance and the energy becomes heat." in text
    assert "The catalog density stays on the 200 MeV ledger." in text
    assert "The fragment share is already inside that density." in text
    assert "The Rankine efficiency multiplies the whole ranked density." in text
    assert "The fragment inventory is heat in solid fuel." in text
    assert "The fragment inventory is not electricity." in text
    assert "The catalog files no collector efficiency for a fragment that leaves the fuel." in text
    assert "The operating electricity on this kilogram is the Rankine product." in text
    assert "A domain scalar is not a plant efficiency." in text
    assert "A domain scalar is not the fragment share." in text
    assert "The nuclear scalar stays out of the fission joules." in text
    assert "The domain scalar sits farther from the share than the nearest seed." in text
    assert "gamma_FSOT does not multiply the plant." in text
    assert "The engine valve stays out of the fission joules." in text
    assert "The nearest seed is `G^2` at " in text
    assert "The next seed is `e/pi` at " in text
    assert "The page figure is a round share." in text
    assert "The nearest seed is not installed as the share." in text
    assert "The checked-in schedule is not rewritten." in text
    partial = (
        report.fragment_share * report.thermal.efficiency * report.thermal.energy_density_j
    )
    assert f"{partial:.6g}" not in text
    assert f"{partial:.6e}" not in text
    for scaled in (
        report.nuclear_scalar * report.thermal.delivered_j,
        report.nuclear_scalar * report.kinetic_inventory_j,
        report.gamma_fsot * report.thermal.delivered_j,
        report.fragment_share * report.thermionic.delivered_j,
    ):
        assert f"{scaled:.6g}" not in text
        assert f"{scaled:.6e}" not in text
    reciprocal = 1.0 / report.valve_beta
    assert f"{reciprocal:.6e}" not in text
    assert f"{reciprocal:.6g}" not in text
    assert "3.328" not in text
    assert "3.703" not in text
    assert "warp_shell_energy_j" not in text
    assert "137.035999177" not in text
    output = default_output_path()
    assert output.parent.parent == project_root() / "runs"
    assert main() == 0
    assert output.read_text(encoding="utf-8") == text
    readme = (project_root() / "README.md").read_text(encoding="utf-8")
    plants = readme.index("## Working plants")
    observables = readme.index("## Fission observables")
    paths = readme.index("## Fission paths")
    tests = readme.index("## Tests")
    assert plants < observables < paths < tests
    assert "py -3.14 -m solcascade.production.fission_paths" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_paths.py"
    ).read_text(encoding="utf-8")
    assert "warp_shell_energy_j" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9"):
        assert needle not in source
    observables_source = (
        project_root() / "src" / "solcascade" / "production" / "fission_observables.py"
    ).read_text(encoding="utf-8")
    assert "The channel search stays closed." in observables_source
    package = (project_root() / "src" / "solcascade" / "production" / "__init__.py").read_text(
        encoding="utf-8"
    )
    assert "fission_paths" not in package
