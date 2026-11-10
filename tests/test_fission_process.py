"""Light-water process figures, from the catalog and the association's own balance."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_process import (
    ANNUAL_SWU,
    HIGHER_BURNUP_GWD,
    RELOAD_TONNES,
    build_process,
    default_output_path,
    enrichment_estimate,
    main,
    render,
)
from solcascade.production.refine import catalog_record, project_root


def test_enrichment_basis_matches_the_filed_joules():
    record = catalog_record("fission-thermal")
    swu, kwh = enrichment_estimate(str(record["input_energy_basis"]))
    assert swu == 200.0
    assert kwh == 50.0
    assert swu * kwh * 3_600_000.0 == pytest.approx(float(record["input_energy_j"]))
    with pytest.raises(RuntimeError):
        enrichment_estimate("no separative work here")


@pytest.fixture(scope="module")
def report():
    return build_process()


def test_heavy_metal_steps_stay_below_the_ranked_kilogram(report):
    assert report.discharge_gwd == 45.0
    assert report.efficiency == pytest.approx(0.33)
    assert report.carnot == pytest.approx(1.0 - 310.0 / 560.0)
    assert report.heavy_electric_j == pytest.approx(report.discharge_j_per_kg * report.efficiency)
    assert report.heavy_carnot_j == pytest.approx(report.discharge_j_per_kg * report.carnot)
    assert report.burnup_step == pytest.approx(HIGHER_BURNUP_GWD / 45.0)
    assert report.higher_electric_j == pytest.approx(report.heavy_electric_j * report.burnup_step)
    assert report.higher_carnot_j == pytest.approx(report.heavy_carnot_j * report.burnup_step)
    assert report.stack_factor == pytest.approx(report.burnup_step * report.carnot_step)
    assert report.higher_over_ranked == pytest.approx(report.higher_discharge_j / report.ranked_density_j)
    assert report.higher_over_ranked < 0.1
    assert report.stack_over_ranked_electric < 0.1
    assert report.swu_per_kg_reload == pytest.approx(ANNUAL_SWU / (RELOAD_TONNES * 1000.0))
    assert report.heu_swu_over_reload == pytest.approx(200.0 / report.swu_per_kg_reload)
    assert report.reload_output_over_input == pytest.approx(
        report.heavy_electric_j / report.enrichment_j_per_kg_reload
    )
    assert report.nameplate_j != pytest.approx(report.fuel_side_j, rel=1.0e-3)
    assert report.fuel_side_j == pytest.approx(RELOAD_TONNES * 1000.0 * report.heavy_electric_j)
    assert report.power_density_w_m3 == 1.0e8
    assert report.discharge_seed.name == "gamma/pi"
    assert report.enrichment_seed.name == "gamma/pi"
    assert report.atomic_seed.name == "gamma/pi"
    assert report.discharge_seed.relative_gap > 1.0
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX


def test_report_keeps_the_stack_off_the_ledger(report):
    text = render()
    assert text == render()
    assert "That power density does not set the joules per kilogram." in text
    assert "A heavy-metal kilogram, a highly enriched kilogram, and a pure U-235 kilogram are three loads." in text
    assert "Those two annual figures are the page's round companions." in text
    assert "The higher burnup band stops at that limit." in text
    assert "That stack is not a filed plant." in text
    assert "The fast row remains the filed route for the actinides a thermal discharge leaves behind." in text
    assert "A domain scalar is not a plant efficiency." in text
    assert "The nuclear scalar stays out of the fission joules." in text
    assert "gamma_FSOT does not multiply the plant." in text
    assert "The engine valve stays out of the fission joules." in text
    assert "No seed is installed as a process fraction." in text
    assert "The checked-in schedule is not rewritten." in text
    for scaled in (
        report.nuclear_scalar * report.heavy_electric_j,
        report.nuclear_scalar * report.higher_carnot_j,
        report.gamma_fsot * report.higher_carnot_j,
        report.discharge_seed.value * report.ranked_density_j,
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
    ceiling = readme.index("## Fission ceiling")
    process = readme.index("## Fission process")
    tests = readme.index("## Tests")
    assert ceiling < process < tests
    assert "py -3.14 -m solcascade.production.fission_process" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_process.py"
    ).read_text(encoding="utf-8")
    assert "warp_shell_energy_j" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9"):
        assert needle not in source
    observables = (
        project_root() / "src" / "solcascade" / "production" / "fission_observables.py"
    ).read_text(encoding="utf-8")
    assert "The channel search stays closed." in observables
    package = (project_root() / "src" / "solcascade" / "production" / "__init__.py").read_text(
        encoding="utf-8"
    )
    assert "fission_process" not in package
