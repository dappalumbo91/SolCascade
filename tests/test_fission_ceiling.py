"""Thermal ceiling tightened from the temperatures written in the catalog source."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_ceiling import (
    build_ceiling,
    default_output_path,
    main,
    render,
    temperature_carnot,
)
from solcascade.production.refine import catalog_record, project_root


def test_temperature_carnot_uses_the_stated_pair():
    hot, cold, carnot = temperature_carnot("near 560 K rejected near 310 K")
    assert hot == 560.0
    assert cold == 310.0
    assert carnot == pytest.approx(1.0 - 310.0 / 560.0)
    with pytest.raises(RuntimeError):
        temperature_carnot("no temperatures here")
    fast = str(catalog_record("fission-fast")["efficiency_ceiling_source"])
    with pytest.raises(RuntimeError):
        temperature_carnot(fast)


@pytest.fixture(scope="module")
def report():
    return build_ceiling()


def test_refined_bound_sits_between_the_plant_and_the_file(report):
    thermal = catalog_record("fission-thermal")
    assert report.hot_k == 560.0
    assert report.cold_k == 310.0
    assert report.demonstrated == pytest.approx(float(thermal["demonstrated_efficiency"]))
    assert report.catalog_ceiling == pytest.approx(float(thermal["efficiency_ceiling"]))
    assert report.carnot == pytest.approx(1.0 - 310.0 / 560.0)
    assert report.demonstrated < report.carnot < report.catalog_ceiling
    assert report.demonstrated_j == pytest.approx(report.density_j * report.demonstrated)
    assert report.carnot_j == pytest.approx(report.density_j * report.carnot)
    assert report.catalog_j == pytest.approx(report.density_j * report.catalog_ceiling)
    assert report.carnot_step_j == pytest.approx(report.carnot_j - report.demonstrated_j)
    assert report.catalog_margin_j == pytest.approx(report.catalog_j - report.carnot_j)
    assert report.catalog_margin_j < report.carnot_step_j
    assert report.demonstrated_seed.name == "gamma^2"
    assert report.demonstrated_seed.relative_gap == pytest.approx(0.009630072, rel=1.0e-4)
    assert report.carnot_seed.name == "phi/pi"
    assert report.fast_seed.name == "phi/pi"
    assert report.fast_ceiling == pytest.approx(0.5)
    assert min(
        report.demonstrated_seed.relative_gap,
        report.carnot_seed.relative_gap,
        report.fast_seed.relative_gap,
    ) > 1.0e-3
    assert report.nuclear_scalar != pytest.approx(report.demonstrated)
    assert report.nuclear_scalar != pytest.approx(report.carnot)
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX


def test_report_keeps_the_seed_out_of_the_ceiling(report):
    text = render()
    assert text == render()
    assert "The catalog ceiling sits above that temperature Carnot." in text
    assert "The refined bound is the temperature Carnot." in text
    assert "That extra is the round-file margin." in text
    assert "The refined electricity is the temperature Carnot." in text
    assert "The fast ceiling source states no temperature pair." in text
    assert "The fast row is not tightened." in text
    assert "A domain scalar is not a plant efficiency." in text
    assert "The nuclear scalar stays out of the fission joules." in text
    assert "gamma_FSOT does not multiply the plant." in text
    assert "The engine valve stays out of the fission joules." in text
    assert "The nearest seed to the demonstrated efficiency is `gamma^2` at " in text
    assert "The nearest seed is not installed as the demonstrated efficiency." in text
    assert "The nearest seed to the temperature Carnot is `phi/pi` at " in text
    assert "The nearest seed is not installed as the temperature Carnot." in text
    assert "The nearest seed is not installed as the fast ceiling." in text
    assert "The closest of these seeds is still wider than one part in a thousand." in text
    assert "The checked-in schedule is not rewritten." in text
    for scaled in (
        report.nuclear_scalar * report.carnot_j,
        report.nuclear_scalar * report.density_j,
        report.gamma_fsot * report.carnot_j,
        report.demonstrated_seed.value * report.density_j,
        report.carnot_seed.value * report.density_j,
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
    balance = readme.index("## Fission balance")
    ceiling = readme.index("## Fission ceiling")
    tests = readme.index("## Tests")
    assert balance < ceiling < tests
    assert "py -3.14 -m solcascade.production.fission_ceiling" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_ceiling.py"
    ).read_text(encoding="utf-8")
    assert "warp_shell_energy_j" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9"):
        assert needle not in source
    package = (project_root() / "src" / "solcascade" / "production" / "__init__.py").read_text(
        encoding="utf-8"
    )
    assert "fission_ceiling" not in package
