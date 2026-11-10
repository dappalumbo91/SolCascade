"""Cooling load, decay-heat removal, and the published small class."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_safety import (
    BUILT_EFFICIENCY,
    SMALL_CLASS_MWE,
    build_safety,
    default_output_path,
    heat_dump_mj_per_kwh,
    main,
    render,
)
from solcascade.production.refine import project_root


def test_heat_dump_follows_the_page_identity():
    demonstrated = heat_dump_mj_per_kwh(0.33)
    built = heat_dump_mj_per_kwh(BUILT_EFFICIENCY)
    assert demonstrated == pytest.approx((1.0 / 0.33 - 1.0) * 3.6)
    assert built == pytest.approx(6.4)
    assert abs(demonstrated - 7.3) < 0.02
    with pytest.raises(RuntimeError):
        heat_dump_mj_per_kwh(0.0)


@pytest.fixture(scope="module")
def report():
    return build_safety()


def test_small_class_scales_the_nameplate(report):
    assert report.efficiency == pytest.approx(0.33)
    assert report.built_efficiency == pytest.approx(0.36)
    assert report.carnot == pytest.approx(1.0 - 310.0 / 560.0)
    assert report.dump_carnot_mj < report.dump_built_mj < report.dump_demonstrated_mj
    assert report.dump_carnot_mj == pytest.approx(4.464)
    assert report.size_ratio == pytest.approx(SMALL_CLASS_MWE / 1000.0)
    assert report.large_volume_m3 == pytest.approx(1000.0e6 / 0.33 / 1.0e8)
    assert report.small_volume_m3 == pytest.approx(report.large_volume_m3 * report.size_ratio)
    assert report.linear_scale == pytest.approx(report.size_ratio ** (1.0 / 3.0))
    assert report.small_nameplate_j == pytest.approx(report.large_nameplate_j * report.size_ratio)
    assert report.power_density_w_m3 == 1.0e8
    assert report.fabrication_cap_percent == 5.0
    assert report.size_seed.name == "G/pi"
    assert report.built_seed.name == "gamma/phi"
    assert min(report.size_seed.relative_gap, report.built_seed.relative_gap) > 1.0e-3
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX


def test_report_keeps_decay_heat_off_the_ledger(report):
    text = render()
    assert text == render()
    assert "Routine cooling for a steam plant is set by its thermal efficiency." in text
    assert "The heat dumped per kilowatt-hour falls as the efficiency rises." in text
    assert "The cooling page states no decay-heat fraction, so no seed is ranked against it." in text
    assert "Emergency core cooling is the page's system when primary cooling fails." in text
    assert "If circulation is maintained, the page says the decay heat ceases to be a problem after a day or two." in text
    assert "Decay heat is a cooling load." in text
    assert "It is not a second electricity product." in text
    assert "The delayed-neutron margin stays in the core." in text
    assert "The core volume scales with the thermal power." in text
    assert "The joules per kilogram of fuel stay on the discharge." in text
    assert "The fast row's filed barrier is the sodium system." in text
    assert "The small light-water class does not retire that barrier." in text
    assert "No seed is installed as a cooling figure or as the size ratio." in text
    assert "The checked-in schedule is not rewritten." in text
    for scaled in (
        report.nuclear_scalar * report.large_nameplate_j,
        report.gamma_fsot * report.small_volume_m3,
        report.size_seed.value * report.large_nameplate_j,
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
    process = readme.index("## Fission process")
    safety = readme.index("## Fission safety")
    tests = readme.index("## Tests")
    assert process < safety < tests
    assert "py -3.14 -m solcascade.production.fission_safety" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_safety.py"
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
    assert "fission_safety" not in package
