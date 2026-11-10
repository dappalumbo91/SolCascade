"""The operating heat is the filed nameplate times the FSOT heat ratio."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_decay_heat import (
    build_decay_heat,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_decay_heat()


def test_operating_heat_meets_the_helium_match(report):
    assert report.anchor_name == "pi-gamma^4"
    assert f"{report.anchor_value:.6g}" == "3.03059"
    assert f"{report.anchor_gap:.6g}" == "9.30911e-05"
    added, closing = report.steps
    assert added.dressing.name == "(pi-gamma^4)*(1-alpha*pi^-4)"
    assert f"{added.dressing.relative_gap:.6g}" == "1.81697e-05"
    assert added.inside_helium == 0
    assert added.inside_review == 30
    assert added.inside_tenth == 6
    assert closing.dressing.name == "((pi-gamma^4)*(1-alpha*pi^-4))*(1-alpha*e^-6)"
    assert f"{closing.dressing.relative_gap:.6g}" == "8.10204e-08"
    assert closing.inside_helium == 1
    assert closing.inside_review == 29
    assert closing.inside_tenth == 10
    assert closing.size.expressions == closing.size.inside == 116
    assert f"{closing.size.median_spacing:.6g}" == "0.000511364"
    assert closing.size.median_spacing > report.helium_mismatch
    assert f"{report.thermal_w:.6g}" == "3.0303e+09"
    assert f"{report.filed_thermal_w:.6g}" == "3.0303e+09"
    assert f"{report.small_thermal_w:.6g}" == "9.09091e+08"
    assert f"{report.filed_small_thermal_w:.6g}" == "9.09091e+08"
    assert f"{report.efficiency:.6g}" == "0.33"
    assert f"{report.filed_efficiency:.6g}" == "0.33"
    heat_gap = abs(report.thermal_w - report.filed_thermal_w) / report.filed_thermal_w
    assert f"{heat_gap:.6g}" == "8.10204e-08"
    assert heat_gap <= report.helium_mismatch
    assert report.thermal_w == 1000.0e6 * report.ratio
    assert report.small_thermal_w == 300.0e6 * report.ratio
    expected = (
        ("one-second fraction", "1.76426e+08", "1.76425e+08", "5.29278e+07", "5.29276e+07", "4.05394e-06"),
        ("one-day fraction", "2.00482e+07", "2.00482e+07", "6.01447e+06", "6.01445e+06", "3.4255e-06"),
        ("two-day fraction", "1.602e+07", "1.602e+07", "4.80601e+06", "4.806e+06", "4.64844e-07"),
    )
    for load, row in zip(report.loads, expected, strict=True):
        label, watts, fit_watts, small, fit_small, gap = row
        assert load.label == label
        assert f"{load.watts:.6g}" == watts
        assert f"{load.fit_watts:.6g}" == fit_watts
        assert f"{load.small_watts:.6g}" == small
        assert f"{load.fit_small_watts:.6g}" == fit_small
        assert f"{load.relative_gap:.6g}" == gap
        assert load.relative_gap <= report.helium_mismatch
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX


def test_report_uses_the_derived_heat_in_the_load(report):
    text = render()
    assert text == render()
    assert "The operating heat is the filed nameplate times the FSOT heat ratio." in text
    assert "The bare value is pi minus gamma to the fourth." in text
    assert "The derived operating heat is the heat in the cooling load." in text
    assert "Every scored gap meets the helium quotient match." in text
    assert "It is not a second electricity product." in text
    assert "The checked-in schedule is not rewritten." in text
    assert "((pi-gamma^4)*(1-alpha*pi^-4))*(1-alpha*e^-6)" in text
    printed = {
        f"{report.thermal_w:.6g}",
        f"{report.filed_thermal_w:.6g}",
        f"{report.small_thermal_w:.6g}",
        f"{report.filed_small_thermal_w:.6g}",
    }
    for load in report.loads:
        for value in (load.watts, load.fit_watts, load.small_watts, load.fit_small_watts):
            printed.add(f"{value:.6g}")
            assert f"{value:.6g}" in text
    for banned in (
        report.nuclear_scalar * report.thermal_w,
        report.thermodynamics_scalar * report.thermal_w,
        report.gamma_fsot * report.thermal_w,
        report.valve_beta * report.thermal_w,
    ):
        for formatted in (f"{banned:.6g}", f"{banned:.6e}"):
            if formatted in printed:
                continue
            assert formatted not in text
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
    load = readme.index("## Fission decay load")
    heat = readme.index("## Fission decay heat")
    tests = readme.index("## Tests")
    assert load < heat < tests
    assert "py -3.14 -m solcascade.production.fission_decay_heat" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_decay_heat.py"
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
    assert "fission_decay_heat" not in package
