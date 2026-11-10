"""The filed light-water plant under the conditions already scored."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_stress import (
    build_stress,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_stress()


def test_filed_conditions_hold(report):
    assert report.held_k == 1.0
    assert f"{report.carnot:.6g}" == "0.446429"
    assert f"{report.catalog_ceiling:.6g}" == "0.45"
    assert report.catalog_ceiling > report.carnot
    assert f"{report.companion_ratio:.6g}" == "0.988644"
    assert f"{report.discharge_over_ranked:.6g}" == "0.047357"
    assert f"{report.stack_factor:.6g}" == "1.80375"
    assert report.fabrication_cap_percent == 5.0
    assert f"{report.reload_tonnes:.6g}" == "24.3"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX
    large, small, built, carnot = report.power
    assert [case.name for case in report.power] == [
        "large nameplate",
        "small class",
        "plants being built",
        "temperature Carnot",
    ]
    assert large.electric_w == 1000.0e6
    assert small.electric_w == 300.0e6
    assert built.electric_w == carnot.electric_w == large.electric_w
    assert large.source == small.source == "FSOT heat ratio"
    assert built.source == "page efficiency"
    assert carnot.source == "page temperatures"
    assert built.gap is None and carnot.gap is None
    for case in report.power:
        assert case.waste_w == case.thermal_w - case.electric_w
        assert case.volume_m3 == case.thermal_w / 1.0e8
        assert case.mj_per_kwh == case.waste_w / case.electric_w * 3.6
        assert case.electric_w <= large.electric_w
    assert large.thermal_w == 1000.0e6 * (small.thermal_w / small.electric_w)
    assert f"{large.thermal_w:.6g}" == "3.0303e+09"
    assert f"{large.waste_w:.6g}" == "2.0303e+09"
    assert f"{large.volume_m3:.6g}" == "30.303"
    assert f"{large.mj_per_kwh:.6g}" == "7.30909"
    assert f"{large.gap:.6g}" == "8.10204e-08"
    assert f"{small.thermal_w:.6g}" == "9.09091e+08"
    assert f"{small.waste_w:.6g}" == "6.09091e+08"
    assert f"{small.volume_m3:.6g}" == "9.09091"
    assert f"{small.gap:.6g}" == "8.10204e-08"
    assert small.gap <= report.helium_mismatch
    assert f"{built.thermal_w:.6g}" == "2.77778e+09"
    assert f"{built.waste_w:.6g}" == "1.77778e+09"
    assert f"{built.volume_m3:.6g}" == "27.7778"
    assert f"{built.mj_per_kwh:.6g}" == "6.4"
    assert f"{carnot.thermal_w:.6g}" == "2.24e+09"
    assert f"{carnot.waste_w:.6g}" == "1.24e+09"
    assert f"{carnot.volume_m3:.6g}" == "22.4"
    assert f"{carnot.mj_per_kwh:.6g}" == "4.464"
    assert carnot.waste_w < built.waste_w < large.waste_w
    assert carnot.volume_m3 < built.volume_m3 < large.volume_m3
    expected = (
        ("large one-second fraction", "1.76426e+08", "0.0868964", "4.05394e-06"),
        ("small one-second fraction", "5.29278e+07", "0.0868964", "4.05394e-06"),
        ("large one-day fraction", "2.00482e+07", "0.0098745", "3.4255e-06"),
        ("small one-day fraction", "6.01447e+06", "0.0098745", "3.4255e-06"),
        ("large two-day fraction", "1.602e+07", "0.00789046", "4.64844e-07"),
        ("small two-day fraction", "4.80601e+06", "0.00789046", "4.64844e-07"),
    )
    assert len(report.shutdown) == len(expected)
    for case, row in zip(report.shutdown, expected, strict=True):
        name, cooling, share, gap = row
        assert case.name == name
        assert f"{case.cooling_w:.6g}" == cooling
        assert f"{case.share:.6g}" == share
        assert f"{case.gap:.6g}" == gap
        assert 0.0 < case.share < 1.0
        assert case.gap <= report.helium_mismatch
        assert case.cooling_w < case.waste_w
    large_shares = [case.share for case in report.shutdown if case.name.startswith("large ")]
    assert large_shares[0] > large_shares[1] > large_shares[2]


def test_report_keeps_the_chain_and_the_bounds(report):
    text = render()
    assert text == render()
    assert "This score runs the light-water plant already on file." in text
    assert "The chain stays at k_eff of 1." in text
    assert "This score leaves that ceiling unrun." in text
    assert "Fission is stopped in the shutdown cases." in text
    assert "This score assigns those systems no watt figure." in text
    assert "That stack is not a stressed case." in text
    assert "Decay heat is a cooling load." in text
    assert "It is not a second electricity product." in text
    assert "This run does not install a seed." in text
    assert "The checked-in schedule is not rewritten." in text
    assert "The fast row's filed barrier is the sodium system." in text
    printed = set()
    for case in report.power:
        for value in (case.electric_w, case.thermal_w, case.waste_w, case.volume_m3, case.mj_per_kwh):
            printed.add(f"{value:.6g}")
    for case in report.shutdown:
        for value in (case.cooling_w, case.waste_w, case.share):
            printed.add(f"{value:.6g}")
    for banned in (
        report.nuclear_scalar * report.power[0].thermal_w,
        report.thermodynamics_scalar * report.power[0].thermal_w,
        report.gamma_fsot * report.power[0].thermal_w,
        report.valve_beta * report.power[0].thermal_w,
        report.gamma_fsot * report.power[0].electric_w,
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
    heat = readme.index("## Fission decay heat")
    stress = readme.index("## Fission stress")
    tests = readme.index("## Tests")
    assert heat < stress < tests
    assert "py -3.14 -m solcascade.production.fission_stress" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_stress.py"
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
    assert "fission_stress" not in package
    assert "import" not in package
