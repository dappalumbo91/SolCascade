"""The remaining Carnot gain outranks the published steam step."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_lever import (
    build_lever,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_lever()


def test_gains_keep_the_named_steps(report):
    assert report.demonstrated == 0.33
    assert report.built == 0.36
    assert report.catalog == 0.45
    assert report.discharge_gwd == 45.0
    assert report.higher_burnup_gwd == 60.0
    assert report.atomic_percent == 6.5
    assert report.fabrication_cap_percent == 5.0
    assert report.reload_tonnes == 24.3
    assert report.demonstrated < report.built < report.carnot < report.catalog
    assert f"{report.carnot:.6g}" == "0.446429"
    assert f"{report.built_step:.6g}" == "1.09091"
    assert f"{report.built_carnot_step:.6g}" == "1.24008"
    assert f"{report.published_gain:.6g}" == "0.0909091"
    assert f"{report.remaining_gain:.6g}" == "0.240079"
    assert f"{report.gain_ratio:.6g}" == "2.64087"
    assert f"{report.burnup_step:.6g}" == "1.33333"
    assert f"{report.burnup_gain:.6g}" == "0.333333"
    assert f"{report.burnup_ratio:.6g}" == "3.66667"
    assert report.published_gain < report.remaining_gain < report.burnup_gain
    assert f"{report.companion:.6g}" == "0.988644"
    assert f"{report.built_companion:.6g}" == "1.07852"
    assert report.companion < 1.0 < report.built_companion
    assert abs(report.published_gain - (report.built_step - 1.0)) / report.published_gain < 1.0e-15
    assert abs(report.remaining_gain - (report.built_carnot_step - 1.0)) / report.remaining_gain < 1.0e-15
    assert abs(report.gain_ratio - report.remaining_gain / report.published_gain) / report.gain_ratio < 1.0e-15
    assert abs(report.burnup_ratio - report.burnup_gain / report.published_gain) / report.burnup_ratio < 1.0e-15
    assert abs(report.built_companion - report.companion * report.built_step) / report.built_companion < 1.0e-15
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_ranks_the_gains(report):
    text = render()
    assert text == render()
    assert "The steam share leaves a larger electricity gain above the plants-being-built efficiency." in text
    assert "The chain stays at k_eff of 1." in text
    assert "The demonstrated efficiency is `0.33`." in text
    assert "The plants-being-built efficiency is `0.36`." in text
    assert "The temperature Carnot is `0.446429`." in text
    assert "The named built step stays `1.09091`." in text
    assert "The published gain is `0.0909091`." in text
    assert "The same quotient is `1/11`." in text
    assert "The published gain is one less than the named built step." in text
    assert "The published gain keeps the built score." in text
    assert "The named Carnot step on the built score stays `1.24008`." in text
    assert "The remaining gain is `0.240079`." in text
    assert "The same quotient is `121/504`." in text
    assert "The remaining gain is one less than the named Carnot step on the built score." in text
    assert "The remaining gain keeps the built score." in text
    assert "The remaining gain over the published gain is `2.64087`." in text
    assert "The same quotient is `1331/504`." in text
    assert "That ratio keeps the two named steps." in text
    assert "The published gain stays under the remaining gain." in text
    assert "The temperature Carnot stays the bound." in text
    assert "The 60 GWd/t step over 45 GWd/t is `1.33333`." in text
    assert "The same quotient is `4/3`." in text
    assert "The burnup step stays the stack score." in text
    assert "The burnup gain is `0.333333`." in text
    assert "The same quotient is `1/3`." in text
    assert "The burnup gain keeps the stack score." in text
    assert "The burnup gain over the published gain is `3.66667`." in text
    assert "The same quotient is `11/3`." in text
    assert "That ratio keeps the burnup step and the named built step." in text
    assert "The remaining gain stays under the burnup gain." in text
    assert "The atomic percent stays `6.5`." in text
    assert "The atomic percent stays above the fabrication cap." in text
    assert "The 60 GWd/t stack stays unrun." in text
    assert "The fabrication cap stays 5 percent." in text
    assert "The filed companion is `0.988644`." in text
    assert "The filed companion stays the reload score." in text
    assert "The fuel-side year at the demonstrated efficiency stays under the nameplate year." in text
    assert "The same reload at the plants-being-built efficiency stands at `1.07852` of the nameplate year." in text
    assert "That quotient is the filed companion times the named built step." in text
    assert "That quotient keeps the reload score and the built score." in text
    assert "The fuel-side year at the plants-being-built efficiency stands above the nameplate year." in text
    assert "The filed reload stays `24.3` t." in text
    assert "The demonstrated plant stays at its filed efficiency." in text
    assert "The catalog ceiling stays unrun." in text
    assert "The fast row's filed barrier is the sodium system." in text
    assert "This score adds no joules." in text
    assert "The published stamp records these gains as rationals of the named steps." in text
    assert "This run does not install a seed as a new efficiency or a new plant row." in text
    assert "The checked-in schedule is not rewritten." in text
    printed = {
        f"{value:.6g}"
        for value in (
            report.demonstrated,
            report.built,
            report.carnot,
            report.catalog,
            report.built_step,
            report.built_carnot_step,
            report.published_gain,
            report.remaining_gain,
            report.gain_ratio,
            report.burnup_step,
            report.burnup_gain,
            report.burnup_ratio,
            report.atomic_percent,
            report.companion,
            report.built_companion,
            report.reload_tonnes,
            report.nuclear_scalar,
            report.thermodynamics_scalar,
            report.gamma_fsot,
            report.valve_beta,
        )
    }
    for banned in (
        report.nuclear_scalar * report.gain_ratio,
        report.thermodynamics_scalar * report.remaining_gain,
        report.gamma_fsot * report.built_companion,
        report.valve_beta * report.burnup_ratio,
        report.nuclear_scalar * report.published_gain,
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
    share = readme.index("## Fission share")
    lever = readme.index("## Fission lever")
    tests = readme.index("## Tests")
    assert share < lever < tests
    assert "py -3.14 -m solcascade.production.fission_lever" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_lever.py"
    ).read_text(encoding="utf-8")
    assert "1-alpha" not in source
    assert "pi^-3*phi^5" not in source
    assert "warp_shell_energy_j" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9"):
        assert needle not in source
        assert needle not in text
    observables = (
        project_root() / "src" / "solcascade" / "production" / "fission_observables.py"
    ).read_text(encoding="utf-8")
    assert "The channel search stays closed." in observables
    package = (project_root() / "src" / "solcascade" / "production" / "__init__.py").read_text(
        encoding="utf-8"
    )
    assert "fission_lever" not in package
    assert "import" not in package
    lean = (project_root() / "lean" / "FissionStamp.lean").read_text(encoding="utf-8")
    assert "installedClosings.length = 79" in lean
    assert "theorem published_gain" in lean
    assert "theorem remaining_gain" in lean
    assert "theorem gain_ratio" in lean
    assert "theorem burnup_step" in lean
    assert "theorem burnup_gain" in lean
    assert "theorem burnup_ratio" in lean
    assert "theorem gain_order" in lean
    assert "fission_lever" not in lean
    log = (project_root() / "data" / "fission_stamp.txt").read_text(encoding="utf-8").splitlines()
    assert log[-4:] == ["closes 79", "pairs 4", "miss 1", "fails 0"]
