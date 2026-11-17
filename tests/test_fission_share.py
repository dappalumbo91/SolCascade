"""The plants-being-built efficiency is the named target under the Carnot."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_share import (
    build_share,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_share()


def test_share_keeps_the_named_steps(report):
    helium = report.helium_mismatch
    assert report.demonstrated == 0.33
    assert report.built == 0.36
    assert report.catalog == 0.45
    assert report.hot_k == 560.0
    assert report.cold_k == 310.0
    assert report.fast_efficiency == 0.4
    assert report.demonstrated < report.built < report.carnot < report.catalog
    assert f"{report.carnot:.6g}" == "0.446429"
    assert report.steam_name == "(((pi^-3*phi^5)*(1+alpha*G))*(1-alpha*e^-4))*(1-alpha*e^-5)"
    assert f"{report.steam_gap:.6g}" == "2.3554e-07"
    assert f"{report.steam_second_gap:.6g}" == "2.50889e-05"
    assert report.steam_second_gap > helium >= report.steam_gap
    assert f"{report.built_step:.6g}" == "1.09091"
    assert f"{report.steam_carnot_step:.6g}" == "1.35281"
    assert f"{report.built_carnot_step:.6g}" == "1.24008"
    assert f"{report.kept_built:.6g}" == "0.8064"
    assert f"{report.kept_demonstrated:.6g}" == "0.7392"
    assert f"{report.closed_share:.6g}" == "0.257669"
    assert f"{report.remaining_share:.6g}" == "0.742331"
    assert f"{report.unused_built:.6g}" == "0.1936"
    assert f"{report.unused_demonstrated:.6g}" == "0.2608"
    assert abs(report.kept_built - 1.0 / report.built_carnot_step) / report.kept_built < 1.0e-15
    assert abs(report.kept_demonstrated - 1.0 / report.steam_carnot_step) / report.kept_demonstrated < 1.0e-15
    assert abs(report.closed_share - (report.built_step - 1.0) / (report.steam_carnot_step - 1.0)) / report.closed_share < 1.0e-15
    assert abs(report.remaining_share - (1.0 - report.closed_share)) / report.remaining_share < 1.0e-15
    assert abs(report.built_step - report.built / report.demonstrated) / report.built_step < 1.0e-15
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_names_the_existing_target(report):
    text = render()
    assert text == render()
    assert "The target inside the open interval is the plants-being-built efficiency." in text
    assert "The steam score names that efficiency." in text
    assert "The chain stays at k_eff of 1." in text
    assert "The demonstrated efficiency is `0.33`." in text
    assert "The plants-being-built efficiency is `0.36`." in text
    assert "The page temperatures stay `560` K and `310` K." in text
    assert "The temperature Carnot is `0.446429`." in text
    assert "The catalog ceiling is `0.45`." in text
    assert "The steam closing is `(((pi^-3*phi^5)*(1+alpha*G))*(1-alpha*e^-4))*(1-alpha*e^-5)`." in text
    assert "The relative gap is `2.3554e-07`." in text
    assert "The next dressing sits at relative gap `2.50889e-05`." in text
    assert "The closing value agrees with `0.36` at six figures." in text
    assert "One expression names that efficiency." in text
    assert "The worked figure stays the plants-being-built efficiency." in text
    assert "The plants-being-built efficiency over the temperature Carnot is `0.8064`." in text
    assert "The same quotient is `504/625`." in text
    assert "The named Carnot step on the built score stays `1.24008`." in text
    assert "That quotient is the reciprocal of the named Carnot step on the built score." in text
    assert "That quotient keeps the built score." in text
    assert "The demonstrated efficiency over the temperature Carnot is `0.7392`." in text
    assert "The same quotient is `462/625`." in text
    assert "The named steam Carnot step stays `1.35281`." in text
    assert "That quotient is the reciprocal of the named steam Carnot step." in text
    assert "That quotient keeps the steam score." in text
    assert "The unused share at the plants-being-built efficiency is `0.1936`." in text
    assert "The same quotient is `121/625`." in text
    assert "The unused share is one minus the built share." in text
    assert "The unused share keeps the built share." in text
    assert "The unused share at the demonstrated efficiency is `0.2608`." in text
    assert "The same quotient is `163/625`." in text
    assert "The unused share is one minus the demonstrated share." in text
    assert "The unused share keeps the demonstrated share." in text
    assert "On the same heat, the electricity ratio is `1.09091`." in text
    assert "The electricity ratio is the named built step." in text
    assert "The named built step stays the built score." in text
    assert "The closed share of the open interval is `0.257669`." in text
    assert "The same quotient is `42/163`." in text
    assert "The closed share is one less than the built step, divided by one less than the steam Carnot step." in text
    assert "The closed share keeps those two steps." in text
    assert "The remaining share of the open interval is `0.742331`." in text
    assert "The same quotient is `121/163`." in text
    assert "The remaining share is one minus the closed share." in text
    assert "The remaining share keeps the closed share." in text
    assert "The temperature Carnot stays the bound." in text
    assert "The fast efficiency stays `0.4`." in text
    assert "The fast efficiency stays the bounds score." in text
    assert "The fast ceiling source states no temperature pair." in text
    assert "The fast row's filed barrier is the sodium system." in text
    assert "The catalog ceiling stays unrun." in text
    assert "The 60 GWd/t stack stays unrun." in text
    assert "The fabrication cap stays 5 percent." in text
    assert "This score adds no joules." in text
    assert "The published stamp records these quotients as rationals of the named steps." in text
    assert "This run does not install a seed as a new efficiency or a new plant row." in text
    assert "The checked-in schedule is not rewritten." in text
    printed = {
        f"{report.demonstrated:.6g}",
        f"{report.built:.6g}",
        f"{report.carnot:.6g}",
        f"{report.catalog:.6g}",
        f"{report.hot_k:.6g}",
        f"{report.cold_k:.6g}",
        f"{report.built_step:.6g}",
        f"{report.steam_carnot_step:.6g}",
        f"{report.built_carnot_step:.6g}",
        f"{report.kept_built:.6g}",
        f"{report.kept_demonstrated:.6g}",
        f"{report.closed_share:.6g}",
        f"{report.remaining_share:.6g}",
        f"{report.unused_built:.6g}",
        f"{report.unused_demonstrated:.6g}",
        f"{report.fast_efficiency:.6g}",
        f"{report.steam_gap:.6g}",
        f"{report.steam_second_gap:.6g}",
        f"{report.alpha:.6g}",
        f"{report.helium_mismatch:.6g}",
        f"{report.nuclear_scalar:.6g}",
        f"{report.thermodynamics_scalar:.6g}",
        f"{report.gamma_fsot:.6g}",
        f"{report.valve_beta:.6g}",
    }
    for banned in (
        report.nuclear_scalar * report.built,
        report.thermodynamics_scalar * report.carnot,
        report.gamma_fsot * report.kept_built,
        report.valve_beta * report.closed_share,
        report.nuclear_scalar * report.built_step,
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
    tails = readme.index("## Fission tails")
    share = readme.index("## Fission share")
    tests = readme.index("## Tests")
    assert tails < share < tests
    assert "py -3.14 -m solcascade.production.fission_share" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_share.py"
    ).read_text(encoding="utf-8")
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
    assert "fission_share" not in package
    assert "import" not in package
    lean = (project_root() / "lean" / "FissionStamp.lean").read_text(encoding="utf-8")
    assert "installedClosings.length = 79" in lean
    assert "theorem kept_built_share" in lean
    assert "theorem kept_demonstrated_share" in lean
    assert "theorem closed_share" in lean
    assert "theorem remaining_share" in lean
    assert "theorem unused_built_share" in lean
    assert "theorem unused_demonstrated_share" in lean
    assert "fission_share" not in lean
    log = (project_root() / "data" / "fission_stamp.txt").read_text(encoding="utf-8").splitlines()
    assert log[-4:] == ["closes 79", "pairs 4", "miss 1", "fails 0"]
