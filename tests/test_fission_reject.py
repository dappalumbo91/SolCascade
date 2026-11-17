"""The cooling page's extra heat keeps the named built step."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_reject import (
    build_reject,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_reject()


def test_extra_heat_keeps_the_named_step(report):
    assert report.demonstrated == 0.33
    assert report.built == 0.36
    assert report.catalog == 0.45
    assert report.demonstrated < report.built < report.carnot < report.catalog
    assert f"{report.carnot:.6g}" == "0.446429"
    assert f"{report.built_step:.6g}" == "1.09091"
    assert f"{report.published_gain:.6g}" == "0.0909091"
    assert f"{report.demonstrated_dump:.6g}" == "7.30909"
    assert f"{report.built_dump:.6g}" == "6.4"
    assert f"{report.dump_ratio:.6g}" == "1.14205"
    assert f"{report.extra_heat:.6g}" == "0.142045"
    assert f"{report.extra_over_gain:.6g}" == "1.5625"
    assert report.built_dump < report.demonstrated_dump
    assert report.extra_heat > report.published_gain
    assert abs(report.dump_ratio - report.demonstrated_dump / report.built_dump) / report.dump_ratio < 1.0e-15
    assert abs(report.extra_heat - (report.dump_ratio - 1.0)) / report.extra_heat < 1.0e-14
    assert abs(report.extra_over_gain - report.extra_heat / report.published_gain) / report.extra_over_gain < 1.0e-14
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_records_the_heat_gap(report):
    text = render()
    assert text == render()
    assert "The cooling a steam plant requires is set by its thermal efficiency." in text
    assert "The chain stays at k_eff of 1." in text
    assert "The demonstrated efficiency is `0.33`." in text
    assert "The plants-being-built efficiency is `0.36`." in text
    assert "The demonstrated dump stays `7.30909` MJ per kilowatt-hour." in text
    assert "The built dump stays `6.4` MJ per kilowatt-hour." in text
    assert "The demonstrated plant dumps `1.14205` times that built dump." in text
    assert "The same quotient is `201/176`." in text
    assert "That ratio keeps the named built step." in text
    assert "The extra heat is `0.142045`." in text
    assert "The same quotient is `25/176`." in text
    assert "The extra heat is one less than the heat ratio." in text
    assert "The extra heat keeps the named built step." in text
    assert "The page's about 14 percent stays that round." in text
    assert "The named built step stays `1.09091`." in text
    assert "The published gain stays `0.0909091`." in text
    assert "The extra heat is `1.5625` times the published gain." in text
    assert "The same quotient is `25/16`." in text
    assert "That multiple keeps the named built step and the published gain." in text
    assert "The extra heat stays above the published gain." in text
    assert "The temperature Carnot is `0.446429`." in text
    assert "The temperature Carnot stays the bound." in text
    assert "The catalog ceiling stays unrun." in text
    assert "The prose about two-thirds stays approximate." in text
    assert "The tower penalty stays the page's 2 to 5 percent range." in text
    assert "The fan share stays the page's about 1 percent." in text
    assert "The page's evaporative figure stays an upper bound." in text
    assert "The 60 GWd/t stack stays unrun." in text
    assert "The fabrication cap stays 5 percent." in text
    assert "The demonstrated plant stays at its filed efficiency." in text
    assert "The fast row's filed barrier is the sodium system." in text
    assert "This score adds no joules." in text
    assert "The published stamp records these quotients as rationals of the named steps." in text
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
            report.published_gain,
            report.demonstrated_dump,
            report.built_dump,
            report.dump_ratio,
            report.extra_heat,
            report.extra_over_gain,
            report.nuclear_scalar,
            report.thermodynamics_scalar,
            report.gamma_fsot,
            report.valve_beta,
        )
    }
    for banned in (
        report.nuclear_scalar * report.extra_heat,
        report.thermodynamics_scalar * report.dump_ratio,
        report.gamma_fsot * report.extra_over_gain,
        report.valve_beta * report.published_gain,
        report.nuclear_scalar * report.built_dump,
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
    column = readme.index("## Fission column")
    reject = readme.index("## Fission reject")
    tests = readme.index("## Tests")
    assert column < reject < tests
    assert "py -3.14 -m solcascade.production.fission_reject" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_reject.py"
    ).read_text(encoding="utf-8")
    assert "1-alpha" not in source
    assert "pi^-3*phi^5" not in source
    assert "warp_shell_energy_j" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9", "3.328", "3.703"):
        assert needle not in source
        assert needle not in text
    observables = (
        project_root() / "src" / "solcascade" / "production" / "fission_observables.py"
    ).read_text(encoding="utf-8")
    assert "The channel search stays closed." in observables
    assert "dialed=bindings_inside and quantum_emitted" in observables
    package = (project_root() / "src" / "solcascade" / "production" / "__init__.py").read_text(
        encoding="utf-8"
    )
    assert "fission_reject" not in package
    assert "import" not in package
    lean = (project_root() / "lean" / "FissionStamp.lean").read_text(encoding="utf-8")
    assert "installedClosings.length = 97" in lean
    assert "theorem reject_ratio" in lean
    assert "theorem reject_extra" in lean
    assert "theorem reject_over_gain" in lean
    assert "fission_reject" not in lean
    log = (project_root() / "data" / "fission_stamp.txt").read_text(encoding="utf-8").splitlines()
    assert log[-4:] == ["closes 97", "pairs 4", "miss 0", "fails 0"]
