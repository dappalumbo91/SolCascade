"""The filed 200 MeV ledger has a name, and no wave row carries it."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_quantum import (
    build_quantum,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_quantum()


def test_the_ledger_closes_off_the_wave_list(report):
    assert report.target == 200.0
    assert report.fragment_share == 0.85
    assert report.added_mev == 30.0
    assert len(report.wave_targets) == 8
    assert max(report.wave_targets[:3]) == 28.3
    assert report.anchor.name == "e^2*gamma^-6"
    assert f"{report.anchor.value:.6g}" == "199.784"
    assert f"{report.anchor.relative_gap:.6g}" == "0.0010809"
    assert report.anchor.review is False
    assert len(report.steps) == 2
    first, closing = report.steps
    assert first.dressing.name == "(e^2*gamma^-6)*(1+alpha*phi^-4)"
    assert f"{first.dressing.value:.6g}" == "199.997"
    assert f"{first.dressing.relative_gap:.6g}" == "1.73855e-05"
    assert f"{first.second_gap:.6g}" == "9.43823e-05"
    assert first.inside_helium == 0
    assert first.expressions == 116
    assert f"{first.median_spacing:.6g}" == "0.000511364"
    assert closing.dressing.name == "((e^2*gamma^-6)*(1+alpha*phi^-4))*(1+alpha*e^-6)"
    assert f"{closing.dressing.value:.6g}" == "200"
    assert f"{closing.dressing.relative_gap:.6g}" == "7.02539e-07"
    assert f"{closing.second_gap:.6g}" == "6.46013e-06"
    assert closing.inside_helium == 1
    assert closing.expressions == 116
    assert f"{closing.median_spacing:.6g}" == "0.000511364"
    assert closing.dressing.relative_gap <= report.helium_mismatch < closing.second_gap
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_names_the_file_quantum(report):
    text = render()
    assert text == render()
    assert "The pinned waves keep their own targets." in text
    assert "The chain stays at k_eff of 1." in text
    assert "The file value is `200` MeV." in text
    assert "No wave result carries that quantum." in text
    assert "The largest binding wave target is the helium-4 row." in text
    assert "The channel search stays closed." in text
    assert "The positive five-atom nearest is `e^2*gamma^-6` at `199.784`." in text
    assert "The anchor gap is `0.0010809`." in text
    assert "The dressing is `(e^2*gamma^-6)*(1+alpha*phi^-4)` at `199.997`." in text
    assert "The relative gap is `1.73855e-05`." in text
    assert "The dressing stays wider than the helium quotient match." in text
    assert "The dressing is `((e^2*gamma^-6)*(1+alpha*phi^-4))*(1+alpha*e^-6)` at `200`." in text
    assert "The relative gap is `7.02539e-07`." in text
    assert "The relative gap meets the helium quotient match." in text
    assert "One dressing of that anchor meets the helium quotient match." in text
    assert "The next dressing sits at relative gap `6.46013e-06`." in text
    assert text.count("The anchor has `116` dressings, median spacing `0.000511364`.") == 2
    assert "The closing value agrees with `200` at six figures." in text
    assert "The closing value names the file quantum." in text
    assert "This expression is not a wave row." in text
    assert "The catalog density stays on the 200 MeV ledger." in text
    assert "The about `30` MeV addendum stays the paths score." in text
    assert "The fragment share stays the round split of this ledger." in text
    assert "The 60 GWd/t stack stays unrun." in text
    assert "The catalog ceiling stays unrun." in text
    assert "The fast row's filed barrier is the sodium system." in text
    assert "This score adds no joules." in text
    assert "The published stamp records this close." in text
    assert "A domain scalar is not a megaelectronvolt." in text
    assert "This run does not install a seed as a new megaelectronvolt or a new plant row." in text
    assert "The checked-in schedule is not rewritten." in text
    printed = {
        f"{value:.6g}"
        for value in (
            report.target,
            report.anchor.value,
            report.ranked_density_j,
            report.fragment_share,
            report.nuclear_scalar,
            report.thermodynamics_scalar,
            report.gamma_fsot,
            report.valve_beta,
        )
    }
    for banned in (
        report.nuclear_scalar * report.target,
        report.thermodynamics_scalar * report.target,
        report.gamma_fsot * report.target,
        report.valve_beta * report.target,
        report.nuclear_scalar * report.ranked_density_j,
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
    reject = readme.index("## Fission reject")
    quantum = readme.index("## Fission quantum")
    tests = readme.index("## Tests")
    assert reject < quantum < tests
    assert "py -3.14 -m solcascade.production.fission_quantum" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_quantum.py"
    ).read_text(encoding="utf-8")
    assert "1-alpha" not in source
    assert "dialed" not in source
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
    assert "No wave result carries that quantum." in observables
    package = (project_root() / "src" / "solcascade" / "production" / "__init__.py").read_text(
        encoding="utf-8"
    )
    assert "fission_quantum" not in package
    assert "import" not in package
    lean = (project_root() / "lean" / "FissionStamp.lean").read_text(encoding="utf-8")
    assert "installedClosings.length = 97" in lean
    assert "((e^2*gamma^-6)*(1+alpha*phi^-4))*(1+alpha*e^-6)" in lean
    assert "fission_quantum" not in lean
    log = (project_root() / "data" / "fission_stamp.txt").read_text(encoding="utf-8").splitlines()
    assert log[-4:] == ["closes 97", "pairs 4", "miss 0", "fails 0"]
