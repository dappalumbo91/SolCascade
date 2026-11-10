"""The 1952 long-weight coefficient."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_fit import build_fit, default_output_path, main, render
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_fit()


def test_long_weight_meets_helium(report):
    figure = report.weight
    assert figure.anchor.name == "pi^-2+G^3"
    assert f"{figure.anchor.value:.6g}" == "0.86981"
    assert f"{figure.anchor.relative_gap:.6g}" == "0.000218531"
    assert f"{figure.target:.6g}" == "0.87"
    assert len(figure.steps) == 2
    first, closing = figure.steps
    assert first.dressing.name == "(pi^-2+G^3)*(1+alpha*pi^-3)"
    assert f"{first.dressing.value:.6g}" == "0.870015"
    assert f"{first.dressing.relative_gap:.6g}" == "1.6768e-05"
    assert f"{first.second_gap:.6g}" == "5.13041e-05"
    assert (first.inside_review, first.inside_tenth, first.inside_helium) == (29, 3, 0)
    assert closing.dressing.name == "((pi^-2+G^3)*(1+alpha*pi^-3))*(1-alpha*e^-6)"
    assert f"{closing.dressing.value:.6g}" == "0.869999"
    assert f"{closing.dressing.relative_gap:.6g}" == "1.32067e-06"
    assert f"{closing.second_gap:.6g}" == "7.07846e-06"
    assert (closing.inside_review, closing.inside_tenth, closing.inside_helium) == (29, 10, 1)
    for step in figure.steps:
        assert step.expressions == 116
        assert f"{step.median_spacing:.6g}" == "0.000511364"
    assert closing.dressing.relative_gap <= report.helium_mismatch
    assert closing.second_gap > report.helium_mismatch
    exponent = report.exponent
    assert exponent.anchor.name == "e^-2+gamma^5"
    assert f"{exponent.anchor.value:.6g}" == "0.199411"
    assert f"{exponent.anchor.relative_gap:.6g}" == "0.00294716"
    assert f"{exponent.target:.6g}" == "0.2"
    assert len(exponent.steps) == 4
    first, second, third, last = exponent.steps
    assert first.dressing.name == "(e^-2+gamma^5)*(1+alpha*phi^-2)"
    assert f"{first.dressing.value:.6g}" == "0.199966"
    assert f"{first.dressing.relative_gap:.6g}" == "0.000168035"
    assert f"{first.second_gap:.6g}" == "0.000523009"
    assert (first.inside_review, first.inside_tenth, first.inside_helium) == (2, 0, 0)
    assert second.dressing.name == "((e^-2+gamma^5)*(1+alpha*phi^-2))*(1+alpha*e^-4)"
    assert f"{second.dressing.value:.6g}" == "0.199993"
    assert f"{second.dressing.relative_gap:.6g}" == "3.44016e-05"
    assert f"{second.second_gap:.6g}" == "6.72764e-05"
    assert (second.inside_review, second.inside_tenth, second.inside_helium) == (30, 3, 0)
    assert third.dressing.name == "(((e^-2+gamma^5)*(1+alpha*phi^-2))*(1+alpha*e^-4))*(1+alpha*pi^-5)"
    assert f"{third.dressing.value:.6g}" == "0.199998"
    assert f"{third.dressing.relative_gap:.6g}" == "1.05564e-05"
    assert f"{third.second_gap:.6g}" == "1.47658e-05"
    assert (third.inside_review, third.inside_tenth, third.inside_helium) == (29, 10, 0)
    assert last.dressing.name == (
        "((((e^-2+gamma^5)*(1+alpha*phi^-2))*(1+alpha*e^-4))*(1+alpha*pi^-5))*(1+alpha*pi^-6)"
    )
    assert f"{last.dressing.value:.6g}" == "0.199999"
    assert f"{last.dressing.relative_gap:.6g}" == "2.96609e-06"
    assert f"{last.second_gap:.6g}" == "7.5317e-06"
    assert (last.inside_review, last.inside_tenth, last.inside_helium) == (30, 10, 1)
    for step in exponent.steps:
        assert step.expressions == 116
        assert f"{step.median_spacing:.6g}" == "0.000511364"
    assert last.dressing.relative_gap <= report.helium_mismatch
    assert last.second_gap > report.helium_mismatch
    assert report.lead == 0.1
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_names_the_long_weight(report):
    text = render()
    assert text == render()
    assert "Equation (2) of the 1952 fit stays the witness." in text
    assert "The lead coefficient is the ten percent band." in text
    assert "The unit-interval dressing of each coefficient drops below 116 values." in text
    assert "The positive dressing from the heat-ratio score names these coefficients." in text
    assert "The filed magnitude is `0.2`." in text
    assert "Four hours over twenty hours is that same number." in text
    assert "This score does not replace equation (2)." in text
    assert "This run does not install a seed as a new decay fraction or a new plant row." in text
    assert "Every closing gap meets the helium quotient match." in text
    assert "The checked-in schedule is not rewritten." in text
    printed = set()
    for figure in (report.weight, report.exponent):
        printed.add(f"{figure.target:.6g}")
        printed.add(f"{figure.anchor.value:.6g}")
        for step in figure.steps:
            printed.add(f"{step.dressing.value:.6g}")
    closing = report.weight.steps[-1].dressing.value
    for banned in (
        report.nuclear_scalar * closing,
        report.thermodynamics_scalar * closing,
        report.gamma_fsot * closing,
        report.valve_beta * closing,
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
    stack = readme.index("## Fission stack")
    fit = readme.index("## Fission fit")
    tests = readme.index("## Tests")
    assert stack < fit < tests
    assert "py -3.14 -m solcascade.production.fission_fit" in readme
    source = (project_root() / "src" / "solcascade" / "production" / "fission_fit.py").read_text(
        encoding="utf-8"
    )
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
    assert "fission_fit" not in package
    decay = (project_root() / "src" / "solcascade" / "production" / "fission_decay.py").read_text(
        encoding="utf-8"
    )
    assert "def decay_fraction(" in decay
