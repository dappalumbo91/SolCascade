"""The filed plant bounds and the 1952 stated bands."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_bounds import build_bounds, default_output_path, main, render
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_bounds()


def _counts(step):
    return (step.inside_review, step.inside_tenth, step.inside_helium)


def _check(figure, anchor, anchor_value, anchor_gap, target, steps):
    assert figure.anchor.name == anchor
    assert f"{figure.anchor.value:.6g}" == anchor_value
    assert f"{figure.anchor.relative_gap:.6g}" == anchor_gap
    assert f"{figure.target:.6g}" == target
    assert len(figure.steps) == len(steps)
    for step, expected in zip(figure.steps, steps, strict=True):
        name, value, gap, second, counts = expected
        assert step.dressing.name == name
        assert f"{step.dressing.value:.6g}" == value
        assert f"{step.dressing.relative_gap:.6g}" == gap
        assert f"{step.second_gap:.6g}" == second
        assert _counts(step) == counts
        assert step.expressions == 116
        assert f"{step.median_spacing:.6g}" == "0.000511364"
    assert figure.steps[-1].inside_helium == 1
    assert figure.steps[-1].dressing.relative_gap <= 4.23431e-06
    assert figure.steps[-1].second_gap > 4.23431e-06


def test_bounds_figures_meet_helium(report):
    _check(
        report.ceiling,
        "sqrt(gamma)*G^6",
        "0.448688",
        "0.00291655",
        "0.45",
        (
            (
                "(sqrt(gamma)*G^6)*(1+alpha*phi^-2)",
                "0.449938",
                "0.000137343",
                "0.000492329",
                (2, 0, 0),
            ),
            (
                "((sqrt(gamma)*G^6)*(1+alpha*phi^-2))*(1+alpha*e^-4)",
                "0.449998",
                "3.70603e-06",
                "6.24391e-05",
                (30, 4, 1),
            ),
        ),
    )
    _check(
        report.headroom,
        "sqrt(phi)+phi^-5",
        "1.36219",
        "0.00106096",
        "1.36364",
        (
            (
                "(sqrt(phi)+phi^-5)*(1+alpha*phi^-4)",
                "1.36364",
                "2.57487e-06",
                "7.44235e-05",
                (14, 2, 1),
            ),
        ),
    )
    _check(
        report.fast_efficiency,
        "phi^-3*G^-6",
        "0.399726",
        "0.000685652",
        "0.4",
        (
            (
                "(phi^-3*G^-6)*(1+alpha*phi^-5)",
                "0.399989",
                "2.81017e-05",
                "5.3217e-05",
                (25, 2, 0),
            ),
            (
                "((phi^-3*G^-6)*(1+alpha*phi^-5))*(1+alpha*pi^-5)",
                "0.399998",
                "4.25634e-06",
                "1.00139e-05",
                (29, 9, 0),
            ),
            (
                "(((phi^-3*G^-6)*(1+alpha*phi^-5))*(1+alpha*pi^-5))*(1+alpha*pi^-6)",
                "0.400001",
                "3.33405e-06",
                "1.18467e-05",
                (30, 10, 1),
            ),
        ),
    )
    _check(
        report.fast_ceiling,
        "sqrt(pi)-sqrt(phi)",
        "0.500434",
        "0.000868403",
        "0.5",
        (
            (
                "(sqrt(pi)-sqrt(phi))*(1-alpha*gamma^4)",
                "0.500029",
                "5.76382e-05",
                "0.000120044",
                (23, 1, 0),
            ),
            (
                "((sqrt(pi)-sqrt(phi))*(1-alpha*gamma^4))*(1-alpha*e^-5)",
                "0.500004",
                "8.46624e-06",
                "1.72806e-05",
                (29, 9, 0),
            ),
            (
                "(((sqrt(pi)-sqrt(phi))*(1-alpha*gamma^4))*(1-alpha*e^-5))*(1-alpha*pi^-6)",
                "0.5",
                "8.7575e-07",
                "9.62224e-06",
                (30, 10, 1),
            ),
        ),
    )
    _check(
        report.fast_headroom,
        "gamma^2+G",
        "1.24914",
        "0.000685186",
        "1.25",
        (
            (
                "(gamma^2+G)*(1+alpha*phi^-5)",
                "1.24997",
                "2.76346e-05",
                "5.36842e-05",
                (25, 2, 0),
            ),
            (
                "((gamma^2+G)*(1+alpha*phi^-5))*(1+alpha*pi^-5)",
                "1.25",
                "3.78923e-06",
                "9.54677e-06",
                (29, 9, 1),
            ),
        ),
    )
    _check(
        report.small_class,
        "phi^-3+gamma^5",
        "0.300143",
        "0.00047754",
        "0.3",
        (
            (
                "(phi^-3+gamma^5)*(1-alpha*gamma^5)",
                "0.300003",
                "9.73714e-06",
                "7.06786e-05",
                (28, 2, 0),
            ),
            (
                "((phi^-3+gamma^5)*(1-alpha*gamma^5))*(1-alpha*pi^-6)",
                "0.300001",
                "2.14664e-06",
                "8.35137e-06",
                (30, 10, 1),
            ),
        ),
    )
    _check(
        report.late_band,
        "e-phi^2",
        "0.100248",
        "0.0024784",
        "0.1",
        (
            (
                "(e-phi^2)*(1-alpha*gamma^2)",
                "0.100004",
                "4.10545e-05",
                "0.000315852",
                (3, 1, 0),
            ),
            (
                "((e-phi^2)*(1-alpha*gamma^2))*(1-alpha*e^-5)",
                "0.0999992",
                "8.11665e-06",
                "1.72075e-05",
                (29, 10, 0),
            ),
            (
                "(((e-phi^2)*(1-alpha*gamma^2))*(1-alpha*e^-5))*(1+alpha*pi^-6)",
                "0.0999999",
                "5.26284e-07",
                "9.97153e-06",
                (30, 10, 1),
            ),
        ),
    )
    _check(
        report.thermal_return,
        "pi^6/sqrt(phi)",
        "755.797",
        "0.00427316",
        "752.582",
        (
            (
                "(pi^6/sqrt(phi))*(1-alpha*gamma)",
                "752.614",
                "4.30112e-05",
                "5.48922e-05",
                (4, 2, 0),
            ),
            (
                "((pi^6/sqrt(phi))*(1-alpha*gamma))*(1-alpha*e^-5)",
                "752.577",
                "6.16005e-06",
                "1.91642e-05",
                (29, 10, 0),
            ),
            (
                "(((pi^6/sqrt(phi))*(1-alpha*gamma))*(1-alpha*e^-5))*(1+alpha*pi^-6)",
                "752.583",
                "1.43033e-06",
                "1.19282e-05",
                (30, 10, 1),
            ),
        ),
    )
    _check(
        report.fast_return,
        "pi^6-e^4",
        "906.791",
        "0.0109992",
        "896.926",
        (
            (
                "(pi^6-e^4)*(1-alpha*G^-5)",
                "896.528",
                "0.000443264",
                "0.000518299",
                (3, 0, 0),
            ),
            (
                "((pi^6-e^4)*(1-alpha*G^-5))*(1+alpha*gamma^5)",
                "896.947",
                "2.41091e-05",
                "3.67763e-05",
                (28, 3, 0),
            ),
            (
                "(((pi^6-e^4)*(1-alpha*G^-5))*(1+alpha*gamma^5))*(1-alpha*pi^-5)",
                "896.926",
                "2.6253e-07",
                "6.02036e-06",
                (29, 10, 1),
            ),
        ),
    )
    _check(
        report.efficiency_step,
        "e^-4+G^-2",
        "1.21022",
        "0.00156789",
        "1.21212",
        (
            (
                "(e^-4+G^-2)*(1+alpha*phi^-3)",
                "1.21231",
                "0.00015208",
                "0.000166696",
                (8, 0, 0),
            ),
            (
                "((e^-4+G^-2)*(1+alpha*phi^-3))*(1-alpha*e^-4)",
                "1.21214",
                "1.84044e-05",
                "7.71545e-05",
                (30, 3, 0),
            ),
            (
                "(((e^-4+G^-2)*(1+alpha*phi^-3))*(1-alpha*e^-4))*(1-alpha*e^-6)",
                "1.21212",
                "3.15745e-07",
                "5.44206e-06",
                (29, 10, 1),
            ),
        ),
    )
    assert report.early_band == report.fast_ceiling.target
    assert report.middle_band == report.small_class.target
    assert f"{report.demonstrated_efficiency:.6g}" == "0.33"
    assert f"{report.carnot:.6g}" == "0.446429"
    assert report.fabrication_cap_percent == 5.0
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_names_the_bounds(report):
    text = render()
    assert text == render()
    assert "The catalog ceiling stays above the temperature Carnot and is not run." in text
    assert "The fast row keeps its sodium barrier." in text
    assert "The demonstrated thermal efficiency stays the heat-ratio score." in text
    assert "Ratios that leave the unit interval use the positive dressing from the heat-ratio score." in text
    assert "The chain stays at k_eff of 1." in text
    assert "The 1952 band under 100 seconds is that same number." in text
    assert "The fast ceiling source states no temperature pair, so this score does not tighten that row." in text
    assert "The 1952 band from 100 seconds to 10,000 seconds is that same number." in text
    assert "The band stays an uncertainty on the fit." in text
    assert "The small class stays the published electric boundary." in text
    assert "Four hours and twenty hours sit in this band." in text
    assert "The equation (2) lead coefficient is that same number." in text
    assert "The 1952 fit remains the witness." in text
    assert "This ratio compares those two filed energies." in text
    assert "The phase-two burnup readings stay at one." in text
    assert "The fabrication cap stays 5 percent." in text
    assert "This run does not install a seed as a new efficiency or a new plant row." in text
    assert "Every closing gap meets the helium quotient match." in text
    assert "The checked-in schedule is not rewritten." in text
    printed = set()
    figures = (
        report.ceiling,
        report.headroom,
        report.fast_efficiency,
        report.fast_ceiling,
        report.fast_headroom,
        report.small_class,
        report.late_band,
        report.thermal_return,
        report.fast_return,
        report.efficiency_step,
    )
    for figure in figures:
        printed.add(f"{figure.target:.6g}")
        printed.add(f"{figure.anchor.value:.6g}")
        for step in figure.steps:
            printed.add(f"{step.dressing.value:.6g}")
    closing = report.ceiling.steps[-1].dressing.value
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
    margin = readme.index("## Fission margin")
    bounds = readme.index("## Fission bounds")
    tests = readme.index("## Tests")
    assert margin < bounds < tests
    assert "py -3.14 -m solcascade.production.fission_bounds" in readme
    source = (project_root() / "src" / "solcascade" / "production" / "fission_bounds.py").read_text(
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
    assert "fission_bounds" not in package
