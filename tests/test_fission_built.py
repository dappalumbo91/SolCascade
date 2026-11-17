"""The plants-being-built steam steps over the demonstrated plant and up to Carnot."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_built import (
    build_built,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_built()


def _assert_figure(figure, anchor, anchor_value, anchor_gap, steps, helium):
    assert figure.anchor.name == anchor
    assert f"{figure.anchor.value:.6g}" == anchor_value
    assert f"{figure.anchor.relative_gap:.6g}" == anchor_gap
    assert figure.anchor.relative_gap > helium
    assert len(figure.steps) == len(steps)
    for step, expected in zip(figure.steps, steps, strict=True):
        name, value, gap, second, counts = expected
        assert step.dressing.name == name
        assert f"{step.dressing.value:.6g}" == value
        assert f"{step.dressing.relative_gap:.6g}" == gap
        assert f"{step.second_gap:.6g}" == second
        assert (step.inside_review, step.inside_tenth, step.inside_helium) == counts
        assert step.expressions == 116
        assert f"{step.median_spacing:.6g}" == "0.000511364"
        assert step.median_spacing > helium
    closing = figure.steps[-1]
    assert closing.second_gap > helium >= closing.dressing.relative_gap
    assert all(step.inside_helium == 0 for step in figure.steps[:-1])


def test_steps_meet_helium(report):
    helium = report.helium_mismatch
    assert report.demonstrated == 0.33
    assert report.built == 0.36
    assert report.catalog == 0.45
    assert f"{report.carnot:.6g}" == "0.446429"
    assert report.demonstrated < report.built < report.carnot < report.catalog
    assert report.built_step.target == report.built / report.demonstrated
    assert report.carnot_step.target == report.carnot / report.built
    assert abs(report.built_difference - report.demonstrated * (report.built_step.target - 1.0)) / report.built_difference < 1.0e-15
    assert abs(report.carnot_difference - report.built * (report.carnot_step.target - 1.0)) / report.carnot_difference < 1.0e-15
    assert report.fast_headroom == 1.25
    assert abs(report.catalog / report.built - report.fast_headroom) / report.fast_headroom < 1.0e-15
    _assert_figure(
        report.built_step,
        "G*G^-2",
        "1.09174",
        "0.000765392",
        (
            (
                "(G*G^-2)*(1-alpha*pi^-2)",
                "1.09094",
                "2.54494e-05",
                "4.52894e-05",
                (24, 2, 0),
            ),
            (
                "((G*G^-2)*(1-alpha*pi^-2))*(1-alpha*pi^-5)",
                "1.09091",
                "1.60279e-06",
                "7.36063e-06",
                (29, 9, 1),
            ),
        ),
        helium,
    )
    assert f"{report.built_step.target:.6g}" == "1.09091"
    assert f"{report.built_step.steps[-1].dressing.value:.6g}" == "1.09091"
    assert len(report.built_step.steps) == 2
    _assert_figure(
        report.carnot_step,
        "e^-3+G^-2",
        "1.24169",
        "0.00130057",
        (
            (
                "(e^-3+G^-2)*(1-alpha*gamma^3)",
                "1.23995",
                "0.000104654",
                "0.000234511",
                (10, 0, 0),
            ),
            (
                "((e^-3+G^-2)*(1-alpha*gamma^3))*(1+alpha*e^-4)",
                "1.24012",
                "2.89874e-05",
                "2.97476e-05",
                (30, 6, 0),
            ),
            (
                "(((e^-3+G^-2)*(1-alpha*gamma^3))*(1+alpha*e^-4))*(1-alpha*pi^-5)",
                "1.24009",
                "5.14072e-06",
                "1.08986e-05",
                (29, 9, 0),
            ),
            (
                "((((e^-3+G^-2)*(1-alpha*gamma^3))*(1+alpha*e^-4))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
                "1.24008",
                "2.44975e-06",
                "1.27312e-05",
                (30, 10, 1),
            ),
        ),
        helium,
    )
    assert len(report.carnot_step.steps) == 4
    assert report.carnot_step.steps[-2].inside_helium == 0
    assert report.carnot_step.steps[-2].dressing.relative_gap > helium
    assert f"{report.carnot_step.target:.6g}" == "1.24008"
    assert f"{report.carnot_step.steps[-1].dressing.value:.6g}" == "1.24008"
    assert f"{report.built_difference:.6g}" == "0.03"
    assert f"{report.carnot_difference:.6g}" == "0.0864286"
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_names_the_built_steps(report):
    text = render()
    assert text == render()
    assert "The demonstrated efficiency stays the heat-ratio score." in text
    assert "The plants-being-built figure stays the page efficiency." in text
    assert "((G*G^-2)*(1-alpha*pi^-2))*(1-alpha*pi^-5)" in text
    assert "The closing value agrees with `1.09091` at six figures." in text
    assert "The operating step stays the plants-being-built efficiency over the demonstrated efficiency." in text
    assert "That difference is the demonstrated efficiency times one less than this step." in text
    assert "That difference keeps the heat-ratio score and this step." in text
    assert "The temperature Carnot stays the page temperatures." in text
    assert "Three factors stay wider than the helium quotient match." in text
    assert "The fourth factor is one more base-piece factor of the same form." in text
    assert "((((e^-3+G^-2)*(1-alpha*gamma^3))*(1+alpha*e^-4))*(1-alpha*pi^-5))*(1-alpha*pi^-6)" in text
    assert "The closing value agrees with `1.24008` at six figures." in text
    assert "The operating step stays the temperature Carnot over the plants-being-built efficiency." in text
    assert "That difference is the plants-being-built efficiency times one less than this step." in text
    assert "That difference keeps the named plants-being-built efficiency and this step." in text
    assert "The catalog ceiling over the plants-being-built figure is the fast headroom." in text
    assert "The fast headroom stays the bounds score." in text
    assert "One expression names that shared number." in text
    assert "The catalog ceiling stays unrun." in text
    assert "The demonstrated plant stays at its filed efficiency." in text
    assert "This score adds no joules." in text
    assert "The fast row's filed barrier is the sodium system." in text
    assert "The 60 GWd/t stack stays unrun." in text
    assert "The fabrication cap stays 5 percent." in text
    assert "The chain stays at k_eff of 1." in text
    assert "Every closing gap meets the helium quotient match." in text
    assert "This run does not install a seed as a new efficiency or a new plant row." in text
    assert "The checked-in schedule is not rewritten." in text
    assert "phi^-5*gamma^2" not in text
    assert "phi^-6*G^-5" not in text
    printed = {
        f"{report.demonstrated:.6g}",
        f"{report.built:.6g}",
        f"{report.carnot:.6g}",
        f"{report.catalog:.6g}",
        f"{report.fast_headroom:.6g}",
        f"{report.built_difference:.6g}",
        f"{report.carnot_difference:.6g}",
        f"{report.alpha:.6g}",
        f"{report.helium_mismatch:.6g}",
        f"{report.nuclear_scalar:.6g}",
        f"{report.thermodynamics_scalar:.6g}",
        f"{report.gamma_fsot:.6g}",
        f"{report.valve_beta:.6g}",
    }
    for figure in (report.built_step, report.carnot_step):
        printed.add(f"{figure.target:.6g}")
        printed.add(f"{figure.anchor.value:.6g}")
        printed.add(f"{figure.anchor.relative_gap:.6g}")
        for step in figure.steps:
            printed.add(f"{step.dressing.value:.6g}")
            printed.add(f"{step.dressing.relative_gap:.6g}")
            printed.add(f"{step.second_gap:.6g}")
    for banned in (
        report.nuclear_scalar * report.built_step.target,
        report.thermodynamics_scalar * report.carnot_step.target,
        report.gamma_fsot * report.built,
        report.valve_beta * report.carnot,
        report.nuclear_scalar * report.catalog,
        report.gamma_fsot * report.fast_headroom,
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
    reference = readme.index("## Fission reference")
    built = readme.index("## Fission built")
    tests = readme.index("## Tests")
    assert reference < built < tests
    assert "py -3.14 -m solcascade.production.fission_built" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_built.py"
    ).read_text(encoding="utf-8")
    assert "warp_shell_energy_j" not in source
    assert "phi^-5*gamma^2" not in source
    assert "phi^-6*G^-5" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9"):
        assert needle not in source
    observables = (
        project_root() / "src" / "solcascade" / "production" / "fission_observables.py"
    ).read_text(encoding="utf-8")
    assert "The channel search stays closed." in observables
    package = (project_root() / "src" / "solcascade" / "production" / "__init__.py").read_text(
        encoding="utf-8"
    )
    assert "fission_built" not in package
    assert "import" not in package
