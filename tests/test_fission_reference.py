"""The Pu-239 kilogram and the Surry core over the scored plant."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_reference import (
    build_reference,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_reference()


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


def test_quotients_meet_helium(report):
    helium = report.helium_mismatch
    assert 0.0 < report.kilogram.target < 1.0
    assert 0.0 < report.core_ratio.target < 1.0
    assert abs(report.delivered_ratio - report.kilogram.target * report.efficiency_step) / report.delivered_ratio < 1.0e-15
    assert abs(report.return_ratio - report.delivered_ratio) / report.delivered_ratio < 1.0e-15
    assert f"{report.efficiency_step:.6g}" == "1.21212"
    assert f"{report.delivered_ratio:.6g}" == "1.1918"
    assert f"{report.return_ratio:.6g}" == "1.1918"
    assert report.core_w == 2.587e9
    assert f"{report.heat_w:.6g}" == "3.0303e+09"
    assert f"{report.core_ratio.target:.6g}" == "0.85371"
    _assert_figure(
        report.kilogram,
        "phi^-4+G^2",
        "0.984891",
        "0.00168549",
        (
            (
                "(phi^-4+G^2)*(1-alpha*phi^-3)",
                "0.983194",
                "4.00854e-05",
                "0.00027973",
                (7, 1, 0),
            ),
            (
                "((phi^-4+G^2)*(1-alpha*phi^-3))*(1+alpha*e^-5)",
                "0.983243",
                "9.08176e-06",
                "1.62404e-05",
                (29, 10, 0),
            ),
            (
                "(((phi^-4+G^2)*(1-alpha*phi^-3))*(1+alpha*e^-5))*(1-alpha*pi^-6)",
                "0.983235",
                "1.49126e-06",
                "9.00674e-06",
                (30, 10, 1),
            ),
        ),
        helium,
    )
    assert f"{report.kilogram.target:.6g}" == "0.983234"
    assert f"{report.kilogram.steps[-1].dressing.value:.6g}" == "0.983235"
    _assert_figure(
        report.core_ratio,
        "sqrt(pi)-G",
        "0.856488",
        "0.00325433",
        (
            (
                "(sqrt(pi)-G)*(1-alpha*phi^-2)",
                "0.854101",
                "0.000457922",
                "0.000815104",
                (3, 0, 0),
            ),
            (
                "((sqrt(pi)-G)*(1-alpha*phi^-2))*(1-alpha*gamma^5)",
                "0.853702",
                "9.87252e-06",
                "5.10678e-05",
                (28, 3, 0),
            ),
            (
                "(((sqrt(pi)-G)*(1-alpha*phi^-2))*(1-alpha*gamma^5))*(1+alpha*pi^-6)",
                "0.853708",
                "2.28217e-06",
                "8.21563e-06",
                (30, 10, 1),
            ),
        ),
        helium,
    )
    assert f"{report.core_ratio.steps[-1].dressing.value:.6g}" == "0.853708"
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_names_the_reference_quotients(report):
    text = render()
    assert text == render()
    assert "Both densities are the 200 MeV ledger." in text
    assert "The unit-interval dressing of the kilogram ratio drops below 116 values." in text
    assert "The positive dressing names it." in text
    assert "(((phi^-4+G^2)*(1-alpha*phi^-3))*(1+alpha*e^-5))*(1-alpha*pi^-6)" in text
    assert "The operating ratio stays the Pu-239 kilogram over the U-235 kilogram." in text
    assert "The two rows use the same enrichment input." in text
    assert "The efficiency step stays the bounds score." in text
    assert "That delivery ratio is the kilogram ratio times the named efficiency step." in text
    assert "The fast return over the thermal return is that same product." in text
    assert "Those two comparisons keep that product." in text
    assert "The fast row stays on its own reference unit." in text
    assert "The fast row's filed barrier is the sodium system." in text
    assert "The UFSAR rates that core at `2.587e+09` W." in text
    assert "The scored operating heat is `3.0303e+09` W." in text
    assert "The scored heat is the nameplate divided by the filed efficiency." in text
    assert "(((sqrt(pi)-G)*(1-alpha*phi^-2))*(1-alpha*gamma^5))*(1+alpha*pi^-6)" in text
    assert "The operating ratio stays the Surry core rating over the scored operating heat." in text
    assert "The Surry core stays the published rating." in text
    assert "The fast power density stays an estimate." in text
    assert "The thermal power density stays the filed order." in text
    assert "The two power densities stay on their own bases." in text
    assert "This score adds no joules." in text
    assert "The catalog ceiling stays unrun." in text
    assert "The 60 GWd/t stack stays unrun." in text
    assert "The fabrication cap stays 5 percent." in text
    assert "The chain stays at k_eff of 1." in text
    assert "Every closing gap meets the helium quotient match." in text
    assert "This run does not install a seed as a new enrichment or a new plant row." in text
    assert "The checked-in schedule is not rewritten." in text
    assert "The closing value agrees with `0.983234` at six figures." not in text
    assert "The closing value agrees with `0.85371` at six figures." not in text
    assert "The fast delivery over the thermal delivery is `1.1918`." in text
    printed = {
        f"{report.efficiency_step:.6g}",
        f"{report.delivered_ratio:.6g}",
        f"{report.return_ratio:.6g}",
        f"{report.core_w:.6g}",
        f"{report.heat_w:.6g}",
        f"{report.alpha:.6g}",
        f"{report.helium_mismatch:.6g}",
        f"{report.nuclear_scalar:.6g}",
        f"{report.thermodynamics_scalar:.6g}",
        f"{report.gamma_fsot:.6g}",
        f"{report.valve_beta:.6g}",
    }
    for figure in (report.kilogram, report.core_ratio):
        printed.add(f"{figure.target:.6g}")
        printed.add(f"{figure.anchor.value:.6g}")
        printed.add(f"{figure.anchor.relative_gap:.6g}")
        for step in figure.steps:
            printed.add(f"{step.dressing.value:.6g}")
            printed.add(f"{step.dressing.relative_gap:.6g}")
            printed.add(f"{step.second_gap:.6g}")
    for banned in (
        report.nuclear_scalar * report.delivered_ratio,
        report.thermodynamics_scalar * report.delivered_ratio,
        report.gamma_fsot * report.kilogram.target,
        report.valve_beta * report.heat_w,
        report.nuclear_scalar * report.core_w,
        report.gamma_fsot * report.core_ratio.target,
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
    filed = readme.index("## Fission filed")
    reference = readme.index("## Fission reference")
    tests = readme.index("## Tests")
    assert filed < reference < tests
    assert "py -3.14 -m solcascade.production.fission_reference" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_reference.py"
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
    assert "fission_reference" not in package
    assert "import" not in package
