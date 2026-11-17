"""The thermionic quotients beside the thermal plant."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_thermionic import (
    build_thermionic,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_thermionic()


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
    assert report.efficiency == 0.05
    assert report.ceiling == 0.3
    _assert_figure(
        report.delivered,
        "e*phi^-6",
        "0.151485",
        "0.000201281",
        (
            (
                "(e*phi^-6)*(1+alpha*pi^-3)",
                "0.15152",
                "3.40224e-05",
                "6.76523e-05",
                (29, 3, 0),
            ),
            (
                "((e*phi^-6)*(1+alpha*pi^-3))*(1-alpha*pi^-5)",
                "0.151517",
                "1.01755e-05",
                "1.51485e-05",
                (29, 10, 0),
            ),
            (
                "(((e*phi^-6)*(1+alpha*pi^-3))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
                "0.151516",
                "2.58504e-06",
                "7.91297e-06",
                (30, 10, 1),
            ),
        ),
        helium,
    )
    assert f"{report.delivered.target:.6g}" == "0.151515"
    assert f"{report.delivered.steps[-1].dressing.value:.6g}" == "0.151516"
    _assert_figure(
        report.enrichment_return,
        "e^5*G^3",
        "114.054",
        "0.000230898",
        (
            (
                "(e^5*G^3)*(1-alpha*pi^-3)",
                "114.027",
                "4.5073e-06",
                "3.90589e-05",
                (29, 3, 0),
            ),
            (
                "((e^5*G^3)*(1-alpha*pi^-3))*(1+alpha*pi^-6)",
                "114.028",
                "3.08309e-06",
                "1.20977e-05",
                (30, 10, 1),
            ),
        ),
        helium,
    )
    assert f"{report.enrichment_return.target:.6g}" == "114.028"
    assert f"{report.enrichment_return.steps[-1].dressing.value:.6g}" == "114.028"
    assert report.enrichment_return.target > 1.0
    _assert_figure(
        report.headroom,
        "gamma^-2+gamma^-2",
        "6.0028",
        "0.000466446",
        (
            (
                "(gamma^-2+gamma^-2)*(1-alpha*gamma^5)",
                "5.99999",
                "1.352e-06",
                "5.95888e-05",
                (28, 2, 1),
            ),
        ),
        helium,
    )
    assert f"{report.headroom.target:.6g}" == "6"
    assert abs(report.headroom.target - 6.0) / 6.0 < 1.0e-15
    assert f"{report.headroom.steps[-1].dressing.value:.6g}" == "5.99999"
    _assert_figure(
        report.ceiling_ratio,
        "gamma+gamma^2",
        "0.910394",
        "0.00143295",
        (
            (
                "(gamma+gamma^2)*(1-alpha*gamma^3)",
                "0.909116",
                "2.75425e-05",
                "0.000292192",
                (8, 1, 0),
            ),
            (
                "((gamma+gamma^2)*(1-alpha*gamma^3))*(1-alpha*pi^-5)",
                "0.909094",
                "3.69578e-06",
                "9.45363e-06",
                (29, 9, 1),
            ),
        ),
        helium,
    )
    assert f"{report.ceiling_ratio.target:.6g}" == "0.909091"
    assert report.delivered.target < report.ceiling_ratio.target < 1.0
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_names_the_printed_quotients(report):
    text = render()
    assert text == render()
    assert "The thermionic efficiency is the fabrication cap." in text
    assert "The thermionic ceiling is the small-class share." in text
    assert "One expression names each shared number." in text
    assert "The unit-interval dressing of the delivered ratio drops below 116 values." in text
    assert "The positive dressing names it." in text
    assert "(((e*phi^-6)*(1+alpha*pi^-3))*(1-alpha*pi^-5))*(1-alpha*pi^-6)" in text
    assert "The operating ratio stays the thermionic delivery over the thermal delivery." in text
    assert "The thermionic row delivers less electricity than the thermal plant." in text
    assert "((e^5*G^3)*(1-alpha*pi^-3))*(1+alpha*pi^-6)" in text
    assert "The closing value agrees with `114.028` at six figures." in text
    assert "The operating return stays the thermionic delivery over the filed enrichment input." in text
    assert "This ratio compares those two filed energies." in text
    assert "The first factor closes that quotient." in text
    assert "(gamma^-2+gamma^-2)*(1-alpha*gamma^5)" in text
    assert "The operating headroom stays the thermionic ceiling over the thermionic efficiency." in text
    assert "The fast row's filed barrier is the sodium system." in text
    assert "The fast row stays on its own reference unit." in text
    assert "((gamma+gamma^2)*(1-alpha*gamma^3))*(1-alpha*pi^-5)" in text
    assert "The operating ceiling ratio stays the thermionic file-ceiling delivery over the thermal delivery." in text
    assert "Even at its file ceiling the thermionic row delivers less than the thermal plant." in text
    assert "This score adds no joules." in text
    assert "The 60 GWd/t stack stays unrun." in text
    assert "The fabrication cap stays 5 percent." in text
    assert "The chain stays at k_eff of 1." in text
    assert "Every closing gap meets the helium quotient match." in text
    assert "The checked-in schedule is not rewritten." in text
    assert "The closing value agrees with `6` at six figures." not in text
    assert "The closing value agrees with `0.151515` at six figures." not in text
    assert "The closing value agrees with `0.909091` at six figures." not in text
    printed = {
        f"{report.efficiency:.6g}",
        f"{report.ceiling:.6g}",
    }
    for figure in (
        report.delivered,
        report.enrichment_return,
        report.headroom,
        report.ceiling_ratio,
    ):
        printed.add(f"{figure.target:.6g}")
        printed.add(f"{figure.anchor.value:.6g}")
        printed.add(f"{figure.anchor.relative_gap:.6g}")
        for step in figure.steps:
            printed.add(f"{step.dressing.value:.6g}")
            printed.add(f"{step.dressing.relative_gap:.6g}")
            printed.add(f"{step.second_gap:.6g}")
    for banned in (
        report.nuclear_scalar * report.enrichment_return.target,
        report.thermodynamics_scalar * report.enrichment_return.target,
        report.gamma_fsot * report.headroom.target,
        report.valve_beta * report.enrichment_return.target,
        report.gamma_fsot * report.delivered.target,
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
    condenser = readme.index("## Fission condenser")
    thermionic = readme.index("## Fission thermionic")
    tests = readme.index("## Tests")
    assert condenser < thermionic < tests
    assert "py -3.14 -m solcascade.production.fission_thermionic" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_thermionic.py"
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
    assert "fission_thermionic" not in package
    assert "import" not in package
