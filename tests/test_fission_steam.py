"""The steam figures above the demonstrated efficiency."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_steam import build_steam, default_output_path, main, render
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_steam()


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
    assert figure.steps[-1].dressing.relative_gap <= report_helium(figure)
    assert figure.steps[-1].second_gap > figure.steps[-1].dressing.relative_gap


def report_helium(figure):
    del figure
    return 4.23431e-06


def test_steam_figures_meet_helium(report):
    _check(
        report.built,
        "pi^-3*phi^5",
        "0.357675",
        "0.00645834",
        "0.36",
        (
            (
                "(pi^-3*phi^5)*(1+alpha*G)",
                "0.360066",
                "0.000182616",
                "0.000375453",
                (5, 0, 0),
            ),
            (
                "((pi^-3*phi^5)*(1+alpha*G))*(1-alpha*e^-4)",
                "0.360018",
                "4.8936e-05",
                "5.27777e-05",
                (30, 3, 0),
            ),
            (
                "(((pi^-3*phi^5)*(1+alpha*G))*(1-alpha*e^-4))*(1-alpha*e^-5)",
                "0.36",
                "2.3554e-07",
                "2.50889e-05",
                (29, 10, 1),
            ),
        ),
    )
    _check(
        report.carnot,
        "phi^-2+gamma^5",
        "0.446041",
        "0.000867497",
        "0.446429",
        (
            (
                "(phi^-2+gamma^5)*(1+alpha*gamma^4)",
                "0.446403",
                "5.81389e-05",
                "0.000119235",
                (23, 1, 0),
            ),
            (
                "((phi^-2+gamma^5)*(1+alpha*gamma^4))*(1+alpha*e^-5)",
                "0.446425",
                "8.97261e-06",
                "1.67712e-05",
                (29, 9, 0),
            ),
            (
                "(((phi^-2+gamma^5)*(1+alpha*gamma^4))*(1+alpha*e^-5))*(1+alpha*pi^-6)",
                "0.446428",
                "1.38225e-06",
                "9.11556e-06",
                (30, 10, 1),
            ),
        ),
    )
    _check(
        report.carnot_step,
        "e^-3+G^-3",
        "1.35104",
        "0.00130947",
        "1.35281",
        (
            (
                "(e^-3+G^-3)*(1+alpha*gamma^3)",
                "1.35294",
                "9.20885e-05",
                "0.000246193",
                (10, 1, 0),
            ),
            (
                "((e^-3+G^-3)*(1+alpha*gamma^3))*(1-alpha*pi^-4)",
                "1.35284",
                "1.71671e-05",
                "4.15795e-05",
                (30, 7, 0),
            ),
            (
                "(((e^-3+G^-3)*(1+alpha*gamma^3))*(1-alpha*pi^-4))*(1-alpha*e^-6)",
                "1.35281",
                "9.21512e-07",
                "6.67931e-06",
                (29, 10, 1),
            ),
        ),
    )
    built = report.built.steps[-1].dressing.value
    carnot = report.carnot.steps[-1].dressing.value
    assert report.demonstrated_efficiency < built < report.carnot.target
    assert report.demonstrated_efficiency < carnot < report.catalog_ceiling
    assert f"{report.demonstrated_efficiency:.6g}" == "0.33"
    assert f"{report.catalog_ceiling:.6g}" == "0.45"
    assert report.fabrication_cap_percent == 5.0
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix
    for figure in (report.built, report.carnot, report.carnot_step):
        assert figure.steps[-1].dressing.relative_gap <= report.helium_mismatch
        assert figure.steps[-1].second_gap > report.helium_mismatch


def test_report_names_the_steam_figures(report):
    text = render()
    assert text == render()
    assert "The temperature Carnot stays the page temperatures." in text
    assert "The catalog ceiling stays above that Carnot and is not run." in text
    assert "The demonstrated efficiency stays the heat-ratio score." in text
    assert "The 60 GWd/t stack stays unrun." in text
    assert "The fabrication cap stays 5 percent." in text
    assert "The catalog ceiling stays `0.45`." in text
    assert "The chain stays at k_eff of 1." in text
    assert "Every closing gap meets the helium quotient match." in text
    assert "This run does not install a seed as a new efficiency or a new plant row." in text
    assert "The checked-in schedule is not rewritten." in text
    printed = {f"{report.catalog_ceiling:.6g}", f"{report.demonstrated_efficiency:.6g}"}
    for figure in (report.built, report.carnot, report.carnot_step):
        printed.add(f"{figure.target:.6g}")
        printed.add(f"{figure.anchor.value:.6g}")
        for step in figure.steps:
            printed.add(f"{step.dressing.value:.6g}")
    closing = report.carnot.steps[-1].dressing.value
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
    reload = readme.index("## Fission reload")
    steam = readme.index("## Fission steam")
    tests = readme.index("## Tests")
    assert reload < steam < tests
    assert "py -3.14 -m solcascade.production.fission_steam" in readme
    source = (project_root() / "src" / "solcascade" / "production" / "fission_steam.py").read_text(
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
    assert "fission_steam" not in package
