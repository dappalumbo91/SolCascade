"""The material-balance year ratios."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_reload import build_reload, default_output_path, main, render
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_reload()


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
    assert figure.steps[-1].dressing.relative_gap <= 4.23431e-06
    assert figure.steps[-1].second_gap > 4.23431e-06


def test_reload_ratios_meet_helium(report):
    _check(
        report.companion,
        "sqrt(gamma)*G^-3",
        "0.988625",
        "1.90905e-05",
        "0.988644",
        (
            (
                "(sqrt(gamma)*G^-3)*(1+alpha*e^-6)",
                "0.988643",
                "1.00247e-06",
                "4.75511e-06",
                (29, 10, 1),
            ),
        ),
    )
    _check(
        report.separative,
        "pi^-2+e^2",
        "7.49038",
        "8.8835e-05",
        "7.48971",
        (
            (
                "(pi^-2+e^2)*(1-alpha*pi^-4)",
                "7.48982",
                "1.39138e-05",
                "3.96614e-05",
                (30, 7, 0),
            ),
            (
                "((pi^-2+e^2)*(1-alpha*pi^-4))*(1-alpha*e^-6)",
                "7.48968",
                "4.17475e-06",
                "6.3233e-06",
                (29, 10, 1),
            ),
        ),
    )
    _check(
        report.highly_enriched,
        "phi^5+gamma^-5",
        "26.6968",
        "0.000242861",
        "26.7033",
        (
            (
                "(phi^5+gamma^-5)*(1+alpha*pi^-3)",
                "26.7031",
                "7.56759e-06",
                "2.69677e-05",
                (29, 2, 0),
            ),
            (
                "((phi^5+gamma^-5)*(1+alpha*pi^-3))*(1+alpha*pi^-6)",
                "26.7033",
                "2.27791e-08",
                "1.05206e-05",
                (30, 10, 1),
            ),
        ),
    )
    _check(
        report.energy_return,
        "pi^6-gamma^-4",
        "952.381",
        "0.000709569",
        "951.705",
        (
            (
                "(pi^6-gamma^-4)*(1-alpha*pi^-2)",
                "951.677",
                "3.03317e-05",
                "5.11006e-05",
                (25, 2, 0),
            ),
            (
                "((pi^6-gamma^-4)*(1-alpha*pi^-2))*(1+alpha*pi^-5)",
                "951.699",
                "6.48635e-06",
                "1.22439e-05",
                (29, 9, 0),
            ),
            (
                "(((pi^6-gamma^-4)*(1-alpha*pi^-2))*(1+alpha*pi^-5))*(1+alpha*pi^-6)",
                "951.707",
                "1.10402e-06",
                "1.16019e-05",
                (30, 10, 1),
            ),
        ),
    )
    assert report.annual_swu == 182_000.0
    assert report.reload_tonnes == 24.3
    assert report.fabrication_cap_percent == 5.0
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_names_the_reload_ratios(report):
    text = render()
    assert text == render()
    assert "The fuel-side year and the nameplate year stay that page pair." in text
    assert "The highly enriched estimate stays the catalog input." in text
    assert "It is not a second plant output." in text
    assert "The page enrichment stays 4.5 percent." in text
    assert "The fabrication cap stays 5 percent." in text
    assert "The 60 GWd/t stack stays unrun." in text
    assert "The chain stays at k_eff of 1." in text
    assert "Every closing gap meets the helium quotient match." in text
    assert "This run does not install a seed as a new enrichment or a new plant row." in text
    assert "The checked-in schedule is not rewritten." in text
    printed = set()
    for figure in (report.companion, report.separative, report.highly_enriched, report.energy_return):
        printed.add(f"{figure.target:.6g}")
        printed.add(f"{figure.anchor.value:.6g}")
        for step in figure.steps:
            printed.add(f"{step.dressing.value:.6g}")
    closing = report.energy_return.steps[-1].dressing.value
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
    discharge = readme.index("## Fission discharge")
    reload = readme.index("## Fission reload")
    tests = readme.index("## Tests")
    assert discharge < reload < tests
    assert "py -3.14 -m solcascade.production.fission_reload" in readme
    source = (project_root() / "src" / "solcascade" / "production" / "fission_reload.py").read_text(
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
    assert "fission_reload" not in package
