"""The heavy-metal discharge as a share of the ranked density."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_discharge import (
    build_discharge,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_discharge()


def _counts(step):
    return (step.inside_review, step.inside_tenth, step.inside_helium)


def test_discharge_share_meets_helium(report):
    assert report.anchor.name == "e^-3-e^-6"
    assert f"{report.anchor.value:.6g}" == "0.0473083"
    assert f"{report.anchor.relative_gap:.6g}" == "0.0010279"
    assert f"{report.target:.6g}" == "0.047357"
    assert report.discharge_gwd == 45.0
    assert report.fabrication_cap_percent == 5.0
    assert f"{report.enriched_percent:.6g}" == "4.5"
    assert f"{report.companion_ratio:.6g}" == "0.988644"
    assert f"{report.built_efficiency:.6g}" == "0.36"
    first, second, closing = report.steps
    assert first.dressing.name == "(e^-3-e^-6)*(1+alpha*phi^-4)"
    assert f"{first.dressing.value:.6g}" == "0.0473587"
    assert f"{first.dressing.relative_gap:.6g}" == "3.56701e-05"
    assert f"{first.second_gap:.6g}" == "4.13308e-05"
    assert _counts(first) == (15, 2, 0)
    assert second.dressing.name == "((e^-3-e^-6)*(1+alpha*phi^-4))*(1-alpha*pi^-5)"
    assert f"{second.dressing.value:.6g}" == "0.0473576"
    assert f"{second.dressing.relative_gap:.6g}" == "1.18232e-05"
    assert f"{second.second_gap:.6g}" == "1.35008e-05"
    assert _counts(second) == (29, 10, 0)
    assert closing.dressing.name == (
        "(((e^-3-e^-6)*(1+alpha*phi^-4))*(1-alpha*pi^-5))*(1-alpha*pi^-6)"
    )
    assert f"{closing.dressing.value:.6g}" == "0.0473572"
    assert f"{closing.dressing.relative_gap:.6g}" == "4.2327e-06"
    assert f"{closing.second_gap:.6g}" == "6.26533e-06"
    assert _counts(closing) == (30, 10, 1)
    assert closing.dressing.relative_gap <= report.helium_mismatch
    assert closing.second_gap > report.helium_mismatch
    for step in report.steps:
        assert step.expressions == 116
        assert f"{step.median_spacing:.6g}" == "0.000511364"
        assert step.median_spacing > report.helium_mismatch
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_names_the_discharge_share(report):
    text = render()
    assert text == render()
    assert "The heavy-metal discharge is the notes figure of 45 GWd per tonne." in text
    assert "The ranked density stays the filed ledger." in text
    assert "The chain stays at k_eff of 1." in text
    assert "The relative gap meets the helium quotient match." in text
    assert "One dressing of that anchor meets the helium quotient match." in text
    assert "The 45 GWd/t discharge stays the material balance." in text
    assert "The fabrication cap stays 5 percent." in text
    assert "The 60 GWd/t stack stays unrun." in text
    assert "The page enrichment stays `4.5` percent." in text
    assert "The fuel-year companion stays the page pair." in text
    assert "The plants-being-built efficiency stays the page figure." in text
    assert "The temperature Carnot stays the page temperatures." in text
    assert "This run does not install a seed as a new burn-up or a new plant row." in text
    assert "The process score's foundational seed stays on that score." in text
    assert "Every closing gap meets the helium quotient match." in text
    assert "The checked-in schedule is not rewritten." in text
    printed = {
        f"{report.target:.6g}",
        f"{report.anchor.value:.6g}",
        f"{report.steps[0].dressing.value:.6g}",
        f"{report.steps[1].dressing.value:.6g}",
        f"{report.steps[2].dressing.value:.6g}",
    }
    for banned in (
        report.nuclear_scalar * report.steps[-1].dressing.value,
        report.thermodynamics_scalar * report.steps[-1].dressing.value,
        report.gamma_fsot * report.steps[-1].dressing.value,
        report.valve_beta * report.steps[-1].dressing.value,
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
    service = readme.index("## Fission service")
    discharge = readme.index("## Fission discharge")
    tests = readme.index("## Tests")
    assert service < discharge < tests
    assert "py -3.14 -m solcascade.production.fission_discharge" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_discharge.py"
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
    assert "fission_discharge" not in package
