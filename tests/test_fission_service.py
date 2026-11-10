"""The residual-heat times the removal score left on the 1952 fit."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_service import (
    build_service,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_service()


def _counts(step):
    return (step.inside_review, step.inside_tenth, step.inside_helium)


def test_service_times_meet_helium(report):
    assert report.service_anchor.name == "sqrt(e)*e^-5"
    assert f"{report.service_anchor.value:.6g}" == "0.011109"
    assert f"{report.service_anchor.relative_gap:.6g}" == "0.00545493"
    first, second, closing = report.service_steps
    assert first.dressing.name == "(sqrt(e)*e^-5)*(1-alpha*sqrt(gamma))"
    assert f"{first.dressing.relative_gap:.6g}" == "0.000119453"
    assert _counts(first) == (5, 0, 0)
    assert second.dressing.name == "((sqrt(e)*e^-5)*(1-alpha*sqrt(gamma)))*(1+alpha*e^-4)"
    assert f"{second.dressing.relative_gap:.6g}" == "1.41863e-05"
    assert _counts(second) == (30, 4, 0)
    assert closing.dressing.name == (
        "(((sqrt(e)*e^-5)*(1-alpha*sqrt(gamma)))*(1+alpha*e^-4))*(1-alpha*e^-6)"
    )
    assert f"{closing.dressing.relative_gap:.6g}" == "3.90227e-06"
    assert f"{closing.second_gap:.6g}" == "6.59578e-06"
    assert _counts(closing) == (29, 10, 1)
    assert f"{report.service_target:.6g}" == "0.0110487"
    assert f"{report.service_w:.6g}" == "3.34809e+07"
    assert f"{report.service_fit_w:.6g}" == "3.3481e+07"
    assert f"{report.service_small_w:.6g}" == "1.00443e+07"
    assert f"{report.service_fit_small_w:.6g}" == "1.00443e+07"
    assert f"{report.service_watt_gap:.6g}" == "3.82125e-06"
    assert report.service_watt_gap <= report.helium_mismatch
    added, cool = report.cooldown_steps
    assert report.cooldown_anchor.name == "pi^-4-pi^-5"
    assert f"{report.cooldown_anchor.relative_gap:.6g}" == "8.2094e-05"
    assert added.dressing.name == "(pi^-4-pi^-5)*(1-alpha*pi^-4)"
    assert f"{added.dressing.relative_gap:.6g}" == "7.17336e-06"
    assert _counts(added) == (30, 7, 0)
    assert cool.dressing.name == "((pi^-4-pi^-5)*(1-alpha*pi^-4))*(1-alpha*pi^-6)"
    assert f"{cool.dressing.relative_gap:.6g}" == "4.17122e-07"
    assert f"{cool.second_gap:.6g}" == "1.09151e-05"
    assert _counts(cool) == (30, 10, 1)
    assert f"{report.cooldown_w:.6g}" == "2.1205e+07"
    assert f"{report.cooldown_fit_w:.6g}" == "2.1205e+07"
    assert f"{report.cooldown_small_w:.6g}" == "6.36149e+06"
    assert f"{report.cooldown_watt_gap:.6g}" == "3.36101e-07"
    assert report.cooldown_watt_gap <= report.helium_mismatch
    assert report.service_steps[-1].dressing.value > report.cooldown_steps[-1].dressing.value > report.derived_duty
    first, second, cross = report.cross_steps
    assert report.cross_anchor.name == "phi^-2*gamma^6"
    assert f"{report.cross_anchor.relative_gap:.6g}" == "0.00836909"
    assert f"{report.cross_target:.6g}" == "0.0142463"
    assert first.dressing.name == "(phi^-2*gamma^6)*(1+alpha*G^-2)"
    assert f"{first.dressing.relative_gap:.6g}" == "0.000255868"
    assert _counts(first) == (2, 0, 0)
    assert second.dressing.name == "((phi^-2*gamma^6)*(1+alpha*G^-2))*(1-alpha*gamma^6)"
    assert f"{second.dressing.relative_gap:.6g}" == "1.40958e-05"
    assert _counts(second) == (29, 2, 0)
    assert cross.dressing.name == (
        "(((phi^-2*gamma^6)*(1+alpha*G^-2))*(1-alpha*gamma^6))*(1+alpha*e^-6)"
    )
    assert f"{cross.dressing.value:.6g}" == "0.0142464"
    assert f"{cross.dressing.relative_gap:.6g}" == "3.99223e-06"
    assert f"{cross.second_gap:.6g}" == "6.50552e-06"
    assert _counts(cross) == (29, 10, 1)
    assert f"{report.fit_hours:.6g}" == "124.798"
    assert f"{report.cross_hours:.6g}" == "124.798"
    assert f"{abs(report.cross_hours - report.fit_hours) / report.fit_hours:.6g}" == "3.99223e-06"
    for step in report.service_steps + report.cooldown_steps + report.cross_steps:
        assert step.expressions == 116
        assert f"{step.median_spacing:.6g}" == "0.000511364"
        assert step.median_spacing > report.helium_mismatch
        assert step.second_gap > report.helium_mismatch or step.inside_helium == 0
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_names_the_two_times(report):
    text = render()
    assert text == render()
    assert "Four hours is the service entry." in text
    assert "Twenty hours is the end of the published cooldown." in text
    assert "The 1952 fit remains the calorimetry witness." in text
    assert "The watt gap meets the helium quotient match." in text
    assert "The twenty-hour fraction stays above the design-point duty fraction." in text
    assert "The four-hour fraction stays above the design-point duty fraction." in text
    assert "This score assigns that hotter rate no figure." in text
    assert "Decay heat is a cooling load." in text
    assert "It is not a second electricity product." in text
    assert "The chain stays at k_eff of 1." in text
    assert "The checked-in schedule is not rewritten." in text
    printed = {
        f"{report.service_w:.6g}",
        f"{report.service_fit_w:.6g}",
        f"{report.service_small_w:.6g}",
        f"{report.service_fit_small_w:.6g}",
        f"{report.cooldown_w:.6g}",
        f"{report.cooldown_fit_w:.6g}",
        f"{report.cooldown_small_w:.6g}",
        f"{report.cooldown_fit_small_w:.6g}",
    }
    for banned in (
        report.nuclear_scalar * report.service_w,
        report.thermodynamics_scalar * report.service_w,
        report.gamma_fsot * report.service_w,
        report.valve_beta * report.service_w,
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
    removal = readme.index("## Fission removal")
    service = readme.index("## Fission service")
    tests = readme.index("## Tests")
    assert removal < service < tests
    assert "py -3.14 -m solcascade.production.fission_service" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_service.py"
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
    assert "fission_service" not in package
