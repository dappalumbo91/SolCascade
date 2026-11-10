"""The condenser shares, the service multiples, and the linear scale."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_condenser import (
    build_condenser,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_condenser()


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
    _assert_figure(
        report.one_second,
        "e^-3+gamma^6",
        "0.0867723",
        "0.00142766",
        (
            (
                "(e^-3+gamma^6)*(1+alpha*gamma^3)",
                "0.0868941",
                "2.62677e-05",
                "0.000292554",
                (8, 1, 0),
            ),
            (
                "((e^-3+gamma^6)*(1+alpha*gamma^3))*(1+alpha*pi^-5)",
                "0.0868962",
                "2.42235e-06",
                "8.17989e-06",
                (29, 9, 1),
            ),
        ),
        helium,
    )
    assert f"{report.one_second.target:.6g}" == "0.0868964"
    _assert_figure(
        report.one_day,
        "pi^-4*sqrt(G)",
        "0.00982517",
        "0.00499559",
        (
            (
                "(pi^-4*sqrt(G))*(1+alpha*G^4)",
                "0.00987564",
                "0.000115427",
                "0.000314073",
                (6, 0, 0),
            ),
            (
                "((pi^-4*sqrt(G))*(1+alpha*G^4))*(1-alpha*e^-4)",
                "0.00987432",
                "1.82437e-05",
                "4.05043e-05",
                (30, 5, 0),
            ),
            (
                "(((pi^-4*sqrt(G))*(1+alpha*G^4))*(1-alpha*e^-4))*(1+alpha*e^-6)",
                "0.0098745",
                "1.55663e-07",
                "5.60193e-06",
                (29, 10, 1),
            ),
        ),
        helium,
    )
    assert f"{report.one_day.target:.6g}" == "0.0098745"
    assert f"{report.one_day.steps[-1].dressing.value:.6g}" == "0.0098745"
    _assert_figure(
        report.two_day,
        "pi^-4*G^3",
        "0.00788929",
        "0.000147628",
        (
            (
                "(pi^-4*G^3)*(1+alpha*e^-4)",
                "0.00789035",
                "1.39919e-05",
                "7.27244e-05",
                (30, 4, 0),
            ),
            (
                "((pi^-4*G^3)*(1+alpha*e^-4))*(1+alpha*e^-6)",
                "0.00789049",
                "4.09621e-06",
                "6.40155e-06",
                (29, 10, 1),
            ),
        ),
        helium,
    )
    assert f"{report.two_day.target:.6g}" == "0.00789046"
    _assert_figure(
        report.four_hour,
        "e+phi^-3",
        "2.95435",
        "0.000365678",
        (
            (
                "(e+phi^-3)*(1+alpha*e^-3)",
                "2.95542",
                "2.49663e-06",
                "4.08412e-05",
                (27, 3, 1),
            ),
        ),
        helium,
    )
    assert f"{report.four_hour.target:.6g}" == "2.95543"
    assert report.four_hour.target > 1.0
    _assert_figure(
        report.twenty_hour,
        "G+sqrt(G)",
        "1.87303",
        "0.000649592",
        (
            (
                "(G+sqrt(G))*(1-alpha*phi^-5)",
                "1.87179",
                "8.83714e-06",
                "9.02645e-05",
                (25, 2, 0),
            ),
            (
                "((G+sqrt(G))*(1-alpha*phi^-5))*(1+alpha*pi^-6)",
                "1.87181",
                "1.24679e-06",
                "9.25103e-06",
                (30, 10, 1),
            ),
        ),
        helium,
    )
    assert f"{report.twenty_hour.target:.6g}" == "1.87181"
    assert f"{report.twenty_hour.steps[-1].dressing.value:.6g}" == "1.87181"
    _assert_figure(
        report.linear_scale,
        "phi^-5+gamma",
        "0.667386",
        "0.00305832",
        (
            (
                "(phi^-5+gamma)*(1+alpha*phi^-2)",
                "0.669246",
                "0.000279506",
                "0.000634441",
                (2, 0, 0),
            ),
            (
                "((phi^-5+gamma)*(1+alpha*phi^-2))*(1+alpha*gamma^6)",
                "0.669426",
                "9.68677e-06",
                "4.42208e-05",
                (28, 3, 0),
            ),
            (
                "(((phi^-5+gamma)*(1+alpha*phi^-2))*(1+alpha*gamma^6))*(1+alpha*pi^-6)",
                "0.669432",
                "2.09642e-06",
                "8.40138e-06",
                (30, 10, 1),
            ),
        ),
        helium,
    )
    assert f"{report.linear_scale.target:.6g}" == "0.669433"
    assert report.size_ratio == 0.3
    assert report.linear_scale.target == report.size_ratio ** (1.0 / 3.0)
    assert f"{report.condenser_w:.6g}" == "2.0303e+09"
    assert f"{report.small_condenser_w:.6g}" == "6.09091e+08"
    assert f"{report.plant_duty_w:.6g}" == "1.13286e+07"
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_names_the_printed_quotients(report):
    text = render()
    assert text == render()
    assert "The small-class shares match the large-class shares." in text
    assert "The operating share stays the one-second load over the running condenser." in text
    assert "((e^-3+gamma^6)*(1+alpha*gamma^3))*(1+alpha*pi^-5)" in text
    assert "The operating share stays the one-day load over the running condenser." in text
    assert "(((pi^-4*sqrt(G))*(1+alpha*G^4))*(1-alpha*e^-4))*(1+alpha*e^-6)" in text
    assert "The closing value agrees with `0.0098745` at six figures." in text
    assert "The operating share stays the two-day load over the running condenser." in text
    assert "((pi^-4*G^3)*(1+alpha*e^-4))*(1+alpha*e^-6)" in text
    assert "The first factor closes that quotient." in text
    assert "(e+phi^-3)*(1+alpha*e^-3)" in text
    assert "The operating multiple stays the four-hour load over the scaled design-point duty." in text
    assert "The operating multiple stays the twenty-hour load over the scaled design-point duty." in text
    assert "((G+sqrt(G))*(1-alpha*phi^-5))*(1+alpha*pi^-6)" in text
    assert "The closing value agrees with `1.87181` at six figures." in text
    assert "The linear scale is the cube root of that share." in text
    assert "The unit-interval dressing of that scale drops below 116 values." in text
    assert "The operating scale stays the cube root of the small-class share." in text
    assert "(((phi^-5+gamma)*(1+alpha*phi^-2))*(1+alpha*gamma^6))*(1+alpha*pi^-6)" in text
    assert "Decay heat is a cooling load." in text
    assert "It is not a second electricity product." in text
    assert "This score assigns that hotter rate no figure." in text
    assert "This score assigns that loss no figure." in text
    assert "The 60 GWd/t stack stays unrun." in text
    assert "The fabrication cap stays 5 percent." in text
    assert "The chain stays at k_eff of 1." in text
    assert "Every closing gap meets the helium quotient match." in text
    assert "The checked-in schedule is not rewritten." in text
    printed = {
        f"{report.condenser_w:.6g}",
        f"{report.small_condenser_w:.6g}",
        f"{report.plant_duty_w:.6g}",
        f"{report.size_ratio:.6g}",
    }
    for figure in (
        report.one_second,
        report.one_day,
        report.two_day,
        report.four_hour,
        report.twenty_hour,
        report.linear_scale,
    ):
        printed.add(f"{figure.target:.6g}")
        printed.add(f"{figure.anchor.value:.6g}")
        printed.add(f"{figure.anchor.relative_gap:.6g}")
        for step in figure.steps:
            printed.add(f"{step.dressing.value:.6g}")
            printed.add(f"{step.dressing.relative_gap:.6g}")
            printed.add(f"{step.second_gap:.6g}")
    for banned in (
        report.nuclear_scalar * report.plant_duty_w,
        report.thermodynamics_scalar * report.plant_duty_w,
        report.gamma_fsot * report.condenser_w,
        report.valve_beta * report.condenser_w,
        report.gamma_fsot * report.plant_duty_w,
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
    fit = readme.index("## Fission fit")
    condenser = readme.index("## Fission condenser")
    tests = readme.index("## Tests")
    assert fit < condenser < tests
    assert "py -3.14 -m solcascade.production.fission_condenser" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_condenser.py"
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
    assert "fission_condenser" not in package
