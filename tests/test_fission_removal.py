"""The fabrication cap and the published residual-heat duty."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_removal import (
    build_removal,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_removal()


def test_cap_and_duty_meet_helium(report):
    assert report.cap_anchor.name == "pi^-3*G^-5"
    assert f"{report.cap_anchor.value:.6g}" == "0.0500213"
    assert f"{report.cap_anchor.relative_gap:.6g}" == "0.000425086"
    added, closing = report.cap_steps
    assert added.dressing.name == "(pi^-3*G^-5)*(1-alpha*phi^-6)"
    assert f"{added.dressing.relative_gap:.6g}" == "1.82457e-05"
    assert (added.inside_review, added.inside_tenth, added.inside_helium) == (28, 3, 0)
    assert closing.dressing.name == "((pi^-3*G^-5)*(1-alpha*phi^-6))*(1-alpha*e^-6)"
    assert f"{closing.dressing.value:.6g}" == "0.05"
    assert f"{closing.dressing.relative_gap:.6g}" == "1.57058e-07"
    assert f"{closing.second_gap:.6g}" == "5.60074e-06"
    assert (closing.inside_review, closing.inside_tenth, closing.inside_helium) == (29, 10, 1)
    assert closing.second_gap > report.helium_mismatch >= closing.dressing.relative_gap
    assert report.duty_anchor.name == "pi^-2*gamma^6"
    assert f"{report.duty_anchor.value:.6g}" == "0.00374739"
    assert f"{report.duty_anchor.relative_gap:.6g}" == "0.00239397"
    first, second, duty = report.duty_steps
    assert first.dressing.name == "(pi^-2*gamma^6)*(1-alpha*gamma^2)"
    assert f"{first.dressing.relative_gap:.6g}" == "4.31661e-05"
    assert (first.inside_review, first.inside_tenth, first.inside_helium) == (4, 1, 0)
    assert second.dressing.name == "((pi^-2*gamma^6)*(1-alpha*gamma^2))*(1+alpha*e^-5)"
    assert f"{second.dressing.relative_gap:.6g}" == "6.00095e-06"
    assert (second.inside_review, second.inside_tenth, second.inside_helium) == (29, 10, 0)
    assert duty.dressing.name == (
        "(((pi^-2*gamma^6)*(1-alpha*gamma^2))*(1+alpha*e^-5))*(1-alpha*pi^-6)"
    )
    assert f"{duty.dressing.relative_gap:.6g}" == "1.58952e-06"
    assert f"{duty.second_gap:.6g}" == "1.20875e-05"
    assert (duty.inside_review, duty.inside_tenth, duty.inside_helium) == (30, 10, 1)
    assert duty.second_gap > report.helium_mismatch >= duty.dressing.relative_gap
    for step in report.cap_steps + report.duty_steps:
        assert step.expressions == 116
        assert f"{step.median_spacing:.6g}" == "0.000511364"
        assert step.median_spacing > report.helium_mismatch
    assert f"{report.duty_w:.6g}" == "9.67135e+06"
    assert f"{report.duty_fraction:.6g}" == "0.00373844"
    assert f"{report.derived_fraction:.6g}" == "0.00373843"
    assert f"{report.surry_duty_w:.6g}" == "9.67133e+06"
    assert f"{report.plant_duty_w:.6g}" == "1.13286e+07"
    assert f"{report.small_duty_w:.6g}" == "3.39858e+06"
    assert report.surry_duty_w == report.derived_fraction * report.core_w
    scale = report.plant_duty_w / report.small_duty_w
    assert abs(scale - 1000.0 / 300.0) / (1000.0 / 300.0) < 1.0e-12
    assert f"{report.service_fraction:.6g}" == "0.0110487"
    assert f"{report.cooldown_fraction:.6g}" == "0.00699764"
    assert f"{report.one_second_over_duty:.6g}" == "15.5735"
    assert f"{report.one_day_over_duty:.6g}" == "1.7697"
    assert f"{report.two_day_over_duty:.6g}" == "1.41412"
    assert f"{report.cross_hours:.6g}" == "124.798"
    assert f"{report.condenser_w:.6g}" == "2.0303e+09"
    assert report.one_second_w < report.condenser_w
    assert report.cross_hours > 48.0
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_cooling_multiples_name_the_duty_ratios(report):
    assert report.one_second_anchor.name == "gamma^-5-gamma^6"
    assert f"{report.one_second_anchor.value:.6g}" == "15.5697"
    assert f"{report.one_second_anchor.relative_gap:.6g}" == "0.00024803"
    second_steps = (
        (
            "(gamma^-5-gamma^6)*(1+alpha*pi^-3)",
            "15.5733",
            "1.27377e-05",
            "2.17974e-05",
            (29, 2, 0),
        ),
        (
            "((gamma^-5-gamma^6)*(1+alpha*pi^-3))*(1+alpha*pi^-6)",
            "15.5734",
            "5.14732e-06",
            "5.35045e-06",
            (29, 10, 0),
        ),
        (
            "(((gamma^-5-gamma^6)*(1+alpha*pi^-3))*(1+alpha*pi^-6))*(1+alpha*pi^-6)",
            "15.5736",
            "2.44306e-06",
            "1.27377e-05",
            (30, 10, 1),
        ),
    )
    assert len(report.one_second_steps) == 3
    for step, expected in zip(report.one_second_steps, second_steps, strict=True):
        name, value, gap, second, counts = expected
        assert step.dressing.name == name
        assert f"{step.dressing.value:.6g}" == value
        assert f"{step.dressing.relative_gap:.6g}" == gap
        assert f"{step.second_gap:.6g}" == second
        assert (step.inside_review, step.inside_tenth, step.inside_helium) == counts
        assert step.expressions == 116
        assert f"{step.median_spacing:.6g}" == "0.000511364"
        assert step.median_spacing > report.helium_mismatch
    closing = report.one_second_steps[-1]
    assert closing.second_gap > report.helium_mismatch >= closing.dressing.relative_gap
    assert f"{closing.dressing.value:.6g}" != f"{report.one_second_over_duty:.6g}"
    assert all(step.inside_helium == 0 for step in report.one_second_steps[:-1])
    assert report.two_day_anchor.name == "G^3+G^5"
    assert f"{report.two_day_anchor.value:.6g}" == "1.41325"
    assert f"{report.two_day_anchor.relative_gap:.6g}" == "0.000620603"
    two_day_steps = (
        (
            "(G^3+G^5)*(1+alpha*phi^-5)",
            "1.41418",
            "3.6991e-05",
            "0.000118315",
            (26, 1, 0),
        ),
        (
            "((G^3+G^5)*(1+alpha*phi^-5))*(1-alpha*e^-5)",
            "1.41411",
            "1.218e-05",
            "1.31441e-05",
            (29, 10, 0),
        ),
        (
            "(((G^3+G^5)*(1+alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6)",
            "1.41412",
            "4.58969e-06",
            "5.90809e-06",
            (30, 10, 0),
        ),
        (
            "((((G^3+G^5)*(1+alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6))*(1+alpha*pi^-6)",
            "1.41413",
            "3.00071e-06",
            "1.21801e-05",
            (30, 10, 1),
        ),
    )
    assert len(report.two_day_steps) == 4
    for step, expected in zip(report.two_day_steps, two_day_steps, strict=True):
        name, value, gap, second, counts = expected
        assert step.dressing.name == name
        assert f"{step.dressing.value:.6g}" == value
        assert f"{step.dressing.relative_gap:.6g}" == gap
        assert f"{step.second_gap:.6g}" == second
        assert (step.inside_review, step.inside_tenth, step.inside_helium) == counts
        assert step.expressions == 116
        assert f"{step.median_spacing:.6g}" == "0.000511364"
    assert report.two_day_steps[-2].dressing.relative_gap > report.helium_mismatch
    late = report.two_day_steps[-1]
    assert late.second_gap > report.helium_mismatch >= late.dressing.relative_gap
    assert f"{late.dressing.value:.6g}" != f"{report.two_day_over_duty:.6g}"
    pair = report.one_day_pair
    assert pair.anchor.name == "sqrt(pi)-e^-6"
    assert f"{pair.anchor.value:.6g}" == "1.76998"
    assert f"{pair.anchor.relative_gap:.6g}" == "0.000154343"
    assert pair.parent.dressing.name == "(sqrt(pi)-e^-6)*(1-alpha*e^-4)"
    assert f"{pair.parent.dressing.value:.6g}" == "1.76974"
    assert f"{pair.parent.dressing.relative_gap:.6g}" == "2.06672e-05"
    assert f"{pair.parent.second_gap:.6g}" == "7.94174e-05"
    assert (pair.parent.inside_review, pair.parent.inside_tenth, pair.parent.inside_helium) == (30, 3, 0)
    assert pair.parent.expressions == 116
    assert f"{pair.parent.median_spacing:.6g}" == "0.000511364"
    assert pair.parent.dressing.relative_gap > report.helium_mismatch
    assert pair.e_name == "((sqrt(pi)-e^-6)*(1-alpha*e^-4))*(1-alpha*e^-6)"
    assert f"{pair.e_value:.6g}" == "1.76971"
    assert f"{pair.e_gap:.6g}" == "2.57849e-06"
    assert pair.pi_name == "((sqrt(pi)-e^-6)*(1-alpha*e^-4))*(1-alpha*pi^-5)"
    assert f"{pair.pi_value:.6g}" == "1.7697"
    assert f"{pair.pi_gap:.6g}" == "3.17932e-06"
    assert pair.pi_value < report.one_day_over_duty < pair.e_value
    assert f"{pair.pi_value:.6g}" == f"{report.one_day_over_duty:.6g}"
    assert pair.pi_value != report.one_day_over_duty
    assert pair.e_gap <= report.helium_mismatch
    assert pair.pi_gap <= report.helium_mismatch
    assert f"{pair.separation:.6g}" == "5.75781e-06"
    assert f"{pair.piece_separation:.6g}" == "5.75769e-06"
    assert report.helium_mismatch < pair.separation < 2.0 * report.helium_mismatch
    assert abs(pair.separation - pair.piece_separation * pair.parent.dressing.value / report.one_day_over_duty) / pair.separation < 1.0e-9
    assert f"{pair.required:.6g}" == "0.00283209"
    assert f"{pair.required_e_gap:.6g}" == "0.142547"
    assert f"{pair.required_pi_gap:.6g}" == "0.133324"
    assert pair.nearest_piece == "pi^-5"
    assert f"{pair.nearest_piece_gap:.6g}" == "0.133324"
    assert pair.required_e_gap > report.helium_mismatch
    assert pair.required_pi_gap > report.helium_mismatch
    assert pair.e_child_name == "(((sqrt(pi)-e^-6)*(1-alpha*e^-4))*(1-alpha*e^-6))*(1-alpha*pi^-6)"
    assert f"{pair.e_child_gap:.6g}" == "5.01195e-06"
    assert pair.pi_child_name == "(((sqrt(pi)-e^-6)*(1-alpha*e^-4))*(1-alpha*pi^-5))*(1+alpha*pi^-6)"
    assert f"{pair.pi_child_gap:.6g}" == "4.41108e-06"
    assert pair.e_child_gap > pair.e_gap
    assert pair.pi_child_gap > pair.pi_gap
    assert pair.e_child_gap > report.helium_mismatch
    assert pair.pi_child_gap > report.helium_mismatch
    assert pair.expressions == 116
    assert f"{pair.median_spacing:.6g}" == "0.000511364"
    assert f"{pair.nested_spacing:.6g}" == "3.7316e-06"
    assert pair.nested_spacing < report.helium_mismatch
    assert closing.dressing.name not in (pair.e_name, pair.pi_name)
    assert late.dressing.name not in (pair.e_name, pair.pi_name)


def test_report_keeps_the_cap_and_the_cold_end(report):
    text = render()
    assert text == render()
    assert "The fabrication cap remains `5` percent U-235." in text
    assert "This score does not install them." in text
    assert "This score assigns that hotter rate no figure." in text
    assert "The condenser remains the early removal path." in text
    assert "A loss of condenser circulation has no second watt capacity on file." in text
    assert "This score assigns that loss no figure." in text
    assert "The third factor closes the one-second multiple." in text
    assert "(((gamma^-5-gamma^6)*(1+alpha*pi^-3))*(1+alpha*pi^-6))*(1+alpha*pi^-6)" in text
    assert "The operating multiple stays the quotient of the named one-second load and the named duty." in text
    assert "The one-day multiple sits between those two values." in text
    assert "((sqrt(pi)-e^-6)*(1-alpha*e^-4))*(1-alpha*e^-6)" in text
    assert "((sqrt(pi)-e^-6)*(1-alpha*e^-4))*(1-alpha*pi^-5)" in text
    assert "The `pi^-5` neighbor prints the same six digits as the quotient." in text
    assert "The two values stay distinct." in text
    assert "Both neighbors stay uninstalled." in text
    assert "The operating multiple stays the quotient of the named one-day load and the named duty." in text
    assert "Both routes stay for a later experiment or a later published figure sharper than their separation." in text
    assert "This score makes no prediction that chooses between them." in text
    assert "Three factors of the two-day multiple stay wider than the helium quotient match." in text
    assert "The fourth factor is one more base-piece factor of the same form." in text
    assert "((((G^3+G^5)*(1+alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6))*(1+alpha*pi^-6)" in text
    assert "The operating multiple stays the quotient of the named two-day load and the named duty." in text
    assert "Decay heat is a cooling load." in text
    assert "It is not a second electricity product." in text
    assert "The chain stays at k_eff of 1." in text
    assert "The checked-in schedule is not rewritten." in text
    assert "((pi^-3*G^-5)*(1-alpha*phi^-6))*(1-alpha*e^-6)" in text
    printed = {
        f"{report.duty_w:.6g}",
        f"{report.surry_duty_w:.6g}",
        f"{report.plant_duty_w:.6g}",
        f"{report.small_duty_w:.6g}",
        f"{report.condenser_w:.6g}",
        f"{report.one_second_w:.6g}",
        f"{report.core_w:.6g}",
        f"{report.one_second_over_duty:.6g}",
        f"{report.one_day_over_duty:.6g}",
        f"{report.two_day_over_duty:.6g}",
        f"{report.one_second_anchor.value:.6g}",
        f"{report.one_second_anchor.relative_gap:.6g}",
        f"{report.two_day_anchor.value:.6g}",
        f"{report.two_day_anchor.relative_gap:.6g}",
    }
    for step in report.one_second_steps + report.two_day_steps:
        printed.add(f"{step.dressing.value:.6g}")
        printed.add(f"{step.dressing.relative_gap:.6g}")
        printed.add(f"{step.second_gap:.6g}")
    pair = report.one_day_pair
    printed.update(
        {
            f"{pair.anchor.value:.6g}",
            f"{pair.anchor.relative_gap:.6g}",
            f"{pair.parent.dressing.value:.6g}",
            f"{pair.parent.dressing.relative_gap:.6g}",
            f"{pair.parent.second_gap:.6g}",
            f"{pair.e_value:.6g}",
            f"{pair.e_gap:.6g}",
            f"{pair.pi_value:.6g}",
            f"{pair.pi_gap:.6g}",
            f"{pair.e_child_gap:.6g}",
            f"{pair.pi_child_gap:.6g}",
            f"{pair.required:.6g}",
            f"{pair.required_e_gap:.6g}",
            f"{pair.required_pi_gap:.6g}",
            f"{pair.nearest_piece_gap:.6g}",
            f"{pair.separation:.6g}",
            f"{pair.piece_separation:.6g}",
            f"{pair.nested_spacing:.6g}",
            f"{pair.median_spacing:.6g}",
        }
    )
    for banned in (
        report.nuclear_scalar * report.plant_duty_w,
        report.thermodynamics_scalar * report.plant_duty_w,
        report.gamma_fsot * report.plant_duty_w,
        report.valve_beta * report.plant_duty_w,
        report.gamma_fsot * report.core_w,
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
    stress = readme.index("## Fission stress")
    removal = readme.index("## Fission removal")
    tests = readme.index("## Tests")
    assert stress < removal < tests
    assert "py -3.14 -m solcascade.production.fission_removal" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_removal.py"
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
    assert "fission_removal" not in package
