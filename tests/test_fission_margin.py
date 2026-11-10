"""The ledger shares and the delayed-neutron counts."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_margin import build_margin, default_output_path, main, render
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_margin()


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


def test_margin_figures_meet_helium(report):
    _check(
        report.balance,
        "e^-3*gamma^-2",
        "0.149431",
        "0.00379417",
        "0.15",
        (
            (
                "(e^-3*gamma^-2)*(1+alpha*gamma)",
                "0.15006",
                "0.000401992",
                "0.000499109",
                (3, 0, 0),
            ),
            (
                "((e^-3*gamma^-2)*(1+alpha*gamma))*(1-alpha*phi^-6)",
                "0.149999",
                "4.83941e-06",
                "3.85317e-05",
                (27, 3, 0),
            ),
            (
                "(((e^-3*gamma^-2)*(1+alpha*gamma))*(1-alpha*phi^-6))*(1+alpha*pi^-6)",
                "0.15",
                "2.75097e-06",
                "1.24298e-05",
                (30, 10, 1),
            ),
        ),
    )
    _check(
        report.fragment,
        "phi^-5+sqrt(gamma)",
        "0.849917",
        "9.75887e-05",
        "0.85",
        (
            (
                "(phi^-5+sqrt(gamma))*(1+alpha*pi^-4)",
                "0.849981",
                "2.26815e-05",
                "3.6054e-05",
                (30, 6, 0),
            ),
            (
                "((phi^-5+sqrt(gamma))*(1+alpha*pi^-4))*(1+alpha*pi^-5)",
                "0.850001",
                "1.164e-06",
                "4.59357e-06",
                (29, 10, 1),
            ),
        ),
    )
    _check(
        report.uranium,
        "pi^-4*G^5",
        "0.00661906",
        "0.00288787",
        "0.0066",
        (
            (
                "(pi^-4*G^5)*(1-alpha*phi^-2)",
                "0.00660061",
                "9.24785e-05",
                "0.000449531",
                (2, 1, 0),
            ),
            (
                "((pi^-4*G^5)*(1-alpha*phi^-2))*(1-alpha*pi^-4)",
                "0.00660012",
                "1.75571e-05",
                "4.11895e-05",
                (30, 6, 0),
            ),
            (
                "(((pi^-4*G^5)*(1-alpha*phi^-2))*(1-alpha*pi^-4))*(1-alpha*e^-6)",
                "0.0066",
                "5.31534e-07",
                "6.28933e-06",
                (29, 10, 1),
            ),
        ),
    )
    _check(
        report.plutonium,
        "pi^-5*G^4",
        "0.00230021",
        "9.08271e-05",
        "0.0023",
        (
            (
                "(pi^-5*G^4)*(1-alpha*pi^-4)",
                "0.00230004",
                "1.59058e-05",
                "4.16535e-05",
                (30, 7, 0),
            ),
            (
                "((pi^-5*G^4)*(1-alpha*pi^-4))*(1-alpha*e^-6)",
                "0.00229999",
                "2.18281e-06",
                "7.9406e-06",
                (29, 10, 1),
            ),
        ),
    )
    _check(
        report.uranium_233,
        "phi^-5*gamma^-2",
        "0.270636",
        "0.00235559",
        "0.27",
        (
            (
                "(phi^-5*gamma^-2)*(1-alpha*gamma^2)",
                "0.269978",
                "8.14533e-05",
                "0.000438316",
                (4, 1, 0),
            ),
            (
                "((phi^-5*gamma^-2)*(1-alpha*gamma^2))*(1+alpha*pi^-4)",
                "0.269998",
                "6.54487e-06",
                "3.22881e-05",
                (30, 8, 0),
            ),
            (
                "(((phi^-5*gamma^-2)*(1-alpha*gamma^2))*(1+alpha*pi^-4))*(1+alpha*pi^-6)",
                "0.27",
                "1.0455e-06",
                "1.15433e-05",
                (30, 10, 1),
            ),
        ),
    )
    named = report.uranium_233.steps[-1].dressing.value / 100.0
    assert abs(named - report.uranium_233_percent / 100.0) / (report.uranium_233_percent / 100.0) <= report.helium_mismatch
    pair = report.proportion_pair
    assert pair.target == report.uranium_233_percent / 100.0
    assert pair.scale == named
    assert pair.parent_name == "((e^-6/G)*(1-alpha*gamma^2))*(1+alpha*e^-4)"
    assert f"{pair.parent_value:.6g}" == "0.00269994"
    assert f"{pair.parent_gap:.6g}" == "2.06483e-05"
    assert pair.parent_gap > report.helium_mismatch
    assert pair.e_name == "(((e^-6/G)*(1-alpha*gamma^2))*(1+alpha*e^-4))*(1+alpha*e^-6)"
    assert f"{pair.e_value:.6g}" == "0.00269999"
    assert f"{pair.e_gap:.6g}" == "2.56039e-06"
    assert pair.pi_name == "(((e^-6/G)*(1-alpha*gamma^2))*(1+alpha*e^-4))*(1+alpha*pi^-5)"
    assert f"{pair.pi_value:.6g}" == "0.00270001"
    assert f"{pair.pi_gap:.6g}" == "3.19719e-06"
    assert pair.e_value < pair.target < pair.pi_value
    assert pair.e_gap <= report.helium_mismatch
    assert pair.pi_gap <= report.helium_mismatch
    assert f"{pair.required:.6g}" == "0.00282962"
    assert f"{pair.required_e_gap:.6g}" == "0.141552"
    assert f"{pair.required_pi_gap:.6g}" == "0.134079"
    assert pair.required_e_gap > report.helium_mismatch
    assert pair.required_pi_gap > report.helium_mismatch
    assert pair.nearest_piece == "pi^-5"
    assert f"{pair.nearest_piece_gap:.6g}" == "0.134079"
    assert pair.nearest_piece_gap > report.helium_mismatch
    assert f"{pair.separation:.6g}" == "5.75758e-06"
    assert f"{pair.piece_separation:.6g}" == "5.75769e-06"
    assert report.helium_mismatch < pair.separation < 2.0 * report.helium_mismatch
    assert abs(pair.separation - pair.piece_separation * pair.parent_value / pair.target) / pair.separation < 1.0e-9
    assert f"{pair.scale_to_low:.6g}" == "3.60589e-06"
    assert f"{pair.scale_to_high:.6g}" == "2.15168e-06"
    assert pair.scale_to_low <= report.helium_mismatch
    assert pair.scale_to_high <= report.helium_mismatch
    assert pair.scale not in (pair.e_value, pair.pi_value)
    assert f"{pair.printed_place:.6g}" == "0.01"
    assert f"{pair.half_step:.6g}" == "0.0185185"
    assert pair.e_gap < pair.half_step
    assert pair.pi_gap < pair.half_step
    assert pair.separation < pair.half_step
    assert pair.half_step > report.helium_mismatch
    assert pair.nested_spacing < report.helium_mismatch
    assert f"{pair.nested_spacing:.6g}" == "3.7316e-06"
    closing_name = report.uranium_233.steps[-1].dressing.name
    assert closing_name not in (pair.e_name, pair.pi_name)
    total = report.fragment.steps[-1].dressing.value + report.balance.steps[-1].dressing.value
    assert f"{abs(total - 1.0):.6g}" == "1.40205e-06"
    assert abs(total - 1.0) <= report.helium_mismatch
    assert report.uranium_percent == 0.66
    assert report.plutonium_percent == 0.23
    assert report.uranium_233_percent == 0.27
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_names_the_margin(report):
    text = render()
    assert text == render()
    assert "This score adds no joules to the Rankine product." in text
    assert "The delayed proportions count neutrons." in text
    assert "They are the control margin, and the control margin stays in the core." in text
    assert "The page shares add to one." in text
    assert "The two closing values add to one at the helium quotient match." in text
    assert "This proportion does not replace the thermal fuel." in text
    assert "The U-233 percent drops below 116 unit-interval dressings, so the positive dressing names it." in text
    assert "The page percent is `0.27`." in text
    assert "The proportion is that percent over 100." in text
    assert "The proportion's own dressing meets the helium quotient match more than once." in text
    assert "Another factor on that proportion moves farther from the page figure." in text
    assert "The parent of the page proportion is `((e^-6/G)*(1-alpha*gamma^2))*(1+alpha*e^-4)` at `0.00269994`." in text
    assert "The `e^-6` neighbor is `(((e^-6/G)*(1-alpha*gamma^2))*(1+alpha*e^-4))*(1+alpha*e^-6)` at `0.00269999`." in text
    assert "The `pi^-5` neighbor is `(((e^-6/G)*(1-alpha*gamma^2))*(1+alpha*e^-4))*(1+alpha*pi^-5)` at `0.00270001`." in text
    assert "The page proportion sits between those two values." in text
    assert "Their separation is `5.75758e-06`." in text
    assert "alpha times the difference of `e^-6` and `pi^-5` is `5.75769e-06`." in text
    assert "The separation equals that product times the parent over the page figure." in text
    assert "The factor the page proportion requires is `0.00282962`." in text
    assert "Its relative gap to `e^-6` is `0.141552`." in text
    assert "Its relative gap to `pi^-5` is `0.134079`." in text
    assert "Both factor gaps stay wider than the helium quotient match." in text
    assert "No base piece meets the required factor." in text
    assert "The page prints the percent to the place `0.01`." in text
    assert "Half of that place over the printed percent is `0.0185185`." in text
    assert "Both neighbors sit inside that half-step." in text
    assert "Their separation sits inside that half-step." in text
    assert "Dressing the factor by another base piece scales the spacing by alpha to `3.7316e-06`." in text
    assert "That spacing is finer than the helium quotient match." in text
    assert "The dressed factor stays unused." in text
    assert "The named percent over 100 meets both neighbors inside the helium quotient match." in text
    assert "Both neighbors stay uninstalled." in text
    assert "The operating proportion stays the printed page figure." in text
    assert "Both routes stay for a later experiment or a later published figure sharper than their separation." in text
    assert "The thermal fuel stays uranium-235." in text
    assert "The fast row's filed barrier is the sodium system." in text
    assert "The chain stays at k_eff of 1." in text
    assert "Every closing gap meets the helium quotient match." in text
    assert "The balance score's foundational seed stays on that score." in text
    assert "This run does not install a seed as a new fuel or a new plant row." in text
    assert "The checked-in schedule is not rewritten." in text
    printed = set()
    for figure in (report.balance, report.fragment, report.uranium, report.plutonium, report.uranium_233):
        printed.add(f"{figure.target:.6g}")
        printed.add(f"{figure.anchor.value:.6g}")
        for step in figure.steps:
            printed.add(f"{step.dressing.value:.6g}")
    pair = report.proportion_pair
    for value in (
        pair.parent_value,
        pair.parent_gap,
        pair.e_value,
        pair.e_gap,
        pair.pi_value,
        pair.pi_gap,
        pair.required,
        pair.required_e_gap,
        pair.required_pi_gap,
        pair.nearest_piece_gap,
        pair.separation,
        pair.piece_separation,
        pair.scale,
        pair.scale_to_low,
        pair.scale_to_high,
        pair.target,
        pair.printed_place,
        pair.half_step,
        pair.nested_spacing,
    ):
        printed.add(f"{value:.6g}")
    closing = report.balance.steps[-1].dressing.value
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
    steam = readme.index("## Fission steam")
    margin = readme.index("## Fission margin")
    tests = readme.index("## Tests")
    assert steam < margin < tests
    assert "py -3.14 -m solcascade.production.fission_margin" in readme
    source = (project_root() / "src" / "solcascade" / "production" / "fission_margin.py").read_text(
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
    assert "fission_margin" not in package
