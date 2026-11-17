"""The filed quotients the ceiling, the year, and the residual-heat table left open."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_filed import (
    build_filed,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


def _kelvin(fahrenheit: float) -> float:
    return (fahrenheit - 32.0) * 5.0 / 9.0 + 273.15


@pytest.fixture(scope="module")
def report():
    return build_filed()


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
    assert report.hot_k == 560.0
    assert report.cold_k == 310.0
    assert report.rejected.target == report.cold_k / report.hot_k
    assert report.rate_kwh == 50.0
    assert report.short_offset_s == 10.0
    assert report.long_offset_s == 2.0e7
    assert report.offset_ratio.target == report.long_offset_s / report.short_offset_s == 2.0e6
    assert report.day_share.target == 86400.0 / (365.0 * 86400.0)
    assert report.cooldown_share.target == 20.0 * 3600.0 / (365.0 * 86400.0)
    assert report.tube_ratio.target == _kelvin(124.0) / _kelvin(140.0)
    assert report.cooldown_ratio.target == _kelvin(140.0) / _kelvin(350.0)
    assert f"{report.margin_j:.6g}" == "2.93214e+11"
    _assert_figure(
        report.margin,
        "phi^-6*gamma^5",
        "0.00357079",
        "0.000177896",
        (
            (
                "(phi^-6*gamma^5)*(1+alpha*e^-4)",
                "0.00357127",
                "4.42644e-05",
                "5.74127e-05",
                (30, 3, 0),
            ),
            (
                "((phi^-6*gamma^5)*(1+alpha*e^-4))*(1+alpha*e^-5)",
                "0.00357145",
                "4.90262e-06",
                "2.04194e-05",
                (29, 10, 0),
            ),
            (
                "(((phi^-6*gamma^5)*(1+alpha*e^-4))*(1+alpha*e^-5))*(1-alpha*pi^-6)",
                "0.00357142",
                "2.68784e-06",
                "1.24931e-05",
                (30, 10, 1),
            ),
        ),
        helium,
    )
    assert f"{report.margin.target:.6g}" == "0.00357143"
    assert f"{report.margin.steps[-1].dressing.value:.6g}" == "0.00357142"
    pair = report.rejected
    assert pair.anchor.name == "gamma*sqrt(G)"
    assert f"{pair.anchor.value:.6g}" == "0.552431"
    assert f"{pair.anchor.relative_gap:.6g}" == "0.00206094"
    assert pair.anchor.relative_gap > helium
    assert pair.first.dressing.name == "(gamma*sqrt(G))*(1+alpha*phi^-3)"
    assert f"{pair.first.dressing.value:.6g}" == "0.553382"
    assert f"{pair.first.dressing.relative_gap:.6g}" == "0.000341817"
    assert f"{pair.first.second_gap:.6g}" == "0.000365368"
    assert (pair.first.inside_review, pair.first.inside_tenth, pair.first.inside_helium) == (5, 0, 0)
    assert pair.first.expressions == 116
    assert f"{pair.first.median_spacing:.6g}" == "0.000511364"
    assert pair.first.dressing.relative_gap > helium
    assert pair.parent.dressing.name == "((gamma*sqrt(G))*(1+alpha*phi^-3))*(1+alpha*e^-3)"
    assert f"{pair.parent.dressing.value:.6g}" == "0.553583"
    assert f"{pair.parent.dressing.relative_gap:.6g}" == "2.13725e-05"
    assert f"{pair.parent.second_gap:.6g}" == "6.47114e-05"
    assert (pair.parent.inside_review, pair.parent.inside_tenth, pair.parent.inside_helium) == (28, 3, 0)
    assert pair.parent.expressions == 116
    assert f"{pair.parent.median_spacing:.6g}" == "0.000511364"
    assert pair.parent.dressing.relative_gap > helium
    assert pair.e_name == "(((gamma*sqrt(G))*(1+alpha*phi^-3))*(1+alpha*e^-3))*(1-alpha*e^-6)"
    assert f"{pair.e_value:.6g}" == "0.553573"
    assert f"{pair.e_gap:.6g}" == "3.28378e-06"
    assert pair.pi_name == "(((gamma*sqrt(G))*(1+alpha*phi^-3))*(1+alpha*e^-3))*(1-alpha*pi^-5)"
    assert f"{pair.pi_value:.6g}" == "0.55357"
    assert f"{pair.pi_gap:.6g}" == "2.47404e-06"
    assert pair.pi_value < pair.target < pair.e_value
    assert f"{pair.pi_value:.6g}" != f"{pair.target:.6g}"
    assert f"{pair.e_value:.6g}" != f"{pair.target:.6g}"
    assert f"{pair.target:.6g}" == "0.553571"
    assert pair.e_gap <= helium
    assert pair.pi_gap <= helium
    assert f"{pair.separation:.6g}" == "5.75782e-06"
    assert f"{pair.piece_separation:.6g}" == "5.75769e-06"
    assert helium < pair.separation < 2.0 * helium
    assert abs(pair.separation - pair.piece_separation * pair.parent.dressing.value / pair.target) / pair.separation < 1.0e-9
    assert f"{pair.required:.6g}" == "0.00292874"
    assert f"{pair.required_e_gap:.6g}" == "0.181537"
    assert f"{pair.required_pi_gap:.6g}" == "0.103748"
    assert pair.nearest_piece == "pi^-5"
    assert f"{pair.nearest_piece_gap:.6g}" == "0.103748"
    assert pair.required_e_gap > helium
    assert pair.required_pi_gap > helium
    assert pair.e_child_name == "((((gamma*sqrt(G))*(1+alpha*phi^-3))*(1+alpha*e^-3))*(1-alpha*e^-6))*(1-alpha*pi^-6)"
    assert f"{pair.e_child_gap:.6g}" == "4.30667e-06"
    assert pair.pi_child_name == "((((gamma*sqrt(G))*(1+alpha*phi^-3))*(1+alpha*e^-3))*(1-alpha*pi^-5))*(1+alpha*pi^-6)"
    assert f"{pair.pi_child_gap:.6g}" == "5.11637e-06"
    assert pair.e_child_gap > pair.e_gap
    assert pair.pi_child_gap > pair.pi_gap
    assert pair.e_child_gap > helium
    assert pair.pi_child_gap > helium
    assert f"{pair.nested_spacing:.6g}" == "3.7316e-06"
    assert pair.nested_spacing < helium
    _assert_figure(
        report.rate,
        "e^4*G",
        "50.01",
        "0.000200539",
        (
            (
                "(e^4*G)*(1-alpha*pi^-3)",
                "49.9983",
                "3.48593e-05",
                "6.68562e-05",
                (29, 3, 0),
            ),
            (
                "((e^4*G)*(1-alpha*pi^-3))*(1+alpha*pi^-5)",
                "49.9994",
                "1.10141e-05",
                "1.43082e-05",
                (29, 10, 0),
            ),
            (
                "(((e^4*G)*(1-alpha*pi^-3))*(1+alpha*pi^-5))*(1+alpha*pi^-6)",
                "49.9998",
                "3.42376e-06",
                "7.07403e-06",
                (30, 10, 1),
            ),
        ),
        helium,
    )
    assert report.rate.target == 50.0
    assert f"{report.rate.target:.6g}" == "50"
    assert f"{report.rate.steps[-1].dressing.value:.6g}" == "49.9998"
    assert report.rate.target > 1.0
    _assert_figure(
        report.day_share,
        "pi^-5*G^2",
        "0.00274163",
        "0.000695214",
        (
            (
                "(pi^-5*G^2)*(1-alpha*phi^-5)",
                "0.00273983",
                "3.67547e-05",
                "4.46765e-05",
                (25, 2, 0),
            ),
            (
                "((pi^-5*G^2)*(1-alpha*phi^-5))*(1-alpha*e^-5)",
                "0.00273969",
                "1.24163e-05",
                "1.29078e-05",
                (29, 10, 0),
            ),
            (
                "(((pi^-5*G^2)*(1-alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6)",
                "0.00273971",
                "4.826e-06",
                "5.67177e-06",
                (30, 10, 0),
            ),
            (
                "((((pi^-5*G^2)*(1-alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6))*(1+alpha*pi^-6)",
                "0.00273973",
                "2.76439e-06",
                "1.24164e-05",
                (30, 10, 1),
            ),
        ),
        helium,
    )
    assert len(report.day_share.steps) == 4
    assert report.day_share.steps[-2].inside_helium == 0
    assert report.day_share.steps[-2].dressing.relative_gap > helium
    assert f"{report.day_share.target:.6g}" == "0.00273973"
    assert f"{report.day_share.steps[-1].dressing.value:.6g}" == "0.00273973"
    _assert_figure(
        report.cooldown_share,
        "e^-6*G",
        "0.00227045",
        "0.00554215",
        (
            (
                "(e^-6*G)*(1+alpha*sqrt(gamma))",
                "0.00228304",
                "2.87349e-05",
                "3.47021e-05",
                (5, 2, 0),
            ),
            (
                "((e^-6*G)*(1+alpha*sqrt(gamma)))*(1+alpha*pi^-5)",
                "0.00228309",
                "4.88952e-06",
                "1.06471e-05",
                (29, 9, 0),
            ),
            (
                "(((e^-6*G)*(1+alpha*sqrt(gamma)))*(1+alpha*pi^-5))*(1+alpha*pi^-6)",
                "0.00228311",
                "2.70086e-06",
                "1.24799e-05",
                (30, 10, 1),
            ),
        ),
        helium,
    )
    assert f"{report.cooldown_share.target:.6g}" == "0.00228311"
    assert f"{report.cooldown_share.steps[-1].dressing.value:.6g}" == "0.00228311"
    _assert_figure(
        report.offset_ratio,
        "pi^6*pi^6",
        "924269",
        "0.537865",
        (
            (
                "(pi^6*pi^6)*(1+alpha*e^5)",
                "1.92527e+06",
                "0.037363",
                "0.209367",
                (0, 0, 0),
            ),
            (
                "((pi^6*pi^6)*(1+alpha*e^5))*(1+alpha*gamma^-3)",
                "1.99833e+06",
                "0.000835986",
                "0.00760584",
                (1, 0, 0),
            ),
            (
                "(((pi^6*pi^6)*(1+alpha*e^5))*(1+alpha*gamma^-3))*(1+alpha*gamma^4)",
                "1.99995e+06",
                "2.66025e-05",
                "9.72281e-05",
                (24, 2, 0),
            ),
            (
                "((((pi^6*pi^6)*(1+alpha*e^5))*(1+alpha*gamma^-3))*(1+alpha*gamma^4))*(1+alpha*pi^-5)",
                "1.99999e+06",
                "2.75712e-06",
                "8.51466e-06",
                (29, 9, 1),
            ),
        ),
        helium,
    )
    assert len(report.offset_ratio.steps) == 4
    assert report.offset_ratio.steps[-2].inside_helium == 0
    assert report.offset_ratio.steps[-2].dressing.relative_gap > helium
    assert f"{report.offset_ratio.target:.6g}" == "2e+06"
    assert f"{report.offset_ratio.steps[-1].dressing.value:.6g}" == "1.99999e+06"
    assert report.offset_ratio.target > 1.0
    _assert_figure(
        report.tube_ratio,
        "phi-G^5",
        "0.973277",
        "4.24131e-05",
        (
            (
                "(phi-G^5)*(1+alpha*e^-5)",
                "0.973325",
                "6.75402e-06",
                "1.85681e-05",
                (29, 10, 0),
            ),
            (
                "((phi-G^5)*(1+alpha*e^-5))*(1-alpha*pi^-6)",
                "0.973318",
                "8.36457e-07",
                "1.13344e-05",
                (30, 10, 1),
            ),
        ),
        helium,
    )
    assert f"{report.tube_ratio.target:.6g}" == "0.973319"
    assert f"{report.tube_ratio.steps[-1].dressing.value:.6g}" == "0.973318"
    assert report.tube_ratio.target < 1.0
    _assert_figure(
        report.cooldown_ratio,
        "gamma^6+G^4",
        "0.740894",
        "0.000350223",
        (
            (
                "(gamma^6+G^4)*(1-alpha*e^-3)",
                "0.740625",
                "1.32177e-05",
                "5.65866e-05",
                (27, 3, 0),
            ),
            (
                "((gamma^6+G^4)*(1-alpha*e^-3))*(1+alpha*e^-6)",
                "0.740639",
                "4.87041e-06",
                "5.62736e-06",
                (29, 10, 0),
            ),
            (
                "(((gamma^6+G^4)*(1-alpha*e^-3))*(1+alpha*e^-6))*(1-alpha*pi^-6)",
                "0.740633",
                "2.72006e-06",
                "1.24609e-05",
                (30, 10, 1),
            ),
        ),
        helium,
    )
    assert f"{report.cooldown_ratio.target:.6g}" == "0.740635"
    assert f"{report.cooldown_ratio.steps[-1].dressing.value:.6g}" == "0.740633"
    assert report.cooldown_ratio.target < 1.0
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_names_the_filed_quotients(report):
    text = render()
    assert text == render()
    assert "The catalog ceiling stays unrun." in text
    assert "(((phi^-6*gamma^5)*(1+alpha*e^-4))*(1+alpha*e^-5))*(1-alpha*pi^-6)" in text
    assert "The operating margin stays the catalog ceiling minus the temperature Carnot." in text
    assert "The extra joules are `2.93214e+11` J." in text
    assert "The extra joules are that margin times the ranked density." in text
    assert "The unit-interval dressing of the rejected fraction drops below 116 values." in text
    assert "The positive dressing compares it." in text
    assert "The rejected fraction sits between those two values." in text
    assert "Two dressings of that parent meet the helium quotient match." in text
    assert "The separation is wider than one helium quotient match." in text
    assert "The separation is narrower than two helium quotient matches." in text
    assert "No base piece meets the required factor." in text
    assert "Another factor on either neighbor moves farther from the fraction." in text
    assert "The dressed factor stays unused." in text
    assert "Both neighbors stay uninstalled." in text
    assert "The operating fraction stays the cold temperature over the hot temperature." in text
    assert "The temperature Carnot is one minus the rejected fraction." in text
    assert "Both routes stay for a later experiment or a later published figure sharper than their separation." in text
    assert "This score makes no prediction that chooses between them." in text
    assert "The separative-work estimate stays on the enrichment basis." in text
    assert "(((e^4*G)*(1-alpha*pi^-3))*(1+alpha*pi^-5))*(1+alpha*pi^-6)" in text
    assert "The operating rate stays 50 kWh per SWU." in text
    assert "Three factors stay wider than the helium quotient match." in text
    assert "The fourth factor is one more base-piece factor of the same form." in text
    assert "((((pi^-5*G^2)*(1-alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6))*(1+alpha*pi^-6)" in text
    assert "The closing value agrees with `0.00273973` at six figures." in text
    assert "The operating share stays one day over the filed year." in text
    assert "Two days over the filed year is twice the one-day share." in text
    assert "Four hours over twenty hours stays the exponent magnitude." in text
    assert "(((e^-6*G)*(1+alpha*sqrt(gamma)))*(1+alpha*pi^-5))*(1+alpha*pi^-6)" in text
    assert "The closing value agrees with `0.00228311` at six figures." in text
    assert "The operating share stays twenty hours over the filed year." in text
    assert "Equation (2) stays the witness." in text
    assert "This score does not replace equation (2)." in text
    assert "((((pi^6*pi^6)*(1+alpha*e^5))*(1+alpha*gamma^-3))*(1+alpha*gamma^4))*(1+alpha*pi^-5)" in text
    assert "The operating ratio stays the long offset over the short offset." in text
    assert "The unit-interval dressing of that ratio drops below 116 values." in text
    assert "The positive dressing names it." in text
    assert "((phi-G^5)*(1+alpha*e^-5))*(1-alpha*pi^-6)" in text
    assert "The operating ratio stays the outlet temperature over the inlet temperature." in text
    assert "The UFSAR does not print a watt rate at the 350°F entry." in text
    assert "This score assigns that hotter rate no figure." in text
    assert "(((gamma^6+G^4)*(1-alpha*e^-3))*(1+alpha*e^-6))*(1-alpha*pi^-6)" in text
    assert "The operating ratio stays the cold-end temperature over the cooldown entry." in text
    assert "The fast ceiling source states no temperature pair." in text
    assert "The fast row's filed barrier is the sodium system." in text
    assert "The 60 GWd/t stack stays unrun." in text
    assert "The fabrication cap stays 5 percent." in text
    assert "The chain stays at k_eff of 1." in text
    assert "Every closing gap meets the helium quotient match." in text
    assert "This run does not install a seed as a new enrichment or a new plant row." in text
    assert "The checked-in schedule is not rewritten." in text
    assert "The closing value agrees with `0.00357143` at six figures." not in text
    assert "The closing value agrees with `50` at six figures." not in text
    assert "The closing value agrees with `2e+06` at six figures." not in text
    assert "The closing value agrees with `0.973319` at six figures." not in text
    assert "The closing value agrees with `0.740635` at six figures." not in text
    assert "same six digits" not in text
    assert "((pi^-4*G^3)*(1+alpha*G^2))*(1-alpha*e^-4)" not in text
    assert "(((e^-5*G^-2)*(1-alpha*gamma))*(1+alpha*e^-3))*(1-alpha*pi^-6)" not in text
    printed = {
        f"{report.hot_k:.6g}",
        f"{report.cold_k:.6g}",
        f"{report.margin_j:.6g}",
        f"{report.rate_kwh:.6g}",
        f"{report.short_offset_s:.6g}",
        f"{report.long_offset_s:.6g}",
        f"{report.alpha:.6g}",
        f"{report.helium_mismatch:.6g}",
        f"{report.nuclear_scalar:.6g}",
        f"{report.thermodynamics_scalar:.6g}",
        f"{report.gamma_fsot:.6g}",
        f"{report.valve_beta:.6g}",
    }
    for figure in (
        report.margin,
        report.rate,
        report.day_share,
        report.cooldown_share,
        report.offset_ratio,
        report.tube_ratio,
        report.cooldown_ratio,
    ):
        printed.add(f"{figure.target:.6g}")
        printed.add(f"{figure.anchor.value:.6g}")
        printed.add(f"{figure.anchor.relative_gap:.6g}")
        for step in figure.steps:
            printed.add(f"{step.dressing.value:.6g}")
            printed.add(f"{step.dressing.relative_gap:.6g}")
            printed.add(f"{step.second_gap:.6g}")
    pair = report.rejected
    printed.update(
        {
            f"{pair.target:.6g}",
            f"{pair.anchor.value:.6g}",
            f"{pair.anchor.relative_gap:.6g}",
            f"{pair.first.dressing.value:.6g}",
            f"{pair.first.dressing.relative_gap:.6g}",
            f"{pair.first.second_gap:.6g}",
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
        }
    )
    for banned in (
        report.nuclear_scalar * report.rate.target,
        report.thermodynamics_scalar * report.rate.target,
        report.gamma_fsot * report.offset_ratio.target,
        report.valve_beta * report.margin_j,
        report.gamma_fsot * report.margin.target,
        report.nuclear_scalar * report.margin_j,
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
    thermionic = readme.index("## Fission thermionic")
    filed = readme.index("## Fission filed")
    tests = readme.index("## Tests")
    assert thermionic < filed < tests
    assert "py -3.14 -m solcascade.production.fission_filed" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_filed.py"
    ).read_text(encoding="utf-8")
    assert "warp_shell_energy_j" not in source
    assert "((pi^-4*G^3)*(1+alpha*G^2))*(1-alpha*e^-4)" not in source
    assert "(((e^-5*G^-2)*(1-alpha*gamma))*(1+alpha*e^-3))*(1-alpha*pi^-6)" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9"):
        assert needle not in source
    observables = (
        project_root() / "src" / "solcascade" / "production" / "fission_observables.py"
    ).read_text(encoding="utf-8")
    assert "The channel search stays closed." in observables
    package = (project_root() / "src" / "solcascade" / "production" / "__init__.py").read_text(
        encoding="utf-8"
    )
    assert "fission_filed" not in package
    assert "import" not in package
