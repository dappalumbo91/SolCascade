"""Name the overview figures the operating column left aside."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_aside import (
    build_aside,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root

_CLOSINGS = {
    "comparison concentrate": "((pi^-4+G^2)*(1-alpha*e^-3))*(1+alpha*pi^-6)",
    "comparison hex": "(((pi^-3+G^5)*(1-alpha*e^-2))*(1+alpha*e^-6))*(1-alpha*pi^-6)",
    "comparison enriched hex": "(sqrt(G)*G^4)*(1+alpha*gamma^5)",
    "comparison oxide": "(((phi^-3+G^5)*(1+alpha*gamma^3))*(1+alpha*e^-4))*(1+alpha*pi^-5)",
    "comparison per terawatt-hour": "((gamma^-2+gamma^-5)*(1-alpha*e^-5))*(1+alpha*pi^-6)",
    "comparison transuranic": "((e^-5*G^-5)*(1+alpha*phi^-3))*(1+alpha*gamma^6)",
    "comparison used oxide": "(pi^-6+sqrt(G))*(1+alpha*e^-6)",
    "red book tails percent": "((gamma^3*G^-3)*(1-alpha*e^-2))*(1-alpha*e^-6)",
    "red book demand": "(((e^5+gamma^-6)*(1-alpha*gamma^2))*(1-alpha*e^-4))*(1-alpha*pi^-6)",
    "red book lower demand": (
        "((((e^5+phi^5)*(1+alpha*phi^-2))*(1+alpha*e^-3))*(1-alpha*e^-5))*(1+alpha*pi^-6)"
    ),
    "red book oxide": (
        "((((e^4+e^5)*(1+alpha*sqrt(pi)))*(1+alpha*phi^-3))*(1+alpha*e^-5))*(1-alpha*pi^-6)"
    ),
    "red book lower oxide": "(((pi^5/phi)*(1-alpha*phi^-5))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
    "candu uranium 238": "(((pi^4+G^-4)*(1-alpha*gamma^6))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
    "candu plutonium": (
        "((((G^5*G^6)*(1-alpha*phi^-3))*(1-alpha*e^-3))*(1+alpha*e^-5))*(1-alpha*pi^-6)"
    ),
    "candu efficiency": "(((phi^-3/sqrt(gamma))*(1-alpha*gamma^2))*(1+alpha*e^-4))*(1-alpha*e^-6)",
    "candu demand": "(((e^5+gamma^-4)*(1-alpha*phi^-2))*(1+alpha*e^-4))*(1-alpha*pi^-5)",
    "candu capacity demand": "((e^5*sqrt(G))*(1-alpha*gamma^6))*(1-alpha*e^-6)",
}

_FOUR_FACTORS = {
    "red book lower demand",
    "red book oxide",
    "candu plutonium",
}

_KEPT = (0.30, 0.0023, 0.23, 0.0038, 0.988, 7.5, 18.6, 17.9, 0.80, 0.90)


@pytest.fixture(scope="module")
def report():
    return build_aside()


def test_aside_figures_meet_helium_without_a_fifth_factor(report):
    assert len(report.figures) == len(_CLOSINGS) == 17
    assert [figure.label for figure in report.figures] == list(_CLOSINGS)
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert f"{report.alpha:.6g}" == "0.00729735"
    unit_named = [figure.label for figure in report.figures if not figure.positive]
    assert unit_named == ["comparison transuranic"]
    for figure in report.figures:
        assert figure.target not in _KEPT
        assert 1 <= len(figure.steps) <= 4
        if figure.label in _FOUR_FACTORS:
            assert len(figure.steps) == 4
        previous = figure.anchor.relative_gap
        for step in figure.steps[:-1]:
            assert step.inside_helium == 0
            assert step.dressing.relative_gap < previous
            assert step.expressions == 116
            assert f"{step.median_spacing:.6g}" == "0.000511364"
            assert step.median_spacing > report.helium_mismatch
            previous = step.dressing.relative_gap
        closing = figure.steps[-1]
        assert closing.dressing.name == _CLOSINGS[figure.label]
        assert closing.inside_helium == 1
        assert closing.dressing.relative_gap <= report.helium_mismatch < closing.second_gap
        assert closing.dressing.relative_gap < previous
        assert closing.expressions == 116
        assert f"{closing.median_spacing:.6g}" == "0.000511364"
    assert report.compare_burnup_step == 65.0 / 45.0
    assert report.candu_burnup_step == 7.5 / 45.0


def test_report_keeps_the_unrun_rows_and_the_shared_names(report):
    text = render()
    assert text == render()
    assert "# Fission aside" in text
    assert "This comparison stays unrun." in text
    assert "That enrichment is the fabrication cap." in text
    assert "The page's 18.6 is not a second target." in text
    assert "The comparison burn-up over the operating discharge is `1.44444`." in text
    assert "The same quotient is `13/9`." in text
    assert "That burn-up stays unrun." in text
    assert "The comparison tails assay stays the named 0.22 percent." in text
    assert "The enrichment note's about 140000 SWU stays approximate." in text
    assert "The mining line stays a range." in text
    assert "The feed multiple stays the named rational `3260/337`." in text
    assert "Those assays stay off this column." in text
    assert "The 0.30 percent assay keeps the named small-class share." in text
    assert "The 0.30 percent fraction is `3/1000`." in text
    assert "The 0.25 percent fraction is `1/400`." in text
    assert "The same quotient is `175/206`." in text
    assert "The same quotient is `160/189`." in text
    assert "The CANDU row stays a different plant." in text
    assert "The same quotient is `494/5`." in text
    assert "The 0.23 percent figure is the delayed plutonium count." in text
    assert "Its fraction keeps the named plutonium margin." in text
    assert "This score does not add a second expression for that count." in text
    assert "The 0.38 percent fraction is `19/5000`." in text
    assert "The 0.6 percent fission products stay one digit." in text
    assert "The same quotient is `1/6`." in text
    assert "The almost 50 GWd/t equivalent stays prose." in text
    assert "The same quotient is `3925/219`." in text
    assert "The page's 17.9 tonnes per terawatt-hour stays the page's round." in text
    assert "That figure is not 90 percent of 157 t." in text
    assert "The page's 80 percent stays that round." in text
    assert "The 90 percent load factor stays the page's round." in text
    assert "The example 190 tU is not a second target." in text
    assert "The example 147 tU stays the page's example." in text
    assert "The example 7.9 TWh stays the page's round." in text
    assert "The catalog ceiling stays unrun." in text
    assert "The catalog ceiling stays above the temperature Carnot." in text
    assert "The 6.5 percent page stays unrun." in text
    assert "The 60 GWd/t stack stays unrun." in text
    assert "The 55 GWd/t step stays unrun." in text
    assert "The 70 GWd/t step stays unrun." in text
    assert "The comparison column stays unrun." in text
    assert "The fabrication cap stays 5 percent." in text
    assert "The fast row's filed barrier is the sodium system." in text
    assert "The demonstrated plant stays at its filed efficiency." in text
    assert "This score adds no joules." in text
    assert "The published stamp records these closes." in text
    assert "The checked-in schedule is not rewritten." in text
    assert "The channel search stays closed." in text
    assert "This run does not install a seed as a new operating plant." in text
    assert "Both routes stay" not in text
    assert text.count("The unit dressing names it.") == 1
    assert text.count("The positive dressing names it.") == 16
    assert text.count("The unit-interval dressing is empty.") == 8
    agreements = 0
    for figure in report.figures:
        closing = figure.steps[-1]
        assert f"The dressing is `{closing.dressing.name}`" in text
        if f"{closing.dressing.value:.6g}" == f"{figure.target:.6g}":
            assert f"The closing value agrees with `{figure.target:.6g}` at six figures." in text
            agreements += 1
        else:
            assert f"The closing print is `{closing.dressing.value:.6g}`." in text
            assert f"The closing value agrees with `{figure.target:.6g}` at six figures." not in text
    assert text.count("The closing value agrees with") == agreements
    printed = {
        f"{value:.6g}"
        for value in (
            report.nuclear_scalar,
            report.thermodynamics_scalar,
            report.gamma_fsot,
            report.valve_beta,
            report.alpha,
            report.helium_mismatch,
        )
    }
    for banned in (
        report.nuclear_scalar * report.valve_beta,
        report.thermodynamics_scalar * report.gamma_fsot,
        1.0 / report.valve_beta,
        1.0 / report.alpha,
    ):
        for formatted in (f"{banned:.6g}", f"{banned:.6e}"):
            if formatted in printed:
                continue
            assert formatted not in text
    for needle in (
        "22.275",
        "3.328",
        "3.703",
        "warp_shell_energy_j",
        "137.035999177",
        "1.602176634",
        "8071.31713",
        "1.17302e-9",
    ):
        assert needle not in text
    output = default_output_path()
    assert output.parent.parent == project_root() / "runs"
    assert main() == 0
    assert output.read_text(encoding="utf-8") == text
    readme = (project_root() / "README.md").read_text(encoding="utf-8")
    quantum = readme.index("## Fission quantum")
    aside = readme.index("## Fission aside")
    tests = readme.index("## Tests")
    assert quantum < aside < tests
    assert "py -3.14 -m solcascade.production.fission_aside" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_aside.py"
    ).read_text(encoding="utf-8")
    assert "dialed" not in source
    assert "pi^-3*phi^5" not in source
    assert "warp_shell_energy_j" not in source
    assert "22.275" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9", "3.328", "3.703"):
        assert needle not in source
    observables = (
        project_root() / "src" / "solcascade" / "production" / "fission_observables.py"
    ).read_text(encoding="utf-8")
    assert "The channel search stays closed." in observables
    assert "dialed=bindings_inside and quantum_emitted" in observables
    package = (project_root() / "src" / "solcascade" / "production" / "__init__.py").read_text(
        encoding="utf-8"
    )
    assert "fission_aside" not in package
    assert "import" not in package
    lean = (project_root() / "lean" / "FissionStamp.lean").read_text(encoding="utf-8")
    assert "installedClosings.length = 96" in lean
    for closing in _CLOSINGS.values():
        assert closing in lean
    for name in (
        "theorem aside_compare_concentrate",
        "theorem aside_compare_hex",
        "theorem aside_compare_enriched_hex",
        "theorem aside_compare_oxide",
        "theorem aside_compare_per_twh",
        "theorem aside_compare_transuranic",
        "theorem aside_compare_used_oxide",
        "theorem aside_compare_burnup",
        "theorem aside_red_fraction",
        "theorem aside_red_tails",
        "theorem aside_red_tails_fraction",
        "theorem aside_red_high_share",
        "theorem aside_red_low_share",
        "theorem aside_candu_uranium",
        "theorem aside_candu_plutonium",
        "theorem aside_candu_plutonium_fraction",
        "theorem aside_candu_efficiency",
        "theorem aside_candu_burnup",
        "theorem aside_candu_per_twh",
    ):
        assert name in lean
    assert "fission_aside" not in lean
    assert "fission_column" not in lean
    assert "fission_reject" not in lean
    assert "fission_quantum" not in lean
    log = (project_root() / "data" / "fission_stamp.txt").read_text(encoding="utf-8").splitlines()
    assert log[-4:] == ["closes 96", "pairs 4", "miss 1", "fails 0"]
    unrun = [line for line in log if " unrun " in line]
    assert [line.split()[0] for line in unrun] == [
        "bounds_ceiling",
        "stack_atomic",
        "stack_burnup",
        "stack_stack",
        "stack_higher",
        "stack_electric",
        "aside_compare_concentrate",
        "aside_compare_hex",
        "aside_compare_enriched_hex",
        "aside_compare_oxide",
        "aside_compare_per_twh",
        "aside_compare_transuranic",
        "aside_compare_used_oxide",
    ]
    named = [line.split()[0] for line in log if " named " in line]
    for ident in (
        "aside_red_tails",
        "aside_red_demand",
        "aside_red_lower_demand",
        "aside_red_oxide",
        "aside_red_lower_oxide",
        "aside_candu_uranium",
        "aside_candu_plutonium",
        "aside_candu_efficiency",
        "aside_candu_demand",
        "aside_candu_capacity",
    ):
        assert ident in named
