"""The operating column leaves no further free quotient."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_column import (
    build_column,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_column()


def test_remaining_quotients_keep_the_named_scores(report):
    assert report.annual_swu == 182_000.0
    assert report.reload_tonnes == 24.3
    assert report.natural_t == 211.0
    assert report.nameplate_twh == 8.76
    assert report.discharge_gwd == 45.0
    assert report.cap_burnup_gwd == 55.0
    assert report.higher_burnup_gwd == 60.0
    assert report.compare_burnup_gwd == 65.0
    assert report.sighted_burnup_gwd == 70.0
    assert report.enriched_percent == 4.5
    assert report.tails_percent == 0.22
    assert report.fabrication_cap_percent == 5.0
    assert report.atomic_percent == 6.5
    assert report.sighted_percent == 6.0
    assert report.fast_ceiling == 0.5
    assert report.discharge_gwd < report.cap_burnup_gwd < report.higher_burnup_gwd
    assert report.higher_burnup_gwd < report.compare_burnup_gwd < report.sighted_burnup_gwd
    assert report.sighted_percent > report.fabrication_cap_percent
    assert f"{report.swu_per_twh:.6g}" == "20776.3"
    assert f"{report.enriched_per_twh:.6g}" == "2.77397"
    assert f"{report.product_share:.6g}" == "0.115166"
    assert f"{report.tails_stream:.6g}" == "0.884834"
    assert f"{report.swu_per_natural_t:.6g}" == "862.559"
    assert f"{report.fissile_part:.6g}" == "0.0049291"
    assert f"{report.cap_step:.6g}" == "1.22222"
    assert f"{report.built_step:.6g}" == "1.09091"
    assert f"{report.compare_multiple:.6g}" == "9.67359"
    assert f"{report.compare_product:.6g}" == "0.103374"
    assert f"{report.sighted_step:.6g}" == "1.55556"
    assert abs(report.product_share * (report.natural_t / report.reload_tonnes) - 1.0) < 1.0e-14
    assert abs(report.tails_stream - (1.0 - report.product_share)) / report.tails_stream < 1.0e-15
    assert abs(report.compare_product * report.compare_multiple - 1.0) < 1.0e-14
    assert abs(report.swu_per_twh - report.annual_swu / report.nameplate_twh) / report.swu_per_twh < 1.0e-15
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_records_the_identities(report):
    text = render()
    assert text == render()
    assert "The operating column leaves no further free quotient." in text
    assert "The chain stays at k_eff of 1." in text
    assert "The nameplate year is `8.76` TWh." in text
    assert "The annual separative work stays `182000` SWU." in text
    assert "The separative work over the nameplate year is `20776.3` SWU per terawatt-hour." in text
    assert "The same quotient is `4550000/219`." in text
    assert (
        "That quotient is the named separative work per kilogram times the enriched kilograms in that terawatt-hour."
        in text
    )
    assert "That quotient keeps the reload score and the feed score." in text
    assert "The enriched uranium in one nameplate terawatt-hour is `2.77397` t." in text
    assert "The same quotient is `405/146`." in text
    assert "That rate keeps the per-terawatt-hour quotient over the feed multiple." in text
    assert "The filed reload stays `24.3` t." in text
    assert "The product share of the natural uranium is `0.115166`." in text
    assert "The same quotient is `243/2110`." in text
    assert "The product share is the reciprocal of the feed multiple." in text
    assert "The product share keeps the feed score." in text
    assert "The tails stream is `0.884834` of the natural uranium." in text
    assert "The same quotient is `1867/2110`." in text
    assert "The tails stream is one minus the product share." in text
    assert "The tails stream keeps the feed score." in text
    assert "The page prints no tails mass." in text
    assert "The separative work per tonne of natural uranium is `862.559` SWU." in text
    assert "The same quotient is `182000/211`." in text
    assert "That quotient keeps the named separative work per kilogram and the feed multiple." in text
    assert "The fissile part of the natural uranium is `0.0049291`." in text
    assert "The same quotient is `26001/5275000`." in text
    assert "That quotient keeps the enriched fraction, the tails fraction, and the feed multiple." in text
    assert "The prose 0.49 percent stays the page's round." in text
    assert "The prose 0.394 percent stays outside the operating column." in text
    assert "The prose 0.6 percent utilization stays approximate." in text
    assert "The page's 55 GWd/t step over that discharge is `1.22222`." in text
    assert "The same quotient is `11/9`." in text
    assert "That burn-up sits on the fabrication cap and stays unrun." in text
    assert "Sixty over that 55 GWd/t stays the named built step `1.09091`." in text
    assert "The same quotient is `12/11`." in text
    assert "That quotient keeps the built score." in text
    assert "The 60 GWd/t stack stays unrun." in text
    assert "The comparison product share is `0.103374`." in text
    assert "The same quotient is `337/3260`." in text
    assert "The comparison product share is the reciprocal of the unrun multiple." in text
    assert "This comparison stays unrun." in text
    assert "The page's 70 GWd/t step over the operating discharge is `1.55556`." in text
    assert "The same quotient is `14/9`." in text
    assert "The page states `6` percent enrichment for that burn-up." in text
    assert "That enrichment stays above the fabrication cap." in text
    assert "That burn-up stays unrun." in text
    assert "The atomic percent stays above the fabrication cap." in text
    assert "The 18-month reload is half the core load." in text
    assert "That half keeps the fast ceiling." in text
    assert "The enrichment note's about 140000 SWU stays approximate." in text
    assert "The mining line stays a range." in text
    assert "The 90 percent load factor stays the page's round." in text
    assert "The Red Book tails assays are 0.30 percent and 0.25 percent." in text
    assert "Those assays stay off this column." in text
    assert "The CANDU row stays a different plant." in text
    assert "The fabrication cap stays 5 percent." in text
    assert "The demonstrated plant stays at its filed efficiency." in text
    assert "The catalog ceiling stays unrun." in text
    assert "The fast row's filed barrier is the sodium system." in text
    assert "This score adds no joules." in text
    assert "The published stamp records these quotients as rationals of the named scores." in text
    assert "This run does not install a seed as a new enrichment or a new plant row." in text
    assert "The checked-in schedule is not rewritten." in text
    assert "186.7" not in text
    assert "22.275" not in text
    printed = {
        f"{value:.6g}"
        for value in (
            report.swu_per_twh,
            report.product_share,
            report.tails_stream,
            report.swu_per_natural_t,
            report.enriched_per_twh,
            report.fissile_part,
            report.cap_step,
            report.sighted_step,
            report.built_step,
            report.compare_product,
            report.compare_multiple,
            report.nameplate_twh,
            report.annual_swu,
            report.reload_tonnes,
            report.natural_t,
            report.nuclear_scalar,
            report.thermodynamics_scalar,
            report.gamma_fsot,
            report.valve_beta,
        )
    }
    for banned in (
        report.nuclear_scalar * report.swu_per_twh,
        report.thermodynamics_scalar * report.product_share,
        report.gamma_fsot * report.enriched_per_twh,
        report.valve_beta * report.fissile_part,
        report.nuclear_scalar * report.tails_stream,
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
    feed = readme.index("## Fission feed")
    column = readme.index("## Fission column")
    tests = readme.index("## Tests")
    assert feed < column < tests
    assert "py -3.14 -m solcascade.production.fission_column" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_column.py"
    ).read_text(encoding="utf-8")
    assert "1-alpha" not in source
    assert "pi^-3*phi^5" not in source
    assert "warp_shell_energy_j" not in source
    assert "186.7" not in source
    assert "22.275" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9", "3.328", "3.703"):
        assert needle not in source
        assert needle not in text
    observables = (
        project_root() / "src" / "solcascade" / "production" / "fission_observables.py"
    ).read_text(encoding="utf-8")
    assert "The channel search stays closed." in observables
    package = (project_root() / "src" / "solcascade" / "production" / "__init__.py").read_text(
        encoding="utf-8"
    )
    assert "fission_column" not in package
    assert "import" not in package
    lean = (project_root() / "lean" / "FissionStamp.lean").read_text(encoding="utf-8")
    assert "installedClosings.length = 96" in lean
    assert "theorem column_swu_year" in lean
    assert "theorem column_product" in lean
    assert "theorem column_tails" in lean
    assert "theorem column_swu_natural" in lean
    assert "theorem column_enriched" in lean
    assert "theorem column_fissile" in lean
    assert "theorem column_cap_step" in lean
    assert "theorem column_cap_keeps_built" in lean
    assert "theorem column_compare_product" in lean
    assert "theorem column_sighted_step" in lean
    assert "fission_column" not in lean
    log = (project_root() / "data" / "fission_stamp.txt").read_text(encoding="utf-8").splitlines()
    assert log[-4:] == ["closes 96", "pairs 4", "miss 1", "fails 0"]
