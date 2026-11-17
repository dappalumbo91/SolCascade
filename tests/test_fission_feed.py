"""The operating uranium column keeps its printed quotients."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_feed import build_feed, default_output_path, main, render
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_feed()


def _figure(report, label):
    return next(figure for figure in report.figures if figure.label == label)


def test_operating_quotients_keep_the_printed_masses(report):
    assert report.reload_tonnes == 24.3
    assert report.natural_t == 211.0
    assert report.concentrate_t == 249.0
    assert report.feed_hex_t == 312.0
    assert report.enriched_hex_t == 35.9
    assert report.fabricated_t == 27.6
    assert report.core_t == 73.0
    assert report.eighteen_month_t == 36.5
    assert report.transuranic_t == 0.280
    assert report.nameplate_twh == 8.76
    assert report.compare_natural_t == 163.0
    assert report.compare_product_t == 16.85
    assert report.compare_burnup_gwd == 65.0
    assert report.enriched_percent == 4.5
    assert report.tails_percent == 0.22
    assert report.discharge_gwd == 45.0
    assert report.higher_burnup_gwd == 60.0
    assert report.fabrication_cap_percent == 5.0
    assert report.built_efficiency == 0.36
    assert report.annual_swu == 182_000.0
    assert report.fast_ceiling == 0.5
    assert report.eighteen_month_t / report.core_t == report.fast_ceiling
    assert report.compare_burnup_gwd > report.higher_burnup_gwd > report.discharge_gwd
    assert report.compare_multiple > report.natural_t / report.reload_tonnes
    assert report.compare_product_t < report.reload_tonnes
    assert report.compare_natural_t < report.natural_t
    assert len(report.figures) == 8
    expected = {
        "feed multiple": (
            "(((pi^2-G^-2)*(1+alpha*phi^-5))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
            "8.68313",
            "8.68314",
            "1.00158e-06",
            "9.49642e-06",
            3,
            True,
        ),
        "natural uranium per terawatt-hour": (
            "((pi^3-phi^4)*(1-alpha*phi^-2))*(1+alpha*pi^-4)",
            "24.0868",
            "24.0867",
            "4.12304e-06",
            "2.98663e-05",
            2,
            False,
        ),
        "concentrate share": (
            "((e^-5+G^2)*(1+alpha*phi^-3))*(1+alpha*pi^-3)",
            "0.84739",
            "0.847387",
            "2.76003e-06",
            "3.17754e-05",
            2,
            False,
        ),
        "feed hexafluoride share": (
            "((pi^-3+G^5)*(1-alpha*phi^-4))*(1-alpha*pi^-6)",
            "0.676282",
            "0.676282",
            "2.53709e-07",
            "1.02443e-05",
            2,
            False,
        ),
        "enriched hexafluoride share": (
            "((pi^-3+G^5)*(1-alpha*pi^-3))*(1+alpha*e^-5)",
            "0.67688",
            "0.676882",
            "2.76098e-06",
            "2.2561e-05",
            2,
            True,
        ),
        "oxide share": (
            "((sqrt(e)-G^3)*(1+alpha*pi^-3))*(1-alpha*pi^-6)",
            "0.880435",
            "0.880433",
            "1.95969e-06",
            "1.24577e-05",
            2,
            True,
        ),
        "core fraction": (
            "(((gamma*gamma)*(1-alpha*e^-2))*(1+alpha*pi^-4))*(1+alpha*pi^-6)",
            "0.332877",
            "0.332876",
            "1.11034e-06",
            "9.38747e-06",
            3,
            True,
        ),
        "transuranic share": (
            "((pi*pi^-5)*(1-alpha*phi))*(1+alpha*e^-6)",
            "0.0101449",
            "0.010145",
            "2.36829e-06",
            "8.1259e-06",
            2,
            False,
        ),
    }
    for label, (name, target, closing, gap, second, steps, review) in expected.items():
        figure = _figure(report, label)
        assert figure.steps[-1].dressing.name == name
        assert f"{figure.target:.6g}" == target
        assert f"{figure.steps[-1].dressing.value:.6g}" == closing
        assert f"{figure.steps[-1].dressing.relative_gap:.6g}" == gap
        assert f"{figure.steps[-1].second_gap:.6g}" == second
        assert len(figure.steps) == steps
        assert figure.steps[-1].inside_helium == 1
        assert all(step.inside_helium == 0 for step in figure.steps[:-1])
        assert figure.steps[-1].expressions == 116
        assert f"{figure.steps[-1].median_spacing:.6g}" == "0.000511364"
        assert (figure.anchor.relative_gap <= 1.0e-3) is review
    assert f"{report.compare_multiple:.6g}" == "9.67359"
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_keeps_the_operating_column(report):
    text = render()
    assert text == render()
    assert "The operating column is 211 t of natural uranium for 24.3 t of enriched uranium at 4.5 percent, 0.22 percent tails, 45 GWd/t, and 33 percent." in text
    assert "The chain stays at k_eff of 1." in text
    assert "The feed multiple is `8.68313`." in text
    assert "The same quotient is `2110/243`." in text
    assert "The unit-interval dressing is empty." in text
    assert "The positive dressing names it." in text
    assert "The five-atom nearest stays inside the review bar." in text
    assert "The five-atom nearest is not the installed close." in text
    assert "The closing print is `8.68314`." in text
    assert "The plants-being-built efficiency leaves these masses in place." in text
    assert "Natural uranium over the nameplate year is `24.0868` t per terawatt-hour." in text
    assert "The same quotient is `5275/219`." in text
    assert "The page rounds that quotient to 24." in text
    assert "The rounded 24 stays the page's hence." in text
    assert "The same quotient is `211/249`." in text
    assert "The same quotient is `211/312`." in text
    assert "The closing value agrees with `0.676282` at six figures." in text
    assert text.count("The closing value agrees with") == 1
    assert "The same quotient is `243/359`." in text
    assert "The feed hexafluoride and the enriched hexafluoride stay different masses." in text
    assert "The same quotient is `81/92`." in text
    assert "The annual reload over that core is `0.332877`." in text
    assert "The same quotient is `243/730`." in text
    assert "The annual reload over the core stays under the burnup gain." in text
    assert "The prose one-third stays approximate." in text
    assert "The 18-month reload is half the core load." in text
    assert "That half keeps the fast ceiling." in text
    assert "The transuranic share of that used fuel is `0.0101449`." in text
    assert "The same quotient is `7/690`." in text
    assert "The unit dressing names it." in text
    assert "The used-fuel uranium oxide is printed as 26 t and does not close the fabricated mass." in text
    assert "The used-fuel fission-product mass stays the printed 1 tonne." in text
    assert "The fast row's filed barrier is the sodium system." in text
    assert "Its feed multiple is `9.67359`." in text
    assert "The same quotient is `3260/337`." in text
    assert "That multiple stays above the operating feed multiple." in text
    assert "That column sits on the fabrication cap." in text
    assert "The `65` GWd/t burn-up stays above the operating discharge." in text
    assert "The 60 GWd/t stack stays unrun." in text
    assert "This comparison stays unrun." in text
    assert "The enrichment note's about 140000 SWU stays approximate." in text
    assert "The mining line stays a range." in text
    assert "The fabrication cap stays 5 percent." in text
    assert "The catalog ceiling stays unrun." in text
    assert "This score adds no joules." in text
    assert "The published stamp records these closes." in text
    assert "The published stamp records the comparison multiple as a rational of the printed masses." in text
    assert "Every closing gap meets the helium quotient match." in text
    assert "This run does not install a seed as a new enrichment or a new plant row." in text
    assert "The checked-in schedule is not rewritten." in text
    output = default_output_path()
    assert output.parent.parent == project_root() / "runs"
    assert main() == 0
    assert output.read_text(encoding="utf-8") == text
    readme = (project_root() / "README.md").read_text(encoding="utf-8")
    lever = readme.index("## Fission lever")
    feed = readme.index("## Fission feed")
    tests = readme.index("## Tests")
    assert lever < feed < tests
    assert "py -3.14 -m solcascade.production.fission_feed" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_feed.py"
    ).read_text(encoding="utf-8")
    assert "phi^-5*gamma^2" not in source
    assert "phi^-6*G^-5" not in source
    assert "pi^-3*phi^5" not in source
    assert "warp_shell_energy_j" not in source
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
    assert "fission_feed" not in package
    assert "import" not in package
    lean = (project_root() / "lean" / "FissionStamp.lean").read_text(encoding="utf-8")
    assert "installedClosings.length = 96" in lean
    assert "theorem feed_multiple" in lean
    assert "theorem feed_per_twh" in lean
    assert "theorem feed_concentrate" in lean
    assert "theorem feed_hex" in lean
    assert "theorem feed_enriched_hex" in lean
    assert "theorem feed_oxide" in lean
    assert "theorem feed_core" in lean
    assert "theorem feed_transuranic" in lean
    assert "theorem feed_compare" in lean
    assert "theorem feed_compare_order" in lean
    assert "(((pi^2-G^-2)*(1+alpha*phi^-5))*(1-alpha*pi^-5))*(1-alpha*pi^-6)" in lean
    assert "((pi*pi^-5)*(1-alpha*phi))*(1+alpha*e^-6)" in lean
    assert "fission_feed" not in lean
    log = (project_root() / "data" / "fission_stamp.txt").read_text(encoding="utf-8").splitlines()
    assert log[-4:] == ["closes 96", "pairs 4", "miss 1", "fails 0"]
