"""The 0.22 percent tails assumption and its fraction over 100."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_tails import (
    build_tails,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_tails()


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


def test_tails_meet_helium(report):
    helium = report.helium_mismatch
    assert report.tails_percent == 0.22
    assert report.percent.target == 0.22
    assert report.fraction.target == report.tails_percent / 100.0
    assert f"{report.fraction.target:.6g}" == "0.0022"
    assert report.enriched_percent == 4.5
    assert report.annual_swu == 182_000.0
    assert report.reload_tonnes == 24.3
    assert report.fabrication_cap_percent == 5.0
    assert report.tails_percent < report.enriched_percent < report.fabrication_cap_percent
    _assert_figure(
        report.percent,
        "phi^-2*gamma",
        "0.220477",
        "0.00216711",
        (
            (
                "(phi^-2*gamma)*(1-alpha*gamma^2)",
                "0.219941",
                "0.000269471",
                "0.00044071",
                (4, 0, 0),
            ),
            (
                "((phi^-2*gamma)*(1-alpha*gamma^2))*(1+alpha*gamma^6)",
                "0.22",
                "3.50354e-07",
                "3.4184e-05",
                (28, 3, 1),
            ),
        ),
        helium,
    )
    assert len(report.percent.steps) == 2
    assert f"{report.percent.steps[-1].dressing.value:.6g}" == "0.22"
    _assert_figure(
        report.fraction,
        "pi^-5-pi^-6",
        "0.0022276",
        "0.0125464",
        (
            (
                "(pi^-5-pi^-6)*(1-alpha*G^-6)",
                "0.00220008",
                "3.50579e-05",
                "0.00036419",
                (4, 1, 0),
            ),
            (
                "((pi^-5-pi^-6)*(1-alpha*G^-6))*(1-alpha*pi^-5)",
                "0.00220002",
                "1.1211e-05",
                "1.4113e-05",
                (29, 10, 0),
            ),
            (
                "(((pi^-5-pi^-6)*(1-alpha*G^-6))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
                "0.00220001",
                "3.62053e-06",
                "6.87749e-06",
                (30, 10, 1),
            ),
        ),
        helium,
    )
    assert len(report.fraction.steps) == 3
    assert f"{report.fraction.steps[-1].dressing.value:.6g}" == "0.00220001"
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX == report.digest_prefix


def test_report_names_the_tails(report):
    text = render()
    assert text == render()
    assert "The material balance is https://world-nuclear.org/information-library/nuclear-fuel-cycle/introduction/nuclear-fuel-cycle-overview." in text
    assert "0.22 percent tails" in text
    assert "FSOT names the tails percent and the tails fraction." in text
    assert "The chain stays at k_eff of 1." in text
    assert "The unit-interval dressing of the tails percent has `115` values." in text
    assert "The positive dressing names it." in text
    assert "The tails percent is `0.22`." in text
    assert "The five-atom nearest is `phi^-2*gamma` at `0.220477`." in text
    assert "The dressing is `((phi^-2*gamma)*(1-alpha*gamma^2))*(1+alpha*gamma^6)` at `0.22`." in text
    assert "The relative gap is `3.50354e-07`." in text
    assert "The closing value agrees with `0.22` at six figures." in text
    assert "The operating percent stays the printed tails percent." in text
    assert "The overview prints no further digit on that percent." in text
    assert "The tails fraction is the tails percent over 100." in text
    assert "The unit dressing names it." in text
    assert "The tails fraction is `0.0022`." in text
    assert "The dressing is `(((pi^-5-pi^-6)*(1-alpha*G^-6))*(1-alpha*pi^-5))*(1-alpha*pi^-6)` at `0.00220001`." in text
    assert "The relative gap is `3.62053e-06`." in text
    assert "The closing value agrees with `0.0022` at six figures." not in text
    assert "The closing print is `0.00220001`." in text
    assert "The operating fraction stays the tails percent over 100." in text
    assert "The enriched percent stays `4.5`." in text
    assert "The enriched percent stays the stack score." in text
    assert "The enriched fraction stays the held pair." in text
    assert "The annual separative work stays `182000` SWU." in text
    assert "The reload stays `24.3` t." in text
    assert "Those figures stay the process score." in text
    assert "The fabrication cap stays 5 percent." in text
    assert "The fabrication page's about 4.8 percent stays approximate." in text
    assert "This score adds no joules." in text
    assert "The catalog ceiling stays unrun." in text
    assert "The 60 GWd/t stack stays unrun." in text
    assert "The fast row's filed barrier is the sodium system." in text
    assert "The published stamp records these two closes." in text
    assert "Every closing gap meets the helium quotient match." in text
    assert "This run does not install a seed as a new enrichment or a new plant row." in text
    assert "The checked-in schedule is not rewritten." in text
    printed = {
        f"{report.tails_percent:.6g}",
        f"{report.enriched_percent:.6g}",
        f"{report.annual_swu:.6g}",
        f"{report.reload_tonnes:.6g}",
        f"{report.fabrication_cap_percent:.6g}",
        f"{report.alpha:.6g}",
        f"{report.helium_mismatch:.6g}",
        f"{report.nuclear_scalar:.6g}",
        f"{report.thermodynamics_scalar:.6g}",
        f"{report.gamma_fsot:.6g}",
        f"{report.valve_beta:.6g}",
    }
    for figure in (report.percent, report.fraction):
        printed.add(f"{figure.target:.6g}")
        printed.add(f"{figure.anchor.value:.6g}")
        printed.add(f"{figure.anchor.relative_gap:.6g}")
        for step in figure.steps:
            printed.add(f"{step.dressing.value:.6g}")
            printed.add(f"{step.dressing.relative_gap:.6g}")
            printed.add(f"{step.second_gap:.6g}")
    for banned in (
        report.nuclear_scalar * report.percent.target,
        report.thermodynamics_scalar * report.percent.target,
        report.gamma_fsot * report.fraction.target,
        report.valve_beta * report.fraction.target,
        report.nuclear_scalar * report.annual_swu,
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
    stamp = readme.index("## Fission stamp")
    tails = readme.index("## Fission tails")
    tests = readme.index("## Tests")
    assert stamp < tails < tests
    assert "py -3.14 -m solcascade.production.fission_tails" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_tails.py"
    ).read_text(encoding="utf-8")
    assert "warp_shell_energy_j" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9"):
        assert needle not in source
        assert needle not in text
    observables = (
        project_root() / "src" / "solcascade" / "production" / "fission_observables.py"
    ).read_text(encoding="utf-8")
    assert "The channel search stays closed." in observables
    package = (project_root() / "src" / "solcascade" / "production" / "__init__.py").read_text(
        encoding="utf-8"
    )
    assert "fission_tails" not in package
    assert "import" not in package
    lean = (project_root() / "lean" / "FissionStamp.lean").read_text(encoding="utf-8")
    assert "installedClosings.length = 70" in lean
    assert "((phi^-2*gamma)*(1-alpha*gamma^2))*(1+alpha*gamma^6)" in lean
    assert "(((pi^-5-pi^-6)*(1-alpha*G^-6))*(1-alpha*pi^-5))*(1-alpha*pi^-6)" in lean
    assert "fission_tails" not in lean
    log = (project_root() / "data" / "fission_stamp.txt").read_text(encoding="utf-8").splitlines()
    assert log[-4:] == ["closes 70", "pairs 4", "miss 1", "fails 0"]
