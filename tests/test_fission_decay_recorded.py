"""Base-piece dressings of the six decay names the first dressing left untouched."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_decay_recorded import (
    build_decay_recorded,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_decay_recorded()


def test_three_recorded_nearests_meet_helium(report):
    assert tuple(figure.label for figure in report.figures) == (
        "one-second fraction",
        "one-day fraction",
        "two-day fraction",
        "two-day limit",
        "one-day drop",
        "two-day drop",
    )
    expected = (
        ("phi^-6/sqrt(G)", "(phi^-6/sqrt(G))*(1-alpha*e^-4)", "3.97292e-06", 30, 4, 1, 116, "0.000511364"),
        ("pi^-4*G^5", "(pi^-4*G^5)*(1-alpha*gamma^5)", "1.0935e-05", 28, 2, 0, 116, "0.000511364"),
        ("pi^-5*phi", "(pi^-5*phi)*(1-alpha*e^-4)", "7.97431e-06", 30, 4, 0, 116, "0.000511364"),
        (
            "1-(gamma^6+sqrt(G))",
            "(1-(gamma^6+sqrt(G)))*(1-alpha*pi^-5)",
            "7.67233e-07",
            29,
            10,
            1,
            116,
            "0.000511364",
        ),
        ("gamma^3*G^6", "(gamma^3*G^6)*(1+alpha*gamma^5)", "4.91253e-05", 28, 1, 0, 116, "0.000511364"),
        (
            "1-(gamma^3/sqrt(G))",
            "(1-(gamma^3/sqrt(G)))*(1+alpha*pi^-5)",
            "1.53457e-06",
            29,
            9,
            1,
            110,
            "0.000433698",
        ),
    )
    for figure, row in zip(report.figures, expected, strict=True):
        anchor, name, gap, review, tenth, helium, expressions, spacing = row
        assert figure.anchor.name == anchor
        assert figure.anchor.review is True
        assert figure.dressing.name == name
        assert f"{figure.dressing.relative_gap:.6g}" == gap
        assert figure.dressing.relative_gap < figure.anchor.relative_gap
        assert figure.inside_review == review
        assert figure.inside_tenth == tenth
        assert figure.inside_helium == helium
        assert figure.size.expressions == figure.size.inside == expressions
        assert f"{figure.size.median_spacing:.6g}" == spacing
        assert figure.dressing.review is True
    helium_labels = tuple(
        figure.label
        for figure in report.figures
        if figure.dressing.relative_gap <= report.helium_mismatch
    )
    assert helium_labels == (
        "one-second fraction",
        "two-day limit",
        "two-day drop",
    )
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX


def test_report_leaves_the_recorded_dressings_uninstalled(report):
    text = render()
    assert text == render()
    assert "These six figures are the decay names the first dressing left untouched." in text
    assert "The two limits already dressed stay on that score." in text
    assert "3 nearest dressings meet the helium quotient match." in text
    assert "3 nearest dressings stay wider than that match." in text
    assert "No seed is installed as a decay-heat fraction." in text
    assert "The dressing is not installed." in text
    assert "It is not a second electricity product." in text
    assert "The checked-in schedule is not rewritten." in text
    assert "`e^-4` dresses the one-second fraction and the two-day fraction with the same sign." in text
    assert "`gamma^5` dresses the one-day fraction and the one-day drop with opposite signs." in text
    assert "`pi^-5` dresses the two-day limit and the two-day drop with opposite signs." in text
    assert "One dressing of this anchor meets the helium quotient match." in text
    assert "No dressing of this anchor meets the helium quotient match." in text
    assert "One dressing of this anchor is inside one part in ten thousand." in text
    printed_loads = {
        f"{watts:.6g}"
        for watts in (report.second_watts, report.day_watts, report.two_watts)
    }
    for figure in report.figures:
        assert figure.dressing.name in text
        for watts in (report.second_watts, report.day_watts, report.two_watts):
            scaled = figure.dressing.value * watts
            for formatted in (f"{scaled:.6g}", f"{scaled:.6e}"):
                if formatted in printed_loads:
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
    dress = readme.index("## Fission decay dress")
    recorded = readme.index("## Fission decay recorded")
    tests = readme.index("## Tests")
    assert dress < recorded < tests
    assert "py -3.14 -m solcascade.production.fission_decay_recorded" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_decay_recorded.py"
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
    assert "fission_decay_recorded" not in package
