"""Alpha dressings of the two decay limits left outside the complement."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_decay_dress import (
    build_decay_dress,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_decay_dress()


def test_both_limits_meet_the_review_bar(report):
    assert tuple(figure.label for figure in report.figures) == (
        "one-second limit",
        "one-day limit",
    )
    second, day = report.figures
    assert second.anchor_name == "pi^-5+phi^-6"
    assert second.dressing.name == "(pi^-5+phi^-6)*(1-alpha*phi^-3)"
    assert f"{second.dressing.relative_gap:.6g}" == "8.50946e-05"
    assert second.dressing.review is True
    assert second.inside_review == 7
    assert second.inside_tenth == 1
    assert day.anchor_name == "e^-3*phi^-4"
    assert day.dressing.name == "(e^-3*phi^-4)*(1+alpha*phi^-2)"
    assert f"{day.dressing.relative_gap:.6g}" == "3.03874e-05"
    assert day.dressing.review is True
    assert day.inside_review == 2
    assert day.inside_tenth == 1
    assert second.size.expressions == day.size.expressions == 116
    assert second.size.inside == day.size.inside == 116
    assert f"{second.size.median_spacing:.6g}" == "0.000511364"
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    for figure in report.figures:
        assert figure.dressing.relative_gap > report.helium_mismatch
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX


def test_report_leaves_the_dressings_uninstalled(report):
    text = render()
    assert text == render()
    assert "The two limits here are the figures the complement left outside the review bar." in text
    assert "No seed is installed as a decay-heat fraction." in text
    assert "Every scored gap stays wider than that match." in text
    assert "The dressing is not installed." in text
    assert "It is not a second electricity product." in text
    assert "The checked-in schedule is not rewritten." in text
    assert "(pi^-5+phi^-6)*(1-alpha*phi^-3)" in text
    assert "(e^-3*phi^-4)*(1+alpha*phi^-2)" in text
    assert "7 dressings of this anchor meet the review bar." in text
    assert "2 dressings of this anchor meet the review bar." in text
    for figure in report.figures:
        scaled = figure.dressing.value * report.second_watts
        assert f"{scaled:.6g}" not in text
        assert f"{scaled:.6e}" not in text
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
    outside = readme.index("## Fission decay outside")
    dress = readme.index("## Fission decay dress")
    tests = readme.index("## Tests")
    assert outside < dress < tests
    assert "py -3.14 -m solcascade.production.fission_decay_dress" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_decay_dress.py"
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
    assert "fission_decay_dress" not in package
