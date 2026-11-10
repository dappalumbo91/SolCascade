"""Complements for the four decay figures left outside the review bar."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_decay_outside import (
    build_decay_outside,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_decay_outside()


def _figure(report, label):
    return next(figure for figure in report.figures if figure.label == label)


def test_two_complements_meet_the_review_bar(report):
    assert tuple(figure.label for figure in report.figures) == (
        "one-second limit",
        "one-day limit",
        "two-day limit",
        "two-day drop",
    )
    assert report.complement_size.expressions == 2274
    assert report.complement_size.inside == 1395
    assert f"{report.complement_size.median_spacing:.6g}" == "0.000333282"
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    second = _figure(report, "one-second limit")
    day = _figure(report, "one-day limit")
    two = _figure(report, "two-day limit")
    drop = _figure(report, "two-day drop")
    assert second.five.name == "pi^-5+phi^-6"
    assert second.complement.name == "1-(pi^-2+G^2)"
    assert second.complement.review is False
    assert day.five.name == "e^-3*phi^-4"
    assert day.complement.review is False
    assert day.alpha_gap > 1.0e-3
    assert f"{day.alpha_gap:.6g}" == "0.00185228"
    assert two.complement.name == "1-(gamma^6+sqrt(G))"
    assert two.complement.review is True
    assert f"{two.complement.relative_gap:.6g}" == "2.46138e-05"
    assert drop.complement.name == "1-(gamma^3/sqrt(G))"
    assert drop.complement.review is True
    assert f"{drop.complement.relative_gap:.6g}" == "2.538e-05"
    for figure in report.figures:
        assert figure.five.review is False
        assert figure.complement.relative_gap > report.helium_mismatch
        assert figure.alpha_gap > report.helium_mismatch
        assert figure.alpha_gap > 1.0e-3


def test_watts_stay_on_the_1952_fit(report):
    assert report.second_watts > report.day_watts > report.two_watts
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX


def test_report_leaves_the_complements_uninstalled(report):
    text = render()
    assert text == render()
    assert "These are the four decay figures the five-atom grammar left outside the review bar." in text
    assert "No seed is installed as a decay-heat fraction." in text
    assert "Every scored gap stays wider than that match." in text
    assert "It is not a second electricity product." in text
    assert "The checked-in schedule is not rewritten." in text
    assert "1-(gamma^6+sqrt(G))" in text
    assert "1-(gamma^3/sqrt(G))" in text
    for figure in report.figures:
        scaled = figure.complement.value * report.second_watts
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
    seed = readme.index("## Fission decay seed")
    outside = readme.index("## Fission decay outside")
    tests = readme.index("## Tests")
    assert seed < outside < tests
    assert "py -3.14 -m solcascade.production.fission_decay_outside" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_decay_outside.py"
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
    assert "fission_decay_outside" not in package
