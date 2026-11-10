"""Depth-two names for the decay-heat figures. None are installed."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_decay import build_decay
from solcascade.production.fission_decay_seed import (
    build_decay_seed,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_decay_seed()


def _figure(report, label):
    return next(figure for figure in report.figures if figure.label == label)


def test_depth_two_tightens_without_reaching_helium(report):
    assert report.ladder_size.expressions == 1224
    assert report.ladder_size.inside == 180
    assert report.five_size.expressions == 7380
    assert report.five_size.inside == 1395
    assert f"{report.ladder_size.median_spacing:.6g}" == "0.032156"
    assert f"{report.five_size.median_spacing:.6g}" == "0.00349101"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    second = _figure(report, "one-second fraction")
    day = _figure(report, "one-day fraction")
    two = _figure(report, "two-day fraction")
    drop = _figure(report, "one-day drop")
    later = _figure(report, "two-day drop")
    assert second.ladder.name == "e^-6+phi^-6"
    assert second.ladder.review is True
    assert second.five.name == "phi^-6/sqrt(G)"
    assert second.five.review is True
    assert day.five.name == "pi^-4*G^5"
    assert day.five.review is True
    assert day.ladder.review is False
    assert two.five.name == "pi^-5*phi"
    assert two.five.review is True
    assert two.ladder.review is False
    assert drop.five.name == "gamma^3*G^6"
    assert drop.five.review is True
    assert later.five.name == "pi^-3+G^3"
    assert later.five.review is False
    assert later.shallow.name == "G^2"
    for figure in report.figures:
        assert figure.five.relative_gap < figure.shallow.relative_gap
        assert figure.five.relative_gap <= figure.ladder.relative_gap
        assert figure.ladder.relative_gap > report.helium_mismatch
        assert figure.five.relative_gap > report.helium_mismatch
        assert figure.shallow.review is False


def test_watts_stay_on_the_1952_fit(report):
    decay = build_decay()
    second, day, two = decay.points
    assert report.second_watts == second.year_watts
    assert report.day_watts == day.year_watts
    assert report.two_watts == two.year_watts
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX


def test_report_leaves_the_expressions_uninstalled(report):
    text = render()
    assert text == render()
    assert "The decay watts stay the 1952 fit." in text
    assert "No seed is installed as a decay-heat fraction." in text
    assert "Every scored gap stays wider than that match." in text
    assert "It is not a second electricity product." in text
    assert "The checked-in schedule is not rewritten." in text
    assert "e^-6+phi^-6" in text
    assert "pi^-5*phi" in text
    assert "pi^-4*G^5" in text
    for figure in report.figures:
        for hit in (figure.ladder, figure.five):
            scaled = hit.value * report.second_watts
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
    decay = readme.index("## Fission decay")
    seed = readme.index("## Fission decay seed")
    tests = readme.index("## Tests")
    assert decay < seed < tests
    assert "py -3.14 -m solcascade.production.fission_decay_seed" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_decay_seed.py"
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
    assert "fission_decay_seed" not in package
