"""Added base-piece factors for the decay names still outside the helium match."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_decay_open import (
    build_decay_open,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_decay_open()


def test_every_open_figure_closes_inside_helium(report):
    assert tuple(figure.label for figure in report.figures) == (
        "one-second limit",
        "one-day fraction",
        "one-day limit",
        "two-day fraction",
        "one-day drop",
    )
    expected = {
        "one-second limit": (
            "(pi^-5+phi^-6)*(1-alpha*phi^-3)",
            "8.50946e-05",
            (
                (
                    "((pi^-5+phi^-6)*(1-alpha*phi^-3))*(1-alpha*pi^-4)",
                    "1.01737e-05",
                    0,
                    30,
                    7,
                ),
                (
                    "(((pi^-5+phi^-6)*(1-alpha*phi^-3))*(1-alpha*pi^-4))*(1-alpha*pi^-6)",
                    "2.58319e-06",
                    1,
                    30,
                    10,
                ),
            ),
        ),
        "one-day fraction": (
            "(pi^-4*G^5)*(1-alpha*gamma^5)",
            "1.0935e-05",
            (
                (
                    "((pi^-4*G^5)*(1-alpha*gamma^5))*(1-alpha*pi^-6)",
                    "3.34448e-06",
                    1,
                    30,
                    10,
                ),
            ),
        ),
        "one-day limit": (
            "(e^-3*phi^-4)*(1+alpha*phi^-2)",
            "3.03874e-05",
            (
                (
                    "((e^-3*phi^-4)*(1+alpha*phi^-2))*(1-alpha*pi^-5)",
                    "6.54067e-06",
                    0,
                    29,
                    9,
                ),
                (
                    "(((e^-3*phi^-4)*(1+alpha*phi^-2))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
                    "1.0498e-06",
                    1,
                    30,
                    10,
                ),
            ),
        ),
        "two-day fraction": (
            "(pi^-5*phi)*(1-alpha*e^-4)",
            "7.97431e-06",
            (
                (
                    "((pi^-5*phi)*(1-alpha*e^-4))*(1-alpha*pi^-6)",
                    "3.83823e-07",
                    1,
                    30,
                    10,
                ),
            ),
        ),
        "one-day drop": (
            "(gamma^3*G^6)*(1+alpha*gamma^5)",
            "4.91253e-05",
            (
                (
                    "((gamma^3*G^6)*(1+alpha*gamma^5))*(1+alpha*e^-5)",
                    "4.14453e-08",
                    1,
                    29,
                    10,
                ),
            ),
        ),
    }
    for figure in report.figures:
        start_name, start_gap, steps = expected[figure.label]
        assert figure.start_name == start_name
        assert f"{figure.start_gap:.6g}" == start_gap
        assert figure.start_gap > report.helium_mismatch
        assert len(figure.steps) == len(steps)
        for step, row in zip(figure.steps, steps, strict=True):
            name, gap, helium_count, review, tenth = row
            assert step.dressing.name == name
            assert f"{step.dressing.relative_gap:.6g}" == gap
            assert step.inside_helium == helium_count
            assert step.inside_review == review
            assert step.inside_tenth == tenth
            assert step.size.expressions == step.size.inside == 116
            assert f"{step.size.median_spacing:.6g}" == "0.000511364"
            assert step.size.median_spacing > report.helium_mismatch
        closing = figure.steps[-1]
        assert closing.dressing.relative_gap <= report.helium_mismatch
        assert closing.inside_helium == 1
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX


def test_report_leaves_the_open_dressings_uninstalled(report):
    text = render()
    assert text == render()
    assert "These five figures were still wider than the helium quotient match." in text
    assert "Every closing gap meets the helium quotient match." in text
    assert "No seed is installed as a decay-heat fraction." in text
    assert "The dressing is not installed." in text
    assert "It is not a second electricity product." in text
    assert "The checked-in schedule is not rewritten." in text
    assert (
        "`pi^-6` closes the one-second limit, the one-day fraction, "
        "the one-day limit, and the two-day fraction with a minus sign."
    ) in text
    assert "`e^-5` closes the one-day drop with a plus sign." in text
    printed_loads = {
        f"{watts:.6g}"
        for watts in (report.second_watts, report.day_watts, report.two_watts)
    }
    for figure in report.figures:
        for step in figure.steps:
            assert step.dressing.name in text
            for watts in (report.second_watts, report.day_watts, report.two_watts):
                scaled = step.dressing.value * watts
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
    recorded = readme.index("## Fission decay recorded")
    opened = readme.index("## Fission decay open")
    tests = readme.index("## Tests")
    assert recorded < opened < tests
    assert "py -3.14 -m solcascade.production.fission_decay_open" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_decay_open.py"
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
    assert "fission_decay_open" not in package
