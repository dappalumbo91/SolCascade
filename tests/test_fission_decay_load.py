"""FSOT fractions times the filed operating heat are the cooling loads."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_decay_load import (
    build_decay_load,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_decay_load()


def test_cooling_loads_are_the_fsot_fractions(report):
    assert tuple(load.label for load in report.loads) == (
        "one-second fraction",
        "one-day fraction",
        "two-day fraction",
    )
    expected = (
        (
            "(phi^-6/sqrt(G))*(1-alpha*e^-4)",
            "3.97292e-06",
            "1.76426e+08",
            "1.76425e+08",
            "5.29278e+07",
            "5.29276e+07",
            "(((pi^-5+phi^-6)*(1-alpha*phi^-3))*(1-alpha*pi^-4))*(1-alpha*pi^-6)",
            "1.78453e+08",
            "1.78452e+08",
            "2.58319e-06",
        ),
        (
            "((pi^-4*G^5)*(1-alpha*gamma^5))*(1-alpha*pi^-6)",
            "3.34448e-06",
            "2.00482e+07",
            "2.00482e+07",
            "6.01447e+06",
            "6.01445e+06",
            "(((e^-3*phi^-4)*(1+alpha*phi^-2))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
            "2.20723e+07",
            "2.20723e+07",
            "1.0498e-06",
        ),
        (
            "((pi^-5*phi)*(1-alpha*e^-4))*(1-alpha*pi^-6)",
            "3.83823e-07",
            "1.602e+07",
            "1.602e+07",
            "4.80601e+06",
            "4.806e+06",
            "(1-(gamma^6+sqrt(G)))*(1-alpha*pi^-5)",
            "1.80415e+07",
            "1.80415e+07",
            "7.67233e-07",
        ),
    )
    for load, row in zip(report.loads, expected, strict=True):
        (
            name,
            gap,
            watts,
            fit_watts,
            small,
            fit_small,
            limit_name,
            limit_watts,
            fit_limit,
            limit_gap,
        ) = row
        assert load.fraction.name == name
        assert f"{load.fraction.relative_gap:.6g}" == gap
        assert f"{load.watts:.6g}" == watts
        assert f"{load.fit_watts:.6g}" == fit_watts
        assert f"{load.small_watts:.6g}" == small
        assert f"{load.fit_small_watts:.6g}" == fit_small
        assert load.limit.name == limit_name
        assert f"{load.limit_watts:.6g}" == limit_watts
        assert f"{load.fit_limit_watts:.6g}" == fit_limit
        assert f"{load.limit.relative_gap:.6g}" == limit_gap
        assert load.watts == load.fraction.value * report.thermal_w
        assert load.small_watts == load.fraction.value * report.small_thermal_w
        assert load.fraction.relative_gap <= report.helium_mismatch
        assert load.limit.relative_gap <= report.helium_mismatch
    assert f"{report.day_over_second:.6g}" == "0.113635"
    assert f"{report.two_over_day:.6g}" == "0.799074"
    assert f"{report.day_ratio_gap:.6g}" == "6.28438e-07"
    assert f"{report.two_ratio_gap:.6g}" == "2.96064e-06"
    assert report.day_ratio_gap <= report.helium_mismatch
    assert report.two_ratio_gap <= report.helium_mismatch
    assert f"{report.thermal_w:.6g}" == "3.0303e+09"
    assert f"{report.small_thermal_w:.6g}" == "9.09091e+08"
    assert f"{report.alpha:.6g}" == "0.00729735"
    assert f"{report.helium_mismatch:.6g}" == "4.23431e-06"
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX


def test_report_keeps_the_fit_as_the_witness(report):
    text = render()
    assert text == render()
    assert "The cooling load is the helium-tight FSOT fraction times the filed operating heat." in text
    assert "The 1952 fit remains the calorimetry witness." in text
    assert "The derived loads are the cooling loads." in text
    assert "It is not a second electricity product." in text
    assert "The checked-in schedule is not rewritten." in text
    assert "gamma_FSOT does not multiply the plant." in text
    printed = set()
    for load in report.loads:
        for value in (
            load.watts,
            load.fit_watts,
            load.small_watts,
            load.fit_small_watts,
            load.limit_watts,
            load.fit_limit_watts,
        ):
            printed.add(f"{value:.6g}")
            assert f"{value:.6g}" in text
    for banned in (
        report.nuclear_scalar * report.thermal_w,
        report.thermodynamics_scalar * report.thermal_w,
        report.gamma_fsot * report.thermal_w,
        report.valve_beta * report.thermal_w,
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
    opened = readme.index("## Fission decay open")
    load = readme.index("## Fission decay load")
    tests = readme.index("## Tests")
    assert opened < load < tests
    assert "py -3.14 -m solcascade.production.fission_decay_load" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_decay_load.py"
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
    assert "fission_decay_load" not in package
