"""The 1952 decay-heat fit on the filed full-power year."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_decay import (
    build_decay,
    decay_fraction,
    default_output_path,
    limit_fraction,
    main,
    render,
    stated_band,
)
from solcascade.production.fission_process import FULL_POWER_HOURS, NAMEPLATE_MWE
from solcascade.production.fission_safety import SMALL_CLASS_MWE
from solcascade.production.refine import project_root


def _brace(time_s: float) -> float:
    return (time_s + 10.0) ** -0.2 - 0.87 * (time_s + 2.0e7) ** -0.2


def test_equation_two_is_the_1952_brace():
    operation_s = FULL_POWER_HOURS * 3600.0
    for cooling_s in (1.0, 86_400.0, 172_800.0):
        expected = 0.1 * (_brace(cooling_s) - _brace(cooling_s + operation_s))
        assert decay_fraction(cooling_s, operation_s) == pytest.approx(expected)
        assert limit_fraction(cooling_s) == pytest.approx(0.1 * _brace(cooling_s))
        assert decay_fraction(cooling_s, 0.0) == pytest.approx(0.0)
    assert stated_band(1.0) == pytest.approx(0.50)
    assert stated_band(86_400.0) == pytest.approx(0.10)
    assert stated_band(172_800.0) == pytest.approx(0.10)
    with pytest.raises(RuntimeError):
        decay_fraction(0.5, operation_s)


@pytest.fixture(scope="module")
def report():
    return build_decay()


def test_loads_fall_and_stay_on_the_heat_ledger(report):
    second, day, two = report.points
    thermal = NAMEPLATE_MWE * 1.0e6 / 0.33
    small = SMALL_CLASS_MWE * 1.0e6 / 0.33
    assert report.thermal_w == pytest.approx(thermal)
    assert report.small_thermal_w == pytest.approx(small)
    assert report.operation_s == pytest.approx(FULL_POWER_HOURS * 3600.0)
    previous = 1.0
    for point in report.points:
        assert point.year_fraction < point.limit_fraction < previous
        assert point.year_watts == pytest.approx(point.year_fraction * thermal)
        assert point.small_watts == pytest.approx(point.year_watts * SMALL_CLASS_MWE / NAMEPLATE_MWE)
        assert point.year_seed.name == "gamma/pi"
        assert point.limit_seed.name == "gamma/pi"
        assert point.year_seed.relative_gap > 1.0e-3
        previous = point.year_fraction
    assert second.band == pytest.approx(0.50)
    assert day.band == pytest.approx(0.10)
    assert two.band == pytest.approx(0.10)
    assert report.day_over_second == pytest.approx(day.year_fraction / second.year_fraction)
    assert report.two_over_day == pytest.approx(two.year_fraction / day.year_fraction)
    assert report.day_seed.name == "gamma/pi"
    assert report.two_seed.name == "G^2"
    assert report.two_seed.relative_gap > 1.0e-3
    assert f"{report.two_seed.relative_gap:.6g}" == "0.0499534"
    assert report.delayed_gap == pytest.approx(abs(day.year_fraction - 0.0066) / 0.0066)
    assert report.delayed_gap > 1.0e-3
    assert f"{report.delayed_gap:.6g}" == "0.00240798"
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX


def test_report_keeps_decay_heat_off_the_electricity(report):
    text = render()
    assert text == render()
    assert "The cooling page states no decay-heat fraction." in text
    assert "It is not a second electricity product." in text
    assert "The delayed-neutron margin stays in the core." in text
    assert "The fit is scored as published, with U-239 and Np-239 left in." in text
    assert "No seed is installed as a decay-heat fraction." in text
    assert "The one-day heat ratio stays off that neutron count." in text
    assert "The closest of these seeds is still wider than one part in a thousand." in text
    assert "The checked-in schedule is not rewritten." in text
    assert "gamma/pi" in text
    assert "G^2" in text
    for scaled in (
        report.nuclear_scalar * report.points[0].year_watts,
        report.gamma_fsot * report.points[0].small_watts,
        report.two_seed.value * report.thermal_w,
        report.points[1].year_watts * 0.33,
    ):
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
    safety = readme.index("## Fission safety")
    decay = readme.index("## Fission decay")
    tests = readme.index("## Tests")
    assert safety < decay < tests
    assert "py -3.14 -m solcascade.production.fission_decay" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_decay.py"
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
    assert "fission_decay" not in package
