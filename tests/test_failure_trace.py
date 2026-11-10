"""Failure trace: beta = 1 ledgers, shell rise, and the FSOT record."""

import math

import pytest

from solcascade.constants import BOSCH_HALE_DT_T_MIN_KEV
from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.refine import SUPPLIES, project_root
from solcascade.production.trace import build, default_output_path, main, render


@pytest.fixture(scope="module")
def trace():
    return build()


def test_cooled_trajectory_sits_at_beta_one_inside_the_fit(trace):
    checked = trace.checked
    assert checked.fusion_j == pytest.approx(0.4977033387900156, abs=1.0e-12)
    assert checked.peak_beta > 1.0
    assert trace.cooled.temperature0_kev == pytest.approx(
        checked.temperature0_kev / checked.peak_beta, rel=1.0e-12
    )
    assert trace.cooled.temperature0_kev >= BOSCH_HALE_DT_T_MIN_KEV
    assert trace.cooled.peak_beta == pytest.approx(1.0, rel=1.0e-9)
    assert trace.cooled.fusion_j < checked.fusion_j
    assert trace.cooled_electric_j < dict(SUPPLIES)["same_b_helmholtz"]
    assert trace.cooled_efold_s > trace.cooled.walk_time_s


def test_density_limit_cuts_fusion_and_the_schedule_stays_under_the_draw(trace):
    helmholtz = dict(SUPPLIES)["same_b_helmholtz"]
    assert trace.store.practical_j > trace.store.helmholtz_stored_j * trace.store.recovery
    assert trace.store.practical_j == pytest.approx(
        trace.store.helmholtz_ledger_j + (1.0 - trace.store.recovery) * trace.store.helmholtz_stored_j
    )
    assert trace.store.unrecovered_store_j == pytest.approx(
        (1.0 - trace.store.recovery) * trace.store.helmholtz_stored_j
    )
    assert trace.store.helmholtz_stored_j > trace.store.helmholtz_ledger_j
    checked_limited = trace.limited[0]
    assert checked_limited.peak_temperature_kev == pytest.approx(trace.checked.peak_temperature_kev)
    assert checked_limited.fusion_j == pytest.approx(
        trace.checked.fusion_j / trace.checked.peak_beta ** 2, rel=1.0e-12
    )
    for row in trace.limited:
        assert row.raw_beta > 1.0
        assert row.density_scale == pytest.approx(1.0 / row.raw_beta)
        assert row.electric_j < helmholtz
        assert row.electric_j < trace.store.practical_j
        assert row.walk_yield_j < row.fusion_j
        assert row.efold_s > row.walk_time_s
        assert row.complete_electric_j == pytest.approx(row.ceiling_j * trace.operating_per_j)
        assert row.complete_electric_j > row.electric_j
    best = max(trace.limited, key=lambda row: row.fusion_j)
    widest = max(trace.limited, key=lambda row: row.complete_electric_j)
    assert best.complete_electric_j < trace.store.practical_j
    assert widest.complete_electric_j > helmholtz
    assert widest.complete_electric_j < trace.store.practical_j
    assert widest.complete_electric_j < trace.store.unrecovered_store_j
    assert widest.cap_complete_electric_j > trace.store.practical_j
    assert widest.electric_j < helmholtz
    assert widest.efold_s / widest.walk_time_s > 1.0e6
    original_complete = trace.operating_per_j * trace.checked.fusion_j / trace.checked.burn_fraction
    assert original_complete < trace.store.unrecovered_store_j
    assert trace.cap_per_j > trace.operating_per_j
    assert trace.cap_per_j < 1.0


def test_shell_rise_is_a_fraction_of_a_kelvin(trace):
    assert 1.0e5 < trace.shell.capacity_j_per_k < 1.0e6
    assert trace.shell.checked_rise_k < 1.0e-3
    assert trace.shell.ceiling_rise_k < 1.0
    assert trace.shell.checked_neutron_j < trace.shell.ceiling_neutron_j


def test_fsot_record_is_the_pin_and_stays_out_of_the_joules(trace):
    assert trace.witness.pin_prefix == PIN_PREFIX
    assert trace.witness.digest_prefix == PIN_PREFIX
    assert trace.witness.valve_beta < 1.0e-15
    assert not math.isclose(trace.witness.valve_beta, trace.checked.peak_beta)
    assert "K" in trace.witness.scalars
    assert "gamma_FSOT" in trace.witness.scalars
    text = render()
    assert "The FSOT scalars are recorded and stay out of the joule product." in text
    assert "Engine valve and plasma beta are different quantities." in text
    assert "The schedule burns a fraction of the ceiling." in text
    assert "The confined schedule remains a net consumer on every beta = 1 row." in text
    assert "Widest complete burn versus the Helmholtz ledger: above." in text
    assert "Widest complete burn versus the practical supply: below." in text
    assert "Widest converter-cap complete burn versus the practical supply: above." in text
    assert "Original-density complete burn versus the unrecovered store: below." in text
    reciprocal = 1.0 / trace.witness.valve_beta
    assert f"{reciprocal:.6e}" not in text
    assert f"{reciprocal:.6g}" not in text
    assert "3.328" not in text
    assert "3.703" not in text
    assert "warp_shell_energy_j" not in text
    assert text == render()


def test_report_writes_under_runs_and_the_readme_orders_the_heading():
    output = default_output_path()
    assert output.parent.parent == project_root() / "runs"
    assert main() == 0
    written = output.read_text(encoding="utf-8")
    assert written.startswith("# Failure trace\n")
    assert "py -3.14 -m solcascade.production.trace" not in written
    readme = (project_root() / "README.md").read_text(encoding="utf-8")
    refinement = readme.index("## Production refinement")
    failure = readme.index("## Failure trace")
    cost = readme.index("## Dominant cost")
    surface = readme.index("## Surface field")
    nature = readme.index("## Nature")
    throat = readme.index("## Throat")
    plants = readme.index("## Working plants")
    observables = readme.index("## Fission observables")
    tests = readme.index("## Tests")
    assert refinement < failure < cost < surface < nature < throat < plants < observables < tests
    assert "py -3.14 -m solcascade.production.trace" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "trace.py"
    ).read_text(encoding="utf-8")
    assert "warp_shell_energy_j" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9"):
        assert needle not in source
