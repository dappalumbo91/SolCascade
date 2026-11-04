"""Phase-three derivations. Phase-one and phase-two tests stay untouched."""

import hashlib
import json
import math
from pathlib import Path

import yaml

from solcascade.constants import (
    BOSCH_HALE_DHE3_T_MAX_KEV,
    BOSCH_HALE_DHE3_T_MIN_KEV,
    NEVINS_SWAIN_PB11_T_MAX_KEV,
    NEVINS_SWAIN_PB11_T_MIN_KEV,
    warp_shell_energy_j,
)
from solcascade.energy_baseline.phase2.gaps import compose_gaps, render_gaps
from solcascade.energy_baseline.phase3.derive import (
    CATALOG_AIR_DENSITY_KG_M3,
    CATALOG_GRAVITY_M_S2,
    CATALOG_WATER_DENSITY_KG_M3,
    CASIMIR_STATED_AREA_M2,
    CASIMIR_STATED_GAP_M,
    FUSION_IDS,
    NEXT_TARGET,
    RESULT_WHEN_NONE,
    annotate,
    apply_coverage_factor,
    blackbody_energy_fraction_above,
    blackbody_energy_tail,
    casimir_energy_magnitude_j,
    catalog_wind_speed_m_s,
    coverage_changed,
    hydraulic_head_m,
    load_overrides,
    moves_by_order,
    phase3_gap_section,
    photovoltaic_coverage_factor,
    solar_blackbody_kt_ev,
    thermal_ion_speed_m_s,
    tidal_range_m,
    walk_clock,
)
from solcascade.plasma.reactivity import dhe3_reactivity_m3_s, pb11_reactivity_m3_s
from solcascade.energy_baseline.phase3.rank_derived import (
    build_gaps,
    build_ranking,
    default_paths,
    render_derivations,
    render_ranking,
)
from solcascade.energy_baseline.rank import rank_directory, score_record, sort_rows
from solcascade.plasma.metrics import bohm_time_s, gyro_walk_time_s

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "src" / "solcascade" / "energy_baseline" / "data"
PHASE3 = ROOT / "src" / "solcascade" / "energy_baseline" / "phase3"
OVERRIDES = PHASE3 / "overrides.yaml"
REPORT = ROOT / "runs" / "fsot-apply" / "report.json"


def _minimal(**overrides):
    row = {
        "id": "sample",
        "name": "Sample",
        "category": "mainstream",
        "governing_equation": "stored joules times a stated efficiency",
        "governing_source": "https://example.invalid/law",
        "efficiency_ceiling": 1.0,
        "efficiency_ceiling_kind": "other",
        "efficiency_ceiling_source": "fixture cap",
        "demonstrated_efficiency": 0.5,
        "demonstrated_year": 2020,
        "demonstrated_source": "https://example.invalid/device",
        "energy_density_j": 100.0,
        "energy_density_basis": "J per reference unit",
        "reference_unit": "one unit",
        "power_density_w": 1.0,
        "power_density_basis": "W",
        "input_energy_j": 40.0,
        "input_energy_basis": "J",
        "chain_reaction": False,
        "self_sustaining": False,
        "continuous_external_drive": True,
        "trl": 3,
        "engineering_barrier": "fixture",
        "fsot_status": "open_gap",
        "fsot_scalar": None,
        "fsot_mapping": None,
        "fsot_note": "open gap",
        "claimed_output_w": None,
        "replication_status": "not_applicable",
        "flag": None,
        "figure_status": "estimate",
        "notes": "fixture",
        "sources": ["https://example.invalid/law"],
    }
    row.update(overrides)
    return row


def _catalog_digest() -> str:
    digest = hashlib.sha256()
    for path in sorted(DATA.glob("*.yaml")):
        digest.update(path.read_bytes())
    return digest.hexdigest()


def test_walk_clock_matches_the_recorded_peak_and_bohm_ratio():
    loaded = yaml.safe_load(OVERRIDES.read_text(encoding="utf-8"))
    clock = walk_clock(loaded["radius_m"], loaded["temperature_keV"], loaded["field_T"])
    direct = gyro_walk_time_s(loaded["radius_m"], loaded["temperature_keV"], loaded["field_T"])
    bohm = bohm_time_s(loaded["radius_m"], loaded["temperature_keV"], loaded["field_T"])
    assert math.isclose(clock.walk_time_s, direct, rel_tol=0.0, abs_tol=0.0)
    assert math.isclose(clock.walk_over_bohm, math.pi / 4.0, rel_tol=1e-12)
    assert math.isclose(clock.walk_time_s / bohm, math.pi / 4.0, rel_tol=1e-12)
    assert clock.gamma_enters_the_time is False
    report = json.loads(REPORT.read_text(encoding="utf-8"))["confinement"]
    assert loaded["radius_m"] == report["radius_m"]
    assert loaded["temperature_keV"] == report["temperature_keV"]
    assert loaded["field_T"] == report["field_T"]
    assert math.isclose(clock.walk_time_s, report["walk_time_s"], rel_tol=1e-9)
    assert report["gamma_enters_the_time"] is False
    assert report["installed_time"] == "gyro_period_random_walk"


def test_confined_yield_moves_scored_dt_rows_by_orders_of_magnitude():
    before = _catalog_digest()
    target = warp_shell_energy_j()
    pulse = yaml.safe_load((ROOT / "configs" / "minimal_pipeline.yaml").read_text(encoding="utf-8"))
    mapping = yaml.safe_load(OVERRIDES.read_text(encoding="utf-8"))
    assert mapping["ion_density_m3"] == pulse["plasma"]["ion_density_m3"]
    assert mapping["ion_temperature_keV"] == pulse["plasma"]["ion_temperature_keV"]
    assert mapping["radius0_m"] == pulse["compression"]["plasma_radius_m"][0]
    assert mapping["radius_m"] == min(pulse["compression"]["plasma_radius_m"])
    assert mapping["deuterium_fraction"] == pulse["plasma"]["deuterium_fraction"]
    loaded = load_overrides(OVERRIDES)
    clock = loaded.clock
    confined = loaded.confined_yield
    derivations = loaded.derivations
    assert clock.gamma_enters_the_time is False
    fusion_ids = [item.mechanism_id for item in derivations if item.mechanism_id in FUSION_IDS]
    assert fusion_ids == list(FUSION_IDS)
    assert math.isclose(
        confined.energy_j,
        confined.power_w * clock.walk_time_s,
        rel_tol=0.0,
        abs_tol=0.0,
    )
    base = {row["id"]: row for row in rank_directory(DATA, target)}
    rows, _clock, ranking, derivations_text = build_ranking(DATA, OVERRIDES, target)
    assert _catalog_digest() == before
    for path in DATA.glob("*.yaml"):
        record = yaml.safe_load(path.read_text(encoding="utf-8"))
        assert record["fsot_status"] == "open_gap"
        assert record["fsot_scalar"] is None
        assert record["fsot_mapping"] is None
    assert len(rows) == len(base) == 38
    scored = ("fusion-dt-magnetic", "fusion-dt-inertial")
    for row in rows:
        prior = base[row["id"]]["output_over_target"]
        assert row["old_coverage"] == prior or (row["old_coverage"] is None and prior is None)
        if row["id"] in scored:
            delivered = base[row["id"]]["delivered_j"]
            factor = confined.energy_j / delivered
            assert math.isclose(row["coverage_factor"], factor, rel_tol=1e-12)
            assert math.isclose(row["new_coverage"], prior * factor, rel_tol=1e-12)
            assert moves_by_order(row["old_coverage"], row["new_coverage"]) is True
            assert abs(math.log10(factor)) >= 1.0
            assert f"{(row['new_coverage'] / row['old_coverage']):.12g}" in ranking
            assert row["derived_formula"] == "confined_dt_yield_j"
        elif row["id"] == "fusion-dt-pulsed-power":
            assert prior is None
            assert math.isclose(row["new_coverage"], confined.energy_j / target, rel_tol=1e-12)
            assert row["coverage_factor"] is None
            assert coverage_changed(prior, row["new_coverage"]) is True
        elif row["id"] == "fusion-aneutronic-pb11":
            assert prior is None
            assert row["coverage_factor"] is None
            assert row["new_coverage"] is not None
            assert row["derived_formula"] == "confined_pb11_yield_j"
            assert coverage_changed(prior, row["new_coverage"]) is True
        elif row["id"] == "fusion-d-he3":
            assert prior is None
            assert row["coverage_factor"] is None
            assert row["new_coverage"] is not None
            assert row["derived_formula"] == "confined_dhe3_yield_j"
            assert coverage_changed(prior, row["new_coverage"]) is True
        elif row["id"] in {"hydroelectric", "tidal", "wind", "photovoltaic"}:
            assert row["coverage_factor"] is not None
            assert moves_by_order(row["old_coverage"], row["new_coverage"]) is True
        elif row["id"] == "casimir-static":
            assert prior == 0.0
            assert row["new_coverage"] == 0.0
            assert row["coverage_factor"] is None
            assert row["derived_side_factor"] is not None
            assert coverage_changed(prior, row["new_coverage"]) is False
        else:
            assert row["coverage_factor"] is None
            assert row["new_coverage"] == prior or (row["new_coverage"] is None and prior is None)
    assert RESULT_WHEN_NONE not in ranking
    assert "confined_dt_yield_j" in ranking
    assert "moves by a factor of" in ranking
    assert NEXT_TARGET in ranking
    assert NEXT_TARGET in derivations_text
    assert "gyro_walk_time_s" in derivations_text
    assert "_fsot_scaled" in derivations_text
    assert ranking == render_ranking(rows, target, clock, confined)
    fresh = load_overrides(OVERRIDES)
    assert derivations_text == render_derivations(fresh.clock, fresh.confined_yield, fresh.derivations)


def test_pulse_units_and_dhe3_burn_move_the_named_rows():
    target = warp_shell_energy_j()
    loaded = load_overrides(OVERRIDES)
    clock = loaded.clock
    confined = loaded.confined_yield
    assert BOSCH_HALE_DHE3_T_MIN_KEV < confined.temperature_kev < BOSCH_HALE_DHE3_T_MAX_KEV
    try:
        dhe3_reactivity_m3_s(BOSCH_HALE_DHE3_T_MIN_KEV * 0.5)
    except ValueError:
        outside = True
    else:
        outside = False
    assert outside is True
    full = math.pi ** 4 / 15.0
    assert math.isclose(blackbody_energy_tail(0.0), full, rel_tol=1e-9)
    assert blackbody_energy_fraction_above(0.0) == 1.0
    assert 0.0 < blackbody_energy_fraction_above(2.0) < 1.0
    gap_ev = clock.temperature_kev * 1000.0
    assert blackbody_energy_fraction_above(gap_ev / solar_blackbody_kt_ev()) == 0.0
    records = {
        path.stem: yaml.safe_load(path.read_text(encoding="utf-8"))
        for path in DATA.glob("*.yaml")
    }
    head = hydraulic_head_m(records["hydroelectric"]["energy_density_j"])
    span = tidal_range_m(records["tidal"]["energy_density_j"])
    assert math.isclose(CATALOG_GRAVITY_M_S2 * head, records["hydroelectric"]["energy_density_j"], rel_tol=1e-12)
    assert math.isclose(head, 100.0, rel_tol=1e-12)
    basin = 0.5 * CATALOG_WATER_DENSITY_KG_M3 * CATALOG_GRAVITY_M_S2 * span ** 2
    assert math.isclose(basin, records["tidal"]["energy_density_j"], rel_tol=1e-12)
    assert math.isclose(span, 8.0, rel_tol=1e-12)
    wind_speed = catalog_wind_speed_m_s(
        records["wind"]["power_density_w"],
        records["wind"]["demonstrated_efficiency"],
    )
    kinetic = records["wind"]["power_density_w"] / records["wind"]["demonstrated_efficiency"]
    assert math.isclose(0.5 * CATALOG_AIR_DENSITY_KG_M3 * wind_speed ** 3, kinetic, rel_tol=1e-12)
    thermal = thermal_ion_speed_m_s(clock.temperature_kev)
    pv_factor = photovoltaic_coverage_factor(gap_ev, records["photovoltaic"]["demonstrated_efficiency"])
    assert pv_factor == 0.0
    exact_cavity = casimir_energy_magnitude_j(CASIMIR_STATED_AREA_M2, CASIMIR_STATED_GAP_M)
    stored_cavity = records["casimir-static"]["energy_density_j"]
    assert abs(exact_cavity - stored_cavity) / stored_cavity < 0.02
    cavity_factor = (CASIMIR_STATED_GAP_M / clock.radius_m) ** 3
    reactivity = dhe3_reactivity_m3_s(confined.temperature_kev)
    partner = (1.0 - confined.deuterium_fraction) * confined.density_m3
    burn = partner * reactivity * clock.walk_time_s
    assert 0.0 < burn < 1.0
    dhe3_joules = records["fusion-d-he3"]["energy_density_j"] * burn
    assert NEVINS_SWAIN_PB11_T_MIN_KEV < confined.temperature_kev < NEVINS_SWAIN_PB11_T_MAX_KEV
    for outside_kev in (
        NEVINS_SWAIN_PB11_T_MIN_KEV * 0.5,
        NEVINS_SWAIN_PB11_T_MAX_KEV + 30.0,
    ):
        try:
            pb11_reactivity_m3_s(outside_kev)
        except ValueError:
            rejected = True
        else:
            rejected = False
        assert rejected is True
    pb11_reactivity = pb11_reactivity_m3_s(confined.temperature_kev)
    pb11_burn = partner * pb11_reactivity * clock.walk_time_s
    assert 0.0 < pb11_burn < 1.0
    pb11_joules = records["fusion-aneutronic-pb11"]["energy_density_j"] * pb11_burn
    rows, _clock, ranking, derivations_text = build_ranking(DATA, OVERRIDES, target)
    by_id = {row["id"]: row for row in rows}
    assert math.isclose(by_id["hydroelectric"]["coverage_factor"], clock.radius_m / head, rel_tol=1e-12)
    assert math.isclose(by_id["tidal"]["coverage_factor"], (clock.radius_m / span) ** 2, rel_tol=1e-12)
    assert math.isclose(by_id["wind"]["coverage_factor"], (thermal / wind_speed) ** 3, rel_tol=1e-12)
    assert by_id["photovoltaic"]["coverage_factor"] == 0.0
    assert by_id["photovoltaic"]["new_coverage"] == 0.0
    assert by_id["hydroelectric"]["coverage_factor"] < 1.0
    assert by_id["tidal"]["coverage_factor"] < 1.0
    assert by_id["wind"]["coverage_factor"] > 1.0
    assert by_id["wind"]["new_coverage"] < 1.0e-10
    assert by_id["wind"]["new_coverage"] > by_id["fission-fast"]["new_coverage"]
    assert math.isclose(by_id["casimir-static"]["derived_side_factor"], cavity_factor, rel_tol=1e-12)
    assert by_id["casimir-static"]["delivered_j"] == 0.0
    assert math.isclose(by_id["fusion-d-he3"]["delivered_j"], dhe3_joules, rel_tol=1e-12)
    assert math.isclose(by_id["fusion-d-he3"]["new_coverage"], dhe3_joules / target, rel_tol=1e-12)
    assert math.isclose(by_id["fusion-aneutronic-pb11"]["delivered_j"], pb11_joules, rel_tol=1e-12)
    assert math.isclose(
        by_id["fusion-aneutronic-pb11"]["new_coverage"],
        pb11_joules / target,
        rel_tol=1e-12,
    )
    assert by_id["fusion-aneutronic-pb11"]["coverage_factor"] is None
    assert by_id["fusion-aneutronic-pb11"]["derived_formula"] == "confined_pb11_yield_j"
    assert f"{(clock.radius_m / head):.12g}" in ranking
    assert f"{((clock.radius_m / span) ** 2):.12g}" in ranking
    assert f"{((thermal / wind_speed) ** 3):.12g}" in ranking
    assert f"{cavity_factor:.12g}" in ranking
    assert "factor of 0 via" in ranking
    assert "falls to 0" in ranking
    assert "moves down" in ranking
    assert "moves up and stays" in ranking
    assert "under the shell" in ranking
    assert "had a blank coverage" in ranking
    assert "Nevins and Swain" in NEXT_TARGET
    assert "Nevins and Swain" in derivations_text
    assert "demonstrated efficiency is 0" in ranking


def test_factor_of_ten_moves_coverage_and_is_reported():
    target = 1.0e45
    base = score_record(_minimal(), target)
    moved = apply_coverage_factor(
        base,
        10.0,
        quantity="confinement time",
        source="fixture-source",
        formula_id="fixture-law",
        value=1.0e-6,
        unit="s",
    )
    assert math.isclose(moved["new_coverage"], base["output_over_target"] * 10.0)
    assert math.isclose(moved["delivered_j"], base["delivered_j"] * 10.0)
    assert moves_by_order(moved["old_coverage"], moved["new_coverage"]) is True
    loaded = load_overrides(OVERRIDES)
    text = render_ranking(sort_rows([moved]), target, loaded.clock, loaded.confined_yield)
    assert "factor of 10" in text
    assert "fixture-law" in text
    assert "fixture-source" in text
    assert RESULT_WHEN_NONE not in text
    blank = score_record(_minimal(id="blank", energy_density_j=None), target)
    untouched = apply_coverage_factor(
        blank,
        10.0,
        quantity="confinement time",
        source="fixture-source",
        formula_id="fixture-law",
    )
    assert untouched["old_coverage"] is None
    assert untouched["new_coverage"] is None
    assert moves_by_order(None, None) is False


def test_phase3_sources_omit_the_shell_literals_and_paths_resolve():
    for path in PHASE3.iterdir():
        if path.suffix not in {".py", ".yaml"}:
            continue
        text = path.read_text(encoding="utf-8")
        assert "3.328" not in text
        assert "3.703" not in text
    catalog, ranking, derivations, gaps = default_paths()
    assert catalog == DATA
    assert ranking == ROOT / "runs" / "energy_baseline" / "phase3_ranking.md"
    assert derivations == ROOT / "runs" / "energy_baseline" / "phase3_derivations.md"
    assert gaps == ROOT / "runs" / "energy_baseline" / "gaps.md"
    assert not (DATA / "overrides.yaml").exists()
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    phase2 = readme.index("## Energy Baseline Phase Two")
    phase3 = readme.index("## Phase Three — Derived Values")
    tests = readme.index("## Tests")
    assert phase2 < phase3 < tests
    assert "py -3.14 -m solcascade.energy_baseline.phase3.rank_derived" in readme


def test_gap_section_is_stable_and_keeps_the_phase_two_table():
    loaded = load_overrides(OVERRIDES)
    first = phase3_gap_section(loaded.clock)
    second = phase3_gap_section(loaded.clock)
    assert first == second
    assert "gyro_walk_time_s" in first
    assert NEXT_TARGET in first
    target = warp_shell_energy_j()
    gaps_a = build_gaps(DATA, target)
    gaps_b = build_gaps(DATA, target)
    assert gaps_a == gaps_b
    assert gaps_a.startswith("# Energy baseline gaps")
    assert "## Phase three" in gaps_a
    assert render_gaps([], target).strip() in compose_gaps([], target)
    phase2_only = render_gaps([], target)
    assert "## Phase three" not in phase2_only
    phase1 = rank_directory(DATA, target)
    annotated = annotate(phase1, load_overrides(OVERRIDES).derivations, target)
    phase3_ids = [row["id"] for row in sort_rows(annotated)]
    assert phase3_ids == [row["id"] for row in sort_rows(list(reversed(annotated)))]
    assert [row["id"] for row in phase1].index("fusion-dt-magnetic") < phase3_ids.index(
        "fusion-dt-magnetic"
    )
