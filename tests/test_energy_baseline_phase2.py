"""Phase-two mapping math. The phase-one tests stay untouched."""

import math
from pathlib import Path

import yaml

from solcascade.constants import warp_shell_energy_j
from solcascade.energy_baseline.phase2.gaps import render_gaps, sort_gaps
from solcascade.energy_baseline.phase2.map_fsot import (
    CAPTURED_DOMAINS,
    EXTRA_OUTPUTS,
    READINGS,
    capture_scalars,
    compare_record,
    default_paths,
    load_engine,
    map_directory,
    render_markdown,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "src" / "solcascade" / "energy_baseline" / "data"
MAP_PY = ROOT / "src" / "solcascade" / "energy_baseline" / "phase2" / "map_fsot.py"
TARGET = 1.0e45


def _record(**overrides):
    row = {
        "id": "fixture",
        "name": "Fixture",
        "category": "mainstream",
        "governing_equation": "coverage scales with one dimensionless efficiency",
        "governing_source": "https://example.invalid/law",
        "efficiency_ceiling": 1.0,
        "efficiency_ceiling_kind": "other",
        "efficiency_ceiling_source": "fixture cap",
        "demonstrated_efficiency": 0.05,
        "demonstrated_year": 2024,
        "demonstrated_source": "https://example.invalid/device",
        "energy_density_j": 1.0e6,
        "energy_density_basis": "J per reference unit",
        "reference_unit": "one unit",
        "power_density_w": 1.0,
        "power_density_basis": "W",
        "input_energy_j": 1.0,
        "input_energy_basis": "J",
        "chain_reaction": False,
        "self_sustaining": False,
        "continuous_external_drive": True,
        "trl": 2,
        "engineering_barrier": "fixture",
        "fsot_status": "open_gap",
        "fsot_scalar": None,
        "fsot_mapping": None,
        "fsot_note": "domain_scalar(Optics)",
        "claimed_output_w": None,
        "replication_status": "not_applicable",
        "flag": None,
        "figure_status": "estimate",
        "notes": "fixture",
        "sources": ["https://example.invalid/law"],
    }
    row.update(overrides)
    return row


def test_map_directory_runs_on_a_minimal_fixture(tmp_path: Path):
    (tmp_path / "fixture.yaml").write_text(
        yaml.safe_dump(_record()),
        encoding="utf-8",
    )
    scalars = {"Optics": 0.2}
    rows = map_directory(tmp_path, scalars, TARGET)
    text = render_markdown(rows, TARGET, scalars, pin=None)
    assert len(rows) == 1
    assert "fixture" in text
    assert "Warp-shell target" in text


def test_in_range_scalar_is_detected():
    # Wave's reading is linear in the file efficiency. The tenfold window
    # on 0.05 with a cap of 1 is [0.5, 1]. 0.8 sits inside it.
    row = compare_record(
        _record(id="wave", fsot_note="domain_scalar(Oceanography)"),
        {"Oceanography": 0.8},
        TARGET,
    )
    assert row["in_range"] is True
    coverage = 1.0e6 * 0.05 / TARGET
    assert math.isclose(row["coverage"], coverage)
    assert math.isclose(row["coverage_if_fixed"], coverage * (0.8 / 0.05))
    assert "In range" in row["gap_note"]


def test_out_of_range_scalar_is_flagged():
    row = compare_record(
        _record(id="wave", fsot_note="domain_scalar(Oceanography)"),
        {"Oceanography": 0.2},
        TARGET,
    )
    assert row["in_range"] is False
    assert row["coverage_if_fixed"] is None
    assert "outside the window" in row["gap_note"]


def test_physical_cap_blocks_a_tenfold_move():
    row = compare_record(
        _record(
            id="wave",
            efficiency_ceiling=0.3,
            fsot_note="domain_scalar(Oceanography)",
        ),
        {"Oceanography": 0.94},
        TARGET,
    )
    assert row["in_range"] is False
    assert row["coverage_if_fixed"] is None
    assert "capped" in row["gap_note"]


def test_gaps_text_is_stable_across_two_runs(tmp_path: Path):
    records = [
        _record(id="b-linear", fsot_note="domain_scalar(Optics)"),
        _record(
            id="a-linear",
            demonstrated_efficiency=0.2,
            efficiency_ceiling=0.4,
            fsot_note="domain_scalar(Optics)",
        ),
    ]
    # The fixture ids are not in READINGS, so give them known kinds by
    # comparing records that use catalog readings through written files
    # of two real shapes: one dimensionless cap and one repeated render.
    scalars = {"Optics": 0.2, "Nuclear_Physics": 0.9363875640629747}
    first = [
        compare_record(records[0], scalars, TARGET),
        compare_record(records[1], scalars, TARGET),
    ]
    # Unknown ids take the unscored path. Build a mixed set from the
    # pure functions so the sort sees upside, a zero, and a missing row.
    hit = compare_record(_record(id="wave"), {"Oceanography": 0.8}, TARGET)
    capped = compare_record(
        _record(id="thermionic", demonstrated_efficiency=0.05, efficiency_ceiling=0.3),
        {"Thermodynamics": 0.94},
        TARGET,
    )
    blank = compare_record(_record(id="unknown-row"), scalars, TARGET)
    rows = [blank, capped, hit, first[0]]
    rendered = [render_gaps(rows, TARGET) for _ in range(2)]
    assert rendered[0] == rendered[1]
    ordered = sort_gaps(rows)
    assert [row["id"] for row in ordered] == [row["id"] for row in sort_gaps(list(reversed(rows)))]
    path_a = tmp_path / "gaps-a.md"
    path_b = tmp_path / "gaps-b.md"
    path_a.write_text(rendered[0], encoding="utf-8")
    path_b.write_text(render_gaps(list(reversed(rows)), TARGET), encoding="utf-8")
    assert path_a.read_bytes() == path_b.read_bytes()


def test_catalog_readings_cover_every_file_and_do_not_hardcode_the_shell():
    stems = {path.stem for path in DATA.glob("*.yaml")}
    assert stems == set(READINGS)
    assert len(stems) == 38
    text = MAP_PY.read_text(encoding="utf-8")
    assert "3.328" not in text
    assert "3.703" not in text
    assert "warp_shell_energy_j" in text
    captured = set(CAPTURED_DOMAINS) | set(EXTRA_OUTPUTS)
    for reading in READINGS.values():
        if reading.candidate is not None:
            assert reading.candidate in captured
    _catalog, mapping, gaps = default_paths()
    assert _catalog == DATA
    assert mapping == ROOT / "runs" / "energy_baseline" / "phase2_mapping.md"
    assert gaps == ROOT / "runs" / "energy_baseline" / "gaps.md"
    rows = map_directory(DATA, {name: 0.01 for name in captured}, warp_shell_energy_j())
    assert len(rows) == 38
    assert all(row["in_range"] is False for row in rows)


def test_live_engine_pin_and_named_scalars():
    module, digest = load_engine()
    assert digest.startswith("aeb2adad6e80f487")
    scalars = capture_scalars(module)
    assert math.isclose(scalars["Nuclear_Physics"], 0.9363875640629747, rel_tol=1e-12)
    assert math.isclose(scalars["Electromagnetism"], 0.9557285700955828, rel_tol=1e-12)
    assert math.isclose(scalars["K"], 0.4201087636498879, rel_tol=1e-12)
    assert math.isclose(scalars["phi^-4"], 0.145898, rel_tol=1e-4)
    condensed = float(module.domain_scalar("Condensed_Matter"))
    assert scalars["gamma_FSOT"] == condensed / scalars["K"]
    assert math.isclose(scalars["gamma_FSOT"], 2.241, rel_tol=1e-3)
