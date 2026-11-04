"""Schema and ranking math for the energy-baseline catalog."""

import math
from pathlib import Path

import pytest
import yaml

from solcascade.constants import warp_shell_energy_j
from solcascade.energy_baseline.rank import (
    REQUIRED_FIELDS,
    default_paths,
    rank_directory,
    score_record,
    sort_rows,
    validate_record,
    warp_target_j,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "src" / "solcascade" / "energy_baseline" / "data"
RANK_PY = ROOT / "src" / "solcascade" / "energy_baseline" / "rank.py"
MECHANISM_FILES = sorted(DATA.glob("*.yaml"))


def _minimal(**overrides):
    row = {
        "id": "sample",
        "name": "Sample",
        "category": "mainstream",
        "governing_equation": "stored joules times a stated efficiency",
        "governing_source": "https://example.invalid/law",
        "efficiency_ceiling": 1.0,
        "efficiency_ceiling_kind": "other",
        "efficiency_ceiling_source": "ideal bound used only by the fixture",
        "demonstrated_efficiency": 0.5,
        "demonstrated_year": 2020,
        "demonstrated_source": "https://example.invalid/device",
        "energy_density_j": 100.0,
        "energy_density_basis": "J per reference unit",
        "reference_unit": "one unit",
        "power_density_w": 1.0,
        "power_density_basis": "W for the fixture",
        "input_energy_j": 40.0,
        "input_energy_basis": "J to start the fixture",
        "chain_reaction": False,
        "self_sustaining": False,
        "continuous_external_drive": True,
        "trl": 3,
        "engineering_barrier": "none in this fixture",
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
    missing = [field for field in REQUIRED_FIELDS if field not in row]
    assert missing == []
    return row


def test_catalog_covers_at_least_one_file_per_ten_mechanisms():
    assert len(MECHANISM_FILES) >= 30


@pytest.mark.parametrize("path", MECHANISM_FILES, ids=lambda path: path.stem)
def test_mechanism_file_matches_schema(path: Path):
    record = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert validate_record(record) == []
    assert record["fsot_status"] == "open_gap"
    assert record["fsot_scalar"] is None
    assert record["fsot_mapping"] is None


def test_delivered_over_target_uses_the_constants_function():
    target = warp_shell_energy_j()
    assert warp_target_j() == target
    row = score_record(
        _minimal(energy_density_j=200.0, demonstrated_efficiency=0.25),
        target,
    )
    assert row["schema_errors"] == []
    assert row["delivered_j"] == 50.0
    assert math.isclose(row["output_over_target"], 50.0 / target)
    assert math.isclose(row["units_to_target"], target / 50.0)
    assert math.isclose(row["output_over_input"], 50.0 / 40.0)


def test_equal_coverage_sorts_by_id():
    target = 1.0e6
    rows = [
        score_record(
            _minimal(id="b", energy_density_j=10.0, demonstrated_efficiency=1.0),
            target,
        ),
        score_record(
            _minimal(id="a", energy_density_j=10.0, demonstrated_efficiency=1.0),
            target,
        ),
        score_record(
            _minimal(id="c", energy_density_j=20.0, demonstrated_efficiency=1.0),
            target,
        ),
    ]
    assert [row["id"] for row in sort_rows(rows)] == ["c", "a", "b"]


def test_missing_energy_density_sorts_last_with_blank_coverage():
    target = warp_shell_energy_j()
    present = score_record(_minimal(id="present"), target)
    blank = score_record(_minimal(id="blank", energy_density_j=None), target)
    assert blank["delivered_j"] is None
    assert blank["output_over_target"] is None
    assert [row["id"] for row in sort_rows([blank, present])] == ["present", "blank"]


def test_missing_required_key_does_not_raise():
    row = score_record({"id": "broken"}, warp_shell_energy_j())
    assert any(error.startswith("missing ") for error in row["schema_errors"])
    assert row["output_over_target"] is None


def test_zero_or_negative_input_is_not_divided(tmp_path: Path):
    target = warp_shell_energy_j()
    zero = score_record(_minimal(input_energy_j=0.0), target)
    negative = score_record(_minimal(id="neg", input_energy_j=-5.0), target)
    assert zero["delivered_j"] == 50.0
    assert zero["output_over_input"] is None
    assert zero["input_ratio_note"]
    assert negative["output_over_input"] is None
    (tmp_path / "zero-in.yaml").write_text(
        yaml.safe_dump(_minimal(id="zero-in", input_energy_j=0.0)),
        encoding="utf-8",
    )
    (tmp_path / "bad.yaml").write_text("id: only\n", encoding="utf-8")
    rows = {row["id"]: row for row in rank_directory(tmp_path, target)}
    assert rows["zero-in"]["output_over_input"] is None
    assert rows["only"]["schema_errors"]
    assert rows["only"]["output_over_target"] is None


def test_rank_source_reads_the_shell_function():
    text = RANK_PY.read_text(encoding="utf-8")
    assert "3.328" not in text
    assert "3.703" not in text
    assert "warp_shell_energy_j" in text


def test_catalog_rank_is_sorted_and_a_flag_does_not_reorder():
    rows = rank_directory(DATA, warp_target_j())
    assert rows
    assert all(row["schema_errors"] == [] for row in rows)
    coverages = [row["output_over_target"] for row in rows]
    numeric = [value for value in coverages if value is not None]
    assert numeric == sorted(numeric, reverse=True)
    if None in coverages:
        first_blank = coverages.index(None)
        assert coverages[first_blank:] == [None] * (len(coverages) - first_blank)
    flagged = score_record(_minimal(id="m", flag="fringe note"), 1.0)
    plain = score_record(_minimal(id="m-plain", flag=None), 1.0)
    assert flagged["output_over_target"] == plain["output_over_target"]
    data_dir, output = default_paths()
    assert data_dir == DATA
    assert output == ROOT / "runs" / "energy_baseline" / "ranking.md"
