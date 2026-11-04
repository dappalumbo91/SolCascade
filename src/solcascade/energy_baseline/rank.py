"""Rank catalogued energy mechanisms against the warp-shell joule target.

The target is ``warp_shell_energy_j()`` in ``constants.py``. This module
does not retype that number and does not derive an FSOT law.

Delivered joules for one reference unit are

    energy_density_j * demonstrated_efficiency

when both values are present. The rank is that product divided by the
warp-shell target, largest first. A missing number stays blank and sorts
after every numeric row. Equal scores sort by mechanism id.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Mapping

import yaml

from solcascade.constants import warp_shell_energy_j

REQUIRED_FIELDS = (
    "id",
    "name",
    "category",
    "governing_equation",
    "governing_source",
    "efficiency_ceiling",
    "efficiency_ceiling_kind",
    "efficiency_ceiling_source",
    "demonstrated_efficiency",
    "demonstrated_year",
    "demonstrated_source",
    "energy_density_j",
    "energy_density_basis",
    "reference_unit",
    "power_density_w",
    "power_density_basis",
    "input_energy_j",
    "input_energy_basis",
    "chain_reaction",
    "self_sustaining",
    "continuous_external_drive",
    "trl",
    "engineering_barrier",
    "fsot_status",
    "fsot_scalar",
    "fsot_mapping",
    "fsot_note",
    "claimed_output_w",
    "replication_status",
    "flag",
    "figure_status",
    "notes",
    "sources",
)

CATEGORIES = frozenset({"mainstream", "cutting_edge", "fringe"})
CEILING_KINDS = frozenset({"carnot", "thermodynamic", "other"})
REPLICATION = frozenset(
    {"not_applicable", "replicated", "unreplicated", "disputed", "single_experiment"}
)
FSOT_STATUS = frozenset({"open_gap", "mapped"})
FIGURE_STATUS = frozenset({"measured", "design", "estimate", "mixed"})
NULLABLE_NUMBERS = frozenset(
    {
        "efficiency_ceiling",
        "demonstrated_efficiency",
        "demonstrated_year",
        "energy_density_j",
        "power_density_w",
        "input_energy_j",
        "claimed_output_w",
    }
)


def warp_target_j() -> float:
    """Joule requirement read from constants. Not a local literal."""
    return float(warp_shell_energy_j())


def _number(value) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def validate_record(record: Mapping) -> list[str]:
    """Return schema errors. An empty list means the file matches the catalog."""
    errors: list[str] = []
    for field in REQUIRED_FIELDS:
        if field not in record:
            errors.append(f"missing {field}")
    if errors:
        return errors
    if record["category"] not in CATEGORIES:
        errors.append("category is not mainstream, cutting_edge, or fringe")
    if record["efficiency_ceiling_kind"] not in CEILING_KINDS:
        errors.append("efficiency_ceiling_kind is not carnot, thermodynamic, or other")
    if record["replication_status"] not in REPLICATION:
        errors.append("replication_status is not a known status")
    if record["fsot_status"] not in FSOT_STATUS:
        errors.append("fsot_status must be open_gap or mapped")
    if record["figure_status"] not in FIGURE_STATUS:
        errors.append("figure_status is not measured, design, estimate, or mixed")
    if not isinstance(record["sources"], list) or not record["sources"]:
        errors.append("sources must be a non-empty list")
    for field in ("chain_reaction", "self_sustaining", "continuous_external_drive"):
        if not isinstance(record[field], bool):
            errors.append(f"{field} must be a boolean")
    trl = record["trl"]
    if isinstance(trl, bool) or not isinstance(trl, int) or not 1 <= trl <= 9:
        errors.append("trl must be an integer from 1 to 9")
    for field in NULLABLE_NUMBERS:
        value = record[field]
        if value is None:
            continue
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            errors.append(f"{field} must be a number or null")
    ceiling = _number(record["efficiency_ceiling"])
    demonstrated = _number(record["demonstrated_efficiency"])
    if ceiling is not None and not 0.0 <= ceiling <= 1.0:
        errors.append("efficiency_ceiling must be in [0, 1] when set")
    if demonstrated is not None and not 0.0 <= demonstrated <= 1.0:
        errors.append("demonstrated_efficiency must be in [0, 1] when set")
    if record["fsot_status"] == "mapped":
        if not record["fsot_scalar"] or not record["fsot_mapping"]:
            errors.append("a mapped FSOT row needs fsot_scalar and fsot_mapping")
    if record["fsot_status"] == "open_gap":
        if record["fsot_scalar"] is not None or record["fsot_mapping"] is not None:
            errors.append("an open gap has no scalar and no mapping")
    return errors


def score_record(record: Mapping, target_j: float) -> dict:
    """Score one mechanism. Missing numbers do not raise."""
    record = dict(record or {})
    schema_errors = validate_record(record)
    gross = _number(record.get("energy_density_j"))
    efficiency = _number(record.get("demonstrated_efficiency"))
    incoming = _number(record.get("input_energy_j"))
    delivered = None
    if gross is not None and efficiency is not None:
        delivered = gross * efficiency
    output_over_input = None
    input_ratio_note = None
    if delivered is not None and incoming is not None:
        if incoming > 0.0:
            output_over_input = delivered / incoming
        else:
            input_ratio_note = "input energy is zero, so output/input is undefined"
    if incoming is not None and incoming < 0.0:
        input_ratio_note = "input energy is negative and was not divided"
        output_over_input = None
    output_over_target = None
    input_over_target = None
    units_to_target = None
    if target_j > 0.0:
        if delivered is not None:
            output_over_target = delivered / target_j
            if delivered > 0.0:
                units_to_target = target_j / delivered
        if incoming is not None:
            input_over_target = incoming / target_j
    return {
        "id": record.get("id"),
        "name": record.get("name"),
        "category": record.get("category"),
        "reference_unit": record.get("reference_unit"),
        "trl": record.get("trl"),
        "energy_density_j": gross,
        "demonstrated_efficiency": efficiency,
        "delivered_j": delivered,
        "input_energy_j": incoming,
        "output_over_input": output_over_input,
        "input_ratio_note": input_ratio_note,
        "output_over_target": output_over_target,
        "input_over_target": input_over_target,
        "units_to_target": units_to_target,
        "chain_reaction": record.get("chain_reaction"),
        "self_sustaining": record.get("self_sustaining"),
        "continuous_external_drive": record.get("continuous_external_drive"),
        "fsot_status": record.get("fsot_status"),
        "replication_status": record.get("replication_status"),
        "claimed_output_w": _number(record.get("claimed_output_w")),
        "flag": record.get("flag"),
        "schema_errors": schema_errors,
    }


def sort_rows(rows: list[dict]) -> list[dict]:
    """Largest target coverage first. Blanks last. Ties break on id."""
    return sorted(
        rows,
        key=lambda row: (
            row["output_over_target"] is None,
            -(row["output_over_target"] if row["output_over_target"] is not None else 0.0),
            "" if row.get("id") is None else str(row["id"]),
        ),
    )


def load_records(directory: Path) -> list[dict]:
    """Load every YAML file. A broken file becomes an error row, not an exception."""
    records = []
    for path in sorted(directory.glob("*.yaml")):
        try:
            loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
        except yaml.YAMLError as exc:
            records.append({"id": path.stem, "name": path.name, "_load_error": str(exc)})
            continue
        if not isinstance(loaded, dict):
            records.append({"id": path.stem, "name": path.name, "_load_error": "file is not a mapping"})
            continue
        records.append(loaded)
    return records


def rank_directory(directory: Path, target_j: float) -> list[dict]:
    rows = []
    for record in load_records(directory):
        if record.get("_load_error"):
            row = score_record({}, target_j)
            row["id"] = record.get("id")
            row["name"] = record.get("name")
            row["schema_errors"] = [record["_load_error"]]
            rows.append(row)
            continue
        rows.append(score_record(record, target_j))
    return sort_rows(rows)


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value).replace("|", "/").replace("\n", " ")


def render_markdown(rows: list[dict], target_j: float) -> str:
    lines = [
        "# Energy baseline ranking",
        "",
        "Phase one catalog. Each row is one mechanism file.",
        "Delivered joules are `energy_density_j * demonstrated_efficiency`",
        "for that file's reference unit. Coverage is delivered joules",
        "divided by the warp-shell target from `warp_shell_energy_j()`.",
        "This table does not derive an FSOT law.",
        "",
        f"Warp-shell target: `{target_j:.6g}` J.",
        "",
        "| Rank | Id | Reference unit | Delivered J | Input J | Output/input | Output/target | Units to target | TRL | FSOT | Flag |",
        "| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for index, row in enumerate(rows, start=1):
        lines.append(
            "| {rank} | {id} | {unit} | {delivered} | {incoming} | {ratio} | {coverage} | {units} | {trl} | {fsot} | {flag} |".format(
                rank=index,
                id=_cell(row.get("id")),
                unit=_cell(row.get("reference_unit")),
                delivered=_cell(row.get("delivered_j")),
                incoming=_cell(row.get("input_energy_j")),
                ratio=_cell(row.get("output_over_input")),
                coverage=_cell(row.get("output_over_target")),
                units=_cell(row.get("units_to_target")),
                trl=_cell(row.get("trl")),
                fsot=_cell(row.get("fsot_status")),
                flag=_cell(row.get("flag")),
            )
        )
    problems = [row for row in rows if row.get("schema_errors")]
    if problems:
        lines.extend(["", "## Schema notes", ""])
        for row in problems:
            lines.append(f"- `{_cell(row.get('id'))}`: {'; '.join(row['schema_errors'])}")
    lines.append("")
    return "\n".join(lines)


def default_paths() -> tuple[Path, Path]:
    package = Path(__file__).resolve().parent
    root = package.parents[2]
    return package / "data", root / "runs" / "energy_baseline" / "ranking.md"


def main(argv: list[str] | None = None) -> int:
    del argv
    data_dir, output = default_paths()
    target = warp_target_j()
    rows = rank_directory(data_dir, target)
    text = render_markdown(rows, target)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text, encoding="utf-8")
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
