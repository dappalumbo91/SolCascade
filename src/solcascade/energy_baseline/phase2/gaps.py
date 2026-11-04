"""Collect phase-one and phase-two open gaps.

The sort is deterministic. Computable upside comes first: the coverage gained
by pushing the dominant dimensionless parameter to its physical cap. That
number is catalog arithmetic. Dimensionful gaps follow, ordered by the
coverage already scored, because this pass does not invent the factor a
missing unit conversion would provide. Rows with no scored lever sort last,
by mechanism id.
"""

from __future__ import annotations

import sys
from pathlib import Path

from solcascade.constants import warp_shell_energy_j
from solcascade.energy_baseline.phase2.map_fsot import (
    capture_scalars,
    default_paths,
    load_engine,
    map_directory,
)


def sort_gaps(rows: list[dict]) -> list[dict]:
    def key(row: dict):
        upside = row.get("upside")
        mechanism_id = "" if row.get("id") is None else str(row["id"])
        if upside is not None:
            return (0, -float(upside), mechanism_id)
        if row.get("kind") == "dimensionful":
            coverage = row.get("coverage")
            missing = coverage is None
            return (1, missing, -(coverage if coverage is not None else 0.0), mechanism_id)
        return (2, mechanism_id)

    return sorted(rows, key=key)


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value).replace("|", "/").replace("\n", " ")


def render_gaps(rows: list[dict], target_j: float) -> str:
    ordered = sort_gaps(rows)
    lines = [
        "# Energy baseline gaps",
        "",
        "Phase one left `fsot_status: open_gap` on every mechanism file.",
        "Those files are not edited here. A phase-two row is `in_range` only",
        "when a captured scalar sits in the tenfold window under the physical cap.",
        "Every other phase-two row stays an open gap.",
        "",
        "Sort: computable upside first (coverage gained by taking the dominant",
        "dimensionless parameter to its cap), then dimensionful gaps by the",
        "coverage already scored, then unscored rows by id. Equal keys break on id.",
        "Dimensionful closures are not given a made-up factor.",
        "",
        f"Warp-shell target: `{target_j:.6g}` J, from `warp_shell_energy_j()`.",
        "",
        "| Mechanism id | Phase 1 | Phase 2 | Dominant parameter | Upside | Gap note |",
        "| --- | --- | --- | --- | ---: | --- |",
    ]
    for row in ordered:
        phase2 = "in_range" if row.get("in_range") else "open_gap"
        lines.append(
            "| {id} | {phase1} | {phase2} | {parameter} | {upside} | {note} |".format(
                id=_cell(row.get("id")),
                phase1=_cell(row.get("phase1_status") or "open_gap"),
                phase2=phase2,
                parameter=_cell(row.get("parameter")),
                upside=_cell(row.get("upside")),
                note=_cell(row.get("gap_note")),
            )
        )
    lines.append("")
    return "\n".join(lines)


def compose_gaps(rows: list[dict], target_j: float) -> str:
    """Phase-two table plus the phase-three section. Both renders match."""
    from solcascade.energy_baseline.phase3.derive import phase3_gap_section

    body = render_gaps(rows, target_j).rstrip("\n")
    extra = phase3_gap_section().rstrip("\n")
    return body + "\n\n" + extra + "\n"


def main(argv: list[str] | None = None) -> int:
    del argv
    catalog, _mapping, output = default_paths()
    module, _digest = load_engine()
    scalars = capture_scalars(module)
    target = float(warp_shell_energy_j())
    text = compose_gaps(map_directory(catalog, scalars, target), target)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text, encoding="utf-8")
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
