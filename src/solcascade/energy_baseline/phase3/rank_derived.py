"""Re-rank the energy catalog with phase-three derived overrides.

Phase-one ``rank.py`` stays on the catalog values. This module loads
``overrides.yaml``, computes the walk time, and writes the phase-three
table, the derivation note, and a gaps file that keeps the phase-two rows.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

from solcascade.constants import warp_shell_energy_j
from solcascade.energy_baseline.phase2.gaps import compose_gaps
from solcascade.energy_baseline.phase2.map_fsot import (
    capture_scalars,
    load_engine,
    map_directory,
)
from solcascade.energy_baseline.phase3.derive import (
    DHE3_FORMULA_ID,
    DHE3_LAW,
    FORMULA_ID,
    PB11_FORMULA_ID,
    PB11_LAW,
    GAP_FORMULA_ID,
    GAP_LAW,
    HEAD_FORMULA_ID,
    HEAD_LAW,
    INPUTS_SOURCE,
    LAW,
    NEXT_TARGET,
    PV_FORMULA_ID,
    PV_LAW,
    RESULT_WHEN_NONE,
    SOURCE,
    TIDAL_FORMULA_ID,
    TIDAL_LAW,
    WIND_FORMULA_ID,
    WIND_LAW,
    YIELD_FORMULA_ID,
    YIELD_LAW,
    YIELD_SOURCE,
    ConfinedYield,
    Derivation,
    WalkClock,
    _dimensionful_status,
    annotate,
    coverage_changed,
    load_overrides,
    moves_by_order,
    order_ratio,
    unmapped_parameters,
)

_LAW_BY_FORMULA = {
    FORMULA_ID: LAW,
    YIELD_FORMULA_ID: YIELD_LAW,
    DHE3_FORMULA_ID: DHE3_LAW,
    PB11_FORMULA_ID: PB11_LAW,
    HEAD_FORMULA_ID: HEAD_LAW,
    TIDAL_FORMULA_ID: TIDAL_LAW,
    WIND_FORMULA_ID: WIND_LAW,
    PV_FORMULA_ID: PV_LAW,
    GAP_FORMULA_ID: GAP_LAW,
}
from solcascade.energy_baseline.rank import rank_directory, sort_rows


def default_paths() -> tuple[Path, Path, Path, Path]:
    package = Path(__file__).resolve().parent
    catalog = package.parent / "data"
    # phase3 -> energy_baseline -> solcascade -> src -> project root
    root = package.parents[3]
    output = root / "runs" / "energy_baseline"
    return (
        catalog,
        output / "phase3_ranking.md",
        output / "phase3_derivations.md",
        output / "gaps.md",
    )


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value).replace("|", "/").replace("\n", " ")


def _direction(ratio: float, new: float) -> str:
    if ratio == 0.0:
        return "The coverage falls to 0."
    if ratio < 1.0:
        return "The coverage moves down."
    if new >= 1.0:
        return "The coverage meets or passes the shell."
    orders = -math.log10(float(new))
    return f"The coverage moves up and stays {orders:.3g} orders under the shell."


def order_sentence(rows: list[dict]) -> str:
    """Exact factor and formula for each order-of-magnitude coverage change."""
    moves = []
    created = []
    for row in rows:
        old = row.get("old_coverage")
        new = row.get("new_coverage")
        if moves_by_order(old, new):
            ratio = order_ratio(old, new)
            moves.append(
                "`{id}` moves by a factor of {factor} via `{formula}` ({source}). {direction}".format(
                    id=_cell(row.get("id")),
                    factor=f"{float(ratio):.12g}",
                    formula=_cell(row.get("derived_formula")),
                    source=_cell(row.get("derived_source")),
                    direction=_direction(float(ratio), float(new)),
                )
            )
            continue
        if old is None and new is not None:
            created.append(
                (
                    "`{id}` had a blank coverage. Derived yield `{value}` J "
                    "gives coverage `{coverage}` via `{formula}`."
                ).format(
                    id=_cell(row.get("id")),
                    value=_cell(row.get("derived_value")),
                    coverage=_cell(new),
                    formula=_cell(row.get("derived_formula")),
                )
            )
    if not moves and not created:
        return RESULT_WHEN_NONE
    parts = moves or [RESULT_WHEN_NONE]
    return " ".join(parts + created)


def _move_line(row: dict) -> str:
    ratio = order_ratio(row.get("old_coverage"), row.get("new_coverage"))
    if ratio is None:
        factor = "undefined"
    else:
        factor = f"{float(ratio):.12g}"
    return (
        "- `{id}`: old coverage `{old}`, new coverage `{new}`, factor `{factor}`, "
        "quantity `{quantity}`, formula `{formula}`, source `{source}`."
    ).format(
        id=_cell(row.get("id")),
        old=_cell(row.get("old_coverage")),
        new=_cell(row.get("new_coverage")),
        factor=factor,
        quantity=_cell(row.get("derived_quantity")),
        formula=_cell(row.get("derived_formula")),
        source=_cell(row.get("derived_source")),
    )


def render_ranking(
    rows: list[dict],
    target_j: float,
    clock: WalkClock,
    confined: ConfinedYield | None = None,
) -> str:
    changed = [row for row in rows if coverage_changed(row.get("old_coverage"), row.get("new_coverage"))]
    lines = [
        "# Energy baseline phase-three ranking",
        "",
        order_sentence(rows),
        "",
        (
            f"Installed confinement time: `{clock.walk_time_s:.6g}` s from "
            f"`{FORMULA_ID}` at radius `{clock.radius_m:.6g}` m, "
            f"temperature `{clock.temperature_kev:.6g}` keV, "
            f"field `{clock.field_t:.6g}` T."
        ),
        f"Formula: {LAW}.",
        (
            f"walk/Bohm is `{clock.walk_over_bohm:.6g}` (π/4). "
            "`gamma_enters_the_time` is false."
        ),
        f"Source: `{SOURCE}`. Inputs: `{INPUTS_SOURCE}`.",
        "",
    ]
    if confined is not None:
        lines.extend(
            [
                (
                    f"Confined DT yield: `{confined.energy_j:.6g}` J from "
                    f"`{YIELD_FORMULA_ID}`."
                ),
                f"Formula: {YIELD_LAW}.",
                (
                    "The three DT rows use this yield as delivered joules. "
                    "A scored catalog yield produces the factor yield/delivered. "
                    "A blank catalog yield receives this coverage directly."
                ),
                f"Source: `{YIELD_SOURCE}`.",
                "",
            ]
        )
    recorded = [row for row in rows if row.get("derived_side_factor") is not None]
    for row in recorded:
        lines.append(
            (
                "`{id}` records `{quantity}` `{factor}` via `{formula}`. "
                "Coverage stays `{coverage}` because the demonstrated efficiency is 0."
            ).format(
                id=_cell(row.get("id")),
                quantity=_cell(row.get("derived_side_quantity")),
                factor=f"{float(row['derived_side_factor']):.12g}",
                formula=_cell(row.get("derived_formula")),
                coverage=_cell(row.get("new_coverage")),
            )
        )
    if recorded:
        lines.append("")
    lines.extend(
        [
            f"Warp-shell target: `{target_j:.6g}` J, from `warp_shell_energy_j()`.",
            "",
            "## Rows that move",
            "",
        ]
    )
    if changed:
        lines.extend(_move_line(row) for row in changed)
    else:
        lines.append("None.")
    lines.extend(
        [
            "",
            "## Full table",
            "",
            "Old coverage is the phase-one score. New coverage is the derived yield",
            "divided by the shell where a yield is installed, and the phase-one",
            "coverage elsewhere. The factor is new/old when both coverages exist.",
            "",
            "| Rank | Id | Old coverage | New coverage | Factor | Derived quantity | Value | Source |",
            "| ---: | --- | ---: | ---: | ---: | --- | ---: | --- |",
        ]
    )
    for index, row in enumerate(rows, start=1):
        value = row.get("derived_value")
        shown = None if value is None else f"{float(value):.6g} {row.get('derived_unit') or ''}".strip()
        lines.append(
            "| {rank} | {id} | {old} | {new} | {factor} | {quantity} | {value} | {source} |".format(
                rank=index,
                id=_cell(row.get("id")),
                old=_cell(row.get("old_coverage")),
                new=_cell(row.get("new_coverage")),
                factor=_cell(row.get("coverage_factor")),
                quantity=_cell(row.get("derived_quantity")),
                value=_cell(shown),
                source=_cell(row.get("derived_source")),
            )
        )
    lines.extend(["", NEXT_TARGET, ""])
    return "\n".join(lines)


def render_derivations(
    clock: WalkClock,
    confined: ConfinedYield,
    derivations: tuple[Derivation, ...],
) -> str:
    """Formulas, steps, and the repo files that were read."""
    applied = ", ".join(
        f"`{item.mechanism_id}`" for item in derivations if item.formula_id == FORMULA_ID
    )
    if applied:
        assign_seconds = f"3. Assign the seconds to {applied}."
    else:
        assign_seconds = (
            "3. The seconds are τ in the confined DT yield, the D-He3 burn, "
            "and the p-B11 burn. No row takes the seconds as its own coverage factor."
        )
    yielded = ", ".join(
        f"`{item.mechanism_id}`"
        for item in derivations
        if item.formula_id == YIELD_FORMULA_ID and item.absolute_joules is not None
    )
    lines = [
        "# Phase-three derivations",
        "",
        "## Installed: confinement time",
        "",
        f"Formula id: `{FORMULA_ID}`.",
        f"Law: {LAW}.",
        "The step is one thermal gyroradius and the clock is one cyclotron period.",
        "The ion mass cancels. The historical Bohm denominator 16 is the comparison",
        "in `bohm_time_s` in the same file. walk/Bohm equals π/4 at every positive",
        "radius, temperature, and field.",
        "",
        f"Source file: `{SOURCE}`.",
        f"Inputs: `{INPUTS_SOURCE}`, stored again in `phase3/overrides.yaml`.",
        "",
        "Steps:",
        "",
        "1. Read `radius_m`, `temperature_keV`, and `field_T` from the mapping file.",
        "2. Call `gyro_walk_time_s`. `gamma_enters_the_time` is false, so the",
        "   Condensed_Matter scalar and `K` stay out of the product.",
        assign_seconds,
        "",
        f"Computed walk time: `{clock.walk_time_s:.6g}` s.",
        f"Computed Bohm time: `{clock.bohm_time_s:.6g}` s.",
        f"walk/Bohm: `{clock.walk_over_bohm:.6g}`.",
        "",
        "## Installed: confined DT yield",
        "",
        f"Formula id: `{YIELD_FORMULA_ID}`.",
        f"Law: {YIELD_LAW}.",
        f"Source: `{YIELD_SOURCE}`.",
        "The power density is the same product `evolve` writes at one knot:",
        "`n_D n_T <σv>(T) Q`. The volume is the sphere at the peak radius.",
        "The time is the walk time, so the yield is that power held for one",
        "confinement time. The schedule integral in the ledger is a different",
        "integral: it follows the prescribed radius schedule.",
        "",
        "Steps:",
        "",
        "1. Read the initial density, temperature, radius, and deuterium fraction",
        "   from the mapping file. They are the checked-in pulse.",
        "2. Compress to the peak radius with the monatomic sphere law in",
        "   `adiabatic_temperature_exponent`: `T = T0 (R0/R)^2` and `n = n0 (R0/R)^3`.",
        "3. Require that compressed temperature to match `temperature_keV`.",
        "4. Call `dt_reactivity_m3_s` and `dt_reaction().q_joule`.",
        f"5. Assign the joules to {yielded}.",
        "",
        f"Compressed density: `{confined.density_m3:.6g}` m^-3.",
        f"Compressed temperature: `{confined.temperature_kev:.6g}` keV.",
        f"Reactivity: `{confined.reactivity_m3_s:.6g}` m^3/s.",
        f"Volume: `{confined.volume_m3:.6g}` m^3.",
        f"Power: `{confined.power_w:.6g}` W.",
        f"Yield: `{confined.energy_j:.6g}` J.",
        "",
        "## Coverage arithmetic",
        "",
        "On a DT row that already has delivered joules, the factor is this yield",
        "divided by those joules. Coverage moves by that factor. On a DT row whose",
        "phase-one coverage is blank, the yield becomes the delivered joules and",
        "the coverage is yield divided by `warp_shell_energy_j()`.",
        "γ_FSOT is recorded on the confinement report and stays out of E.",
        "A coverage that falls from a positive value to 0 is reported as a factor",
        "of 0. A blank coverage that receives a yield has no factor.",
        "",
        "## Installed: pulse units, D-He3, and p-B11",
        "",
        "The confinement radius is the length already inside `τ = a² / (2 D)`.",
        "The ion thermal speed is the speed already inside that step,",
        "`v = sqrt(kT / m)`, with `m` the mean deuterium and tritium ion mass.",
        "The electron-volt is that same `kT`. Head scales linearly, tidal range",
        "scales as the square, wind speed scales as the cube, and the Casimir",
        "cavity scales as the inverse cube of the gap. The photovoltaic factor is",
        "the fraction of 6000 K blackbody energy above that `kT`. These are unit",
        "substitutions of the pulse. A dam, a tide, a rotor, a solar cell, and a",
        "Casimir cavity are not the plasma.",
        "",
        "D-He3 uses `dhe3_reactivity_m3_s` at the compressed temperature. The burn",
        "fraction is `n_He3 <σv> τ` with `n_He3` equal to half the compressed ion",
        "density, the pulse's stoichiometric split. Joules per reference kilogram",
        "are the catalog complete-burn density times that fraction. The checked-in",
        "plasma is deuterium-tritium. The D-He3 fit is evaluated at that peak.",
        "",
        "p-B11 uses `pb11_reactivity_m3_s` at the same temperature. That function",
        "is the Nevins and Swain low-temperature branch, Nucl. Fusion 40, 865",
        "(2000), reprinted as equations (9)-(15) by Tentori and Belloni, Nucl.",
        "Fusion 63, 086001 (2023), with the Nevins column of their table 1. The",
        "interval is 1 keV to 70 keV. The burn fraction is `n_B <σv> τ` with",
        "`n_B` equal to half the compressed ion density, one boron per proton",
        "on the catalog kilogram. The checked-in plasma is deuterium-tritium.",
        "The reactivity is p-11B. The high-temperature Padé starts at 50 keV",
        "and is not evaluated.",
        "",
    ]
    for item in derivations:
        if item.formula_id in {FORMULA_ID, YIELD_FORMULA_ID}:
            continue
        law = _LAW_BY_FORMULA.get(item.formula_id, "")
        lines.extend(
            [
                f"### `{item.mechanism_id}`",
                "",
                f"Formula id: `{item.formula_id}`.",
                f"Law: {law}.",
                f"Source: `{item.source}`.",
                item.detail,
                "",
            ]
        )
    lines.extend(
        [
            "## Other catalog levers",
            "",
            "Phase two already scored the dimensionless efficiencies against captured",
            "domain scalars. This pass leaves those efficiencies on the catalog values.",
            "Fission burnup on the ranked unit is already 1. The first-collision tally",
            "in the hybrid ledger is a fissions-per-source surrogate for one pulse, and",
            "the literature burnup note in the fission files is a phase-one citation.",
            "Neither one is installed here as a burnup of that ranked kilogram.",
            "",
            "Dimensionful parameters:",
            "",
        ]
    )
    for mechanism_id, parameter in unmapped_parameters():
        lines.append(
            f"- `{mechanism_id}`: {parameter} ({_dimensionful_status(mechanism_id)})"
        )
    lines.extend(
        [
            "",
            "## Inspected and left unused",
            "",
            "### FSOT-2.1-Lean `_fsot_scaled`",
            "",
            "File: `FSOT-2.1-Lean/scripts/tier_gap_fill_lib.py`, function `_fsot_scaled`.",
            "The expression is `measured * (1 + abs(scalar) * factor)`, with default",
            "factor `0.001`. `scripts/tier71_fusion_lab_lib.py` calls it on measured",
            "triple products with factor `1e-12`, and on reaction energetics with",
            "factor `1e-6`. The relative move is `|S|` times that factor. The panel",
            "relays a measured confinement time from its public anchors and then",
            "scales the triple product. The walk clock is the installed time.",
            "",
            "### Engine chemistry",
            "",
            "File: `FSOT-2.1-Lean/vendor/fsot_compute.py`, pin prefix `aeb2adad6e80f487`.",
            "The file was read and left unedited. `chemistry_bond_energies` returns",
            "kJ/mol for H-H, C-H, C-C, C=C, C≡C, O-H, N-H, N≡N, O=O, and F-F.",
            "`chemistry_bond_lengths` returns angstroms, `chemistry_ionization`",
            "returns eV, `chemistry_radii` returns angstroms, and `chemistry_molecular`",
            "returns kelvin, debye, and pH. The bond table has no C=O entry. Those",
            "comparisons stay inside the engine. They were not written onto a",
            "mechanism's ranked reference unit.",
            "",
            "### FSOT-Genetics",
            "",
            "File: `FSOT-Genetics/formulas/FSOT_PROTEIN_DERIVATIONS.md`, formulas",
            "F01 through F09. They are amino-acid trinary, contact, and chemistry",
            "amplitude scalars. That repository is on pin D1D38A. This pass keeps",
            "the live engine on pin AEB2AD and does not import the Genetics vendor.",
            "",
            "### FSOT-2.1-Neural",
            "",
            "No checkout is present at `Desktop/FSOT-2.1-Neural`,",
            "`fsot-xval/siblings/FSOT-2.1-Neural`, or",
            "`Desktop/fsot neuron family/FSOT-2.1-Neural`. The public repository",
            "is the codon neural substrate. No seconds law was read from a local tree.",
            "",
            "### FSOT-2.1-Cpp",
            "",
            "The C++ repository states the same scalar definitions the apply script",
            "already checks in `include/fsot/core.hpp`. The confinement time installed",
            "here is `gyro_walk_time_s`.",
            "",
            "## Next derivation target",
            "",
            NEXT_TARGET,
            "",
        ]
    )
    return "\n".join(lines)


def build_ranking(catalog: Path, overrides: Path, target_j: float) -> tuple[list[dict], WalkClock, str, str]:
    loaded = load_overrides(overrides, catalog)
    rows = sort_rows(annotate(rank_directory(catalog, target_j), loaded.derivations, target_j))
    ranking = render_ranking(rows, target_j, loaded.clock, loaded.confined_yield)
    derivations_text = render_derivations(loaded.clock, loaded.confined_yield, loaded.derivations)
    return rows, loaded.clock, ranking, derivations_text


def build_gaps(catalog: Path, target_j: float) -> str:
    module, _digest = load_engine()
    scalars = capture_scalars(module)
    return compose_gaps(map_directory(catalog, scalars, target_j), target_j)


def main(argv: list[str] | None = None) -> int:
    del argv
    catalog, ranking_path, derivations_path, gaps_path = default_paths()
    target = float(warp_shell_energy_j())
    _rows, _clock, ranking, derivations_text = build_ranking(
        catalog,
        default_paths_overrides(),
        target,
    )
    gaps_text = build_gaps(catalog, target)
    for path, text in (
        (ranking_path, ranking),
        (derivations_path, derivations_text),
        (gaps_path, gaps_text),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    sys.stdout.write(ranking)
    return 0


def default_paths_overrides() -> Path:
    return Path(__file__).resolve().parent / "overrides.yaml"


if __name__ == "__main__":
    raise SystemExit(main())
