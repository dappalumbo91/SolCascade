"""Compare the converters already filed on one U-235 fission kilogram.

Delivered joules stay the catalog product of energy density and demonstrated
efficiency. The fragment share is the cited page's round split of the 200 MeV
ledger. It is heat inside that ledger. A foundational seed is recorded beside
the share and is not multiplied into the joules. The checked-in pulse is not
rewritten.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from solcascade.constants import FISSION_RECOVERABLE_ENERGY_MEV
from solcascade.energy_baseline.phase2.map_fsot import (
    PIN_PREFIX,
    capture_scalars,
    load_engine,
)
from solcascade.production.refine import catalog_record, project_root

THERMAL_ID = "fission-thermal"
THERMIONIC_ID = "thermionic"
FAST_ID = "fission-fast"
LEDGER_MEV = FISSION_RECOVERABLE_ENERGY_MEV
ADDED_MEV = 30.0
FRAGMENT_SHARE = (LEDGER_MEV - ADDED_MEV) / LEDGER_MEV
_ATOMS = ("pi", "e", "phi", "gamma", "G")


@dataclass(frozen=True)
class Converter:
    """One catalog converter on the shared U-235 kilogram."""

    mechanism_id: str
    name: str
    trl: int
    reference_unit: str
    efficiency: float
    ceiling: float
    energy_density_j: float
    delivered_j: float
    ceiling_delivered_j: float
    input_j: float
    output_over_input: float


@dataclass(frozen=True)
class SeedGap:
    """One foundational product beside the page share."""

    name: str
    value: float
    relative_gap: float


@dataclass(frozen=True)
class PathsReport:
    thermal: Converter
    thermionic: Converter
    same_ids: tuple[str, ...]
    fast_reference_unit: str
    governing_source: str
    fragment_share: float
    added_mev: float
    ledger_mev: float
    kinetic_inventory_j: float
    nearest: SeedGap
    next_seed: SeedGap
    seeds_named: int
    seeds_ranked: int
    nuclear_scalar: float
    thermodynamics_scalar: float
    scalar_gap: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_paths.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _converter(mechanism_id: str) -> tuple[Converter, dict]:
    record = catalog_record(mechanism_id)
    if record.get("fsot_status") != "open_gap" or record.get("fsot_scalar") is not None:
        raise RuntimeError(f"{mechanism_id} carries a scalar mapping")
    if record.get("chain_reaction") is not True:
        raise RuntimeError(f"{mechanism_id} is not the fission chain")
    if record.get("self_sustaining") is not True:
        raise RuntimeError(f"{mechanism_id} is not a self-sustaining chain")
    if record.get("continuous_external_drive") is not False:
        raise RuntimeError(f"{mechanism_id} asks for a continuous external drive")
    efficiency = float(record["demonstrated_efficiency"])
    ceiling = float(record["efficiency_ceiling"])
    density = float(record["energy_density_j"])
    incoming = float(record["input_energy_j"])
    if not 0.0 < efficiency <= 1.0:
        raise RuntimeError(f"{mechanism_id} efficiency is outside (0, 1]")
    if ceiling < efficiency or ceiling > 1.0:
        raise RuntimeError(f"{mechanism_id} ceiling is below the demonstrated efficiency")
    if density <= 0.0 or incoming <= 0.0:
        raise RuntimeError(f"{mechanism_id} joules are not positive")
    delivered = density * efficiency
    return (
        Converter(
            mechanism_id=mechanism_id,
            name=str(record["name"]),
            trl=int(record["trl"]),
            reference_unit=str(record["reference_unit"]),
            efficiency=efficiency,
            ceiling=ceiling,
            energy_density_j=density,
            delivered_j=delivered,
            ceiling_delivered_j=density * ceiling,
            input_j=incoming,
            output_over_input=delivered / incoming,
        ),
        record,
    )


def _ids_at_density(density: float) -> tuple[str, ...]:
    data = project_root() / "src" / "solcascade" / "energy_baseline" / "data"
    found: list[str] = []
    for path in sorted(data.glob("*.yaml")):
        record = catalog_record(path.stem)
        value = record.get("energy_density_j")
        if value is None:
            continue
        if float(value) == density:
            found.append(str(record["id"]))
    return tuple(found)


def _no_collector_key(record: dict) -> None:
    for key in record:
        if "collector" in str(key).lower():
            raise RuntimeError(f"{record.get('id')} files a collector field")


def _share_seeds(module) -> tuple[tuple[str, float], ...]:
    """Products of the five foundational seeds. No free coefficient."""
    values = {
        "pi": module.PI,
        "e": module.E,
        "phi": module.PHI,
        "gamma": module.GAMMA,
        "G": module.G_CAT,
    }
    named: list[tuple[str, float]] = []
    for name in _ATOMS:
        value = values[name]
        named.append((name, float(value)))
        named.append((f"{name}^2", float(value**2)))
        named.append((f"{name}^-1", float(1 / value)))
    for index, left_name in enumerate(_ATOMS):
        left = values[left_name]
        for right_name in _ATOMS[index + 1 :]:
            right = values[right_name]
            named.append((f"{left_name}*{right_name}", float(left * right)))
            named.append((f"{left_name}/{right_name}", float(left / right)))
            named.append((f"{right_name}/{left_name}", float(right / left)))
    return tuple(named)


def _nearest_seeds(module, share: float) -> tuple[SeedGap, SeedGap, int, int]:
    named = _share_seeds(module)
    ranked: list[SeedGap] = []
    for name, value in named:
        if not math.isfinite(value):
            raise RuntimeError(f"seed {name} is not finite")
        if not 0.0 < value < 1.0:
            continue
        gap = abs(value - share) / share
        ranked.append(SeedGap(name=name, value=value, relative_gap=gap))
    ranked.sort(key=lambda item: (item.relative_gap, item.name))
    if len(ranked) < 2:
        raise RuntimeError("the share comparison needs two seeds inside (0, 1)")
    if ranked[0].relative_gap <= 0.0:
        raise RuntimeError("a seed reproduces the page share")
    return ranked[0], ranked[1], len(named), len(ranked)


def _witness() -> tuple[str, float, dict[str, float], SeedGap, SeedGap, int, int]:
    module, digest = load_engine()
    if not str(digest).startswith(PIN_PREFIX):
        raise RuntimeError(f"engine pin is {str(digest)[:16]}")
    scalars = capture_scalars(module)
    valve = float(module.BETA)
    if valve <= 0.0 or valve >= 1.0:
        raise RuntimeError("engine valve is outside the recorded range")
    nearest, following, named_count, ranked_count = _nearest_seeds(module, FRAGMENT_SHARE)
    return str(digest)[:16], valve, scalars, nearest, following, named_count, ranked_count


def build_paths() -> PathsReport:
    if LEDGER_MEV <= ADDED_MEV or ADDED_MEV <= 0.0:
        raise RuntimeError("the page addendum does not sit inside the ledger")
    if not 0.0 < FRAGMENT_SHARE < 1.0:
        raise RuntimeError("the fragment share is outside (0, 1)")
    thermal, thermal_record = _converter(THERMAL_ID)
    thermionic, thermionic_record = _converter(THERMIONIC_ID)
    _no_collector_key(thermal_record)
    _no_collector_key(thermionic_record)
    if "U-235" not in thermal.reference_unit or "U-235" not in thermionic.reference_unit:
        raise RuntimeError("the pair is not on the U-235 kilogram")
    if thermal.energy_density_j != thermionic.energy_density_j:
        raise RuntimeError("the pair does not share one energy density")
    if thermal.input_j != thermionic.input_j:
        raise RuntimeError("the enrichment estimates differ")
    if thermionic.delivered_j >= thermal.delivered_j:
        raise RuntimeError("the thermionic row meets the thermal delivery")
    if thermionic.ceiling_delivered_j >= thermal.delivered_j:
        raise RuntimeError("the thermionic ceiling meets the thermal delivery")
    if thermionic.ceiling >= thermal.efficiency:
        raise RuntimeError("the thermionic ceiling reaches the thermal efficiency")
    same_ids = _ids_at_density(thermal.energy_density_j)
    if same_ids != (THERMAL_ID, THERMIONIC_ID):
        raise RuntimeError("the shared density is not exactly the thermal and thermionic rows")
    fast = catalog_record(FAST_ID)
    _no_collector_key(fast)
    fast_unit = str(fast["reference_unit"])
    if "Pu-239" not in fast_unit:
        raise RuntimeError("the fast reference unit is not Pu-239")
    if float(fast["energy_density_j"]) == thermal.energy_density_j:
        raise RuntimeError("the fast row shares the U-235 density")
    source = str(thermal_record["governing_source"])
    if "physics-of-nuclear-energy" not in source:
        raise RuntimeError("the thermal governing source is not the physics page")
    if "TOPAZ" not in str(thermionic_record.get("notes") or ""):
        raise RuntimeError("the thermionic notes do not name TOPAZ")
    digest, valve, scalars, nearest, following, named_count, ranked_count = _witness()
    nuclear = float(scalars["Nuclear_Physics"])
    thermo = float(scalars["Thermodynamics"])
    if not math.isfinite(nuclear) or not math.isfinite(thermo):
        raise RuntimeError("a witness scalar is not finite")
    scalar_gap = abs(nuclear - FRAGMENT_SHARE) / FRAGMENT_SHARE
    if scalar_gap <= nearest.relative_gap:
        raise RuntimeError("the domain scalar is at least as close as the nearest seed")
    inventory = FRAGMENT_SHARE * thermal.energy_density_j
    if inventory >= thermal.energy_density_j:
        raise RuntimeError("the fragment inventory meets the whole density")
    return PathsReport(
        thermal=thermal,
        thermionic=thermionic,
        same_ids=same_ids,
        fast_reference_unit=fast_unit,
        governing_source=source,
        fragment_share=FRAGMENT_SHARE,
        added_mev=ADDED_MEV,
        ledger_mev=LEDGER_MEV,
        kinetic_inventory_j=inventory,
        nearest=nearest,
        next_seed=following,
        seeds_named=named_count,
        seeds_ranked=ranked_count,
        nuclear_scalar=nuclear,
        thermodynamics_scalar=thermo,
        scalar_gap=scalar_gap,
        gamma_fsot=float(scalars["gamma_FSOT"]),
        valve_beta=valve,
        pin_prefix=PIN_PREFIX,
        digest_prefix=digest,
    )


def _row(converter: Converter) -> str:
    return (
        f"| {converter.name} | {converter.trl} | {converter.efficiency:.6g} | "
        f"{converter.ceiling:.6g} | {converter.delivered_j:.6g} | "
        f"{converter.ceiling_delivered_j:.6g} | {converter.output_over_input:.6g} |"
    )


def render() -> str:
    report = build_paths()
    thermal = report.thermal
    thermionic = report.thermionic
    nearest = report.nearest
    following = report.next_seed
    lines = [
        "# Fission paths",
        "",
        "The thermal plant and the thermionic row are the two catalog converters",
        "on one kilogram of U-235 fissioned completely.",
        "Delivered joules are the catalog product of that shared density and each",
        "row's demonstrated efficiency.",
        "",
        "| Converter | TRL | Efficiency | Ceiling | Delivered (J) | At the ceiling (J) | Output / enrichment estimate |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        _row(thermal),
        _row(thermionic),
        "",
        f"The thermal reference unit is {thermal.reference_unit}.",
        f"The thermionic reference unit is {thermionic.reference_unit}.",
        f"The shared density is `{_num(thermal.energy_density_j)}` J per kilogram.",
        "The thermionic row is the other filed converter on this same U-235 kilogram.",
        "The thermionic notes name TOPAZ as a flown thermionic converter.",
        f"The thermionic row delivers `{_num(thermionic.delivered_j)}` J.",
        "The thermionic row delivers less electricity than the thermal plant.",
        f"The delivered ratio, thermionic over thermal, is `{_num(thermionic.delivered_j / thermal.delivered_j)}`.",
        "Even at its file ceiling the thermionic row delivers less than the thermal plant.",
        "Both rows use the same enrichment estimate.",
        f"The fast reference unit is {report.fast_reference_unit}.",
        "The fast row stays on that unit.",
        "The refinement inside the operating plants is the file efficiency ceiling.",
        "",
        "## Fragment inventory",
        "",
        f"The thermal governing source is {report.governing_source}.",
        "That page states that about 85 percent of the released energy starts as fission-fragment kinetic energy.",
        "In solid fuel the fragments stop in a microscopic distance and the energy becomes heat.",
        f"Neutron, gamma, and delayed releases add about `{_num(report.added_mev)}` MeV on a total of about `{_num(report.ledger_mev)}` MeV.",
        f"Subtracting that addendum from the ledger leaves the share `{_num(report.fragment_share)}`.",
        "The catalog density stays on the 200 MeV ledger.",
        f"The fragment kinetic inventory on the ranked kilogram is `{_num(report.kinetic_inventory_j)}` J.",
        "The fragment share is already inside that density.",
        "The Rankine efficiency multiplies the whole ranked density.",
        "The fragment inventory is heat in solid fuel.",
        "The fragment inventory is not electricity.",
        "The catalog files no collector efficiency for a fragment that leaves the fuel.",
        "The operating electricity on this kilogram is the Rankine product.",
        "",
        "## FSOT record",
        "",
        f"Pin `{report.pin_prefix}` matches the loaded engine prefix `{report.digest_prefix}`.",
        "The scalars below are dimensionless.",
        "A domain scalar is not a plant efficiency.",
        "A domain scalar is not the fragment share.",
        "The nuclear scalar stays out of the fission joules.",
        f"Nuclear_Physics: `{_num(report.nuclear_scalar)}`.",
        f"Thermodynamics: `{_num(report.thermodynamics_scalar)}`.",
        f"The domain scalar's relative gap to the share is `{_num(report.scalar_gap)}`.",
        "The domain scalar sits farther from the share than the nearest seed.",
        f"gamma_FSOT: `{_num(report.gamma_fsot)}`.",
        "gamma_FSOT does not multiply the plant.",
        "The engine valve stays out of the fission joules.",
        f"Engine valve: `{_num(report.valve_beta)}`.",
        "",
        "The comparison uses each of pi, e, phi, gamma, and G, its square, its reciprocal,",
        "and each pair's product and two quotients.",
        f"`{report.seeds_named}` products are named and `{report.seeds_ranked}` lie strictly between zero and one.",
        f"The nearest seed is `{nearest.name}` at `{_num(nearest.value)}`.",
        f"Its relative gap to the share is `{_num(nearest.relative_gap)}`.",
        f"The next seed is `{following.name}` at `{_num(following.value)}`.",
        f"Its relative gap to the share is `{_num(following.relative_gap)}`.",
        "The page figure is a round share.",
        "The nearest seed is not installed as the share.",
        "",
        "The checked-in schedule is not rewritten.",
        "",
    ]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    del argv
    text = render()
    output = default_output_path()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text, encoding="utf-8")
    print(text, end="" if text.endswith("\n") else "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
