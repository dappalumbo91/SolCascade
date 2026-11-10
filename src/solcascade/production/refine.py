"""Recompute the production ledger on the checked-in pulse.

The four coil draws are the supplies locked in ``scripts/apply_fsot.py``.
Fusion energy is the schedule integral of Bosch-Hale power. The gyro-walk
time is reported beside that integral. It is not a multiplier, and gamma
does not enter it.

The baseline store stays closed on this fuel load: the complete-burn ceiling
is below that draw. Hotter trajectories stay inside the Bosch-Hale window.
A larger inventory is more fuel at the same density and the same burn
fraction. It is the mass a JET pulse or a Pu-239 kilogram would require.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import yaml

from solcascade.assemble import build_coils, build_pulse
from solcascade.constants import (
    BOSCH_HALE_DT_T_MAX_KEV,
    BOSCH_HALE_DT_T_MIN_KEV,
    MASS_EXCESS_DEUTERIUM_KEV,
    MASS_EXCESS_TRITIUM_KEV,
    dt_reaction,
    ion_mass_kg,
)
from solcascade.conversion import SupplyLedger, convert
from solcascade.interfaces import NullProvider, apply_provider, collect_slots
from solcascade.load import load_problem
from solcascade.plasma.metrics import evolve, gyro_walk_time_s
from solcascade.plasma.reactivity import dt_reactivity_m3_s, pb11_reactivity_m3_s

# Draws locked by scripts/apply_fsot.py _check_reference. This module does not
# rerun the coil solver.
SUPPLIES: tuple[tuple[str, float], ...] = (
    ("baseline", 1423874.990828079),
    ("same_b_flat", 2106.2160000000003),
    ("same_b_optimal", 2020.4381637452752),
    ("same_b_helmholtz", 827.5009597961755),
)

PEAK_TEMPERATURES_KEV: tuple[float, ...] = (10.0, 20.0, 50.0, 100.0)


@dataclass(frozen=True)
class Problem:
    config: dict
    slots: dict[str, float]
    provenance: dict[str, str]
    coils: object
    schedule: object


@dataclass(frozen=True)
class Fuel:
    pairs: float
    ceiling_j: float
    mass_kg: float
    deuterium_kg: float
    tritium_kg: float


@dataclass(frozen=True)
class PulseYield:
    temperature0_kev: float
    peak_temperature_kev: float
    peak_radius_m: float
    peak_density_m3: float
    peak_field_t: float
    peak_beta: float
    peak_power_w: float
    fusion_j: float
    burn_fraction: float
    walk_time_s: float
    walk_yield_j: float
    walk_burn_fraction: float


@dataclass(frozen=True)
class TopologyScore:
    name: str
    supply_j: float
    fusion_j: float
    burn_fraction: float
    electric_j: float
    net_j: float
    ceiling_over_supply: float
    closed: bool


@dataclass(frozen=True)
class NeutronAccount:
    label: str
    fusion_j: float
    fuel_mass_kg: float
    burn_fraction: float
    tritium_loaded_kg: float
    tritium_consumed_kg: float
    neutrons: float
    neutron_energy_j: float
    neutron_energy_kev: float


def project_root() -> Path:
    package = Path(__file__).resolve().parent
    return package.parents[2]


def default_config_path() -> Path:
    return project_root() / "configs" / "minimal_pipeline.yaml"


def default_output_path() -> Path:
    return project_root() / "runs" / "production" / "refinement.md"


def catalog_record(mechanism_id: str) -> dict:
    path = (
        project_root()
        / "src"
        / "solcascade"
        / "energy_baseline"
        / "data"
        / f"{mechanism_id}.yaml"
    )
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(loaded, dict) or loaded.get("id") != mechanism_id:
        raise ValueError(f"catalog record {mechanism_id} is missing")
    return loaded


def load_checked_in(path: Path | None = None) -> Problem:
    config_path = default_config_path() if path is None else Path(path)
    config, materials = load_problem(config_path)
    slots, provenance = apply_provider(collect_slots(config, materials), NullProvider())
    coils = build_coils(config, slots)
    schedule = build_pulse(config, slots, provenance)
    return Problem(
        config=config,
        slots=slots,
        provenance=provenance,
        coils=coils,
        schedule=schedule,
    )


def fuel_of(problem: Problem) -> Fuel:
    radius = float(problem.schedule.plasma_radius_m[0])
    density = float(problem.slots["plasma.ion_density_m3"])
    fraction = float(problem.slots["plasma.deuterium_fraction"])
    if not 0.0 < fraction < 1.0:
        raise ValueError("deuterium fraction must be between 0 and 1")
    volume = (4.0 / 3.0) * math.pi * radius ** 3
    pairs = fraction * density * volume
    reaction = dt_reaction()
    deuterium = ion_mass_kg(2, MASS_EXCESS_DEUTERIUM_KEV, 1)
    tritium = ion_mass_kg(3, MASS_EXCESS_TRITIUM_KEV, 1)
    return Fuel(
        pairs=pairs,
        ceiling_j=pairs * float(reaction.q_joule),
        mass_kg=pairs * (deuterium + tritium),
        deuterium_kg=pairs * deuterium,
        tritium_kg=pairs * tritium,
    )


def electric_per_fusion_j(problem: Problem) -> float:
    """Joules of electricity per fusion joule on the checked-in conversion."""
    conversion = problem.config["conversion"]
    if conversion["shell_channel"] != "none":
        raise ValueError("production score expects shell_channel none")
    reaction = dt_reaction()
    supply = SupplyLedger(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    ledger = convert(
        fusion_j=1.0,
        alpha_fraction=float(reaction.alpha_fraction),
        neutron_deposited_fraction=0.0,
        alpha_escape_fraction=float(conversion["alpha_escape_fraction"]),
        alpha_wall_fraction=float(conversion["alpha_wall_fraction"]),
        alpha_direct_efficiency=float(conversion["alpha_direct_efficiency"]),
        shell_channel="none",
        carnot_fraction=float(conversion["carnot_fraction"]),
        hot_k=600.0,
        cold_k=float(conversion["cold_temperature_K"]),
        supply=supply,
    )
    return float(ledger.electric_gross_j)


def schedule_yield(problem: Problem, temperature0_kev: float, fuel: Fuel) -> PulseYield:
    """Integrate one imposed radius schedule at one initial temperature."""
    if temperature0_kev < BOSCH_HALE_DT_T_MIN_KEV:
        raise ValueError("initial temperature is below the Bosch-Hale window")
    history = evolve(
        problem.coils,
        problem.schedule,
        float(problem.slots["plasma.ion_density_m3"]),
        float(temperature0_kev),
        float(problem.slots["plasma.deuterium_fraction"]),
        str(problem.config["compression"]["flux_conservation"]),
    )
    peak = max(
        history.points,
        key=lambda point: (point.ion_temperature_kev, point.fusion_power_w),
    )
    if peak.ion_temperature_kev > BOSCH_HALE_DT_T_MAX_KEV:
        raise ValueError("peak temperature is above the Bosch-Hale window")
    # Touch the reactivity at both ends so a trajectory outside the fit fails
    # in this function rather than inside the integrator's private call.
    dt_reactivity_m3_s(float(history.points[0].ion_temperature_kev))
    dt_reactivity_m3_s(float(peak.ion_temperature_kev))
    fusion = history.integrate("fusion_power_w")
    walk = gyro_walk_time_s(
        peak.radius_m, peak.ion_temperature_kev, peak.internal_field_t
    )
    walk_yield = peak.fusion_power_w * walk
    ceiling = fuel.ceiling_j
    return PulseYield(
        temperature0_kev=float(temperature0_kev),
        peak_temperature_kev=float(peak.ion_temperature_kev),
        peak_radius_m=float(peak.radius_m),
        peak_density_m3=float(peak.ion_density_m3),
        peak_field_t=float(peak.internal_field_t),
        peak_beta=float(peak.beta),
        peak_power_w=float(peak.fusion_power_w),
        fusion_j=fusion,
        burn_fraction=fusion / ceiling,
        walk_time_s=float(walk),
        walk_yield_j=float(walk_yield),
        walk_burn_fraction=walk_yield / ceiling,
    )


def temperature_scan(problem: Problem, fuel: Fuel) -> tuple[PulseYield, ...]:
    """Checked-in trajectory, then hotter peaks still inside the fit window."""
    radius0 = float(problem.schedule.plasma_radius_m[0])
    radius_peak = float(min(problem.schedule.plasma_radius_m))
    ratio = (radius0 / radius_peak) ** 2
    checked = float(problem.slots["plasma.ion_temperature_keV"])
    rows = [schedule_yield(problem, checked, fuel)]
    for peak in PEAK_TEMPERATURES_KEV:
        initial = peak / ratio
        if initial <= checked:
            continue
        rows.append(schedule_yield(problem, initial, fuel))
    return tuple(rows)


def score_topologies(
    fusion_j: float,
    fuel: Fuel,
    electric_per_j: float,
) -> tuple[TopologyScore, ...]:
    scores = []
    for name, supply in SUPPLIES:
        electric = fusion_j * electric_per_j
        scores.append(
            TopologyScore(
                name=name,
                supply_j=float(supply),
                fusion_j=float(fusion_j),
                burn_fraction=fusion_j / fuel.ceiling_j,
                electric_j=electric,
                net_j=electric - float(supply),
                ceiling_over_supply=fuel.ceiling_j / float(supply),
                closed=fuel.ceiling_j < float(supply),
            )
        )
    return tuple(scores)


def inventory_mass_kg(target_j: float, row: PulseYield, fuel: Fuel) -> float:
    """Fuel mass that releases ``target_j`` at this row's burn fraction."""
    if row.fusion_j <= 0.0:
        raise ValueError("inventory scale needs a positive yield")
    return fuel.mass_kg * (float(target_j) / row.fusion_j)


def neutron_account(
    label: str,
    fusion_j: float,
    fuel_mass_kg: float,
    burn_fraction: float,
    fuel: Fuel,
) -> NeutronAccount:
    reaction = dt_reaction()
    if fuel.mass_kg <= 0.0:
        raise ValueError("fuel mass must be positive")
    scale = float(fuel_mass_kg) / fuel.mass_kg
    loaded = fuel.tritium_kg * scale
    neutrons = float(fusion_j) / float(reaction.q_joule)
    consumed = neutrons * ion_mass_kg(3, MASS_EXCESS_TRITIUM_KEV, 1)
    return NeutronAccount(
        label=label,
        fusion_j=float(fusion_j),
        fuel_mass_kg=float(fuel_mass_kg),
        burn_fraction=float(burn_fraction),
        tritium_loaded_kg=loaded,
        tritium_consumed_kg=consumed,
        neutrons=neutrons,
        neutron_energy_j=neutrons * float(reaction.neutron_energy_joule),
        neutron_energy_kev=float(reaction.neutron_energy_kev),
    )


def pb11_burn_fraction(row: PulseYield, problem: Problem) -> float:
    """Partner-ion burn at one peak. The low-temperature fit rejects itself outside 1–70 keV."""
    fraction = float(problem.slots["plasma.deuterium_fraction"])
    partner = (1.0 - fraction) * row.peak_density_m3
    reactivity = pb11_reactivity_m3_s(row.peak_temperature_kev)
    return partner * reactivity * row.walk_time_s


def _num(value: float) -> str:
    return f"{value:.6g}"


def _joule(value: float) -> str:
    return f"{value:.6g}"


def render(problem: Problem | None = None) -> str:
    problem = load_checked_in() if problem is None else problem
    fuel = fuel_of(problem)
    per_joule = electric_per_fusion_j(problem)
    scan = temperature_scan(problem, fuel)
    checked = scan[0]
    hottest = max(scan, key=lambda row: row.fusion_j)
    walk_at_20 = next(row for row in scan if math.isclose(row.peak_temperature_kev, 20.0, rel_tol=1e-6))
    scores = score_topologies(checked.fusion_j, fuel, per_joule)
    helmholtz = next(row for row in scores if row.name == "same_b_helmholtz")
    jet = catalog_record("fusion-dt-magnetic")
    fission = catalog_record("fission-fast")
    jet_j = float(jet["energy_density_j"])
    fission_j = float(fission["energy_density_j"])
    hot_jet_mass = inventory_mass_kg(jet_j, hottest, fuel)
    hot_fission_mass = inventory_mass_kg(fission_j, hottest, fuel)
    checked_jet_mass = inventory_mass_kg(jet_j, checked, fuel)
    checked_fission_mass = inventory_mass_kg(fission_j, checked, fuel)
    reaction = dt_reaction()
    accounts = (
        neutron_account("checked-in schedule", checked.fusion_j, fuel.mass_kg, checked.burn_fraction, fuel),
        neutron_account("hottest schedule", hottest.fusion_j, fuel.mass_kg, hottest.burn_fraction, fuel),
        neutron_account("JET mass at the hottest burn fraction", jet_j, hot_jet_mass, hottest.burn_fraction, fuel),
        neutron_account(
            "Pu-239 kilogram at the hottest burn fraction",
            fission_j,
            hot_fission_mass,
            hottest.burn_fraction,
            fuel,
        ),
    )
    pb11 = pb11_burn_fraction(checked, problem)
    dwell = helmholtz.supply_j / hottest.peak_power_w
    lines = [
        "# Production refinement",
        "",
        "The coil draws are the four supplies locked in `scripts/apply_fsot.py`.",
        "Fusion energy is the integral of Bosch-Hale power over the checked-in",
        "radius schedule. Electricity is the checked-in direct converter:",
        "escaped alphas only, `shell_channel` none. The warp-shell rest energy",
        "is not this score.",
        "",
        f"The compressed sphere holds `{_num(fuel.mass_kg)}` kg of deuterium-tritium,",
        f"`{_num(fuel.tritium_kg)}` kg of it tritium. Complete burn releases",
        f"`{_joule(fuel.ceiling_j)}` J.",
        "",
        "## Topologies",
        "",
        "Checked-in schedule. A closed row has a coil draw above the complete-burn ceiling.",
        "",
        "| Topology | Coil supply (J) | Fusion (J) | Burn fraction | Electric (J) | Net (J) | Ceiling / supply |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in scores:
        lines.append(
            f"| {row.name} | {row.supply_j:.6g} | {row.fusion_j:.6g} | "
            f"{row.burn_fraction:.6g} | {row.electric_j:.6g} | {row.net_j:.6g} | "
            f"{row.ceiling_over_supply:.6g} |"
        )
    lines.extend(
        [
            "",
            "The baseline store stays closed on this fuel load. The flat, optimal,",
            "and Helmholtz same-B draws sit under the ceiling. At the checked-in",
            "yield every net is negative.",
            "",
            "## Schedule on the Bosch-Hale curve",
            "",
            "Initial temperature sets the adiabatic trajectory. Every temperature",
            f"on a row stays inside {BOSCH_HALE_DT_T_MIN_KEV:g}–{BOSCH_HALE_DT_T_MAX_KEV:g} keV.",
            "The walk yield is peak power times `gyro_walk_time_s` at that peak.",
            "It is a different integral from the schedule.",
            "",
            "| T0 (keV) | Peak (keV) | Fusion (J) | Burn fraction | Peak beta | Walk time (s) | Walk yield (J) |",
            "| ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in scan:
        lines.append(
            f"| {row.temperature0_kev:.6g} | {row.peak_temperature_kev:.6g} | "
            f"{row.fusion_j:.6g} | {row.burn_fraction:.6g} | {row.peak_beta:.6g} | "
            f"{row.walk_time_s:.6g} | {row.walk_yield_j:.6g} |"
        )
    lines.extend(
        [
            "",
            f"The hottest schedule in the window releases `{_joule(hottest.fusion_j)}` J",
            f"at a peak of `{_num(hottest.peak_temperature_kev)}` keV. Its burn fraction",
            f"is `{_num(hottest.burn_fraction)}`. The Helmholtz draw is",
            f"`{_joule(helmholtz.supply_j)}` J, and the burn fraction that matches it",
            f"on this ceiling is `{_num(helmholtz.supply_j / fuel.ceiling_j)}`.",
            f"The schedule remains short of that draw by a factor of",
            f"`{_num(helmholtz.supply_j / hottest.fusion_j)}`.",
            "",
            f"Matching the Helmholtz draw from the hottest peak power would take",
            f"`{_num(dwell)}` s at that power. The gyro-walk time at that peak is",
            f"`{_num(hottest.walk_time_s)}` s, and the walk yield is",
            f"`{_joule(hottest.walk_yield_j)}` J. The walk yield at 100 keV is below",
            f"the walk yield at 20 keV (`{_joule(walk_at_20.walk_yield_j)}` J).",
            "Heating shortens the walk. The confinement time does not close the",
            "Helmholtz gap.",
            "",
            f"Peak beta on the checked-in schedule is `{_num(checked.peak_beta)}`.",
            "It rises with temperature. Beta above 1 means the thermal pressure",
            "exceeds the magnetic pressure of the frozen-flux field. The schedule",
            "integral is the prescribed compression.",
            "",
            "## Inventory",
            "",
            "Mass below is more fuel at the same density, temperature, and burn",
            "fraction. The joule target is the catalog energy of that reference unit.",
            "",
            f"- One JET pulse, `{jet['reference_unit']}`: `{_joule(jet_j)}` J.",
            f"  Checked-in burn fraction needs `{_num(checked_jet_mass)}` kg.",
            f"  The hottest burn fraction needs `{_num(hot_jet_mass)}` kg.",
            f"- One Pu-239 kilogram, `{fission['reference_unit']}`:",
            f"  `{_joule(fission_j)}` J of fission energy, before the catalog efficiency.",
            f"  Checked-in burn fraction needs `{_num(checked_fission_mass)}` kg of DT.",
            f"  The hottest burn fraction needs `{_num(hot_fission_mass)}` kg of DT.",
            "",
            "A denser sphere would raise the schedule burn fraction in proportion",
            "to density, and it would raise beta by the same factor. The checked-in",
            "peak is already above beta 1, so that direction moves farther from",
            "magnetic balance.",
            "",
            "## Tritium and DT neutrons",
            "",
            "Each DT reaction consumes one tritium nucleus and emits one neutron.",
            f"The neutron energy in this ledger is `{_num(reaction.neutron_energy_kev)}` keV,",
            f"`{_num(reaction.neutron_energy_kev / 1000.0)}` MeV, from the two-body split of Q.",
            "Loaded tritium scales with fuel mass. Consumed tritium scales with",
            "the number of reactions.",
            "",
            "| Account | Fuel (kg) | Tritium loaded (kg) | Tritium consumed (kg) | Neutrons | Neutron energy (J) |",
            "| --- | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for account in accounts:
        lines.append(
            f"| {account.label} | {account.fuel_mass_kg:.6g} | "
            f"{account.tritium_loaded_kg:.6g} | {account.tritium_consumed_kg:.6g} | "
            f"{account.neutrons:.6g} | {account.neutron_energy_j:.6g} |"
        )
    lines.extend(
        [
            "",
            f"The p-B11 burn fraction at the checked-in peak, one walk at",
            f"`{_num(checked.peak_temperature_kev)}` keV, is `{pb11:.6g}`.",
            "That branch does not carry this tritium account.",
            "",
        ]
    )
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
