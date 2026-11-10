"""Trace the checked-in net through pressure, burn, conversion, and the store.

Net electrical energy is electric output minus coil supply. The four draws
stay the supplies locked in ``scripts/apply_fsot.py``. Plasma beta is thermal
pressure over magnetic pressure. The FSOT valve is recorded beside that ratio
and stays out of the joule product.

Two confined ledgers use the installed laws. Cooling divides the initial
temperature by the checked-in peak beta and calls ``evolve``. Cutting the
density by that same beta scales the undepleted schedule integral by beta
squared and the fuel ceiling by beta. Neither ledger is written back into
the checked-in YAML.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from solcascade.assemble import build_shell_from_slots
from solcascade.constants import dt_reaction
from solcascade.conversion import SupplyLedger, carnot_efficiency, convert
from solcascade.energy_baseline.phase2.map_fsot import (
    PIN_PREFIX,
    capture_scalars,
    load_engine,
)
from solcascade.geometry.coils import Coil, CoilSet
from solcascade.plasma.reactivity import dt_reactivity_m3_s
from solcascade.production.refine import (
    SUPPLIES,
    Problem,
    PulseYield,
    electric_per_fusion_j,
    fuel_of,
    load_checked_in,
    project_root,
    schedule_yield,
    temperature_scan,
)

# Geometry of scripts/apply_fsot.py helmholtz_config. The offset is the
# clearance past the shell outer radius used for that locked draw.
HELMHOLTZ_OFFSET_M = 0.021
HELMHOLTZ_TURNS = 12.0
HELMHOLTZ_WIRE_M = 0.008


@dataclass(frozen=True)
class ShellHeat:
    capacity_j_per_k: float
    checked_neutron_j: float
    checked_rise_k: float
    ceiling_neutron_j: float
    ceiling_rise_k: float


@dataclass(frozen=True)
class FieldStore:
    confining_field_t: float
    baseline_first_knot_j: float
    baseline_peak_j: float
    helmholtz_stored_j: float
    helmholtz_ledger_j: float
    recovery: float
    practical_j: float
    unrecovered_store_j: float


@dataclass(frozen=True)
class LimitedRow:
    """One temperature held fixed, density lowered until peak beta is 1."""

    peak_temperature_kev: float
    raw_beta: float
    density_scale: float
    fusion_j: float
    burn_fraction: float
    electric_j: float
    walk_time_s: float
    walk_yield_j: float
    efold_s: float
    ceiling_j: float
    complete_electric_j: float
    cap_complete_electric_j: float


@dataclass(frozen=True)
class FsotRecord:
    pin_prefix: str
    digest_prefix: str
    valve_beta: float
    scalars: dict[str, float]


@dataclass(frozen=True)
class FailureTrace:
    operating_per_j: float
    cap_per_j: float
    alpha_fraction: float
    neutron_fraction: float
    shell_channel: str
    shell: ShellHeat
    store: FieldStore
    checked: PulseYield
    cooled: PulseYield
    cooled_electric_j: float
    cooled_efold_s: float
    limited: tuple[LimitedRow, ...]
    witness: FsotRecord
    baseline_supply_j: float
    schedule_duration_s: float


def default_output_path():
    return project_root() / "runs" / "production" / "failure_trace.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _shell_capacity_j_per_k(problem: Problem) -> float:
    shell = build_shell_from_slots(problem.config, problem.slots)
    capacity = 0.0
    for cell in shell.cells:
        density = float(problem.slots[f"materials.{cell.material}.density_kg_m3"])
        heat = float(problem.slots[f"materials.{cell.material}.specific_heat_J_kgK"])
        capacity += density * heat * cell.volume_m3
    if capacity <= 0.0:
        raise ValueError("shell heat capacity must be positive")
    return capacity


def _efold_s(density_m3: float, temperature_kev: float, deuterium_fraction: float) -> float:
    """Time for one burn of the partner ions if the peak state is held."""
    partner = (1.0 - deuterium_fraction) * density_m3
    rate = partner * dt_reactivity_m3_s(temperature_kev)
    if rate <= 0.0:
        raise ValueError("peak burn rate must be positive")
    return 1.0 / rate


def _minimum_energy_j(coils: CoilSet, field_t: float) -> float:
    """Minimum of (1/2) I^T L I at one on-axis field. Same station as apply_fsot."""
    inductance = coils.inductance_matrix_h()
    count = len(coils.coils)
    weights = np.zeros(count, dtype=float)
    for index in range(count):
        probe = np.zeros(count, dtype=float)
        probe[index] = 1.0
        weights[index] = coils.on_axis_field(0.0, probe)
    solved = np.linalg.solve(inductance, weights)
    denominator = float(weights @ solved)
    if denominator == 0.0:
        raise ValueError("coils produce no on-axis field")
    currents = (field_t / denominator) * solved
    return 0.5 * float(currents @ inductance @ currents)


def _field_store(problem: Problem) -> FieldStore:
    coils = problem.coils
    base = coils.base_currents_a()
    scales = problem.schedule.current_scales(base.size)
    first = base * scales[0]
    peak_scale = float(np.max(problem.schedule.coil_current_scale))
    field = float(coils.on_axis_field(0.0, first))
    shell = build_shell_from_slots(problem.config, problem.slots)
    radius = shell.outer_radius_m + HELMHOLTZ_OFFSET_M
    spacing = 0.5 * radius
    pair = CoilSet(
        (
            Coil(
                "helmholtz_minus",
                radius,
                -spacing,
                HELMHOLTZ_TURNS,
                1.0,
                HELMHOLTZ_WIRE_M,
                "copper",
            ),
            Coil(
                "helmholtz_plus",
                radius,
                spacing,
                HELMHOLTZ_TURNS,
                1.0,
                HELMHOLTZ_WIRE_M,
                "copper",
            ),
        )
    )
    stored = _minimum_energy_j(pair, field)
    recovery = float(problem.slots["conversion.coil_energy_recovery_fraction"])
    ledger = float(dict(SUPPLIES)["same_b_helmholtz"])
    return FieldStore(
        confining_field_t=field,
        baseline_first_knot_j=float(coils.magnetic_energy_j(first)),
        baseline_peak_j=float(coils.magnetic_energy_j(base * peak_scale)),
        helmholtz_stored_j=stored,
        helmholtz_ledger_j=ledger,
        recovery=recovery,
        practical_j=ledger + (1.0 - recovery) * stored,
        unrecovered_store_j=(1.0 - recovery) * stored,
    )


def _conversion_cap_per_j(problem: Problem) -> float:
    """``convert`` with every alpha collected and every neutron joule through Carnot.

    The hot temperature is the configured shell temperature. This cap is a
    bound of the installed function. The checked-in ledger keeps
    ``shell_channel`` none.
    """
    reaction = dt_reaction()
    hot = float(problem.slots["heat.initial_temperature_K"])
    cold = float(problem.slots["conversion.cold_temperature_K"])
    supply = SupplyLedger(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    ledger = convert(
        fusion_j=1.0,
        alpha_fraction=float(reaction.alpha_fraction),
        neutron_deposited_fraction=1.0,
        alpha_escape_fraction=1.0,
        alpha_wall_fraction=0.0,
        alpha_direct_efficiency=1.0,
        shell_channel="carnot_fraction",
        carnot_fraction=1.0,
        hot_k=hot,
        cold_k=cold,
        supply=supply,
    )
    expected = float(reaction.alpha_fraction) + float(reaction.neutron_fraction) * carnot_efficiency(
        hot, cold
    )
    if not math.isclose(ledger.electric_gross_j, expected, rel_tol=1.0e-12, abs_tol=0.0):
        raise ValueError("converter cap drifted from the alpha plus Carnot split")
    return float(ledger.electric_gross_j)


def _limit_row(
    row: PulseYield,
    fuel_ceiling_j: float,
    operating_per_j: float,
    cap_per_j: float,
    deuterium_fraction: float,
) -> LimitedRow:
    if row.peak_beta <= 0.0:
        raise ValueError("peak beta must be positive")
    scale = 1.0 / row.peak_beta
    ceiling = fuel_ceiling_j * scale
    fusion = row.fusion_j * scale * scale
    return LimitedRow(
        peak_temperature_kev=row.peak_temperature_kev,
        raw_beta=row.peak_beta,
        density_scale=scale,
        fusion_j=fusion,
        burn_fraction=row.burn_fraction * scale,
        electric_j=fusion * operating_per_j,
        walk_time_s=row.walk_time_s,
        walk_yield_j=row.walk_yield_j * scale * scale,
        efold_s=_efold_s(row.peak_density_m3 * scale, row.peak_temperature_kev, deuterium_fraction),
        ceiling_j=ceiling,
        complete_electric_j=ceiling * operating_per_j,
        cap_complete_electric_j=ceiling * cap_per_j,
    )


def _witness() -> FsotRecord:
    module, digest = load_engine()
    if not str(digest).startswith(PIN_PREFIX):
        raise RuntimeError(f"engine pin is {str(digest)[:16]}")
    scalars = capture_scalars(module)
    valve = float(module.BETA)
    if valve <= 0.0 or valve >= 1.0:
        raise RuntimeError("engine valve is outside the recorded range")
    return FsotRecord(
        pin_prefix=PIN_PREFIX,
        digest_prefix=str(digest)[:16],
        valve_beta=valve,
        scalars=scalars,
    )


def _side(left: float, right: float) -> str:
    if left > right:
        return "above"
    if left < right:
        return "below"
    return "level with"


def build(problem: Problem | None = None) -> FailureTrace:
    problem = load_checked_in() if problem is None else problem
    fuel = fuel_of(problem)
    operating = electric_per_fusion_j(problem)
    cap = _conversion_cap_per_j(problem)
    reaction = dt_reaction()
    fraction = float(problem.slots["plasma.deuterium_fraction"])
    scan = temperature_scan(problem, fuel)
    checked = scan[0]
    if checked.peak_beta <= 1.0:
        raise ValueError("the checked-in peak is already at or under beta 1")
    cooled_t0 = checked.temperature0_kev / checked.peak_beta
    cooled = schedule_yield(problem, cooled_t0, fuel)
    shell_capacity = _shell_capacity_j_per_k(problem)
    checked_neutrons = checked.fusion_j * float(reaction.neutron_fraction)
    ceiling_neutrons = fuel.ceiling_j * float(reaction.neutron_fraction)
    limited = tuple(
        _limit_row(row, fuel.ceiling_j, operating, cap, fraction) for row in scan
    )
    return FailureTrace(
        operating_per_j=operating,
        cap_per_j=cap,
        alpha_fraction=float(reaction.alpha_fraction),
        neutron_fraction=float(reaction.neutron_fraction),
        shell_channel=str(problem.config["conversion"]["shell_channel"]),
        shell=ShellHeat(
            capacity_j_per_k=shell_capacity,
            checked_neutron_j=checked_neutrons,
            checked_rise_k=checked_neutrons / shell_capacity,
            ceiling_neutron_j=ceiling_neutrons,
            ceiling_rise_k=ceiling_neutrons / shell_capacity,
        ),
        store=_field_store(problem),
        checked=checked,
        cooled=cooled,
        cooled_electric_j=cooled.fusion_j * operating,
        cooled_efold_s=_efold_s(cooled.peak_density_m3, cooled.peak_temperature_kev, fraction),
        limited=limited,
        witness=_witness(),
        baseline_supply_j=float(dict(SUPPLIES)["baseline"]),
        schedule_duration_s=float(problem.schedule.duration_s),
    )


def render(problem: Problem | None = None) -> str:
    trace = build(problem)
    store = trace.store
    checked_electric = trace.checked.fusion_j * trace.operating_per_j
    best = max(trace.limited, key=lambda row: row.fusion_j)
    widest = max(trace.limited, key=lambda row: row.complete_electric_j)
    burn_factor = widest.complete_electric_j / widest.electric_j
    store_factor = store.practical_j / widest.complete_electric_j
    converter_factor = trace.cap_per_j / trace.operating_per_j
    original_complete = trace.checked.fusion_j / trace.checked.burn_fraction * trace.operating_per_j
    hold_factor = widest.efold_s / widest.walk_time_s
    lines = [
        "# Failure trace",
        "",
        "Net electrical energy is electric output minus coil supply. The chain",
        "below is that difference on the checked-in pulse: pressure balance,",
        "the burn integral, the converter, and the stored field. The warp-shell",
        "rest energy is outside this score.",
        "",
        "## Checked-in product",
        "",
        f"The operating converter returns `{_num(trace.operating_per_j)}` J of electricity",
        "per fusion joule. That product is the alpha share of Q times the YAML",
        "escape fraction times the YAML direct efficiency.",
        f"The alpha share is `{_num(trace.alpha_fraction)}`. The neutron share is",
        f"`{_num(trace.neutron_fraction)}`. `shell_channel` is `{trace.shell_channel}`,",
        "so the neutron share is heat on this ledger.",
        "",
        f"Checked-in fusion is `{_num(trace.checked.fusion_j)}` J. Electricity is",
        f"`{_num(checked_electric)}` J. The baseline supply is `{_num(trace.baseline_supply_j)}` J,",
        f"and the net is `{_num(checked_electric - trace.baseline_supply_j)}` J.",
        f"The Helmholtz ledger draw is `{_num(store.helmholtz_ledger_j)}` J, and the",
        f"net against that draw is `{_num(checked_electric - store.helmholtz_ledger_j)}` J.",
        "",
        "## Pressure",
        "",
        f"Peak beta on the checked-in schedule is `{_num(trace.checked.peak_beta)}`.",
        "Beta is thermal pressure divided by magnetic pressure. Above 1, the",
        "thermal pressure exceeds the frozen-flux field that the schedule prescribes.",
        "",
        f"Cooling at fixed density divides the initial temperature by that beta.",
        f"The evolved initial temperature is `{_num(trace.cooled.temperature0_kev)}` keV",
        f"and the evolved peak beta is `{_num(trace.cooled.peak_beta)}`. The evolved",
        f"schedule releases `{_num(trace.cooled.fusion_j)}` J, electric",
        f"`{_num(trace.cooled_electric_j)}` J. The partner-ion e-folding time at that",
        f"peak is `{_num(trace.cooled_efold_s)}` s. The gyro-walk there is",
        f"`{_num(trace.cooled.walk_time_s)}` s.",
        "",
        "The density-limited rows keep each temperature trajectory and divide",
        "the density by the raw peak beta. Fusion then scales as the square of",
        "that factor because the integrator is undepleted, and the ceiling",
        "scales as the density. Peak beta on each of those rows is 1.",
        "",
        "| Peak (keV) | Raw beta | Density scale | Fusion (J) | Electric (J) | Walk yield (J) | E-fold (s) | Complete electric (J) |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in trace.limited:
        lines.append(
            f"| {row.peak_temperature_kev:.6g} | {row.raw_beta:.6g} | "
            f"{row.density_scale:.6g} | {row.fusion_j:.6g} | {row.electric_j:.6g} | "
            f"{row.walk_yield_j:.6g} | {row.efold_s:.6g} | {row.complete_electric_j:.6g} |"
        )
    lines.extend(
        [
            "",
            f"The strongest density-limited schedule releases `{_num(best.fusion_j)}` J",
            f"at `{_num(best.peak_temperature_kev)}` keV, electric `{_num(best.electric_j)}` J.",
            f"Its walk yield is `{_num(best.walk_yield_j)}` J. Its e-folding time is",
            f"`{_num(best.efold_s)}` s beside a walk of `{_num(best.walk_time_s)}` s and a",
            f"schedule duration of `{_num(trace.schedule_duration_s)}` s.",
            "The schedule burns a fraction of the ceiling.",
            "",
            f"Complete burn of that reduced density would deliver `{_num(best.complete_electric_j)}` J",
            "at the checked-in converter. The schedule on the same row is short of",
            f"that ceiling by `{_num(best.complete_electric_j / best.electric_j)}`.",
            f"The Helmholtz practical supply is `{_num(store.practical_j)}` J, and that",
            f"ceiling is short of that supply by `{_num(store.practical_j / best.complete_electric_j)}`.",
            "",
            "## Shell heat",
            "",
            f"The built shell stores `{_num(trace.shell.capacity_j_per_k)}` J per kelvin.",
            f"Checked-in neutron energy is `{_num(trace.shell.checked_neutron_j)}` J and",
            f"would raise that shell by `{_num(trace.shell.checked_rise_k)}` K if it all",
            "stopped in the shell. Complete burn of the original sphere would raise",
            f"it by `{_num(trace.shell.ceiling_rise_k)}` K. The configured shell temperature",
            "is the heat boundary the converter reads. The rises above are the pulse.",
            "",
            "## Converter bound",
            "",
            f"The installed `convert` function, with every alpha collected and every",
            "neutron joule passed through the Carnot factor of the configured shell",
            f"and the configured cold sink, returns `{_num(trace.cap_per_j)}` J per fusion",
            f"joule. The operating product is `{_num(trace.operating_per_j)}` J. The ratio",
            f"of those two products is `{_num(converter_factor)}`.",
            "",
            f"On the strongest density-limited schedule the bound times the integral is",
            f"`{_num(best.fusion_j * trace.cap_per_j)}` J. On the widest density-limited",
            f"ceiling the bound times complete burn is `{_num(widest.cap_complete_electric_j)}` J.",
            "Widest converter-cap complete burn versus the practical supply: "
            f"{_side(widest.cap_complete_electric_j, store.practical_j)}.",
            "The checked-in ledger keeps the operating product. The bound is the",
            "parameter cap of `convert()`. The fusion column is the schedule integral.",
            "",
            "## Stored field",
            "",
            f"The confining on-axis field at the first current knot is `{_num(store.confining_field_t)}` T.",
            f"The four-coil store is `{_num(store.baseline_first_knot_j)}` J at that knot",
            f"and `{_num(store.baseline_peak_j)}` J at the peak scale. The baseline supply",
            "on the locked ledger is the draw of that ramp.",
            "",
            f"The Helmholtz pair at the same on-axis field stores `{_num(store.helmholtz_stored_j)}` J.",
            f"Its locked ledger draw is `{_num(store.helmholtz_ledger_j)}` J. At recovery",
            f"`{_num(store.recovery)}` the practical supply is the ledger plus the",
            f"unrecovered store, `{_num(store.practical_j)}` J. A zero-resistivity winding",
            "sets the ohmic ledger to zero and leaves the unrecovered store,",
            f"`{_num(store.unrecovered_store_j)}` J.",
            "",
            f"Original-density complete burn at the checked-in converter is `{_num(original_complete)}` J.",
            f"Original-density complete burn versus the Helmholtz ledger: {_side(original_complete, store.helmholtz_ledger_j)}.",
            f"Original-density complete burn versus the unrecovered store: {_side(original_complete, store.unrecovered_store_j)}.",
            f"The widest density-limited ceiling is the `{_num(widest.peak_temperature_kev)}` keV row,",
            f"`{_num(widest.complete_electric_j)}` J at the checked-in converter. Its schedule",
            f"electricity is `{_num(widest.electric_j)}` J and its e-folding time is",
            f"`{_num(widest.efold_s)}` s.",
            f"Widest complete burn versus the Helmholtz ledger: {_side(widest.complete_electric_j, store.helmholtz_ledger_j)}.",
            f"Widest complete burn versus the practical supply: {_side(widest.complete_electric_j, store.practical_j)}.",
            f"Widest complete burn versus the unrecovered store: {_side(widest.complete_electric_j, store.unrecovered_store_j)}.",
            "",
            "## Field scaling",
            "",
            "At beta 1 the density scales with the square of B. Undepleted schedule",
            "fusion scales with the square of density, so with B to the fourth,",
            "while magnetic energy and ohmic energy on a fixed winding and a fixed",
            "pulse length scale with B squared. The undepleted ratio therefore",
            "improves with B squared while the burn fraction stays below 1.",
            "At complete burn the ceiling scales with density, so with B squared,",
            "and electric energy over the store stays at the ratio for that temperature.",
            "The highest complete-burn ratio in this table is the "
            f"`{_num(widest.peak_temperature_kev)}` keV row, `{_num(widest.complete_electric_j / store.practical_j)}`.",
            "Hotter rows sit lower because beta is higher and the beta = 1 ceiling",
            "is smaller. Raising the field while holding density and temperature",
            "multiplies the store by the present beta and leaves the fusion integral",
            "where it is.",
            "",
            "## FSOT record",
            "",
            f"Pin `{trace.witness.pin_prefix}` matches the loaded engine prefix",
            f"`{trace.witness.digest_prefix}`. The scalars below are dimensionless.",
            "The FSOT scalars are recorded and stay out of the joule product.",
            "",
            f"Engine valve: `{trace.witness.valve_beta:.6e}`.",
            f"Plasma beta on the checked-in peak: `{_num(trace.checked.peak_beta)}`.",
            "Engine valve and plasma beta are different quantities.",
            "",
            "| Scalar | Value |",
            "| --- | ---: |",
        ]
    )
    for name, value in trace.witness.scalars.items():
        lines.append(f"| {name} | {value:.6g} |")
    lines.extend(
        [
            "",
            "K is the scalar rest unit. gamma_FSOT is the orifice witness.",
            "phi^-4 is the filament-fraction ceiling recorded on the coil scan.",
            "The enthalpy identity in the FSOT physics notes uses a power that is",
            "the logarithm of one chemistry leaf. This trace assigns no power to",
            "the fusion yield.",
            "",
            "## Standing",
            "",
            "The confined schedule remains a net consumer on every beta = 1 row.",
            f"On the widest ceiling the burn factor is `{_num(burn_factor)}`, the store",
            f"factor is `{_num(store_factor)}`, and the converter factor is",
            f"`{_num(converter_factor)}`. The burn factor is complete-burn electricity",
            "divided by the schedule electricity. The store factor is the practical",
            "supply divided by that complete-burn electricity. The converter factor",
            "is the Carnot-and-alpha cap divided by the operating product. Holding",
            f"that peak for one e-folding takes `{_num(hold_factor)}` walks. The walk",
            "is the confinement time already installed. The FSOT record does not",
            "move reactivity, density, field, or the converter product.",
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
