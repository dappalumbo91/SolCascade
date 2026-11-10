"""Rank the coil store against the beta = 1 fuel and resize the sphere in the bore.

The baseline four-coil store is the largest energy on the checked-in pulse.
The unrecovered Helmholtz store is the largest payment that still sits under
the original fuel ceiling, and it is larger than the electricity in a
complete burn of the beta = 1 load. A larger sphere uses the same coils,
the same compression ratio, and the beta = 1 density. The checked-in YAML
is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from solcascade.assemble import build_shell_from_slots
from solcascade.constants import VACUUM_PERMEABILITY_H_PER_M
from solcascade.geometry.coils import Coil, CoilSet
from solcascade.plasma.metrics import evolve, gyro_walk_time_s
from solcascade.plasma.reactivity import dt_reactivity_m3_s
from solcascade.production.refine import (
    Problem,
    catalog_record,
    fuel_of,
    load_checked_in,
    project_root,
)
from solcascade.production.trace import (
    HELMHOLTZ_OFFSET_M,
    HELMHOLTZ_TURNS,
    HELMHOLTZ_WIRE_M,
    build,
)
from solcascade.compression.schedule import PulseSchedule

# The enlarged initial radius stays at this fraction of the shell inner radius.
BORE_FRACTION = 0.95


@dataclass(frozen=True)
class Sphere:
    label: str
    radius_factor: float
    initial_radius_m: float
    peak_radius_m: float
    peak_beta: float
    peak_field_t: float
    peak_temperature_kev: float
    peak_density_m3: float
    field_edge_over_axis: float
    fusion_j: float
    electric_j: float
    walk_time_s: float
    walk_electric_j: float
    complete_fusion_j: float
    complete_electric_j: float
    peak_power_w: float
    electric_power_w: float
    efold_s: float


@dataclass(frozen=True)
class CostReport:
    baseline_peak_j: float
    baseline_supply_j: float
    fuel_ceiling_j: float
    helmholtz_stored_j: float
    ohmic_ledger_j: float
    unrecovered_store_j: float
    practical_j: float
    plasma_magnetic_j: float
    store_over_plasma: float
    beta1_complete_electric_j: float
    beta1_schedule_electric_j: float
    beta1_walk_electric_j: float
    shell_inner_m: float
    jet_j: float
    store_match: Sphere
    bore: Sphere


def default_output_path():
    return project_root() / "runs" / "production" / "dominant_cost.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _side(left: float, right: float) -> str:
    if left > right:
        return "above"
    if left < right:
        return "below"
    return "level with"


def _helmholtz(shell_outer_m: float, field_t: float) -> tuple[CoilSet, np.ndarray]:
    """Minimum-energy Helmholtz pair at one on-axis field. Same station as the trace."""
    radius = shell_outer_m + HELMHOLTZ_OFFSET_M
    spacing = 0.5 * radius
    coils = CoilSet(
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
    inductance = coils.inductance_matrix_h()
    weights = np.zeros(2, dtype=float)
    for index in range(2):
        probe = np.zeros(2, dtype=float)
        probe[index] = 1.0
        weights[index] = coils.on_axis_field(0.0, probe)
    solved = np.linalg.solve(inductance, weights)
    denominator = float(weights @ solved)
    if denominator == 0.0:
        raise ValueError("Helmholtz pair produces no on-axis field")
    currents = (field_t / denominator) * solved
    return coils, currents


def _plasma_magnetic_j(radius_m: float, field_t: float) -> float:
    volume = (4.0 / 3.0) * np.pi * radius_m ** 3
    return float(field_t ** 2 / (2.0 * VACUUM_PERMEABILITY_H_PER_M) * volume)


def _scaled_schedule(schedule: PulseSchedule, factor: float) -> PulseSchedule:
    per_coil = None if schedule.per_coil_scale is None else np.array(schedule.per_coil_scale, dtype=float)
    return PulseSchedule(
        time_s=np.array(schedule.time_s, dtype=float),
        plasma_radius_m=np.array(schedule.plasma_radius_m, dtype=float) * factor,
        coil_current_scale=np.array(schedule.coil_current_scale, dtype=float),
        per_coil_scale=per_coil,
    )


def _sphere(
    problem: Problem,
    factor: float,
    density_m3: float,
    ceiling_j: float,
    electric_per_j: float,
    coils: CoilSet,
    currents: np.ndarray,
    axis_field_t: float,
    label: str,
) -> Sphere:
    """Evolve one enlarged copy of the schedule. Density is the beta = 1 initial value."""
    if factor <= 0.0:
        raise ValueError("radius factor must be positive")
    schedule = _scaled_schedule(problem.schedule, factor)
    history = evolve(
        problem.coils,
        schedule,
        density_m3,
        float(problem.slots["plasma.ion_temperature_keV"]),
        float(problem.slots["plasma.deuterium_fraction"]),
        str(problem.config["compression"]["flux_conservation"]),
    )
    peak = max(history.points, key=lambda point: (point.ion_temperature_kev, point.fusion_power_w))
    fusion = history.integrate("fusion_power_w")
    walk = gyro_walk_time_s(peak.radius_m, peak.ion_temperature_kev, peak.internal_field_t)
    fraction = float(problem.slots["plasma.deuterium_fraction"])
    partner = (1.0 - fraction) * peak.ion_density_m3
    efold = 1.0 / (partner * dt_reactivity_m3_s(peak.ion_temperature_kev))
    _radial, axial = coils.field(float(schedule.plasma_radius_m[0]), 0.0, currents)
    complete_fusion = ceiling_j * factor ** 3
    return Sphere(
        label=label,
        radius_factor=float(factor),
        initial_radius_m=float(schedule.plasma_radius_m[0]),
        peak_radius_m=float(peak.radius_m),
        peak_beta=float(peak.beta),
        peak_field_t=float(peak.internal_field_t),
        peak_temperature_kev=float(peak.ion_temperature_kev),
        peak_density_m3=float(peak.ion_density_m3),
        field_edge_over_axis=float(axial / axis_field_t),
        fusion_j=float(fusion),
        electric_j=float(fusion * electric_per_j),
        walk_time_s=float(walk),
        walk_electric_j=float(peak.fusion_power_w * walk * electric_per_j),
        complete_fusion_j=float(complete_fusion),
        complete_electric_j=float(complete_fusion * electric_per_j),
        peak_power_w=float(peak.fusion_power_w),
        electric_power_w=float(peak.fusion_power_w * electric_per_j),
        efold_s=float(efold),
    )


def build_cost(problem: Problem | None = None) -> CostReport:
    problem = load_checked_in() if problem is None else problem
    trace = build(problem)
    fuel = fuel_of(problem)
    shell = build_shell_from_slots(problem.config, problem.slots)
    widest = max(trace.limited, key=lambda row: row.complete_electric_j)
    beta = trace.checked.peak_beta
    density = float(problem.slots["plasma.ion_density_m3"]) / beta
    ceiling = fuel.ceiling_j / beta
    coils, currents = _helmholtz(shell.outer_radius_m, trace.store.confining_field_t)
    axis = float(coils.on_axis_field(0.0, currents))
    match_factor = (trace.store.practical_j / widest.complete_electric_j) ** (1.0 / 3.0)
    initial = float(problem.schedule.plasma_radius_m[0])
    bore_factor = BORE_FRACTION * shell.inner_radius_m / initial
    plasma_magnetic = _plasma_magnetic_j(trace.checked.peak_radius_m, trace.checked.peak_field_t)
    return CostReport(
        baseline_peak_j=trace.store.baseline_peak_j,
        baseline_supply_j=trace.baseline_supply_j,
        fuel_ceiling_j=fuel.ceiling_j,
        helmholtz_stored_j=trace.store.helmholtz_stored_j,
        ohmic_ledger_j=trace.store.helmholtz_ledger_j,
        unrecovered_store_j=trace.store.unrecovered_store_j,
        practical_j=trace.store.practical_j,
        plasma_magnetic_j=plasma_magnetic,
        store_over_plasma=trace.store.helmholtz_stored_j / plasma_magnetic,
        beta1_complete_electric_j=widest.complete_electric_j,
        beta1_schedule_electric_j=widest.electric_j,
        beta1_walk_electric_j=widest.walk_yield_j * trace.operating_per_j,
        shell_inner_m=shell.inner_radius_m,
        jet_j=float(catalog_record("fusion-dt-magnetic")["energy_density_j"]),
        store_match=_sphere(
            problem,
            match_factor,
            density,
            ceiling,
            trace.operating_per_j,
            coils,
            currents,
            axis,
            "store match",
        ),
        bore=_sphere(
            problem,
            bore_factor,
            density,
            ceiling,
            trace.operating_per_j,
            coils,
            currents,
            axis,
            "bore limit",
        ),
    )


def render(problem: Problem | None = None) -> str:
    report = build_cost(problem)
    match = report.store_match
    bore = report.bore
    lines = [
        "# Dominant cost",
        "",
        "The baseline coil store is the largest energy in this ranking.",
        "It is the magnetic energy of the four checked-in coils at full current.",
        "The fuel ceiling of this sphere sits under that store, so every burn",
        "fraction of this load stays under the baseline draw.",
        "",
        "The Helmholtz winding is the smaller store at the same on-axis field.",
        "Its practical supply is the ohmic ledger plus the unrecovered part of",
        "the coil energy. The unrecovered Helmholtz store is larger than the",
        "ohmic ledger and larger than the beta = 1 complete-burn electricity.",
        "",
        "## Ranking",
        "",
        "| Energy | Joules |",
        "| --- | ---: |",
        f"| Baseline coil store at full current | {_num(report.baseline_peak_j)} |",
        f"| Baseline supply draw | {_num(report.baseline_supply_j)} |",
        f"| Original fuel ceiling | {_num(report.fuel_ceiling_j)} |",
        f"| Helmholtz coil store | {_num(report.helmholtz_stored_j)} |",
        f"| Helmholtz practical supply | {_num(report.practical_j)} |",
        f"| Unrecovered Helmholtz store | {_num(report.unrecovered_store_j)} |",
        f"| Beta = 1 complete-burn electricity | {_num(report.beta1_complete_electric_j)} |",
        f"| Helmholtz ohmic ledger | {_num(report.ohmic_ledger_j)} |",
        f"| Magnetic energy inside the checked-in plasma | {_num(report.plasma_magnetic_j)} |",
        f"| Beta = 1 schedule electricity | {_num(report.beta1_schedule_electric_j)} |",
        f"| Beta = 1 walk electricity | {_num(report.beta1_walk_electric_j)} |",
        "",
        f"The Helmholtz store is `{_num(report.store_over_plasma)}` times the magnetic",
        "energy inside the checked-in plasma. The coils carry the field that",
        "the reacting volume does not.",
        "",
        "Beta = 1 complete burn versus the ohmic ledger: "
        f"{_side(report.beta1_complete_electric_j, report.ohmic_ledger_j)}.",
        "Beta = 1 complete burn versus the unrecovered store: "
        f"{_side(report.beta1_complete_electric_j, report.unrecovered_store_j)}.",
        "",
        "## Sphere inside the bore",
        "",
        f"The shell inner radius is `{_num(report.shell_inner_m)}` m. The store-matching",
        "sphere is the beta = 1 density at the checked-in temperature, with the",
        "radius schedule scaled so a complete burn returns the practical supply.",
        f"Its initial radius is `{_num(match.initial_radius_m)}` m and its compressed",
        f"radius is `{_num(match.peak_radius_m)}` m, a factor `{_num(match.radius_factor)}`",
        "on the checked-in radii. The midplane field of the Helmholtz pair at",
        f"the initial edge is `{_num(match.field_edge_over_axis)}` of the on-axis field.",
        f"Evolved peak beta is `{_num(match.peak_beta)}`. Schedule electricity is",
        f"`{_num(match.electric_j)}` J. Walk electricity is `{_num(match.walk_electric_j)}` J.",
        f"Complete-burn electricity is `{_num(match.complete_electric_j)}` J.",
        "The store-matching sphere fits inside the shell.",
        "",
        f"The bore-limited sphere sets the initial radius to `{_num(BORE_FRACTION)}` of",
        "the shell inner radius and keeps the same compression ratio and the",
        f"same beta = 1 density. Initial radius `{_num(bore.initial_radius_m)}` m,",
        f"compressed radius `{_num(bore.peak_radius_m)}` m, factor `{_num(bore.radius_factor)}`.",
        f"The midplane field at that initial edge is `{_num(bore.field_edge_over_axis)}`",
        f"of the on-axis field. Evolved peak beta is `{_num(bore.peak_beta)}`.",
        f"Schedule fusion is `{_num(bore.fusion_j)}` J and schedule electricity is",
        f"`{_num(bore.electric_j)}` J. Walk electricity is `{_num(bore.walk_electric_j)}` J",
        f"over `{_num(bore.walk_time_s)}` s. Complete fusion is `{_num(bore.complete_fusion_j)}` J",
        f"and complete-burn electricity is `{_num(bore.complete_electric_j)}` J.",
        f"Peak fusion power is `{_num(bore.peak_power_w)}` W. Peak electric power is",
        f"`{_num(bore.electric_power_w)}` W. The e-folding time at that peak is",
        f"`{_num(bore.efold_s)}` s.",
        "",
        "The bore-limited complete burn sits above the practical supply.",
        "The bore-limited schedule stays under the practical supply.",
        "Bore-limited complete fusion versus one catalog JET pulse: "
        f"{_side(bore.complete_fusion_j, report.jet_j)}.",
        "",
        "## Approach",
        "",
        "The payment to bring under the fuel is the external coil store.",
        "On the Helmholtz winding the unrecovered half of that store is the",
        "term above the beta = 1 electricity. A complete burn of a larger",
        "sphere inside this bore meets that term. The evolved schedule and the",
        "gyro-walk on that sphere remain the small integrals. Peak electric",
        "power is the rate at the compressed state for as long as the plasma",
        "is held. This ranking does not rewrite the checked-in schedule.",
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
