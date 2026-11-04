"""Minimal working pipeline.

Coil field, plasma pulse, shell deposition, transient heat, stress, and the
electrical ledger. OpenMC, MOOSE, Cardinal, TokaMaker, Griffin, and Bison
are written as inputs or tables and executed only when the local install
can run them.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Mapping

import numpy as np

from solcascade.assemble import build_coils, build_pulse, build_shell_from_slots, validate_geometry
from solcascade.constants import VACUUM_PERMEABILITY_H_PER_M, dt_reaction
from solcascade.conversion import convert, supply_energy
from solcascade.coupling.cardinal import write_skeleton
from solcascade.coupling.moose import run_moose, write_heat_input
from solcascade.coupling.neams import write_tables
from solcascade.coupling.tokamaker import tokamaker_status, write_script
from solcascade.geometry.coils import CoilSet
from solcascade.hybrid.blanket import fission_power_w, run_hybrid, scale_fertile_material
from solcascade.interfaces import (
    MappingProvider,
    NullProvider,
    ScalarProvider,
    apply_provider,
    collect_slots,
    overrides_only,
    reconcile_energy_density,
)
from solcascade.load import load_problem
from solcascade.materials.library import Material, load_materials
from solcascade.mechanics.heat import solve_transient
from solcascade.mechanics.stress import evaluate_stress
from solcascade.plasma.equilibrium import equilibrium_summary
from solcascade.plasma.metrics import evolve
from solcascade.transport.analytical import first_collision
from solcascade.transport.openmc_backend import (
    neutron_source_energy_ev,
    openmc_status,
    run_openmc,
    write_model_script,
)


def run_pipeline(
    config_path: Path,
    workdir: Path | None = None,
    provider: ScalarProvider | None = None,
    provider_overrides: Mapping[str, float] | None = None,
) -> dict:
    config, raw_materials = load_problem(Path(config_path))
    slots = collect_slots(config, raw_materials)
    active = _provider(config, provider, provider_overrides)
    resolved, provenance = apply_provider(slots, active)
    resolved, provenance = reconcile_energy_density(resolved, provenance)
    library = scale_fertile_material(config, resolved, load_materials(raw_materials, resolved))
    coils = build_coils(config, resolved)
    shell = build_shell_from_slots(config, resolved)
    schedule = build_pulse(config, resolved, provenance)
    validate_geometry(coils, shell, schedule)
    history = evolve(
        coils,
        schedule,
        ion_density_m3=resolved["plasma.ion_density_m3"],
        ion_temperature_kev=resolved["plasma.ion_temperature_keV"],
        deuterium_fraction=resolved["plasma.deuterium_fraction"],
        flux_conservation=str(config["compression"]["flux_conservation"]),
        electron_temperature_equals_ion=bool(
            config["plasma"].get("electron_temperature_equals_ion", True)
        ),
    )
    deposition = first_collision(shell, library)
    reaction = dt_reaction()
    duration = schedule.duration_s
    fusion_j = history.integrate("fusion_power_w")
    neutron_j = history.integrate("neutron_power_w")
    alpha_j = history.integrate("alpha_power_w")
    destination = Path(workdir) if workdir else Path(config["output"]["workdir"])
    # OpenMC execution changes the working directory. Later reads of the
    # generated model must not depend on that temporary directory.
    destination = destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    neutron_power, neutron_fraction, openmc_info, transport_backend = _neutron_power(
        destination,
        config,
        shell,
        library,
        history,
        deposition,
        neutron_j,
        duration,
        resolved["heat.initial_temperature_K"],
    )
    hybrid = run_hybrid(
        config, resolved, shell, library, deposition, neutron_j, openmc_info, destination
    )
    alpha_wall_power = np.zeros(len(shell.cells), dtype=float)
    if duration > 0.0 and shell.cells:
        alpha_wall_power[0] = (
            alpha_j * resolved["conversion.alpha_wall_fraction"] / duration
        )
    source_w = neutron_power + alpha_wall_power
    if not (isinstance(openmc_info, dict) and openmc_info.get("heating_applied")):
        source_w = source_w + fission_power_w(
            shell, deposition, hybrid["fission_thermal_J"], duration
        )
    heat = solve_transient(
        shell,
        library,
        source_w,
        duration_s=duration,
        steps=int(config["heat"]["time_steps"]),
        initial_temperature_k=resolved["heat.initial_temperature_K"],
        sink_temperature_k=resolved["heat.sink_temperature_K"],
        heat_transfer_coefficient=resolved["heat.outer_heat_transfer_coefficient_W_m2K"],
    )
    pressure_pa, pressure_source = _shell_pressure(coils, shell, history)
    stress = evaluate_stress(
        shell,
        library,
        heat.temperature_k,
        stress_free_temperature_k=resolved["mechanics.stress_free_temperature_K"],
        internal_pressure_pa=pressure_pa,
        pressure_source=pressure_source,
    )
    supply = _supply(coils, library, schedule, resolved)
    ledger = convert(
        fusion_j=fusion_j,
        alpha_fraction=reaction.alpha_fraction,
        neutron_deposited_fraction=neutron_fraction,
        alpha_escape_fraction=resolved["conversion.alpha_escape_fraction"],
        alpha_wall_fraction=resolved["conversion.alpha_wall_fraction"],
        alpha_direct_efficiency=resolved["conversion.alpha_direct_efficiency"],
        shell_channel=str(config["conversion"]["shell_channel"]),
        carnot_fraction=resolved["conversion.carnot_fraction"],
        hot_k=float(np.max(heat.temperature_k)),
        cold_k=resolved["conversion.cold_temperature_K"],
        supply=supply,
        fission_j=hybrid["fission_thermal_J"],
        thermal_efficiency=hybrid["thermal_efficiency"],
    )
    _write_history(destination / "plasma_history.csv", history)
    _write_radial(destination / "radial_map.csv", shell, deposition, source_w, heat, stress)
    _write_field(destination, config, coils, schedule)
    moose_info = _write_moose(
        destination, config, resolved, shell, library, heat.source_w_m3, duration
    )
    write_skeleton(
        destination / "cardinal_skeleton.i",
        particles=int(config["transport"]["openmc"]["particles"]),
        batches=int(config["transport"]["openmc"]["batches"]),
        source_rate_per_s=_average_rate(history),
        heat_input_name="heat_conduction.i",
    )
    write_script(
        destination / "tokamaker_equilibrium.py",
        coils,
        config["equilibrium"].get("tokamaker_mesh"),
    )
    mises = _mises_per_cell(shell, stress)
    neams = write_tables(
        destination, shell, heat.source_w_m3, heat.temperature_k, mises
    )
    equilibrium = equilibrium_summary(history, str(config["equilibrium"]["solver"]))
    equilibrium["tokamaker"] = tokamaker_status()
    summary = {
        "name": config["name"],
        "description": str(config.get("description", "")).strip(),
        "provider": active.name(),
        "overrides": overrides_only(provenance),
        "reaction": {
            "q_keV": reaction.q_kev,
            "neutron_energy_keV": reaction.neutron_energy_kev,
            "alpha_energy_keV": reaction.alpha_energy_kev,
            "neutron_fraction": reaction.neutron_fraction,
            "alpha_fraction": reaction.alpha_fraction,
        },
        "plasma": {
            "flux_conservation": history.flux_conservation,
            "adiabatic_index": history.adiabatic_index,
            "duration_s": duration,
            "peak_fusion_power_W": float(np.max(history.column("fusion_power_w"))),
            "peak_ion_temperature_keV": float(np.max(history.column("ion_temperature_kev"))),
            "peak_ion_density_m3": float(np.max(history.column("ion_density_m3"))),
            "peak_internal_field_T": float(np.max(history.column("internal_field_t"))),
            "peak_vacuum_field_T": float(np.max(history.column("vacuum_field_t"))),
            "peak_triple_product_ratio": float(np.max(history.column("triple_product_ratio"))),
            "fusion_energy_J": fusion_j,
            "neutron_energy_J": neutron_j,
            "alpha_energy_J": alpha_j,
        },
        "transport": {
            "backend": transport_backend,
            "approximation": deposition.approximation,
            "deposited_fraction": neutron_fraction,
            "leakage_uncollided_fraction": deposition.leakage_uncollided_fraction,
            "unresolved_fraction": deposition.unresolved_fraction,
            "openmc": openmc_info,
        },
        "heat": {
            "reported_by": heat.mode,
            "peak_temperature_K": float(np.max(heat.temperature_k)),
            "sink_temperature_K": resolved["heat.sink_temperature_K"],
            "moose": moose_info,
        },
        "stress": {
            "internal_pressure_Pa": stress.internal_pressure_pa,
            "pressure_source": stress.pressure_source,
            "layers": [
                {
                    "layer": layer.layer,
                    "structural": layer.structural,
                    "max_von_mises_Pa": None
                    if not layer.structural
                    else float(np.max(layer.von_mises_pa)),
                    "yield_strength_Pa": layer.yield_strength_pa,
                    "exceeds_yield": layer.exceeds_yield,
                }
                for layer in stress.layers
            ],
        },
        "equilibrium": equilibrium,
        "neams": {"executed": False, "tables": neams},
        "hybrid": hybrid,
        "conversion": _ledger_dict(ledger),
    }
    (destination / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    (destination / "report.md").write_text(_report(summary), encoding="utf-8")
    summary["workdir"] = str(destination)
    return summary


def _provider(config, provider, overrides):
    if provider is not None and overrides:
        raise ValueError("pass either a provider or override values, not both")
    if provider is not None:
        return provider
    values = dict(config.get("provider", {}).get("overrides") or {})
    if overrides:
        values.update(overrides)
    if not values:
        return NullProvider()
    return MappingProvider(values)


def _shell_pressure(coils: CoilSet, shell, history):
    """Conservative inner-wall pressure: the larger of plasma and magnetic pressure."""
    peak_plasma = max(point.pressure_pa for point in history.points)
    strongest = int(np.argmax(history.column("vacuum_field_t")))
    base_field = coils.on_axis_field(0.0, coils.base_currents_a())
    scale = 0.0 if base_field == 0.0 else history.points[strongest].vacuum_field_t / base_field
    currents = coils.base_currents_a() * scale
    _radial, axial = coils.field(shell.inner_radius_m, 0.0, currents)
    magnetic = axial * axial / (2.0 * VACUUM_PERMEABILITY_H_PER_M)
    if magnetic >= peak_plasma:
        return float(magnetic), "vacuum_magnetic_pressure_at_inner_wall"
    return float(peak_plasma), "peak_plasma_pressure"


def _supply(coils: CoilSet, library: dict[str, Material], schedule, resolved):
    temperature = resolved["conversion.conductor_temperature_K"]
    resistivity = []
    for coil in coils.coils:
        material = library[coil.material]
        resistivity.append(material.resistivity_at(temperature))
    resistivity_array = np.array(resistivity, dtype=float)
    inductance = coils.inductance_matrix_h()
    base = coils.base_currents_a()
    magnetic = []
    ohmic = []
    for scale in schedule.current_scales(base.size):
        currents = base * scale
        magnetic.append(coils.magnetic_energy_j(currents, inductance))
        ohmic.append(coils.ohmic_power_w(currents, resistivity_array))
    return supply_energy(
        schedule.time_s,
        np.array(magnetic),
        np.array(ohmic),
        resolved["conversion.coil_energy_recovery_fraction"],
    )


def _neutron_power(
    destination, config, shell, library, history, deposition, neutron_j, duration, temperature_k
):
    """Analytical cell power, replaced by OpenMC heating when a tally is read."""
    from solcascade.constants import JOULE_PER_EV

    analytical = deposition.power_w(neutron_j / duration if duration else 0.0)
    info = _write_openmc(
        destination, config, shell, library, history, neutron_j, duration, temperature_k
    )
    heating = (info.get("heating") or {}) if isinstance(info, dict) else {}
    values = heating.get("heating_eV_per_source")
    if not values:
        return analytical, float(np.sum(deposition.energy_fraction)), info, deposition.backend
    rate = _average_rate(history)
    power = np.zeros(len(shell.cells), dtype=float)
    if len(values) != len(shell.layers):
        info = dict(info)
        info["heating_applied"] = False
        info["heating_reason"] = "tally length does not match shell layers"
        return analytical, float(np.sum(deposition.energy_fraction)), info, deposition.backend
    for layer, heating_ev in zip(shell.layers, values, strict=True):
        indexes = [cell.index for cell in shell.cells if cell.layer == layer.name]
        layer_power = max(float(heating_ev), 0.0) * JOULE_PER_EV * rate
        weights = np.array([deposition.energy_fraction[index] for index in indexes], dtype=float)
        if float(weights.sum()) <= 0.0:
            weights = np.array([shell.cells[index].volume_m3 for index in indexes], dtype=float)
        weights = weights / float(weights.sum())
        for index, weight in zip(indexes, weights, strict=True):
            power[index] = layer_power * float(weight)
    heating_j = float(power.sum()) * duration
    fraction = 0.0 if neutron_j <= 0.0 else heating_j / neutron_j
    info = dict(info)
    info["heating_applied"] = True
    return power, fraction, info, "openmc"


def _write_openmc(
    destination, config, shell, library, history, neutron_j, duration, temperature_k
):
    transport = config["transport"]
    rate = _average_rate(history)
    script = destination / "openmc_model.py"
    write_model_script(
        script,
        shell,
        library,
        neutron_energy_ev=neutron_source_energy_ev(),
        particles=int(transport["openmc"]["particles"]),
        batches=int(transport["openmc"]["batches"]),
        photon_transport=bool(transport["openmc"]["photon_transport"]),
        seed=int(transport["openmc"]["seed"]),
        source_rate_per_s=rate,
        temperature_k=float(temperature_k),
    )
    backend = str(transport["backend"])
    execute = bool(transport.get("execute", False))
    status = openmc_status()
    info = {"script": script.name, "requested": backend, "executed": False, "status": status}
    runnable = status.get("importable") and "reason" not in status
    if backend == "analytical" or not execute:
        info["reason"] = "analytical deposition was used"
        return info
    if backend == "openmc" and not runnable:
        raise RuntimeError(status.get("reason", "OpenMC cannot run"))
    if backend == "auto" and not runnable:
        info["reason"] = status.get("reason", "OpenMC cannot run")
        return info
    # Export XML by executing the generated script's builder only when the
    # import works. The statepoint is optional; failure falls back in auto.
    try:
        import openmc

        previous = Path.cwd()
        os.chdir(destination)
        try:
            namespace: dict = {}
            exec(script.read_text(encoding="utf-8"), namespace)  # noqa: S102
            run = run_openmc(
                destination,
                threads=transport["openmc"].get("threads"),
                executable=transport["openmc"].get("executable"),
            )
        finally:
            os.chdir(previous)
        info.update(run)
        if run.get("executed"):
            heating = _read_heating(openmc, destination)
            info["heating"] = heating
            for key in (
                "fissions_per_source",
                "fissions_per_source_std",
                "tbr_per_source",
                "tbr_per_source_std",
                "heating_over_source",
                "heating_over_source_std",
            ):
                if key in heating:
                    info[key] = heating[key]
    except Exception as exc:
        if backend == "openmc":
            raise
        info["executed"] = False
        info["reason"] = str(exc)
    return info


def _score_total(tally, score: str) -> tuple[float, float | None]:
    """Sum a score across cells. The std is the std of that sum, not a sum of stds."""
    mean = float(np.sum(tally.get_values(scores=[score])))
    try:
        reduced = tally.summation()
        std = float(np.sum(reduced.get_values(scores=[score], value="std_dev")))
    except Exception:
        std = None
    return mean, std


def _named_tally(state, name: str, scores: tuple[str, ...]) -> tuple[float, float | None] | None:
    try:
        tally = state.get_tally(name=name)
        available = list(tally.scores)
    except Exception:
        return None
    for score in scores:
        if score not in available:
            continue
        try:
            return _score_total(tally, score)
        except Exception:
            continue
    return None


def _read_heating(openmc, destination: Path) -> dict:
    statepoints = sorted(destination.glob("statepoint.*.h5"))
    if not statepoints:
        return {"read": False, "reason": "no statepoint"}
    with openmc.StatePoint(statepoints[-1]) as state:
        heating_tally = state.get_tally(name="shell-heating")
        heating_bins = [
            float(value)
            for value in heating_tally.get_values(scores=["heating"]).reshape(-1)
        ]
        heating_std = _score_total(heating_tally, "heating")[1]
        # OpenMC 0.16 accepts the score name H3-production and stores it as (n,Xt).
        fissions = _named_tally(state, "fission-rate", ("fission",))
        tritium = _named_tally(state, "tritium-production", ("H3-production", "(n,Xt)"))
    source_ev = neutron_source_energy_ev()
    result = {
        "read": True,
        "heating_eV_per_source": heating_bins,
        "statepoint": statepoints[-1].name,
    }
    if source_ev > 0.0:
        result["heating_over_source"] = float(sum(heating_bins) / source_ev)
        if heating_std is not None:
            result["heating_over_source_std"] = float(heating_std / source_ev)
    if fissions is not None:
        result["fissions_per_source"] = fissions[0]
        if fissions[1] is not None:
            result["fissions_per_source_std"] = fissions[1]
    if tritium is not None:
        result["tbr_per_source"] = tritium[0]
        if tritium[1] is not None:
            result["tbr_per_source_std"] = tritium[1]
    return result


def _write_moose(destination, config, resolved, shell, library, source_w_m3, duration):
    path = destination / "heat_conduction.i"
    write_heat_input(
        path,
        shell,
        library,
        source_w_m3,
        duration_s=duration,
        steps=int(config["heat"]["time_steps"]),
        initial_temperature_k=resolved["heat.initial_temperature_K"],
        sink_temperature_k=resolved["heat.sink_temperature_K"],
        heat_transfer_coefficient=resolved["heat.outer_heat_transfer_coefficient_W_m2K"],
    )
    info = {"input": path.name, "executed": False}
    if not bool(config["heat"].get("execute_moose", False)):
        info["reason"] = "execute_moose is false; the native transient is the reported solve"
        return info
    info.update(run_moose(config["heat"].get("moose_executable"), path, destination))
    return info


def _average_rate(history) -> float:
    total = history.integrate("neutron_rate_per_s")
    duration = history.points[-1].time_s - history.points[0].time_s
    if duration <= 0.0:
        return 0.0
    return total / duration


def _write_history(path: Path, history) -> None:
    columns = [
        "time_s",
        "radius_m",
        "ion_density_m3",
        "ion_temperature_kev",
        "vacuum_field_t",
        "internal_field_t",
        "beta",
        "fusion_power_w",
        "alpha_power_w",
        "neutron_power_w",
        "neutron_rate_per_s",
        "tau_e_s",
        "triple_product_ratio",
        "pressure_pa",
    ]
    lines = [",".join(columns)]
    for point in history.points:
        lines.append(",".join(repr(float(getattr(point, name))) for name in columns))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_radial(path, shell, deposition, source_w, heat, stress) -> None:
    mises = _mises_per_cell(shell, stress)
    header = (
        "cell,layer,material,r_inner_m,r_outer_m,deposited_fraction,"
        "source_W,source_W_m3,temperature_K,von_mises_Pa"
    )
    lines = [header]
    for cell in shell.cells:
        stress_text = "nan" if np.isnan(mises[cell.index]) else repr(float(mises[cell.index]))
        lines.append(
            ",".join(
                [
                    str(cell.index),
                    cell.layer,
                    cell.material,
                    repr(cell.inner_radius_m),
                    repr(cell.outer_radius_m),
                    repr(float(deposition.energy_fraction[cell.index])),
                    repr(float(source_w[cell.index])),
                    repr(float(heat.source_w_m3[cell.index])),
                    repr(float(heat.temperature_k[cell.index])),
                    stress_text,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_field(destination, config, coils: CoilSet, schedule) -> None:
    field = config["field"]
    base = coils.base_currents_a()
    rows = schedule.current_scales(base.size)
    currents = base * max(rows, key=lambda row: abs(coils.on_axis_field(0.0, base * row)))
    scale = float(np.max(np.abs(rows)))
    radial = np.linspace(0.0, float(field["grid_r_max_m"]), int(field["grid_nr"]))
    axial = np.linspace(
        -float(field["grid_z_max_m"]), float(field["grid_z_max_m"]), int(field["grid_nz"])
    )
    b_r, b_z = coils.field_grid(
        radial,
        axial,
        currents,
        method=str(field["method"]),
        array_backend=str(field["array_backend"]),
    )
    np.savez(
        destination / "field_grid.npz",
        radial_m=radial,
        axial_m=axial,
        B_r_T=b_r,
        B_z_T=b_z,
        current_scale=scale,
    )


def _mises_per_cell(shell, stress) -> np.ndarray:
    values = np.full(len(shell.cells), np.nan)
    for layer in stress.layers:
        indexes = [cell.index for cell in shell.cells if cell.layer == layer.layer]
        if layer.structural:
            values[indexes] = layer.von_mises_pa
    return values


def _ledger_dict(ledger) -> dict:
    supply = ledger.supply
    return {
        "fusion_J": ledger.fusion_j,
        "alpha_J": ledger.alpha_j,
        "neutron_J": ledger.neutron_j,
        "alpha_escaped_J": ledger.alpha_escaped_j,
        "alpha_wall_J": ledger.alpha_wall_j,
        "alpha_retained_J": ledger.alpha_retained_j,
        "neutron_deposited_J": ledger.neutron_deposited_j,
        "neutron_not_deposited_J": ledger.neutron_not_deposited_j,
        "shell_heat_J": ledger.shell_heat_j,
        "electric_alpha_J": ledger.electric_alpha_j,
        "electric_shell_J": ledger.electric_shell_j,
        "fission_thermal_J": ledger.fission_j,
        "electric_fission_J": ledger.electric_fission_j,
        "thermal_efficiency": ledger.thermal_efficiency,
        "electric_gross_J": ledger.electric_gross_j,
        "thermal_ceiling_J": ledger.thermal_ceiling_j,
        "carnot_efficiency": ledger.carnot_efficiency,
        "shell_channel": ledger.shell_channel,
        "ohmic_J": supply.ohmic_j,
        "magnetic_invested_J": supply.magnetic_invested_j,
        "magnetic_returned_J": supply.magnetic_returned_j,
        "magnetic_unrecovered_J": supply.magnetic_unrecovered_j,
        "magnetic_energy_change_J": supply.magnetic_energy_change_j,
        "supply_drawn_J": supply.drawn_j,
        "net_electrical_J": ledger.net_electrical_j,
        "fusion_to_electric_efficiency": ledger.fusion_to_electric_efficiency,
    }


def _report(summary: dict) -> str:
    plasma = summary["plasma"]
    conversion = summary["conversion"]
    transport = summary["transport"]
    lines = [
        f"# {summary['name']}",
        "",
        summary["description"],
        "",
        "This pulse is a zero-dimensional compression schedule through a",
        "spherical shell. It is not a claim that the device ignites or that",
        "the net electricity is a design point.",
        "",
        f"- Fusion energy: {conversion['fusion_J']:.6g} J",
        f"- Peak fusion power: {plasma['peak_fusion_power_W']:.6g} W",
        f"- Peak ion temperature: {plasma['peak_ion_temperature_keV']:.4g} keV",
        f"- Peak triple product / 3e21 benchmark: {plasma['peak_triple_product_ratio']:.4g}",
        f"- Neutron energy deposited in the shell: {transport['deposited_fraction']:.4g} of the neutron energy",
        f"- Direct electric energy from escaped alphas: {conversion['electric_alpha_J']:.6g} J",
        f"- Electric energy counted from shell heat: {conversion['electric_shell_J']:.6g} J",
        f"- Fission thermal energy: {conversion['fission_thermal_J']:.6g} J",
        f"- Electric energy from fission: {conversion['electric_fission_J']:.6g} J",
        f"- Blanket multiplication M: {summary['hybrid']['energy_multiplication_M']:.4g}",
        f"- Carnot ceiling of that heat, not counted unless the shell channel says so: {conversion['thermal_ceiling_J']:.6g} J",
        f"- Energy drawn from the coil supply: {conversion['supply_drawn_J']:.6g} J",
        f"- Net electrical energy: {conversion['net_electrical_J']:.6g} J",
        "",
        f"Transport backend: {transport['backend']} ({transport['approximation']}).",
        f"Heat reported by: {summary['heat']['reported_by']}.",
        f"OpenMC: {transport['openmc'].get('reason', transport['openmc'].get('executed'))}.",
        f"MOOSE: {summary['heat']['moose'].get('reason', summary['heat']['moose'].get('executed'))}.",
        "",
    ]
    return "\n".join(lines)
