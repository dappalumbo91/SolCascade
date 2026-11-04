"""Build coils and a shell from a resolved slot map."""

from __future__ import annotations

import math
from typing import Mapping

from solcascade.compression.schedule import PulseSchedule, build_schedule
from solcascade.geometry.coils import Coil, CoilSet
from solcascade.geometry.shell import Shell, build_shell


def build_coils(config: Mapping, slots: Mapping[str, float]) -> CoilSet:
    coils = []
    for spec in config["geometry"]["coils"]:
        name = spec["name"]
        coils.append(
            Coil(
                name=name,
                radius_m=slots[f"coils.{name}.radius_m"],
                z_m=slots[f"coils.{name}.z_m"],
                turns=slots[f"coils.{name}.turns"],
                current_a=slots[f"coils.{name}.current_a"],
                wire_radius_m=slots[f"coils.{name}.wire_radius_m"],
                material=str(spec["material"]),
            )
        )
    return CoilSet(tuple(coils))


def build_shell_from_slots(config: Mapping, slots: Mapping[str, float]) -> Shell:
    specs = []
    for layer in config["geometry"]["shell"]["layers"]:
        specs.append(
            {
                "name": layer["name"],
                "material": layer["material"],
                "thickness_m": slots[f"shell.layers.{layer['name']}.thickness_m"],
            }
        )
    from solcascade.hybrid.blanket import fertile_layer_spec

    extra = fertile_layer_spec(config, slots)
    if extra is not None:
        specs.append(extra)
    return build_shell(
        inner_radius_m=slots["shell.inner_radius_m"],
        layer_specs=specs,
        radial_cells=int(config["geometry"]["shell"]["radial_cells"]),
    )


def build_pulse(
    config: Mapping,
    slots: Mapping[str, float],
    provenance: Mapping[str, str] | None = None,
) -> PulseSchedule:
    names = [str(coil["name"]) for coil in config["geometry"]["coils"]]
    knots = build_schedule(
        dict(slots),
        coil_names=names,
        provenance=None if provenance is None else dict(provenance),
    )
    return knots.substep(int(config["compression"]["substeps_per_interval"]))


def validate_geometry(coils: CoilSet, shell: Shell, schedule: PulseSchedule) -> None:
    if np_max(schedule.plasma_radius_m) >= shell.inner_radius_m:
        raise ValueError("plasma radius reaches the inner shell wall")
    for coil in coils.coils:
        reach = math.hypot(coil.radius_m, coil.z_m)
        clearance = shell.outer_radius_m + coil.wire_radius_m
        if reach <= clearance:
            raise ValueError(
                f"coil {coil.name} at spherical radius {reach:.4f} m intersects the shell"
            )


def np_max(values) -> float:
    return float(max(values))
