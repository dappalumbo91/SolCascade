"""Fertile zone wrapped outside the existing shell.

The lithium-lead absorber stays the tritium breeding layer. U-238 or Th-232
is an added outer zone. ``enrichment`` is the weight fraction of that fertile
nuclide, not a U-235 enrichment. Density enters as g/cm3 and is converted to
kg/m3 at this boundary.

The analytical fission count is a first-collision tally. It does not build the
cascade that published hybrid reactors quote as an energy multiplication of
10 to 90. Those figures are literature anchors. They are not installed here.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Mapping

import numpy as np

from solcascade.constants import (
    AVOGADRO_PER_MOL,
    MOLAR_MASS_TH232_G,
    MOLAR_MASS_U238_G,
    dt_reaction,
    fission_recoverable_energy_j,
)
from solcascade.geometry.shell import Shell
from solcascade.hybrid.backends import write_fission_backend
from solcascade.hybrid.depletion import deplete_pulses
from solcascade.materials.library import Material
from solcascade.transport.analytical import Deposition

FUEL_TYPES = ("u238", "th232", "none")
FISSION_BACKENDS = ("serpent", "mcnp", "griffin", "none")
_FUEL_MATERIAL = {"u238": "uranium_238", "th232": "thorium_232"}
_FUEL_NUCLIDE = {"u238": "U238", "th232": "Th232"}
_MOLAR_MASS = {"u238": MOLAR_MASS_U238_G, "th232": MOLAR_MASS_TH232_G}

LITERATURE_NOTE = (
    "Published hybrid energy multiplications of about 10 to 90, including "
    "FFHER near 13 and FDS-EM near 90, are scale anchors from the literature. "
    "They are not the multiplication computed for this pulse."
)


@dataclass(frozen=True)
class Blanket:
    fuel_type: str
    thickness_m: float
    density_kg_m3: float
    enrichment: float
    neutron_yield_scale: float
    thermal_efficiency: float
    backend: str
    execute: bool
    pulses: int
    interpulse_s: float

    @property
    def active(self) -> bool:
        return self.fuel_type != "none"

    @property
    def material_name(self) -> str | None:
        if not self.active:
            return None
        return _FUEL_MATERIAL[self.fuel_type]


def blanket_from_config(config: Mapping, slots: Mapping[str, float]) -> Blanket:
    """Validate the hybrid YAML. A missing blanket section is pure fusion."""
    raw = config.get("blanket")
    fission = config.get("fission") or {}
    backend = str(fission.get("backend", "none"))
    if backend not in FISSION_BACKENDS:
        raise ValueError(
            "fission.backend must be one of " + ", ".join(FISSION_BACKENDS)
        )
    pulses = int(fission.get("pulses", 0))
    interpulse_s = float(fission.get("interpulse_s", 0.0))
    execute = bool(fission.get("execute", False))
    if pulses < 0:
        raise ValueError("fission.pulses must be >= 0")
    if interpulse_s < 0.0:
        raise ValueError("fission.interpulse_s must be >= 0")
    efficiency = float(slots.get("conversion.thermal_efficiency", 0.0))
    if "conversion.thermal_efficiency" not in slots and "thermal_efficiency" in (
        config.get("conversion") or {}
    ):
        efficiency = float(config["conversion"]["thermal_efficiency"])
    if not 0.0 <= efficiency <= 1.0:
        raise ValueError("conversion.thermal_efficiency must be in [0, 1]")
    if not raw:
        return Blanket(
            "none", 0.0, 0.0, 1.0, 1.0, efficiency, backend, execute, pulses, interpulse_s
        )
    fuel = str(raw.get("fuel_type", "none"))
    if fuel not in FUEL_TYPES:
        raise ValueError("blanket.fuel_type must be one of " + ", ".join(FUEL_TYPES))
    if fuel == "none":
        return Blanket(
            "none", 0.0, 0.0, 1.0, 1.0, efficiency, backend, execute, pulses, interpulse_s
        )
    thickness = float(slots["blanket.thickness_m"])
    density_g_cm3 = float(slots["blanket.density_g_cm3"])
    enrichment = float(slots["blanket.enrichment"])
    scale = float(slots.get("blanket.neutron_yield_scale", 1.0))
    if thickness <= 0.0:
        raise ValueError("blanket.thickness_m must be positive")
    if density_g_cm3 <= 0.0:
        raise ValueError("blanket.density_g_cm3 must be positive")
    if not 0.0 < enrichment <= 1.0:
        raise ValueError("blanket.enrichment must be in (0, 1]")
    if scale <= 0.0:
        raise ValueError("blanket.neutron_yield_scale must be positive")
    if "thermal_efficiency" not in (config.get("conversion") or {}):
        raise ValueError("a fertile blanket requires conversion.thermal_efficiency")
    return Blanket(
        fuel,
        thickness,
        density_g_cm3 * 1000.0,
        enrichment,
        scale,
        efficiency,
        backend,
        execute,
        pulses,
        interpulse_s,
    )


def fertile_layer_spec(config: Mapping, slots: Mapping[str, float]) -> dict | None:
    blanket = blanket_from_config(config, slots)
    if not blanket.active:
        return None
    return {
        "name": "fertile",
        "material": blanket.material_name,
        "thickness_m": blanket.thickness_m,
    }


def scale_fertile_material(
    config: Mapping, slots: Mapping[str, float], library: dict[str, Material]
) -> dict[str, Material]:
    """Copy the fertile record at the zone density and enrichment.

    Macroscopic totals scale with fertile number density. Reaction fractions
    stay at the 14.1 MeV point values. Enrichment below 1 dilutes the fertile
    nuclide with Pb-208 so the weight fractions sum to 1. The surrogate total
    does not add a lead cross section.
    """
    blanket = blanket_from_config(config, slots)
    if not blanket.active or blanket.material_name is None:
        return library
    reference = library[blanket.material_name]
    if reference.one_group_macroscopic_total_per_m is None:
        raise ValueError(f"{reference.name} has no one-group macroscopic total")
    if reference.density_kg_m3 <= 0.0:
        raise ValueError(f"{reference.name} density must be positive")
    number_scale = (blanket.enrichment * blanket.density_kg_m3) / reference.density_kg_m3
    nuclide = _FUEL_NUCLIDE[blanket.fuel_type]
    if blanket.enrichment == 1.0:
        composition = ((nuclide, 1.0),)
    else:
        composition = ((nuclide, blanket.enrichment), ("Pb208", 1.0 - blanket.enrichment))
    scaled = replace(
        reference,
        density_kg_m3=blanket.density_kg_m3,
        composition=composition,
        one_group_macroscopic_total_per_m=reference.one_group_macroscopic_total_per_m
        * number_scale,
    )
    updated = dict(library)
    updated[reference.name] = scaled
    return updated


def _source_neutrons(neutron_j: float) -> float:
    energy = dt_reaction().neutron_energy_joule
    if neutron_j <= 0.0 or energy <= 0.0:
        return 0.0
    return neutron_j / energy


def _rates_from_deposition(deposition: Deposition) -> tuple[float, float, float]:
    return (
        float(np.sum(deposition.fission_per_source)),
        float(np.sum(deposition.capture_per_source)),
        float(np.sum(deposition.tritium_per_source)),
    )


def multiplication(deposition: Deposition, fissions_per_source: float) -> float:
    """Kinetic kerma plus fission Q, per source neutron.

    The yield scale is not in this ratio. Fission joules on the ledger do
    carry that scale, so joules can move while this M stays put.
    """
    kinetic = float(np.sum(deposition.energy_fraction))
    neutron_j = dt_reaction().neutron_energy_joule
    if neutron_j <= 0.0:
        return kinetic
    return kinetic + fissions_per_source * fission_recoverable_energy_j() / neutron_j


def _fertile_atoms(shell: Shell, material: Material, enrichment: float, molar_mass_g: float) -> float:
    volume = sum(cell.volume_m3 for cell in shell.cells if cell.material == material.name)
    mass_kg = material.density_kg_m3 * volume * enrichment
    return mass_kg / (molar_mass_g * 1.0e-3) * AVOGADRO_PER_MOL


def fission_power_w(shell: Shell, deposition: Deposition, fission_j: float, duration_s: float) -> np.ndarray:
    power = np.zeros(len(shell.cells), dtype=float)
    if duration_s <= 0.0 or fission_j <= 0.0:
        return power
    weights = np.array(deposition.fission_per_source, dtype=float)
    total = float(weights.sum())
    if total <= 0.0:
        return power
    power[:] = weights / total * (fission_j / duration_s)
    return power


def run_hybrid(
    config: Mapping,
    slots: Mapping[str, float],
    shell: Shell,
    library: dict[str, Material],
    deposition: Deposition,
    neutron_j: float,
    openmc_info: Mapping,
    workdir: Path,
) -> dict:
    blanket = blanket_from_config(config, slots)
    fissions, captures, tbr = _rates_from_deposition(deposition)
    approximation = "first_collision"
    heating_ratio = None
    if blanket.active:
        fission_tally = openmc_info.get("fissions_per_source")
        if fission_tally is not None:
            fissions = float(fission_tally)
            approximation = "openmc"
        tbr_tally = openmc_info.get("tbr_per_source")
        if tbr_tally is not None:
            tbr = float(tbr_tally)
            approximation = "openmc"
    if openmc_info.get("heating_over_source") is not None:
        heating_ratio = float(openmc_info["heating_over_source"])
    source = _source_neutrons(neutron_j) * blanket.neutron_yield_scale
    fission_j = fissions * source * fission_recoverable_energy_j()
    # Capture and fission used by depletion are the fertile-zone reactions.
    # The LiPb tritium tally is not a capture on the fertile parent.
    fertile_captures = captures
    if blanket.active and blanket.material_name is not None:
        indexes = [
            cell.index for cell in shell.cells if cell.material == blanket.material_name
        ]
        fertile_captures = float(sum(deposition.capture_per_source[index] for index in indexes))
        # The OpenMC fission score, when present, already replaced ``fissions``.
        # Capture stays on the first-collision surrogate. No capture tally is read.
    inventory = None
    material = None
    if blanket.active and blanket.material_name is not None:
        material = library[blanket.material_name]
        if blanket.pulses > 0:
            atoms = _fertile_atoms(
                shell, material, blanket.enrichment, _MOLAR_MASS[blanket.fuel_type]
            )
            one_pulse_neutrons = _source_neutrons(neutron_j) * blanket.neutron_yield_scale
            inventory = deplete_pulses(
                blanket.fuel_type,
                atoms,
                fissions * one_pulse_neutrons,
                fertile_captures * one_pulse_neutrons,
                blanket.pulses,
                blanket.interpulse_s,
            )
            inventory["note"] = (
                "First-collision (n,gamma) at 14.1 MeV is a small capture. "
                "Pu-239 and U-233 breeding in a real blanket is carried by "
                "slowed neutrons. This inventory does not add that cascade. "
                "interpulse_s is a chosen gap, not a measured repetition rate."
            )
            (workdir / "isotopic_inventory.json").write_text(
                json.dumps(inventory, indent=2), encoding="utf-8"
            )
    neutron_mev = dt_reaction().neutron_energy_kev / 1.0e3
    if blanket.active:
        backend = write_fission_backend(
            workdir,
            blanket.backend,
            blanket.execute,
            shell,
            material,
            neutron_mev,
            blanket.fuel_type,
        )
    else:
        backend = {
            "backend": blanket.backend,
            "executed": False,
            "reason": "no fertile blanket",
        }
    return {
        "fuel_type": blanket.fuel_type,
        "enrichment": blanket.enrichment if blanket.active else None,
        "enrichment_meaning": "weight fraction of the named fertile nuclide",
        "density_kg_m3": blanket.density_kg_m3 if blanket.active else None,
        "thickness_m": blanket.thickness_m if blanket.active else None,
        "approximation": approximation if blanket.active else "none",
        "fissions_per_source_neutron": fissions if blanket.active else 0.0,
        "captures_per_source_neutron": fertile_captures if blanket.active else 0.0,
        "tbr_per_source_neutron": tbr,
        "energy_multiplication_M": multiplication(deposition, fissions if blanket.active else 0.0),
        "neutron_yield_scale": blanket.neutron_yield_scale,
        "fission_thermal_J": fission_j if blanket.active else 0.0,
        "thermal_efficiency": blanket.thermal_efficiency if blanket.active else 0.0,
        "openmc_heating_over_source": heating_ratio,
        "literature_note": LITERATURE_NOTE,
        "depletion": None
        if inventory is None
        else {
            "pulses": inventory["pulses"],
            "interpulse_s": inventory["interpulse_s"],
            "fissioned": inventory["fissioned"],
            "inventory": inventory["inventory"],
            "file": "isotopic_inventory.json",
        },
        "fission_backend": backend,
    }
