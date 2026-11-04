"""Numeric slots an external scalar provider may replace.

Every measured or chosen physical number that the pipeline uses is copied
into a flat slot map before any physics runs. A provider may replace a slot
by key. Geometry code then reads the resolved map, so a later FSOT module
does not have to edit the solvers.

Non-numeric choices (backend names, element symbols, flux-conservation law)
stay in the YAML. They are not slots.
"""

from __future__ import annotations

from typing import Mapping, Protocol

from solcascade.constants import JOULE_PER_EV

MATERIAL_SCALAR_FIELDS = (
    "density_kg_m3",
    "thermal_conductivity_W_mK",
    "specific_heat_J_kgK",
    "youngs_modulus_Pa",
    "poissons_ratio",
    "thermal_expansion_per_K",
    "yield_strength_Pa",
    "resistivity_ohm_m",
    "resistivity_temperature_coeff_per_K",
    "resistivity_reference_K",
    "one_group_macroscopic_total_per_m",
    "one_group_kerma_fraction",
    "one_group_fission_fraction",
    "one_group_capture_fraction",
    "one_group_tritium_fraction",
)


class ScalarProvider(Protocol):
    """Return an injected float for ``key``, or None to keep the YAML value."""

    def value(self, key: str) -> float | None: ...

    def name(self) -> str: ...


class NullProvider:
    """No external scalars. Every slot keeps its YAML value."""

    def value(self, key: str) -> None:
        return None

    def name(self) -> str:
        return "null"


class MappingProvider:
    """Inject a dict of already-computed scalars.

    The dict is data. This class does not evaluate FSOT formulas. Unknown
    keys raise, so a typo cannot disappear into an unused override.
    """

    def __init__(self, values: Mapping[str, float], allowed: set[str] | None = None):
        self._values = {str(key): float(raw) for key, raw in values.items()}
        self._allowed = allowed

    def bind(self, allowed: set[str]) -> None:
        unknown = sorted(set(self._values) - allowed)
        if unknown:
            raise KeyError(
                "Scalar provider keys are not slots in this problem: "
                + ", ".join(unknown)
            )
        self._allowed = allowed

    def value(self, key: str) -> float | None:
        if self._allowed is not None and key not in self._allowed and key in self._values:
            raise KeyError(key)
        if key not in self._values:
            return None
        return self._values[key]

    def name(self) -> str:
        return "mapping"


def collect_slots(config: Mapping, materials: Mapping) -> dict[str, float]:
    """Flatten the numeric inputs of one problem into stable keys."""
    slots: dict[str, float] = {}
    for coil in config["geometry"]["coils"]:
        name = coil["name"]
        for field in ("radius_m", "z_m", "turns", "current_a", "wire_radius_m"):
            slots[f"coils.{name}.{field}"] = float(coil[field])
    shell = config["geometry"]["shell"]
    slots["shell.inner_radius_m"] = float(shell["inner_radius_m"])
    for layer in shell["layers"]:
        slots[f"shell.layers.{layer['name']}.thickness_m"] = float(layer["thickness_m"])
    for mat_name, mat in materials.items():
        for field in MATERIAL_SCALAR_FIELDS:
            if field in mat and mat[field] is not None:
                slots[f"materials.{mat_name}.{field}"] = float(mat[field])
    plasma = config["plasma"]
    slots["plasma.ion_density_m3"] = float(plasma["ion_density_m3"])
    slots["plasma.ion_temperature_keV"] = float(plasma["ion_temperature_keV"])
    slots["plasma.deuterium_fraction"] = float(plasma["deuterium_fraction"])
    # Te = Ti, monatomic: thermal energy density is 3 n kT.
    slots["plasma.thermal_energy_density_J_m3"] = (
        3.0
        * slots["plasma.ion_density_m3"]
        * slots["plasma.ion_temperature_keV"]
        * JOULE_PER_EV
        * 1.0e3
    )
    compression = config["compression"]
    for index, time_s in enumerate(compression["time_s"]):
        slots[f"compression.time_s.{index}"] = float(time_s)
    for index, radius in enumerate(compression["plasma_radius_m"]):
        slots[f"compression.plasma_radius_m.{index}"] = float(radius)
    shared_scale = [float(scale) for scale in compression["coil_current_scale"]]
    for index, scale in enumerate(shared_scale):
        slots[f"compression.coil_current_scale.{index}"] = scale
    coil_scales = compression.get("coil_scales") or {}
    for coil in config["geometry"]["coils"]:
        name = coil["name"]
        series = coil_scales.get(name, shared_scale)
        if len(series) != len(shared_scale):
            raise ValueError(f"coil {name} scale length does not match compression.time_s")
        for index, scale in enumerate(series):
            slots[f"compression.coils.{name}.scale.{index}"] = float(scale)
    heat = config["heat"]
    slots["heat.sink_temperature_K"] = float(heat["sink_temperature_K"])
    slots["heat.outer_heat_transfer_coefficient_W_m2K"] = float(
        heat["outer_heat_transfer_coefficient_W_m2K"]
    )
    slots["heat.initial_temperature_K"] = float(heat["initial_temperature_K"])
    slots["mechanics.stress_free_temperature_K"] = float(
        config["mechanics"]["stress_free_temperature_K"]
    )
    conversion = config["conversion"]
    for field in (
        "alpha_escape_fraction",
        "alpha_wall_fraction",
        "alpha_direct_efficiency",
        "carnot_fraction",
        "cold_temperature_K",
        "coil_energy_recovery_fraction",
        "conductor_temperature_K",
        "thermal_efficiency",
    ):
        if field in conversion and conversion[field] is not None:
            slots[f"conversion.{field}"] = float(conversion[field])
    blanket = config.get("blanket") or {}
    for field in ("thickness_m", "density_g_cm3", "enrichment", "neutron_yield_scale"):
        if field in blanket and blanket[field] is not None:
            slots[f"blanket.{field}"] = float(blanket[field])
    return slots


def apply_provider(
    slots: Mapping[str, float], provider: ScalarProvider
) -> tuple[dict[str, float], dict[str, str]]:
    """Replace slots the provider fills. Return values and a provenance map."""
    if isinstance(provider, MappingProvider):
        provider.bind(set(slots))
    resolved: dict[str, float] = {}
    provenance: dict[str, str] = {}
    for key, configured in slots.items():
        injected = provider.value(key)
        if injected is None:
            resolved[key] = float(configured)
            provenance[key] = "config"
        else:
            resolved[key] = float(injected)
            provenance[key] = provider.name()
    return resolved, provenance


def overrides_only(provenance: Mapping[str, str]) -> dict[str, str]:
    return {key: origin for key, origin in provenance.items() if origin != "config"}


def reconcile_energy_density(
    slots: dict[str, float], provenance: dict[str, str]
) -> tuple[dict[str, float], dict[str, str]]:
    """If an injected energy density is present, it sets the ion temperature.

    Density stays as resolved. Temperature is then whatever that energy
    density and density imply for Te = Ti monatomic fuel. If the energy
    density was not injected, it is rewritten from the resolved density and
    temperature so the ledger stays consistent.
    """
    density_key = "plasma.ion_density_m3"
    temperature_key = "plasma.ion_temperature_keV"
    energy_key = "plasma.thermal_energy_density_J_m3"
    joule_per_kev = JOULE_PER_EV * 1.0e3
    density = slots[density_key]
    if density <= 0.0:
        raise ValueError("ion density must be positive")
    if provenance[energy_key] != "config":
        energy = slots[energy_key]
        if energy <= 0.0:
            raise ValueError("thermal energy density must be positive")
        slots[temperature_key] = energy / (3.0 * density * joule_per_kev)
        provenance[temperature_key] = provenance[energy_key] + "->temperature"
    else:
        slots[energy_key] = 3.0 * density * slots[temperature_key] * joule_per_kev
    return slots, provenance
