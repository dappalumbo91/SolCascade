"""Material records.

Scalar properties are overwritten from the resolved slot map. Compositions
stay on the record because they are lists of element symbols, not scalars.
One-group removal data are inputs to the analytical transport surrogate.
They are not a nuclear-data library.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping


@dataclass(frozen=True)
class Material:
    name: str
    density_kg_m3: float
    thermal_conductivity_W_mK: float
    specific_heat_J_kgK: float
    youngs_modulus_Pa: float
    poissons_ratio: float
    thermal_expansion_per_K: float
    yield_strength_Pa: float
    resistivity_ohm_m: float | None
    resistivity_temperature_coeff_per_K: float | None
    resistivity_reference_K: float | None
    one_group_macroscopic_total_per_m: float | None
    one_group_kerma_fraction: float | None
    composition: tuple[tuple[str, float], ...]
    note: str
    one_group_fission_fraction: float | None = None
    one_group_capture_fraction: float | None = None
    one_group_tritium_fraction: float | None = None

    @property
    def structural(self) -> bool:
        """Liquids are entered with zero modulus and carry no elastic stress."""
        return self.youngs_modulus_Pa > 0.0

    def resistivity_at(self, temperature_k: float) -> float:
        if self.resistivity_ohm_m is None:
            raise ValueError(f"material {self.name} has no electrical resistivity")
        coeff = self.resistivity_temperature_coeff_per_K or 0.0
        reference = self.resistivity_reference_K
        if reference is None or coeff == 0.0:
            return self.resistivity_ohm_m
        return self.resistivity_ohm_m * (1.0 + coeff * (temperature_k - reference))


def _optional(slots: Mapping[str, float], key: str) -> float | None:
    if key not in slots:
        return None
    return float(slots[key])


def _symbol(entry: Mapping) -> str:
    if "nuclide" in entry:
        return str(entry["nuclide"])
    return str(entry["element"])


def material_from_slots(name: str, raw: Mapping, slots: Mapping[str, float]) -> Material:
    prefix = f"materials.{name}."
    composition = tuple(
        (_symbol(entry), float(entry["weight_fraction"])) for entry in raw["composition"]
    )
    weight = sum(fraction for _element, fraction in composition)
    if abs(weight - 1.0) > 1.0e-6:
        raise ValueError(f"material {name} weight fractions sum to {weight}, not 1")
    return Material(
        name=name,
        density_kg_m3=slots[prefix + "density_kg_m3"],
        thermal_conductivity_W_mK=slots[prefix + "thermal_conductivity_W_mK"],
        specific_heat_J_kgK=slots[prefix + "specific_heat_J_kgK"],
        youngs_modulus_Pa=slots[prefix + "youngs_modulus_Pa"],
        poissons_ratio=slots[prefix + "poissons_ratio"],
        thermal_expansion_per_K=slots[prefix + "thermal_expansion_per_K"],
        yield_strength_Pa=slots[prefix + "yield_strength_Pa"],
        resistivity_ohm_m=_optional(slots, prefix + "resistivity_ohm_m"),
        resistivity_temperature_coeff_per_K=_optional(
            slots, prefix + "resistivity_temperature_coeff_per_K"
        ),
        resistivity_reference_K=_optional(slots, prefix + "resistivity_reference_K"),
        one_group_macroscopic_total_per_m=_optional(
            slots, prefix + "one_group_macroscopic_total_per_m"
        ),
        one_group_kerma_fraction=_optional(slots, prefix + "one_group_kerma_fraction"),
        composition=composition,
        note=str(raw.get("note", "")),
        one_group_fission_fraction=_optional(slots, prefix + "one_group_fission_fraction"),
        one_group_capture_fraction=_optional(slots, prefix + "one_group_capture_fraction"),
        one_group_tritium_fraction=_optional(slots, prefix + "one_group_tritium_fraction"),
    )


def load_materials(raw_library: Mapping, slots: Mapping[str, float]) -> dict[str, Material]:
    return {
        name: material_from_slots(name, raw, slots)
        for name, raw in raw_library.items()
    }
