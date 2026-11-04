"""First-collision energy deposition in a spherical shell.

A point isotropic source sits at the center, so every ray is radial.
The uncollided collision probability in a cell is

    exp(-tau_inner) - exp(-tau_outer).

A configured kerma fraction of that collision energy is deposited in the
cell. The rest of the collision energy, and the uncollided leakage, are
reported as not deposited. Scattered neutrons are not transported again.
Optional fission, capture, and tritium fractions split those same
collisions. They do not add a second flight. This is a surrogate. The
OpenMC backend is the Monte Carlo replacement.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from solcascade.geometry.shell import Shell
from solcascade.materials.library import Material


@dataclass(frozen=True)
class Deposition:
    energy_fraction: np.ndarray
    collision_fraction: np.ndarray
    leakage_uncollided_fraction: float
    unresolved_fraction: float
    macroscopic_total_per_m: np.ndarray
    fission_per_source: np.ndarray
    capture_per_source: np.ndarray
    tritium_per_source: np.ndarray
    backend: str
    approximation: str

    def power_w(self, neutron_power_w: float) -> np.ndarray:
        return self.energy_fraction * neutron_power_w


def _channel(material: Material, value: float | None, label: str) -> float:
    if value is None:
        return 0.0
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"material {material.name} {label} must be in [0, 1]")
    return value


def first_collision(shell: Shell, library: dict[str, Material]) -> Deposition:
    fractions = np.zeros(len(shell.cells), dtype=float)
    collisions = np.zeros(len(shell.cells), dtype=float)
    fissions = np.zeros(len(shell.cells), dtype=float)
    captures = np.zeros(len(shell.cells), dtype=float)
    tritium = np.zeros(len(shell.cells), dtype=float)
    sigma = np.zeros(len(shell.cells), dtype=float)
    tau = 0.0
    for cell in shell.cells:
        material = library[cell.material]
        if material.one_group_macroscopic_total_per_m is None:
            raise ValueError(
                f"material {material.name} has no one-group macroscopic total"
            )
        if material.one_group_kerma_fraction is None:
            raise ValueError(f"material {material.name} has no one-group kerma fraction")
        macro = material.one_group_macroscopic_total_per_m
        if macro < 0.0:
            raise ValueError(f"material {material.name} macroscopic total is negative")
        kerma = material.one_group_kerma_fraction
        if not 0.0 <= kerma <= 1.0:
            raise ValueError(f"material {material.name} kerma fraction must be in [0, 1]")
        fission = _channel(material, material.one_group_fission_fraction, "fission fraction")
        capture = _channel(material, material.one_group_capture_fraction, "capture fraction")
        tritium_fraction = _channel(
            material, material.one_group_tritium_fraction, "tritium fraction"
        )
        if fission + capture + tritium_fraction > 1.0 + 1.0e-12:
            raise ValueError(
                f"material {material.name} reaction fractions sum to more than 1"
            )
        dtau = macro * cell.thickness_m
        collision = math.exp(-tau) * (-math.expm1(-dtau))
        collisions[cell.index] = collision
        fractions[cell.index] = collision * kerma
        fissions[cell.index] = collision * fission
        captures[cell.index] = collision * capture
        tritium[cell.index] = collision * tritium_fraction
        sigma[cell.index] = macro
        tau += dtau
    leakage = math.exp(-tau)
    unresolved = float(collisions.sum() - fractions.sum() + leakage)
    return Deposition(
        energy_fraction=fractions,
        collision_fraction=collisions,
        leakage_uncollided_fraction=leakage,
        unresolved_fraction=unresolved,
        macroscopic_total_per_m=sigma,
        fission_per_source=fissions,
        capture_per_source=captures,
        tritium_per_source=tritium,
        backend="analytical",
        approximation="first_collision",
    )
