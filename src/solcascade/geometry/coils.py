"""Axisymmetric circular coils and their vacuum magnetic field.

The off-axis field of one filamentary loop is the standard elliptic-integral
expression. On the axis it reduces to the closed form

    B_z = mu_0 * I * R^2 / (2 * (R^2 + z^2)^(3/2)),  B_r = 0.

``mu_0`` is the derived CODATA value from :mod:`solcascade.constants`.
A segmented Biot-Savart sum is available for checks and for an optional
CuPy grid. Self-inductance uses the low-frequency round-wire formula.
Mutual inductance uses the Neumann coaxial-loop formula.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.special import ellipe, ellipk

from solcascade.constants import VACUUM_PERMEABILITY_H_PER_M

_K2_CAP = 1.0 - 1.0e-14
_AXIS_RADIUS_FRACTION = 1.0e-8


@dataclass(frozen=True)
class Coil:
    name: str
    radius_m: float
    z_m: float
    turns: float
    current_a: float
    wire_radius_m: float
    material: str

    def __post_init__(self) -> None:
        if self.radius_m <= 0.0:
            raise ValueError(f"Coil {self.name} radius must be positive")
        if self.turns <= 0.0:
            raise ValueError(f"Coil {self.name} turns must be positive")
        if self.wire_radius_m <= 0.0 or self.wire_radius_m >= self.radius_m:
            raise ValueError(f"Coil {self.name} wire radius is not physical")


def _elliptic_pair(k2: float) -> tuple[float, float]:
    m = min(max(k2, 0.0), _K2_CAP)
    return float(ellipk(m)), float(ellipe(m))


def loop_field(radial_m: float, axial_m: float, radius_m: float, current_a: float) -> tuple[float, float]:
    """Return (B_r, B_z) in tesla at one cylindrical point, for one filament."""
    radial = float(radial_m)
    axial = float(axial_m)
    radius = float(radius_m)
    current = float(current_a)
    mu = VACUUM_PERMEABILITY_H_PER_M
    if radial < 0.0:
        raise ValueError("cylindrical radius must be >= 0")
    axis_limit = _AXIS_RADIUS_FRACTION * radius
    if radial <= axis_limit:
        denom = (radius * radius + axial * axial) ** 1.5
        return 0.0, mu * current * radius * radius / (2.0 * denom)
    distance_sq = (radius + radial) ** 2 + axial * axial
    k2 = 4.0 * radius * radial / distance_sq
    complete_k, complete_e = _elliptic_pair(k2)
    root = math.sqrt(distance_sq)
    common = mu * current / (2.0 * math.pi * root)
    radial_denom = (radius - radial) ** 2 + axial * axial
    if radial_denom <= (1.0e-12 * radius) ** 2:
        return math.nan, math.nan
    b_z = common * (
        complete_k + (radius * radius - radial * radial - axial * axial) * complete_e / radial_denom
    )
    b_r = common * (axial / radial) * (
        -complete_k + (radius * radius + radial * radial + axial * axial) * complete_e / radial_denom
    )
    return b_r, b_z


def loop_self_inductance_h(radius_m: float, wire_radius_m: float) -> float:
    """Low-frequency self-inductance of one circular turn of round wire, henry."""
    ratio = 8.0 * radius_m / wire_radius_m
    # Internal inductance 1/4 for uniform current density: ln(8R/a) - 2 + 1/4.
    return VACUUM_PERMEABILITY_H_PER_M * radius_m * (math.log(ratio) - 1.75)


def coaxial_mutual_inductance_h(
    radius_a_m: float, radius_b_m: float, separation_m: float
) -> float:
    """Mutual inductance of two coaxial filamentary loops, henry."""
    k2 = (4.0 * radius_a_m * radius_b_m) / (
        (radius_a_m + radius_b_m) ** 2 + separation_m * separation_m
    )
    if k2 >= _K2_CAP:
        raise ValueError("coaxial loops are too close for the filament mutual-inductance formula")
    complete_k, complete_e = _elliptic_pair(k2)
    k = math.sqrt(k2)
    bracket = (2.0 / k - k) * complete_k - (2.0 / k) * complete_e
    return VACUUM_PERMEABILITY_H_PER_M * math.sqrt(radius_a_m * radius_b_m) * bracket


@dataclass(frozen=True)
class CoilSet:
    coils: tuple[Coil, ...]

    def __post_init__(self) -> None:
        names = [coil.name for coil in self.coils]
        if len(names) != len(set(names)):
            raise ValueError("coil names must be unique")
        if not self.coils:
            raise ValueError("at least one coil is required")

    def base_currents_a(self) -> np.ndarray:
        return np.array([coil.current_a for coil in self.coils], dtype=float)

    def field(self, radial_m: float, axial_m: float, currents_a: np.ndarray) -> tuple[float, float]:
        """Superpose ampere-turns. ``currents_a`` is the conductor current of each coil."""
        if len(currents_a) != len(self.coils):
            raise ValueError("current vector length does not match the coil set")
        b_r = 0.0
        b_z = 0.0
        for coil, current in zip(self.coils, currents_a, strict=True):
            point_r, point_z = loop_field(
                radial_m, axial_m - coil.z_m, coil.radius_m, float(current) * coil.turns
            )
            b_r += point_r
            b_z += point_z
        return b_r, b_z

    def on_axis_field(self, axial_m: float, currents_a: np.ndarray) -> float:
        _b_r, b_z = self.field(0.0, axial_m, currents_a)
        return b_z

    def inductance_matrix_h(self) -> np.ndarray:
        """Inductance matrix in the conductor-current basis, including turn counts."""
        count = len(self.coils)
        matrix = np.zeros((count, count), dtype=float)
        for i, left in enumerate(self.coils):
            matrix[i, i] = (left.turns ** 2) * loop_self_inductance_h(
                left.radius_m, left.wire_radius_m
            )
            for j in range(i + 1, count):
                right = self.coils[j]
                mutual = left.turns * right.turns * coaxial_mutual_inductance_h(
                    left.radius_m, right.radius_m, left.z_m - right.z_m
                )
                matrix[i, j] = mutual
                matrix[j, i] = mutual
        return matrix

    def magnetic_energy_j(self, currents_a: np.ndarray, inductance_h: np.ndarray | None = None) -> float:
        matrix = self.inductance_matrix_h() if inductance_h is None else inductance_h
        current = np.asarray(currents_a, dtype=float)
        return 0.5 * float(current @ matrix @ current)

    def ohmic_power_w(
        self,
        currents_a: np.ndarray,
        resistivity_ohm_m: np.ndarray,
    ) -> float:
        power = 0.0
        for coil, current, resistivity in zip(self.coils, currents_a, resistivity_ohm_m, strict=True):
            length = coil.turns * 2.0 * math.pi * coil.radius_m
            area = math.pi * coil.wire_radius_m ** 2
            resistance = float(resistivity) * length / area
            power += float(current) ** 2 * resistance
        return power

    def field_grid(
        self,
        radial_m: np.ndarray,
        axial_m: np.ndarray,
        currents_a: np.ndarray,
        method: str = "elliptic",
        array_backend: str = "numpy",
        segments: int = 96,
    ) -> tuple[np.ndarray, np.ndarray]:
        """B_r and B_z on a tensor grid. Arrays are shaped (n_z, n_r)."""
        if method == "elliptic":
            return self._grid_elliptic(radial_m, axial_m, currents_a)
        if method == "biot_savart":
            return self._grid_biot_savart(radial_m, axial_m, currents_a, array_backend, segments)
        raise ValueError(f"unknown field method {method!r}")

    def _grid_elliptic(
        self, radial_m: np.ndarray, axial_m: np.ndarray, currents_a: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        b_r = np.zeros((axial_m.size, radial_m.size), dtype=float)
        b_z = np.zeros_like(b_r)
        for iz, axial in enumerate(axial_m):
            for ir, radial in enumerate(radial_m):
                point_r, point_z = self.field(float(radial), float(axial), currents_a)
                b_r[iz, ir] = point_r
                b_z[iz, ir] = point_z
        return b_r, b_z

    def _grid_biot_savart(
        self,
        radial_m: np.ndarray,
        axial_m: np.ndarray,
        currents_a: np.ndarray,
        array_backend: str,
        segments: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        xp = _array_module(array_backend)
        phi = xp.linspace(0.0, 2.0 * math.pi, int(segments), endpoint=False)
        dphi = 2.0 * math.pi / float(segments)
        rr, zz = xp.meshgrid(xp.asarray(radial_m, dtype=float), xp.asarray(axial_m, dtype=float))
        b_r = xp.zeros_like(rr)
        b_z = xp.zeros_like(rr)
        # Grid points sit in the x-z plane: (r, 0, z).
        mu_factor = VACUUM_PERMEABILITY_H_PER_M / (4.0 * math.pi)
        for coil, current in zip(self.coils, currents_a, strict=True):
            ampere_turns = float(current) * coil.turns
            if ampere_turns == 0.0:
                continue
            cx = coil.radius_m * xp.cos(phi)
            cy = coil.radius_m * xp.sin(phi)
            cz = xp.full_like(phi, coil.z_m)
            dlx = -coil.radius_m * xp.sin(phi) * dphi
            dly = coil.radius_m * xp.cos(phi) * dphi
            dlz = xp.zeros_like(phi)
            # Broadcast segments against the grid. Shape (n_seg, n_z, n_r).
            rx = rr[None, :, :] - cx[:, None, None]
            ry = 0.0 - cy[:, None, None]
            rz = zz[None, :, :] - cz[:, None, None]
            inv_r3 = (rx * rx + ry * ry + rz * rz) ** -1.5
            # dl × r
            cross_x = dly[:, None, None] * rz - dlz[:, None, None] * ry
            cross_y = dlz[:, None, None] * rx - dlx[:, None, None] * rz
            cross_z = dlx[:, None, None] * ry - dly[:, None, None] * rx
            scale = mu_factor * ampere_turns
            bx = scale * xp.sum(cross_x * inv_r3, axis=0)
            by = scale * xp.sum(cross_y * inv_r3, axis=0)
            bz = scale * xp.sum(cross_z * inv_r3, axis=0)
            # Cylindrical radial component in the x-z plane is the cartesian x component.
            b_r = b_r + bx
            b_z = b_z + bz
            del by
        return _to_numpy(b_r), _to_numpy(b_z)


def _array_module(name: str):
    if name == "numpy":
        return np
    if name == "cupy":
        import cupy  # type: ignore

        return cupy
    if name == "auto":
        try:
            import cupy  # type: ignore

            return cupy
        except ImportError:
            return np
    raise ValueError(f"unknown array backend {name!r}")


def _to_numpy(array: np.ndarray) -> np.ndarray:
    if isinstance(array, np.ndarray):
        return array
    return np.asarray(array.get())
