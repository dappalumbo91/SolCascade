"""Spherical radial heat conduction with a volumetric source.

The control volumes are the shell cells. Heat flow between cells uses the
series conduction resistance of the two half-cells. The inner face is
adiabatic. The outer face couples to a sink through conduction to the
surface plus Newton's law of cooling.

Backward Euler advances a pulse. With conduction removed, the update is the
local adiabatic rise q dt / (rho cp). A long pulse with conduction approaches
the steady solution of the same spatial operator.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from solcascade.geometry.shell import Shell
from solcascade.materials.library import Material


@dataclass(frozen=True)
class HeatResult:
    temperature_k: np.ndarray
    source_w_m3: np.ndarray
    time_s: float
    steps: int
    mode: str
    adiabatic_rise_k: np.ndarray


def solve_transient(
    shell: Shell,
    library: dict[str, Material],
    source_w: np.ndarray,
    duration_s: float,
    steps: int,
    initial_temperature_k: float,
    sink_temperature_k: float,
    heat_transfer_coefficient: float,
) -> HeatResult:
    if steps < 1:
        raise ValueError("heat solve needs at least one time step")
    if duration_s < 0.0:
        raise ValueError("heat duration must be >= 0")
    faces = shell.face_radii_m()
    count = len(shell.cells)
    if source_w.shape != (count,):
        raise ValueError("source vector does not match the shell cells")
    conductivity = np.array(
        [library[cell.material].thermal_conductivity_W_mK for cell in shell.cells],
        dtype=float,
    )
    density = np.array(
        [library[cell.material].density_kg_m3 for cell in shell.cells], dtype=float
    )
    specific_heat = np.array(
        [library[cell.material].specific_heat_J_kgK for cell in shell.cells],
        dtype=float,
    )
    volume = np.array([cell.volume_m3 for cell in shell.cells], dtype=float)
    source_density = np.divide(
        source_w, volume, out=np.zeros(count), where=volume > 0.0
    )
    conductance = _internal_conductance(faces, conductivity)
    outer = _outer_conductance(faces, conductivity[-1], heat_transfer_coefficient)
    capacity = density * specific_heat * volume
    if np.any(capacity <= 0.0):
        raise ValueError("density and specific heat must be positive")
    temperature = np.full(count, float(initial_temperature_k), dtype=float)
    dt = duration_s / steps if steps else 0.0
    if duration_s == 0.0 or dt == 0.0:
        rise = np.zeros(count)
    else:
        for _step in range(steps):
            temperature = _backward_euler(
                temperature, conductance, outer, capacity, source_w, sink_temperature_k, dt
            )
        rise = source_density * duration_s / (density * specific_heat)
    return HeatResult(
        temperature_k=temperature,
        source_w_m3=source_density,
        time_s=float(duration_s),
        steps=int(steps),
        mode="pulse_transient",
        adiabatic_rise_k=rise,
    )


def steady_temperature(
    shell: Shell,
    library: dict[str, Material],
    source_w: np.ndarray,
    sink_temperature_k: float,
    heat_transfer_coefficient: float,
) -> np.ndarray:
    """Steady state of the same spatial operator, for tests and checks."""
    faces = shell.face_radii_m()
    conductivity = np.array(
        [library[cell.material].thermal_conductivity_W_mK for cell in shell.cells],
        dtype=float,
    )
    conductance = _internal_conductance(faces, conductivity)
    outer = _outer_conductance(faces, conductivity[-1], heat_transfer_coefficient)
    count = len(shell.cells)
    matrix, rhs = _spatial_system(count, conductance, outer, sink_temperature_k)
    rhs = rhs + source_w
    return np.linalg.solve(matrix, rhs)


def _backward_euler(
    temperature: np.ndarray,
    conductance: np.ndarray,
    outer: float,
    capacity: np.ndarray,
    source_w: np.ndarray,
    sink: float,
    dt: float,
) -> np.ndarray:
    count = temperature.size
    matrix, rhs = _spatial_system(count, conductance, outer, sink)
    diagonal = capacity / dt
    matrix = matrix + np.diag(diagonal)
    rhs = rhs + source_w + diagonal * temperature
    return np.linalg.solve(matrix, rhs)


def _spatial_system(
    count: int, conductance: np.ndarray, outer: float, sink: float
) -> tuple[np.ndarray, np.ndarray]:
    """Matrix for H_left - H_right, moved to the left-hand side.

    The equation is written as M T = q_watts + boundary, so conduction
    leaving a hot cell is a positive diagonal.
    """
    matrix = np.zeros((count, count), dtype=float)
    rhs = np.zeros(count, dtype=float)
    for index in range(count - 1):
        link = conductance[index]
        matrix[index, index] += link
        matrix[index + 1, index + 1] += link
        matrix[index, index + 1] -= link
        matrix[index + 1, index] -= link
    matrix[-1, -1] += outer
    rhs[-1] += outer * sink
    return matrix, rhs


def _internal_conductance(faces: np.ndarray, conductivity: np.ndarray) -> np.ndarray:
    centers = 0.5 * (faces[:-1] + faces[1:])
    links = np.zeros(conductivity.size - 1, dtype=float)
    for index in range(links.size):
        face = faces[index + 1]
        area = 4.0 * math.pi * face * face
        left = conductivity[index]
        right = conductivity[index + 1]
        if left <= 0.0 or right <= 0.0 or area <= 0.0:
            continue
        gap_left = face - centers[index]
        gap_right = centers[index + 1] - face
        resistance = gap_left / (left * area) + gap_right / (right * area)
        links[index] = 1.0 / resistance
    return links


def _outer_conductance(faces: np.ndarray, conductivity: float, coefficient: float) -> float:
    if coefficient <= 0.0 or conductivity <= 0.0:
        return 0.0
    face = float(faces[-1])
    center = 0.5 * (float(faces[-2]) + face)
    area = 4.0 * math.pi * face * face
    gap = face - center
    resistance = gap / (conductivity * area) + 1.0 / (coefficient * area)
    return 1.0 / resistance
