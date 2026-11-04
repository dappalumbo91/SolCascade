"""Piecewise-linear plasma radius and coil-current scale.

The schedule is an input. The solver does not search it for yield.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PulseSchedule:
    time_s: np.ndarray
    plasma_radius_m: np.ndarray
    coil_current_scale: np.ndarray
    per_coil_scale: np.ndarray | None = None

    def __post_init__(self) -> None:
        count = self.time_s.size
        if count < 2:
            raise ValueError("compression schedule needs at least two times")
        if self.coil_current_scale.ndim != 1:
            raise ValueError("shared coil current scale is one series")
        if self.plasma_radius_m.size != count or self.coil_current_scale.size != count:
            raise ValueError("compression columns must share one time base")
        if np.any(np.diff(self.time_s) <= 0.0):
            raise ValueError("compression times must be strictly increasing")
        if np.any(self.plasma_radius_m <= 0.0):
            raise ValueError("plasma radius must stay positive")
        if np.any(self.coil_current_scale < 0.0):
            raise ValueError("coil current scale must be >= 0")
        if self.per_coil_scale is not None:
            if self.per_coil_scale.ndim != 2 or self.per_coil_scale.shape[0] != count:
                raise ValueError("per-coil scales need one row per time and one column per coil")
            if self.per_coil_scale.shape[1] < 1:
                raise ValueError("per-coil scales need at least one coil")
            if np.any(self.per_coil_scale < 0.0):
                raise ValueError("per-coil current scale must be >= 0")

    def current_scales(self, n_coils: int) -> np.ndarray:
        """Scale of each coil at each knot. Shape ``(n_times, n_coils)``."""
        if n_coils < 1:
            raise ValueError("coil count must be positive")
        if self.per_coil_scale is None:
            return np.repeat(self.coil_current_scale[:, None], n_coils, axis=1)
        if self.per_coil_scale.shape[1] != n_coils:
            raise ValueError("per-coil scale columns do not match the coil set")
        return self.per_coil_scale

    @property
    def duration_s(self) -> float:
        return float(self.time_s[-1] - self.time_s[0])

    def substep(self, substeps_per_interval: int) -> "PulseSchedule":
        if substeps_per_interval < 1:
            raise ValueError("substeps_per_interval must be >= 1")
        pieces_t = []
        pieces_r = []
        pieces_s = []
        for index in range(self.time_s.size - 1):
            grid = np.linspace(
                self.time_s[index],
                self.time_s[index + 1],
                substeps_per_interval + 1,
                endpoint=True,
            )
            if index > 0:
                grid = grid[1:]
            pieces_t.append(grid)
            pieces_r.append(np.interp(grid, self.time_s, self.plasma_radius_m))
            pieces_s.append(np.interp(grid, self.time_s, self.coil_current_scale))
        per_coil = None
        if self.per_coil_scale is not None:
            columns = []
            for column in range(self.per_coil_scale.shape[1]):
                pieces = []
                for index in range(self.time_s.size - 1):
                    grid = np.linspace(
                        self.time_s[index],
                        self.time_s[index + 1],
                        substeps_per_interval + 1,
                        endpoint=True,
                    )
                    if index > 0:
                        grid = grid[1:]
                    pieces.append(np.interp(grid, self.time_s, self.per_coil_scale[:, column]))
                columns.append(np.concatenate(pieces))
            per_coil = np.column_stack(columns)
        return PulseSchedule(
            time_s=np.concatenate(pieces_t),
            plasma_radius_m=np.concatenate(pieces_r),
            coil_current_scale=np.concatenate(pieces_s),
            per_coil_scale=per_coil,
        )


def build_schedule(
    slots: dict[str, float],
    prefix: str = "compression",
    coil_names: list[str] | None = None,
    provenance: dict[str, str] | None = None,
) -> PulseSchedule:
    times = _indexed(slots, f"{prefix}.time_s")
    radii = _indexed(slots, f"{prefix}.plasma_radius_m")
    scales = _indexed(slots, f"{prefix}.coil_current_scale")
    per_coil = _per_coil_matrix(slots, prefix, coil_names or [], len(times), scales, provenance)
    return PulseSchedule(
        time_s=np.array(times, dtype=float),
        plasma_radius_m=np.array(radii, dtype=float),
        coil_current_scale=np.array(scales, dtype=float),
        per_coil_scale=per_coil,
    )


def _per_coil_matrix(
    slots: dict[str, float],
    prefix: str,
    coil_names: list[str],
    count: int,
    shared: list[float],
    provenance: dict[str, str] | None,
) -> np.ndarray | None:
    """Independent waveforms. Absent keys keep the shared scale."""
    if not coil_names:
        return None
    rows: list[list[float]] = []
    any_key = False
    for name in coil_names:
        column: list[float] = []
        for index, shared_scale in enumerate(shared):
            key = f"{prefix}.coils.{name}.scale.{index}"
            shared_key = f"{prefix}.coil_current_scale.{index}"
            own_injected = (
                provenance is not None and provenance.get(key, "config") != "config"
            )
            shared_injected = (
                provenance is not None
                and provenance.get(shared_key, "config") != "config"
            )
            if key in slots and (own_injected or not shared_injected):
                any_key = True
                column.append(float(slots[key]))
            else:
                column.append(float(shared_scale))
        if len(column) != count:
            raise ValueError(f"{name} scale length does not match the schedule")
        rows.append(column)
    if not any_key:
        return None
    return np.array(rows, dtype=float).T


def _indexed(slots: dict[str, float], prefix: str) -> list[float]:
    values: list[tuple[int, float]] = []
    head = prefix + "."
    for key, value in slots.items():
        if key.startswith(head):
            values.append((int(key[len(head) :]), float(value)))
    if not values:
        raise KeyError(prefix)
    values.sort()
    indexes = [index for index, _value in values]
    if indexes != list(range(len(indexes))):
        raise ValueError(f"{prefix} indexes are not a contiguous range starting at 0")
    return [value for _index, value in values]
