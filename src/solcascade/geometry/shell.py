"""Spherical shell layers and the radial cells used by transport and heat."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ShellLayer:
    name: str
    material: str
    inner_radius_m: float
    outer_radius_m: float

    @property
    def thickness_m(self) -> float:
        return self.outer_radius_m - self.inner_radius_m


@dataclass(frozen=True)
class ShellCell:
    index: int
    layer: str
    material: str
    inner_radius_m: float
    outer_radius_m: float

    @property
    def center_radius_m(self) -> float:
        return 0.5 * (self.inner_radius_m + self.outer_radius_m)

    @property
    def thickness_m(self) -> float:
        return self.outer_radius_m - self.inner_radius_m

    @property
    def volume_m3(self) -> float:
        return (4.0 / 3.0) * math.pi * (
            self.outer_radius_m ** 3 - self.inner_radius_m ** 3
        )


@dataclass(frozen=True)
class Shell:
    inner_radius_m: float
    layers: tuple[ShellLayer, ...]
    cells: tuple[ShellCell, ...]

    @property
    def outer_radius_m(self) -> float:
        return self.layers[-1].outer_radius_m

    def face_radii_m(self) -> np.ndarray:
        faces = [self.cells[0].inner_radius_m]
        faces.extend(cell.outer_radius_m for cell in self.cells)
        return np.array(faces, dtype=float)

    def contains_radius(self, radius_m: float) -> bool:
        return self.inner_radius_m < radius_m < self.outer_radius_m


def build_shell(
    inner_radius_m: float,
    layer_specs: list[dict],
    radial_cells: int,
) -> Shell:
    if inner_radius_m <= 0.0:
        raise ValueError("shell inner radius must be positive")
    if radial_cells < len(layer_specs):
        raise ValueError("radial cell count must be at least the number of layers")
    layers: list[ShellLayer] = []
    cursor = float(inner_radius_m)
    for spec in layer_specs:
        thickness = float(spec["thickness_m"])
        if thickness <= 0.0:
            raise ValueError(f"layer {spec['name']} thickness must be positive")
        layers.append(
            ShellLayer(
                name=str(spec["name"]),
                material=str(spec["material"]),
                inner_radius_m=cursor,
                outer_radius_m=cursor + thickness,
            )
        )
        cursor += thickness
    # Distribute cells in proportion to thickness, at least one per layer.
    thicknesses = np.array([layer.thickness_m for layer in layers], dtype=float)
    raw = thicknesses / thicknesses.sum() * radial_cells
    counts = np.maximum(1, np.floor(raw).astype(int))
    # Fix the sum without dropping a layer below one cell.
    while int(counts.sum()) < radial_cells:
        counts[int(np.argmax(raw - counts))] += 1
    while int(counts.sum()) > radial_cells and np.any(counts > 1):
        candidates = np.where(counts > 1)[0]
        counts[candidates[int(np.argmax(counts[candidates] - raw[candidates]))]] -= 1
    cells: list[ShellCell] = []
    index = 0
    for layer, count in zip(layers, counts, strict=True):
        edges = np.linspace(layer.inner_radius_m, layer.outer_radius_m, int(count) + 1)
        for start, stop in zip(edges[:-1], edges[1:], strict=True):
            cells.append(
                ShellCell(
                    index=index,
                    layer=layer.name,
                    material=layer.material,
                    inner_radius_m=float(start),
                    outer_radius_m=float(stop),
                )
            )
            index += 1
    return Shell(
        inner_radius_m=float(inner_radius_m),
        layers=tuple(layers),
        cells=tuple(cells),
    )
