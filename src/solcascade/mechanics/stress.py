"""Elastic stress in the solid layers of a spherical shell.

Pressure stress is the Lamé solution for internal pressure and zero external
pressure. Tensile stress is positive, so the inner radial stress equals the
negative of the internal pressure.

Thermal stress uses the free hollow-sphere solution on each solid layer
separately, with that layer's modulus, Poisson ratio, and expansion
coefficient. A liquid layer (zero modulus) is skipped. Treating each layer
as free of traction at its own faces omits constraint from the layers next
to it, so the thermal stress is a per-layer estimate.

Uniform temperature above the stress-free temperature produces zero thermal
stress. Pressure and thermal stresses are superposed linearly.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from solcascade.geometry.shell import Shell, ShellLayer
from solcascade.materials.library import Material


@dataclass(frozen=True)
class LayerStress:
    layer: str
    material: str
    structural: bool
    radius_m: np.ndarray
    sigma_rr_pa: np.ndarray
    sigma_theta_pa: np.ndarray
    von_mises_pa: np.ndarray
    yield_strength_pa: float
    exceeds_yield: bool


@dataclass(frozen=True)
class StressResult:
    layers: tuple[LayerStress, ...]
    internal_pressure_pa: float
    pressure_source: str


def evaluate_stress(
    shell: Shell,
    library: dict[str, Material],
    temperature_k: np.ndarray,
    stress_free_temperature_k: float,
    internal_pressure_pa: float,
    pressure_source: str,
) -> StressResult:
    if temperature_k.shape != (len(shell.cells),):
        raise ValueError("temperature vector does not match the shell cells")
    results: list[LayerStress] = []
    for layer in shell.layers:
        material = library[layer.material]
        indexes = [cell.index for cell in shell.cells if cell.layer == layer.name]
        radii = np.array(
            [shell.cells[index].center_radius_m for index in indexes], dtype=float
        )
        temperatures = temperature_k[indexes]
        if not material.structural:
            empty = np.full(radii.shape, np.nan)
            results.append(
                LayerStress(
                    layer=layer.name,
                    material=material.name,
                    structural=False,
                    radius_m=radii,
                    sigma_rr_pa=empty,
                    sigma_theta_pa=empty,
                    von_mises_pa=empty,
                    yield_strength_pa=material.yield_strength_Pa,
                    exceeds_yield=False,
                )
            )
            continue
        radial, hoop = _superposed(
            radii,
            temperatures,
            layer,
            material,
            stress_free_temperature_k,
            internal_pressure_pa,
            shell.inner_radius_m,
            shell.outer_radius_m,
        )
        mises = np.abs(hoop - radial)
        results.append(
            LayerStress(
                layer=layer.name,
                material=material.name,
                structural=True,
                radius_m=radii,
                sigma_rr_pa=radial,
                sigma_theta_pa=hoop,
                von_mises_pa=mises,
                yield_strength_pa=material.yield_strength_Pa,
                exceeds_yield=bool(np.nanmax(mises) > material.yield_strength_Pa),
            )
        )
    return StressResult(
        layers=tuple(results),
        internal_pressure_pa=float(internal_pressure_pa),
        pressure_source=pressure_source,
    )


def lame_stress(
    radius_m: np.ndarray, inner_m: float, outer_m: float, pressure_pa: float
) -> tuple[np.ndarray, np.ndarray]:
    """Lamé radial and hoop stress. Internal pressure, zero external pressure."""
    radius = np.asarray(radius_m, dtype=float)
    a3 = inner_m ** 3
    b3 = outer_m ** 3
    r3 = radius ** 3
    denom = b3 - a3
    if denom <= 0.0:
        raise ValueError("shell radii are not ordered")
    sigma_rr = pressure_pa * a3 * (r3 - b3) / (r3 * denom)
    sigma_theta = pressure_pa * a3 * (2.0 * r3 + b3) / (2.0 * r3 * denom)
    return sigma_rr, sigma_theta


def thermal_stress_free_sphere(
    radius_m: np.ndarray,
    temperature_delta_k: np.ndarray,
    inner_m: float,
    outer_m: float,
    youngs_pa: float,
    poisson: float,
    expansion_per_k: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Free-boundary thermal stress. ``temperature_delta_k`` is T - T_ref."""
    radius = np.asarray(radius_m, dtype=float)
    delta = np.asarray(temperature_delta_k, dtype=float)
    if radius.shape != delta.shape:
        raise ValueError("radius and temperature arrays differ")
    if not 0.0 <= poisson < 0.5:
        raise ValueError("Poisson ratio must be in [0, 0.5)")
    samples_r = np.concatenate(([inner_m], radius, [outer_m]))
    # Temperature at the layer faces takes the nearest sample.
    samples_t = np.concatenate(([delta[0]], delta, [delta[-1]]))
    order = np.argsort(samples_r, kind="mergesort")
    samples_r = samples_r[order]
    samples_t = samples_t[order]
    samples_r, samples_t = _unique_sorted(samples_r, samples_t)
    cumulative = _piecewise_linear_r2_integral(samples_r, samples_t)
    total = float(cumulative[-1])
    a3 = inner_m ** 3
    b3 = outer_m ** 3
    factor = youngs_pa * expansion_per_k / (1.0 - poisson)
    integral = np.interp(radius, samples_r, cumulative)
    r3 = radius ** 3
    sigma_rr = (2.0 * factor) * (1.0 / r3) * (
        (r3 - a3) / (b3 - a3) * total - integral
    )
    sigma_theta = factor * (1.0 / r3) * (
        (2.0 * r3 + a3) / (b3 - a3) * total + integral - delta * r3
    )
    return sigma_rr, sigma_theta


def _superposed(
    radii: np.ndarray,
    temperature_k: np.ndarray,
    layer: ShellLayer,
    material: Material,
    stress_free_k: float,
    pressure_pa: float,
    shell_inner_m: float,
    shell_outer_m: float,
) -> tuple[np.ndarray, np.ndarray]:
    pressure_rr, pressure_theta = lame_stress(
        radii, shell_inner_m, shell_outer_m, pressure_pa
    )
    thermal_rr, thermal_theta = thermal_stress_free_sphere(
        radii,
        temperature_k - stress_free_k,
        layer.inner_radius_m,
        layer.outer_radius_m,
        material.youngs_modulus_Pa,
        material.poissons_ratio,
        material.thermal_expansion_per_K,
    )
    return pressure_rr + thermal_rr, pressure_theta + thermal_theta


def _unique_sorted(radius: np.ndarray, temperature: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    kept_r = [float(radius[0])]
    kept_t = [float(temperature[0])]
    for r_value, t_value in zip(radius[1:], temperature[1:], strict=True):
        if abs(float(r_value) - kept_r[-1]) <= 1.0e-15 * max(1.0, abs(kept_r[-1])):
            kept_t[-1] = 0.5 * (kept_t[-1] + float(t_value))
            continue
        kept_r.append(float(r_value))
        kept_t.append(float(t_value))
    return np.array(kept_r), np.array(kept_t)


def _piecewise_linear_r2_integral(radius: np.ndarray, temperature: np.ndarray) -> np.ndarray:
    """Running integral of T(r) r^2 for piecewise-linear T. Exact on each piece."""
    cumulative = np.zeros(radius.shape, dtype=float)
    for index in range(radius.size - 1):
        r0 = float(radius[index])
        r1 = float(radius[index + 1])
        t0 = float(temperature[index])
        t1 = float(temperature[index + 1])
        dr = r1 - r0
        slope = 0.0 if dr == 0.0 else (t1 - t0) / dr
        intercept = t0 - slope * r0

        def antiderivative(radius_value: float) -> float:
            return intercept * (radius_value ** 3 / 3.0) + slope * (radius_value ** 4 / 4.0)

        cumulative[index + 1] = cumulative[index] + antiderivative(r1) - antiderivative(r0)
    return cumulative
