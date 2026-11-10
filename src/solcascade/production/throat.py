"""Design a bore throat stronger than the frozen internal field.

The checked-in coils freeze the internal field at the first current knot.
A particle mirrors only when the strongest field on its way to the wall
exceeds the field it occupies. This module places one symmetric pair in
the bore, on the checked-in wire, and solves its ampere-turns with the
coil field and the coil inductance. The checked-in YAML is not rewritten.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from solcascade.assemble import build_shell_from_slots
from solcascade.geometry.coils import Coil, CoilSet
from solcascade.production.nature import isotropic_loss_fraction
from solcascade.production.refine import (
    Problem,
    fuel_of,
    load_checked_in,
    project_root,
    schedule_yield,
)
from solcascade.production.surface import build_surface

# Final on-axis samples from the compressed plasma surface to the inner wall.
PATH_SAMPLES = 2001
# Coarse samples used only to rank candidate windings.
SEARCH_SAMPLES = 41
HALF_LOSS = 0.5
TENTH_LOSS = 0.1


@dataclass(frozen=True)
class ThroatScale:
    """One ampere-turn solution on the designed bore pair."""

    label: str
    mirror_ratio: float
    loss_fraction: float
    ampere_turns: float
    turns: float
    current_a: float
    stored_j: float
    pulsed_supply_j: float
    persistent_draw_j: float
    throat_field_t: float
    internal_field_t: float
    trapped_electric_j: float


@dataclass(frozen=True)
class ExternalScore:
    """A pair that stays outside the shell, scored with the same ratio."""

    vacuum_rise: float
    asymptotic_ratio: float
    radius_m: float
    z_m: float
    ampere_turns: float
    stored_j: float
    pulsed_supply_j: float
    steady_reaches: bool


@dataclass(frozen=True)
class ThroatReport:
    wire_m: float
    conductor_a: float
    radius_m: float
    z_m: float
    compression: float
    vacuum_rise: float
    bare_internal_t: float
    bare_ratio: float
    recovery: float
    walk_electric_j: float
    checked_fusion_j: float
    threshold: ThroatScale
    design: ThroatScale
    tenth: ThroatScale
    external: ExternalScore


def default_output_path():
    return project_root() / "runs" / "production" / "throat.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _side(left: float, right: float) -> str:
    if left > right:
        return "above"
    if left < right:
        return "below"
    return "level with"


def mirror_ratio_for_loss(loss_fraction: float) -> float:
    """Invert :func:`isotropic_loss_fraction`.

    Half the distribution is trapped when the loss fraction is one half,
    and that mirror ratio is 4/3. A loss fraction of 1 is mirror ratio 1.
    """
    if not 0.0 < loss_fraction <= 1.0:
        raise ValueError("loss fraction is between 0 and 1")
    if loss_fraction == 1.0:
        return 1.0
    trapped = 1.0 - loss_fraction
    return 1.0 / (1.0 - trapped * trapped)


def _pair(radius_m: float, z_m: float, wire_m: float, turns: float, current_a: float) -> CoilSet:
    return CoilSet(
        (
            Coil("throat_minus", radius_m, -z_m, turns, current_a, wire_m, "copper"),
            Coil("throat_plus", radius_m, z_m, turns, current_a, wire_m, "copper"),
        )
    )


def _conductor(ampere_turns: float, conductor_a: float) -> tuple[float, float]:
    """Split ampere-turns onto the checked-in conductor.

    One turn carries the current when that current stays at or below the
    checked-in conductor current. A larger solution keeps that current and
    adds turns.
    """
    if ampere_turns <= 0.0:
        raise ValueError("ampere-turns must be positive")
    if ampere_turns <= conductor_a:
        return 1.0, float(ampere_turns)
    return float(ampere_turns / conductor_a), float(conductor_a)


def _axis(coils: CoilSet, currents: np.ndarray, stations: np.ndarray) -> np.ndarray:
    return np.array(
        [float(coils.on_axis_field(float(station), currents)) for station in stations],
        dtype=float,
    )


def _steady_ampere_turns(
    throat_mid_t: float,
    throat_path_t: np.ndarray,
    checked_mid_t: float,
    checked_path_t: np.ndarray,
    compression: float,
    target: float,
) -> float:
    """Smallest positive ampere-turns whose mirror ratio meets ``target``.

    The throat is on at freeze-in, so its midplane field is multiplied by
    the compression. A non-positive denominator means that station never
    reaches the target.
    """
    if target < 1.0:
        raise ValueError("a throat target is at least 1")
    best = math.inf
    for throat_t, checked_t in zip(throat_path_t, checked_path_t, strict=True):
        numer = target * compression * checked_mid_t - checked_t
        denom = throat_t - target * compression * throat_mid_t
        if numer <= 0.0:
            return 0.0
        if denom <= 0.0:
            continue
        best = min(best, numer / denom)
    return best


def _pulsed_ampere_turns(
    throat_path_t: np.ndarray,
    checked_path_t: np.ndarray,
    internal_t: float,
    target: float,
) -> float:
    """Ampere-turns for a throat that is off at freeze-in."""
    need = target * internal_t
    best = math.inf
    for throat_t, checked_t in zip(throat_path_t, checked_path_t, strict=True):
        gap = need - checked_t
        if gap <= 0.0:
            return 0.0
        if throat_t <= 0.0:
            continue
        best = min(best, gap / throat_t)
    return best


def _bore_corner(wire_m: float, plasma_m: float, wall_m: float) -> tuple[float, float]:
    """Smallest centerline the wire allows, against the inner wall.

    The coil constructor requires the centerline radius to exceed the wire
    radius. The conductor stays strictly outside the initial plasma sphere
    and strictly inside the inner wall.
    """
    radius = math.nextafter(wire_m, math.inf)
    reach = math.nextafter(wall_m - wire_m, 0.0)
    remain = reach * reach - radius * radius
    if remain <= 0.0:
        raise RuntimeError("the checked-in wire does not fit in the bore")
    z_m = math.sqrt(remain)
    while math.hypot(radius, z_m) + wire_m >= wall_m:
        z_m = math.nextafter(z_m, 0.0)
    if math.hypot(radius, z_m) - wire_m <= plasma_m or z_m < wire_m:
        raise RuntimeError("the bore corner is not an open path")
    return radius, z_m


def _external_corner(wire_m: float, clearance_m: float) -> tuple[float, float]:
    """Smallest loop just outside the shell clearance, near the axis."""
    radius = math.nextafter(wire_m, math.inf)
    remain = clearance_m * clearance_m - radius * radius
    z_m = math.sqrt(remain) if remain > 0.0 else 0.0
    while math.hypot(radius, z_m) <= clearance_m:
        z_m = math.nextafter(z_m, math.inf)
    return radius, z_m


def _bore_feasible(radius_m: float, z_m: float, wire_m: float, plasma_m: float, wall_m: float) -> bool:
    if radius_m <= wire_m or z_m < wire_m:
        return False
    reach = math.hypot(radius_m, z_m)
    return reach - wire_m > plasma_m and reach + wire_m < wall_m


def _stations(path: np.ndarray, z_m: float, plasma_m: float, wall_m: float) -> np.ndarray:
    stations = np.unique(np.concatenate([path, np.array([z_m], dtype=float)]))
    return stations[(stations >= plasma_m) & (stations <= wall_m)]


def _scale(
    label: str,
    radius_m: float,
    z_m: float,
    wire_m: float,
    conductor_a: float,
    ampere_turns: float,
    checked_mid_t: float,
    checked_path_t: np.ndarray,
    path: np.ndarray,
    compression: float,
    recovery: float,
    walk_electric_j: float,
) -> ThroatScale:
    turns, current = _conductor(ampere_turns, conductor_a)
    winding = _pair(radius_m, z_m, wire_m, turns, current)
    currents = np.array([current, current], dtype=float)
    stored = float(winding.magnetic_energy_j(currents))
    throat_path = _axis(winding, currents, path)
    throat_mid = float(winding.on_axis_field(0.0, currents))
    internal = (checked_mid_t + throat_mid) * compression
    # The occupied field is the frozen internal field. The throat is the
    # strongest vacuum field on the escape path, including the checked-in coils.
    throat_field = float(np.max(checked_path_t + throat_path))
    if internal <= 0.0:
        raise RuntimeError("internal field is not positive")
    ratio = throat_field / internal
    # Same numerical tie as axial_peak: a rise this close to 1 has no trapped set.
    if ratio <= 1.0 + 1.0e-9:
        ratio = 1.0
        throat_field = float(internal)
    occupied = ratio if ratio > 1.0 else 1.0
    loss = isotropic_loss_fraction(occupied)
    trapped = 1.0 - loss
    return ThroatScale(
        label=label,
        mirror_ratio=float(ratio),
        loss_fraction=float(loss),
        ampere_turns=float(ampere_turns),
        turns=float(turns),
        current_a=float(current),
        stored_j=stored,
        pulsed_supply_j=float((1.0 - recovery) * stored),
        persistent_draw_j=0.0,
        throat_field_t=throat_field,
        internal_field_t=float(internal),
        trapped_electric_j=float(trapped ** 2 * walk_electric_j),
    )


def _bore_energy(
    radius_m: float,
    z_m: float,
    wire_m: float,
    target: float,
    path: np.ndarray,
    checked_path_t: np.ndarray,
    checked_mid_t: float,
    compression: float,
) -> tuple[float, float]:
    """Return steady ampere-turns and stored energy for one bore station."""
    unit = _pair(radius_m, z_m, wire_m, 1.0, 1.0)
    probe = np.ones(2, dtype=float)
    stations = _stations(path, z_m, float(path[0]), float(path[-1]))
    throat_path = _axis(unit, probe, stations)
    checked_on = np.interp(stations, path, checked_path_t)
    throat_mid = float(unit.on_axis_field(0.0, probe))
    ampere_turns = _steady_ampere_turns(
        throat_mid, throat_path, checked_mid_t, checked_on, compression, target
    )
    if not math.isfinite(ampere_turns) or ampere_turns <= 0.0:
        return ampere_turns, math.inf
    stored = float(unit.magnetic_energy_j(probe * ampere_turns))
    return ampere_turns, stored


def _search_bore(
    wire_m: float,
    plasma_m: float,
    wall_m: float,
    target: float,
    path: np.ndarray,
    checked_path_t: np.ndarray,
    checked_mid_t: float,
    compression: float,
) -> tuple[float, float, float]:
    """Lowest steady store in the bore. Returns radius, z, and ampere-turns."""
    search = np.linspace(float(path[0]), float(path[-1]), SEARCH_SAMPLES)
    corner_r, corner_z = _bore_corner(wire_m, plasma_m, wall_m)
    best_energy = math.inf
    best: tuple[float, float, float] | None = None

    def consider(radius_m: float, z_m: float, stations: np.ndarray) -> None:
        nonlocal best, best_energy
        if not _bore_feasible(radius_m, z_m, wire_m, plasma_m, wall_m):
            return
        try:
            ampere_turns, stored = _bore_energy(
                radius_m,
                z_m,
                wire_m,
                target,
                stations,
                np.interp(stations, path, checked_path_t),
                checked_mid_t,
                compression,
            )
        except ValueError:
            return
        if stored < best_energy:
            best_energy = stored
            best = (radius_m, z_m, ampere_turns)

    consider(corner_r, corner_z, path)
    outer = wall_m - wire_m
    radii = np.geomspace(corner_r * (1.0 + 1.0e-3), outer * 0.98, 18)
    checked_search = np.interp(search, path, checked_path_t)
    for radius_m in radii:
        inner = plasma_m + wire_m
        z_lo = wire_m
        if radius_m < inner:
            remain = inner * inner - radius_m * radius_m
            if remain <= 0.0:
                continue
            z_lo = max(z_lo, math.sqrt(remain))
        remain = outer * outer - radius_m * radius_m
        if remain <= 0.0:
            continue
        z_hi = math.sqrt(remain)
        if z_hi <= z_lo:
            continue
        for z_m in np.linspace(z_lo, math.nextafter(z_hi, 0.0), 12):
            if not _bore_feasible(float(radius_m), float(z_m), wire_m, plasma_m, wall_m):
                continue
            try:
                coarse = _bore_energy(
                    float(radius_m),
                    float(z_m),
                    wire_m,
                    target,
                    search,
                    checked_search,
                    checked_mid_t,
                    compression,
                )[1]
            except ValueError:
                continue
            if coarse < best_energy:
                consider(float(radius_m), float(z_m), path)
    if best is None:
        raise RuntimeError("no bore pair reaches the half-trapped ratio")
    return best


def _unit_fields(
    radius_m: float, z_m: float, wire_m: float, path: np.ndarray
) -> tuple[CoilSet, float, np.ndarray]:
    unit = _pair(radius_m, z_m, wire_m, 1.0, 1.0)
    probe = np.ones(2, dtype=float)
    throat_mid = float(unit.on_axis_field(0.0, probe))
    throat_path = _axis(unit, probe, path)
    return unit, throat_mid, throat_path


def _search_external(
    wire_m: float,
    clearance_m: float,
    target: float,
    path: np.ndarray,
    checked_path_t: np.ndarray,
    checked_mid_t: float,
    compression: float,
    recovery: float,
) -> ExternalScore:
    """Strongest external vacuum rise, and the cheapest pulsed store."""
    search = np.linspace(float(path[0]), float(path[-1]), 81)
    checked_search = np.interp(search, path, checked_path_t)
    bare_internal = checked_mid_t * compression
    best_rise = 0.0
    best_rise_at: tuple[float, float] | None = None
    best_pulsed = math.inf
    best_pulsed_at: tuple[float, float, float] | None = None

    def rank(radius_m: float, z_m: float, stations: np.ndarray, checked_on: np.ndarray) -> None:
        nonlocal best_rise, best_rise_at, best_pulsed, best_pulsed_at
        if radius_m <= wire_m or z_m < wire_m or math.hypot(radius_m, z_m) <= clearance_m:
            return
        try:
            unit, throat_mid, throat_path = _unit_fields(radius_m, z_m, wire_m, stations)
        except ValueError:
            return
        if throat_mid <= 0.0:
            return
        rise = float(np.max(throat_path) / throat_mid)
        if rise > best_rise:
            best_rise = rise
            best_rise_at = (radius_m, z_m)
        ampere_turns = _pulsed_ampere_turns(throat_path, checked_on, bare_internal, target)
        if not math.isfinite(ampere_turns) or ampere_turns <= 0.0:
            return
        stored = float(unit.magnetic_energy_j(np.ones(2) * ampere_turns))
        if stored < best_pulsed:
            best_pulsed = stored
            best_pulsed_at = (radius_m, z_m, ampere_turns)

    corner_r, corner_z = _external_corner(wire_m, clearance_m)
    rank(corner_r, corner_z, search, checked_search)
    for radius_m in np.geomspace(wire_m * 2.0, 1.2, 14):
        for z_m in np.geomspace(wire_m * 2.0, 1.2, 14):
            rank(float(radius_m), float(z_m), search, checked_search)
    if best_pulsed_at is None or best_rise_at is None:
        raise RuntimeError("no external pair produces an on-axis field")

    def pulsed_energy(point: np.ndarray) -> float:
        radius_m = float(point[0])
        z_m = float(point[1])
        if radius_m <= wire_m or z_m < wire_m or math.hypot(radius_m, z_m) <= clearance_m:
            return 1.0e30
        try:
            unit, _throat_mid, throat_path = _unit_fields(radius_m, z_m, wire_m, search)
        except ValueError:
            return 1.0e30
        ampere_turns = _pulsed_ampere_turns(throat_path, checked_search, bare_internal, target)
        if not math.isfinite(ampere_turns) or ampere_turns <= 0.0:
            return 1.0e30
        return float(unit.magnetic_energy_j(np.ones(2) * ampere_turns))

    start = np.array([best_pulsed_at[0], best_pulsed_at[1]], dtype=float)
    polished = minimize(
        pulsed_energy,
        start,
        method="Nelder-Mead",
        options={"xatol": 1.0e-5, "fatol": 1.0e-2, "maxiter": 80, "disp": False},
    )
    if math.isfinite(float(polished.fun)) and float(polished.fun) < 1.0e29:
        radius_m = float(polished.x[0])
        z_m = float(polished.x[1])
    else:
        radius_m, z_m = start
    rank(radius_m, z_m, search, checked_search)
    # The reported external store is the fine-path CoilSet energy.
    unit, throat_mid, throat_path = _unit_fields(radius_m, z_m, wire_m, path)
    ampere_turns = _pulsed_ampere_turns(throat_path, checked_path_t, bare_internal, target)
    if not math.isfinite(ampere_turns):
        raise RuntimeError("the external pulsed pair missed the fine path")
    stored = float(unit.magnetic_energy_j(np.ones(2) * ampere_turns))
    rise_r, rise_z = best_rise_at
    _unit, rise_mid, rise_path = _unit_fields(rise_r, rise_z, wire_m, path)
    vacuum_rise = float(np.max(rise_path) / rise_mid)
    asymptotic = vacuum_rise / compression
    steady = _steady_ampere_turns(
        rise_mid, rise_path, checked_mid_t, checked_path_t, compression, target
    )
    return ExternalScore(
        vacuum_rise=vacuum_rise,
        asymptotic_ratio=float(asymptotic),
        radius_m=float(radius_m),
        z_m=float(z_m),
        ampere_turns=float(ampere_turns),
        stored_j=stored,
        pulsed_supply_j=float((1.0 - recovery) * stored),
        steady_reaches=bool(math.isfinite(steady)),
    )


def _checked_wire(problem: Problem) -> tuple[float, float]:
    wires = {float(coil.wire_radius_m) for coil in problem.coils.coils}
    currents = {float(coil.current_a) for coil in problem.coils.coils}
    if len(wires) != 1 or len(currents) != 1:
        raise RuntimeError("the checked-in coils do not share one wire and one current")
    wire = wires.pop()
    current = currents.pop()
    if wire <= 0.0 or current <= 0.0:
        raise RuntimeError("the checked-in wire is not a conductor")
    return wire, current


def build_throat(problem: Problem | None = None) -> ThroatReport:
    problem = load_checked_in() if problem is None else problem
    surface = build_surface(problem)
    shell = build_shell_from_slots(problem.config, problem.slots)
    coils = problem.coils
    base = coils.base_currents_a()
    scales = problem.schedule.current_scales(base.size)
    radii = np.asarray(problem.schedule.plasma_radius_m, dtype=float)
    peak_index = int(np.argmin(radii))
    first = base * scales[0]
    peak = base * scales[peak_index]
    wire, conductor = _checked_wire(problem)
    plasma0 = float(radii[0])
    plasma_peak = float(radii[peak_index])
    wall = float(shell.inner_radius_m)
    compression = (plasma0 / plasma_peak) ** 2
    if compression <= 1.0:
        raise RuntimeError("the schedule does not compress the plasma")
    path = np.linspace(plasma_peak, wall, PATH_SAMPLES)
    checked_mid = float(coils.on_axis_field(0.0, first))
    checked_on_escape = _axis(coils, peak, path)
    bare_internal = checked_mid * compression
    if abs(bare_internal - surface.plasma.field_t) / bare_internal > 1.0e-9:
        raise RuntimeError("throat internal field disagrees with the surface score")
    bare_ratio = float(np.max(checked_on_escape) / bare_internal)
    half = mirror_ratio_for_loss(HALF_LOSS)
    radius_m, z_m, _found = _search_bore(
        wire,
        plasma0,
        wall,
        half,
        path,
        checked_on_escape,
        checked_mid,
        compression,
    )
    # The loop peaks in its own plane. The score samples that plane.
    score_path = _stations(path, z_m, plasma_peak, wall)
    checked_path = _axis(coils, peak, score_path)
    unit, throat_mid, throat_path = _unit_fields(radius_m, z_m, wire, score_path)
    vacuum_rise = float(np.max(throat_path) / throat_mid)
    if vacuum_rise <= half * compression:
        raise RuntimeError("the bore pair's vacuum rise does not clear the compression")
    recovery = float(problem.slots["conversion.coil_energy_recovery_fraction"])
    if not 0.0 <= recovery <= 1.0:
        raise ValueError("coil energy recovery is outside 0 to 1")
    walk = float(surface.beta1_walk_electric_j)
    checked = schedule_yield(
        problem,
        float(problem.slots["plasma.ion_temperature_keV"]),
        fuel_of(problem),
    )
    scales_on_pair = []
    for label, target in (
        ("threshold", 1.0),
        ("half trapped", half),
        ("loss 0.1", mirror_ratio_for_loss(TENTH_LOSS)),
    ):
        solved = _steady_ampere_turns(
            throat_mid, throat_path, checked_mid, checked_path, compression, target
        )
        if not math.isfinite(solved) or solved <= 0.0:
            raise RuntimeError(f"the bore pair has no steady solution for {label}")
        scales_on_pair.append(
            _scale(
                label,
                radius_m,
                z_m,
                wire,
                conductor,
                solved,
                checked_mid,
                checked_path,
                score_path,
                compression,
                recovery,
                walk,
            )
        )
    threshold, design, tenth = scales_on_pair
    if design.mirror_ratio <= 1.0:
        raise RuntimeError("the designed throat is not stronger than the internal field")
    external = _search_external(
        wire,
        float(shell.outer_radius_m) + wire,
        half,
        path,
        checked_on_escape,
        checked_mid,
        compression,
        recovery,
    )
    return ThroatReport(
        wire_m=wire,
        conductor_a=conductor,
        radius_m=float(radius_m),
        z_m=float(z_m),
        compression=float(compression),
        vacuum_rise=vacuum_rise,
        bare_internal_t=float(bare_internal),
        bare_ratio=bare_ratio,
        recovery=recovery,
        walk_electric_j=walk,
        checked_fusion_j=float(checked.fusion_j),
        threshold=threshold,
        design=design,
        tenth=tenth,
        external=external,
    )


def _row(scale: ThroatScale) -> str:
    return (
        f"| {scale.label} | {scale.mirror_ratio:.6g} | {scale.loss_fraction:.6g} | "
        f"{scale.ampere_turns:.6g} | {scale.stored_j:.6g} | {scale.pulsed_supply_j:.6g} | "
        f"{scale.persistent_draw_j:.6g} | {scale.throat_field_t:.6g} | "
        f"{scale.internal_field_t:.6g} | {scale.trapped_electric_j:.6g} |"
    )


def render(problem: Problem | None = None) -> str:
    report = build_throat(problem)
    design = report.design
    external = report.external
    lines = [
        "# Throat",
        "",
        "The gyro-walk remains the cross-field clock. A particle mirrors",
        "when the strongest field on the path to the wall exceeds the field",
        "it occupies. The isotropic share that misses that throat is",
        "`1 - sqrt(1 - 1/R)`. Half of an isotropic distribution is trapped",
        "when the missing share is one half, and the mirror ratio is then 4/3.",
        "",
        f"The frozen internal field on the checked-in pulse is `{_num(report.bare_internal_t)}` T.",
        f"The checked-in coils at peak scale reach `{_num(report.bare_ratio)}` of that field.",
        "The designed pair is a steady bore winding on the checked-in wire.",
        "The checked-in coil list is unchanged.",
        "",
        f"The centerline radius is `{_num(report.radius_m)}` m and the pair sits at",
        f"`±{_num(report.z_m)}` m. That centerline is the smallest radius the",
        f"checked-in wire, `{_num(report.wire_m)}` m, allows, and the conductor lies",
        "against the inner wall. The bore pair is stronger than the internal field.",
        f"Its vacuum rise is `{_num(report.vacuum_rise)}`. The compression is",
        f"`{_num(report.compression)}`. The pair is present at freeze-in.",
        "",
        "The threshold row is mirror ratio 1 on the same pair. "
        "A path whose strongest field is the occupied field has no trapped set. "
        "The half-trapped row is the design. The one-tenth loss row is the same pair "
        "at the mirror ratio where one tenth of an isotropic distribution is in the loss cone.",
        "",
        "| Scale | Mirror ratio | Loss | Ampere-turns | Stored (J) | Pulsed supply (J) | Persistent draw (J) | Throat (T) | Internal (T) | Trapped electricity (J) |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        _row(report.threshold),
        _row(design),
        _row(report.tenth),
        "",
        (
            f"The designed conductor is `{_num(design.turns)}` "
            f"{'turn' if design.turns == 1.0 else 'turns'} at `{_num(design.current_a)}` A,"
        ),
        f"on the checked-in conductor class of `{_num(report.conductor_a)}` A.",
        f"The throat field is `{_num(design.throat_field_t)}` T and the internal field",
        f"on this current is `{_num(design.internal_field_t)}` T.",
        "",
        "A persistent throat current draws no magnetic energy from the supply.",
        f"The pulsed supply at recovery `{_num(report.recovery)}` is",
        f"`{_num(design.pulsed_supply_j)}` J.",
        "Trapped particles stay on the gyro-walk.",
        "Trapped electricity versus the gyro-walk electricity: "
        f"{_side(design.trapped_electric_j, report.walk_electric_j)}.",
        "Trapped electricity versus the pulsed supply: "
        f"{_side(design.trapped_electric_j, design.pulsed_supply_j)}.",
        "",
        "A steady pair outside the shell stays below the internal field.",
        f"The strongest vacuum rise outside the shell is `{_num(external.vacuum_rise)}`.",
        "At large current that pair's mirror ratio approaches the vacuum rise",
        f"divided by the compression, `{_num(external.asymptotic_ratio)}`.",
        f"A pulsed pair outside the shell, radius `{_num(external.radius_m)}` m at",
        f"`±{_num(external.z_m)}` m, stores `{_num(external.stored_j)}` J to reach",
        "the half-trapped ratio.",
        "Bore stored energy versus that pulsed store: "
        f"{_side(design.stored_j, external.stored_j)}.",
        "",
        f"Checked-in fusion remains `{_num(report.checked_fusion_j)}` J.",
        "The checked-in schedule is not rewritten.",
        "",
    ]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    del argv
    text = render()
    output = default_output_path()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text, encoding="utf-8")
    print(text, end="" if text.endswith("\n") else "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
