"""Call the live FSOT 2.1 engine and run the SolCascade coil ledgers.

The engine stays in FSOT-2.1-Lean. This file does not copy vendor/fsot_compute.py.
Pin AEB2AD is the SHA-256 prefix of that file. The C++ port states the same
K and beta formulas in include/fsot/core.hpp.

The fusion panels in the Lean repo rescale a published number,
measured * (1 + |S| * epsilon). Those products are not written into slots.
Domain scalars are dimensionless. They are recorded here and are not
converted into keV, tesla, or amperes.

What does change the ledger is the coil current that produces the
frozen-flux field. Under frozen_axial_flux the internal field is fixed by
the current at the first time knot. Later current raises stored energy and
does not raise that field. The minimum of (1/2) I^T L I at that on-axis
field is the current distribution this script runs.

Coil location is scored on the shell clearance sphere. Engine numbers are
candidate values of |z| / cylindrical radius. A midplane angle scan is the
Maxwell minimum on that same sphere. Neither result is written back into
the engine. The ledger confinement time is one thermal gyroradius per cyclotron
period. Bohm's 16 is the comparison. The domain scalar is the orifice
witness and is not a second coefficient.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
ENGINE_PATH = Path(r"C:\Users\damia\Desktop\FSOT-2.1-Lean\vendor\fsot_compute.py")
CPP_CORE = Path(r"C:\Users\damia\Desktop\FSOT-2.1-Cpp\include\fsot\core.hpp")
PIN_PREFIX = "aeb2adad6e80f487"
CONFIG = ROOT / "configs" / "minimal_pipeline.yaml"
MATERIALS = ROOT / "configs" / "materials.yaml"
OUT = ROOT / "runs" / "fsot-apply"
TURNS = 12.0
WIRE_M = 0.008
CLEARANCE_M = 1.0e-4

sys.path.insert(0, str(ROOT / "src"))

from solcascade.assemble import build_coils, build_shell_from_slots  # noqa: E402
from solcascade.constants import (  # noqa: E402
    ELEMENTARY_CHARGE_C,
    JOULE_PER_EV,
    MASS_EXCESS_DEUTERIUM_KEV,
    MASS_EXCESS_TRITIUM_KEV,
    ion_mass_kg,
)
from solcascade.geometry.coils import Coil, CoilSet  # noqa: E402
from solcascade.interfaces import NullProvider, apply_provider, collect_slots  # noqa: E402
from solcascade.load import load_problem, load_yaml  # noqa: E402
from solcascade.materials.library import load_materials  # noqa: E402
from solcascade.pipeline import run_pipeline  # noqa: E402
from solcascade.plasma.metrics import bohm_time_s, gyro_walk_time_s  # noqa: E402


def load_engine():
    digest = hashlib.sha256(ENGINE_PATH.read_bytes()).hexdigest()
    if not digest.startswith(PIN_PREFIX):
        raise SystemExit(f"engine pin is {digest[:16]}, expected {PIN_PREFIX}")
    spec = importlib.util.spec_from_file_location("fsot_compute_live", ENGINE_PATH)
    if spec is None or spec.loader is None:
        raise SystemExit(f"cannot load {ENGINE_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module, digest


def engine_record(module, digest: str) -> dict:
    names = (
        "Particle_Physics",
        "High_Energy_Physics",
        "Electromagnetism",
        "Nuclear_Physics",
        "Thermodynamics",
        "Condensed_Matter",
    )
    scalars = {name: float(module.domain_scalar(name)) for name in names}
    cpp = CPP_CORE.read_text(encoding="utf-8")
    k_line = "PHI * (GAMMA / E) * sqrt(R(2)) / log(PI) * (R(1) - R(1) / ipow(PI, 4))"
    beta_line = "R(1) / exp(pow(PI, PI) + (E - R(1)))"
    return {
        "pin_prefix": digest[:16],
        "sha256": digest,
        "engine": str(ENGINE_PATH),
        "cpp_core": str(CPP_CORE),
        "cpp_K_formula_present": k_line in cpp,
        "cpp_beta_formula_present": beta_line in cpp,
        "K": float(module.K),
        "beta": float(module.BETA),
        "phi": float(module.PHI),
        "phi_to_minus_4": float(module.PHI ** (-4)),
        "edge_theta_C_eff_P_var": float(module.C_EFF * module.P_VAR),
        "domain_scalar": scalars,
        "nuclear_equals_thermodynamics": scalars["Nuclear_Physics"] == scalars["Thermodynamics"],
        "not_injected": (
            "Fusion-panel rescales and these dimensionless scalars are not "
            "slot values. Bare engine ALPHA is not installed as mu_0."
        ),
    }


def field_per_amp(coils: CoilSet) -> np.ndarray:
    count = len(coils.coils)
    weights = np.zeros(count, dtype=float)
    for index in range(count):
        current = np.zeros(count, dtype=float)
        current[index] = 1.0
        weights[index] = coils.on_axis_field(0.0, current)
    return weights


def minimum_energy_currents(coils: CoilSet, field_t: float) -> tuple[np.ndarray, float]:
    """Currents minimizing (1/2) I^T L I at a fixed on-axis field."""
    inductance = coils.inductance_matrix_h()
    weights = field_per_amp(coils)
    solved = np.linalg.solve(inductance, weights)
    denominator = float(weights @ solved)
    if denominator == 0.0:
        raise SystemExit("coils produce no on-axis field")
    lagrange = field_t / denominator
    currents = lagrange * solved
    energy = 0.5 * float(currents @ inductance @ currents)
    return currents, energy


def write_case(
    name: str,
    raw: dict,
    currents: dict[str, float],
    scales: list[float],
    coil_scales: dict[str, list[float]] | None = None,
) -> Path:
    case = json.loads(json.dumps(raw))
    case["name"] = name
    case["materials_file"] = str(MATERIALS)
    case["transport"]["execute"] = False
    case["compression"]["coil_current_scale"] = scales
    if coil_scales is not None:
        case["compression"]["coil_scales"] = coil_scales
    by_name = {coil["name"]: coil for coil in case["geometry"]["coils"]}
    for coil_name, current in currents.items():
        by_name[coil_name]["current_a"] = float(current)
    path = OUT / "configs" / f"{name}.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(case, sort_keys=False), encoding="utf-8")
    return path


def helmholtz_config(raw: dict, field_t: float) -> tuple[dict, np.ndarray, float]:
    """Two coils just outside the shell, spaced at half a radius."""
    shell = raw["geometry"]["shell"]
    outer = float(shell["inner_radius_m"]) + sum(
        float(layer["thickness_m"]) for layer in shell["layers"]
    )
    radius = outer + 0.021
    spacing = 0.5 * radius
    spec = {
        "turns": 12,
        "wire_radius_m": 0.008,
        "material": "copper",
    }
    coils = CoilSet(
        (
            Coil("helmholtz_minus", radius, -spacing, spec["turns"], 1.0, spec["wire_radius_m"], "copper"),
            Coil("helmholtz_plus", radius, spacing, spec["turns"], 1.0, spec["wire_radius_m"], "copper"),
        )
    )
    currents, energy = minimum_energy_currents(coils, field_t)
    case = json.loads(json.dumps(raw))
    case["geometry"]["coils"] = [
        {
            "name": "helmholtz_minus",
            "radius_m": radius,
            "z_m": -spacing,
            "current_a": float(currents[0]),
            **spec,
        },
        {
            "name": "helmholtz_plus",
            "radius_m": radius,
            "z_m": spacing,
            "current_a": float(currents[1]),
            **spec,
        },
    ]
    return case, currents, energy


def clearance_winding(sphere_m: float, aspect: float, turns: float, wire_m: float) -> CoilSet:
    """Coils on the clearance sphere. aspect is |z| / cylindrical radius."""
    if aspect < 0.0:
        raise ValueError("aspect must be >= 0")
    if aspect == 0.0:
        return CoilSet((Coil("midplane", sphere_m, 0.0, turns, 1.0, wire_m, "copper"),))
    rho = sphere_m / math.sqrt(1.0 + aspect * aspect)
    axial = aspect * rho
    return CoilSet(
        (
            Coil("pair_minus", rho, -axial, turns, 1.0, wire_m, "copper"),
            Coil("pair_plus", rho, axial, turns, 1.0, wire_m, "copper"),
        )
    )


def score_winding(
    winding: CoilSet,
    field_t: float,
    resistivity: float,
    duration_s: float,
    recovery: float,
) -> dict:
    currents, stored = minimum_energy_currents(winding, field_t)
    ohmic_w = winding.ohmic_power_w(currents, np.full(currents.size, resistivity))
    ohmic_j = ohmic_w * duration_s
    practical = ohmic_j + (1.0 - recovery) * stored
    achieved = winding.on_axis_field(0.0, currents)
    return {
        "currents_A": [float(current) for current in currents],
        "radius_m": [coil.radius_m for coil in winding.coils],
        "z_m": [coil.z_m for coil in winding.coils],
        "stored_J": stored,
        "ohmic_J": ohmic_j,
        "practical_J": practical,
        "achieved_field_T": achieved,
        "field_error": abs(achieved - field_t) / field_t,
    }


def geometry_config(raw: dict, score: dict, names: list[str]) -> dict:
    case = json.loads(json.dumps(raw))
    case["geometry"]["coils"] = [
        {
            "name": name,
            "radius_m": radius,
            "z_m": axial,
            "turns": TURNS,
            "current_a": current,
            "wire_radius_m": WIRE_M,
            "material": "copper",
        }
        for name, radius, axial, current in zip(
            names, score["radius_m"], score["z_m"], score["currents_A"], strict=True
        )
    ]
    return case


def confinement_check(module, radius_m: float, temperature_kev: float, field_t: float) -> dict:
    """Ledger time is the gyro-period walk. Bohm's 16 is the comparison."""
    bohm = bohm_time_s(radius_m, temperature_kev, field_t)
    walk = gyro_walk_time_s(radius_m, temperature_kev, field_t)
    mass = 0.5 * (
        ion_mass_kg(2, MASS_EXCESS_DEUTERIUM_KEV, 1)
        + ion_mass_kg(3, MASS_EXCESS_TRITIUM_KEV, 1)
    )
    thermal = temperature_kev * JOULE_PER_EV * 1.0e3
    speed = math.sqrt(thermal / mass)
    rho = mass * speed / (ELEMENTARY_CHARGE_C * field_t)
    # Transit-time clock, χ = ρ² v / a. Same a²/(2D) definition. Not installed.
    chi_gyro = rho * rho * speed / radius_m
    gyro = radius_m * radius_m / (2.0 * chi_gyro)
    parent_s = float(module.domain_scalar("Condensed_Matter"))
    rest_k = float(module.K)
    return {
        "radius_m": radius_m,
        "temperature_keV": temperature_kev,
        "field_T": field_t,
        "law": "D = kT / (4 π e B), one thermal gyroradius per cyclotron period",
        "walk_time_s": walk,
        "bohm_time_s": bohm,
        "bohm_denominator": 16.0,
        "walk_denominator": 4.0 * math.pi,
        "walk_over_bohm": walk / bohm,
        "plasma_parent_orifice": "Condensed_Matter",
        "parent_scalar_S": parent_s,
        "rest_unit_K": rest_k,
        "gamma_fsot": parent_s / rest_k,
        "gamma_enters_the_time": False,
        "ion_mass_kg": mass,
        "thermal_speed_m_s": speed,
        "ion_gyroradius_m": rho,
        "radius_over_gyroradius": radius_m / rho,
        "transit_clock_time_s": gyro,
        "transit_clock_over_walk": gyro / walk,
        "installed_time": "gyro_period_random_walk",
        "note": (
            "The step is the thermal gyroradius and the clock is one cyclotron "
            "period. D = ρ² / (2 τ_c) is kT / (4 π e B). The historical 16 is "
            "the coefficient this replaces. The walk is π/4 of that Bohm time "
            "at every state, so n T τ still scales with n and B. The "
            "Condensed_Matter scalar is the orifice witness. It is not a "
            "second coefficient, and β is not a time. A transit-time clock "
            "would let temperature remain in the triple product. That clock "
            "is reported and is not the ledger time. Fusion joules follow the "
            "prescribed schedule."
        ),
    }


def _legal_pair(winding: CoilSet) -> bool:
    if len(winding.coils) == 2 and abs(winding.coils[0].z_m) <= WIRE_M:
        return False
    return True


def score_aspect(
    sphere_m: float,
    aspect: float,
    field_t: float,
    resistivity: float,
    duration_s: float,
    recovery: float,
    label: str,
    origin: str,
) -> dict:
    """One clearance-sphere winding. Overlap and a non-physical bore are skipped."""
    try:
        winding = clearance_winding(sphere_m, aspect, TURNS, WIRE_M)
    except ValueError as exc:
        return {
            "label": label,
            "origin": origin,
            "aspect_abs_z_over_rho": aspect,
            "skipped": str(exc),
        }
    if not _legal_pair(winding):
        return {
            "label": label,
            "origin": origin,
            "aspect_abs_z_over_rho": aspect,
            "skipped": "coil bundles overlap",
        }
    score = score_winding(winding, field_t, resistivity, duration_s, recovery)
    if score["field_error"] > 1.0e-8:
        raise SystemExit(f"{label} missed the on-axis field by {score['field_error']}")
    return {
        "label": label,
        "origin": origin,
        "aspect_abs_z_over_rho": aspect,
        "n_coils": len(winding.coils),
        "names": [coil.name for coil in winding.coils],
        "skipped": None,
        **score,
    }


def _same_geometry(left: dict, right: dict) -> bool:
    if left.get("n_coils") != right.get("n_coils"):
        return False
    left_points = sorted(zip(left["radius_m"], left["z_m"], strict=True))
    right_points = sorted(zip(right["radius_m"], right["z_m"], strict=True))
    return all(
        abs(radius_a - radius_b) <= 1.0e-9 and abs(axial_a - axial_b) <= 1.0e-9
        for (radius_a, axial_a), (radius_b, axial_b) in zip(left_points, right_points, strict=True)
    )


def _with_practical(row: dict, energy_j: float, recovery: float) -> dict:
    """Ledger draw plus the unrecovered part of the energy already stored at t=0."""
    extra = (1.0 - recovery) * energy_j
    total = row["supply_drawn_J"] + extra
    enriched = dict(row)
    enriched["energy_at_first_knot_J"] = energy_j
    enriched["unrecovered_first_knot_J"] = extra
    enriched["practical_supply_J"] = total
    enriched["fusion_over_practical_supply"] = row["fusion_J"] / total
    return enriched


def _check_reference(row: dict) -> None:
    # Bohm-era ratio. The walk is π/4 of that time, and the pulse cap does not bind.
    bohm_ratio = 3.091993244213298e-06
    walk_ratio = bohm_ratio * math.pi / 4.0
    expected = {
        "baseline": (0.4977033387900156, 1423874.990828079, walk_ratio),
        "same_b_flat": (0.4977033387900156, 2106.2160000000003, walk_ratio),
        "same_b_optimal": (0.4977033387900156, 2020.4381637452752, walk_ratio),
        "same_b_helmholtz": (0.4977033387900156, 827.5009597961755, walk_ratio),
    }
    if row["name"] not in expected:
        return
    fusion, supply, ratio = expected[row["name"]]
    if not math.isclose(row["fusion_J"], fusion, rel_tol=0.0, abs_tol=1.0e-12):
        raise SystemExit(f"{row['name']} fusion moved to {row['fusion_J']}")
    if not math.isclose(row["supply_drawn_J"], supply, rel_tol=0.0, abs_tol=1.0e-6):
        raise SystemExit(f"{row['name']} supply moved to {row['supply_drawn_J']}")
    if not math.isclose(row["peak_triple_product_ratio"], ratio, rel_tol=1.0e-9, abs_tol=0.0):
        raise SystemExit(f"{row['name']} triple-product ratio moved to {row['peak_triple_product_ratio']}")


def _staged_schedule(
    coils: CoilSet,
    field_t: float,
    optimal_map: dict[str, float],
    shared: list[float],
) -> tuple[dict[str, float], dict[str, list[float]], float, float]:
    """Inners hold the confining field. Outers turn on only at the flat top."""
    inners = [coil for coil in coils.coils if coil.name.startswith("mirror_inner")]
    outers = [coil for coil in coils.coils if coil.name.startswith("mirror_outer")]
    if len(inners) != 2 or len(outers) != 2:
        raise SystemExit("staged case expects two inner and two outer mirror coils")
    inner_currents, _stored = minimum_energy_currents(CoilSet(tuple(inners)), field_t)
    currents = {coil.name: float(current) for coil, current in zip(inners, inner_currents, strict=True)}
    for coil in outers:
        currents[coil.name] = float(optimal_map[coil.name])
    peak = max(shared)
    outer_scale = [1.0 if value == peak else 0.0 for value in shared]
    held = [1.0] * len(shared)
    coil_scales = {
        coil.name: list(held) if coil.name.startswith("mirror_inner") else list(outer_scale)
        for coil in coils.coils
    }
    vector = np.array(
        [
            currents[coil.name] * (1.0 if coil.name.startswith("mirror_inner") else 0.0)
            for coil in coils.coils
        ]
    )
    return currents, coil_scales, coils.magnetic_energy_j(vector), coils.on_axis_field(0.0, vector)


def run_case(name: str, config_path: Path, record: dict) -> dict:
    workdir = OUT / name
    summary = run_pipeline(config_path, workdir=workdir)
    (workdir / "fsot_engine.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    conversion = summary["conversion"]
    plasma = summary["plasma"]
    return {
        "name": name,
        "fusion_J": conversion["fusion_J"],
        "ohmic_J": conversion["ohmic_J"],
        "magnetic_unrecovered_J": conversion["magnetic_unrecovered_J"],
        "supply_drawn_J": conversion["supply_drawn_J"],
        "net_electrical_J": conversion["net_electrical_J"],
        "peak_ion_temperature_keV": plasma["peak_ion_temperature_keV"],
        "peak_internal_field_T": plasma["peak_internal_field_T"],
        "peak_vacuum_field_T": plasma["peak_vacuum_field_T"],
        "peak_triple_product_ratio": plasma["peak_triple_product_ratio"],
        "fusion_over_supply": conversion["fusion_J"] / conversion["supply_drawn_J"],
    }


def main() -> int:
    module, digest = load_engine()
    record = engine_record(module, digest)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "fsot_engine.json").write_text(json.dumps(record, indent=2), encoding="utf-8")

    raw = load_yaml(CONFIG)
    config, materials = load_problem(CONFIG)
    slots, _provenance = apply_provider(collect_slots(config, materials), NullProvider())
    coils = build_coils(config, slots)
    base = coils.base_currents_a()
    initial_scale = float(raw["compression"]["coil_current_scale"][0])
    confining = base * initial_scale
    field_t = coils.on_axis_field(0.0, confining)
    inductance = coils.inductance_matrix_h()
    energy_confining = coils.magnetic_energy_j(confining, inductance)
    energy_peak = coils.magnetic_energy_j(base, inductance)
    optimal, energy_optimal = minimum_energy_currents(coils, field_t)
    names = [coil.name for coil in coils.coils]
    optimal_map = {name: float(current) for name, current in zip(names, optimal, strict=True)}

    helm_raw, _helm_currents, helm_energy = helmholtz_config(raw, field_t)
    flat = [initial_scale] * len(raw["compression"]["coil_current_scale"])
    ones = [1.0] * len(flat)

    cases = [
        ("baseline", raw, {name: float(current) for name, current in zip(names, base, strict=True)}, list(raw["compression"]["coil_current_scale"])),
        ("same_b_flat", raw, {name: float(current) for name, current in zip(names, base, strict=True)}, flat),
        ("same_b_optimal", raw, optimal_map, ones),
        ("same_b_helmholtz", helm_raw, {coil["name"]: float(coil["current_a"]) for coil in helm_raw["geometry"]["coils"]}, ones),
    ]
    rows = []
    for name, body, currents, scales in cases:
        path = write_case(name, body, currents, scales)
        row = run_case(name, path, record)
        _check_reference(row)
        rows.append(row)

    recovery = float(raw["conversion"]["coil_energy_recovery_fraction"])
    duration = float(raw["compression"]["time_s"][-1] - raw["compression"]["time_s"][0])
    library = load_materials(materials, slots)
    resistivity = library["copper"].resistivity_at(
        float(raw["conversion"]["conductor_temperature_K"])
    )
    shell = build_shell_from_slots(config, slots)
    sphere = shell.outer_radius_m + WIRE_M + CLEARANCE_M
    phi = float(module.PHI)
    seed_specs = (
        ("midplane", "clearance_sphere_midplane", 0.0),
        ("classical_half", "classical_helmholtz_aspect", 0.5),
        ("phi_inverse", "engine_number_as_aspect", 1.0 / phi),
        ("phi_to_minus_2", "engine_number_as_aspect", phi ** -2),
        ("phi_to_minus_4", "engine_number_as_aspect", phi ** -4),
        (
            "S_electromagnetism",
            "engine_number_as_aspect",
            float(module.domain_scalar("Electromagnetism")),
        ),
        (
            "S_nuclear",
            "engine_number_as_aspect",
            float(module.domain_scalar("Nuclear_Physics")),
        ),
        ("edge_theta", "engine_number_as_aspect", float(module.C_EFF * module.P_VAR)),
    )
    seeds = [
        score_aspect(sphere, aspect, field_t, resistivity, duration, recovery, label, origin)
        for label, origin, aspect in seed_specs
    ]
    scan = [
        score_aspect(
            sphere,
            0.0 if degree == 0 else math.tan(math.radians(degree)),
            field_t,
            resistivity,
            duration,
            recovery,
            f"scan_{degree:02d}deg",
            "maxwell_scan",
        )
        for degree in range(0, 76)
    ]
    runnable = [item for item in seeds if item.get("skipped") is None]
    scanned = [item for item in scan if item.get("skipped") is None]
    if not runnable or not scanned:
        raise SystemExit("clearance sphere produced no legal winding")
    best_seed = min(runnable, key=lambda item: item["practical_J"])
    best_scan = min(scanned, key=lambda item: item["practical_J"])

    def practical_at_degrees(degree: float) -> float:
        aspect = 0.0 if degree <= 0.0 else math.tan(math.radians(degree))
        scored = score_aspect(
            sphere, aspect, field_t, resistivity, duration, recovery,
            "stationary", "maxwell_stationary",
        )
        if scored.get("skipped") is not None:
            return math.inf
        return float(scored["practical_J"])

    # Pairs start once the bundles clear the wire. The midplane loop is scored separately.
    low, high = 2.0, 70.0
    best_degree = low
    best_cost = practical_at_degrees(low)
    degree = low
    while degree <= high:
        cost = practical_at_degrees(degree)
        if cost < best_cost:
            best_degree = degree
            best_cost = cost
        degree += 0.5
    left = max(low, best_degree - 0.5)
    right = min(high, best_degree + 0.5)
    for _ in range(40):
        span = right - left
        if span < 1.0e-4:
            break
        first = left + span / 3.0
        second = right - span / 3.0
        if practical_at_degrees(first) < practical_at_degrees(second):
            right = second
        else:
            left = first
    stationary_degree = 0.5 * (left + right)
    stationary = score_aspect(
        sphere,
        math.tan(math.radians(stationary_degree)),
        field_t,
        resistivity,
        duration,
        recovery,
        f"stationary_{stationary_degree:.4f}deg",
        "maxwell_stationary",
    )
    best_pair = min(
        (item for item in scanned if item["n_coils"] == 2),
        key=lambda item: item["practical_J"],
        default=None,
    )
    first_knot = {
        "baseline": energy_confining,
        "same_b_flat": energy_confining,
        "same_b_optimal": energy_optimal,
        "same_b_helmholtz": helm_energy,
    }
    ran = []
    for chosen in (best_seed, best_scan, stationary):
        if any(_same_geometry(chosen, earlier) for earlier in ran):
            continue
        body = geometry_config(raw, chosen, chosen["names"])
        currents = {
            name: current
            for name, current in zip(chosen["names"], chosen["currents_A"], strict=True)
        }
        case_name = "topology_" + chosen["label"]
        path = write_case(case_name, body, currents, ones)
        row = run_case(case_name, path, record)
        ohmic_gap = abs(row["ohmic_J"] - chosen["ohmic_J"]) / chosen["ohmic_J"]
        if ohmic_gap > 1.0e-6:
            raise SystemExit(f"{case_name} ledger ohmic disagrees with the score by {ohmic_gap}")
        if not math.isclose(row["fusion_J"], rows[0]["fusion_J"], rel_tol=0.0, abs_tol=1.0e-12):
            raise SystemExit(f"{case_name} moved the fusion yield")
        rows.append(row)
        first_knot[case_name] = chosen["stored_J"]
        ran.append(chosen)

    staged_currents, coil_scales, staged_energy, staged_field = _staged_schedule(
        coils, field_t, optimal_map, list(raw["compression"]["coil_current_scale"])
    )
    if abs(staged_field - field_t) / field_t > 1.0e-8:
        raise SystemExit(f"staged first knot field is {staged_field}, target {field_t}")
    staged_path = write_case(
        "staged_outers_after_knot", raw, staged_currents, ones, coil_scales
    )
    staged_row = run_case("staged_outers_after_knot", staged_path, record)
    if not math.isclose(
        staged_row["peak_internal_field_T"],
        rows[0]["peak_internal_field_T"],
        rel_tol=0.0,
        abs_tol=1.0e-9,
    ):
        raise SystemExit("staged waveform moved the frozen internal field")
    if not math.isclose(staged_row["fusion_J"], rows[0]["fusion_J"], rel_tol=0.0, abs_tol=1.0e-12):
        raise SystemExit("staged waveform moved the fusion yield")
    rows.append(staged_row)
    first_knot["staged_outers_after_knot"] = staged_energy

    rows = [_with_practical(row, first_knot[row["name"]], recovery) for row in rows]
    baseline = rows[0]
    confinement = confinement_check(
        module,
        radius_m=min(float(radius) for radius in raw["compression"]["plasma_radius_m"]),
        temperature_kev=baseline["peak_ion_temperature_keV"],
        field_t=baseline["peak_internal_field_T"],
    )
    report = {
        "confining_on_axis_field_T": field_t,
        "initial_scale": initial_scale,
        "stored_J_at_confining_current": energy_confining,
        "stored_J_at_yaml_peak_current": energy_peak,
        "stored_J_optimal_same_field": energy_optimal,
        "stored_J_helmholtz_same_field": helm_energy,
        "optimal_current_A": optimal_map,
        "unrecovered_if_confining_field_is_later_dumped_J": (1.0 - recovery) * energy_confining,
        "unrecovered_if_optimal_field_is_later_dumped_J": (1.0 - recovery) * energy_optimal,
        "ledger": rows,
        "topology": {
            "objective": (
                "minimum of ohmic*duration + (1-recovery)*stored energy "
                "at the baseline on-axis field"
            ),
            "sphere_radius_m": sphere,
            "shell_outer_radius_m": shell.outer_radius_m,
            "wire_radius_m": WIRE_M,
            "turns": TURNS,
            "duration_s": duration,
            "recovery": recovery,
            "resistivity_ohm_m": resistivity,
            "aspect_definition": "|z| / cylindrical radius, coils on the clearance sphere",
            "seed_family": seeds,
            "maxwell_scan_minimum": best_scan,
            "stationary_minimum": stationary,
            "stationary_degrees_from_midplane": stationary_degree,
            "best_pair_on_scan": best_pair,
            "pipeline_cases": [item["label"] for item in ran],
            "scan_matches_seed_winner": _same_geometry(best_seed, best_scan),
            "fitted_to_ledger": False,
            "written_back_into_engine": False,
            "note": (
                "Engine numbers are candidate aspects. They are not a fit to "
                "15 kJ or 47 kJ. The degree scan is the Maxwell minimum on "
                "this sphere. Conversion channels stay idle."
            ),
        },
        "staged": {
            "name": "staged_outers_after_knot",
            "inner_current_A": {
                name: current
                for name, current in staged_currents.items()
                if name.startswith("mirror_inner")
            },
            "outer_current_A": {
                name: current
                for name, current in staged_currents.items()
                if name.startswith("mirror_outer")
            },
            "outer_scale": coil_scales["mirror_outer_minus"],
            "first_knot_field_T": staged_field,
            "energy_at_first_knot_J": staged_energy,
            "note": (
                "Independent scales let the outers turn on after the first knot. "
                "Under frozen_axial_flux that adds supply and does not raise "
                "the internal field."
            ),
        },
        "confinement": confinement,
        "confinement_law": (
            "The ledger time is one thermal gyroradius per cyclotron period, "
            "D = kT / (4 π e B). That is π/4 of the historical Bohm time, so "
            "n T τ still scales with n and B. The Condensed_Matter scalar "
            "witnesses the orifice and does not multiply the time."
        ),
    }
    (OUT / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    compact = {
        "confining_on_axis_field_T": field_t,
        "ledger": [
            {
                "name": row["name"],
                "fusion_J": row["fusion_J"],
                "supply_drawn_J": row["supply_drawn_J"],
                "practical_supply_J": row["practical_supply_J"],
                "fusion_over_supply": row["fusion_over_supply"],
                "fusion_over_practical_supply": row["fusion_over_practical_supply"],
                "peak_triple_product_ratio": row["peak_triple_product_ratio"],
            }
            for row in rows
        ],
        "topology_pipeline": [item["label"] for item in ran],
        "topology_best_seed_J": best_seed["practical_J"],
        "topology_best_scan_J": best_scan["practical_J"],
        "installed_time": confinement["installed_time"],
        "walk_time_s": confinement["walk_time_s"],
        "walk_over_bohm": confinement["walk_over_bohm"],
        "stationary_degrees": stationary_degree,
        "radius_over_gyroradius": confinement["radius_over_gyroradius"],
    }
    print(json.dumps(compact, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
