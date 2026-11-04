"""Dimensionful derivations for the energy-baseline ranking.

The confinement clock is ``gyro_walk_time_s``. The joule map holds the
SolCascade peak state for that one time:

    E = n_D n_T <σv>(T) Q V τ

Density and temperature follow the same compression as ``plasma.metrics``.
``gamma_FSOT`` stays out of both expressions. A catalog row is replaced
only where this yield is the same reaction the file scores.

The open catalog units take quantities the same peak already has. The
radius inside the walk is the length. The ion thermal speed inside that
step is the speed. That same ``kT`` is the electron-volt. Each factor uses
the power the phase-two reading already records. D-He3 uses the Bosch-Hale
fit on its published interval. p-B11 uses the Nevins-Swain low-temperature
branch on 1 keV to 70 keV.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import yaml

from solcascade.constants import (
    BOLTZMANN_J_K,
    EV_PER_JOULE,
    JOULE_PER_EV,
    MASS_EXCESS_DEUTERIUM_KEV,
    MASS_EXCESS_TRITIUM_KEV,
    MONATOMIC_ADIABATIC_INDEX,
    REDUCED_PLANCK_J_S,
    SPEED_OF_LIGHT_M_S,
    dt_reaction,
    ion_mass_kg,
)
from solcascade.energy_baseline.rank import load_records
from solcascade.plasma.metrics import (
    adiabatic_temperature_exponent,
    bohm_time_s,
    gyro_walk_time_s,
)
from solcascade.plasma.reactivity import (
    dhe3_reactivity_m3_s,
    dt_reactivity_m3_s,
    pb11_reactivity_m3_s,
)
from solcascade.units import EV_PER_KEV

FUSION_IDS = (
    "fusion-dt-magnetic",
    "fusion-dt-inertial",
    "fusion-dt-pulsed-power",
    "fusion-d-he3",
    "fusion-aneutronic-pb11",
)
DT_YIELD_IDS = (
    "fusion-dt-magnetic",
    "fusion-dt-inertial",
    "fusion-dt-pulsed-power",
)

FORMULA_ID = "gyro_walk_time_s"
LAW = "D = kT / (4 π e B), τ = a² / (2 D) = 2 π a² e B / kT"
SOURCE = "src/solcascade/plasma/metrics.py gyro_walk_time_s"
INPUTS_SOURCE = "runs/fsot-apply/report.json confinement"
YIELD_FORMULA_ID = "confined_dt_yield_j"
YIELD_LAW = "E = n_D n_T <σv>(T) Q V τ"
YIELD_SOURCE = (
    "src/solcascade/plasma/metrics.py power density times gyro_walk_time_s"
)
ORDER_OF_MAGNITUDE = 10.0
RESULT_WHEN_NONE = (
    "No derived substitution moves any coverage by an order of magnitude."
)
NEXT_TARGET = (
    "Next derivation target: the boundary-drive frequency of casimir-dynamical. "
    "Breakdown field, acoustic drive, and radiator temperature difference "
    "remain unmapped. The static Casimir coverage stays 0 because the "
    "demonstrated efficiency is 0. The Nevins and Swain low-temperature "
    "p-B11 reactivity is installed on 1 keV to 70 keV. The high-temperature "
    "branch of that paper starts at 50 keV and is not evaluated at this peak."
)

# Catalog inversion constants. They recover the head, the range, and the
# wind speed from the equations the YAML already states. They are not a fit
# to the shell. g is the 9.81 m/s^2 in the hydroelectric governing equation.
# Water density 1000 kg/m^3 is the density in (1/2) rho g R^2. Air density
# 1.225 kg/m^3 is ISA sea-level air in (1/2) rho v^3. The photovoltaic
# ceiling paper, Shockley and Queisser, J. Appl. Phys. 32, 510 (1961),
# takes the sun as a 6000 K blackbody. The Casimir reference unit is
# 10 cm^2 at 1 micrometer; the file stores that energy rounded.
CATALOG_GRAVITY_M_S2 = 9.81
CATALOG_WATER_DENSITY_KG_M3 = 1000.0
CATALOG_AIR_DENSITY_KG_M3 = 1.225
SOLAR_BLACKBODY_K = 6000.0
CASIMIR_STATED_AREA_M2 = 1.0e-3
CASIMIR_STATED_GAP_M = 1.0e-6

HEAD_FORMULA_ID = "confinement_radius_over_head"
HEAD_LAW = "E ∝ h, factor = a / h, with h = energy_density / g"
HEAD_SOURCE = "hydroelectric.yaml potential energy m g h; a is gyro_walk radius"
TIDAL_FORMULA_ID = "confinement_radius_over_range_squared"
TIDAL_LAW = "E ∝ R², factor = (a / R)², with (1/2) ρ g R² = energy_density"
TIDAL_SOURCE = "tidal.yaml basin potential; a is gyro_walk radius"
WIND_FORMULA_ID = "thermal_speed_over_wind_cubed"
WIND_LAW = "P ∝ v³, factor = (v_thermal / v_catalog)³, with (1/2) ρ v³ = power / Cp"
WIND_SOURCE = "wind.yaml Betz flux; v_thermal is the ion speed inside the walk"
PV_FORMULA_ID = "pulse_kt_detailed_balance"
PV_LAW = (
    "fraction of 6000 K blackbody energy above E = kT; "
    "a zero tail gives coverage factor 0"
)
PV_SOURCE = "photovoltaic.yaml Shockley-Queisser ceiling; E is the pulse kT"
GAP_FORMULA_ID = "confinement_radius_casimir_gap"
GAP_LAW = "E ∝ d^{-3}; coverage is left unchanged because demonstrated efficiency is 0"
GAP_SOURCE = "casimir-static.yaml E = -(π² ħ c A) / (720 d³); d derived is the walk radius"
DHE3_FORMULA_ID = "confined_dhe3_yield_j"
DHE3_LAW = "f = n_He3 <σv>(T) τ, joules per reference kilogram = energy_density * f"
DHE3_SOURCE = "src/solcascade/plasma/reactivity.py dhe3_reactivity_m3_s times gyro_walk_time_s"
DHE3_ID = "fusion-d-he3"
PB11_FORMULA_ID = "confined_pb11_yield_j"
PB11_LAW = "f = n_B <σv>(T) τ, joules per reference kilogram = energy_density * f"
PB11_SOURCE = "src/solcascade/plasma/reactivity.py pb11_reactivity_m3_s times gyro_walk_time_s"
PB11_ID = "fusion-aneutronic-pb11"
BURN_FORMULA_IDS = frozenset({DHE3_FORMULA_ID, PB11_FORMULA_ID})


@dataclass(frozen=True)
class WalkClock:
    radius_m: float
    temperature_kev: float
    field_t: float
    walk_time_s: float
    bohm_time_s: float
    walk_over_bohm: float
    gamma_enters_the_time: bool


@dataclass(frozen=True)
class ConfinedYield:
    energy_j: float
    power_w: float
    density_m3: float
    temperature_kev: float
    reactivity_m3_s: float
    volume_m3: float
    deuterium_fraction: float
    radius0_m: float
    density0_m3: float
    temperature0_kev: float


@dataclass(frozen=True)
class Derivation:
    mechanism_id: str
    quantity: str
    unit: str
    value: float
    coverage_factor: float | None
    absolute_joules: float | None
    formula_id: str
    source: str
    detail: str = ""
    side_factor: float | None = None
    side_quantity: str = ""


@dataclass(frozen=True)
class LoadedOverrides:
    clock: WalkClock
    confined_yield: ConfinedYield
    derivations: tuple[Derivation, ...]


def default_overrides_path() -> Path:
    return Path(__file__).resolve().parent / "overrides.yaml"


def walk_clock(radius_m: float, temperature_kev: float, field_t: float) -> WalkClock:
    """Evaluate the installed clock. The seconds are a function of a, T, and B."""
    walk = float(gyro_walk_time_s(radius_m, temperature_kev, field_t))
    bohm = float(bohm_time_s(radius_m, temperature_kev, field_t))
    return WalkClock(
        radius_m=float(radius_m),
        temperature_kev=float(temperature_kev),
        field_t=float(field_t),
        walk_time_s=walk,
        bohm_time_s=bohm,
        walk_over_bohm=walk / bohm,
        gamma_enters_the_time=False,
    )


def _positive(name: str, value) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be a positive finite number")
    return number


def confined_dt_yield(
    clock: WalkClock,
    *,
    density0_m3: float,
    temperature0_kev: float,
    radius0_m: float,
    deuterium_fraction: float,
) -> ConfinedYield:
    """Joules from one walk time at the compressed peak.

    Density scales as the inverse volume and temperature scales as the
    monatomic sphere in ``evolve``. The reactivity is the project's DT fit.
    """
    if not 0.0 < deuterium_fraction < 1.0:
        raise ValueError("deuterium fraction must be between 0 and 1")
    ratio = float(radius0_m) / clock.radius_m
    thermal_power = -adiabatic_temperature_exponent(MONATOMIC_ADIABATIC_INDEX)
    temperature = float(temperature0_kev) * ratio ** thermal_power
    if not math.isclose(temperature, clock.temperature_kev, rel_tol=1e-9, abs_tol=0.0):
        raise ValueError("peak temperature does not match the adiabatic compression")
    density = float(density0_m3) * ratio ** 3
    n_d = deuterium_fraction * density
    n_t = (1.0 - deuterium_fraction) * density
    reactivity = float(dt_reactivity_m3_s(temperature))
    volume = (4.0 / 3.0) * math.pi * clock.radius_m ** 3
    power = n_d * n_t * reactivity * float(dt_reaction().q_joule) * volume
    return ConfinedYield(
        energy_j=power * clock.walk_time_s,
        power_w=power,
        density_m3=density,
        temperature_kev=temperature,
        reactivity_m3_s=reactivity,
        volume_m3=volume,
        deuterium_fraction=float(deuterium_fraction),
        radius0_m=float(radius0_m),
        density0_m3=float(density0_m3),
        temperature0_kev=float(temperature0_kev),
    )


def default_catalog_path() -> Path:
    return Path(__file__).resolve().parent.parent / "data"


def thermal_ion_speed_m_s(temperature_kev: float) -> float:
    """Ion thermal speed inside the walk step, v = sqrt(kT / m).

    m is the mean of the deuterium and tritium ion masses, the same mean
    ``confinement_check`` uses. The mass cancels in the walk time and remains
    in the speed.
    """
    mass = 0.5 * (
        ion_mass_kg(2, MASS_EXCESS_DEUTERIUM_KEV, 1)
        + ion_mass_kg(3, MASS_EXCESS_TRITIUM_KEV, 1)
    )
    thermal = float(temperature_kev) * JOULE_PER_EV * EV_PER_KEV
    return math.sqrt(thermal / mass)


def hydraulic_head_m(energy_density_j: float) -> float:
    """Head recovered from m g h for a one-kilogram reference."""
    return float(energy_density_j) / CATALOG_GRAVITY_M_S2


def tidal_range_m(energy_density_j: float) -> float:
    """Range recovered from (1/2) rho g R^2."""
    scale = 0.5 * CATALOG_WATER_DENSITY_KG_M3 * CATALOG_GRAVITY_M_S2
    return math.sqrt(float(energy_density_j) / scale)


def catalog_wind_speed_m_s(power_density_w: float, power_coefficient: float) -> float:
    """Speed recovered from (1/2) rho v^3 = electrical power / Cp."""
    if power_coefficient <= 0.0:
        raise ValueError("wind power coefficient must be positive")
    kinetic = float(power_density_w) / float(power_coefficient)
    return (2.0 * kinetic / CATALOG_AIR_DENSITY_KG_M3) ** (1.0 / 3.0)


def blackbody_energy_tail(x: float) -> float:
    """Integral from x to infinity of t^3 / (exp(t) - 1).

    The Debye series is exact. A term whose exponential underflows in
    float64 is zero, and every later term is smaller, so the sum stops.
    """
    if x < 0.0 or not math.isfinite(x):
        raise ValueError("x must be a non-negative finite number")
    total = 0.0
    for m in range(1, 4097):
        decay = math.exp(-m * x)
        if decay == 0.0:
            break
        inv = 1.0 / m
        poly = inv * (x ** 3 + inv * (3.0 * x ** 2 + inv * (6.0 * x + inv * 6.0)))
        total += decay * poly
    return total


def blackbody_energy_fraction_above(x: float) -> float:
    """Share of blackbody energy in photons above x = E / kT."""
    if x == 0.0:
        return 1.0
    full = math.pi ** 4 / 15.0
    return blackbody_energy_tail(x) / full


def solar_blackbody_kt_ev() -> float:
    """kT of the 6000 K solar blackbody in the cited ceiling paper, in eV."""
    return BOLTZMANN_J_K * SOLAR_BLACKBODY_K * EV_PER_JOULE


def photovoltaic_coverage_factor(gap_ev: float, demonstrated_efficiency: float) -> float:
    """Coverage factor from the above-gap share of the 6000 K blackbody.

    A zero tail means the junction absorbs none of that sunlight, so the
    factor is 0. A positive tail replaces the catalog electrical yield with
    the insolation above the gap: factor = fraction / demonstrated efficiency.
    """
    efficiency = float(demonstrated_efficiency)
    if efficiency < 0.0:
        raise ValueError("demonstrated efficiency must be non-negative")
    kt_ev = solar_blackbody_kt_ev()
    fraction = blackbody_energy_fraction_above(float(gap_ev) / kt_ev)
    if fraction == 0.0:
        return 0.0
    if efficiency == 0.0:
        raise ValueError("a positive tail cannot scale a zero catalog efficiency")
    return fraction / efficiency


def casimir_energy_magnitude_j(area_m2: float, gap_m: float) -> float:
    """Magnitude of E = -(π² ħ c A) / (720 d³)."""
    return (
        math.pi ** 2
        * REDUCED_PLANCK_J_S
        * SPEED_OF_LIGHT_M_S
        * float(area_m2)
        / (720.0 * float(gap_m) ** 3)
    )


def _catalog_by_id(catalog: Path) -> dict[str, dict]:
    found: dict[str, dict] = {}
    for record in load_records(catalog):
        if record.get("_load_error"):
            raise ValueError(f"{record.get('id')}: {record['_load_error']}")
        found[str(record.get("id"))] = record
    return found


def _required_number(record: dict, key: str) -> float:
    value = record.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{record.get('id')} is missing a numeric {key}")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{record.get('id')} {key} is not finite")
    return number


def _dhe3_derivation(clock: WalkClock, confined: ConfinedYield, record: dict) -> Derivation:
    """Joules per reference kilogram at the burn fraction of one walk time.

    The partner density is the pulse's stoichiometric split. The checked-in
    gas is deuterium-tritium. The reactivity is the D-He3 fit. The catalog
    joules are a complete burn of one kilogram, so the pulse volume is not
    the unit of the row.
    """
    gross = _required_number(record, "energy_density_j")
    if gross <= 0.0:
        raise ValueError("D-He3 energy density must be positive")
    partner = (1.0 - confined.deuterium_fraction) * confined.density_m3
    reactivity = float(dhe3_reactivity_m3_s(confined.temperature_kev))
    burn = partner * reactivity * clock.walk_time_s
    if not 0.0 < burn <= 1.0:
        raise ValueError("D-He3 burn fraction is outside the constant-density interval")
    joules = gross * burn
    detail = (
        f"Partner density `{partner:.6g}` m^-3, reactivity `{reactivity:.6g}` m^3/s, "
        f"burn fraction `{burn:.6g}`, yield `{joules:.6g}` J per reference kilogram. "
        "The 50/50 split is the pulse deuterium fraction. It matches a stoichiometric "
        "number density. The fit is evaluated only inside 0.5–190 keV."
    )
    return Derivation(
        mechanism_id=DHE3_ID,
        quantity="confined D-He3 yield per reference kilogram",
        unit="J",
        value=joules,
        coverage_factor=None,
        absolute_joules=joules,
        formula_id=DHE3_FORMULA_ID,
        source=DHE3_SOURCE,
        detail=detail,
    )


def _pb11_derivation(clock: WalkClock, confined: ConfinedYield, record: dict) -> Derivation:
    """Joules per reference kilogram at the p-11B burn fraction of one walk time.

    The partner density follows the same half-and-half number split as D-He3.
    The catalog kilogram is stoichiometric proton plus boron-11. The checked-in
    plasma is deuterium-tritium. The reactivity is the Nevins-Swain p-11B fit.
    """
    gross = _required_number(record, "energy_density_j")
    if gross <= 0.0:
        raise ValueError("p-B11 energy density must be positive")
    partner = (1.0 - confined.deuterium_fraction) * confined.density_m3
    reactivity = float(pb11_reactivity_m3_s(confined.temperature_kev))
    burn = partner * reactivity * clock.walk_time_s
    if not 0.0 < burn <= 1.0:
        raise ValueError("p-B11 burn fraction is outside the constant-density interval")
    joules = gross * burn
    detail = (
        f"Partner density `{partner:.6g}` m^-3, reactivity `{reactivity:.6g}` m^3/s, "
        f"burn fraction `{burn:.6g}`, yield `{joules:.6g}` J per reference kilogram. "
        "The 50/50 split is the pulse deuterium fraction. It matches one boron "
        "per proton on the stoichiometric kilogram. The checked-in plasma is "
        "deuterium-tritium. The reactivity is p-11B. The fit is the Nevins-Swain "
        "low-temperature branch, evaluated only inside 1 keV to 70 keV."
    )
    return Derivation(
        mechanism_id=PB11_ID,
        quantity="confined p-B11 yield per reference kilogram",
        unit="J",
        value=joules,
        coverage_factor=None,
        absolute_joules=joules,
        formula_id=PB11_FORMULA_ID,
        source=PB11_SOURCE,
        detail=detail,
    )


def _pulse_unit_derivations(
    clock: WalkClock,
    records: dict[str, dict],
) -> list[Derivation]:
    """Head, range, gap, wind speed, and bandgap from the peak state."""
    length = clock.radius_m
    speed = thermal_ion_speed_m_s(clock.temperature_kev)
    gap_ev = clock.temperature_kev * EV_PER_KEV
    hydro = records["hydroelectric"]
    tidal = records["tidal"]
    wind = records["wind"]
    photovoltaic = records["photovoltaic"]
    casimir = records["casimir-static"]
    head = hydraulic_head_m(_required_number(hydro, "energy_density_j"))
    tidal_range = tidal_range_m(_required_number(tidal, "energy_density_j"))
    wind_speed = catalog_wind_speed_m_s(
        _required_number(wind, "power_density_w"),
        _required_number(wind, "demonstrated_efficiency"),
    )
    if head <= 0.0 or tidal_range <= 0.0 or wind_speed <= 0.0:
        raise ValueError("catalog head, range, and wind speed must be positive")
    head_factor = length / head
    range_factor = (length / tidal_range) ** 2
    wind_factor = (speed / wind_speed) ** 3
    pv_factor = photovoltaic_coverage_factor(
        gap_ev,
        _required_number(photovoltaic, "demonstrated_efficiency"),
    )
    exact_cavity = casimir_energy_magnitude_j(CASIMIR_STATED_AREA_M2, CASIMIR_STATED_GAP_M)
    cavity_factor = (CASIMIR_STATED_GAP_M / length) ** 3
    file_cavity = _required_number(casimir, "energy_density_j")
    return [
        Derivation(
            mechanism_id="hydroelectric",
            quantity="hydraulic head",
            unit="m",
            value=length,
            coverage_factor=head_factor,
            absolute_joules=None,
            formula_id=HEAD_FORMULA_ID,
            source=HEAD_SOURCE,
            detail=(
                f"Catalog head `{head:.6g}` m from energy density over g. "
                f"Derived head `{length:.6g}` m. Factor `{head_factor:.12g}`. "
                "The length is the confinement radius. The dam is not the plasma."
            ),
        ),
        Derivation(
            mechanism_id="tidal",
            quantity="tidal range",
            unit="m",
            value=length,
            coverage_factor=range_factor,
            absolute_joules=None,
            formula_id=TIDAL_FORMULA_ID,
            source=TIDAL_SOURCE,
            detail=(
                f"Catalog range `{tidal_range:.6g}` m from (1/2) rho g R^2. "
                f"Derived range `{length:.6g}` m. Factor `{range_factor:.12g}`. "
                "The length is the confinement radius. The basin is not the plasma."
            ),
        ),
        Derivation(
            mechanism_id="wind",
            quantity="wind speed",
            unit="m/s",
            value=speed,
            coverage_factor=wind_factor,
            absolute_joules=None,
            formula_id=WIND_FORMULA_ID,
            source=WIND_SOURCE,
            detail=(
                f"Catalog speed `{wind_speed:.6g}` m/s from (1/2) rho v^3 = power / Cp. "
                f"Thermal speed `{speed:.6g}` m/s. Factor `{wind_factor:.12g}`. "
                "The speed is the ion thermal speed. The rotor is not the plasma."
            ),
        ),
        Derivation(
            mechanism_id="photovoltaic",
            quantity="bandgap",
            unit="eV",
            value=gap_ev,
            coverage_factor=pv_factor,
            absolute_joules=None,
            formula_id=PV_FORMULA_ID,
            source=PV_SOURCE,
            detail=(
                f"Pulse kT `{gap_ev:.6g}` eV. Solar kT `{solar_blackbody_kt_ev():.6g}` eV. "
                f"Above-gap fraction gives coverage factor `{pv_factor:.12g}`. "
                "The electron-volt is the pulse kT. It is not a silicon gap."
            ),
        ),
        Derivation(
            mechanism_id="casimir-static",
            quantity="plate separation",
            unit="m",
            value=length,
            coverage_factor=None,
            absolute_joules=None,
            formula_id=GAP_FORMULA_ID,
            source=GAP_SOURCE,
            detail=(
                f"Stated geometry `{CASIMIR_STATED_AREA_M2:.6g}` m^2 at "
                f"`{CASIMIR_STATED_GAP_M:.6g}` m gives cavity magnitude "
                f"`{exact_cavity:.6g}` J. The file stores `{file_cavity:.6g}` J. "
                f"Derived gap `{length:.6g}` m. Cavity-energy factor "
                f"`{cavity_factor:.12g}`. Demonstrated efficiency is 0, so "
                "delivered joules stay 0 and the coverage stays 0."
            ),
            side_factor=cavity_factor,
            side_quantity="cavity energy factor",
        ),
    ]


def load_overrides(path: Path | None = None, catalog: Path | None = None) -> LoadedOverrides:
    """Load the mapping file, the walk time, the yields, and the pulse-unit factors."""
    path = default_overrides_path() if path is None else Path(path)
    catalog_path = default_catalog_path() if catalog is None else Path(catalog)
    records = _catalog_by_id(catalog_path)
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(loaded, dict):
        raise ValueError("overrides file is not a mapping")
    if loaded.get("formula_id") != FORMULA_ID:
        raise ValueError(f"formula_id must be {FORMULA_ID}")
    if loaded.get("yield_formula_id") != YIELD_FORMULA_ID:
        raise ValueError(f"yield_formula_id must be {YIELD_FORMULA_ID}")
    if loaded.get("gamma_enters_the_time") is not False:
        raise ValueError("gamma_enters_the_time must be false")
    applies = tuple(loaded.get("applies_to") or ())
    if applies != FUSION_IDS:
        raise ValueError("applies_to must be the five fusion mechanism ids")
    yield_ids = tuple(loaded.get("yield_applies_to") or ())
    if yield_ids != DT_YIELD_IDS:
        raise ValueError("yield_applies_to must be the three DT mechanism ids")
    clock = walk_clock(
        _positive("radius_m", loaded["radius_m"]),
        _positive("temperature_keV", loaded["temperature_keV"]),
        _positive("field_T", loaded["field_T"]),
    )
    fraction = float(loaded["deuterium_fraction"])
    confined = confined_dt_yield(
        clock,
        density0_m3=_positive("ion_density_m3", loaded["ion_density_m3"]),
        temperature0_kev=_positive("ion_temperature_keV", loaded["ion_temperature_keV"]),
        radius0_m=_positive("radius0_m", loaded["radius0_m"]),
        deuterium_fraction=fraction,
    )
    source = str(loaded.get("source") or SOURCE)
    derivations = []
    for mechanism_id in applies:
        if mechanism_id in yield_ids:
            derivations.append(
                Derivation(
                    mechanism_id=mechanism_id,
                    quantity="confined DT yield",
                    unit="J",
                    value=confined.energy_j,
                    coverage_factor=None,
                    absolute_joules=confined.energy_j,
                    formula_id=YIELD_FORMULA_ID,
                    source=YIELD_SOURCE,
                )
            )
            continue
        if mechanism_id == DHE3_ID:
            derivations.append(_dhe3_derivation(clock, confined, records[DHE3_ID]))
            continue
        if mechanism_id == PB11_ID:
            derivations.append(_pb11_derivation(clock, confined, records[PB11_ID]))
            continue
        derivations.append(
            Derivation(
                mechanism_id=mechanism_id,
                quantity="confinement time",
                unit="s",
                value=clock.walk_time_s,
                coverage_factor=1.0,
                absolute_joules=None,
                formula_id=FORMULA_ID,
                source=source,
            )
        )
    derivations.extend(_pulse_unit_derivations(clock, records))
    return LoadedOverrides(clock, confined, tuple(derivations))


def apply_coverage_factor(
    row: dict,
    factor: float,
    *,
    quantity: str,
    source: str,
    formula_id: str,
    value: float | None = None,
    unit: str | None = None,
) -> dict:
    """Scale one scored row. A blank coverage stays blank."""
    updated = dict(row)
    old = row.get("output_over_target")
    updated["old_coverage"] = old
    updated["coverage_factor"] = float(factor)
    updated["derived_quantity"] = quantity
    updated["derived_source"] = source
    updated["derived_formula"] = formula_id
    updated["derived_value"] = value
    updated["derived_unit"] = unit
    updated["derived_side_factor"] = None
    updated["derived_side_quantity"] = ""
    if old is None:
        updated["new_coverage"] = None
        updated["output_over_target"] = None
        return updated
    new = float(old) * float(factor)
    updated["new_coverage"] = new
    updated["output_over_target"] = new
    delivered = row.get("delivered_j")
    if isinstance(delivered, (int, float)) and not isinstance(delivered, bool):
        updated["delivered_j"] = float(delivered) * float(factor)
    units = row.get("units_to_target")
    if isinstance(units, (int, float)) and not isinstance(units, bool) and factor != 0.0:
        updated["units_to_target"] = float(units) / float(factor)
    ratio = row.get("output_over_input")
    if isinstance(ratio, (int, float)) and not isinstance(ratio, bool):
        updated["output_over_input"] = float(ratio) * float(factor)
    return updated


def apply_absolute_joules(
    row: dict,
    joules: float,
    *,
    target_j: float,
    quantity: str,
    source: str,
    formula_id: str,
) -> dict:
    """Replace delivered joules with a derived yield. A blank row gains a coverage."""
    updated = dict(row)
    old = row.get("output_over_target")
    delivered = row.get("delivered_j")
    updated["old_coverage"] = old
    updated["derived_quantity"] = quantity
    updated["derived_source"] = source
    updated["derived_formula"] = formula_id
    updated["derived_value"] = float(joules)
    updated["derived_unit"] = "J"
    updated["derived_side_factor"] = None
    updated["derived_side_quantity"] = ""
    updated["delivered_j"] = float(joules)
    if target_j > 0.0:
        coverage = float(joules) / float(target_j)
        updated["new_coverage"] = coverage
        updated["output_over_target"] = coverage
        if joules > 0.0:
            updated["units_to_target"] = float(target_j) / float(joules)
    else:
        updated["new_coverage"] = None
        updated["output_over_target"] = None
    incoming = row.get("input_energy_j")
    if isinstance(incoming, (int, float)) and not isinstance(incoming, bool) and incoming > 0.0:
        updated["output_over_input"] = float(joules) / float(incoming)
    if (
        isinstance(delivered, (int, float))
        and not isinstance(delivered, bool)
        and delivered > 0.0
        and old is not None
    ):
        updated["coverage_factor"] = float(joules) / float(delivered)
    else:
        updated["coverage_factor"] = None
    return updated


def _blank_annotation(row: dict) -> dict:
    updated = dict(row)
    coverage = row.get("output_over_target")
    updated["old_coverage"] = coverage
    updated["new_coverage"] = coverage
    updated["coverage_factor"] = None
    updated["derived_quantity"] = None
    updated["derived_source"] = None
    updated["derived_formula"] = None
    updated["derived_value"] = None
    updated["derived_unit"] = None
    updated["derived_side_factor"] = None
    updated["derived_side_quantity"] = ""
    return updated


def apply_recorded_only(
    row: dict,
    *,
    quantity: str,
    source: str,
    formula_id: str,
    value: float,
    unit: str,
    side_factor: float | None = None,
    side_quantity: str = "",
) -> dict:
    """Record a derived quantity and leave the phase-one coverage where it is."""
    updated = dict(row)
    coverage = row.get("output_over_target")
    updated["old_coverage"] = coverage
    updated["new_coverage"] = coverage
    updated["output_over_target"] = coverage
    updated["coverage_factor"] = None
    updated["derived_quantity"] = quantity
    updated["derived_source"] = source
    updated["derived_formula"] = formula_id
    updated["derived_value"] = float(value)
    updated["derived_unit"] = unit
    updated["derived_side_factor"] = side_factor
    updated["derived_side_quantity"] = side_quantity
    return updated


def annotate(rows: list[dict], derivations: tuple[Derivation, ...], target_j: float) -> list[dict]:
    """Attach a derivation where one exists. Other rows keep phase-one coverage."""
    by_id = {item.mechanism_id: item for item in derivations}
    annotated: list[dict] = []
    for row in rows:
        item = by_id.get(row.get("id"))
        if item is None:
            annotated.append(_blank_annotation(row))
            continue
        if item.absolute_joules is not None:
            annotated.append(
                apply_absolute_joules(
                    row,
                    item.absolute_joules,
                    target_j=target_j,
                    quantity=item.quantity,
                    source=item.source,
                    formula_id=item.formula_id,
                )
            )
            continue
        if item.coverage_factor is None:
            annotated.append(
                apply_recorded_only(
                    row,
                    quantity=item.quantity,
                    source=item.source,
                    formula_id=item.formula_id,
                    value=item.value,
                    unit=item.unit,
                    side_factor=item.side_factor,
                    side_quantity=item.side_quantity,
                )
            )
            continue
        annotated.append(
            apply_coverage_factor(
                row,
                float(item.coverage_factor),
                quantity=item.quantity,
                source=item.source,
                formula_id=item.formula_id,
                value=item.value,
                unit=item.unit,
            )
        )
    return annotated


def coverage_changed(old, new) -> bool:
    if old is None and new is None:
        return False
    if old is None or new is None:
        return True
    return not math.isclose(float(old), float(new), rel_tol=1e-12, abs_tol=0.0)


def order_ratio(old, new) -> float | None:
    """New coverage over old coverage. A fall from a positive coverage to 0 is 0."""
    if old is None or new is None:
        return None
    old_f = float(old)
    new_f = float(new)
    if old_f == 0.0 or old_f < 0.0 or new_f < 0.0:
        return None
    if new_f == 0.0:
        return 0.0
    return new_f / old_f


def moves_by_order(old, new) -> bool:
    ratio = order_ratio(old, new)
    if ratio is None or ratio < 0.0:
        return False
    if ratio == 0.0:
        return True
    return abs(math.log10(ratio)) >= math.log10(ORDER_OF_MAGNITUDE) - 1e-12


def unmapped_parameters() -> list[tuple[str, str]]:
    """Dimensionful catalog levers, in phase-two reading order."""
    from solcascade.energy_baseline.phase2.map_fsot import READINGS

    found: list[tuple[str, str]] = []
    for mechanism_id, reading in READINGS.items():
        if reading.kind == "dimensionful":
            found.append((mechanism_id, reading.parameter))
    return found


def phase3_gap_section(clock: WalkClock | None = None) -> str:
    """Markdown appended to the phase-two gap table. Deterministic."""
    loaded = load_overrides()
    clock = loaded.clock if clock is None else clock
    confined = loaded.confined_yield
    lines = [
        "## Phase three",
        "",
        f"Confinement time for the SolCascade pulse is `{clock.walk_time_s:.6g}` s.",
        f"Formula `{FORMULA_ID}`: {LAW}.",
        (
            f"Inputs: radius `{clock.radius_m:.6g}` m, "
            f"temperature `{clock.temperature_kev:.6g}` keV, "
            f"field `{clock.field_t:.6g}` T."
        ),
        (
            f"Source: `{SOURCE}`. The inputs are stored in "
            f"`phase3/overrides.yaml` and match `{INPUTS_SOURCE}`."
        ),
        (
            f"walk/Bohm is `{clock.walk_over_bohm:.6g}`, equal to π/4. "
            "`gamma_enters_the_time` is false."
        ),
        "",
        (
            f"Confined DT yield: `{confined.energy_j:.6g}` J from `{YIELD_FORMULA_ID}`: "
            f"{YIELD_LAW}."
        ),
        (
            "n and T are the adiabatic compression of the checked-in pulse from "
            f"`radius0` `{confined.radius0_m:.6g}` m, "
            f"`n0` `{confined.density0_m3:.6g}` m^-3, "
            f"`T0` `{confined.temperature0_kev:.6g}` keV "
            "to the peak radius. Q is `dt_reaction().q_joule`. "
            "`<σv>` is `dt_reactivity_m3_s`. τ is the walk time. "
            "The three DT rows take this yield in place of the catalog delivered joules."
        ),
        (
            "That substitution moves the scored DT coverages down, by the ratio of "
            "this yield to the catalog yield. D-He3 multiplies the catalog joules "
            "per reference kilogram by the burn fraction "
            "`f = n_He3 <σv>(T) τ` at the same peak. p-B11 multiplies its "
            "catalog joules by `f = n_B <σv>(T) τ` on the Nevins-Swain "
            "low-temperature branch."
        ),
        (
            "Hydraulic head, tidal range, and wind speed scale coverage by the "
            "power already recorded for those rows: linear in head, squared in "
            "range, cubed in speed. The length is the confinement radius. The "
            "speed is the ion thermal speed. The photovoltaic bandgap is that "
            "kT in electron-volts, scored as the above-gap fraction of the "
            "6000 K blackbody in the cited ceiling paper. The Casimir plate "
            "separation is the same radius. Cavity energy scales as the inverse "
            "cube. Demonstrated efficiency is 0, so that coverage stays 0."
        ),
        "",
        "Installed factors from this peak:",
        "",
    ]
    for item in loaded.derivations:
        if item.coverage_factor is not None and item.formula_id != FORMULA_ID:
            lines.append(
                f"- `{item.mechanism_id}`: coverage factor `{item.coverage_factor:.12g}` "
                f"via `{item.formula_id}`."
            )
        elif item.side_factor is not None:
            lines.append(
                f"- `{item.mechanism_id}`: cavity-energy factor `{item.side_factor:.12g}` "
                f"via `{item.formula_id}`. Coverage is unchanged."
            )
        elif item.formula_id in BURN_FORMULA_IDS:
            lines.append(f"- `{item.mechanism_id}`: {item.detail}")
    lines.extend(
        [
            "",
            "Dimensionful parameters:",
            "",
        ]
    )
    for mechanism_id, parameter in unmapped_parameters():
        lines.append(f"- `{mechanism_id}`: {parameter} ({_dimensionful_status(mechanism_id)})")
    lines.extend(["", NEXT_TARGET, ""])
    return "\n".join(lines)


def _dimensionful_status(mechanism_id: str) -> str:
    if mechanism_id in DT_YIELD_IDS:
        return "DT yield installed; see the phase-three ranking for the factor"
    if mechanism_id == DHE3_ID:
        return "D-He3 burn fraction installed; the phase-one coverage was blank"
    if mechanism_id == PB11_ID:
        return "p-B11 burn fraction installed; the phase-one coverage was blank"
    if mechanism_id == "casimir-static":
        return "plate gap installed; coverage stays 0 because demonstrated efficiency is 0"
    if mechanism_id in {"hydroelectric", "tidal", "wind", "photovoltaic"}:
        return "pulse unit installed; see the phase-three ranking for the factor"
    return "no derivation installed"
