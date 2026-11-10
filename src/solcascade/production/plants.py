"""Score the operating fission plants and leave the FSOT scalars beside them.

The thermal and fast rows are the catalog plants that already run. Delivered
joules are the catalog product of energy density and demonstrated efficiency.
The nuclear scalar is recorded and is not multiplied into those joules. The
checked-in pulse is not rewritten.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import (
    PIN_PREFIX,
    READINGS,
    capture_scalars,
    load_engine,
)
from solcascade.production.refine import catalog_record, project_root

THERMAL_ID = "fission-thermal"
FAST_ID = "fission-fast"
_DISCHARGE_GWD = re.compile(r"(\d+(?:\.\d+)?) GWd")
_SECONDS_PER_DAY = 86400.0
_WATT_PER_GIGAWATT = 1.0e9
_KILOGRAMS_PER_TONNE = 1000.0


@dataclass(frozen=True)
class Plant:
    """One catalog plant on its own reference unit."""

    mechanism_id: str
    name: str
    trl: int
    reference_unit: str
    efficiency: float
    ceiling: float
    headroom: float
    energy_density_j: float
    delivered_j: float
    ceiling_delivered_j: float
    input_j: float
    output_over_input: float
    self_sustaining: bool
    continuous_external_drive: bool


@dataclass(frozen=True)
class PlantsReport:
    thermal: Plant
    fast: Plant
    efficiency_step: float
    discharge_gwd_per_tonne: float
    discharge_j_per_kg: float
    discharge_over_ranked: float
    burnup_demonstrated: float
    burnup_cap: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "plants.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _plant(mechanism_id: str) -> tuple[Plant, str]:
    record = catalog_record(mechanism_id)
    if record.get("fsot_status") != "open_gap" or record.get("fsot_scalar") is not None:
        raise RuntimeError(f"{mechanism_id} carries a scalar mapping")
    if record.get("self_sustaining") is not True:
        raise RuntimeError(f"{mechanism_id} is not a self-sustaining chain")
    if record.get("continuous_external_drive") is not False:
        raise RuntimeError(f"{mechanism_id} asks for a continuous external drive")
    efficiency = float(record["demonstrated_efficiency"])
    ceiling = float(record["efficiency_ceiling"])
    density = float(record["energy_density_j"])
    incoming = float(record["input_energy_j"])
    if not 0.0 < efficiency <= 1.0:
        raise RuntimeError(f"{mechanism_id} efficiency is outside (0, 1]")
    if ceiling < efficiency or ceiling > 1.0:
        raise RuntimeError(f"{mechanism_id} ceiling is below the demonstrated efficiency")
    if density <= 0.0 or incoming <= 0.0:
        raise RuntimeError(f"{mechanism_id} joules are not positive")
    delivered = density * efficiency
    return (
        Plant(
            mechanism_id=mechanism_id,
            name=str(record["name"]),
            trl=int(record["trl"]),
            reference_unit=str(record["reference_unit"]),
            efficiency=efficiency,
            ceiling=ceiling,
            headroom=ceiling / efficiency,
            energy_density_j=density,
            delivered_j=delivered,
            ceiling_delivered_j=density * ceiling,
            input_j=incoming,
            output_over_input=delivered / incoming,
            self_sustaining=True,
            continuous_external_drive=False,
        ),
        str(record.get("notes") or ""),
    )


def discharge_j_per_kg(notes: str) -> tuple[float, float]:
    """Convert the thermal note's gigawatt-day discharge to joules per kilogram.

    The note cites the discharge. A gigawatt-day is 1e9 watts times 86400
    seconds, and a tonne is 1000 kilograms.
    """
    match = _DISCHARGE_GWD.search(notes)
    if match is None or "heavy metal" not in notes:
        raise RuntimeError("the thermal notes do not cite a heavy-metal discharge")
    gigawatt_days = float(match.group(1))
    if gigawatt_days <= 0.0:
        raise RuntimeError("the discharge citation is not positive")
    per_kg = gigawatt_days * _WATT_PER_GIGAWATT * _SECONDS_PER_DAY / _KILOGRAMS_PER_TONNE
    return gigawatt_days, per_kg


def _witness() -> tuple[str, float, dict[str, float]]:
    module, digest = load_engine()
    if not str(digest).startswith(PIN_PREFIX):
        raise RuntimeError(f"engine pin is {str(digest)[:16]}")
    scalars = capture_scalars(module)
    valve = float(module.BETA)
    if valve <= 0.0 or valve >= 1.0:
        raise RuntimeError("engine valve is outside the recorded range")
    return str(digest)[:16], valve, scalars


def build_plants() -> PlantsReport:
    thermal, notes = _plant(THERMAL_ID)
    fast, _fast_notes = _plant(FAST_ID)
    gigawatt_days, discharge = discharge_j_per_kg(notes)
    thermal_reading = READINGS[THERMAL_ID]
    fast_reading = READINGS[FAST_ID]
    if thermal_reading.demonstrated != fast_reading.demonstrated:
        raise RuntimeError("the two fission burnup readings disagree")
    if thermal_reading.physical_max != fast_reading.physical_max:
        raise RuntimeError("the two fission burnup caps disagree")
    demonstrated = float(thermal_reading.demonstrated)
    cap = float(thermal_reading.physical_max)
    if demonstrated <= 0.0 or cap < demonstrated:
        raise RuntimeError("the phase-two burnup reading is not a closed unit")
    digest, valve, scalars = _witness()
    nuclear = float(scalars["Nuclear_Physics"])
    thermo = float(scalars["Thermodynamics"])
    if not math.isfinite(nuclear) or not math.isfinite(thermo):
        raise RuntimeError("a witness scalar is not finite")
    return PlantsReport(
        thermal=thermal,
        fast=fast,
        efficiency_step=fast.efficiency / thermal.efficiency,
        discharge_gwd_per_tonne=gigawatt_days,
        discharge_j_per_kg=discharge,
        discharge_over_ranked=discharge / thermal.energy_density_j,
        burnup_demonstrated=demonstrated,
        burnup_cap=cap,
        nuclear_scalar=nuclear,
        thermodynamics_scalar=thermo,
        gamma_fsot=float(scalars["gamma_FSOT"]),
        valve_beta=valve,
        pin_prefix=PIN_PREFIX,
        digest_prefix=digest,
    )


def _row(plant: Plant) -> str:
    return (
        f"| {plant.name} | {plant.trl} | {plant.efficiency:.6g} | {plant.ceiling:.6g} | "
        f"{plant.delivered_j:.6g} | {plant.ceiling_delivered_j:.6g} | "
        f"{plant.output_over_input:.6g} |"
    )


def render() -> str:
    report = build_plants()
    thermal = report.thermal
    fast = report.fast
    lines = [
        "# Working plants",
        "",
        "The thermal-reactor row is an operating chain. It is self-sustaining.",
        "The catalog sets continuous external drive to false.",
        "The fast-reactor row is the operating fast plant on the same fission energy.",
        "The chain sustains once critical.",
        "",
        "Delivered joules are the catalog product of energy density and",
        "demonstrated efficiency, on each row's own reference unit.",
        "",
        "| Plant | TRL | Efficiency | Ceiling | Delivered (J) | At the ceiling (J) | Output / enrichment estimate |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        _row(thermal),
        _row(fast),
        "",
        f"The thermal reference unit is {thermal.reference_unit}.",
        f"The fast reference unit is {fast.reference_unit}.",
        f"Thermal headroom, ceiling over demonstrated efficiency, is `{_num(thermal.headroom)}`.",
        f"Fast headroom is `{_num(fast.headroom)}`.",
        f"The fast demonstrated efficiency over the thermal efficiency is `{_num(report.efficiency_step)}`.",
        "",
        f"The thermal notes cite `{_num(report.discharge_gwd_per_tonne)}` GWd per tonne of heavy metal.",
        f"That discharge is `{_num(report.discharge_j_per_kg)}` J per kilogram of heavy metal.",
        f"The ranked thermal density is `{_num(thermal.energy_density_j)}` J per kilogram of pure U-235.",
        "Per kilogram, the discharge citation stands at "
        f"`{_num(report.discharge_over_ranked)}` of that ranked density.",
        "A heavy-metal kilogram and a pure U-235 kilogram are different loads.",
        "The enrichment estimate belongs to one kilogram of highly enriched uranium",
        "and is scored against the ranked kilogram. The fast row uses that same",
        "enrichment estimate.",
        "",
        "## FSOT record",
        "",
        f"Pin `{report.pin_prefix}` matches the loaded engine prefix `{report.digest_prefix}`.",
        "The scalars below are dimensionless.",
        "A domain scalar is not a plant efficiency.",
        "The nuclear scalar stays out of the fission joules.",
        f"Nuclear_Physics: `{_num(report.nuclear_scalar)}`.",
        f"Thermodynamics: `{_num(report.thermodynamics_scalar)}`.",
        f"gamma_FSOT: `{_num(report.gamma_fsot)}`.",
        "gamma_FSOT does not multiply the plant.",
        "The engine valve stays out of the fission joules.",
        f"Engine valve: `{_num(report.valve_beta)}`.",
        "",
        f"The phase-two burnup on the ranked unit is `{_num(report.burnup_demonstrated)}`",
        f"with a cap of `{_num(report.burnup_cap)}`.",
        "The tenfold window is empty.",
        "",
        "The refinement inside these plants is the file efficiency ceiling.",
        "The fast row is the operating plant on the actinides a thermal reactor",
        "leaves behind. The discharge citation is the gap between a heavy-metal",
        "tonne and the ranked pure kilogram. The 200 MeV ledger quantum, the",
        "efficiencies, and that discharge stay on the catalog figures.",
        "",
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
