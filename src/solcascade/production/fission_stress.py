"""Run the filed light-water plant across the conditions already on file.

The electric boundary stays the published nameplate. The demonstrated
operating heat is the FSOT heat ratio times that nameplate. The steam
figures already scored, and the shutdown times already scored, are the
other conditions. A domain scalar is not multiplied into the watts.
The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.production.fission_ceiling import build_ceiling
from solcascade.production.fission_decay_heat import HeatLoad, build_decay_heat
from solcascade.production.fission_process import (
    FABRICATION_CAP_PERCENT,
    NAMEPLATE_MWE,
    RELOAD_TONNES,
    build_process,
)
from solcascade.production.fission_safety import SMALL_CLASS_MWE, build_safety, heat_dump_mj_per_kwh
from solcascade.production.plants import FAST_ID, THERMAL_ID
from solcascade.production.refine import catalog_record, project_root

_REVIEW_GAP = 1.0e-3
_LOAD_ORDER = ("one-second fraction", "one-day fraction", "two-day fraction")
_HELD_K = 1.0


@dataclass(frozen=True)
class PowerCase:
    name: str
    electric_w: float
    thermal_w: float
    waste_w: float
    volume_m3: float
    mj_per_kwh: float
    source: str
    gap: float | None


@dataclass(frozen=True)
class ShutdownCase:
    name: str
    cooling_w: float
    waste_w: float
    share: float
    gap: float


@dataclass(frozen=True)
class StressReport:
    power: tuple[PowerCase, ...]
    shutdown: tuple[ShutdownCase, ...]
    catalog_ceiling: float
    carnot: float
    reload_tonnes: float
    nameplate_j: float
    fuel_side_j: float
    companion_ratio: float
    discharge_over_ranked: float
    stack_factor: float
    fabrication_cap_percent: float
    held_k: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_stress.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _gap(derived: float, witness: float) -> float:
    if witness == 0.0:
        raise RuntimeError("the witness is zero")
    return abs(derived - witness) / abs(witness)


def _power_case(
    name: str,
    electric_w: float,
    thermal_w: float,
    density: float,
    source: str,
    gap: float | None,
) -> PowerCase:
    if electric_w <= 0.0 or thermal_w <= electric_w:
        raise RuntimeError(f"{name} does not reject heat")
    waste = thermal_w - electric_w
    return PowerCase(
        name=name,
        electric_w=electric_w,
        thermal_w=thermal_w,
        waste_w=waste,
        volume_m3=thermal_w / density,
        mj_per_kwh=waste / electric_w * 3.6,
        source=source,
        gap=gap,
    )


def _shutdown_cases(
    loads: tuple[HeatLoad, ...],
    large_waste: float,
    small_waste: float,
) -> tuple[ShutdownCase, ...]:
    if tuple(load.label for load in loads) != _LOAD_ORDER:
        raise RuntimeError("the shutdown order changed")
    cases: list[ShutdownCase] = []
    for load in loads:
        for scale, cooling, waste in (
            ("large", load.watts, large_waste),
            ("small", load.small_watts, small_waste),
        ):
            if not 0.0 < cooling < waste:
                raise RuntimeError(f"the {scale} {load.label} is outside the condenser duty")
            cases.append(
                ShutdownCase(
                    name=f"{scale} {load.label}",
                    cooling_w=cooling,
                    waste_w=waste,
                    share=cooling / waste,
                    gap=load.relative_gap,
                )
            )
    return tuple(cases)


def build_stress() -> StressReport:
    process = build_process()
    safety = build_safety()
    ceiling = build_ceiling()
    heat = build_decay_heat()
    thermal = catalog_record(THERMAL_ID)
    fast = catalog_record(FAST_ID)
    governing = str(thermal["governing_equation"])
    if "k_eff is held at 1" not in governing:
        raise RuntimeError("the governing equation does not hold k_eff at 1")
    if "Sodium" not in str(fast.get("engineering_barrier") or ""):
        raise RuntimeError("the fast barrier does not name sodium")
    if heat.pin_prefix != process.pin_prefix or heat.digest_prefix != heat.pin_prefix:
        raise RuntimeError("the stress pin does not match the plant pin")
    if ceiling.carnot != process.carnot or ceiling.catalog_ceiling <= ceiling.carnot:
        raise RuntimeError("the file ceiling is not past the temperature Carnot")
    if safety.fabrication_cap_percent != FABRICATION_CAP_PERCENT:
        raise RuntimeError("the fabrication cap changed")
    for seed in (
        process.discharge_seed,
        process.enrichment_seed,
        process.atomic_seed,
        safety.size_seed,
        safety.built_seed,
        ceiling.demonstrated_seed,
        ceiling.carnot_seed,
        ceiling.fast_seed,
    ):
        if seed.relative_gap <= _REVIEW_GAP:
            raise RuntimeError(f"{seed.name} meets a filed figure closely enough to review")
    density = process.power_density_w_m3
    large_w = NAMEPLATE_MWE * 1.0e6
    small_w = SMALL_CLASS_MWE * 1.0e6
    heat_gap = _gap(heat.thermal_w, heat.filed_thermal_w)
    small_gap = _gap(heat.small_thermal_w, heat.filed_small_thermal_w)
    if heat_gap > heat.helium_mismatch or small_gap > heat.helium_mismatch:
        raise RuntimeError("the FSOT operating heat misses the helium quotient match")
    if abs(heat_gap - small_gap) > 1.0e-15:
        raise RuntimeError("the two FSOT heats do not share one gap")
    large = _power_case(
        "large nameplate",
        large_w,
        heat.thermal_w,
        density,
        "FSOT heat ratio",
        heat_gap,
    )
    small = _power_case(
        "small class",
        small_w,
        heat.small_thermal_w,
        density,
        "FSOT heat ratio",
        small_gap,
    )
    built = _power_case(
        "plants being built",
        large_w,
        large_w / safety.built_efficiency,
        density,
        "page efficiency",
        None,
    )
    carnot = _power_case(
        "temperature Carnot",
        large_w,
        large_w / process.carnot,
        density,
        "page temperatures",
        None,
    )
    if _gap(built.mj_per_kwh, heat_dump_mj_per_kwh(safety.built_efficiency)) > 1.0e-12:
        raise RuntimeError("the built condenser figure left the cooling identity")
    if _gap(carnot.mj_per_kwh, heat_dump_mj_per_kwh(process.carnot)) > 1.0e-12:
        raise RuntimeError("the Carnot condenser figure left the cooling identity")
    demonstrated_gap = _gap(large.mj_per_kwh, safety.dump_demonstrated_mj)
    if demonstrated_gap > heat.helium_mismatch:
        raise RuntimeError("the FSOT condenser figure misses the helium quotient match")
    if not carnot.waste_w < built.waste_w < large.waste_w:
        raise RuntimeError("the condenser duty does not fall as the steam figure rises")
    if not carnot.volume_m3 < built.volume_m3 < large.volume_m3:
        raise RuntimeError("the core volume does not fall as the steam figure rises")
    if small.electric_w / large.electric_w != SMALL_CLASS_MWE / NAMEPLATE_MWE:
        raise RuntimeError("the small class left the published size ratio")
    for case in (large, small):
        if case.gap is None or case.gap > heat.helium_mismatch:
            raise RuntimeError(f"{case.name} is missing a helium-tight heat")
    shutdown = _shutdown_cases(heat.loads, large.waste_w, small.waste_w)
    shares = [case.share for case in shutdown if case.name.startswith("large ")]
    if not shares[0] > shares[1] > shares[2]:
        raise RuntimeError("the shutdown load does not fall with time")
    for large_case, small_case in zip(shutdown[0::2], shutdown[1::2], strict=True):
        if _gap(large_case.share, small_case.share) > 1.0e-12:
            raise RuntimeError("the small-class share left the large-unit share")
        if large_case.gap > heat.helium_mismatch:
            raise RuntimeError(f"{large_case.name} misses the helium quotient match")
    if process.fuel_side_j <= 0.0 or process.nameplate_j <= 0.0:
        raise RuntimeError("an annual figure is not positive")
    return StressReport(
        power=(large, small, built, carnot),
        shutdown=shutdown,
        catalog_ceiling=ceiling.catalog_ceiling,
        carnot=process.carnot,
        reload_tonnes=RELOAD_TONNES,
        nameplate_j=process.nameplate_j,
        fuel_side_j=process.fuel_side_j,
        companion_ratio=process.fuel_side_j / process.nameplate_j,
        discharge_over_ranked=process.discharge_over_ranked,
        stack_factor=process.stack_factor,
        fabrication_cap_percent=safety.fabrication_cap_percent,
        held_k=_HELD_K,
        helium_mismatch=heat.helium_mismatch,
        nuclear_scalar=heat.nuclear_scalar,
        thermodynamics_scalar=heat.thermodynamics_scalar,
        gamma_fsot=heat.gamma_fsot,
        valve_beta=heat.valve_beta,
        pin_prefix=heat.pin_prefix,
        digest_prefix=heat.digest_prefix,
    )


def _power_lines(case: PowerCase) -> list[str]:
    lines = [
        f"The {case.name} electric boundary is `{_num(case.electric_w)}` W.",
        f"The operating heat is `{_num(case.thermal_w)}` W.",
        f"The condenser duty is `{_num(case.waste_w)}` W.",
        f"The heat dumped per kilowatt-hour is `{_num(case.mj_per_kwh)}` MJ.",
        f"The core volume is `{_num(case.volume_m3)}` cubic metres.",
        f"The heat source is the {case.source}.",
    ]
    if case.gap is None:
        lines.append("This steam figure stays a page witness.")
        return lines
    lines.append(f"The relative gap is `{_num(case.gap)}`.")
    lines.append("The relative gap meets the helium quotient match.")
    return lines


def _shutdown_lines(case: ShutdownCase) -> list[str]:
    return [
        f"The {case.name} cooling load is `{_num(case.cooling_w)}` W.",
        f"That load is `{_num(case.share)}` of the running condenser duty `{_num(case.waste_w)}` W.",
        f"The cooling-load gap against the 1952 fit is `{_num(case.gap)}`.",
    ]


def render() -> str:
    report = build_stress()
    lines = [
        "# Fission stress",
        "",
        "This score runs the light-water plant already on file.",
        "The electric boundary stays the published nameplate.",
        "The demonstrated operating heat is the FSOT heat ratio times that nameplate.",
        "The chain stays at k_eff of 1.",
        "The core volume is the operating heat divided by the filed power density.",
        "That volume is not a measured vessel.",
        "Every FSOT gap in these cases meets the helium quotient match.",
        "",
        "## Power",
        "",
        "Four power cases use the conditions already scored.",
    ]
    for case in report.power:
        lines.extend(_power_lines(case))
        lines.append("")
    lines.extend(
        [
            "The condenser duty falls as the steam figure rises.",
            "The core volume falls with that duty at the same electric boundary.",
            f"The temperature Carnot is `{_num(report.carnot)}`.",
            f"The catalog ceiling stays `{_num(report.catalog_ceiling)}`.",
            "The catalog ceiling stays above the temperature Carnot.",
            "This score leaves that ceiling unrun.",
            "The temperature Carnot is the last steam case.",
            "",
            "## Shutdown",
            "",
            "Fission is stopped in the shutdown cases.",
            "The cooling load is the FSOT fraction times the FSOT operating heat.",
            "The share is that load against the running full-power condenser duty.",
            "The comparison assumes the condenser circulation the cooling page requires.",
            "A loss of that circulation has no second capacity on file.",
            "This score assigns those systems no watt figure.",
            "",
        ]
    )
    for case in report.shutdown:
        lines.extend(_shutdown_lines(case))
    lines.extend(
        [
            "",
            "Every scored share is below one.",
            "The one-day share is below the one-second share.",
            "The two-day share is below the one-day share.",
            "The small-class share matches the large-unit share.",
            "Decay heat is a cooling load.",
            "It is not a second electricity product.",
            "",
            "## Fuel",
            "",
            f"The page reload is `{_num(report.reload_tonnes)}` tonnes of enriched uranium.",
            f"The fuel-side year is `{_num(report.fuel_side_j)}` J.",
            f"The nameplate year is `{_num(report.nameplate_j)}` J.",
            f"The fuel-side year over the nameplate year is `{_num(report.companion_ratio)}`.",
            "The fuel-side year and the nameplate year stay the page's round companions.",
            "This run keeps both figures.",
            f"The discharge stands at `{_num(report.discharge_over_ranked)}` of the ranked pure U-235 density.",
            "The operating fuel load is that discharge.",
            f"The 60 GWd/t stack factor is `{_num(report.stack_factor)}`.",
            "That stack is not a stressed case.",
            f"The fabrication cap remains `{_num(report.fabrication_cap_percent)}` percent U-235.",
            "The fast row stays on its own reference unit.",
            "The fast row's filed barrier is the sodium system.",
            "The small light-water class does not retire that barrier.",
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
            f"The helium quotient match is `{_num(report.helium_mismatch)}`.",
            "The FSOT figures in this run are the heat ratio and the three cooling fractions.",
            "The page efficiencies, the size ratio, and the discharge ratio stay on their earlier scores.",
            "This run does not install a seed.",
            "",
            "The checked-in schedule is not rewritten.",
            "",
        ]
    )
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
