"""Score the cooling failure and the published small class.

Routine heat rejection is set by the steam-cycle efficiency. Decay heat after
shutdown is the nuclear failure the cooling page names, and that page states
no fraction for it. The small class is the association's published unit under
the large nameplate. A foundational seed is not multiplied into the joules.
The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_paths import SeedGap, _nearest_seeds
from solcascade.production.fission_process import (
    FABRICATION_CAP_PERCENT,
    FULL_POWER_HOURS,
    NAMEPLATE_MWE,
    build_process,
)
from solcascade.production.plants import FAST_ID
from solcascade.production.refine import catalog_record, project_root

COOLING_SOURCE = (
    "https://world-nuclear.org/information-library/current-and-future-generation/"
    "cooling-power-plants"
)
SMALL_CLASS_SOURCE = (
    "https://world-nuclear.org/information-library/nuclear-power-reactors/"
    "small-modular-reactors/small-modular-reactors"
)
# The cooling page's worked figure for plants now being built, inside its 34–36 percent band.
BUILT_EFFICIENCY = 0.36
# The small-class boundary on that page: typically less than 300 MWe, and the
# maturity note's "300 MWe or less."
SMALL_CLASS_MWE = 300.0
_MJ_PER_KWH = 3.6
_REVIEW_GAP = 1.0e-3


@dataclass(frozen=True)
class SafetyReport:
    efficiency: float
    built_efficiency: float
    carnot: float
    dump_demonstrated_mj: float
    dump_built_mj: float
    dump_carnot_mj: float
    nameplate_mwe: float
    small_mwe: float
    size_ratio: float
    large_volume_m3: float
    small_volume_m3: float
    linear_scale: float
    large_nameplate_j: float
    small_nameplate_j: float
    power_density_w_m3: float
    fabrication_cap_percent: float
    size_seed: SeedGap
    built_seed: SeedGap
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_safety.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def heat_dump_mj_per_kwh(efficiency: float) -> float:
    """Megajoules of waste heat per kilowatt-hour of electricity.

    The cooling page's identity: a steam plant dumps ``1/efficiency - 1``
    kilowatt-hours of heat for each kilowatt-hour delivered.
    """
    if not 0.0 < efficiency < 1.0:
        raise RuntimeError("the efficiency is outside (0, 1)")
    return (1.0 / efficiency - 1.0) * _MJ_PER_KWH


def build_safety() -> SafetyReport:
    process = build_process()
    if not process.efficiency < BUILT_EFFICIENCY < process.carnot:
        raise RuntimeError("the plants-being-built figure is not between the plant and the Carnot")
    if not 0.0 < SMALL_CLASS_MWE < NAMEPLATE_MWE:
        raise RuntimeError("the small class is not under the nameplate")
    fast = catalog_record(FAST_ID)
    if "Sodium" not in str(fast.get("engineering_barrier") or ""):
        raise RuntimeError("the fast barrier does not name sodium")
    dump_demonstrated = heat_dump_mj_per_kwh(process.efficiency)
    dump_built = heat_dump_mj_per_kwh(BUILT_EFFICIENCY)
    dump_carnot = heat_dump_mj_per_kwh(process.carnot)
    if not dump_carnot < dump_built < dump_demonstrated:
        raise RuntimeError("the heat dump does not fall as the efficiency rises")
    size_ratio = SMALL_CLASS_MWE / NAMEPLATE_MWE
    large_thermal_w = NAMEPLATE_MWE * 1.0e6 / process.efficiency
    large_volume = large_thermal_w / process.power_density_w_m3
    small_volume = large_volume * size_ratio
    hours = FULL_POWER_HOURS
    large_j = process.nameplate_j
    small_j = large_j * size_ratio
    if abs(large_j - NAMEPLATE_MWE * 1.0e6 * hours * 3600.0) > 1.0:
        raise RuntimeError("the nameplate year does not match the process score")
    module, _digest = load_engine()
    size_seed = _nearest_seeds(module, size_ratio)[0]
    built_seed = _nearest_seeds(module, BUILT_EFFICIENCY)[0]
    for seed in (size_seed, built_seed):
        if seed.relative_gap <= _REVIEW_GAP:
            raise RuntimeError(f"{seed.name} meets a safety figure closely enough to review")
    return SafetyReport(
        efficiency=process.efficiency,
        built_efficiency=BUILT_EFFICIENCY,
        carnot=process.carnot,
        dump_demonstrated_mj=dump_demonstrated,
        dump_built_mj=dump_built,
        dump_carnot_mj=dump_carnot,
        nameplate_mwe=NAMEPLATE_MWE,
        small_mwe=SMALL_CLASS_MWE,
        size_ratio=size_ratio,
        large_volume_m3=large_volume,
        small_volume_m3=small_volume,
        linear_scale=size_ratio ** (1.0 / 3.0),
        large_nameplate_j=large_j,
        small_nameplate_j=small_j,
        power_density_w_m3=process.power_density_w_m3,
        fabrication_cap_percent=FABRICATION_CAP_PERCENT,
        size_seed=size_seed,
        built_seed=built_seed,
        nuclear_scalar=process.nuclear_scalar,
        thermodynamics_scalar=process.thermodynamics_scalar,
        gamma_fsot=process.gamma_fsot,
        valve_beta=process.valve_beta,
        pin_prefix=process.pin_prefix,
        digest_prefix=process.digest_prefix,
    )


def render() -> str:
    report = build_safety()
    lines = [
        "# Fission safety",
        "",
        f"The cooling page is {COOLING_SOURCE}.",
        "Routine cooling for a steam plant is set by its thermal efficiency.",
        f"At the demonstrated efficiency the plant dumps `{_num(report.dump_demonstrated_mj)}` MJ of heat per kilowatt-hour.",
        "The page rounds that 33 percent case to 7.3 MJ.",
        f"Plants now being built are stated near `{_num(report.built_efficiency)}`.",
        f"At that figure the dump is `{_num(report.dump_built_mj)}` MJ per kilowatt-hour.",
        "The page rounds that 36 percent case to 6.4 MJ.",
        f"At the temperature Carnot the dump is `{_num(report.dump_carnot_mj)}` MJ per kilowatt-hour.",
        "The heat dumped per kilowatt-hour falls as the efficiency rises.",
        "",
        "## After shutdown",
        "",
        "Fission stops at shutdown, and decay heat continues.",
        "The cooling page states no decay-heat fraction, so no seed is ranked against it.",
        "The routine path is the main steam bypassing the turbine and dumping heat into the condenser.",
        "After the pressure drops, a residual heat removal system carries the load.",
        "Emergency core cooling is the page's system when primary cooling fails.",
        "If circulation is maintained, the page says the decay heat ceases to be a problem after a day or two.",
        "For large units, loss of power after an emergency shutdown is the decay-heat case the page flags.",
        "Decay heat is a cooling load.",
        "It is not a second electricity product.",
        "The delayed-neutron margin stays in the core.",
        "Loss of cooling is a heat-removal problem.",
        "",
        "## Smaller unit",
        "",
        f"The small-class page is {SMALL_CLASS_SOURCE}.",
        f"The large nameplate in the fuel-cycle balance is `{_num(report.nameplate_mwe)}` MWe.",
        f"The small class is `{_num(report.small_mwe)}` MWe or less.",
        f"That boundary is `{_num(report.size_ratio)}` of the nameplate.",
        "Near-term small designs are light-water reactors on the same low-enriched fuel.",
        f"The fabrication cap remains `{_num(report.fabrication_cap_percent)}` percent U-235.",
        f"At the filed power density `{_num(report.power_density_w_m3)}` W per cubic metre, the large core volume is `{_num(report.large_volume_m3)}` cubic metres.",
        f"The same density at the small-class boundary is `{_num(report.small_volume_m3)}` cubic metres.",
        "The core volume scales with the thermal power.",
        f"If the core keeps its shape, the linear scale is `{_num(report.linear_scale)}`.",
        f"A full-power year at the small-class boundary is `{_num(report.small_nameplate_j)}` J.",
        f"The large nameplate year is `{_num(report.large_nameplate_j)}` J.",
        "The joules per kilogram of fuel stay on the discharge.",
        "The small class is the published unit for a grid that does not host the large nameplate.",
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
        "",
        f"The nearest seed to the size ratio is `{report.size_seed.name}` at `{_num(report.size_seed.value)}`.",
        f"Its relative gap is `{_num(report.size_seed.relative_gap)}`.",
        f"The nearest seed to the plants-being-built efficiency is `{report.built_seed.name}` at `{_num(report.built_seed.value)}`.",
        f"Its relative gap is `{_num(report.built_seed.relative_gap)}`.",
        "The closest of these seeds is still wider than one part in a thousand.",
        "No seed is installed as a cooling figure or as the size ratio.",
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
