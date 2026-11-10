"""Walk the light-water figures from enrichment through discharge.

The catalog scores a pure U-235 kilogram, a highly enriched enrichment
estimate, and a heavy-metal discharge. The World Nuclear Association's
fuel-cycle overview and fuel-fabrication pages use the same 45 GWd/t,
33 percent, and 50 kWh/SWU. Their other round figures stay attributed to
those pages. A foundational seed is not multiplied into the joules. The
checked-in pulse is not rewritten.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_ceiling import temperature_carnot
from solcascade.production.fission_paths import SeedGap, _nearest_seeds
from solcascade.production.plants import THERMAL_ID, build_plants, discharge_j_per_kg
from solcascade.production.refine import catalog_record, project_root

MATERIAL_BALANCE_SOURCE = (
    "https://world-nuclear.org/information-library/nuclear-fuel-cycle/introduction/"
    "nuclear-fuel-cycle-overview"
)
FABRICATION_SOURCE = (
    "https://world-nuclear.org/information-library/nuclear-fuel-cycle/"
    "conversion-enrichment-and-fabrication/fuel-fabrication"
)
# Typical annual balance on the overview, at the catalog's 45 GWd/t and 33 percent.
ENRICHED_PERCENT = 4.5
ANNUAL_SWU = 182_000.0
RELOAD_TONNES = 24.3
NAMEPLATE_MWE = 1000.0
FULL_POWER_HOURS = 8760.0
# Fabrication page band, anchored by the overview's 60 GWd/t atomic figure.
HIGHER_BURNUP_GWD = 60.0
ATOMIC_PERCENT = 6.5
FABRICATION_CAP_PERCENT = 5.0
_JOULE_PER_KWH = 3_600_000.0
_REVIEW_GAP = 1.0e-3
_SWU = re.compile(r"(\d+(?:\.\d+)?) SWU")
_KWH = re.compile(r"(\d+(?:\.\d+)?) kWh/SWU")


@dataclass(frozen=True)
class ProcessReport:
    swu_per_kg_heu: float
    kwh_per_swu: float
    enrichment_j_per_kg_heu: float
    discharge_gwd: float
    discharge_j_per_kg: float
    ranked_density_j: float
    discharge_over_ranked: float
    efficiency: float
    carnot: float
    heavy_electric_j: float
    heavy_carnot_j: float
    higher_discharge_j: float
    higher_electric_j: float
    higher_carnot_j: float
    burnup_step: float
    carnot_step: float
    stack_factor: float
    higher_over_ranked: float
    stack_over_ranked_electric: float
    swu_per_kg_reload: float
    enrichment_j_per_kg_reload: float
    heu_swu_over_reload: float
    reload_output_over_input: float
    nameplate_j: float
    fuel_side_j: float
    power_density_w_m3: float
    discharge_seed: SeedGap
    enrichment_seed: SeedGap
    atomic_seed: SeedGap
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_process.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def enrichment_estimate(basis: str) -> tuple[float, float]:
    swu = _SWU.search(basis)
    kwh = _KWH.search(basis)
    if swu is None or kwh is None:
        raise RuntimeError("the input basis does not state SWU and kWh per SWU")
    return float(swu.group(1)), float(kwh.group(1))


def build_process() -> ProcessReport:
    plants = build_plants()
    thermal = plants.thermal
    record = catalog_record(THERMAL_ID)
    if "k_eff" not in str(record["governing_equation"]):
        raise RuntimeError("the governing equation does not hold k_eff")
    swu, kwh = enrichment_estimate(str(record["input_energy_basis"]))
    enrichment_j = swu * kwh * _JOULE_PER_KWH
    if abs(enrichment_j - thermal.input_j) > 1.0:
        raise RuntimeError("the SWU product does not match the filed enrichment joules")
    gigawatt_days, discharge = discharge_j_per_kg(str(record.get("notes") or ""))
    if gigawatt_days != 45.0:
        raise RuntimeError("the notes discharge is not the 45 GWd material balance")
    if abs(thermal.efficiency - 0.33) > 1.0e-12:
        raise RuntimeError("the plant efficiency is not the 33 percent material balance")
    if kwh != 50.0:
        raise RuntimeError("the filed kWh per SWU is not the material-balance rate")
    hot_source = str(record["efficiency_ceiling_source"])
    _hot, _cold, carnot = temperature_carnot(hot_source)
    if not thermal.efficiency < carnot < thermal.ceiling:
        raise RuntimeError("the temperature Carnot is not the refined ceiling")
    density = thermal.energy_density_j
    higher_discharge = discharge * (HIGHER_BURNUP_GWD / gigawatt_days)
    heavy_electric = discharge * thermal.efficiency
    heavy_carnot = discharge * carnot
    higher_electric = higher_discharge * thermal.efficiency
    higher_carnot = higher_discharge * carnot
    burnup_step = HIGHER_BURNUP_GWD / gigawatt_days
    carnot_step = carnot / thermal.efficiency
    reload_kg = RELOAD_TONNES * 1000.0
    swu_per_kg = ANNUAL_SWU / reload_kg
    enrichment_per_kg = swu_per_kg * kwh * _JOULE_PER_KWH
    nameplate_j = NAMEPLATE_MWE * 1.0e6 * FULL_POWER_HOURS * 3600.0
    fuel_side_j = reload_kg * heavy_electric
    power_density = float(record["power_density_w"])
    if power_density != 1.0e8:
        raise RuntimeError("the core power density is not the filed 100 W per cubic centimetre")
    module, _digest = load_engine()
    fractions = (
        discharge / density,
        ENRICHED_PERCENT / 100.0,
        ATOMIC_PERCENT / 100.0,
    )
    seeds = tuple(_nearest_seeds(module, fraction)[0] for fraction in fractions)
    for seed in seeds:
        if seed.relative_gap <= _REVIEW_GAP:
            raise RuntimeError(f"{seed.name} meets a process fraction closely enough to review")
    ranked_electric = density * thermal.efficiency
    return ProcessReport(
        swu_per_kg_heu=swu,
        kwh_per_swu=kwh,
        enrichment_j_per_kg_heu=enrichment_j,
        discharge_gwd=gigawatt_days,
        discharge_j_per_kg=discharge,
        ranked_density_j=density,
        discharge_over_ranked=discharge / density,
        efficiency=thermal.efficiency,
        carnot=carnot,
        heavy_electric_j=heavy_electric,
        heavy_carnot_j=heavy_carnot,
        higher_discharge_j=higher_discharge,
        higher_electric_j=higher_electric,
        higher_carnot_j=higher_carnot,
        burnup_step=burnup_step,
        carnot_step=carnot_step,
        stack_factor=burnup_step * carnot_step,
        higher_over_ranked=higher_discharge / density,
        stack_over_ranked_electric=higher_carnot / ranked_electric,
        swu_per_kg_reload=swu_per_kg,
        enrichment_j_per_kg_reload=enrichment_per_kg,
        heu_swu_over_reload=swu / swu_per_kg,
        reload_output_over_input=heavy_electric / enrichment_per_kg,
        nameplate_j=nameplate_j,
        fuel_side_j=fuel_side_j,
        power_density_w_m3=power_density,
        discharge_seed=seeds[0],
        enrichment_seed=seeds[1],
        atomic_seed=seeds[2],
        nuclear_scalar=plants.nuclear_scalar,
        thermodynamics_scalar=plants.thermodynamics_scalar,
        gamma_fsot=plants.gamma_fsot,
        valve_beta=plants.valve_beta,
        pin_prefix=plants.pin_prefix,
        digest_prefix=plants.digest_prefix,
    )


def render() -> str:
    report = build_process()
    lines = [
        "# Fission process",
        "",
        "The operating implementation is a light-water core. Water is the moderator and the coolant.",
        "The chain is held at k_eff of 1.",
        f"The core power density is `{_num(report.power_density_w_m3)}` W per cubic metre.",
        "That power density does not set the joules per kilogram.",
        "The fragments stop in the solid fuel, and the Rankine cycle takes the heat.",
        "",
        "## Three loads",
        "",
        f"The ranked density is `{_num(report.ranked_density_j)}` J per kilogram of pure U-235.",
        f"The enrichment basis is `{_num(report.swu_per_kg_heu)}` SWU at `{_num(report.kwh_per_swu)}` kWh per SWU.",
        f"That product is `{_num(report.enrichment_j_per_kg_heu)}` J per kilogram of highly enriched uranium.",
        f"The notes cite `{_num(report.discharge_gwd)}` GWd per tonne of heavy metal.",
        f"That discharge is `{_num(report.discharge_j_per_kg)}` J per kilogram of heavy metal.",
        f"Per kilogram, the discharge stands at `{_num(report.discharge_over_ranked)}` of the ranked density.",
        "A heavy-metal kilogram, a highly enriched kilogram, and a pure U-235 kilogram are three loads.",
        "",
        "## Normal reload",
        "",
        f"The material balance is {MATERIAL_BALANCE_SOURCE}.",
        "That page's typical year is a 1000 MWe light-water reactor at 4.5 percent enrichment, 45 GWd/t, and 33 percent.",
        f"The overview spends `{_num(ANNUAL_SWU)}` SWU on `{_num(RELOAD_TONNES)}` tonnes of enriched uranium.",
        f"That is `{_num(report.swu_per_kg_reload)}` SWU per kilogram.",
        f"At the same `{_num(report.kwh_per_swu)}` kWh per SWU, the enrichment electricity is `{_num(report.enrichment_j_per_kg_reload)}` J per kilogram.",
        f"The highly enriched estimate is `{_num(report.heu_swu_over_reload)}` times that separative work per kilogram.",
        f"Discharge electricity over that reload enrichment is `{_num(report.reload_output_over_input)}`.",
        f"The nameplate year at `{_num(NAMEPLATE_MWE)}` MWe for `{_num(FULL_POWER_HOURS)}` hours is `{_num(report.nameplate_j)}` J.",
        f"The fuel-side year, reload mass times discharge electricity, is `{_num(report.fuel_side_j)}` J.",
        "Those two annual figures are the page's round companions.",
        "",
        "## Where the output can still move",
        "",
        f"Discharge electricity at the demonstrated efficiency is `{_num(report.heavy_electric_j)}` J per kilogram of heavy metal.",
        f"The temperature Carnot on that same discharge is `{_num(report.heavy_carnot_j)}` J.",
        f"The Carnot step on the discharge is `{_num(report.carnot_step)}`.",
        f"The fabrication page is {FABRICATION_SOURCE}.",
        "That page carries burn-up from about 40 GWd/t past 60 GWd/t as enrichment rises toward 5 percent U-235.",
        "The overview anchors 60 GWd/t at about 6.5 percent of the uranium atoms.",
        "The normal fabrication limit is 5 percent U-235.",
        "The higher burnup band stops at that limit.",
        f"The 60 GWd/t discharge is `{_num(report.higher_discharge_j)}` J per kilogram.",
        f"At the demonstrated efficiency that discharge delivers `{_num(report.higher_electric_j)}` J.",
        f"The burnup step, 60 over 45, is `{_num(report.burnup_step)}`.",
        f"The temperature Carnot on the 60 GWd/t discharge is `{_num(report.higher_carnot_j)}` J.",
        f"The two steps together are `{_num(report.stack_factor)}`.",
        "That stack is not a filed plant.",
        f"The 60 GWd/t discharge stands at `{_num(report.higher_over_ranked)}` of the ranked density.",
        f"The stacked electricity stands at `{_num(report.stack_over_ranked_electric)}` of the ranked electricity at the demonstrated efficiency.",
        "The fast row remains the filed route for the actinides a thermal discharge leaves behind.",
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
        f"The nearest seed to the discharge ratio is `{report.discharge_seed.name}` at `{_num(report.discharge_seed.value)}`.",
        f"Its relative gap is `{_num(report.discharge_seed.relative_gap)}`.",
        f"The nearest seed to the 4.5 percent enrichment is `{report.enrichment_seed.name}` at `{_num(report.enrichment_seed.value)}`.",
        f"Its relative gap is `{_num(report.enrichment_seed.relative_gap)}`.",
        f"The nearest seed to the 6.5 percent atomic burn-up is `{report.atomic_seed.name}` at `{_num(report.atomic_seed.value)}`.",
        f"Its relative gap is `{_num(report.atomic_seed.relative_gap)}`.",
        "No seed is installed as a process fraction.",
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
