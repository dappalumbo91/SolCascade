"""Tighten the thermal plant ceiling from the temperatures the file already states.

The catalog ceiling stays the round figure in the YAML. The refined bound is
the Carnot value of the hot and cold temperatures written in that same
source. A foundational seed is ranked beside each efficiency and is not
multiplied into the joules. The checked-in pulse is not rewritten.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_paths import SeedGap, _nearest_seeds
from solcascade.production.plants import FAST_ID, THERMAL_ID, build_plants
from solcascade.production.refine import catalog_record, project_root

_TEMPERATURES = re.compile(r"near (\d+(?:\.\d+)?) K rejected near (\d+(?:\.\d+)?) K")
_REVIEW_GAP = 1.0e-3


@dataclass(frozen=True)
class CeilingReport:
    hot_k: float
    cold_k: float
    demonstrated: float
    carnot: float
    catalog_ceiling: float
    density_j: float
    demonstrated_j: float
    carnot_j: float
    catalog_j: float
    carnot_step_j: float
    catalog_margin_j: float
    fast_reference_unit: str
    fast_ceiling: float
    demonstrated_seed: SeedGap
    carnot_seed: SeedGap
    fast_seed: SeedGap
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_ceiling.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def temperature_carnot(source: str) -> tuple[float, float, float]:
    """Hot kelvin, cold kelvin, and the Carnot bound 1 - cold/hot."""
    match = _TEMPERATURES.search(source)
    if match is None:
        raise RuntimeError("the ceiling source does not state a rejected temperature pair")
    hot = float(match.group(1))
    cold = float(match.group(2))
    if cold <= 0.0 or hot <= cold:
        raise RuntimeError("the stated temperatures are not a heat engine")
    return hot, cold, 1.0 - cold / hot


def build_ceiling() -> CeilingReport:
    plants = build_plants()
    thermal = plants.thermal
    fast = plants.fast
    thermal_record = catalog_record(THERMAL_ID)
    fast_record = catalog_record(FAST_ID)
    if thermal_record.get("efficiency_ceiling_kind") != "carnot":
        raise RuntimeError("the thermal ceiling is not marked carnot")
    if fast_record.get("efficiency_ceiling_kind") != "carnot":
        raise RuntimeError("the fast ceiling is not marked carnot")
    hot, cold, carnot = temperature_carnot(str(thermal_record["efficiency_ceiling_source"]))
    if _TEMPERATURES.search(str(fast_record["efficiency_ceiling_source"])) is not None:
        raise RuntimeError("the fast ceiling source states a temperature pair")
    if carnot <= thermal.efficiency:
        raise RuntimeError("the temperature Carnot does not clear the demonstrated efficiency")
    if thermal.ceiling <= carnot:
        raise RuntimeError("the catalog ceiling is already at or below the temperature Carnot")
    if not 0.0 < carnot < 1.0:
        raise RuntimeError("the temperature Carnot is outside (0, 1)")
    density = thermal.energy_density_j
    demonstrated_j = thermal.delivered_j
    carnot_j = density * carnot
    catalog_j = thermal.ceiling_delivered_j
    if not demonstrated_j < carnot_j < catalog_j:
        raise RuntimeError("the refined bound does not sit between the demonstrated plant and the catalog ceiling")
    module, _digest = load_engine()
    demonstrated_seed, _next, _named, _ranked = _nearest_seeds(module, thermal.efficiency)
    carnot_seed, _carnot_next, _named_again, _ranked_again = _nearest_seeds(module, carnot)
    fast_seed, _fast_next, _named_fast, _ranked_fast = _nearest_seeds(module, fast.ceiling)
    for seed in (demonstrated_seed, carnot_seed, fast_seed):
        if seed.relative_gap <= _REVIEW_GAP:
            raise RuntimeError(f"{seed.name} meets an efficiency closely enough to review")
    return CeilingReport(
        hot_k=hot,
        cold_k=cold,
        demonstrated=thermal.efficiency,
        carnot=carnot,
        catalog_ceiling=thermal.ceiling,
        density_j=density,
        demonstrated_j=demonstrated_j,
        carnot_j=carnot_j,
        catalog_j=catalog_j,
        carnot_step_j=carnot_j - demonstrated_j,
        catalog_margin_j=catalog_j - carnot_j,
        fast_reference_unit=fast.reference_unit,
        fast_ceiling=fast.ceiling,
        demonstrated_seed=demonstrated_seed,
        carnot_seed=carnot_seed,
        fast_seed=fast_seed,
        nuclear_scalar=plants.nuclear_scalar,
        thermodynamics_scalar=plants.thermodynamics_scalar,
        gamma_fsot=plants.gamma_fsot,
        valve_beta=plants.valve_beta,
        pin_prefix=plants.pin_prefix,
        digest_prefix=plants.digest_prefix,
    )


def render() -> str:
    report = build_ceiling()
    lines = [
        "# Fission ceiling",
        "",
        f"The thermal ceiling source states a steam cycle near `{_num(report.hot_k)}` K rejected near `{_num(report.cold_k)}` K.",
        f"The Carnot bound of those temperatures is `{_num(report.carnot)}`.",
        f"The catalog ceiling stays `{_num(report.catalog_ceiling)}`.",
        "The catalog ceiling sits above that temperature Carnot.",
        "The refined bound is the temperature Carnot.",
        "",
        "| Bound | Efficiency | Delivered (J) |",
        "| --- | ---: | ---: |",
        f"| Demonstrated | {report.demonstrated:.6g} | {report.demonstrated_j:.6g} |",
        f"| Temperature Carnot | {report.carnot:.6g} | {report.carnot_j:.6g} |",
        f"| Catalog ceiling | {report.catalog_ceiling:.6g} | {report.catalog_j:.6g} |",
        "",
        f"The step from the demonstrated efficiency to the temperature Carnot is `{_num(report.carnot_step_j)}` J.",
        f"The catalog ceiling adds `{_num(report.catalog_margin_j)}` J beyond the temperature Carnot.",
        "That extra is the round-file margin.",
        "The refined electricity is the temperature Carnot.",
        "",
        f"The fast reference unit is {report.fast_reference_unit}.",
        "The fast ceiling source states no temperature pair.",
        f"The fast catalog ceiling stays `{_num(report.fast_ceiling)}`.",
        "The fast row is not tightened.",
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
        f"The nearest seed to the demonstrated efficiency is `{report.demonstrated_seed.name}` at `{_num(report.demonstrated_seed.value)}`.",
        f"Its relative gap to that efficiency is `{_num(report.demonstrated_seed.relative_gap)}`.",
        "The nearest seed is not installed as the demonstrated efficiency.",
        f"The nearest seed to the temperature Carnot is `{report.carnot_seed.name}` at `{_num(report.carnot_seed.value)}`.",
        f"Its relative gap to that Carnot bound is `{_num(report.carnot_seed.relative_gap)}`.",
        "The nearest seed is not installed as the temperature Carnot.",
        f"The nearest seed to the fast ceiling is `{report.fast_seed.name}` at `{_num(report.fast_seed.value)}`.",
        f"Its relative gap to that ceiling is `{_num(report.fast_seed.relative_gap)}`.",
        "The nearest seed is not installed as the fast ceiling.",
        "The closest of these seeds is still wider than one part in a thousand.",
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
