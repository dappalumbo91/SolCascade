"""Score the neutron, gamma, and delayed balance on the U-235 kilogram.

The balance is the cited page's 30 MeV addendum. It already sits inside the
200 MeV ledger, so the Rankine product already includes it. The same
foundational products are ranked beside that share and beside the delayed
neutron count. A loose seed stays a witness. The checked-in pulse is not
rewritten.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_observables import _adopted_alpha
from solcascade.production.fission_paths import (
    ADDED_MEV,
    LEDGER_MEV,
    SeedGap,
    _nearest_seeds,
    build_paths,
)
from solcascade.production.refine import catalog_record, project_root

BALANCE_SHARE = ADDED_MEV / LEDGER_MEV
_DELAYED_PERCENT = (
    ("U-235", 0.66),
    ("U-233", 0.27),
    ("Pu-239", 0.23),
)


@dataclass(frozen=True)
class DelayedProportion:
    """A page count of delayed neutrons, not an energy share."""

    nucleus: str
    percent: float
    proportion: float


@dataclass(frozen=True)
class BalanceReport:
    density_j: float
    efficiency: float
    ceiling: float
    delivered_j: float
    balance_share: float
    balance_inventory_j: float
    balance_electricity_j: float
    ceiling_headroom_j: float
    upper_bound_increment_j: float
    governing_source: str
    fragment_inventory_j: float
    delayed: tuple[DelayedProportion, ...]
    balance_nearest: SeedGap
    complement: SeedGap
    delayed_nearest: SeedGap
    fine_structure: float
    fine_structure_gap: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_balance.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def delayed_proportions() -> tuple[DelayedProportion, ...]:
    rows = []
    for nucleus, percent in _DELAYED_PERCENT:
        if not 0.0 < percent < 1.0:
            raise RuntimeError(f"the {nucleus} delayed percent is not a small count")
        rows.append(DelayedProportion(nucleus, percent, percent / 100.0))
    return tuple(rows)


def _no_shield_key() -> None:
    data = project_root() / "src" / "solcascade" / "energy_baseline" / "data"
    for path in sorted(data.glob("*.yaml")):
        record = catalog_record(path.stem)
        for key in record:
            if "shield" in str(key).lower():
                raise RuntimeError(f"{record.get('id')} files a shield field")


def build_balance() -> BalanceReport:
    if not 0.0 < BALANCE_SHARE < 1.0:
        raise RuntimeError("the balance share is outside (0, 1)")
    paths = build_paths()
    thermal = paths.thermal
    fragment_share = float(paths.fragment_share)
    if fragment_share + BALANCE_SHARE != 1.0:
        raise RuntimeError("the fragment share and the balance share do not add to one")
    density = thermal.energy_density_j
    inventory = BALANCE_SHARE * density
    fragment_inventory = paths.kinetic_inventory_j
    if abs(inventory + fragment_inventory - density) > 1.0e-6 * density:
        raise RuntimeError("the fragment inventory and the balance do not add to the density")
    electricity = thermal.efficiency * inventory
    if electricity >= thermal.delivered_j:
        raise RuntimeError("the balance electricity meets the Rankine product")
    headroom = (thermal.ceiling - thermal.efficiency) * density
    upper = (1.0 - thermal.efficiency) * inventory
    if headroom <= upper:
        raise RuntimeError("the file-ceiling headroom does not clear the balance bound")
    if inventory >= density:
        raise RuntimeError("the balance inventory meets the whole density")
    _no_shield_key()
    delayed = delayed_proportions()
    thermal_delayed = delayed[0]
    if thermal_delayed.nucleus != "U-235":
        raise RuntimeError("the thermal delayed proportion is not U-235")
    module, _digest = load_engine()
    balance_nearest, _following, _named, _ranked = _nearest_seeds(module, BALANCE_SHARE)
    delayed_nearest, _delayed_next, _named_again, _ranked_again = _nearest_seeds(
        module, thermal_delayed.proportion
    )
    if balance_nearest.relative_gap <= 0.0 or delayed_nearest.relative_gap <= 0.0:
        raise RuntimeError("a seed reproduces a page figure")
    complement_value = 1.0 - paths.nearest.value
    complement_gap = abs(complement_value - BALANCE_SHARE) / BALANCE_SHARE
    if complement_gap <= 0.0:
        raise RuntimeError("the complement reproduces the balance")
    if complement_gap >= balance_nearest.relative_gap:
        raise RuntimeError("the complement is not closer than the balance grammar")
    alpha = float(_adopted_alpha(module))
    if not math.isfinite(alpha) or alpha <= 0.0:
        raise RuntimeError("the fine-structure leaf is not a positive finite number")
    alpha_gap = abs(alpha - thermal_delayed.proportion) / thermal_delayed.proportion
    if alpha_gap <= 0.0:
        raise RuntimeError("the fine-structure leaf reproduces the delayed proportion")
    return BalanceReport(
        density_j=density,
        efficiency=thermal.efficiency,
        ceiling=thermal.ceiling,
        delivered_j=thermal.delivered_j,
        balance_share=BALANCE_SHARE,
        balance_inventory_j=inventory,
        balance_electricity_j=electricity,
        ceiling_headroom_j=headroom,
        upper_bound_increment_j=upper,
        governing_source=paths.governing_source,
        fragment_inventory_j=fragment_inventory,
        delayed=delayed,
        balance_nearest=balance_nearest,
        complement=SeedGap(
            name=f"1-{paths.nearest.name}",
            value=complement_value,
            relative_gap=complement_gap,
        ),
        delayed_nearest=delayed_nearest,
        fine_structure=alpha,
        fine_structure_gap=alpha_gap,
        nuclear_scalar=paths.nuclear_scalar,
        thermodynamics_scalar=paths.thermodynamics_scalar,
        gamma_fsot=paths.gamma_fsot,
        valve_beta=paths.valve_beta,
        pin_prefix=paths.pin_prefix,
        digest_prefix=paths.digest_prefix,
    )


def _delayed_row(row: DelayedProportion) -> str:
    return f"| {row.nucleus} | {row.percent:.6g} |"


def render() -> str:
    report = build_balance()
    nearest = report.balance_nearest
    complement = report.complement
    delayed = report.delayed_nearest
    lines = [
        "# Fission balance",
        "",
        f"The thermal governing source is {report.governing_source}.",
        "That page assigns the balance, after the fragment kinetic energy, to gamma rays and to the kinetic energy of the neutrons.",
        f"Neutron, gamma, and delayed releases add about `{_num(ADDED_MEV)}` MeV on a total of about `{_num(LEDGER_MEV)}` MeV.",
        f"That balance is the share `{_num(report.balance_share)}`.",
        "The catalog density stays on the 200 MeV ledger.",
        f"The balance inventory on the ranked kilogram is `{_num(report.balance_inventory_j)}` J.",
        f"The fragment inventory on the ranked kilogram is `{_num(report.fragment_inventory_j)}` J.",
        "The fragment inventory and the balance inventory add to the ranked density.",
        "The balance is already inside the ranked density.",
        "The Rankine product already includes it.",
        f"Its electricity at the demonstrated efficiency is `{_num(report.balance_electricity_j)}` J, already inside the Rankine product of `{_num(report.delivered_j)}` J.",
        "Scoring that electricity again would count the same joules twice.",
        "",
        "## Shield and control",
        "",
        "The catalog files no shield-recovery efficiency.",
        "Neutrons and gamma rays that stop outside the coolant are a shield load.",
        "The shield stays a shield.",
        "The page gives these delayed-neutron proportions.",
        "",
        "| Nucleus | Delayed neutrons (percent) |",
        "| --- | ---: |",
        *(_delayed_row(row) for row in report.delayed),
        "",
        "Those proportions count neutrons.",
        "They are not an energy share of the ledger.",
        "Delayed neutrons enable control of the chain.",
        "The control margin stays in the core.",
        "",
        "## Ceiling",
        "",
        f"The file-ceiling headroom on the whole ledger is `{_num(report.ceiling_headroom_j)}` J.",
        "An upper bound that converts the balance at one, and leaves the rest of the ledger at the demonstrated efficiency, "
        f"adds `{_num(report.upper_bound_increment_j)}` J.",
        "That bound is not a filed plant.",
        "The file-ceiling headroom is larger than that bound.",
        "The refinement inside the operating plant is the file efficiency ceiling.",
        "",
        "## FSOT record",
        "",
        f"Pin `{report.pin_prefix}` matches the loaded engine prefix `{report.digest_prefix}`.",
        "The scalars below are dimensionless.",
        "A domain scalar is not a plant efficiency.",
        "A domain scalar is not the balance share.",
        "The nuclear scalar stays out of the fission joules.",
        f"Nuclear_Physics: `{_num(report.nuclear_scalar)}`.",
        f"Thermodynamics: `{_num(report.thermodynamics_scalar)}`.",
        f"gamma_FSOT: `{_num(report.gamma_fsot)}`.",
        "gamma_FSOT does not multiply the plant.",
        "The engine valve stays out of the fission joules.",
        f"Engine valve: `{_num(report.valve_beta)}`.",
        "",
        "The balance uses the same foundational products as the fragment share.",
        f"The balance grammar's nearest seed is `{nearest.name}` at `{_num(nearest.value)}`.",
        f"Its relative gap to the balance is `{_num(nearest.relative_gap)}`.",
        f"The complement of the fragment nearest is `{complement.name}` at `{_num(complement.value)}`.",
        f"Its relative gap to the balance is `{_num(complement.relative_gap)}`.",
        "The complement sits closer than the balance grammar.",
        "The complement is not installed as the balance.",
        f"The nearest seed to the U-235 delayed proportion is `{delayed.name}` at `{_num(delayed.value)}`.",
        f"Its relative gap to that proportion is `{_num(delayed.relative_gap)}`.",
        f"The fine-structure leaf is `{_num(report.fine_structure)}`.",
        f"Its relative gap to the U-235 delayed proportion is `{_num(report.fine_structure_gap)}`.",
        "The fine-structure leaf is not the delayed-neutron proportion.",
        "The page figures stay round figures.",
        "No seed is installed as the balance or as the delayed proportion.",
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
