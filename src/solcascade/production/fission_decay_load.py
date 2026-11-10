"""Derive the decay cooling loads from the helium-tight FSOT fractions.

The fraction is the recorded expression. The operating heat is the filed
nameplate divided by the filed efficiency. The product is the cooling
load. The 1952 fit stays the calorimetry witness. The checked-in schedule
is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay import build_decay
from solcascade.production.fission_decay_open import build_decay_open
from solcascade.production.fission_decay_recorded import build_decay_recorded
from solcascade.production.fission_decay_seed import NamedGap
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.refine import project_root

_LOADS = ("one-second fraction", "one-day fraction", "two-day fraction")
_LIMITS = ("one-second limit", "one-day limit", "two-day limit")


@dataclass(frozen=True)
class DerivedLoad:
    label: str
    fraction: NamedGap
    fit_fraction: float
    watts: float
    fit_watts: float
    small_watts: float
    fit_small_watts: float
    limit: NamedGap
    fit_limit_fraction: float
    limit_watts: float
    fit_limit_watts: float


@dataclass(frozen=True)
class LoadReport:
    loads: tuple[DerivedLoad, ...]
    day_over_second: float
    two_over_day: float
    fit_day_over_second: float
    fit_two_over_day: float
    day_ratio_gap: float
    two_ratio_gap: float
    thermal_w: float
    small_thermal_w: float
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_decay_load.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _leaf(label: str, recorded, opened) -> NamedGap:
    if label in opened:
        return opened[label].steps[-1].dressing
    return recorded[label].dressing


def _gap(derived: float, witness: float) -> float:
    if witness == 0.0:
        raise RuntimeError("the witness load is zero")
    return abs(derived - witness) / abs(witness)


def build_decay_load() -> LoadReport:
    decay = build_decay()
    recorded = {figure.label: figure for figure in build_decay_recorded().figures}
    opened = {figure.label: figure for figure in build_decay_open().figures}
    module, _digest = load_engine()
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    if len(decay.points) != len(_LOADS):
        raise RuntimeError("the cooling schedule changed")
    loads: list[DerivedLoad] = []
    for label, limit_label, point in zip(_LOADS, _LIMITS, decay.points, strict=True):
        fraction = _leaf(label, recorded, opened)
        limit = _leaf(limit_label, recorded, opened)
        if fraction.relative_gap > helium or limit.relative_gap > helium:
            raise RuntimeError(f"the {label} is outside the helium match")
        watts = fraction.value * decay.thermal_w
        small = fraction.value * decay.small_thermal_w
        limit_watts = limit.value * decay.thermal_w
        if _gap(watts, point.year_watts) > helium:
            raise RuntimeError(f"the {label} load misses the helium match")
        if _gap(small, point.small_watts) > helium:
            raise RuntimeError(f"the {label} small-class load misses the helium match")
        if _gap(limit_watts, point.limit_watts) > helium:
            raise RuntimeError(f"the {limit_label} misses the helium match")
        loads.append(
            DerivedLoad(
                label=label,
                fraction=fraction,
                fit_fraction=point.year_fraction,
                watts=watts,
                fit_watts=point.year_watts,
                small_watts=small,
                fit_small_watts=point.small_watts,
                limit=limit,
                fit_limit_fraction=point.limit_fraction,
                limit_watts=limit_watts,
                fit_limit_watts=point.limit_watts,
            )
        )
    second, day, two = loads
    day_over_second = day.watts / second.watts
    two_over_day = two.watts / day.watts
    day_ratio_gap = _gap(day_over_second, decay.day_over_second)
    two_ratio_gap = _gap(two_over_day, decay.two_over_day)
    if day_ratio_gap > helium or two_ratio_gap > helium:
        raise RuntimeError("a derived load ratio misses the helium match")
    if not 0.0 < day_over_second < two_over_day < 1.0:
        raise RuntimeError("the derived early drop is not steeper than the later drop")
    return LoadReport(
        loads=tuple(loads),
        day_over_second=day_over_second,
        two_over_day=two_over_day,
        fit_day_over_second=decay.day_over_second,
        fit_two_over_day=decay.two_over_day,
        day_ratio_gap=day_ratio_gap,
        two_ratio_gap=two_ratio_gap,
        thermal_w=decay.thermal_w,
        small_thermal_w=decay.small_thermal_w,
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=decay.nuclear_scalar,
        thermodynamics_scalar=decay.thermodynamics_scalar,
        gamma_fsot=decay.gamma_fsot,
        valve_beta=decay.valve_beta,
        pin_prefix=decay.pin_prefix,
        digest_prefix=decay.digest_prefix,
    )


def _load_lines(load: DerivedLoad) -> list[str]:
    return [
        f"The {load.label} is `{load.fraction.name}`.",
        f"The FSOT fraction is `{_num(load.fraction.value)}`.",
        f"The 1952 fraction is `{_num(load.fit_fraction)}`.",
        f"The fraction gap is `{_num(load.fraction.relative_gap)}`.",
        f"The derived cooling load is `{_num(load.watts)}` W.",
        f"The 1952 fit load is `{_num(load.fit_watts)}` W.",
        f"The derived small-class load is `{_num(load.small_watts)}` W.",
        f"The 1952 small-class load is `{_num(load.fit_small_watts)}` W.",
        f"The long-operation fraction is `{load.limit.name}`.",
        f"The derived limit load is `{_num(load.limit_watts)}` W.",
        f"The 1952 limit load is `{_num(load.fit_limit_watts)}` W.",
        f"The limit gap is `{_num(load.limit.relative_gap)}`.",
    ]


def render() -> str:
    report = build_decay_load()
    lines = [
        "# Fission decay load",
        "",
        "The cooling load is the helium-tight FSOT fraction times the filed operating heat.",
        "The filed operating heat is the nameplate divided by the filed efficiency.",
        "That heat is not an FSOT scalar.",
        "The 1952 fit remains the calorimetry witness.",
        "Every scored gap meets the helium quotient match.",
        f"alpha in the fraction is the fine-structure leaf `{_num(report.alpha)}`.",
        f"The helium quotient match is `{_num(report.helium_mismatch)}`.",
        f"Operating heat on the large nameplate is `{_num(report.thermal_w)}` W.",
        f"Operating heat on the small-class boundary is `{_num(report.small_thermal_w)}` W.",
        "",
        "## Loads",
        "",
    ]
    for load in report.loads:
        lines.extend(_load_lines(load))
        lines.append("")
    lines.extend(
        [
            "## Ratios",
            "",
            f"The derived one-day load over the derived one-second load is `{_num(report.day_over_second)}`.",
            f"The 1952 ratio is `{_num(report.fit_day_over_second)}`.",
            f"The relative gap is `{_num(report.day_ratio_gap)}`.",
            f"The derived two-day load over the derived one-day load is `{_num(report.two_over_day)}`.",
            f"The 1952 ratio is `{_num(report.fit_two_over_day)}`.",
            f"The relative gap is `{_num(report.two_ratio_gap)}`.",
            "",
            "## Ledger",
            "",
            "The derived loads are the cooling loads.",
            "Decay heat is a cooling load.",
            "It is not a second electricity product.",
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
