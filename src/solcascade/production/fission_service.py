"""Name the decay fractions at the two times the residual-heat section prints.

Four hours is the service entry. Twenty hours is the end of the published
cooldown, where the cold-end duty applies. The crossing is that duty's
share of the filed year. The 1952 fit remains the witness. The hotter
entry rate stays unprinted. The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay import decay_fraction
from solcascade.production.fission_decay_dress import _base_pieces, _dressings
from solcascade.production.fission_decay_heat import build_decay_heat
from solcascade.production.fission_decay_seed import NamedGap, _five_pieces, _nearest, _size
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.fission_process import FULL_POWER_HOURS
from solcascade.production.fission_removal import _cross_hours, build_removal
from solcascade.production.refine import project_root

_SERVICE_S = 4.0 * 3600.0
_COOLDOWN_S = 20.0 * 3600.0
_SERVICE_ANCHOR = "sqrt(e)*e^-5"
_SERVICE_STEPS = (
    "(sqrt(e)*e^-5)*(1-alpha*sqrt(gamma))",
    "((sqrt(e)*e^-5)*(1-alpha*sqrt(gamma)))*(1+alpha*e^-4)",
    "(((sqrt(e)*e^-5)*(1-alpha*sqrt(gamma)))*(1+alpha*e^-4))*(1-alpha*e^-6)",
)
_COOLDOWN_ANCHOR = "pi^-4-pi^-5"
_COOLDOWN_STEPS = (
    "(pi^-4-pi^-5)*(1-alpha*pi^-4)",
    "((pi^-4-pi^-5)*(1-alpha*pi^-4))*(1-alpha*pi^-6)",
)
_CROSS_ANCHOR = "phi^-2*gamma^6"
_CROSS_STEPS = (
    "(phi^-2*gamma^6)*(1+alpha*G^-2)",
    "((phi^-2*gamma^6)*(1+alpha*G^-2))*(1-alpha*gamma^6)",
    "(((phi^-2*gamma^6)*(1+alpha*G^-2))*(1-alpha*gamma^6))*(1+alpha*e^-6)",
)
_COUNT = 116
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4


@dataclass(frozen=True)
class ServiceStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class ServiceReport:
    service_anchor: NamedGap
    service_steps: tuple[ServiceStep, ...]
    service_target: float
    service_w: float
    service_fit_w: float
    service_small_w: float
    service_fit_small_w: float
    service_watt_gap: float
    cooldown_anchor: NamedGap
    cooldown_steps: tuple[ServiceStep, ...]
    cooldown_target: float
    cooldown_w: float
    cooldown_fit_w: float
    cooldown_small_w: float
    cooldown_fit_small_w: float
    cooldown_watt_gap: float
    cross_anchor: NamedGap
    cross_steps: tuple[ServiceStep, ...]
    cross_target: float
    cross_hours: float
    fit_hours: float
    duty_fraction: float
    derived_duty: float
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_service.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _gap(derived: float, witness: float) -> float:
    if witness == 0.0:
        raise RuntimeError("the witness is zero")
    return abs(derived - witness) / abs(witness)


def _second_gap(rows, target: float, best: float) -> float:
    gaps = sorted(_gap(value, target) for _name, value in rows)
    if len(gaps) < 2 or abs(gaps[0] - best) > 1.0e-12:
        raise RuntimeError("the nearest dressing is missing from its own list")
    return gaps[1]


def _walk(target: float, anchor: NamedGap, expected: tuple[str, ...], alpha: float, helium: float, base) -> tuple[ServiceStep, ...]:
    name = anchor.name
    value = anchor.value
    steps: list[ServiceStep] = []
    for step_name in expected:
        rows = _dressings(name, value, alpha, base)
        hit = _nearest(rows, target)
        if hit.name != step_name:
            raise RuntimeError("the service dressing changed")
        if hit.relative_gap >= _gap(value, target):
            raise RuntimeError("the service dressing did not tighten")
        size = _size(rows)
        if size.expressions != _COUNT or size.inside != _COUNT:
            raise RuntimeError("the service dressing count changed")
        if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
            raise RuntimeError("the service spacing is not the base-piece spacing")
        inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
        if inside_helium > 1:
            raise RuntimeError("several service dressings meet the helium standard")
        steps.append(
            ServiceStep(
                dressing=hit,
                inside_review=sum(_gap(number, target) <= _REVIEW_GAP for _row, number in rows),
                inside_tenth=sum(_gap(number, target) <= _TENTH for _row, number in rows),
                inside_helium=inside_helium,
                second_gap=_second_gap(rows, target, hit.relative_gap),
                expressions=size.expressions,
                median_spacing=size.median_spacing,
            )
        )
        name = hit.name
        value = hit.value
    closing = steps[-1]
    if closing.inside_helium != 1 or closing.dressing.relative_gap > helium:
        raise RuntimeError("the closing service figure misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError("a second service dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier service figure already meets helium")
    return tuple(steps)


def _load(value: float, target: float, thermal_w: float, filed_w: float, small_w: float, filed_small_w: float):
    watts = value * thermal_w
    fit = target * filed_w
    small = value * small_w
    fit_small = target * filed_small_w
    return watts, fit, small, fit_small, _gap(watts, fit)


def build_service() -> ServiceReport:
    heat = build_decay_heat()
    removal = build_removal()
    module, digest = load_engine()
    if digest[:16] != heat.pin_prefix:
        raise RuntimeError("the service pin does not match the plant pin")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    five = _five_pieces(module)
    base = _base_pieces(module)
    operation_s = FULL_POWER_HOURS * 3600.0
    service_target = decay_fraction(_SERVICE_S, operation_s)
    cooldown_target = decay_fraction(_COOLDOWN_S, operation_s)
    fit_hours = _cross_hours(removal.duty_fraction, operation_s)
    cross_target = fit_hours * 3600.0 / operation_s
    service_anchor = _nearest(five, service_target)
    cooldown_anchor = _nearest(five, cooldown_target)
    cross_anchor = _nearest(five, cross_target)
    if service_anchor.name != _SERVICE_ANCHOR:
        raise RuntimeError("the four-hour anchor changed")
    if cooldown_anchor.name != _COOLDOWN_ANCHOR:
        raise RuntimeError("the twenty-hour anchor changed")
    if cross_anchor.name != _CROSS_ANCHOR:
        raise RuntimeError("the crossing anchor changed")
    service_steps = _walk(service_target, service_anchor, _SERVICE_STEPS, alpha, helium, base)
    cooldown_steps = _walk(cooldown_target, cooldown_anchor, _COOLDOWN_STEPS, alpha, helium, base)
    cross_steps = _walk(cross_target, cross_anchor, _CROSS_STEPS, alpha, helium, base)
    service_value = service_steps[-1].dressing.value
    cooldown_value = cooldown_steps[-1].dressing.value
    cross_value = cross_steps[-1].dressing.value
    if not service_value > cooldown_value > removal.derived_fraction:
        raise RuntimeError("an FSOT service fraction fell through the design-point duty")
    service_w, service_fit, service_small, service_fit_small, service_gap = _load(
        service_value,
        service_target,
        heat.thermal_w,
        heat.filed_thermal_w,
        heat.small_thermal_w,
        heat.filed_small_thermal_w,
    )
    cooldown_w, cooldown_fit, cooldown_small, cooldown_fit_small, cooldown_gap = _load(
        cooldown_value,
        cooldown_target,
        heat.thermal_w,
        heat.filed_thermal_w,
        heat.small_thermal_w,
        heat.filed_small_thermal_w,
    )
    cross_hours = cross_value * FULL_POWER_HOURS
    if _gap(cross_hours, fit_hours) > helium:
        raise RuntimeError("the crossing time misses the helium quotient match")
    return ServiceReport(
        service_anchor=service_anchor,
        service_steps=service_steps,
        service_target=service_target,
        service_w=service_w,
        service_fit_w=service_fit,
        service_small_w=service_small,
        service_fit_small_w=service_fit_small,
        service_watt_gap=service_gap,
        cooldown_anchor=cooldown_anchor,
        cooldown_steps=cooldown_steps,
        cooldown_target=cooldown_target,
        cooldown_w=cooldown_w,
        cooldown_fit_w=cooldown_fit,
        cooldown_small_w=cooldown_small,
        cooldown_fit_small_w=cooldown_fit_small,
        cooldown_watt_gap=cooldown_gap,
        cross_anchor=cross_anchor,
        cross_steps=cross_steps,
        cross_target=cross_target,
        cross_hours=cross_hours,
        fit_hours=fit_hours,
        duty_fraction=removal.duty_fraction,
        derived_duty=removal.derived_fraction,
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=heat.nuclear_scalar,
        thermodynamics_scalar=heat.thermodynamics_scalar,
        gamma_fsot=heat.gamma_fsot,
        valve_beta=heat.valve_beta,
        pin_prefix=heat.pin_prefix,
        digest_prefix=heat.digest_prefix,
    )


def _step_lines(step: ServiceStep, closing: bool) -> list[str]:
    lines = [
        f"The dressing is `{step.dressing.name}` at `{_num(step.dressing.value)}`.",
        f"The relative gap is `{_num(step.dressing.relative_gap)}`.",
    ]
    if closing:
        lines.append("The relative gap meets the helium quotient match.")
        lines.append("One dressing of that anchor meets the helium quotient match.")
        lines.append(f"The next dressing sits at relative gap `{_num(step.second_gap)}`.")
    else:
        lines.append("The dressing stays wider than the helium quotient match.")
    lines.append(
        f"The anchor has `{step.expressions}` dressings, median spacing `{_num(step.median_spacing)}`."
    )
    return lines


def _watt_lines(label: str, watts: float, fit: float, small: float, fit_small: float, gap: float, helium: float) -> list[str]:
    lines = [
        f"The {label} load is `{_num(watts)}` W.",
        f"The 1952 load is `{_num(fit)}` W.",
        f"The small-class load is `{_num(small)}` W.",
        f"The 1952 small-class load is `{_num(fit_small)}` W.",
        f"The watt gap is `{_num(gap)}`.",
    ]
    if gap <= helium:
        lines.append("The watt gap meets the helium quotient match.")
    else:
        lines.append("The watt gap stays wider than the helium quotient match.")
    return lines


def render() -> str:
    report = build_service()
    lines = [
        "# Fission service",
        "",
        "The residual-heat section names two times the removal score left on the 1952 fit.",
        "Four hours is the service entry.",
        "Twenty hours is the end of the published cooldown.",
        "The crossing is the share of the filed year at which the one-year fraction falls through the design-point duty.",
        "The 1952 fit remains the calorimetry witness.",
        "The chain stays at k_eff of 1.",
        "",
        "## Four hours",
        "",
        f"The 1952 fraction is `{_num(report.service_target)}`.",
        f"The five-atom nearest is `{report.service_anchor.name}` at `{_num(report.service_anchor.value)}`.",
        f"The anchor gap is `{_num(report.service_anchor.relative_gap)}`.",
    ]
    for step in report.service_steps[:-1]:
        lines.extend(_step_lines(step, False))
    lines.extend(_step_lines(report.service_steps[-1], True))
    lines.extend(
        _watt_lines(
            "four-hour",
            report.service_w,
            report.service_fit_w,
            report.service_small_w,
            report.service_fit_small_w,
            report.service_watt_gap,
            report.helium_mismatch,
        )
    )
    lines.extend(
        [
            "",
            "## Twenty hours",
            "",
            f"The 1952 fraction is `{_num(report.cooldown_target)}`.",
            f"The five-atom nearest is `{report.cooldown_anchor.name}` at `{_num(report.cooldown_anchor.value)}`.",
            f"The anchor gap is `{_num(report.cooldown_anchor.relative_gap)}`.",
        ]
    )
    for step in report.cooldown_steps[:-1]:
        lines.extend(_step_lines(step, False))
    lines.extend(_step_lines(report.cooldown_steps[-1], True))
    lines.extend(
        _watt_lines(
            "twenty-hour",
            report.cooldown_w,
            report.cooldown_fit_w,
            report.cooldown_small_w,
            report.cooldown_fit_small_w,
            report.cooldown_watt_gap,
            report.helium_mismatch,
        )
    )
    lines.extend(
        [
            "The last factor is the same base piece that closes the design-point duty, with the same minus sign.",
            "The twenty-hour fraction stays above the design-point duty fraction.",
            "",
            "## Crossing",
            "",
            f"The 1952 crossing is `{_num(report.fit_hours)}` hours.",
            f"That time is `{_num(report.cross_target)}` of the filed year.",
            f"The five-atom nearest is `{report.cross_anchor.name}` at `{_num(report.cross_anchor.value)}`.",
            f"The anchor gap is `{_num(report.cross_anchor.relative_gap)}`.",
        ]
    )
    for step in report.cross_steps[:-1]:
        lines.extend(_step_lines(step, False))
    lines.extend(_step_lines(report.cross_steps[-1], True))
    lines.extend(
        [
            f"The derived crossing is `{_num(report.cross_hours)}` hours.",
            f"The relative gap is `{_num(_gap(report.cross_hours, report.fit_hours))}`.",
            "The relative gap meets the helium quotient match.",
            "The four-hour fraction stays above the design-point duty fraction.",
            "The UFSAR prints no watt rate at the 350°F entry.",
            "This score assigns that hotter rate no figure.",
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
            f"alpha in the dressing is the fine-structure leaf `{_num(report.alpha)}`.",
            f"The helium quotient match is `{_num(report.helium_mismatch)}`.",
            "Every closing gap meets the helium quotient match.",
            "The design-point duty stays the removal score.",
            "This run does not install a seed as a new enrichment or a new plant row.",
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
