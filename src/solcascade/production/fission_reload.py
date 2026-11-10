"""Name the material-balance year ratios the process score left on the page.

The positive dressing is the heat-ratio grammar. It keeps every positive
value. The page enrichment stays 4.5 percent. The checked-in pulse is not
rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay_dress import _base_pieces
from solcascade.production.fission_decay_heat import _dressings, _nearest_row, _positive_nearest, _positive_size
from solcascade.production.fission_decay_seed import NamedGap, _five_pieces
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.fission_process import (
    ANNUAL_SWU,
    FABRICATION_CAP_PERCENT,
    RELOAD_TONNES,
    build_process,
)
from solcascade.production.refine import project_root

_COUNT = 116
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4
_COMPANION_ANCHOR = "sqrt(gamma)*G^-3"
_COMPANION_STEPS = ("(sqrt(gamma)*G^-3)*(1+alpha*e^-6)",)
_SWU_ANCHOR = "pi^-2+e^2"
_SWU_STEPS = (
    "(pi^-2+e^2)*(1-alpha*pi^-4)",
    "((pi^-2+e^2)*(1-alpha*pi^-4))*(1-alpha*e^-6)",
)
_HEU_ANCHOR = "phi^5+gamma^-5"
_HEU_STEPS = (
    "(phi^5+gamma^-5)*(1+alpha*pi^-3)",
    "((phi^5+gamma^-5)*(1+alpha*pi^-3))*(1+alpha*pi^-6)",
)
_RETURN_ANCHOR = "pi^6-gamma^-4"
_RETURN_STEPS = (
    "(pi^6-gamma^-4)*(1-alpha*pi^-2)",
    "((pi^6-gamma^-4)*(1-alpha*pi^-2))*(1+alpha*pi^-5)",
    "(((pi^6-gamma^-4)*(1-alpha*pi^-2))*(1+alpha*pi^-5))*(1+alpha*pi^-6)",
)


@dataclass(frozen=True)
class ReloadStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class ReloadFigure:
    anchor: NamedGap
    steps: tuple[ReloadStep, ...]
    target: float


@dataclass(frozen=True)
class ReloadReport:
    companion: ReloadFigure
    separative: ReloadFigure
    highly_enriched: ReloadFigure
    energy_return: ReloadFigure
    annual_swu: float
    reload_tonnes: float
    fabrication_cap_percent: float
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_reload.md"


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


def _walk(target: float, anchor_name: str, expected: tuple[str, ...], alpha: float, helium: float, five, base) -> ReloadFigure:
    name, value, anchor_gap = _positive_nearest(five, target)
    if name != anchor_name:
        raise RuntimeError("the reload anchor changed")
    anchor = NamedGap(name, value, anchor_gap, anchor_gap <= _REVIEW_GAP)
    steps: list[ReloadStep] = []
    for step_name in expected:
        rows = _dressings(name, value, alpha, base)
        hit_name, hit_value, hit_gap = _nearest_row(rows, target)
        if hit_name != step_name:
            raise RuntimeError("the reload dressing changed")
        if hit_gap >= _gap(value, target):
            raise RuntimeError("the reload dressing did not tighten")
        size = _positive_size(rows)
        if size.expressions != _COUNT or size.inside != _COUNT:
            raise RuntimeError("the reload dressing count changed")
        if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
            raise RuntimeError("the reload spacing is not the base-piece spacing")
        inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
        if inside_helium > 1:
            raise RuntimeError("several reload dressings meet the helium standard")
        steps.append(
            ReloadStep(
                dressing=NamedGap(hit_name, hit_value, hit_gap, hit_gap <= _REVIEW_GAP),
                inside_review=sum(_gap(number, target) <= _REVIEW_GAP for _row, number in rows),
                inside_tenth=sum(_gap(number, target) <= _TENTH for _row, number in rows),
                inside_helium=inside_helium,
                second_gap=_second_gap(rows, target, hit_gap),
                expressions=size.expressions,
                median_spacing=size.median_spacing,
            )
        )
        name = hit_name
        value = hit_value
    closing = steps[-1]
    if closing.inside_helium != 1 or closing.dressing.relative_gap > helium:
        raise RuntimeError("the closing reload figure misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError("a second reload dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier reload figure already meets helium")
    return ReloadFigure(anchor, tuple(steps), target)


def build_reload() -> ReloadReport:
    process = build_process()
    module, digest = load_engine()
    if digest[:16] != process.pin_prefix:
        raise RuntimeError("the reload pin does not match the plant pin")
    if FABRICATION_CAP_PERCENT != 5.0:
        raise RuntimeError("the fabrication cap left 5 percent")
    reload_kg = RELOAD_TONNES * 1000.0
    if _gap(process.swu_per_kg_reload, ANNUAL_SWU / reload_kg) > 1.0e-12:
        raise RuntimeError("the separative work per kilogram left the page quotient")
    if _gap(process.heu_swu_over_reload, process.swu_per_kg_heu / process.swu_per_kg_reload) > 1.0e-12:
        raise RuntimeError("the highly enriched multiple left the page quotient")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    five = _five_pieces(module)
    base = _base_pieces(module)
    companion = process.fuel_side_j / process.nameplate_j
    return ReloadReport(
        companion=_walk(companion, _COMPANION_ANCHOR, _COMPANION_STEPS, alpha, helium, five, base),
        separative=_walk(process.swu_per_kg_reload, _SWU_ANCHOR, _SWU_STEPS, alpha, helium, five, base),
        highly_enriched=_walk(process.heu_swu_over_reload, _HEU_ANCHOR, _HEU_STEPS, alpha, helium, five, base),
        energy_return=_walk(
            process.reload_output_over_input, _RETURN_ANCHOR, _RETURN_STEPS, alpha, helium, five, base
        ),
        annual_swu=ANNUAL_SWU,
        reload_tonnes=RELOAD_TONNES,
        fabrication_cap_percent=FABRICATION_CAP_PERCENT,
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=process.nuclear_scalar,
        thermodynamics_scalar=process.thermodynamics_scalar,
        gamma_fsot=process.gamma_fsot,
        valve_beta=process.valve_beta,
        pin_prefix=process.pin_prefix,
        digest_prefix=digest[:16],
    )


def _step_lines(step: ReloadStep, closing: bool) -> list[str]:
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


def _figure_lines(title: str, intro: list[str], figure: ReloadFigure) -> list[str]:
    lines = ["", f"## {title}", "", *intro]
    lines.append(f"The positive five-atom nearest is `{figure.anchor.name}` at `{_num(figure.anchor.value)}`.")
    lines.append(f"The anchor gap is `{_num(figure.anchor.relative_gap)}`.")
    for step in figure.steps[:-1]:
        lines.extend(_step_lines(step, False))
    lines.extend(_step_lines(figure.steps[-1], True))
    return lines


def render() -> str:
    report = build_reload()
    lines = [
        "# Fission reload",
        "",
        "The material-balance year spends the page separative work on the page reload mass.",
        "The fuel-side year and the nameplate year stay that page pair.",
        "The highly enriched estimate stays the catalog input.",
        "Ratios that leave the unit interval use the positive dressing from the heat-ratio score.",
        "The chain stays at k_eff of 1.",
    ]
    lines.extend(
        _figure_lines(
            "Fuel-year companion",
            [f"The filed companion is `{_num(report.companion.target)}`."],
            report.companion,
        )
    )
    lines.extend(
        _figure_lines(
            "Separative work",
            [
                f"The overview spends `{_num(report.annual_swu)}` SWU on `{_num(report.reload_tonnes)}` tonnes.",
                f"That quotient is `{_num(report.separative.target)}` SWU per kilogram.",
            ],
            report.separative,
        )
    )
    lines.extend(
        _figure_lines(
            "Highly enriched multiple",
            [
                f"The catalog estimate is `{_num(report.highly_enriched.target)}` times the reload separative work per kilogram.",
            ],
            report.highly_enriched,
        )
    )
    lines.extend(
        _figure_lines(
            "Reload return",
            [
                f"Discharge electricity over the reload enrichment electricity is `{_num(report.energy_return.target)}`.",
                "This ratio compares those two filed energies.",
                "It is not a second plant output.",
            ],
            report.energy_return,
        )
    )
    lines.extend(
        [
            "The page enrichment stays 4.5 percent.",
            "The fabrication cap stays 5 percent.",
            "The 60 GWd/t stack stays unrun.",
            "This run does not install a seed as a new enrichment or a new plant row.",
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
