"""Name the steam figures the ceiling score left on page seeds.

The demonstrated efficiency stays the heat-ratio score. The temperature
Carnot stays the page temperatures. The catalog ceiling stays above that
Carnot. The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay_dress import _base_pieces
from solcascade.production.fission_decay_heat import _dressings, _nearest_row, _positive_nearest, _positive_size
from solcascade.production.fission_decay_seed import NamedGap, _five_pieces
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.fission_process import FABRICATION_CAP_PERCENT, build_process
from solcascade.production.fission_safety import BUILT_EFFICIENCY
from solcascade.production.plants import build_plants
from solcascade.production.refine import project_root

_COUNT = 116
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4
_BUILT_ANCHOR = "pi^-3*phi^5"
_BUILT_STEPS = (
    "(pi^-3*phi^5)*(1+alpha*G)",
    "((pi^-3*phi^5)*(1+alpha*G))*(1-alpha*e^-4)",
    "(((pi^-3*phi^5)*(1+alpha*G))*(1-alpha*e^-4))*(1-alpha*e^-5)",
)
_CARNOT_ANCHOR = "phi^-2+gamma^5"
_CARNOT_STEPS = (
    "(phi^-2+gamma^5)*(1+alpha*gamma^4)",
    "((phi^-2+gamma^5)*(1+alpha*gamma^4))*(1+alpha*e^-5)",
    "(((phi^-2+gamma^5)*(1+alpha*gamma^4))*(1+alpha*e^-5))*(1+alpha*pi^-6)",
)
_STEP_ANCHOR = "e^-3+G^-3"
_STEP_STEPS = (
    "(e^-3+G^-3)*(1+alpha*gamma^3)",
    "((e^-3+G^-3)*(1+alpha*gamma^3))*(1-alpha*pi^-4)",
    "(((e^-3+G^-3)*(1+alpha*gamma^3))*(1-alpha*pi^-4))*(1-alpha*e^-6)",
)


@dataclass(frozen=True)
class SteamStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class SteamFigure:
    anchor: NamedGap
    steps: tuple[SteamStep, ...]
    target: float


@dataclass(frozen=True)
class SteamReport:
    built: SteamFigure
    carnot: SteamFigure
    carnot_step: SteamFigure
    demonstrated_efficiency: float
    catalog_ceiling: float
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
    return project_root() / "runs" / "production" / "fission_steam.md"


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


def _walk(target: float, anchor_name: str, expected: tuple[str, ...], alpha: float, helium: float, five, base) -> SteamFigure:
    name, value, anchor_gap = _positive_nearest(five, target)
    if name != anchor_name:
        raise RuntimeError("the steam anchor changed")
    anchor = NamedGap(name, value, anchor_gap, anchor_gap <= _REVIEW_GAP)
    steps: list[SteamStep] = []
    for step_name in expected:
        rows = _dressings(name, value, alpha, base)
        hit_name, hit_value, hit_gap = _nearest_row(rows, target)
        if hit_name != step_name:
            raise RuntimeError("the steam dressing changed")
        if hit_gap >= _gap(value, target):
            raise RuntimeError("the steam dressing did not tighten")
        size = _positive_size(rows)
        if size.expressions != _COUNT or size.inside != _COUNT:
            raise RuntimeError("the steam dressing count changed")
        if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
            raise RuntimeError("the steam spacing is not the base-piece spacing")
        inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
        if inside_helium > 1:
            raise RuntimeError("several steam dressings meet the helium standard")
        steps.append(
            SteamStep(
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
        raise RuntimeError("the closing steam figure misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError("a second steam dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier steam figure already meets helium")
    return SteamFigure(anchor, tuple(steps), target)


def build_steam() -> SteamReport:
    process = build_process()
    plants = build_plants()
    module, digest = load_engine()
    if digest[:16] != process.pin_prefix:
        raise RuntimeError("the steam pin does not match the plant pin")
    if FABRICATION_CAP_PERCENT != 5.0:
        raise RuntimeError("the fabrication cap left 5 percent")
    ceiling = plants.thermal.ceiling
    if not process.efficiency < process.carnot < ceiling:
        raise RuntimeError("the temperature Carnot is not under the catalog ceiling")
    if _gap(process.carnot_step, process.carnot / process.efficiency) > 1.0e-12:
        raise RuntimeError("the Carnot step left the temperature quotient")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    five = _five_pieces(module)
    base = _base_pieces(module)
    built = _walk(BUILT_EFFICIENCY, _BUILT_ANCHOR, _BUILT_STEPS, alpha, helium, five, base)
    carnot = _walk(process.carnot, _CARNOT_ANCHOR, _CARNOT_STEPS, alpha, helium, five, base)
    step = _walk(process.carnot_step, _STEP_ANCHOR, _STEP_STEPS, alpha, helium, five, base)
    built_value = built.steps[-1].dressing.value
    carnot_value = carnot.steps[-1].dressing.value
    if not process.efficiency < built_value < process.carnot:
        raise RuntimeError("the named built efficiency left the open interval under the Carnot")
    if not process.efficiency < carnot_value < ceiling:
        raise RuntimeError("the named Carnot left the open interval under the catalog ceiling")
    return SteamReport(
        built=built,
        carnot=carnot,
        carnot_step=step,
        demonstrated_efficiency=process.efficiency,
        catalog_ceiling=ceiling,
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


def _step_lines(step: SteamStep, closing: bool) -> list[str]:
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


def _figure_lines(title: str, intro: list[str], figure: SteamFigure) -> list[str]:
    lines = ["", f"## {title}", "", *intro]
    lines.append(f"The positive five-atom nearest is `{figure.anchor.name}` at `{_num(figure.anchor.value)}`.")
    lines.append(f"The anchor gap is `{_num(figure.anchor.relative_gap)}`.")
    for step in figure.steps[:-1]:
        lines.extend(_step_lines(step, False))
    lines.extend(_step_lines(figure.steps[-1], True))
    return lines


def render() -> str:
    report = build_steam()
    lines = [
        "# Fission steam",
        "",
        "Plants now being built use the page efficiency.",
        "The temperature Carnot stays the page temperatures.",
        "The catalog ceiling stays above that Carnot and is not run.",
        "The demonstrated efficiency stays the heat-ratio score.",
        "Ratios that leave the unit interval use the positive dressing from the heat-ratio score.",
        "The chain stays at k_eff of 1.",
    ]
    lines.extend(
        _figure_lines(
            "Plants being built",
            [f"The page efficiency is `{_num(report.built.target)}`."],
            report.built,
        )
    )
    lines.extend(
        _figure_lines(
            "Temperature Carnot",
            [f"The page temperatures give `{_num(report.carnot.target)}`."],
            report.carnot,
        )
    )
    lines.extend(
        _figure_lines(
            "Carnot step",
            [
                f"The temperature Carnot divided by the demonstrated efficiency is `{_num(report.carnot_step.target)}`.",
                "The 60 GWd/t stack stays unrun.",
            ],
            report.carnot_step,
        )
    )
    lines.extend(
        [
            "The fabrication cap stays 5 percent.",
            f"The catalog ceiling stays `{_num(report.catalog_ceiling)}`.",
            "This run does not install a seed as a new efficiency or a new plant row.",
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
