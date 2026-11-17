"""Name the two steam steps the plants-being-built figure still left open.

The built step is that page efficiency over the demonstrated efficiency.
The upper step is the temperature Carnot over the page efficiency.
The efficiency differences stay products of those named figures.
The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_ceiling import build_ceiling
from solcascade.production.fission_decay_dress import _base_pieces
from solcascade.production.fission_decay_heat import _dressings as positive_dressings
from solcascade.production.fission_decay_heat import _nearest_row, _positive_nearest, _positive_size
from solcascade.production.fission_decay_seed import NamedGap, _five_pieces
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.fission_safety import BUILT_EFFICIENCY
from solcascade.production.plants import FAST_ID, THERMAL_ID, build_plants
from solcascade.production.refine import catalog_record, project_root

_COUNT = 116
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4
_BUILT_ANCHOR = "G*G^-2"
_BUILT_STEPS = (
    "(G*G^-2)*(1-alpha*pi^-2)",
    "((G*G^-2)*(1-alpha*pi^-2))*(1-alpha*pi^-5)",
)
_CARNOT_ANCHOR = "e^-3+G^-2"
_CARNOT_STEPS = (
    "(e^-3+G^-2)*(1-alpha*gamma^3)",
    "((e^-3+G^-2)*(1-alpha*gamma^3))*(1+alpha*e^-4)",
    "(((e^-3+G^-2)*(1-alpha*gamma^3))*(1+alpha*e^-4))*(1-alpha*pi^-5)",
    "((((e^-3+G^-2)*(1-alpha*gamma^3))*(1+alpha*e^-4))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
)


@dataclass(frozen=True)
class BuiltStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class BuiltFigure:
    anchor: NamedGap
    steps: tuple[BuiltStep, ...]
    target: float


@dataclass(frozen=True)
class BuiltReport:
    built_step: BuiltFigure
    carnot_step: BuiltFigure
    demonstrated: float
    built: float
    carnot: float
    catalog: float
    fast_headroom: float
    built_difference: float
    carnot_difference: float
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_built.md"


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


def _record(rows, target: float, hit: NamedGap, helium: float) -> BuiltStep:
    size = _positive_size(rows)
    if size.expressions != _COUNT or size.inside != _COUNT:
        raise RuntimeError("the built dressing count changed")
    if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
        raise RuntimeError("the built spacing is not the base-piece spacing")
    inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
    if inside_helium > 1:
        raise RuntimeError("several built dressings meet the helium standard")
    return BuiltStep(
        dressing=hit,
        inside_review=sum(_gap(number, target) <= _REVIEW_GAP for _row, number in rows),
        inside_tenth=sum(_gap(number, target) <= _TENTH for _row, number in rows),
        inside_helium=inside_helium,
        second_gap=_second_gap(rows, target, hit.relative_gap),
        expressions=size.expressions,
        median_spacing=size.median_spacing,
    )


def _walk(target, anchor_name, expected, alpha, helium, five, base) -> BuiltFigure:
    if not target > 1.0:
        raise RuntimeError("the built step left the ratios above one")
    name, value, anchor_gap = _positive_nearest(five, target)
    if name != anchor_name:
        raise RuntimeError(f"the built anchor changed to {name}")
    if anchor_gap <= helium:
        raise RuntimeError("the built anchor already meets helium")
    anchor = NamedGap(name, value, anchor_gap, anchor_gap <= _REVIEW_GAP)
    steps: list[BuiltStep] = []
    for step_name in expected:
        rows = positive_dressings(name, value, alpha, base)
        hit_name, hit_value, hit_gap = _nearest_row(rows, target)
        if hit_name != step_name:
            raise RuntimeError(f"the built dressing changed to {hit_name}")
        if hit_gap >= _gap(value, target):
            raise RuntimeError("the built dressing did not tighten")
        recorded = _record(
            rows, target, NamedGap(hit_name, hit_value, hit_gap, hit_gap <= _REVIEW_GAP), helium
        )
        steps.append(recorded)
        name, value = hit_name, hit_value
    closing = steps[-1]
    if closing.inside_helium != 1 or closing.dressing.relative_gap > helium:
        raise RuntimeError("the closing built figure misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError("a second built dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier built figure already meets helium")
    return BuiltFigure(anchor, tuple(steps), target)


def build_built() -> BuiltReport:
    plants = build_plants()
    ceiling = build_ceiling()
    module, digest = load_engine()
    if digest[:16] != plants.pin_prefix or plants.digest_prefix != plants.pin_prefix:
        raise RuntimeError("the built pin does not match the plant pin")
    if ceiling.pin_prefix != plants.pin_prefix:
        raise RuntimeError("the ceiling pin does not match the plant pin")
    thermal = plants.thermal
    fast = plants.fast
    record = catalog_record(THERMAL_ID)
    if "k_eff is held at 1" not in str(record["governing_equation"]):
        raise RuntimeError("the governing equation does not hold k_eff at 1")
    if "Sodium" not in str(catalog_record(FAST_ID).get("engineering_barrier") or ""):
        raise RuntimeError("the fast barrier does not name sodium")
    if BUILT_EFFICIENCY != 0.36:
        raise RuntimeError("the plants-being-built figure left 0.36")
    if _gap(thermal.efficiency, 0.33) > 1.0e-12:
        raise RuntimeError("the demonstrated efficiency left 0.33")
    if ceiling.hot_k != 560.0 or ceiling.cold_k != 310.0:
        raise RuntimeError("the ceiling temperatures left 560 K and 310 K")
    if _gap(ceiling.carnot, 1.0 - ceiling.cold_k / ceiling.hot_k) > 0.0:
        raise RuntimeError("the temperature Carnot left one minus the rejected fraction")
    if thermal.ceiling != 0.45:
        raise RuntimeError("the catalog ceiling left 0.45")
    demonstrated = thermal.efficiency
    built = BUILT_EFFICIENCY
    carnot = ceiling.carnot
    if not demonstrated < built < carnot < thermal.ceiling:
        raise RuntimeError("the plants-being-built figure left the open interval under the Carnot")
    built_step = built / demonstrated
    carnot_step = carnot / built
    built_difference = built - demonstrated
    carnot_difference = carnot - built
    if _gap(built_difference, demonstrated * (built_step - 1.0)) > 1.0e-15:
        raise RuntimeError("the built difference left the demonstrated efficiency times the step")
    if _gap(carnot_difference, built * (carnot_step - 1.0)) > 1.0e-15:
        raise RuntimeError("the Carnot difference left the built efficiency times the step")
    if _gap(thermal.ceiling / built, fast.headroom) > 1.0e-15:
        raise RuntimeError("the catalog ceiling over the built figure left the fast headroom")
    if _gap(fast.headroom, fast.ceiling / fast.efficiency) > 1.0e-15:
        raise RuntimeError("the fast headroom left the ceiling quotient")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    five = _five_pieces(module)
    base = _base_pieces(module)
    carnot_figure = _walk(carnot_step, _CARNOT_ANCHOR, _CARNOT_STEPS, alpha, helium, five, base)
    if len(carnot_figure.steps) != 4 or carnot_figure.steps[-2].inside_helium != 0:
        raise RuntimeError("the Carnot step over the built figure met helium before the fourth factor")
    return BuiltReport(
        built_step=_walk(built_step, _BUILT_ANCHOR, _BUILT_STEPS, alpha, helium, five, base),
        carnot_step=carnot_figure,
        demonstrated=demonstrated,
        built=built,
        carnot=carnot,
        catalog=thermal.ceiling,
        fast_headroom=fast.headroom,
        built_difference=built_difference,
        carnot_difference=carnot_difference,
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=plants.nuclear_scalar,
        thermodynamics_scalar=plants.thermodynamics_scalar,
        gamma_fsot=plants.gamma_fsot,
        valve_beta=plants.valve_beta,
        pin_prefix=plants.pin_prefix,
        digest_prefix=plants.digest_prefix,
    )


def _step_lines(step: BuiltStep, closing: bool) -> list[str]:
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


def _figure_lines(figure: BuiltFigure, lead: str, role: str, fourth: bool) -> list[str]:
    lines = [
        lead,
        f"The five-atom nearest is `{figure.anchor.name}` at `{_num(figure.anchor.value)}`.",
        f"The anchor gap is `{_num(figure.anchor.relative_gap)}`.",
    ]
    for step in figure.steps[:-1]:
        lines.extend(_step_lines(step, False))
    if fourth:
        lines.append("Three factors stay wider than the helium quotient match.")
        lines.append("The fourth factor is one more base-piece factor of the same form.")
    lines.extend(_step_lines(figure.steps[-1], True))
    if _num(figure.steps[-1].dressing.value) == _num(figure.target):
        lines.append(f"The closing value agrees with `{_num(figure.target)}` at six figures.")
    lines.append("The closing value names that quotient.")
    lines.append(role)
    return lines


def render() -> str:
    report = build_built()
    lines = [
        "# Fission built",
        "",
        "The cooling page states plants now being built at its worked efficiency.",
        "The demonstrated plant stays at its filed efficiency.",
        "The temperature Carnot stays the bound from the page temperatures.",
        "FSOT names the two steps those three figures leave open.",
        "The chain stays at k_eff of 1.",
        "",
        "## Built step",
        "",
        "Ratios above one use the positive dressing from the heat-ratio score.",
        f"The demonstrated efficiency is `{_num(report.demonstrated)}`.",
        f"The plants-being-built efficiency is `{_num(report.built)}`.",
        "The demonstrated efficiency stays the heat-ratio score.",
        "The plants-being-built figure stays the page efficiency.",
    ]
    lines.extend(
        _figure_lines(
            report.built_step,
            (
                "The plants-being-built efficiency over the demonstrated efficiency is "
                f"`{_num(report.built_step.target)}`."
            ),
            "The operating step stays the plants-being-built efficiency over the demonstrated efficiency.",
            False,
        )
    )
    lines.extend(
        [
            f"The efficiency difference is `{_num(report.built_difference)}`.",
            "That difference is the demonstrated efficiency times one less than this step.",
            "That difference keeps the heat-ratio score and this step.",
            "",
            "## Carnot step",
            "",
            f"The temperature Carnot is `{_num(report.carnot)}`.",
            "The temperature Carnot stays the page temperatures.",
        ]
    )
    lines.extend(
        _figure_lines(
            report.carnot_step,
            (
                "The temperature Carnot over the plants-being-built efficiency is "
                f"`{_num(report.carnot_step.target)}`."
            ),
            "The operating step stays the temperature Carnot over the plants-being-built efficiency.",
            True,
        )
    )
    lines.extend(
        [
            f"The remaining difference up to the temperature Carnot is `{_num(report.carnot_difference)}`.",
            "That difference is the plants-being-built efficiency times one less than this step.",
            "That difference keeps the named plants-being-built efficiency and this step.",
            f"The catalog ceiling is `{_num(report.catalog)}`.",
            f"The fast headroom is `{_num(report.fast_headroom)}`.",
            "The catalog ceiling over the plants-being-built figure is the fast headroom.",
            "The fast headroom stays the bounds score.",
            "One expression names that shared number.",
            "The catalog ceiling stays unrun.",
            "The demonstrated plant stays at its filed efficiency.",
            "This score adds no joules.",
            "The fast row's filed barrier is the sodium system.",
            "The 60 GWd/t stack stays unrun.",
            "The fabrication cap stays 5 percent.",
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
            "This run does not install a seed as a new efficiency or a new plant row.",
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
