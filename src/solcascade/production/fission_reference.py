"""Name the two filed comparisons the plant rows still left open.

The Pu-239 kilogram over the U-235 kilogram is the 200 MeV ledger.
The Surry core over the scored operating heat is the published rating
over the nameplate heat. The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay import build_decay
from solcascade.production.fission_decay_dress import _base_pieces
from solcascade.production.fission_decay_dress import _dressings as unit_dressings
from solcascade.production.fission_decay_heat import _dressings as positive_dressings
from solcascade.production.fission_decay_heat import _nearest_row, _positive_nearest, _positive_size
from solcascade.production.fission_decay_seed import NamedGap, _five_pieces
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.fission_process import NAMEPLATE_MWE
from solcascade.production.fission_removal import _CORE_W
from solcascade.production.plants import FAST_ID, THERMAL_ID, build_plants
from solcascade.production.refine import catalog_record, project_root

_COUNT = 116
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4
_KILOGRAM_ANCHOR = "phi^-4+G^2"
_KILOGRAM_STEPS = (
    "(phi^-4+G^2)*(1-alpha*phi^-3)",
    "((phi^-4+G^2)*(1-alpha*phi^-3))*(1+alpha*e^-5)",
    "(((phi^-4+G^2)*(1-alpha*phi^-3))*(1+alpha*e^-5))*(1-alpha*pi^-6)",
)
_CORE_ANCHOR = "sqrt(pi)-G"
_CORE_STEPS = (
    "(sqrt(pi)-G)*(1-alpha*phi^-2)",
    "((sqrt(pi)-G)*(1-alpha*phi^-2))*(1-alpha*gamma^5)",
    "(((sqrt(pi)-G)*(1-alpha*phi^-2))*(1-alpha*gamma^5))*(1+alpha*pi^-6)",
)


@dataclass(frozen=True)
class ReferenceStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class ReferenceFigure:
    anchor: NamedGap
    steps: tuple[ReferenceStep, ...]
    target: float


@dataclass(frozen=True)
class ReferenceReport:
    kilogram: ReferenceFigure
    core_ratio: ReferenceFigure
    efficiency_step: float
    delivered_ratio: float
    return_ratio: float
    core_w: float
    heat_w: float
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_reference.md"


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


def _record(rows, target: float, hit: NamedGap, helium: float) -> ReferenceStep:
    size = _positive_size(rows)
    if size.expressions != _COUNT or size.inside != _COUNT:
        raise RuntimeError("the reference dressing count changed")
    if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
        raise RuntimeError("the reference spacing is not the base-piece spacing")
    inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
    if inside_helium > 1:
        raise RuntimeError("several reference dressings meet the helium standard")
    return ReferenceStep(
        dressing=hit,
        inside_review=sum(_gap(number, target) <= _REVIEW_GAP for _row, number in rows),
        inside_tenth=sum(_gap(number, target) <= _TENTH for _row, number in rows),
        inside_helium=inside_helium,
        second_gap=_second_gap(rows, target, hit.relative_gap),
        expressions=size.expressions,
        median_spacing=size.median_spacing,
    )


def _walk(target, anchor_name, expected, alpha, helium, five, base) -> ReferenceFigure:
    name, value, anchor_gap = _positive_nearest(five, target)
    if name != anchor_name:
        raise RuntimeError(f"the reference anchor changed to {name}")
    if anchor_gap <= helium:
        raise RuntimeError("the reference anchor already meets helium")
    if len(unit_dressings(name, value, alpha, base)) >= _COUNT:
        raise RuntimeError("the unit grammar still has 116 dressings")
    anchor = NamedGap(name, value, anchor_gap, anchor_gap <= _REVIEW_GAP)
    steps: list[ReferenceStep] = []
    for step_name in expected:
        rows = positive_dressings(name, value, alpha, base)
        hit_name, hit_value, hit_gap = _nearest_row(rows, target)
        if hit_name != step_name:
            raise RuntimeError(f"the reference dressing changed to {hit_name}")
        if hit_gap >= _gap(value, target):
            raise RuntimeError("the reference dressing did not tighten")
        recorded = _record(
            rows, target, NamedGap(hit_name, hit_value, hit_gap, hit_gap <= _REVIEW_GAP), helium
        )
        steps.append(recorded)
        name, value = hit_name, hit_value
    closing = steps[-1]
    if closing.inside_helium != 1 or closing.dressing.relative_gap > helium:
        raise RuntimeError("the closing reference figure misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError("a second reference dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier reference figure already meets helium")
    return ReferenceFigure(anchor, tuple(steps), target)


def build_reference() -> ReferenceReport:
    plants = build_plants()
    decay = build_decay()
    module, digest = load_engine()
    if digest[:16] != plants.pin_prefix or plants.digest_prefix != plants.pin_prefix:
        raise RuntimeError("the reference pin does not match the plant pin")
    if decay.pin_prefix != plants.pin_prefix:
        raise RuntimeError("the decay pin does not match the plant pin")
    thermal = plants.thermal
    fast = plants.fast
    if thermal.mechanism_id != THERMAL_ID or fast.mechanism_id != FAST_ID:
        raise RuntimeError("the reference rows left the filed plants")
    if "U-235" not in thermal.reference_unit or "Pu-239" not in fast.reference_unit:
        raise RuntimeError("the reference units left U-235 and Pu-239")
    thermal_record = catalog_record(THERMAL_ID)
    if "k_eff is held at 1" not in str(thermal_record["governing_equation"]):
        raise RuntimeError("the governing equation does not hold k_eff at 1")
    fast_record = catalog_record(FAST_ID)
    if "Sodium" not in str(fast_record.get("engineering_barrier") or ""):
        raise RuntimeError("the fast barrier does not name sodium")
    if "estimate" not in str(fast_record.get("power_density_basis") or ""):
        raise RuntimeError("the fast power density left the estimate basis")
    if "about" not in str(thermal_record.get("power_density_basis") or ""):
        raise RuntimeError("the thermal power density left the stated order")
    if thermal.input_j != fast.input_j:
        raise RuntimeError("the two rows left the shared enrichment input")
    if _gap(plants.efficiency_step, fast.efficiency / thermal.efficiency) > 1.0e-15:
        raise RuntimeError("the efficiency step left the fast-over-thermal quotient")
    kilogram = fast.energy_density_j / thermal.energy_density_j
    delivered = fast.delivered_j / thermal.delivered_j
    returned = fast.output_over_input / thermal.output_over_input
    if not 0.0 < kilogram < 1.0:
        raise RuntimeError("the kilogram ratio left the unit interval")
    if _gap(delivered, kilogram * plants.efficiency_step) > 1.0e-15:
        raise RuntimeError("the delivery ratio left the kilogram ratio times the efficiency step")
    if _gap(returned, delivered) > 1.0e-15:
        raise RuntimeError("the return ratio left the delivery ratio")
    heat = decay.thermal_w
    if _gap(heat, NAMEPLATE_MWE * 1.0e6 / thermal.efficiency) > 1.0e-15:
        raise RuntimeError("the scored heat left the nameplate over the filed efficiency")
    core_ratio = _CORE_W / heat
    if not 0.0 < core_ratio < 1.0:
        raise RuntimeError("the Surry core ratio left the unit interval")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    five = _five_pieces(module)
    base = _base_pieces(module)
    return ReferenceReport(
        kilogram=_walk(kilogram, _KILOGRAM_ANCHOR, _KILOGRAM_STEPS, alpha, helium, five, base),
        core_ratio=_walk(core_ratio, _CORE_ANCHOR, _CORE_STEPS, alpha, helium, five, base),
        efficiency_step=plants.efficiency_step,
        delivered_ratio=delivered,
        return_ratio=returned,
        core_w=_CORE_W,
        heat_w=heat,
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=plants.nuclear_scalar,
        thermodynamics_scalar=plants.thermodynamics_scalar,
        gamma_fsot=plants.gamma_fsot,
        valve_beta=plants.valve_beta,
        pin_prefix=plants.pin_prefix,
        digest_prefix=plants.digest_prefix,
    )


def _step_lines(step: ReferenceStep, closing: bool) -> list[str]:
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


def _figure_lines(figure: ReferenceFigure, lead: str, role: str) -> list[str]:
    lines = [
        lead,
        f"The five-atom nearest is `{figure.anchor.name}` at `{_num(figure.anchor.value)}`.",
        f"The anchor gap is `{_num(figure.anchor.relative_gap)}`.",
    ]
    for step in figure.steps[:-1]:
        lines.extend(_step_lines(step, False))
    lines.extend(_step_lines(figure.steps[-1], True))
    if _num(figure.steps[-1].dressing.value) == _num(figure.target):
        lines.append(f"The closing value agrees with `{_num(figure.target)}` at six figures.")
    lines.append("The closing value names that quotient.")
    lines.append(role)
    return lines


def render() -> str:
    report = build_reference()
    lines = [
        "# Fission reference",
        "",
        "The plant rows leave two comparisons without their own names.",
        "One is the Pu-239 reference kilogram over the U-235 reference kilogram.",
        "The other is the Surry core rating over the scored operating heat.",
        "Both densities are the 200 MeV ledger.",
        "The chain stays at k_eff of 1.",
        "",
        "## Kilogram",
        "",
        "The unit-interval dressing of the kilogram ratio drops below 116 values.",
        "The positive dressing names it.",
    ]
    lines.extend(
        _figure_lines(
            report.kilogram,
            f"The Pu-239 kilogram over the U-235 kilogram is `{_num(report.kilogram.target)}`.",
            "The operating ratio stays the Pu-239 kilogram over the U-235 kilogram.",
        )
    )
    lines.extend(
        [
            "The two rows use the same enrichment input.",
            f"The named efficiency step is `{_num(report.efficiency_step)}`.",
            "The efficiency step stays the bounds score.",
            f"The fast delivery over the thermal delivery is `{_num(report.delivered_ratio)}`.",
            "That delivery ratio is the kilogram ratio times the named efficiency step.",
            "The fast return over the thermal return is that same product.",
            "Those two comparisons keep that product.",
            "The fast row stays on its own reference unit.",
            "The fast row's filed barrier is the sodium system.",
            "",
            "## Surry core",
            "",
            f"The UFSAR rates that core at `{_num(report.core_w)}` W.",
            f"The scored operating heat is `{_num(report.heat_w)}` W.",
            "The scored heat is the nameplate divided by the filed efficiency.",
            "The unit-interval dressing of that ratio drops below 116 values.",
            "The positive dressing names it.",
        ]
    )
    lines.extend(
        _figure_lines(
            report.core_ratio,
            (
                "The Surry core over the scored operating heat is "
                f"`{_num(report.core_ratio.target)}`."
            ),
            "The operating ratio stays the Surry core rating over the scored operating heat.",
        )
    )
    lines.extend(
        [
            "The Surry core stays the published rating.",
            "The fast power density stays an estimate.",
            "The thermal power density stays the filed order.",
            "The two power densities stay on their own bases.",
            "This score adds no joules.",
            "The catalog ceiling stays unrun.",
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
