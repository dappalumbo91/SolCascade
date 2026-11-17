"""Name the thermionic quotients the paths score left beside the thermal plant.

The thermionic efficiency is the fabrication cap. The thermionic ceiling is
the small-class share. One expression names each of those shared numbers.
The delivered ratio, the enrichment return, the headroom, and the
file-ceiling ratio are the four quotients still open. The checked-in pulse
is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay_dress import _base_pieces
from solcascade.production.fission_decay_dress import _dressings as unit_dressings
from solcascade.production.fission_decay_heat import _dressings as positive_dressings
from solcascade.production.fission_decay_heat import _nearest_row, _positive_nearest, _positive_size
from solcascade.production.fission_decay_seed import NamedGap, _five_pieces, _nearest, _size
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.fission_paths import FAST_ID, THERMAL_ID, THERMIONIC_ID, build_paths
from solcascade.production.fission_process import FABRICATION_CAP_PERCENT, NAMEPLATE_MWE
from solcascade.production.fission_safety import SMALL_CLASS_MWE
from solcascade.production.refine import catalog_record, project_root

_COUNT = 116
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4
_DELIVERED_ANCHOR = "e*phi^-6"
_DELIVERED_STEPS = (
    "(e*phi^-6)*(1+alpha*pi^-3)",
    "((e*phi^-6)*(1+alpha*pi^-3))*(1-alpha*pi^-5)",
    "(((e*phi^-6)*(1+alpha*pi^-3))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
)
_RETURN_ANCHOR = "e^5*G^3"
_RETURN_STEPS = (
    "(e^5*G^3)*(1-alpha*pi^-3)",
    "((e^5*G^3)*(1-alpha*pi^-3))*(1+alpha*pi^-6)",
)
_HEAD_ANCHOR = "gamma^-2+gamma^-2"
_HEAD_STEPS = ("(gamma^-2+gamma^-2)*(1-alpha*gamma^5)",)
_CEILING_ANCHOR = "gamma+gamma^2"
_CEILING_STEPS = (
    "(gamma+gamma^2)*(1-alpha*gamma^3)",
    "((gamma+gamma^2)*(1-alpha*gamma^3))*(1-alpha*pi^-5)",
)


@dataclass(frozen=True)
class ThermionicStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class ThermionicFigure:
    anchor: NamedGap
    steps: tuple[ThermionicStep, ...]
    target: float


@dataclass(frozen=True)
class ThermionicReport:
    delivered: ThermionicFigure
    enrichment_return: ThermionicFigure
    headroom: ThermionicFigure
    ceiling_ratio: ThermionicFigure
    efficiency: float
    ceiling: float
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_thermionic.md"


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


def _walk(target, anchor_name, expected, alpha, helium, five, base, positive: bool) -> ThermionicFigure:
    if positive:
        name, value, anchor_gap = _positive_nearest(five, target)
        anchor = NamedGap(name, value, anchor_gap, anchor_gap <= _REVIEW_GAP)
    else:
        anchor = _nearest(five, target)
        name, value = anchor.name, anchor.value
    if anchor.name != anchor_name:
        raise RuntimeError(f"the thermionic anchor changed to {anchor.name}")
    if anchor.relative_gap <= helium:
        raise RuntimeError("the thermionic anchor already meets helium")
    if positive and not target > 1.0:
        dropped = len(unit_dressings(anchor.name, anchor.value, alpha, base))
        if dropped >= _COUNT:
            raise RuntimeError("the unit grammar still has 116 dressings")
    steps: list[ThermionicStep] = []
    for step_name in expected:
        rows = (positive_dressings if positive else unit_dressings)(name, value, alpha, base)
        if positive:
            hit_name, hit_value, hit_gap = _nearest_row(rows, target)
            hit = NamedGap(hit_name, hit_value, hit_gap, hit_gap <= _REVIEW_GAP)
        else:
            hit = _nearest(rows, target)
        if hit.name != step_name:
            raise RuntimeError(f"the thermionic dressing changed to {hit.name}")
        if hit.relative_gap >= _gap(value, target):
            raise RuntimeError("the thermionic dressing did not tighten")
        size = (_positive_size if positive else _size)(rows)
        if size.expressions != _COUNT or size.inside != _COUNT:
            raise RuntimeError("the thermionic dressing count changed")
        if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
            raise RuntimeError("the thermionic spacing is not the base-piece spacing")
        inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
        if inside_helium > 1:
            raise RuntimeError("several thermionic dressings meet the helium standard")
        steps.append(
            ThermionicStep(
                dressing=hit,
                inside_review=sum(_gap(number, target) <= _REVIEW_GAP for _row, number in rows),
                inside_tenth=sum(_gap(number, target) <= _TENTH for _row, number in rows),
                inside_helium=inside_helium,
                second_gap=_second_gap(rows, target, hit.relative_gap),
                expressions=size.expressions,
                median_spacing=size.median_spacing,
            )
        )
        name, value = hit.name, hit.value
    closing = steps[-1]
    if closing.inside_helium != 1 or closing.dressing.relative_gap > helium:
        raise RuntimeError("the closing thermionic figure misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError("a second thermionic dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier thermionic figure already meets helium")
    return ThermionicFigure(anchor, tuple(steps), target)


def build_thermionic() -> ThermionicReport:
    paths = build_paths()
    module, digest = load_engine()
    if digest[:16] != paths.pin_prefix or paths.digest_prefix != paths.pin_prefix:
        raise RuntimeError("the thermionic pin does not match the plant pin")
    thermal = paths.thermal
    row = paths.thermionic
    thermal_record = catalog_record(THERMAL_ID)
    if "k_eff is held at 1" not in str(thermal_record["governing_equation"]):
        raise RuntimeError("the governing equation does not hold k_eff at 1")
    fast = catalog_record(FAST_ID)
    if "Sodium" not in str(fast.get("engineering_barrier") or ""):
        raise RuntimeError("the fast barrier does not name sodium")
    notes = catalog_record(THERMIONIC_ID).get("notes") or ""
    if "TOPAZ" not in str(notes):
        raise RuntimeError("the thermionic notes do not name TOPAZ")
    cap = FABRICATION_CAP_PERCENT / 100.0
    small = SMALL_CLASS_MWE / NAMEPLATE_MWE
    if row.efficiency != cap:
        raise RuntimeError("the thermionic efficiency left the fabrication cap")
    if row.ceiling != small:
        raise RuntimeError("the thermionic ceiling left the small-class share")
    delivered = row.delivered_j / thermal.delivered_j
    if _gap(delivered, row.efficiency / thermal.efficiency) > 1.0e-12:
        raise RuntimeError("the delivered ratio left the efficiency quotient")
    headroom = row.ceiling / row.efficiency
    if _gap(headroom, 6.0) > 1.0e-15:
        raise RuntimeError("the thermionic headroom left six")
    ceiling_ratio = row.ceiling_delivered_j / thermal.delivered_j
    if _gap(ceiling_ratio, row.ceiling / thermal.efficiency) > 1.0e-12:
        raise RuntimeError("the ceiling ratio left the efficiency quotient")
    if not 0.0 < delivered < ceiling_ratio < 1.0:
        raise RuntimeError("the thermionic ceiling meets the thermal delivery")
    if _gap(row.output_over_input, thermal.output_over_input * delivered) > 1.0e-12:
        raise RuntimeError("the thermionic return left the delivered ratio of the thermal return")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    five = _five_pieces(module)
    base = _base_pieces(module)
    return ThermionicReport(
        delivered=_walk(delivered, _DELIVERED_ANCHOR, _DELIVERED_STEPS, alpha, helium, five, base, True),
        enrichment_return=_walk(
            row.output_over_input, _RETURN_ANCHOR, _RETURN_STEPS, alpha, helium, five, base, True
        ),
        headroom=_walk(headroom, _HEAD_ANCHOR, _HEAD_STEPS, alpha, helium, five, base, True),
        ceiling_ratio=_walk(
            ceiling_ratio, _CEILING_ANCHOR, _CEILING_STEPS, alpha, helium, five, base, True
        ),
        efficiency=row.efficiency,
        ceiling=row.ceiling,
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=paths.nuclear_scalar,
        thermodynamics_scalar=paths.thermodynamics_scalar,
        gamma_fsot=paths.gamma_fsot,
        valve_beta=paths.valve_beta,
        pin_prefix=paths.pin_prefix,
        digest_prefix=paths.digest_prefix,
    )


def _step_lines(step: ThermionicStep, closing: bool) -> list[str]:
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


def _figure_lines(figure: ThermionicFigure, lead: str, role: str) -> list[str]:
    lines = [
        lead,
        f"The five-atom nearest is `{figure.anchor.name}` at `{_num(figure.anchor.value)}`.",
        f"The anchor gap is `{_num(figure.anchor.relative_gap)}`.",
    ]
    if len(figure.steps) == 1:
        lines.append("The first factor closes that quotient.")
    for step in figure.steps[:-1]:
        lines.extend(_step_lines(step, False))
    lines.extend(_step_lines(figure.steps[-1], True))
    closing = figure.steps[-1].dressing.value
    if _num(closing) == _num(figure.target):
        lines.append(f"The closing value agrees with `{_num(figure.target)}` at six figures.")
    lines.append("The closing value names that quotient.")
    lines.append(role)
    return lines


def render() -> str:
    report = build_thermionic()
    lines = [
        "# Fission thermionic",
        "",
        "The paths score leaves the thermionic row beside the thermal plant on one U-235 kilogram.",
        "The thermionic notes name TOPAZ as a flown thermionic converter.",
        "The thermionic efficiency is the fabrication cap.",
        "The thermionic ceiling is the small-class share.",
        "One expression names each shared number.",
        "FSOT names the four quotients those shared numbers still leave open.",
        "The chain stays at k_eff of 1.",
        "",
        "## Delivery",
        "",
        "The unit-interval dressing of the delivered ratio drops below 116 values.",
        "The positive dressing names it.",
    ]
    lines.extend(
        _figure_lines(
            report.delivered,
            f"The thermionic delivery over the thermal delivery is `{_num(report.delivered.target)}`.",
            "The operating ratio stays the thermionic delivery over the thermal delivery.",
        )
    )
    lines.extend(
        [
            "The thermionic row delivers less electricity than the thermal plant.",
            "",
            "## Return",
            "",
            "Ratios that leave the unit interval use the positive dressing from the heat-ratio score.",
        ]
    )
    lines.extend(
        _figure_lines(
            report.enrichment_return,
            (
                "The thermionic delivery over the filed enrichment input is "
                f"`{_num(report.enrichment_return.target)}`."
            ),
            "The operating return stays the thermionic delivery over the filed enrichment input.",
        )
    )
    lines.extend(
        [
            "This ratio compares those two filed energies.",
            "",
            "## Headroom",
            "",
        ]
    )
    lines.extend(
        _figure_lines(
            report.headroom,
            (
                "The thermionic ceiling divided by the thermionic efficiency is "
                f"`{_num(report.headroom.target)}`."
            ),
            "The operating headroom stays the thermionic ceiling over the thermionic efficiency.",
        )
    )
    lines.extend(
        [
            "The fast row's filed barrier is the sodium system.",
            "The fast row stays on its own reference unit.",
            "",
            "## Ceiling",
            "",
            "The unit-interval dressing of the ceiling ratio drops below 116 values.",
            "The positive dressing names it.",
        ]
    )
    lines.extend(
        _figure_lines(
            report.ceiling_ratio,
            (
                "The thermionic file-ceiling delivery over the thermal delivery is "
                f"`{_num(report.ceiling_ratio.target)}`."
            ),
            "The operating ceiling ratio stays the thermionic file-ceiling delivery over the thermal delivery.",
        )
    )
    lines.extend(
        [
            "Even at its file ceiling the thermionic row delivers less than the thermal plant.",
            "This score adds no joules.",
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
