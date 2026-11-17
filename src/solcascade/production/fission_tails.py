"""Name the tails percent and the tails fraction of the 1000 MWe example.

The overview states 0.22 percent tails with the 4.5 percent enrichment.
The percent uses the positive dressing because the unit list has 115 values.
The fraction keeps the unit dressing. The checked-in pulse is not rewritten.
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
from solcascade.production.fission_process import (
    ANNUAL_SWU,
    ENRICHED_PERCENT,
    FABRICATION_CAP_PERCENT,
    MATERIAL_BALANCE_SOURCE,
    RELOAD_TONNES,
)
from solcascade.production.plants import FAST_ID, THERMAL_ID, build_plants
from solcascade.production.refine import catalog_record, project_root

_COUNT = 116
_UNIT_SHORT = 115
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4
TAILS_PERCENT = 0.22
_PERCENT_ANCHOR = "phi^-2*gamma"
_PERCENT_STEPS = (
    "(phi^-2*gamma)*(1-alpha*gamma^2)",
    "((phi^-2*gamma)*(1-alpha*gamma^2))*(1+alpha*gamma^6)",
)
_FRACTION_ANCHOR = "pi^-5-pi^-6"
_FRACTION_STEPS = (
    "(pi^-5-pi^-6)*(1-alpha*G^-6)",
    "((pi^-5-pi^-6)*(1-alpha*G^-6))*(1-alpha*pi^-5)",
    "(((pi^-5-pi^-6)*(1-alpha*G^-6))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
)


@dataclass(frozen=True)
class TailsStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class TailsFigure:
    anchor: NamedGap
    steps: tuple[TailsStep, ...]
    target: float


@dataclass(frozen=True)
class TailsReport:
    percent: TailsFigure
    fraction: TailsFigure
    tails_percent: float
    enriched_percent: float
    annual_swu: float
    reload_tonnes: float
    fabrication_cap_percent: float
    source: str
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_tails.md"


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


def _record(rows, target: float, hit: NamedGap, helium: float, size) -> TailsStep:
    if size.expressions != _COUNT or size.inside != _COUNT:
        raise RuntimeError("the tails dressing count changed")
    if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
        raise RuntimeError("the tails spacing is not the base-piece spacing")
    inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
    if inside_helium > 1:
        raise RuntimeError("several tails dressings meet the helium standard")
    return TailsStep(
        dressing=hit,
        inside_review=sum(_gap(number, target) <= _REVIEW_GAP for _row, number in rows),
        inside_tenth=sum(_gap(number, target) <= _TENTH for _row, number in rows),
        inside_helium=inside_helium,
        second_gap=_second_gap(rows, target, hit.relative_gap),
        expressions=size.expressions,
        median_spacing=size.median_spacing,
    )


def _walk(target, anchor_name, expected, alpha, helium, five, base, positive: bool) -> TailsFigure:
    if positive:
        name, value, anchor_gap = _positive_nearest(five, target)
    else:
        anchor = _nearest(five, target)
        name, value, anchor_gap = anchor.name, anchor.value, anchor.relative_gap
    if name != anchor_name:
        raise RuntimeError(f"the tails anchor changed to {name}")
    if anchor_gap <= helium:
        raise RuntimeError("the tails anchor already meets helium")
    anchor_gap_named = NamedGap(name, value, anchor_gap, anchor_gap <= _REVIEW_GAP)
    steps: list[TailsStep] = []
    for step_name in expected:
        unit_rows = unit_dressings(name, value, alpha, base)
        if positive:
            if len(unit_rows) != _UNIT_SHORT:
                raise RuntimeError("the tails percent left 115 unit dressings")
            rows = positive_dressings(name, value, alpha, base)
            hit_name, hit_value, hit_gap = _nearest_row(rows, target)
            hit = NamedGap(hit_name, hit_value, hit_gap, hit_gap <= _REVIEW_GAP)
            size = _positive_size(rows)
        else:
            if len(unit_rows) != _COUNT:
                raise RuntimeError("the tails fraction left 116 unit dressings")
            rows = unit_rows
            hit = _nearest(rows, target)
            size = _size(rows)
        if hit.name != step_name:
            raise RuntimeError(f"the tails dressing changed to {hit.name}")
        if hit.relative_gap >= _gap(value, target):
            raise RuntimeError("the tails dressing did not tighten")
        recorded = _record(rows, target, hit, helium, size)
        steps.append(recorded)
        name, value = hit.name, hit.value
    closing = steps[-1]
    if closing.inside_helium != 1 or closing.dressing.relative_gap > helium:
        raise RuntimeError("the closing tails figure misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError("a second tails dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier tails figure already meets helium")
    return TailsFigure(anchor_gap_named, tuple(steps), target)


def build_tails() -> TailsReport:
    plants = build_plants()
    module, digest = load_engine()
    if digest[:16] != plants.pin_prefix or plants.digest_prefix != plants.pin_prefix:
        raise RuntimeError("the tails pin does not match the plant pin")
    record = catalog_record(THERMAL_ID)
    if "k_eff is held at 1" not in str(record["governing_equation"]):
        raise RuntimeError("the governing equation does not hold k_eff at 1")
    if "Sodium" not in str(catalog_record(FAST_ID).get("engineering_barrier") or ""):
        raise RuntimeError("the fast barrier does not name sodium")
    if TAILS_PERCENT != 0.22:
        raise RuntimeError("the tails percent left 0.22")
    if not TAILS_PERCENT < ENRICHED_PERCENT < FABRICATION_CAP_PERCENT:
        raise RuntimeError("the tails percent left its place under the enriched product and the cap")
    if ANNUAL_SWU != 182_000.0 or RELOAD_TONNES != 24.3:
        raise RuntimeError("the 1000 MWe example left 182000 SWU and 24.3 t")
    if FABRICATION_CAP_PERCENT != 5.0:
        raise RuntimeError("the fabrication cap left 5 percent")
    fraction = TAILS_PERCENT / 100.0
    if not 0.0 < fraction < 1.0:
        raise RuntimeError("the tails fraction left the unit interval")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    five = _five_pieces(module)
    base = _base_pieces(module)
    percent = _walk(TAILS_PERCENT, _PERCENT_ANCHOR, _PERCENT_STEPS, alpha, helium, five, base, True)
    fraction_figure = _walk(
        fraction, _FRACTION_ANCHOR, _FRACTION_STEPS, alpha, helium, five, base, False
    )
    if len(percent.steps) != 2 or len(fraction_figure.steps) != 3:
        raise RuntimeError("the tails closes left their factor counts")
    if _num(percent.steps[-1].dressing.value) != _num(percent.target):
        raise RuntimeError("the tails percent left six-figure agreement")
    if _num(fraction_figure.steps[-1].dressing.value) == _num(fraction_figure.target):
        raise RuntimeError("the tails fraction print met the page figure at six figures")
    return TailsReport(
        percent=percent,
        fraction=fraction_figure,
        tails_percent=TAILS_PERCENT,
        enriched_percent=ENRICHED_PERCENT,
        annual_swu=ANNUAL_SWU,
        reload_tonnes=RELOAD_TONNES,
        fabrication_cap_percent=FABRICATION_CAP_PERCENT,
        source=MATERIAL_BALANCE_SOURCE,
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=plants.nuclear_scalar,
        thermodynamics_scalar=plants.thermodynamics_scalar,
        gamma_fsot=plants.gamma_fsot,
        valve_beta=plants.valve_beta,
        pin_prefix=plants.pin_prefix,
        digest_prefix=plants.digest_prefix,
    )


def _step_lines(step: TailsStep, closing: bool) -> list[str]:
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


def _figure_lines(figure: TailsFigure, lead: str, role: str) -> list[str]:
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
    report = build_tails()
    lines = [
        "# Fission tails",
        "",
        f"The material balance is {report.source}.",
        "That page's typical year is a 1000 MWe light-water reactor at 4.5 percent enrichment, 0.22 percent tails, 45 GWd/t, and 33 percent.",
        "FSOT names the tails percent and the tails fraction.",
        "The chain stays at k_eff of 1.",
        "",
        "## Tails percent",
        "",
        f"The unit-interval dressing of the tails percent has `{_UNIT_SHORT}` values.",
        "The positive dressing names it.",
    ]
    lines.extend(
        _figure_lines(
            report.percent,
            f"The tails percent is `{_num(report.percent.target)}`.",
            "The operating percent stays the printed tails percent.",
        )
    )
    lines.extend(
        [
            "The overview prints no further digit on that percent.",
            "",
            "## Tails fraction",
            "",
            "The tails fraction is the tails percent over 100.",
            "The unit dressing names it.",
        ]
    )
    lines.extend(
        _figure_lines(
            report.fraction,
            f"The tails fraction is `{_num(report.fraction.target)}`.",
            "The operating fraction stays the tails percent over 100.",
        )
    )
    lines.extend(
        [
            f"The closing print is `{_num(report.fraction.steps[-1].dressing.value)}`.",
            "The overview prints the fraction as the tails percent over 100.",
            "",
            "## Page figures",
            "",
            f"The enriched percent stays `{_num(report.enriched_percent)}`.",
            "The enriched percent stays the stack score.",
            "The enriched fraction stays the held pair.",
            f"The annual separative work stays `{_num(report.annual_swu)}` SWU.",
            f"The reload stays `{_num(report.reload_tonnes)}` t.",
            "Those figures stay the process score.",
            "The fabrication cap stays 5 percent.",
            "The fabrication page's about 4.8 percent stays approximate.",
            "This score adds no joules.",
            "The catalog ceiling stays unrun.",
            "The 60 GWd/t stack stays unrun.",
            "The fast row's filed barrier is the sodium system.",
            "The published stamp records these two closes.",
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
