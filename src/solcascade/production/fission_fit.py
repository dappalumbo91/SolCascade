"""Name the 1952 long-weight and exponent coefficients.

The lead coefficient is the ten percent band. Equation (2) is not replaced.
The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay import _EXPONENT, _LEAD, _LONG_WEIGHT, stated_band
from solcascade.production.fission_decay_dress import _base_pieces
from solcascade.production.fission_decay_dress import _dressings as unit_dressings
from solcascade.production.fission_decay_heat import _dressings as positive_dressings
from solcascade.production.fission_decay_heat import _nearest_row, _positive_nearest, _positive_size
from solcascade.production.fission_decay_seed import NamedGap, _five_pieces, _size
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.fission_process import build_process
from solcascade.production.fission_service import _COOLDOWN_S, _SERVICE_S
from solcascade.production.refine import project_root

_COUNT = 116
_SPACING = "0.000511364"
_REVIEW_GAP = 1.0e-3
_TENTH = 1.0e-4
_WEIGHT_ANCHOR = "pi^-2+G^3"
_WEIGHT_STEPS = (
    "(pi^-2+G^3)*(1+alpha*pi^-3)",
    "((pi^-2+G^3)*(1+alpha*pi^-3))*(1-alpha*e^-6)",
)
_EXPONENT_ANCHOR = "e^-2+gamma^5"
_EXPONENT_STEPS = (
    "(e^-2+gamma^5)*(1+alpha*phi^-2)",
    "((e^-2+gamma^5)*(1+alpha*phi^-2))*(1+alpha*e^-4)",
    "(((e^-2+gamma^5)*(1+alpha*phi^-2))*(1+alpha*e^-4))*(1+alpha*pi^-5)",
    "((((e^-2+gamma^5)*(1+alpha*phi^-2))*(1+alpha*e^-4))*(1+alpha*pi^-5))*(1+alpha*pi^-6)",
)


@dataclass(frozen=True)
class FitStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    second_gap: float
    expressions: int
    median_spacing: float


@dataclass(frozen=True)
class FitFigure:
    anchor: NamedGap
    steps: tuple[FitStep, ...]
    target: float


@dataclass(frozen=True)
class FitReport:
    weight: FitFigure
    exponent: FitFigure
    lead: float
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_fit.md"


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


def _walk(target, anchor_name, expected, alpha, helium, five, base) -> FitFigure:
    name, value, anchor_gap = _positive_nearest(five, target)
    anchor = NamedGap(name, value, anchor_gap, anchor_gap <= _REVIEW_GAP)
    if anchor.name != anchor_name:
        raise RuntimeError("the fit anchor changed")
    steps: list[FitStep] = []
    for step_name in expected:
        rows = positive_dressings(name, value, alpha, base)
        hit_name, hit_value, hit_gap = _nearest_row(rows, target)
        hit = NamedGap(hit_name, hit_value, hit_gap, hit_gap <= _REVIEW_GAP)
        if hit.name != step_name:
            raise RuntimeError("the fit dressing changed")
        if hit.relative_gap >= _gap(value, target):
            raise RuntimeError("the fit dressing did not tighten")
        size = _positive_size(rows)
        if size.expressions != _COUNT or size.inside != _COUNT:
            raise RuntimeError("the fit dressing count changed")
        if _num(size.median_spacing) != _SPACING or size.median_spacing <= helium:
            raise RuntimeError("the fit spacing is not the base-piece spacing")
        inside_helium = sum(_gap(number, target) <= helium for _row, number in rows)
        if inside_helium > 1:
            raise RuntimeError("several fit dressings meet the helium standard")
        steps.append(
            FitStep(
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
        raise RuntimeError("the closing fit figure misses the helium standard")
    if closing.second_gap <= helium:
        raise RuntimeError("a second fit dressing also meets helium")
    if any(step.inside_helium != 0 for step in steps[:-1]):
        raise RuntimeError("an earlier fit figure already meets helium")
    return FitFigure(anchor, tuple(steps), target)


def build_fit() -> FitReport:
    process = build_process()
    module, digest = load_engine()
    if digest[:16] != process.pin_prefix:
        raise RuntimeError("the fit pin does not match the plant pin")
    if _LEAD != stated_band(1.0e4):
        raise RuntimeError("the equation (2) lead left the ten percent band")
    if _EXPONENT >= 0.0:
        raise RuntimeError("the equation (2) exponent left its minus sign")
    if _SERVICE_S / _COOLDOWN_S != abs(_EXPONENT):
        raise RuntimeError("four hours over twenty hours left the exponent magnitude")
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    five = _five_pieces(module)
    base = _base_pieces(module)
    weight = _walk(_LONG_WEIGHT, _WEIGHT_ANCHOR, _WEIGHT_STEPS, alpha, helium, five, base)
    exponent = _walk(abs(_EXPONENT), _EXPONENT_ANCHOR, _EXPONENT_STEPS, alpha, helium, five, base)
    for figure, label in ((weight, "long weight"), (exponent, "exponent")):
        unit_rows = unit_dressings(figure.anchor.name, figure.anchor.value, alpha, base)
        if _size(unit_rows).expressions == _COUNT:
            raise RuntimeError(f"the {label} kept the unit-interval count")
    if _gap(weight.steps[-1].dressing.value, _LONG_WEIGHT) > helium:
        raise RuntimeError("the named long weight left equation (2)")
    if _gap(exponent.steps[-1].dressing.value, abs(_EXPONENT)) > helium:
        raise RuntimeError("the named exponent left equation (2)")
    return FitReport(
        weight=weight,
        exponent=exponent,
        lead=_LEAD,
        alpha=alpha,
        helium_mismatch=helium,
        nuclear_scalar=process.nuclear_scalar,
        thermodynamics_scalar=process.thermodynamics_scalar,
        gamma_fsot=process.gamma_fsot,
        valve_beta=process.valve_beta,
        pin_prefix=process.pin_prefix,
        digest_prefix=digest[:16],
    )


def _step_lines(step: FitStep, closing: bool) -> list[str]:
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


def _figure_lines(title: str, intro: list[str], figure: FitFigure) -> list[str]:
    lines = ["", f"## {title}", "", *intro]
    lines.append(
        f"The positive five-atom nearest is `{figure.anchor.name}` at `{_num(figure.anchor.value)}`."
    )
    lines.append(f"The anchor gap is `{_num(figure.anchor.relative_gap)}`.")
    for step in figure.steps[:-1]:
        lines.extend(_step_lines(step, False))
    lines.extend(_step_lines(figure.steps[-1], True))
    return lines


def render() -> str:
    report = build_fit()
    lines = [
        "# Fission fit",
        "",
        "Equation (2) of the 1952 fit stays the witness.",
        "The lead coefficient is the ten percent band.",
        "The unit-interval dressing of each coefficient drops below 116 values.",
        "The positive dressing from the heat-ratio score names these coefficients.",
    ]
    lines.extend(
        _figure_lines(
            "Long weight",
            [f"The filed coefficient is `{_num(report.weight.target)}`."],
            report.weight,
        )
    )
    lines.extend(
        _figure_lines(
            "Exponent",
            [
                f"The filed magnitude is `{_num(report.exponent.target)}`.",
                "Four hours over twenty hours is that same number.",
            ],
            report.exponent,
        )
    )
    lines.extend(
        [
            "This score does not replace equation (2).",
            "This run does not install a seed as a new decay fraction or a new plant row.",
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
