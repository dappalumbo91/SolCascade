"""Score the decay figures the five-atom grammar left outside the review bar.

The added operation is the complement already used on the balance share:
one minus a five-atom piece. The fine-structure leaf is scored beside it.
Neither is installed, and neither multiplies the watts. The 1952 fit stays
the cooling load.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay_seed import (
    GrammarSize,
    NamedGap,
    _REVIEW_GAP,
    _five_pieces,
    _nearest,
    _size,
    build_decay_seed,
)
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.refine import project_root

_OUTSIDE = (
    "one-second limit",
    "one-day limit",
    "two-day limit",
    "two-day drop",
)


@dataclass(frozen=True)
class OutsideFigure:
    label: str
    target: float
    five: NamedGap
    complement: NamedGap
    alpha_gap: float


@dataclass(frozen=True)
class OutsideReport:
    figures: tuple[OutsideFigure, ...]
    complement_size: GrammarSize
    alpha: float
    helium_mismatch: float
    second_watts: float
    day_watts: float
    two_watts: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_decay_outside.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _complements(module) -> tuple[tuple[str, float], ...]:
    rows = []
    for name, value in _five_pieces(module):
        number = float(value)
        if 0.0 < number < 1.0:
            rows.append((f"1-({name})", 1.0 - number))
    if not rows:
        raise RuntimeError("the complement grammar is empty")
    return tuple(rows)


def build_decay_outside() -> OutsideReport:
    seed = build_decay_seed()
    labels = tuple(figure.label for figure in seed.figures if not figure.five.review)
    if labels != _OUTSIDE:
        raise RuntimeError("the outside set changed")
    module, _digest = load_engine()
    complements = _complements(module)
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    if not 0.0 < alpha < 1.0:
        raise RuntimeError("the fine-structure leaf is outside (0, 1)")
    figures: list[OutsideFigure] = []
    for figure in seed.figures:
        if figure.five.review:
            continue
        complement = _nearest(complements, figure.target)
        alpha_gap = abs(alpha - figure.target) / figure.target
        if complement.relative_gap <= helium or alpha_gap <= helium:
            raise RuntimeError(f"the {figure.label} meets the helium standard")
        figures.append(
            OutsideFigure(
                label=figure.label,
                target=figure.target,
                five=figure.five,
                complement=complement,
                alpha_gap=alpha_gap,
            )
        )
    return OutsideReport(
        figures=tuple(figures),
        complement_size=_size(complements),
        alpha=alpha,
        helium_mismatch=helium,
        second_watts=seed.second_watts,
        day_watts=seed.day_watts,
        two_watts=seed.two_watts,
        nuclear_scalar=seed.nuclear_scalar,
        thermodynamics_scalar=seed.thermodynamics_scalar,
        gamma_fsot=seed.gamma_fsot,
        valve_beta=seed.valve_beta,
        pin_prefix=seed.pin_prefix,
        digest_prefix=seed.digest_prefix,
    )


def _figure_lines(figure: OutsideFigure) -> list[str]:
    lines = [
        f"The {figure.label} is `{_num(figure.target)}`.",
        f"The five-atom nearest remains `{figure.five.name}` at `{_num(figure.five.value)}`, relative gap `{_num(figure.five.relative_gap)}`.",
        f"The complement nearest is `{figure.complement.name}` at `{_num(figure.complement.value)}`, relative gap `{_num(figure.complement.relative_gap)}`.",
        f"The fine-structure leaf gap is `{_num(figure.alpha_gap)}`.",
    ]
    if figure.complement.review:
        lines.append("The complement meets the review bar.")
    else:
        lines.append("The complement stays outside the review bar.")
    if figure.alpha_gap <= _REVIEW_GAP:
        lines.append("The fine-structure leaf meets the review bar.")
    else:
        lines.append("The fine-structure leaf stays outside the review bar.")
    lines.append("The expression is not installed.")
    return lines


def render() -> str:
    report = build_decay_outside()
    lines = [
        "# Fission decay outside",
        "",
        "These are the four decay figures the five-atom grammar left outside the review bar.",
        "The added operation is the complement of a five-atom piece.",
        f"The complement grammar has `{report.complement_size.expressions}` expressions, `{report.complement_size.inside}` inside (0, 1).",
        f"Its median relative spacing is `{_num(report.complement_size.median_spacing)}`.",
        f"The fine-structure leaf is `{_num(report.alpha)}`.",
        f"The helium quotient match is `{_num(report.helium_mismatch)}`.",
        "Every scored gap stays wider than that match.",
        "",
        "## Figures",
        "",
    ]
    for figure in report.figures:
        lines.extend(_figure_lines(figure))
        lines.append("")
    lines.extend(
        [
            "## Ledger",
            "",
            "No seed is installed as a decay-heat fraction.",
            "The 1952 cooling loads are unchanged.",
            f"The one-second load remains `{_num(report.second_watts)}` W.",
            f"The one-day load remains `{_num(report.day_watts)}` W.",
            f"The two-day load remains `{_num(report.two_watts)}` W.",
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
