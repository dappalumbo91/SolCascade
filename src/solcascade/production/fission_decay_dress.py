"""Dress the two decay limits the complement left outside.

The bare value is the recorded five-atom nearest. The factor is the
fine-structure leaf times one base piece, with either sign. A base piece is
an atom, its square root, or a power from 2 through 6. The 1952 fit stays
the cooling load. A dressing is not installed.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay_outside import build_decay_outside
from solcascade.production.fission_decay_seed import (
    GrammarSize,
    NamedGap,
    _REVIEW_GAP,
    _five_pieces,
    _nearest,
    _size,
)
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.refine import project_root

_REMAINING = ("one-second limit", "one-day limit")
_TENTH = 1.0e-4


@dataclass(frozen=True)
class DressFigure:
    label: str
    target: float
    anchor_name: str
    anchor_value: float
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    size: GrammarSize


@dataclass(frozen=True)
class DressReport:
    figures: tuple[DressFigure, ...]
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
    return project_root() / "runs" / "production" / "fission_decay_dress.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _base_pieces(module) -> tuple[tuple[str, float], ...]:
    rows = []
    for name, value in _five_pieces(module):
        if any(token in name for token in ("+", "*", "/")):
            continue
        if "-" in name.replace("^-", ""):
            continue
        rows.append((name, float(value)))
    if len(rows) != 60:
        raise RuntimeError("the base piece count changed")
    return tuple(rows)


def _dressings(anchor_name: str, anchor: float, alpha: float, base) -> tuple[tuple[str, float], ...]:
    rows = []
    for name, value in base:
        for sign in (1.0, -1.0):
            number = anchor * (1.0 + sign * alpha * value)
            if not 0.0 < number < 1.0:
                continue
            mark = "-" if sign < 0 else "+"
            rows.append((f"({anchor_name})*(1{mark}alpha*{name})", number))
    if not rows:
        raise RuntimeError(f"the {anchor_name} dressing is empty")
    return tuple(rows)


def build_decay_dress() -> DressReport:
    outside = build_decay_outside()
    remaining = tuple(figure for figure in outside.figures if not figure.complement.review)
    if tuple(figure.label for figure in remaining) != _REMAINING:
        raise RuntimeError("the remaining outside set changed")
    module, _digest = load_engine()
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    base = _base_pieces(module)
    figures: list[DressFigure] = []
    for figure in remaining:
        rows = _dressings(figure.five.name, figure.five.value, alpha, base)
        hit = _nearest(rows, figure.target)
        if hit.relative_gap <= helium:
            raise RuntimeError(f"the {figure.label} dressing meets the helium standard")
        size = _size(rows)
        inside_review = sum(
            abs(value - figure.target) / figure.target <= _REVIEW_GAP for _, value in rows
        )
        inside_tenth = sum(
            abs(value - figure.target) / figure.target <= _TENTH for _, value in rows
        )
        figures.append(
            DressFigure(
                label=figure.label,
                target=figure.target,
                anchor_name=figure.five.name,
                anchor_value=figure.five.value,
                dressing=hit,
                inside_review=inside_review,
                inside_tenth=inside_tenth,
                size=size,
            )
        )
    return DressReport(
        figures=tuple(figures),
        alpha=alpha,
        helium_mismatch=helium,
        second_watts=outside.second_watts,
        day_watts=outside.day_watts,
        two_watts=outside.two_watts,
        nuclear_scalar=outside.nuclear_scalar,
        thermodynamics_scalar=outside.thermodynamics_scalar,
        gamma_fsot=outside.gamma_fsot,
        valve_beta=outside.valve_beta,
        pin_prefix=outside.pin_prefix,
        digest_prefix=outside.digest_prefix,
    )


def _count_sentence(count: int, singular: str, plural: str) -> str:
    if count == 1:
        return f"One dressing of this anchor {singular}."
    return f"{count} dressings of this anchor {plural}."


def _figure_lines(figure: DressFigure) -> list[str]:
    lines = [
        f"The {figure.label} is `{_num(figure.target)}`.",
        f"Its five-atom anchor is `{figure.anchor_name}` at `{_num(figure.anchor_value)}`.",
        f"The nearest dressing is `{figure.dressing.name}` at `{_num(figure.dressing.value)}`.",
        f"Its relative gap is `{_num(figure.dressing.relative_gap)}`.",
    ]
    if figure.dressing.review:
        lines.append("The dressing meets the review bar.")
    else:
        lines.append("The dressing stays outside the review bar.")
    lines.append("The dressing stays wider than the helium quotient match.")
    lines.append(_count_sentence(figure.inside_review, "meets the review bar", "meet the review bar"))
    lines.append(
        _count_sentence(
            figure.inside_tenth,
            "is inside one part in ten thousand",
            "are inside one part in ten thousand",
        )
    )
    lines.append("The dressing is not installed.")
    return lines


def render() -> str:
    report = build_decay_dress()
    first, second = report.figures
    if first.size.expressions != second.size.expressions or first.size.inside != second.size.inside:
        raise RuntimeError("the two anchors grew different dressings")
    lines = [
        "# Fission decay dress",
        "",
        "The two limits here are the figures the complement left outside the review bar.",
        "The bare value is the recorded five-atom nearest.",
        "The factor is the fine-structure leaf times one base piece, with either sign.",
        "A base piece is an atom, its square root, or a power from 2 through 6 in either sign.",
        f"alpha in the dressing is the fine-structure leaf `{_num(report.alpha)}`.",
        f"Each anchor has `{first.size.expressions}` dressings inside (0, 1), `{first.size.inside}` distinct.",
        f"Their median relative spacing is `{_num(first.size.median_spacing)}`.",
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
