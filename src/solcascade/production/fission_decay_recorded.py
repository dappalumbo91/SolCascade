"""Dress the six decay names the first dressing left untouched.

The bare value is the recorded nearest already inside the review bar. The
factor is the fine-structure leaf times one base piece, with either sign.
A unique match inside the helium quotient standard is recorded. The 1952
fit stays the cooling load.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay_dress import (
    _base_pieces,
    _dressings,
    build_decay_dress,
)
from solcascade.production.fission_decay_outside import build_decay_outside
from solcascade.production.fission_decay_seed import (
    GrammarSize,
    NamedGap,
    _REVIEW_GAP,
    _nearest,
    _size,
    build_decay_seed,
)
from solcascade.production.fission_observables import _adopted_alpha, _helium_mismatch
from solcascade.production.refine import project_root

_DRESSED = ("one-second limit", "one-day limit")
_FIGURES = (
    "one-second fraction",
    "one-day fraction",
    "two-day fraction",
    "two-day limit",
    "one-day drop",
    "two-day drop",
)
_HELIUM = (
    "one-second fraction",
    "two-day limit",
    "two-day drop",
)
_TENTH = 1.0e-4


@dataclass(frozen=True)
class RecordedFigure:
    label: str
    target: float
    anchor: NamedGap
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    size: GrammarSize


@dataclass(frozen=True)
class RecordedReport:
    figures: tuple[RecordedFigure, ...]
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
    return project_root() / "runs" / "production" / "fission_decay_recorded.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _recorded_anchor(figure, outside) -> NamedGap:
    candidates = [figure.five]
    outside_figure = outside.get(figure.label)
    if outside_figure is not None:
        candidates.append(outside_figure.complement)
    anchor = min(candidates, key=lambda item: (item.relative_gap, item.name))
    if not anchor.review:
        raise RuntimeError(f"the {figure.label} has no recorded nearest inside the review bar")
    return anchor


def build_decay_recorded() -> RecordedReport:
    seed = build_decay_seed()
    outside = {figure.label: figure for figure in build_decay_outside().figures}
    dress = build_decay_dress()
    if tuple(figure.label for figure in dress.figures) != _DRESSED:
        raise RuntimeError("the dressed set changed")
    module, _digest = load_engine()
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    base = _base_pieces(module)
    dressed = set(_DRESSED)
    figures: list[RecordedFigure] = []
    for figure in seed.figures:
        if figure.label in dressed:
            continue
        anchor = _recorded_anchor(figure, outside)
        rows = _dressings(anchor.name, anchor.value, alpha, base)
        hit = _nearest(rows, figure.target)
        if hit.relative_gap >= anchor.relative_gap:
            raise RuntimeError(f"the {figure.label} dressing did not tighten")
        if hit.relative_gap > _REVIEW_GAP:
            raise RuntimeError(f"the {figure.label} dressing leaves the review bar")
        size = _size(rows)
        inside_review = sum(
            abs(value - figure.target) / figure.target <= _REVIEW_GAP for _, value in rows
        )
        inside_tenth = sum(
            abs(value - figure.target) / figure.target <= _TENTH for _, value in rows
        )
        inside_helium = sum(
            abs(value - figure.target) / figure.target <= helium for _, value in rows
        )
        if inside_helium > 1:
            raise RuntimeError(f"several {figure.label} dressings meet the helium standard")
        if (hit.relative_gap <= helium) != (inside_helium == 1):
            raise RuntimeError(f"the {figure.label} helium count does not match its nearest")
        figures.append(
            RecordedFigure(
                label=figure.label,
                target=figure.target,
                anchor=anchor,
                dressing=hit,
                inside_review=inside_review,
                inside_tenth=inside_tenth,
                inside_helium=inside_helium,
                size=size,
            )
        )
    if tuple(figure.label for figure in figures) != _FIGURES:
        raise RuntimeError("the recorded figure set changed")
    met = tuple(figure.label for figure in figures if figure.dressing.relative_gap <= helium)
    if met != _HELIUM:
        raise RuntimeError("the helium set changed")
    return RecordedReport(
        figures=tuple(figures),
        alpha=alpha,
        helium_mismatch=helium,
        second_watts=dress.second_watts,
        day_watts=dress.day_watts,
        two_watts=dress.two_watts,
        nuclear_scalar=dress.nuclear_scalar,
        thermodynamics_scalar=dress.thermodynamics_scalar,
        gamma_fsot=dress.gamma_fsot,
        valve_beta=dress.valve_beta,
        pin_prefix=dress.pin_prefix,
        digest_prefix=dress.digest_prefix,
    )


def _count_sentence(count: int, singular: str, plural: str) -> str:
    if count == 1:
        return f"One dressing of this anchor {singular}."
    return f"{count} dressings of this anchor {plural}."


def _helium_sentence(count: int) -> str:
    if count == 0:
        return "No dressing of this anchor meets the helium quotient match."
    if count == 1:
        return "One dressing of this anchor meets the helium quotient match."
    raise RuntimeError("several dressings meet the helium standard")


def _piece(name: str) -> tuple[str, str]:
    if "1-alpha*" in name:
        sign = "minus"
    elif "1+alpha*" in name:
        sign = "plus"
    else:
        raise RuntimeError("the dressing sign is missing")
    piece = name.split("alpha*")[-1]
    if piece.endswith(")"):
        piece = piece[:-1]
    if not piece:
        raise RuntimeError("the dressing piece is empty")
    return sign, piece


def _shared_lines(figures: tuple[RecordedFigure, ...]) -> list[str]:
    groups: dict[str, list[tuple[str, str]]] = {}
    for figure in figures:
        sign, piece = _piece(figure.dressing.name)
        groups.setdefault(piece, []).append((figure.label, sign))
    lines = []
    for piece, rows in groups.items():
        if len(rows) < 2:
            continue
        labels = " and ".join(f"the {label}" for label, _sign in rows)
        signs = {sign for _label, sign in rows}
        if len(signs) == 1:
            lines.append(f"`{piece}` dresses {labels} with the same sign.")
        else:
            lines.append(f"`{piece}` dresses {labels} with opposite signs.")
    if len(lines) != 3:
        raise RuntimeError("the shared factor pieces changed")
    return lines


def _figure_lines(figure: RecordedFigure, helium: float) -> list[str]:
    lines = [
        f"The {figure.label} is `{_num(figure.target)}`.",
        f"Its recorded nearest is `{figure.anchor.name}` at `{_num(figure.anchor.value)}`.",
        f"The nearest dressing is `{figure.dressing.name}` at `{_num(figure.dressing.value)}`.",
        f"Its relative gap is `{_num(figure.dressing.relative_gap)}`.",
        (
            f"This anchor has `{figure.size.expressions}` dressings inside (0, 1), "
            f"`{figure.size.inside}` distinct, median spacing `{_num(figure.size.median_spacing)}`."
        ),
    ]
    if figure.dressing.relative_gap <= helium:
        lines.append("The dressing meets the helium quotient match.")
    else:
        lines.append("The dressing stays wider than the helium quotient match.")
    lines.append(_helium_sentence(figure.inside_helium))
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
    report = build_decay_recorded()
    met = sum(figure.dressing.relative_gap <= report.helium_mismatch for figure in report.figures)
    wider = len(report.figures) - met
    lines = [
        "# Fission decay recorded",
        "",
        "These six figures are the decay names the first dressing left untouched.",
        "The two limits already dressed stay on that score.",
        "The bare value is the recorded nearest already inside the review bar.",
        "The factor is the fine-structure leaf times one base piece, with either sign.",
        "A base piece is an atom, its square root, or a power from 2 through 6 in either sign.",
        f"alpha in the dressing is the fine-structure leaf `{_num(report.alpha)}`.",
        f"The helium quotient match is `{_num(report.helium_mismatch)}`.",
        f"{met} nearest dressings meet the helium quotient match.",
        f"{wider} nearest dressings stay wider than that match.",
        "",
        "## Figures",
        "",
    ]
    for figure in report.figures:
        lines.extend(_figure_lines(figure, report.helium_mismatch))
        lines.append("")
    lines.extend(_shared_lines(report.figures))
    lines.extend(
        [
            "",
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
