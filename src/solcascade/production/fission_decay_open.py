"""Close the decay names that stayed wider than the helium match.

The bare value is the recorded dressing. Each added factor is the
fine-structure leaf times one base piece, with either sign. The spacing
stays the base-piece spacing. A unique closing match is recorded. The
1952 fit stays the cooling load.
"""

from __future__ import annotations

from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay_dress import (
    _base_pieces,
    _dressings,
    build_decay_dress,
)
from solcascade.production.fission_decay_recorded import build_decay_recorded
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

_OPEN = (
    "one-second limit",
    "one-day fraction",
    "one-day limit",
    "two-day fraction",
    "one-day drop",
)
_ADDED = {
    "one-second limit": 2,
    "one-day fraction": 1,
    "one-day limit": 2,
    "two-day fraction": 1,
    "one-day drop": 1,
}
_MAX_ADDED = 2
_SPACING = "0.000511364"
_COUNT = 116
_TENTH = 1.0e-4


@dataclass(frozen=True)
class OpenStep:
    dressing: NamedGap
    inside_review: int
    inside_tenth: int
    inside_helium: int
    size: GrammarSize


@dataclass(frozen=True)
class OpenFigure:
    label: str
    target: float
    start_name: str
    start_value: float
    start_gap: float
    steps: tuple[OpenStep, ...]


@dataclass(frozen=True)
class OpenReport:
    figures: tuple[OpenFigure, ...]
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
    return project_root() / "runs" / "production" / "fission_decay_open.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _current(label: str, recorded, dress):
    if label in dress:
        figure = dress[label]
    else:
        figure = recorded[label]
    return figure.target, figure.dressing


def build_decay_open() -> OpenReport:
    seed = build_decay_seed()
    recorded = {figure.label: figure for figure in build_decay_recorded().figures}
    dress = {figure.label: figure for figure in build_decay_dress().figures}
    module, _digest = load_engine()
    alpha = float(_adopted_alpha(module))
    helium = _helium_mismatch(module)
    base = _base_pieces(module)
    labels = []
    for figure in seed.figures:
        _target, dressing = _current(figure.label, recorded, dress)
        if dressing.relative_gap > helium:
            labels.append(figure.label)
    if tuple(labels) != _OPEN:
        raise RuntimeError("the open set changed")
    figures: list[OpenFigure] = []
    for label in labels:
        target, dressing = _current(label, recorded, dress)
        name = dressing.name
        value = dressing.value
        gap = dressing.relative_gap
        steps: list[OpenStep] = []
        for _added in range(_MAX_ADDED):
            rows = _dressings(name, value, alpha, base)
            hit = _nearest(rows, target)
            if hit.relative_gap >= gap:
                raise RuntimeError(f"the {label} dressing did not tighten")
            size = _size(rows)
            if size.expressions != _COUNT or size.inside != _COUNT:
                raise RuntimeError(f"the {label} dressing count changed")
            if _num(size.median_spacing) != _SPACING:
                raise RuntimeError(f"the {label} spacing changed")
            if size.median_spacing <= helium:
                raise RuntimeError(f"the {label} spacing is finer than the helium match")
            inside_review = sum(
                abs(number - target) / target <= _REVIEW_GAP for _, number in rows
            )
            inside_tenth = sum(abs(number - target) / target <= _TENTH for _, number in rows)
            inside_helium = sum(abs(number - target) / target <= helium for _, number in rows)
            if inside_helium > 1:
                raise RuntimeError(f"several {label} dressings meet the helium standard")
            steps.append(
                OpenStep(
                    dressing=hit,
                    inside_review=inside_review,
                    inside_tenth=inside_tenth,
                    inside_helium=inside_helium,
                    size=size,
                )
            )
            name = hit.name
            value = hit.value
            gap = hit.relative_gap
            if gap <= helium:
                break
        if len(steps) != _ADDED[label]:
            raise RuntimeError(f"the {label} added-factor count changed")
        closing = steps[-1]
        if closing.dressing.relative_gap > helium or closing.inside_helium != 1:
            raise RuntimeError(f"the {label} closing match is not unique")
        if any(step.inside_helium != 0 for step in steps[:-1]):
            raise RuntimeError(f"an earlier {label} dressing already meets helium")
        figures.append(
            OpenFigure(
                label=label,
                target=target,
                start_name=dressing.name,
                start_value=dressing.value,
                start_gap=dressing.relative_gap,
                steps=tuple(steps),
            )
        )
    ledger = build_decay_dress()
    return OpenReport(
        figures=tuple(figures),
        alpha=alpha,
        helium_mismatch=helium,
        second_watts=ledger.second_watts,
        day_watts=ledger.day_watts,
        two_watts=ledger.two_watts,
        nuclear_scalar=ledger.nuclear_scalar,
        thermodynamics_scalar=ledger.thermodynamics_scalar,
        gamma_fsot=ledger.gamma_fsot,
        valve_beta=ledger.valve_beta,
        pin_prefix=ledger.pin_prefix,
        digest_prefix=ledger.digest_prefix,
    )


def _count_sentence(count: int, singular: str, plural: str) -> str:
    if count == 1:
        return f"One dressing of the closing anchor {singular}."
    return f"{count} dressings of the closing anchor {plural}."


def _closing_piece(name: str) -> tuple[str, str]:
    body = name[:-1] if name.endswith(")") else name
    minus = "*(1-alpha*"
    plus = "*(1+alpha*"
    if minus in body:
        return "minus", body[body.rfind(minus) + len(minus) :]
    if plus in body:
        return "plus", body[body.rfind(plus) + len(plus) :]
    raise RuntimeError("the closing sign is missing")


def _closing_lines(figures: tuple[OpenFigure, ...]) -> list[str]:
    groups: dict[tuple[str, str], list[str]] = {}
    for figure in figures:
        sign, piece = _closing_piece(figure.steps[-1].dressing.name)
        if not piece:
            raise RuntimeError("the closing piece is empty")
        groups.setdefault((piece, sign), []).append(figure.label)
    lines = []
    for (piece, sign), labels in groups.items():
        if len(labels) == 1:
            text = f"the {labels[0]}"
        elif len(labels) == 2:
            text = f"the {labels[0]} and the {labels[1]}"
        else:
            head = ", ".join(f"the {label}" for label in labels[:-1])
            text = f"{head}, and the {labels[-1]}"
        word = "a minus sign" if sign == "minus" else "a plus sign"
        lines.append(f"`{piece}` closes {text} with {word}.")
    if len(lines) != 2:
        raise RuntimeError("the closing pieces changed")
    return lines


def _figure_lines(figure: OpenFigure) -> list[str]:
    lines = [
        f"The {figure.label} is `{_num(figure.target)}`.",
        (
            f"Its recorded dressing is `{figure.start_name}` at `{_num(figure.start_value)}`, "
            f"relative gap `{_num(figure.start_gap)}`."
        ),
    ]
    early = figure.steps[:-1]
    closing = figure.steps[-1]
    for step in early:
        lines.append(
            f"The added dressing is `{step.dressing.name}` at `{_num(step.dressing.value)}`, "
            f"relative gap `{_num(step.dressing.relative_gap)}`."
        )
        lines.append("The added dressing stays wider than the helium quotient match.")
    lines.extend(
        [
            (
                f"The closing dressing is `{closing.dressing.name}` at "
                f"`{_num(closing.dressing.value)}`, relative gap "
                f"`{_num(closing.dressing.relative_gap)}`."
            ),
            "The closing dressing meets the helium quotient match.",
            "One dressing of the closing anchor meets the helium quotient match.",
            _count_sentence(
                closing.inside_review,
                "meets the review bar",
                "meet the review bar",
            ),
            _count_sentence(
                closing.inside_tenth,
                "is inside one part in ten thousand",
                "are inside one part in ten thousand",
            ),
            (
                f"The closing anchor has `{closing.size.expressions}` dressings inside (0, 1), "
                f"`{closing.size.inside}` distinct, median spacing "
                f"`{_num(closing.size.median_spacing)}`."
            ),
            "The dressing is not installed.",
        ]
    )
    return lines


def render() -> str:
    report = build_decay_open()
    lines = [
        "# Fission decay open",
        "",
        "These five figures were still wider than the helium quotient match.",
        "The bare value is the recorded dressing.",
        "Each added factor is the fine-structure leaf times one base piece, with either sign.",
        "A base piece is an atom, its square root, or a power from 2 through 6 in either sign.",
        f"alpha in the dressing is the fine-structure leaf `{_num(report.alpha)}`.",
        f"The helium quotient match is `{_num(report.helium_mismatch)}`.",
        f"Every added anchor has `{_COUNT}` dressings inside (0, 1), median spacing `{_SPACING}`.",
        "That spacing stays wider than the helium quotient match.",
        "Every closing gap meets the helium quotient match.",
        "",
        "## Figures",
        "",
    ]
    for figure in report.figures:
        lines.extend(_figure_lines(figure))
        lines.append("")
    lines.extend(_closing_lines(report.figures))
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
