"""Tighten the decay-heat names with the depth-two seed grammar.

The watts stay the 1952 fit. The e-and-phi ladder is the depth-two grammar
already used for the deuteron. The five-atom copy uses the same operations
on pi, e, phi, gamma, and G. A gap inside one part in a thousand is recorded
for review. A gap inside the helium quotient match would stop the score. No
expression is installed, and none is multiplied into the watts.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from solcascade.energy_baseline.phase2.map_fsot import load_engine
from solcascade.production.fission_decay import build_decay
from solcascade.production.fission_observables import _e_phi_pieces, _helium_mismatch
from solcascade.production.fission_paths import _nearest_seeds
from solcascade.production.refine import project_root

_REVIEW_GAP = 1.0e-3
_ROUND = 12


@dataclass(frozen=True)
class NamedGap:
    name: str
    value: float
    relative_gap: float
    review: bool


@dataclass(frozen=True)
class FigureSeed:
    label: str
    target: float
    shallow: NamedGap
    ladder: NamedGap
    five: NamedGap


@dataclass(frozen=True)
class GrammarSize:
    expressions: int
    inside: int
    median_spacing: float


@dataclass(frozen=True)
class SeedReport:
    figures: tuple[FigureSeed, ...]
    ladder_size: GrammarSize
    five_size: GrammarSize
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
    return project_root() / "runs" / "production" / "fission_decay_seed.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _five_pieces(module) -> tuple[tuple[str, object], ...]:
    """Same operations as the e-and-phi ladder, on the five fission atoms."""
    atoms = (
        ("pi", module.PI),
        ("e", module.E),
        ("phi", module.PHI),
        ("gamma", module.GAMMA),
        ("G", module.G_CAT),
    )
    values: list[tuple[str, object]] = []
    for name, value in atoms:
        values.append((name, value))
        values.append((f"sqrt({name})", module.sqrt(value)))
        for power in range(2, 7):
            values.append((f"{name}^{power}", value**power))
            values.append((f"{name}^-{power}", value ** (-power)))
    pieces = list(values)
    for index, (left_name, left) in enumerate(values):
        for right_name, right in values[index:]:
            pieces.append((f"{left_name}+{right_name}", left + right))
            pieces.append((f"{left_name}*{right_name}", left * right))
            if right != 0:
                pieces.append((f"{left_name}/{right_name}", left / right))
            pieces.append((f"{left_name}-{right_name}", left - right))
    return tuple(pieces)


def _nearest(pieces, target: float) -> NamedGap:
    best: tuple[float, str, float] | None = None
    for name, value in pieces:
        number = float(value)
        if not math.isfinite(number) or not 0.0 < number < 1.0:
            continue
        gap = abs(number - target) / target
        if best is None or (gap, name) < (best[0], best[1]):
            best = (gap, name, number)
    if best is None:
        raise RuntimeError("the grammar has no seed inside (0, 1)")
    gap, name, number = best
    return NamedGap(name, number, gap, gap <= _REVIEW_GAP)


def _size(pieces) -> GrammarSize:
    vals = sorted({round(float(value), _ROUND) for _, value in pieces if 0.0 < float(value) < 1.0})
    spacings = [(right - left) / right for left, right in zip(vals, vals[1:]) if right > left]
    if not spacings:
        raise RuntimeError("the grammar has no interior spacing")
    spacings.sort()
    return GrammarSize(len(pieces), len(vals), spacings[len(spacings) // 2])


def _shallow(module, target: float) -> NamedGap:
    seed = _nearest_seeds(module, target)[0]
    return NamedGap(seed.name, seed.value, seed.relative_gap, seed.relative_gap <= _REVIEW_GAP)


def build_decay_seed() -> SeedReport:
    decay = build_decay()
    module, _digest = load_engine()
    ladder = _e_phi_pieces(module)
    five = _five_pieces(module)
    second, day, two = decay.points
    targets = (
        ("one-second fraction", second.year_fraction),
        ("one-second limit", second.limit_fraction),
        ("one-day fraction", day.year_fraction),
        ("one-day limit", day.limit_fraction),
        ("two-day fraction", two.year_fraction),
        ("two-day limit", two.limit_fraction),
        ("one-day drop", decay.day_over_second),
        ("two-day drop", decay.two_over_day),
    )
    helium = _helium_mismatch(module)
    figures: list[FigureSeed] = []
    for label, target in targets:
        shallow = _shallow(module, target)
        ladder_hit = _nearest(ladder, target)
        five_hit = _nearest(five, target)
        if not five_hit.relative_gap <= ladder_hit.relative_gap <= shallow.relative_gap:
            raise RuntimeError(f"the {label} grammar is out of order")
        if not five_hit.relative_gap < shallow.relative_gap:
            raise RuntimeError(f"the {label} did not tighten")
        for hit in (ladder_hit, five_hit):
            if hit.relative_gap <= helium:
                raise RuntimeError(f"{hit.name} meets the helium standard")
        figures.append(FigureSeed(label, target, shallow, ladder_hit, five_hit))
    return SeedReport(
        figures=tuple(figures),
        ladder_size=_size(ladder),
        five_size=_size(five),
        helium_mismatch=helium,
        second_watts=second.year_watts,
        day_watts=day.year_watts,
        two_watts=two.year_watts,
        nuclear_scalar=decay.nuclear_scalar,
        thermodynamics_scalar=decay.thermodynamics_scalar,
        gamma_fsot=decay.gamma_fsot,
        valve_beta=decay.valve_beta,
        pin_prefix=decay.pin_prefix,
        digest_prefix=decay.digest_prefix,
    )


def _figure_lines(figure: FigureSeed) -> list[str]:
    lines = [
        f"The {figure.label} is `{_num(figure.target)}`.",
        f"The shallow nearest is `{figure.shallow.name}` at `{_num(figure.shallow.value)}`, relative gap `{_num(figure.shallow.relative_gap)}`.",
        f"The e-and-phi ladder nearest is `{figure.ladder.name}` at `{_num(figure.ladder.value)}`, relative gap `{_num(figure.ladder.relative_gap)}`.",
        f"The five-atom nearest is `{figure.five.name}` at `{_num(figure.five.value)}`, relative gap `{_num(figure.five.relative_gap)}`.",
    ]
    if figure.ladder.review:
        lines.append("The e-and-phi expression meets the review bar.")
    else:
        lines.append("The e-and-phi expression stays outside the review bar.")
    if figure.five.review:
        lines.append("The five-atom expression meets the review bar.")
    else:
        lines.append("The five-atom expression stays outside the review bar.")
    lines.append("The expression is not installed.")
    return lines


def render() -> str:
    report = build_decay_seed()
    lines = [
        "# Fission decay seed",
        "",
        "The decay watts stay the 1952 fit.",
        "The shallow grammar is the one in the decay score.",
        "The e-and-phi ladder is the depth-two grammar from the deuteron check.",
        "The five-atom copy uses those same operations on pi, e, phi, gamma, and G.",
        f"The e-and-phi ladder has `{report.ladder_size.expressions}` expressions, `{report.ladder_size.inside}` inside (0, 1).",
        f"Its median relative spacing is `{_num(report.ladder_size.median_spacing)}`.",
        f"The five-atom copy has `{report.five_size.expressions}` expressions, `{report.five_size.inside}` inside (0, 1).",
        f"Its median relative spacing is `{_num(report.five_size.median_spacing)}`.",
        "A review-bar gap is one part in a thousand or closer.",
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
