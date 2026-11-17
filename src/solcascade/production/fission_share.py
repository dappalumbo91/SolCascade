"""Record the steam target inside the open interval under the Carnot.

The plants-being-built efficiency is that target. The steam score already
names it. The kept shares of the temperature Carnot stay reciprocals of the
named steps. The checked-in pulse is not rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from solcascade.production.fission_built import build_built
from solcascade.production.fission_ceiling import build_ceiling
from solcascade.production.fission_safety import BUILT_EFFICIENCY
from solcascade.production.fission_steam import build_steam
from solcascade.production.plants import FAST_ID, THERMAL_ID, build_plants
from solcascade.production.refine import catalog_record, project_root

_RATIONAL_GAP = 1.0e-15
_BUILT_SHARE = Fraction(504, 625)
_DEMONSTRATED_SHARE = Fraction(462, 625)
_CLOSED_SHARE = Fraction(42, 163)
_REMAINING_SHARE = Fraction(121, 163)
_UNUSED_BUILT = Fraction(121, 625)
_UNUSED_DEMONSTRATED = Fraction(163, 625)


@dataclass(frozen=True)
class ShareReport:
    steam_name: str
    steam_gap: float
    steam_second_gap: float
    demonstrated: float
    built: float
    carnot: float
    catalog: float
    hot_k: float
    cold_k: float
    built_step: float
    steam_carnot_step: float
    built_carnot_step: float
    kept_built: float
    kept_demonstrated: float
    closed_share: float
    remaining_share: float
    unused_built: float
    unused_demonstrated: float
    fast_efficiency: float
    alpha: float
    helium_mismatch: float
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_share.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _gap(derived: float, witness: float) -> float:
    if witness == 0.0:
        raise RuntimeError("the witness is zero")
    return abs(derived - witness) / abs(witness)


def _require_rational(value: float, rational: Fraction, label: str) -> None:
    if _gap(value, float(rational)) > _RATIONAL_GAP:
        raise RuntimeError(f"the {label} left {rational.numerator}/{rational.denominator}")


def build_share() -> ShareReport:
    plants = build_plants()
    ceiling = build_ceiling()
    steam = build_steam()
    built = build_built()
    if steam.pin_prefix != built.pin_prefix or steam.pin_prefix != plants.pin_prefix:
        raise RuntimeError("the share pin does not match the plant pin")
    if steam.digest_prefix != built.digest_prefix:
        raise RuntimeError("the share digest does not match the built digest")
    if ceiling.pin_prefix != plants.pin_prefix:
        raise RuntimeError("the ceiling pin does not match the plant pin")
    record = catalog_record(THERMAL_ID)
    if "k_eff is held at 1" not in str(record["governing_equation"]):
        raise RuntimeError("the governing equation does not hold k_eff at 1")
    fast_record = catalog_record(FAST_ID)
    if "Sodium" not in str(fast_record.get("engineering_barrier") or ""):
        raise RuntimeError("the fast barrier does not name sodium")
    fast_source = str(fast_record.get("efficiency_ceiling_source") or "")
    if "rejected near" in fast_source:
        raise RuntimeError("the fast ceiling source states a temperature pair")
    if ceiling.hot_k != 560.0 or ceiling.cold_k != 310.0:
        raise RuntimeError("the ceiling temperatures left 560 K and 310 K")
    if _gap(ceiling.carnot, 1.0 - ceiling.cold_k / ceiling.hot_k) > 0.0:
        raise RuntimeError("the temperature Carnot left one minus the rejected fraction")
    if BUILT_EFFICIENCY != 0.36 or built.built != BUILT_EFFICIENCY:
        raise RuntimeError("the plants-being-built figure left 0.36")
    if steam.built.target != BUILT_EFFICIENCY:
        raise RuntimeError("the steam target left the plants-being-built efficiency")
    if plants.fast.efficiency == BUILT_EFFICIENCY:
        raise RuntimeError("the fast efficiency collided with the plants-being-built figure")
    if plants.fast.efficiency != 0.4:
        raise RuntimeError("the fast efficiency left 0.4")
    demonstrated = built.demonstrated
    built_efficiency = built.built
    carnot = built.carnot
    if _gap(carnot, ceiling.carnot) > 0.0 or _gap(carnot, steam.carnot.target) > 1.0e-15:
        raise RuntimeError("the temperature Carnot left the ceiling score")
    if _gap(demonstrated, steam.demonstrated_efficiency) > 0.0:
        raise RuntimeError("the demonstrated efficiency left the steam score")
    if not demonstrated < built_efficiency < carnot < built.catalog:
        raise RuntimeError("the plants-being-built figure left the open interval under the Carnot")
    if built.catalog != 0.45 or steam.catalog_ceiling != built.catalog:
        raise RuntimeError("the catalog ceiling left 0.45")
    closing = steam.built.steps[-1]
    if len(steam.built.steps) != 3 or closing.inside_helium != 1:
        raise RuntimeError("the steam closing left the single helium dressing")
    if closing.dressing.relative_gap > steam.helium_mismatch or closing.second_gap <= steam.helium_mismatch:
        raise RuntimeError("the steam closing left the helium quotient match")
    if any(step.inside_helium != 0 for step in steam.built.steps[:-1]):
        raise RuntimeError("an earlier steam dressing already meets helium")
    if _num(closing.dressing.value) != _num(steam.built.target):
        raise RuntimeError("the steam closing left six-figure agreement")
    kept_built = built_efficiency / carnot
    kept_demonstrated = demonstrated / carnot
    if _gap(kept_built, 1.0 / built.carnot_step.target) > _RATIONAL_GAP:
        raise RuntimeError("the built share left the reciprocal of the built Carnot step")
    if _gap(kept_demonstrated, 1.0 / steam.carnot_step.target) > _RATIONAL_GAP:
        raise RuntimeError("the demonstrated share left the reciprocal of the steam Carnot step")
    closed_share = (built_efficiency - demonstrated) / (carnot - demonstrated)
    steam_span = steam.carnot_step.target - 1.0
    if _gap(closed_share, (built.built_step.target - 1.0) / steam_span) > _RATIONAL_GAP:
        raise RuntimeError("the closed share left the two named steps")
    remaining_share = (carnot - built_efficiency) / (carnot - demonstrated)
    if _gap(remaining_share, 1.0 - closed_share) > _RATIONAL_GAP:
        raise RuntimeError("the remaining share left the complement of the closed share")
    unused_built = (carnot - built_efficiency) / carnot
    unused_demonstrated = (carnot - demonstrated) / carnot
    if _gap(unused_built, 1.0 - kept_built) > _RATIONAL_GAP:
        raise RuntimeError("the unused built share left one minus the built share")
    if _gap(unused_demonstrated, 1.0 - kept_demonstrated) > _RATIONAL_GAP:
        raise RuntimeError("the unused demonstrated share left one minus the demonstrated share")
    if Fraction(9, 25) / Fraction(25, 56) != _BUILT_SHARE:
        raise RuntimeError("the built share left 504/625")
    if Fraction(33, 100) / Fraction(25, 56) != _DEMONSTRATED_SHARE:
        raise RuntimeError("the demonstrated share left 462/625")
    page_closed = (Fraction(9, 25) - Fraction(33, 100)) / (Fraction(25, 56) - Fraction(33, 100))
    if page_closed != _CLOSED_SHARE or Fraction(1) - _CLOSED_SHARE != _REMAINING_SHARE:
        raise RuntimeError("the closed share left 42/163")
    if Fraction(1) - _BUILT_SHARE != _UNUSED_BUILT:
        raise RuntimeError("the unused built share left 121/625")
    if Fraction(1) - _DEMONSTRATED_SHARE != _UNUSED_DEMONSTRATED:
        raise RuntimeError("the unused demonstrated share left 163/625")
    _require_rational(kept_built, _BUILT_SHARE, "built share")
    _require_rational(kept_demonstrated, _DEMONSTRATED_SHARE, "demonstrated share")
    _require_rational(closed_share, _CLOSED_SHARE, "closed share")
    _require_rational(remaining_share, _REMAINING_SHARE, "remaining share")
    _require_rational(unused_built, _UNUSED_BUILT, "unused built share")
    _require_rational(unused_demonstrated, _UNUSED_DEMONSTRATED, "unused demonstrated share")
    if _gap(built.built_step.target, built_efficiency / demonstrated) > _RATIONAL_GAP:
        raise RuntimeError("the built step left the plants-being-built efficiency over the demonstrated efficiency")
    return ShareReport(
        steam_name=closing.dressing.name,
        steam_gap=closing.dressing.relative_gap,
        steam_second_gap=closing.second_gap,
        demonstrated=demonstrated,
        built=built_efficiency,
        carnot=carnot,
        catalog=built.catalog,
        hot_k=ceiling.hot_k,
        cold_k=ceiling.cold_k,
        built_step=built.built_step.target,
        steam_carnot_step=steam.carnot_step.target,
        built_carnot_step=built.carnot_step.target,
        kept_built=kept_built,
        kept_demonstrated=kept_demonstrated,
        closed_share=closed_share,
        remaining_share=remaining_share,
        unused_built=unused_built,
        unused_demonstrated=unused_demonstrated,
        fast_efficiency=plants.fast.efficiency,
        alpha=steam.alpha,
        helium_mismatch=steam.helium_mismatch,
        nuclear_scalar=steam.nuclear_scalar,
        thermodynamics_scalar=steam.thermodynamics_scalar,
        gamma_fsot=steam.gamma_fsot,
        valve_beta=steam.valve_beta,
        pin_prefix=steam.pin_prefix,
        digest_prefix=steam.digest_prefix,
    )


def render() -> str:
    report = build_share()
    lines = [
        "# Fission share",
        "",
        "The target inside the open interval is the plants-being-built efficiency.",
        "The steam score names that efficiency.",
        "The chain stays at k_eff of 1.",
        "",
        "## Named target",
        "",
        f"The demonstrated efficiency is `{_num(report.demonstrated)}`.",
        f"The plants-being-built efficiency is `{_num(report.built)}`.",
        f"The page temperatures stay `{_num(report.hot_k)}` K and `{_num(report.cold_k)}` K.",
        f"The temperature Carnot is `{_num(report.carnot)}`.",
        f"The catalog ceiling is `{_num(report.catalog)}`.",
        "The demonstrated efficiency stays under the plants-being-built efficiency.",
        "The plants-being-built efficiency stays under the temperature Carnot.",
        "The temperature Carnot stays under the catalog ceiling.",
        f"The steam closing is `{report.steam_name}`.",
        f"The relative gap is `{_num(report.steam_gap)}`.",
        "The relative gap meets the helium quotient match.",
        "One dressing of that anchor meets the helium quotient match.",
        f"The next dressing sits at relative gap `{_num(report.steam_second_gap)}`.",
        "The closing value agrees with `0.36` at six figures.",
        "One expression names that efficiency.",
        "The worked figure stays the plants-being-built efficiency.",
        "",
        "## Kept shares",
        "",
        f"The plants-being-built efficiency over the temperature Carnot is `{_num(report.kept_built)}`.",
        "The same quotient is `504/625`.",
        f"The named Carnot step on the built score stays `{_num(report.built_carnot_step)}`.",
        "That quotient is the reciprocal of the named Carnot step on the built score.",
        "That quotient keeps the built score.",
        f"The demonstrated efficiency over the temperature Carnot is `{_num(report.kept_demonstrated)}`.",
        "The same quotient is `462/625`.",
        f"The named steam Carnot step stays `{_num(report.steam_carnot_step)}`.",
        "That quotient is the reciprocal of the named steam Carnot step.",
        "That quotient keeps the steam score.",
        f"The unused share at the plants-being-built efficiency is `{_num(report.unused_built)}`.",
        "The same quotient is `121/625`.",
        "The unused share is one minus the built share.",
        "The unused share keeps the built share.",
        f"The unused share at the demonstrated efficiency is `{_num(report.unused_demonstrated)}`.",
        "The same quotient is `163/625`.",
        "The unused share is one minus the demonstrated share.",
        "The unused share keeps the demonstrated share.",
        "",
        "## Closed interval",
        "",
        f"On the same heat, the electricity ratio is `{_num(report.built_step)}`.",
        "The electricity ratio is the named built step.",
        "The named built step stays the built score.",
        f"The closed share of the open interval is `{_num(report.closed_share)}`.",
        "The same quotient is `42/163`.",
        "The closed share is one less than the built step, divided by one less than the steam Carnot step.",
        "The closed share keeps those two steps.",
        f"The remaining share of the open interval is `{_num(report.remaining_share)}`.",
        "The same quotient is `121/163`.",
        "The remaining share is one minus the closed share.",
        "The remaining share keeps the closed share.",
        "The temperature Carnot stays the bound.",
        f"The fast efficiency stays `{_num(report.fast_efficiency)}`.",
        "The fast efficiency stays the bounds score.",
        "The fast ceiling source states no temperature pair.",
        "The fast row's filed barrier is the sodium system.",
        "The catalog ceiling stays unrun.",
        "The 60 GWd/t stack stays unrun.",
        "The fabrication cap stays 5 percent.",
        "This score adds no joules.",
        "The published stamp records these quotients as rationals of the named steps.",
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
        "This run does not install a seed as a new efficiency or a new plant row.",
        "",
        "The checked-in schedule is not rewritten.",
        "",
    ]
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
