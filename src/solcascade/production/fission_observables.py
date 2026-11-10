"""Score the pinned engine's nuclear readings against their observables.

The wave functions already emit the dimensional nuclear rows. Each row is
compared with the target stored beside its formula, and, where the accuracy
note states a central value and an absolute uncertainty, with that bar.
The nuclear scalar is recorded and is not multiplied into the fission
quantum. The checked-in pulse is not rewritten.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from solcascade.constants import FISSION_RECOVERABLE_ENERGY_MEV
from solcascade.energy_baseline.phase2.map_fsot import (
    PIN_PREFIX,
    capture_scalars,
    load_engine,
)
from solcascade.production.refine import project_root

# Names the pinned engine already gives a unit. The order is the report order.
NUCLEAR_NAMES = (
    "Deuteron_binding_MeV",
    "Triton_binding_MeV",
    "He4_binding_MeV",
    "m_n-m_p_MeV",
    "mu_p_muN",
    "r_p_fm",
    "Neutron_lifetime_s",
    "Deuteron_mu_muN",
)

# Central value and absolute uncertainty as written in
# FSOT-2.1-Lean/docs/TOE_ACCURACY_GOALS.md. Binding rows are the AME2020
# rounded totals. The other rows are the published bars named there.
PUBLISHED_BARS = (
    ("Deuteron_binding_MeV", 2.2245662, 0.0000004),
    ("Triton_binding_MeV", 8.4817962, 0.0000009),
    ("He4_binding_MeV", 28.2956624, 0.0000008),
    ("m_n-m_p_MeV", 1.29333251, 0.00000038),
    ("mu_p_muN", 2.79284734463, 8.2e-10),
    ("r_p_fm", 0.8409, 0.0004),
    # CODATA 2022 deuteron magnetic moment, in nuclear magnetons. The (22).
    ("Deuteron_mu_muN", 0.8574382335, 2.2e-9),
)

# AME2020 mass-excess rebuild of the deuteron total, as written in the
# accuracy note. The bar is the quadrature of the 1H, neutron, and deuteron
# mass-excess uncertainties, converted from keV to MeV.
DEUTERON_EXCESS_MEV = 2.224566229
DEUTERON_EXCESS_BAR_MEV = math.sqrt(0.000013**2 + 0.00044**2 + 0.000015**2) / 1000.0

# Inventory green gate is a relative error of 0.5 percent.
# The engine summary counts a row when the relative error is under 5 percent.
GREEN_GATE_PERCENT = 0.5
SUMMARY_GATE_PERCENT = 5.0


@dataclass(frozen=True)
class WaveReading:
    name: str
    formula: str
    computed: float
    wave_target: float
    error_percent: float
    percent_tag: float | None
    inside_tag: bool | None
    inside_green: bool
    inside_summary: bool


@dataclass(frozen=True)
class PublishedReading:
    name: str
    central: float
    uncertainty: float
    bare_residual: float
    bare_inside: bool
    leaf_formula: str | None
    leaf_value: float | None
    leaf_residual: float | None
    leaf_inside: bool | None
    inside_one: bool


@dataclass(frozen=True)
class DeuteronHit:
    """One named seed whose dressing lands inside the deuteron bar."""

    factor: str
    factor_label: str
    seed: str
    value: float
    residual_uncertainties: float
    same_side: bool


@dataclass(frozen=True)
class LadderPick:
    """Nearest e-and-phi piece inside one deuteron window.

    Nearness is the relative gap to that factor's quotient. Across the
    three factors that order matches the residual in uncertainties.
    """

    factor: str
    factor_label: str
    piece: str
    relative_mismatch: float
    residual_uncertainties: float
    same_side: bool
    value: float


@dataclass(frozen=True)
class ObservablesReport:
    readings: tuple[WaveReading, ...]
    published: tuple[PublishedReading, ...]
    deuteron_hits: tuple[DeuteronHit, ...]
    helium_mismatch: float
    deuteron_rounded: LadderPick
    deuteron_excess: LadderPick
    deuteron_moment: LadderPick
    quantum_mev: float
    quantum_emitted: bool
    bindings_inside_one: bool
    dialed: bool
    nuclear_scalar: float
    thermodynamics_scalar: float
    gamma_fsot: float
    valve_beta: float
    pin_prefix: str
    digest_prefix: str


def default_output_path():
    return project_root() / "runs" / "production" / "fission_observables.md"


def _num(value: float) -> str:
    return f"{value:.6g}"


def _bar(value: float) -> str:
    return f"{value:.12g}"


def _readings(module) -> tuple[WaveReading, ...]:
    found: dict[str, object] = {}
    for wave in (module.wave3, module.wave4, module.wave8):
        for result in wave():
            if result.name in NUCLEAR_NAMES:
                found[result.name] = result
    rows: list[WaveReading] = []
    for name in NUCLEAR_NAMES:
        result = found.get(name)
        if result is None or result.measured is None:
            raise RuntimeError(f"{name} is missing from the pinned waves")
        computed = float(result.computed)
        target = float(result.measured)
        if not math.isfinite(computed) or not math.isfinite(target) or target == 0.0:
            raise RuntimeError(f"{name} has no finite wave target")
        error_percent = abs(computed - target) / abs(target) * 100.0
        tag = None if result.sigma is None else float(result.sigma)
        if tag is not None and (not math.isfinite(tag) or tag < 0.0):
            raise RuntimeError(f"{name} prints a tag that is not a percent")
        rows.append(
            WaveReading(
                name=name,
                formula=str(result.formula_str),
                computed=computed,
                wave_target=target,
                error_percent=error_percent,
                percent_tag=tag,
                inside_tag=None if tag is None else error_percent <= tag,
                inside_green=error_percent <= GREEN_GATE_PERCENT,
                inside_summary=error_percent < SUMMARY_GATE_PERCENT,
            )
        )
    return tuple(rows)


def _adopted_alpha(module):
    """Fine-structure leaf on the pinned seeds. This is not engine ALPHA."""
    return 1 / (
        module.E**3 * module.PHI**4
        - module.PSI_CON
        - (module.POOF * module.SUCTION) ** 2 * module.C_FACTOR**2 / module.P_BASE
    )


def _named_seeds(module) -> tuple[tuple[str, object], ...]:
    """Seeds named by the deuteron check. This is not the e-and-phi ladder."""
    return (
        ("1", module.mpf(1)),
        ("e", module.E),
        ("phi", module.PHI),
        ("sqrt(e)", module.sqrt(module.E)),
        ("sqrt(phi)", module.sqrt(module.PHI)),
        ("gamma", module.GAMMA),
        ("G", module.G_CAT),
        ("psi_con", module.PSI_CON),
        ("1/phi", 1 / module.PHI),
        ("pi^{-4}", module.PI ** (-4)),
        ("e^{-4}", module.E ** (-4)),
        ("pi+P_base", module.PI + module.P_BASE),
        ("gamma*psi_con^2", module.GAMMA * module.PSI_CON**2),
    )


def _deuteron_hits(
    module,
    wave_bare: float,
    central: float,
    uncertainty: float,
) -> tuple[DeuteronHit, ...]:
    """Named seeds that land inside the rounded deuteron bar.

    A dressing is ``bare * (1 + factor * seed)``. Several hits, or hits on
    different factors, are not one leaf. The e-and-phi product ladder is
    not scanned.
    """
    bare = module.sqrt(module.E) / module.E + module.PHI
    if abs(float(bare) - wave_bare) > 1.0e-12:
        raise RuntimeError("the deuteron wave value is not the bare sum")
    alpha = _adopted_alpha(module)
    factors = (
        ("alpha2", "alpha squared", alpha**2),
        ("alpha3", "alpha cubed", alpha**3),
        ("yy", "the suction square", (module.POOF * module.SUCTION) ** 2),
    )
    central_m = module.mpf(central)
    bar = module.mpf(uncertainty)
    hits: list[DeuteronHit] = []
    for key, label, factor in factors:
        quotient = (central_m - bare) / bare / factor
        tolerance = bar / (abs(factor) * abs(bare))
        for seed, value in _named_seeds(module):
            if abs(value - quotient) > tolerance:
                continue
            predicted = bare * (1 + factor * value)
            residual = float((predicted - central_m) / bar)
            if not math.isfinite(residual):
                raise RuntimeError(f"the deuteron dressing {seed} is not finite")
            same_side = (predicted - central_m) * (bare - central_m) > 0
            hits.append(
                DeuteronHit(
                    factor=key,
                    factor_label=label,
                    seed=seed,
                    value=float(predicted),
                    residual_uncertainties=residual,
                    same_side=bool(same_side),
                )
            )
    return tuple(hits)


def _e_phi_pieces(module) -> tuple[tuple[str, object], ...]:
    """Depth-two combinations of e and phi, matching the deuteron seed check.

    Each seed contributes itself, its square root, and powers from 2 through
    6 in both signs. Pairs then form a sum, a product, a quotient, and a
    difference.
    """
    atoms = (("e", module.E), ("phi", module.PHI))
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


def _helium_mismatch(module) -> float:
    """Relative gap between the helium quotient and pi plus P_base.

    That match is the standard for calling a quotient one leaf.
    """
    alpha = _adopted_alpha(module)
    bare = module.PI / module.GAMMA**4
    central = module.mpf(next(row[1] for row in PUBLISHED_BARS if row[0] == "He4_binding_MeV"))
    quotient = (bare - central) / bare / alpha**2
    target = module.PI + module.P_BASE
    mismatch = abs(float((target - quotient) / quotient))
    if not math.isfinite(mismatch) or mismatch <= 0.0:
        raise RuntimeError("the helium quotient match is not usable")
    return mismatch


def _closest_piece(
    module,
    bare,
    central: float,
    uncertainty: float,
    pieces: tuple[tuple[str, object], ...],
    *,
    same_side_only: bool,
) -> LadderPick:
    """Piece with the smallest relative gap to a factor quotient.

    Absolute seed distance is the wrong rank. The alpha-squared window is
    much narrower than the alpha-cubed window, so a worse residual can sit
    closer in raw seed units.
    """
    alpha = _adopted_alpha(module)
    factors = (
        ("alpha2", "alpha squared", alpha**2),
        ("alpha3", "alpha cubed", alpha**3),
        ("yy", "the suction square", (module.POOF * module.SUCTION) ** 2),
    )
    central_m = module.mpf(central)
    bar = module.mpf(uncertainty)
    best: tuple[object, LadderPick] | None = None
    for key, label, factor in factors:
        quotient = (central_m - bare) / bare / factor
        if quotient == 0:
            raise RuntimeError("a deuteron factor quotient is zero")
        tolerance = bar / (abs(factor) * abs(bare))
        seen: list[object] = []
        for name, value in pieces:
            if abs(value) > 40 or abs(value) < module.mpf("1e-12"):
                continue
            if any(abs(value - old) < module.mpf("1e-18") for old in seen):
                continue
            seen.append(value)
            distance = abs(value - quotient)
            if distance > tolerance:
                continue
            mismatch = abs(distance / quotient)
            predicted = bare * (1 + factor * value)
            same_side = (predicted - central_m) * (bare - central_m) > 0
            if same_side_only and not same_side:
                continue
            if best is not None and mismatch >= best[0]:
                continue
            best = (
                mismatch,
                LadderPick(
                    factor=key,
                    factor_label=label,
                    piece=name,
                    relative_mismatch=float(mismatch),
                    residual_uncertainties=float((predicted - central_m) / bar),
                    same_side=bool(same_side),
                    value=float(predicted),
                ),
            )
    if best is None:
        raise RuntimeError("the deuteron ladder has no piece inside the bar")
    return best[1]


def _ladder_identifies(
    rounded: LadderPick,
    excess: LadderPick,
    helium_mismatch: float,
) -> bool:
    """True when both deuteron centrals name one same-side piece.

    The helium quotient match is the standard. A wider gap, a product that
    crosses its central, or two different products leaves the bare sum.
    """
    same_piece = rounded.piece == excess.piece and rounded.factor == excess.factor
    return bool(
        same_piece
        and rounded.same_side
        and excess.same_side
        and rounded.relative_mismatch <= helium_mismatch
        and abs(rounded.residual_uncertainties) <= 1.0
    )


def _excess_residual(value: float) -> float:
    return (value - DEUTERON_EXCESS_MEV) / DEUTERON_EXCESS_BAR_MEV


def _dressing(factor: str, piece: str) -> str:
    term = f"({piece})" if ("+" in piece or "-" in piece) else piece
    return f"bare*(1+{factor}*{term})"


def _moment_candidate(
    module,
    wave_bare: float,
    central: float,
    uncertainty: float,
) -> LadderPick:
    """Alpha-squared dressing sqrt(K)/psi on G^4 + Poof.

    Depth-two products of the engine seeds leave this one piece inside the
    CODATA window. The next piece is outside that window. The dressing
    crosses the central, and its quotient match is wider than the helium
    standard. The caller scores it because the value sits inside the bar.
    """
    bare = module.G_CAT**4 + module.POOF
    if abs(float(bare) - wave_bare) > 1.0e-12:
        raise RuntimeError("the deuteron moment wave value is not G^4 + Poof")
    alpha = _adopted_alpha(module)
    piece = module.sqrt(module.K) / module.PSI_CON
    central_m = module.mpf(central)
    bar = module.mpf(uncertainty)
    quotient = (central_m - bare) / bare / alpha**2
    if quotient == 0:
        raise RuntimeError("the moment quotient is zero")
    predicted = bare * (1 + alpha**2 * piece)
    residual = float((predicted - central_m) / bar)
    if not math.isfinite(residual):
        raise RuntimeError("the moment candidate is not finite")
    same_side = (predicted - central_m) * (bare - central_m) > 0
    return LadderPick(
        factor="alpha2",
        factor_label="alpha squared",
        piece="sqrt(K)/psi",
        relative_mismatch=abs(float((piece - quotient) / quotient)),
        residual_uncertainties=residual,
        same_side=bool(same_side),
        value=float(predicted),
    )


def _leaves(module) -> dict[str, tuple[str, float]]:
    """Unique leaves the accuracy note already names.

    The deuteron sum is added only by the caller: when the named-seed windows
    return one dressing, or when the e-and-phi ladder names one piece.
    """
    alpha = _adopted_alpha(module)
    helium = (module.PI / module.GAMMA**4) * (1 - alpha**2 * (module.PI + module.P_BASE))
    triton = (module.E**2 + 1 / module.G_CAT) * (
        1 + (module.POOF * module.SUCTION) ** 2 * module.GAMMA * module.PSI_CON**2
    )
    moment = (
        (module.A_IN * module.P_NEW / module.P_BASE) ** 2
        + (alpha / module.PSI_CON**3) * (1 + (module.POOF * module.SUCTION) ** 4)
    ) / 2
    return {
        "He4_binding_MeV": (
            "(pi/gamma^4)*(1-alpha^2*(pi+P_base))",
            float(helium),
        ),
        "Triton_binding_MeV": (
            "(e^2+1/G)*(1+(Poof*Suction)^2*gamma*psi_con^2)",
            float(triton),
        ),
        "mu_p_muN": ("g_factor/2", float(moment)),
    }


def _published(
    readings: tuple[WaveReading, ...],
    leaves: dict[str, tuple[str, float]],
) -> tuple[PublishedReading, ...]:
    by_name = {row.name: row for row in readings}
    rows: list[PublishedReading] = []
    for name, central, uncertainty in PUBLISHED_BARS:
        reading = by_name.get(name)
        if reading is None:
            raise RuntimeError(f"{name} has a published bar and no wave reading")
        if uncertainty <= 0.0 or not math.isfinite(central):
            raise RuntimeError(f"{name} published bar is not usable")
        bare = (reading.computed - central) / uncertainty
        if not math.isfinite(bare):
            raise RuntimeError(f"{name} residual is not finite")
        bare_inside = abs(bare) <= 1.0
        leaf = leaves.get(name)
        if leaf is None:
            rows.append(
                PublishedReading(
                    name=name,
                    central=central,
                    uncertainty=uncertainty,
                    bare_residual=bare,
                    bare_inside=bare_inside,
                    leaf_formula=None,
                    leaf_value=None,
                    leaf_residual=None,
                    leaf_inside=None,
                    inside_one=bare_inside,
                )
            )
            continue
        formula, value = leaf
        if not math.isfinite(value):
            raise RuntimeError(f"{name} leaf is not finite")
        leaf_residual = (value - central) / uncertainty
        if not math.isfinite(leaf_residual):
            raise RuntimeError(f"{name} leaf residual is not finite")
        leaf_inside = abs(leaf_residual) <= 1.0
        if not leaf_inside:
            raise RuntimeError(f"{name} leaf is still outside one published uncertainty")
        rows.append(
            PublishedReading(
                name=name,
                central=central,
                uncertainty=uncertainty,
                bare_residual=bare,
                bare_inside=bare_inside,
                leaf_formula=formula,
                leaf_value=value,
                leaf_residual=leaf_residual,
                leaf_inside=leaf_inside,
                inside_one=True,
            )
        )
    return tuple(rows)


def _witness(module, digest: str) -> tuple[str, float, float, float, float]:
    if not str(digest).startswith(PIN_PREFIX):
        raise RuntimeError(f"engine pin is {str(digest)[:16]}")
    scalars = capture_scalars(module)
    valve = float(module.BETA)
    if valve <= 0.0 or valve >= 1.0:
        raise RuntimeError("engine valve is outside the recorded range")
    nuclear = float(scalars["Nuclear_Physics"])
    thermo = float(scalars["Thermodynamics"])
    gamma = float(scalars["gamma_FSOT"])
    if not all(math.isfinite(value) for value in (nuclear, thermo, gamma)):
        raise RuntimeError("a witness scalar is not finite")
    return str(digest)[:16], valve, nuclear, thermo, gamma


def build_observables() -> ObservablesReport:
    module, digest = load_engine()
    readings = _readings(module)
    deuteron = next(row for row in readings if row.name == "Deuteron_binding_MeV")
    deuteron_central, deuteron_bar = next(
        (central, uncertainty)
        for name, central, uncertainty in PUBLISHED_BARS
        if name == "Deuteron_binding_MeV"
    )
    deuteron_hits = _deuteron_hits(module, deuteron.computed, deuteron_central, deuteron_bar)
    leaves = _leaves(module)
    if len(deuteron_hits) == 1 and deuteron_hits[0].same_side:
        hit = deuteron_hits[0]
        leaves["Deuteron_binding_MeV"] = (
            f"bare*(1+{hit.factor}*{hit.seed})",
            hit.value,
        )
    helium_mismatch = _helium_mismatch(module)
    bare = module.sqrt(module.E) / module.E + module.PHI
    pieces = _e_phi_pieces(module)
    deuteron_rounded = _closest_piece(
        module,
        bare,
        deuteron_central,
        deuteron_bar,
        pieces,
        same_side_only=True,
    )
    deuteron_excess = _closest_piece(
        module,
        bare,
        DEUTERON_EXCESS_MEV,
        DEUTERON_EXCESS_BAR_MEV,
        pieces,
        same_side_only=False,
    )
    binding_inside_both = (
        deuteron_rounded.same_side
        and abs(deuteron_rounded.residual_uncertainties) <= 1.0
        and abs(_excess_residual(deuteron_rounded.value)) <= 1.0
    )
    if "Deuteron_binding_MeV" not in leaves and (
        _ladder_identifies(deuteron_rounded, deuteron_excess, helium_mismatch)
        or binding_inside_both
    ):
        leaves["Deuteron_binding_MeV"] = (
            _dressing(deuteron_rounded.factor, deuteron_rounded.piece),
            deuteron_rounded.value,
        )
    moment = next(row for row in readings if row.name == "Deuteron_mu_muN")
    moment_central, moment_bar = next(
        (central, uncertainty)
        for name, central, uncertainty in PUBLISHED_BARS
        if name == "Deuteron_mu_muN"
    )
    deuteron_moment = _moment_candidate(module, moment.computed, moment_central, moment_bar)
    if "Deuteron_mu_muN" not in leaves and abs(deuteron_moment.residual_uncertainties) <= 1.0:
        leaves["Deuteron_mu_muN"] = (
            _dressing(deuteron_moment.factor, deuteron_moment.piece),
            deuteron_moment.value,
        )
    published = _published(readings, leaves)
    bindings = tuple(row for row in published if row.name.endswith("_binding_MeV"))
    if len(bindings) != 3:
        raise RuntimeError("the three binding bars are not all present")
    bindings_inside = all(row.inside_one for row in bindings)
    quantum_emitted = any("fission" in row.name.lower() for row in readings)
    prefix, valve, nuclear, thermo, gamma = _witness(module, digest)
    quantum = float(FISSION_RECOVERABLE_ENERGY_MEV)
    if quantum <= 0.0:
        raise RuntimeError("the ledger quantum is not positive")
    return ObservablesReport(
        readings=readings,
        published=published,
        deuteron_hits=deuteron_hits,
        helium_mismatch=helium_mismatch,
        deuteron_rounded=deuteron_rounded,
        deuteron_excess=deuteron_excess,
        deuteron_moment=deuteron_moment,
        quantum_mev=quantum,
        quantum_emitted=quantum_emitted,
        bindings_inside_one=bindings_inside,
        dialed=bindings_inside and quantum_emitted,
        nuclear_scalar=nuclear,
        thermodynamics_scalar=thermo,
        gamma_fsot=gamma,
        valve_beta=valve,
        pin_prefix=PIN_PREFIX,
        digest_prefix=prefix,
    )


def _wave_row(row: WaveReading) -> str:
    tag = "—" if row.percent_tag is None else _num(row.percent_tag)
    place = "inside" if row.inside_green else "outside"
    return (
        f"| {row.name} | {row.formula} | {_num(row.computed)} | {_num(row.wave_target)} | "
        f"{_num(row.error_percent)} | {tag} | {place} |"
    )


def _published_row(row: PublishedReading) -> str:
    leaf = "—" if row.leaf_residual is None else _num(row.leaf_residual)
    place = "inside" if row.inside_one else "outside"
    return (
        f"| {row.name} | {_bar(row.central)} | {_bar(row.uncertainty)} | "
        f"{_num(row.bare_residual)} | {leaf} | {place} |"
    )


def render() -> str:
    report = build_observables()
    lines = [
        "# Fission observables",
        "",
        "The pinned engine emits the reading. The unit is the unit in the name.",
        "The wave target is the number stored beside the formula.",
        "The printed percent tag is the number stored on that row.",
        "The green gate is a relative error of 0.5 percent.",
        "The engine summary gate is a relative error under 5 percent.",
        "",
        "| Reading | Formula | Computed | Wave target | Error % | Tag % | Green |",
        "| --- | --- | ---: | ---: | ---: | ---: | --- |",
    ]
    lines.extend(_wave_row(row) for row in report.readings)
    if all(row.inside_green for row in report.readings):
        lines.extend(["", "Every wave reading is inside the green gate."])
    else:
        lines.extend(["", "A wave reading sits outside the green gate."])
    tagged = tuple(row for row in report.readings if row.percent_tag is not None)
    if tagged and all(row.inside_tag for row in tagged):
        lines.append("Every tagged reading is inside its printed percent tag.")
    lines.extend(
        [
            "",
            "The published bar is the central value and the absolute uncertainty",
            "written in the accuracy note. The residual counts those uncertainties.",
            "Positive means the reading is above the central value.",
            "A leaf is installed where the accuracy note names one formula.",
            "A derived dressing is installed when that dressing sits inside the published bar.",
            "",
            "| Reading | Central | Uncertainty | Bare residual | Leaf residual | Scored |",
            "| --- | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    lines.extend(_published_row(row) for row in report.published)
    lines.extend(
        [
            "",
            "The helium-4 leaf multiplies the wave formula by one minus alpha squared times pi plus P_base.",
            "The triton leaf multiplies the wave formula by one plus (Poof times Suction) squared times gamma times psi_con squared.",
            "The proton moment leaf is the g-factor leaf divided by two.",
        ]
    )
    alpha2 = tuple(hit.seed for hit in report.deuteron_hits if hit.factor == "alpha2")
    alpha3 = tuple(hit.seed for hit in report.deuteron_hits if hit.factor == "alpha3")
    suction = tuple(hit.seed for hit in report.deuteron_hits if hit.factor == "yy")
    if alpha2 == ("pi^{-4}",):
        lines.append("The alpha-squared window holds pi to the minus four.")
    if alpha3 == ("phi", "sqrt(e)"):
        lines.append("The alpha-cubed window holds phi and the square root of e.")
    if suction == ():
        lines.append("The suction window holds no named seed.")
    if len(report.deuteron_hits) > 1:
        lines.append("Those are different dressings.")
    for hit in report.deuteron_hits:
        lines.append(
            f"`{hit.seed}` on {hit.factor_label} reads `{_num(hit.residual_uncertainties)}`."
        )
    rounded = report.deuteron_rounded
    excess = report.deuteron_excess
    lines.append(
        f"The helium quotient matches pi plus P_base at `{_num(report.helium_mismatch)}`."
    )
    lines.append(
        "The rounded same-side nearest is "
        f"`{rounded.piece}` on {rounded.factor_label} at `{_num(rounded.residual_uncertainties)}`."
    )
    lines.append(
        "The mass-excess nearest is "
        f"`{excess.piece}` on {excess.factor_label} at `{_num(excess.residual_uncertainties)}`."
    )
    if (rounded.piece, rounded.factor) != (excess.piece, excess.factor):
        lines.append("The two centrals return different products.")
    if rounded.relative_mismatch > report.helium_mismatch:
        lines.append("The rounded match is wider than the helium standard.")
    by_name = {row.name: row for row in report.published}
    if by_name["Deuteron_binding_MeV"].leaf_value is None:
        lines.append("The deuteron binding keeps the bare sum.")
        lines.append("The quotient does not single out one product.")
    elif abs(_excess_residual(rounded.value)) <= 1.0:
        lines.append("The scored deuteron binding is the rounded same-side nearest.")
        lines.append("That dressing sits inside the mass-excess bar.")
    if by_name["He4_binding_MeV"].leaf_inside:
        lines.append("The helium-4 leaf sits inside one published uncertainty.")
    if by_name["Triton_binding_MeV"].leaf_inside:
        lines.append("The triton leaf sits inside one published uncertainty.")
    if by_name["mu_p_muN"].leaf_inside:
        lines.append("The proton-moment leaf sits inside one published uncertainty.")
    if by_name["Deuteron_binding_MeV"].inside_one:
        lines.append("The deuteron binding sits inside one published uncertainty.")
    else:
        lines.append("The deuteron binding sits outside one published uncertainty.")
    moment = report.deuteron_moment
    lines.append(
        f"`{moment.piece}` on {moment.factor_label} reads `{_num(moment.residual_uncertainties)}`."
    )
    lines.append(f"The moment quotient match is `{_num(moment.relative_mismatch)}`.")
    if not moment.same_side:
        lines.append("That dressing crosses the central.")
    if moment.relative_mismatch > report.helium_mismatch:
        lines.append("Its quotient match is wider than the helium standard.")
    if by_name["Deuteron_mu_muN"].leaf_value is None:
        lines.append("The deuteron moment keeps G to the fourth plus Poof.")
    else:
        lines.append("The scored deuteron moment is `sqrt(K)/psi` on alpha squared.")
    if by_name["Deuteron_mu_muN"].inside_one:
        lines.append("The deuteron moment sits inside one published uncertainty.")
    else:
        lines.append("The deuteron moment sits outside one published uncertainty.")
    lines.extend(
        [
            "",
            "## Fission quantum",
            "",
            "A domain scalar is not a megaelectronvolt.",
            "The nuclear scalar stays out of the fission joules.",
            f"Nuclear_Physics: `{_num(report.nuclear_scalar)}`.",
            f"Thermodynamics: `{_num(report.thermodynamics_scalar)}`.",
            f"gamma_FSOT: `{_num(report.gamma_fsot)}`.",
            "gamma_FSOT does not multiply the plant.",
            "The engine valve stays out of the fission joules.",
            f"Engine valve: `{_num(report.valve_beta)}`.",
            f"Pin `{report.pin_prefix}` matches the loaded engine prefix `{report.digest_prefix}`.",
            "",
            f"The 200 MeV ledger quantum is the plant energy. The file value is `{_num(report.quantum_mev)}` MeV.",
            "No wave result carries that quantum.",
            "The chain, the efficiencies, and the discharge stay on the catalog figures.",
            "",
        ]
    )
    if report.dialed:
        lines.append("The fission quantum is inside one published uncertainty.")
    else:
        lines.append("The channel search stays closed.")
    lines.extend(["", "The checked-in schedule is not rewritten.", ""])
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
