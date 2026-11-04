"""Score FSOT scalars against the phase-one energy catalog.

The catalog files are read and not written. ``vendor/fsot_compute.py`` is
read and not written. A scalar is in range only when it already lies in the
interval that would multiply this row's coverage by ten and still respect
the physical cap. This module does not build a seed product to enter that
interval, and it does not write a mapping back into a mechanism file.
"""

from __future__ import annotations

import hashlib
import importlib.util
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import yaml

from solcascade.constants import warp_shell_energy_j
from solcascade.energy_baseline.rank import score_record

ENGINE_PATH = Path(r"C:\Users\damia\Desktop\FSOT-2.1-Lean\vendor\fsot_compute.py")
PIN_PREFIX = "aeb2adad6e80f487"
ORDER_OF_MAGNITUDE = 10.0

# Domains named by the catalog notes, plus the stable outputs the prompt lists.
CAPTURED_DOMAINS = (
    "Particle_Physics",
    "High_Energy_Physics",
    "Electromagnetism",
    "Nuclear_Physics",
    "Thermodynamics",
    "Condensed_Matter",
    "Chemistry",
    "Quantum_Mechanics",
    "Optics",
    "Acoustics",
    "Fluid_Dynamics",
    "Atmospheric_Physics",
    "Oceanography",
    "Geophysics",
)
EXTRA_OUTPUTS = ("K", "phi^-4", "gamma_FSOT")

_DOMAIN_IN_NOTE = re.compile(r"domain_scalar\(([A-Za-z0-9_]+)\)")
_ENGINE: tuple[object, str] | None = None


@dataclass(frozen=True)
class ParameterReading:
    """Which free parameter dominates one mechanism's coverage.

    ``use_record_efficiency`` reads the demonstrated value from the YAML.
    ``use_record_ceiling`` reads the physical cap from the YAML. An explicit
    ``physical_max`` is used when the ceiling field answers a different
    question than the scored reference unit. ``candidate`` overrides the
    domain named in the phase-one note.
    """

    parameter: str
    kind: str
    use_record_efficiency: bool = False
    use_record_ceiling: bool = False
    demonstrated: float | None = None
    physical_max: float | None = None
    candidate: str | None = None
    aside: str | None = None


def _linear(parameter: str, *, candidate: str | None = None, aside: str | None = None, physical_max: float | None = None) -> ParameterReading:
    return ParameterReading(
        parameter=parameter,
        kind="dimensionless_linear",
        use_record_efficiency=True,
        use_record_ceiling=physical_max is None,
        physical_max=physical_max,
        candidate=candidate,
        aside=aside,
    )


def _dim(parameter: str, *, candidate: str | None = None, aside: str | None = None) -> ParameterReading:
    return ParameterReading(
        parameter=parameter,
        kind="dimensionful",
        candidate=candidate,
        aside=aside,
    )


def _open(parameter: str, *, aside: str, candidate: str | None = None) -> ParameterReading:
    return ParameterReading(
        parameter=parameter,
        kind="unscored",
        candidate=candidate,
        aside=aside,
    )


# Dominant parameter for each phase-one id. The YAML remains the source of
# the joules and of any efficiency that multiplies them.
READINGS: dict[str, ParameterReading] = {
    "fission-thermal": ParameterReading(
        parameter="burnup fraction",
        kind="dimensionless_linear",
        demonstrated=1.0,
        physical_max=1.0,
        candidate="Nuclear_Physics",
        aside="The reference unit is already every nucleus fissioning. Rankine efficiency in the file runs from 0.33 to a 0.45 ceiling and cannot supply a factor of ten either.",
    ),
    "fission-fast": ParameterReading(
        parameter="burnup fraction",
        kind="dimensionless_linear",
        demonstrated=1.0,
        physical_max=1.0,
        candidate="Nuclear_Physics",
        aside="The reference unit is already every nucleus fissioning.",
    ),
    "fusion-dt-magnetic": _dim(
        "confinement time",
        candidate="Nuclear_Physics",
        aside="The scored joules are one JET yield. A domain scalar is not a confinement time and is not applied to the ledger clock.",
    ),
    "fusion-dt-inertial": _dim(
        "confinement time",
        candidate="Nuclear_Physics",
        aside="The scored joules are one NIF yield. The nuclear scalar is not a pulse width.",
    ),
    "fusion-dt-pulsed-power": _dim(
        "confinement time",
        candidate="Nuclear_Physics",
        aside="The catalog leaves the shot yield blank.",
    ),
    "fusion-d-he3": _dim(
        "confinement time",
        candidate="Nuclear_Physics",
        aside="The file holds the reaction energy of a complete burn and no demonstrated device efficiency.",
    ),
    "fusion-aneutronic-pb11": _dim(
        "confinement time",
        candidate="Nuclear_Physics",
        aside="The file holds the reaction energy of a complete burn and no demonstrated device efficiency.",
    ),
    "chemical-combustion": _linear("plant efficiency"),
    "chemical-battery": _linear("round-trip efficiency"),
    "chemical-fuel-cell": _linear("stack efficiency"),
    "photovoltaic": _dim(
        "bandgap",
        aside="Fill factor sits inside the scored 0.267. A factor of ten on that efficiency would pass both 1 and the 0.33 single-junction ceiling in the file. Bandgap is in electron-volts.",
    ),
    "thermoelectric": _linear("device efficiency"),
    "thermionic": _linear("converter efficiency"),
    "mhd-generator": _linear(
        "channel efficiency",
        candidate="Electromagnetism",
        aside="The phase-one note names Fluid_Dynamics. The Faraday emf is electromagnetic, so the candidate is that scalar.",
    ),
    "betavoltaic": _linear(
        "conversion efficiency",
        aside="The catalog did not score an isotope inventory, so coverage is blank.",
    ),
    "rtg": _linear("converter efficiency"),
    "hydroelectric": _dim(
        "hydraulic head",
        aside="Turbine efficiency is 0.90 against a ceiling of 1, so the head is the lever that can still move the joules. The head is in meters.",
    ),
    "wind": _dim(
        "wind speed",
        aside="Power scales as speed cubed. The Betz cap on the power coefficient is 16/27 against a scored 0.45, which is also short of a factor of ten.",
    ),
    "geothermal": _linear("cycle efficiency"),
    "tidal": _dim(
        "tidal range",
        aside="Basin energy scales as range squared. Turbine efficiency is 0.80 against a ceiling of 1. The range is in meters.",
    ),
    "wave": _linear("absorption efficiency"),
    "biomass": _linear(
        "boiler efficiency",
        physical_max=1.0,
        aside="The YAML ceiling 0.06 is the photosynthetic bound. It does not cap joules per kilogram of wood, so the cap used here is 1.",
    ),
    "osmotic": _linear(
        "membrane efficiency",
        aside="Demonstrated electrical efficiency is blank, so the mixing resource is not a coverage.",
    ),
    "electrostatic-generator": _dim(
        "breakdown field",
        candidate="Electromagnetism",
        aside="The file already scores air at 3 MV/m. The 0.50 efficiency allowance can only double.",
    ),
    "casimir-static": _dim(
        "plate separation",
        aside="Cavity energy scales as the inverse cube of the gap. The ranked efficiency is 0 because no net-work cycle is demonstrated.",
    ),
    "casimir-dynamical": _dim(
        "boundary-drive frequency",
        aside="The experiment reports photons, not a joule inventory.",
    ),
    "quantum-vacuum-claims": _open(
        "no device parameter",
        aside="The row is a claim class with no accepted device and no joules.",
    ),
    "lenr-cold-fusion": _open(
        "excess-heat fraction",
        aside="No joule figure was copied, and the replication status is disputed.",
    ),
    "sonoluminescence": _dim(
        "acoustic drive",
        aside="The flash is a measured optical result. No electrical joules are scored.",
    ),
    "triboelectric": _linear("conversion efficiency"),
    "piezoelectric": _linear("conversion efficiency"),
    "atmospheric-electricity": _linear("collection efficiency"),
    "mhd-propulsion": _open(
        "none",
        candidate="Electromagnetism",
        aside="The row is a propulsion load. It scores no delivered joules.",
    ),
    "thermoelectric-antarctica-2025": _linear(
        "device efficiency",
        physical_max=1.0,
        aside="The scored 38 kWh is already electrical, so the file efficiency is 1. The 0.20 ceiling is the general thermoelectric bound and is not a second factor on those joules.",
    ),
    "radiative-cooling-stirling-2025": _dim(
        "radiator temperature difference",
        aside="The scored 0.40 W/m2 is already mechanical output. The paper's higher potential is not substituted.",
    ),
    "earth-rotation-generator-2025": _open(
        "claimed voltage",
        candidate="Electromagnetism",
        aside="Joules were left blank. The reported effect is disputed.",
    ),
    "thermoelectric-ionogel-2026": _open(
        "device efficiency",
        aside="The paper reports a power density. Electrical efficiency is blank, so there is no coverage to move.",
    ),
    "quantum-battery-2026": _open(
        "no device parameter",
        aside="The row is a theory paper. No device joules were measured.",
    ),
}


def load_engine():
    """Load the pinned engine the way ``scripts/apply_fsot.py`` loads it."""
    global _ENGINE
    if _ENGINE is not None:
        return _ENGINE
    digest = hashlib.sha256(ENGINE_PATH.read_bytes()).hexdigest()
    if not digest.startswith(PIN_PREFIX):
        raise SystemExit(f"engine pin is {digest[:16]}, expected {PIN_PREFIX}")
    spec = importlib.util.spec_from_file_location("fsot_compute_live", ENGINE_PATH)
    if spec is None or spec.loader is None:
        raise SystemExit(f"cannot load {ENGINE_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    _ENGINE = (module, digest)
    return _ENGINE


def capture_scalars(module) -> dict[str, float]:
    """Dimensionless outputs only. No rescaling into keV, tesla, or seconds."""
    scalars = {name: float(module.domain_scalar(name)) for name in CAPTURED_DOMAINS}
    condensed = scalars["Condensed_Matter"]
    k_value = float(module.K)
    scalars["K"] = k_value
    scalars["phi^-4"] = float(module.PHI ** (-4))
    # gamma_FSOT is the orifice witness S/K already used on the coil ledger.
    # It is not the engine's Euler-Mascheroni constant.
    scalars["gamma_FSOT"] = condensed / k_value
    return scalars


def candidate_name(record: dict, reading: ParameterReading) -> str:
    if reading.candidate:
        return reading.candidate
    match = _DOMAIN_IN_NOTE.search(str(record.get("fsot_note") or ""))
    if match:
        return match.group(1)
    return "unspecified"


def _with_aside(text: str, aside: str | None) -> str:
    if aside:
        return f"{text} {aside}"
    return text


def compare_record(record: dict, scalars: dict[str, float], target_j: float) -> dict:
    """One mechanism against one captured scalar. Missing numbers do not raise."""
    record = dict(record or {})
    reading = READINGS.get(record.get("id"))
    if reading is None:
        reading = _open(
            "unspecified",
            aside="No dominant-parameter reading is recorded for this mechanism.",
        )
    scored = score_record(record, target_j)
    coverage = scored["output_over_target"]
    name = candidate_name(record, reading)
    scalar = scalars.get(name)
    row = {
        "id": record.get("id"),
        "parameter": reading.parameter,
        "kind": reading.kind,
        "candidate": name,
        "scalar_value": scalar,
        "in_range": False,
        "coverage": coverage,
        "coverage_if_fixed": None,
        "upside": None,
        "phase1_status": record.get("fsot_status") or "open_gap",
        "gap_note": "",
    }
    if reading.kind == "dimensionful":
        row["gap_note"] = _with_aside(
            f"The dominant parameter is {reading.parameter}, which is dimensionful. "
            f"{name} is dimensionless, and the engine has no zero-parameter conversion into that unit.",
            reading.aside,
        )
        return row
    if reading.kind == "unscored":
        row["gap_note"] = _with_aside(
            "No scored coverage can be moved.",
            reading.aside,
        )
        return row
    if coverage is None or coverage <= 0.0:
        row["gap_note"] = _with_aside(
            "The catalog row has no positive delivered joules, so a coverage change is undefined.",
            reading.aside,
        )
        return row
    demonstrated = reading.demonstrated
    if reading.use_record_efficiency:
        demonstrated = scored["demonstrated_efficiency"]
    physical_max = reading.physical_max
    if reading.use_record_ceiling:
        ceiling = record.get("efficiency_ceiling")
        physical_max = None if ceiling is None else float(ceiling)
    if demonstrated is None or demonstrated <= 0.0 or physical_max is None:
        row["gap_note"] = _with_aside(
            f"The demonstrated {reading.parameter} is missing or not positive, so a ratio to the scalar is undefined.",
            reading.aside,
        )
        return row
    if physical_max >= demonstrated:
        row["upside"] = coverage * (physical_max / demonstrated) - coverage
    else:
        row["upside"] = 0.0
    if scalar is None:
        row["gap_note"] = _with_aside(
            f"The captured outputs have no value named {name}.",
            reading.aside,
        )
        return row
    needed_low = ORDER_OF_MAGNITUDE * demonstrated
    if needed_low > physical_max:
        row["gap_note"] = _with_aside(
            f"A tenfold coverage gain needs {reading.parameter} >= {needed_low:.6g}, "
            f"and a real device is capped at {physical_max:.6g}.",
            reading.aside,
        )
        return row
    if not (needed_low <= float(scalar) <= physical_max):
        row["gap_note"] = _with_aside(
            f"{name} at {float(scalar):.6g} sits outside the window "
            f"[{needed_low:.6g}, {physical_max:.6g}] that would raise coverage by one order.",
            reading.aside,
        )
        return row
    fixed = coverage * (float(scalar) / demonstrated)
    orders_short = -math.log10(fixed) if fixed > 0.0 else None
    row["in_range"] = True
    row["coverage_if_fixed"] = fixed
    orders_text = f"{orders_short:.6g}" if orders_short is not None else ""
    row["gap_note"] = _with_aside(
        f"In range for a one-order move. Fixing {reading.parameter} at {name} "
        f"{float(scalar):.6g} sets coverage to {fixed:.6g}. "
        f"The shell gap that remains is {orders_text} orders.",
        reading.aside,
    )
    return row


def load_mechanism_records(directory: Path) -> list[dict]:
    records = []
    for path in sorted(directory.glob("*.yaml")):
        try:
            loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
        except yaml.YAMLError as exc:
            records.append({"id": path.stem, "_load_error": str(exc)})
            continue
        if not isinstance(loaded, dict):
            records.append({"id": path.stem, "_load_error": "file is not a mapping"})
            continue
        records.append(loaded)
    return records


def map_directory(directory: Path, scalars: dict[str, float], target_j: float) -> list[dict]:
    rows = []
    for record in load_mechanism_records(directory):
        if record.get("_load_error"):
            rows.append(
                {
                    "id": record.get("id"),
                    "parameter": "unspecified",
                    "kind": "unscored",
                    "candidate": "unspecified",
                    "scalar_value": None,
                    "in_range": False,
                    "coverage": None,
                    "coverage_if_fixed": None,
                    "upside": None,
                    "phase1_status": "open_gap",
                    "gap_note": f"The mechanism file did not load: {record['_load_error']}",
                }
            )
            continue
        rows.append(compare_record(record, scalars, target_j))
    return rows


def result_sentence(rows: list[dict]) -> str:
    hits = [row for row in rows if row["in_range"]]
    if not hits:
        return (
            "No captured FSOT scalar lands in a window that moves any row "
            "by one order of magnitude toward the shell."
        )
    largest = max(hits, key=lambda row: row["coverage_if_fixed"] / row["coverage"])
    factor = largest["coverage_if_fixed"] / largest["coverage"]
    return (
        f"{len(hits)} row(s) land in a one-order window. "
        f"The largest factor is {factor:.6g} on {largest['id']}, "
        f"from coverage {largest['coverage']:.6g} to {largest['coverage_if_fixed']:.6g}. "
        "That factor is the arithmetic of fixing the named parameter at the scalar. "
        "It is not a closure of the shell gap."
    )


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value).replace("|", "/").replace("\n", " ")


def render_markdown(rows: list[dict], target_j: float, scalars: dict[str, float], pin: str | None) -> str:
    lines = [
        "# Energy baseline phase two",
        "",
        "Observational comparison. The phase-one catalog is not edited.",
        "A row is in range when the candidate scalar is at least ten times",
        "the demonstrated parameter and still at or under the physical cap.",
        "Coverage if fixed is the catalog coverage times scalar / demonstrated.",
        "No seed product is formed.",
        "",
        result_sentence(rows),
        "",
        f"Warp-shell target: `{target_j:.6g}` J, from `warp_shell_energy_j()`.",
        f"Engine pin: `{pin}`." if pin else "Engine pin: scalars were supplied by the caller.",
        "",
        "Captured dimensionless outputs:",
        "",
    ]
    for name in list(CAPTURED_DOMAINS) + list(EXTRA_OUTPUTS):
        if name in scalars:
            lines.append(f"- `{name}`: `{scalars[name]:.6g}`")
    lines.extend(
        [
            "",
            "| Mechanism id | Dominant free parameter | FSOT scalar candidate | Scalar value | In range | Coverage if fixed | Gap note |",
            "| --- | --- | --- | ---: | --- | ---: | --- |",
        ]
    )
    for row in rows:
        lines.append(
            "| {id} | {parameter} | {candidate} | {scalar} | {flag} | {fixed} | {note} |".format(
                id=_cell(row.get("id")),
                parameter=_cell(row.get("parameter")),
                candidate=_cell(row.get("candidate")),
                scalar=_cell(row.get("scalar_value")),
                flag=_cell(row.get("in_range")),
                fixed=_cell(row.get("coverage_if_fixed")),
                note=_cell(row.get("gap_note")),
            )
        )
    lines.append("")
    return "\n".join(lines)


def default_paths() -> tuple[Path, Path, Path]:
    package = Path(__file__).resolve().parent
    catalog = package.parent / "data"
    # phase2 -> energy_baseline -> solcascade -> src -> project root
    root = package.parents[3]
    output = root / "runs" / "energy_baseline" / "phase2_mapping.md"
    return catalog, output, root / "runs" / "energy_baseline" / "gaps.md"


def main(argv: list[str] | None = None) -> int:
    del argv
    catalog, output, _gaps = default_paths()
    _module, digest = load_engine()
    scalars = capture_scalars(_module)
    target = float(warp_shell_energy_j())
    rows = map_directory(catalog, scalars, target)
    text = render_markdown(rows, target, scalars, digest[:16])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text, encoding="utf-8")
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
