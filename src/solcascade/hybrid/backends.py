"""Input writers for Serpent 2, MCNP6, and the Griffin table path.

Serpent 2 is issued by VTT for non-commercial research after registration.
It is not redistributed here. MCNP6 is distributed by RSICC and is export
controlled. It is not an open-source package and it is not downloaded.
Griffin remains the MOOSE reactor-physics table written by coupling.neams.
A binary is launched only when it is already on PATH and the config asks
for execution. A missing binary leaves ``executed`` false and the Python
Bateman inventory still runs.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from solcascade.geometry.shell import Shell
from solcascade.materials.library import Material
from solcascade.units import meters_to_centimeters


def _binary(names: tuple[str, ...]) -> str | None:
    for name in names:
        found = shutil.which(name)
        if found:
            return found
    return None


def _zaid(symbol: str) -> str:
    table = {"U238": "92238", "Th232": "90232", "Pb208": "82208"}
    if symbol not in table:
        raise ValueError(f"no ZAID is defined for {symbol}")
    return table[symbol]


def write_serpent_input(
    path: Path,
    shell: Shell,
    material: Material,
    neutron_energy_mev: float,
) -> None:
    """Fertile shell only. Inner wall compositions stay in the OpenMC model."""
    fertile = [layer for layer in shell.layers if layer.material == material.name]
    if len(fertile) != 1:
        raise ValueError("Serpent input expects one fertile layer")
    layer = fertile[0]
    lines = [
        "% SolCascade fertile blanket. Serpent is not bundled.",
        "% ACE libraries are a local file. The .09c suffix is a placeholder.",
        "% The tungsten, lithium-lead, and steel driver is the OpenMC model.",
        f"% Source energy {neutron_energy_mev:.6g} MeV is the DT neutron from constants.py.",
        f"set title \"solcascade {material.name}\"",
        "",
        f"mat fuel {-material.density_kg_m3 / 1000.0:.6g}",
    ]
    for symbol, fraction in material.composition:
        lines.append(f"{_zaid(symbol)}.09c {fraction:.8g}")
    inner_cm = meters_to_centimeters(layer.inner_radius_m)
    outer_cm = meters_to_centimeters(layer.outer_radius_m)
    lines.extend(
        [
            "",
            f"surf 1 sph 0.0 0.0 0.0 {inner_cm:.8g}",
            f"surf 2 sph 0.0 0.0 0.0 {outer_cm:.8g}",
            "",
            "cell 1 0 fuel 1 -2",
            "cell 2 0 outside 2",
            "cell 3 0 void -1",
            "",
            "set pop 1000 20 20",
            f"src 1 n se {neutron_energy_mev:.8g} sp 0.0 0.0 0.0",
            "det 1 dr -6 void 1",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def write_mcnp_input(
    path: Path,
    shell: Shell,
    material: Material,
    neutron_energy_mev: float,
) -> None:
    fertile = [layer for layer in shell.layers if layer.material == material.name]
    if len(fertile) != 1:
        raise ValueError("MCNP input expects one fertile layer")
    layer = fertile[0]
    inner_cm = meters_to_centimeters(layer.inner_radius_m)
    outer_cm = meters_to_centimeters(layer.outer_radius_m)
    composition = " ".join(
        f"{_zaid(symbol)} {fraction:.8g}" for symbol, fraction in material.composition
    )
    text = "\n".join(
        [
            f"SolCascade {material.name} blanket, DT neutron point source",
            "c MCNP6 is not bundled. This file is an input, not a run.",
            "c The inner driver remains the OpenMC model.",
            f"1 1 {-material.density_kg_m3 / 1000.0:.6g} 1 -2",
            "2 0 -1",
            "3 0 2",
            "",
            f"1 so {inner_cm:.8g}",
            f"2 so {outer_cm:.8g}",
            "",
            "mode n",
            "imp:n 1 1 0",
            f"m1 {composition}",
            f"sdef pos=0 0 0 erg={neutron_energy_mev:.8g} par=n",
            "nps 1000",
            "f4:n 1",
            "fm4 1.0 1 -6",
            "",
        ]
    )
    path.write_text(text, encoding="utf-8")


def griffin_record(fuel_type: str) -> dict:
    return {
        "backend": "griffin",
        "executed": False,
        "binary": None,
        "reason": (
            "Griffin stays the table handoff in griffin_interface.json. "
            "SolCascade does not launch a Griffin binary."
        ),
        "fuel_type": fuel_type,
        "table": "griffin_interface.json",
    }


def _run_if_present(
    binary: str | None,
    workdir: Path,
    input_name: str,
    execute: bool,
    missing_reason: str,
) -> dict:
    record = {
        "executed": False,
        "binary": binary,
        "input": input_name,
    }
    if not execute:
        record["reason"] = "fission.execute is false; the input was written"
        return record
    if binary is None:
        record["reason"] = missing_reason
        return record
    try:
        completed = subprocess.run(
            [binary, input_name],
            cwd=workdir,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
    except Exception as exc:
        record["reason"] = str(exc)
        return record
    record["returncode"] = completed.returncode
    record["executed"] = completed.returncode == 0
    if completed.returncode != 0:
        record["reason"] = (completed.stderr or completed.stdout or "non-zero exit")[:500]
    return record


def write_fission_backend(
    workdir: Path,
    backend: str,
    execute: bool,
    shell: Shell,
    material: Material | None,
    neutron_energy_mev: float,
    fuel_type: str,
) -> dict:
    if backend == "none":
        return {"backend": "none", "executed": False, "reason": "fission.backend is none"}
    if backend == "griffin":
        return griffin_record(fuel_type)
    if material is None:
        raise ValueError(f"fission backend {backend} needs a fertile material")
    if backend == "serpent":
        path = workdir / "blanket_serpent.txt"
        write_serpent_input(path, shell, material, neutron_energy_mev)
        record = _run_if_present(
            _binary(("sss2", "serpent")),
            workdir,
            path.name,
            execute,
            "Serpent 2 binary (sss2 or serpent) is not on PATH",
        )
        record["backend"] = "serpent"
        record["license_note"] = (
            "Serpent 2 is available from VTT for non-commercial research. "
            "It is not redistributed and it was not downloaded."
        )
        return record
    if backend == "mcnp":
        path = workdir / "blanket_mcnp.txt"
        write_mcnp_input(path, shell, material, neutron_energy_mev)
        record = _run_if_present(
            _binary(("mcnp6", "mcnp")),
            workdir,
            path.name,
            execute,
            "MCNP6 binary is not on PATH",
        )
        record["backend"] = "mcnp"
        record["license_note"] = (
            "MCNP6 is distributed by RSICC and is export controlled. "
            "It is not open source and it was not downloaded."
        )
        return record
    raise ValueError(f"unknown fission backend {backend!r}")
