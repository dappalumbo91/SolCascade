"""OpenMC model writer and optional runner.

The generated script is ordinary OpenMC Python API: spherical shells in
centimeters, a monoenergetic isotropic DT neutron source, and a cell heating
tally. Running it needs an OpenMC install and a continuous-energy library in
``OPENMC_CROSS_SECTIONS``. No account is required. The download of public
nuclear data is documented in the README and is not performed here.

Stable OpenMC parallelizes with OpenMP and MPI. This module passes a thread
count through. It does not pass a CUDA flag. Use a GPU build only if the
installed ``openmc`` executable was built with one.
"""

from __future__ import annotations

import os
import shutil
import textwrap
from pathlib import Path

from solcascade.constants import dt_reaction
from solcascade.geometry.shell import Shell
from solcascade.materials.library import Material
from solcascade.units import meters_to_centimeters


def openmc_status() -> dict:
    try:
        import openmc
    except ImportError as exc:
        return {"importable": False, "reason": f"openmc is not installed ({exc})"}
    cross_sections = os.environ.get("OPENMC_CROSS_SECTIONS")
    executable = shutil.which("openmc")
    status = {
        "importable": True,
        "version": getattr(openmc, "__version__", "unknown"),
        "cross_sections": cross_sections,
        "executable": executable,
    }
    if not cross_sections:
        status["reason"] = "OPENMC_CROSS_SECTIONS is not set"
    elif not executable:
        status["reason"] = "the openmc executable is not on PATH"
    return status


def write_model_script(
    path: Path,
    shell: Shell,
    library: dict[str, Material],
    neutron_energy_ev: float,
    particles: int,
    batches: int,
    photon_transport: bool,
    seed: int,
    source_rate_per_s: float,
    temperature_k: float,
) -> None:
    """Write a standalone OpenMC build script. It imports openmc when run.

    ``temperature_k`` is the shell temperature used to select cross sections.
    The heat solve has not run yet, so this is the shell's starting temperature.
    """
    if float(temperature_k) <= 0.0:
        raise ValueError("OpenMC material temperature must be positive")
    lines = [
        "\"\"\"Generated SolCascade OpenMC model. Run with OpenMC on PYTHONPATH.\"\"\"",
        "import openmc",
        "",
        "materials = []",
    ]
    for material_name in dict.fromkeys(cell.material for cell in shell.cells):
        material = library[material_name]
        density_g_cm3 = material.density_kg_m3 / 1000.0
        lines.append(f"mat_{material_name} = openmc.Material(name={material_name!r})")
        for symbol, fraction in material.composition:
            if any(character.isdigit() for character in symbol):
                adder = "add_nuclide"
            else:
                adder = "add_element"
            lines.append(
                f"mat_{material_name}.{adder}({symbol!r}, {fraction!r}, percent_type='wo')"
            )
        lines.append(f"mat_{material_name}.set_density('g/cm3', {density_g_cm3!r})")
        lines.append(f"mat_{material_name}.temperature = {float(temperature_k)!r}")
        lines.append(f"materials.append(mat_{material_name})")
    lines.append("mats = openmc.Materials(materials)")
    lines.append("")
    # Surfaces in centimeters. The outermost sphere is a vacuum boundary.
    radii = [shell.inner_radius_m]
    radii.extend(layer.outer_radius_m for layer in shell.layers)
    for index, radius_m in enumerate(radii):
        boundary = "vacuum" if index in {0, len(radii) - 1} else "transmission"
        # The inner cavity surface should not be vacuum: neutrons must cross
        # it into the shell. Only the outer surface is a vacuum boundary.
        if index != len(radii) - 1:
            boundary = "transmission"
        lines.append(
            "surf_{0} = openmc.Sphere(r={1!r}, boundary_type={2!r})".format(
                index, meters_to_centimeters(radius_m), boundary
            )
        )
    lines.append("cells = []")
    lines.append("cavity = openmc.Cell(name='cavity', region=-surf_0)")
    lines.append("cells.append(cavity)")
    layer_cells = []
    for index, layer in enumerate(shell.layers):
        name = f"cell_{layer.name}"
        layer_cells.append(name)
        lines.append(
            f"{name} = openmc.Cell(name={layer.name!r}, fill=mat_{layer.material}, "
            f"region=+surf_{index} & -surf_{index + 1})"
        )
        lines.append(f"cells.append({name})")
    lines.append("geometry = openmc.Geometry(openmc.Universe(cells=cells))")
    lines.append("source = openmc.IndependentSource()")
    lines.append("source.space = openmc.stats.Point((0.0, 0.0, 0.0))")
    lines.append("source.angle = openmc.stats.Isotropic()")
    lines.append(
        f"source.energy = openmc.stats.Discrete([{neutron_energy_ev!r}], [1.0])"
    )
    lines.append("source.particle = 'neutron'")
    lines.append("settings = openmc.Settings()")
    lines.append("settings.run_mode = 'fixed source'")
    lines.append(f"settings.particles = {int(particles)}")
    lines.append(f"settings.batches = {int(batches)}")
    lines.append(f"settings.seed = {int(seed)}")
    lines.append(f"settings.photon_transport = {bool(photon_transport)}")
    lines.append("settings.source = source")
    lines.append(
        f"# Source rate used by SolCascade when converting the tally: {float(source_rate_per_s)!r} 1/s"
    )
    cell_list = ", ".join(layer_cells)
    lines.append(f"cell_filter = openmc.CellFilter([{cell_list}])")
    lines.append("tally = openmc.Tally(name='shell-heating')")
    lines.append("tally.filters = [cell_filter]")
    lines.append("tally.scores = ['heating', 'flux']")
    lines.append("tally_fission = openmc.Tally(name='fission-rate')")
    lines.append("tally_fission.filters = [cell_filter]")
    lines.append("tally_fission.scores = ['fission']")
    lines.append("tally_tritium = openmc.Tally(name='tritium-production')")
    lines.append("tally_tritium.filters = [cell_filter]")
    lines.append("tally_tritium.scores = ['H3-production']")
    lines.append("tallies = openmc.Tallies([tally, tally_fission, tally_tritium])")
    lines.append("model = openmc.Model(geometry, mats, settings, tallies)")
    lines.append("model.export_to_model_xml()")
    lines.append("")
    text = "\n".join(lines) + "\n"
    path.write_text(text, encoding="utf-8")


def generated_script_text_mentions(shell: Shell) -> str:
    """Used by tests that do not execute OpenMC. Returns the outer radius in cm."""
    return repr(meters_to_centimeters(shell.outer_radius_m))


def run_openmc(workdir: Path, threads: int | None, executable: str | None) -> dict:
    """Execute ``openmc`` in ``workdir``. The caller has already exported XML."""
    status = openmc_status()
    if not status["importable"]:
        return {"executed": False, **status}
    if "reason" in status and status.get("cross_sections") in {None, ""}:
        return {"executed": False, **status}
    import openmc

    try:
        openmc.run(
            threads=threads,
            cwd=str(workdir),
            openmc_exec=executable or "openmc",
            output=True,
        )
    except Exception as exc:
        return {"executed": False, "reason": str(exc), **status}
    return {"executed": True, **status}


def neutron_source_energy_ev() -> float:
    return dt_reaction().neutron_energy_ev


def script_header_comment() -> str:
    return textwrap.dedent(
        """
        The heating tally is eV per source particle. Multiply by the source
        rate and by the joule-per-eV constant to get watts. SolCascade does
        that conversion in the runner when a statepoint is present.
        """
    ).strip()
