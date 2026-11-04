"""Griffin and Bison interfaces.

Griffin and Bison are MOOSE applications distributed through the INL NEAMS
program. SolCascade does not download them and does not require an account.
The files written here are the tables those applications would read: power
density for Griffin and temperature for a Bison-style structural pass.

The open-source structural numbers in the pipeline come from the Python
Lamé and thermal-stress solution. A MOOSE solid-mechanics input is not
generated, because the 1D spherical eigenstrain syntax is too easy to get
wrong without a local MOOSE build to check it. The tested stress result is
the Python one.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from solcascade.geometry.shell import Shell


def write_tables(
    directory: Path,
    shell: Shell,
    source_w_m3: np.ndarray,
    temperature_k: np.ndarray,
    von_mises_pa: np.ndarray,
) -> dict:
    directory.mkdir(parents=True, exist_ok=True)
    rows = []
    for cell, source, temperature, mises in zip(
        shell.cells, source_w_m3, temperature_k, von_mises_pa, strict=True
    ):
        rows.append(
            {
                "cell": cell.index,
                "layer": cell.layer,
                "material": cell.material,
                "r_inner_m": cell.inner_radius_m,
                "r_outer_m": cell.outer_radius_m,
                "power_density_W_m3": float(source),
                "temperature_K": float(temperature),
                "von_mises_Pa": None if np.isnan(mises) else float(mises),
            }
        )
    griffin = {
        "application": "Griffin",
        "role": "neutron and photon transport source for a coupled MOOSE run",
        "status": "table_only",
        "note": (
            "One-group removal results from the Python surrogate, or OpenMC "
            "when that backend ran. Griffin DFEM needs its own cross-section "
            "library and is not executed by SolCascade."
        ),
        "cells": rows,
    }
    bison = {
        "application": "Bison",
        "role": "thermomechanical response of the solid shell",
        "status": "table_only",
        "note": (
            "Temperature and von Mises stress from the Python solvers. "
            "Liquid layers have a null stress. Bison is not executed."
        ),
        "cells": rows,
    }
    griffin_path = directory / "griffin_interface.json"
    bison_path = directory / "bison_interface.json"
    griffin_path.write_text(json.dumps(griffin, indent=2), encoding="utf-8")
    bison_path.write_text(json.dumps(bison, indent=2), encoding="utf-8")
    return {"griffin": str(griffin_path), "bison": str(bison_path), "executed": False}
