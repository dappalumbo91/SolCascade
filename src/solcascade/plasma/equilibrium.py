"""Vacuum coil equilibrium, with an optional TokaMaker handoff.

The default equilibrium is the vacuum field of the coil set plus the
zero-dimensional plasma pressure. That is not a Grad-Shafranov solution.
When the Open FUSION Toolkit Python module is installed and a mesh path is
configured, :mod:`solcascade.coupling.tokamaker` can run TokaMaker. This
module only reports which description is in use.
"""

from __future__ import annotations

from solcascade.plasma.metrics import PlasmaHistory


def equilibrium_summary(history: PlasmaHistory, requested: str) -> dict:
    peak = max(history.points, key=lambda point: point.fusion_power_w)
    return {
        "requested": requested,
        "solver": "vacuum_coil_field",
        "grad_shafranov": False,
        "peak_vacuum_field_T": max(point.vacuum_field_t for point in history.points),
        "peak_internal_field_T": max(point.internal_field_t for point in history.points),
        "peak_beta": peak.beta,
        "peak_pressure_Pa": peak.pressure_pa,
        "note": (
            "Vacuum field of the axisymmetric coils plus a uniform-sphere "
            "pressure. A TokaMaker Grad-Shafranov solve is not part of this "
            "result unless the coupling script is run separately."
        ),
    }
