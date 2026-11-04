"""Bateman inventory for one fertile chain.

U-238 captures to U-239, which beta-decays to Np-239 and then Pu-239.
Th-232 captures to Th-233, which beta-decays to Pa-233 and then U-233.
The pulse is treated as a burst: it is short next to these half-lives, so
decay runs only in the gap between pulses. Fission removes a parent nucleus.
Decay conserves the heavy-atom count. Fissioned nuclei leave the heavy inventory.
"""

from __future__ import annotations

import math

from solcascade.constants import (
    HALF_LIFE_NP239_S,
    HALF_LIFE_PA233_S,
    HALF_LIFE_TH233_S,
    HALF_LIFE_U239_S,
    decay_constant_per_s,
)

Chain = dict[str, float]

U238_CHAIN = ("U238", "U239", "Np239", "Pu239")
TH232_CHAIN = ("Th232", "Th233", "Pa233", "U233")

_LAMBDAS = {
    "u238": (
        decay_constant_per_s(HALF_LIFE_U239_S),
        decay_constant_per_s(HALF_LIFE_NP239_S),
    ),
    "th232": (
        decay_constant_per_s(HALF_LIFE_TH233_S),
        decay_constant_per_s(HALF_LIFE_PA233_S),
    ),
}


def chain_names(fuel_type: str) -> tuple[str, str, str, str]:
    if fuel_type == "u238":
        return U238_CHAIN
    if fuel_type == "th232":
        return TH232_CHAIN
    raise ValueError(f"no depletion chain for fuel {fuel_type!r}")


def _decay_daughters(
    first: float,
    second: float,
    third: float,
    lambda_first: float,
    lambda_second: float,
    interpulse_s: float,
) -> tuple[float, float, float]:
    """Linear chain A -> B -> C with C stable. The sum is unchanged."""
    if interpulse_s < 0.0:
        raise ValueError("interpulse interval must be >= 0")
    if interpulse_s == 0.0:
        return first, second, third
    decay_first = math.exp(-lambda_first * interpulse_s)
    decay_second = math.exp(-lambda_second * interpulse_s)
    next_first = first * decay_first
    gap = lambda_second - lambda_first
    if abs(gap) <= 1.0e-12 * max(lambda_first, lambda_second, 1.0):
        next_second = (second + first * lambda_first * interpulse_s) * decay_first
    else:
        next_second = second * decay_second + first * lambda_first / gap * (
            decay_first - decay_second
        )
    next_third = third + (first + second - next_first - next_second)
    return next_first, next_second, next_third


def deplete_pulses(
    fuel_type: str,
    parent_atoms: float,
    fissions_per_pulse: float,
    captures_per_pulse: float,
    pulses: int,
    interpulse_s: float,
) -> dict:
    """Advance a fertile inventory. The returned fissioned count is cumulative."""
    if parent_atoms < 0.0:
        raise ValueError("parent inventory must be >= 0")
    if fissions_per_pulse < 0.0 or captures_per_pulse < 0.0:
        raise ValueError("reaction counts must be >= 0")
    if pulses < 0:
        raise ValueError("pulse count must be >= 0")
    names = chain_names(fuel_type)
    lambda_first, lambda_second = _LAMBDAS[fuel_type]
    parent, first, second, third = float(parent_atoms), 0.0, 0.0, 0.0
    fissioned = 0.0
    history = []
    initial = parent_atoms
    for index in range(pulses):
        removed = fissions_per_pulse + captures_per_pulse
        if removed > parent + 1.0e-6:
            raise ValueError("capture plus fission exceeds the parent inventory")
        parent -= removed
        first += captures_per_pulse
        fissioned += fissions_per_pulse
        first, second, third = _decay_daughters(
            first, second, third, lambda_first, lambda_second, interpulse_s
        )
        inventory = {
            names[0]: parent,
            names[1]: first,
            names[2]: second,
            names[3]: third,
        }
        history.append(
            {
                "pulse": index + 1,
                "atoms": inventory,
                "fissioned_cumulative": fissioned,
                "heavy_plus_fissioned": parent + first + second + third + fissioned,
            }
        )
    final = {
        names[0]: parent,
        names[1]: first,
        names[2]: second,
        names[3]: third,
    }
    return {
        "fuel_type": fuel_type,
        "pulses": pulses,
        "interpulse_s": interpulse_s,
        "initial_parent_atoms": initial,
        "fissioned": fissioned,
        "inventory": final,
        "history": history,
        "conserved_atoms": parent + first + second + third + fissioned,
    }
