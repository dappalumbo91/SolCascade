"""Constants are derived once, and only in the constants module."""

import math
from pathlib import Path

from solcascade.constants import (
    FINE_STRUCTURE_CONSTANT,
    JOULE_PER_EV,
    PLANCK_H_J_S,
    SPEED_OF_LIGHT_M_S,
    VACUUM_PERMEABILITY_H_PER_M,
    VACUUM_PERMITTIVITY_F_PER_M,
    dt_q_kev,
    dt_reaction,
)

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "solcascade"


def test_vacuum_permeability_matches_the_definition():
    expected = (
        2.0 * FINE_STRUCTURE_CONSTANT * PLANCK_H_J_S
        / (JOULE_PER_EV ** 2 * SPEED_OF_LIGHT_M_S)
    )
    assert VACUUM_PERMEABILITY_H_PER_M == expected
    # mu0 * eps0 * c^2 is 1 by definition of eps0. Left-associated float
    # multiplication misses that identity by one unit in the last place.
    product = (
        VACUUM_PERMEABILITY_H_PER_M
        * VACUUM_PERMITTIVITY_F_PER_M
        * SPEED_OF_LIGHT_M_S ** 2
    )
    assert math.isclose(product, 1.0, rel_tol=1.0e-15, abs_tol=0.0)


def test_dt_partition_sums_to_q_and_favors_the_neutron():
    reaction = dt_reaction()
    assert reaction.q_kev == dt_q_kev()
    assert reaction.neutron_energy_kev + reaction.alpha_energy_kev == reaction.q_kev
    assert 17.5 < reaction.q_kev / 1.0e3 < 17.7
    assert reaction.neutron_energy_kev > reaction.alpha_energy_kev
    assert abs(reaction.neutron_fraction + reaction.alpha_fraction - 1.0) < 1.0e-12


def test_measured_digits_live_only_in_the_constants_module():
    needles = ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9")
    offenders = []
    for path in SRC.rglob("*.py"):
        if path.name == "constants.py":
            continue
        text = path.read_text(encoding="utf-8")
        for needle in needles:
            if needle in text:
                offenders.append(f"{path.relative_to(ROOT)} contains {needle}")
    assert offenders == []
