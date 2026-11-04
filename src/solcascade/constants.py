"""Shared physical constants and published coefficients.

This is the only module that stores CODATA values, AME2020 mass excesses,
and published reactivity-fit coefficients. Other modules import the names
below. Derived quantities are computed here so they are not retyped.

Sources
-------
CODATA 2022: Mohr, Newell, Taylor, and Tiesinga, Rev. Mod. Phys. 97, 025002 (2025).
SI exact constants: 2019 redefinition.
AME2020: Wang et al., Chin. Phys. C 45, 030003 (2021).
Bosch-Hale DT and D-He3 reactivity: H.-S. Bosch and G. M. Hale,
Nucl. Fusion 32, 611 (1992).
Nevins-Swain low-temperature p-11B reactivity: W. M. Nevins and R. Swain,
Nucl. Fusion 40, 865 (2000), as reprinted in equations (9)-(15) and the
Nevins column of table 1 by A. Tentori and F. Belloni, Nucl. Fusion 63,
086001 (2023).
Bohm coefficient: the historical comparison D = kT / (16 e B). The ledger
walk uses 4 π from the cyclotron period, computed in plasma.metrics.
Triple-product benchmark: order of magnitude quoted by Wesson, Tokamaks, and
Freidberg, Plasma Physics and Fusion Energy, near 10–20 keV for DT.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# Exact by the 2019 SI definition.
SPEED_OF_LIGHT_M_S = 299_792_458.0
ELEMENTARY_CHARGE_C = 1.602_176_634e-19
PLANCK_H_J_S = 6.626_070_15e-34
BOLTZMANN_J_K = 1.380_649e-23
AVOGADRO_PER_MOL = 6.022_140_76e23

# CODATA 2022 measured inverse fine-structure constant.
ALPHA_INVERSE = 137.035_999_177
FINE_STRUCTURE_CONSTANT = 1.0 / ALPHA_INVERSE

# mu_0 = 2 h alpha / (e^2 c). Not the pre-2019 exact 4 pi e-7.
VACUUM_PERMEABILITY_H_PER_M = (
    2.0 * FINE_STRUCTURE_CONSTANT * PLANCK_H_J_S
    / (ELEMENTARY_CHARGE_C ** 2 * SPEED_OF_LIGHT_M_S)
)
VACUUM_PERMITTIVITY_F_PER_M = 1.0 / (
    VACUUM_PERMEABILITY_H_PER_M * SPEED_OF_LIGHT_M_S ** 2
)

REDUCED_PLANCK_J_S = PLANCK_H_J_S / (2.0 * math.pi)
JOULE_PER_EV = ELEMENTARY_CHARGE_C
EV_PER_JOULE = 1.0 / JOULE_PER_EV

# AME2020 atomic mass excesses, keV, and that evaluation's keV-per-u factor.
# Q(D+T) is computed from the excesses. Electron count matches: the deuteron
# and triton atoms contribute two electrons, and the helium atom contributes
# two, so the atomic excesses give the correct reaction Q.
AME2020_U_KEV = 931_494.10242
MASS_EXCESS_NEUTRON_KEV = 8_071.31713
MASS_EXCESS_DEUTERIUM_KEV = 13_135.72176
MASS_EXCESS_TRITIUM_KEV = 14_949.80993
MASS_EXCESS_HELIUM4_KEV = 2_424.91561

# CODATA 2022 electron mass. Used only to remove atomic electrons from the
# helium rest mass before the two-body kinetic split. Electron binding
# (~79 eV) is negligible next to the keV excesses and is not applied.
ELECTRON_MASS_KG = 9.109_383_7139e-31
ELECTRON_REST_ENERGY_J = ELECTRON_MASS_KG * SPEED_OF_LIGHT_M_S ** 2
ELECTRON_REST_ENERGY_KEV = ELECTRON_REST_ENERGY_J * EV_PER_JOULE / 1.0e3

# Monatomic plasma: three translational degrees of freedom, gamma = 1 + 2/f.
MONATOMIC_TRANSLATIONAL_DEGREES = 3.0
MONATOMIC_ADIABATIC_INDEX = 1.0 + 2.0 / MONATOMIC_TRANSLATIONAL_DEGREES

# Historical Bohm diffusion coefficient in D_B = kT / (BOHM_DENOMINATOR * e * B).
BOHM_DENOMINATOR = 16.0

# Order-of-magnitude DT ignition benchmark, keV s / m^3, not a sharp threshold.
DT_TRIPLE_PRODUCT_BENCHMARK_KEV_S_M3 = 3.0e21

# Bosch and Hale (1992) DT fit. Valid on the temperature interval below.
# Result of the published formula is cm^3/s. Coefficients are not CODATA.
BOSCH_HALE_DT_BG_SQRT_KEV = 34.3827
BOSCH_HALE_DT_REDUCED_MASS_ENERGY_KEV = 1_124_656.0
BOSCH_HALE_DT_C1 = 1.17302e-9
BOSCH_HALE_DT_C2 = 1.51361e-2
BOSCH_HALE_DT_C3 = 7.51886e-2
BOSCH_HALE_DT_C4 = 4.60643e-3
BOSCH_HALE_DT_C5 = 1.35000e-2
BOSCH_HALE_DT_C6 = -1.06750e-4
BOSCH_HALE_DT_C7 = 1.36600e-5
BOSCH_HALE_DT_T_MIN_KEV = 0.2
BOSCH_HALE_DT_T_MAX_KEV = 100.0

# Bosch and Hale (1992) 3He(d,p)4He fit. Same functional form as DT.
# B_G is stored the way the DT row stores it: the formula squares it.
# The reactivity interval is 0.5–190 keV. The wider cross-section energy
# window is not used. Coefficients are the published fit, not CODATA.
# Result of the published formula is cm^3/s.
BOSCH_HALE_DHE3_BG_SQRT_KEV = 68.7508
BOSCH_HALE_DHE3_REDUCED_MASS_ENERGY_KEV = 1_124_572.0
BOSCH_HALE_DHE3_C1 = 5.51036e-10
BOSCH_HALE_DHE3_C2 = 6.41918e-3
BOSCH_HALE_DHE3_C3 = -2.02896e-3
BOSCH_HALE_DHE3_C4 = -1.91080e-5
BOSCH_HALE_DHE3_C5 = 1.35776e-4
BOSCH_HALE_DHE3_C6 = 0.0
BOSCH_HALE_DHE3_C7 = 0.0
BOSCH_HALE_DHE3_T_MIN_KEV = 0.5
BOSCH_HALE_DHE3_T_MAX_KEV = 190.0

# Nevins and Swain (2000) low-temperature p-11B reactivity.
# Tentori and Belloni reprint that branch as equations (9)-(15). The
# S-factor coefficients are the Nevins column of their table 1, for
# E <= 0.400 MeV. E_G and mu c^2 are the values stated with those equations.
# The 2000 abstract gives this branch for 1 keV to 70 keV. The same paper's
# high-temperature Padé starts at 50 keV and is not stored here.
# C0, C1, C2, and A_L are MeV b. E_L and delta E_L are keV.
NEVINS_SWAIN_PB11_C0_MEV_B = 197.0
NEVINS_SWAIN_PB11_C1_MEV_B = 0.240
NEVINS_SWAIN_PB11_C2_MEV_B = 2.31e-4
NEVINS_SWAIN_PB11_AL_MEV_B = 1.82e4
NEVINS_SWAIN_PB11_EL_KEV = 148.0
NEVINS_SWAIN_PB11_DEL_KEV = 2.35
NEVINS_SWAIN_PB11_EG_MEV = 22.589
NEVINS_SWAIN_PB11_MUC2_MEV = 859.526
NEVINS_SWAIN_PB11_T_MIN_KEV = 1.0
NEVINS_SWAIN_PB11_T_MAX_KEV = 70.0

# Conventional recoverable energy per fission in the hybrid ledger.
# This is the engineering quantum near 200 MeV. It is not an isotope Q from
# ENDF. The local ENDF/B-VIII.0 files store a recoverable curve; this ledger
# does not evaluate that curve.
FISSION_RECOVERABLE_ENERGY_MEV = 200.0

# ENSDF half-lives for the fertile breeding chains. ENSDF April 2014, as
# listed in the IAEA comparison against Nubase2012, and named in the NNDC
# ENSDF index for 233Th beta decay (21.83 min) and 233Pa beta decay
# (26.975 d). Pu-239 and U-233 are the chain ends. Their half-lives are
# 24.11 kyr and 159.2 kyr, so the Bateman step treats them as stable.
HALF_LIFE_U239_S = 23.45 * 60.0
HALF_LIFE_NP239_S = 2.356 * 86400.0
HALF_LIFE_TH233_S = 21.83 * 60.0
HALF_LIFE_PA233_S = 26.975 * 86400.0

# Nuclide molar masses used to turn the metal density into a number density
# when the ENDF/B-VIII.0 point cross sections were written as macroscopic
# totals in materials.yaml. The depletion inventory uses the same masses.
MOLAR_MASS_U238_G = 238.0507882
MOLAR_MASS_TH232_G = 232.0380553

# Positive-mass warp shell scored in the local warp ledger: compactness 1/2,
# outer radius 110 m, mass 3.703e28 kg. The ledger prints the rest energy
# as 3.328e45 J. This mass is the source number. The joule requirement is
# the product with c^2, so the catalog does not retype that product.
WARP_SHELL_MASS_KG = 3.703e28


def fission_recoverable_energy_j() -> float:
    """Ledger energy deposited as heat by one fission."""
    return FISSION_RECOVERABLE_ENERGY_MEV * 1.0e6 * JOULE_PER_EV


def decay_constant_per_s(half_life_s: float) -> float:
    if half_life_s <= 0.0:
        raise ValueError("half-life must be positive")
    return math.log(2.0) / half_life_s


def warp_shell_energy_j() -> float:
    """Rest energy of the scored 110 m warp shell, in joules."""
    return WARP_SHELL_MASS_KG * SPEED_OF_LIGHT_M_S ** 2


def dt_q_kev() -> float:
    """DT reaction Q from AME2020 mass excesses, in keV."""
    return (
        MASS_EXCESS_DEUTERIUM_KEV
        + MASS_EXCESS_TRITIUM_KEV
        - MASS_EXCESS_HELIUM4_KEV
        - MASS_EXCESS_NEUTRON_KEV
    )


def _rest_energy_kev(mass_number: int, excess_kev: float) -> float:
    return mass_number * AME2020_U_KEV + excess_kev


def ion_mass_kg(mass_number: int, atomic_excess_kev: float, electrons: int) -> float:
    """Ion mass from an AME2020 atomic excess, with electron rest mass removed."""
    joule = _rest_energy_kev(mass_number, atomic_excess_kev) * JOULE_PER_EV * 1.0e3
    return (joule - electrons * ELECTRON_REST_ENERGY_J) / SPEED_OF_LIGHT_M_S ** 2


@dataclass(frozen=True)
class DTReaction:
    """Cold-target energy partition of D + T -> alpha + n."""

    q_kev: float
    q_joule: float
    neutron_energy_kev: float
    alpha_energy_kev: float
    neutron_energy_joule: float
    alpha_energy_joule: float
    neutron_fraction: float
    alpha_fraction: float

    @property
    def neutron_energy_ev(self) -> float:
        return self.neutron_energy_kev * 1.0e3


def dt_reaction() -> DTReaction:
    """Split Q by two-body kinematics with the reactants at rest.

    The neutron carries the alpha rest-mass fraction of Q. Rest masses are
    the AME2020 atomic masses, with two electron rest energies removed from
    helium so the alpha is a bare nucleus. The neutron excess is already the
    bare neutron.
    """
    q_kev = dt_q_kev()
    neutron_kev = _rest_energy_kev(1, MASS_EXCESS_NEUTRON_KEV)
    helium_atomic_kev = _rest_energy_kev(4, MASS_EXCESS_HELIUM4_KEV)
    alpha_kev = helium_atomic_kev - 2.0 * ELECTRON_REST_ENERGY_KEV
    total = alpha_kev + neutron_kev
    neutron_fraction = alpha_kev / total
    alpha_fraction = neutron_kev / total
    neutron_energy_kev = q_kev * neutron_fraction
    alpha_energy_kev = q_kev * alpha_fraction
    joule_per_kev = JOULE_PER_EV * 1.0e3
    return DTReaction(
        q_kev=q_kev,
        q_joule=q_kev * joule_per_kev,
        neutron_energy_kev=neutron_energy_kev,
        alpha_energy_kev=alpha_energy_kev,
        neutron_energy_joule=neutron_energy_kev * joule_per_kev,
        alpha_energy_joule=alpha_energy_kev * joule_per_kev,
        neutron_fraction=neutron_fraction,
        alpha_fraction=alpha_fraction,
    )
