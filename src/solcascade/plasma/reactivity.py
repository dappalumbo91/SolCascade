"""Published thermal reactivities.

Each fit is used only inside the temperature interval published with it.
Outside that interval the function raises rather than extrapolating.
Bosch-Hale covers DT and D-He3. The p-11B row is the Nevins-Swain
low-temperature branch, 1 keV to 70 keV. The high-temperature Padé is
not evaluated.
"""

from __future__ import annotations

import math

from solcascade.constants import (
    BOSCH_HALE_DHE3_BG_SQRT_KEV,
    BOSCH_HALE_DHE3_C1,
    BOSCH_HALE_DHE3_C2,
    BOSCH_HALE_DHE3_C3,
    BOSCH_HALE_DHE3_C4,
    BOSCH_HALE_DHE3_C5,
    BOSCH_HALE_DHE3_C6,
    BOSCH_HALE_DHE3_C7,
    BOSCH_HALE_DHE3_REDUCED_MASS_ENERGY_KEV,
    BOSCH_HALE_DHE3_T_MAX_KEV,
    BOSCH_HALE_DHE3_T_MIN_KEV,
    BOSCH_HALE_DT_BG_SQRT_KEV,
    BOSCH_HALE_DT_C1,
    BOSCH_HALE_DT_C2,
    BOSCH_HALE_DT_C3,
    BOSCH_HALE_DT_C4,
    BOSCH_HALE_DT_C5,
    BOSCH_HALE_DT_C6,
    BOSCH_HALE_DT_C7,
    BOSCH_HALE_DT_REDUCED_MASS_ENERGY_KEV,
    BOSCH_HALE_DT_T_MAX_KEV,
    BOSCH_HALE_DT_T_MIN_KEV,
    JOULE_PER_EV,
    NEVINS_SWAIN_PB11_AL_MEV_B,
    NEVINS_SWAIN_PB11_C0_MEV_B,
    NEVINS_SWAIN_PB11_C1_MEV_B,
    NEVINS_SWAIN_PB11_C2_MEV_B,
    NEVINS_SWAIN_PB11_DEL_KEV,
    NEVINS_SWAIN_PB11_EG_MEV,
    NEVINS_SWAIN_PB11_EL_KEV,
    NEVINS_SWAIN_PB11_MUC2_MEV,
    NEVINS_SWAIN_PB11_T_MAX_KEV,
    NEVINS_SWAIN_PB11_T_MIN_KEV,
    SPEED_OF_LIGHT_M_S,
)
from solcascade.units import BARN_M2, CUBIC_METERS_PER_CUBIC_CENTIMETER, EV_PER_KEV


def bosch_hale_reactivity_m3_s(
    temperature_kev: float,
    *,
    bg_sqrt_kev: float,
    reduced_mass_energy_kev: float,
    c1: float,
    c2: float,
    c3: float,
    c4: float,
    c5: float,
    c6: float,
    c7: float,
    t_min_kev: float,
    t_max_kev: float,
    reaction: str,
) -> float:
    """Thermal <sigma v> in m^3/s for one Bosch-Hale coefficient row."""
    temperature = float(temperature_kev)
    if temperature < t_min_kev or temperature > t_max_kev:
        raise ValueError(
            f"{reaction} reactivity is outside the Bosch-Hale interval "
            f"[{t_min_kev}, {t_max_kev}] keV: "
            f"{temperature} keV"
        )
    numerator = temperature * (c2 + temperature * (c4 + temperature * c6))
    denominator = 1.0 + temperature * (c3 + temperature * (c5 + temperature * c7))
    theta = temperature / (1.0 - numerator / denominator)
    if theta <= 0.0:
        raise ValueError("Bosch-Hale theta is not positive")
    xi = (bg_sqrt_kev ** 2 / (4.0 * theta)) ** (1.0 / 3.0)
    cm3_s = (
        c1
        * theta
        * math.sqrt(xi / (reduced_mass_energy_kev * temperature ** 3))
        * math.exp(-3.0 * xi)
    )
    return cm3_s * CUBIC_METERS_PER_CUBIC_CENTIMETER


def dt_reactivity_m3_s(temperature_kev: float) -> float:
    """Thermal <sigma v> for DT in m^3/s."""
    return bosch_hale_reactivity_m3_s(
        temperature_kev,
        bg_sqrt_kev=BOSCH_HALE_DT_BG_SQRT_KEV,
        reduced_mass_energy_kev=BOSCH_HALE_DT_REDUCED_MASS_ENERGY_KEV,
        c1=BOSCH_HALE_DT_C1,
        c2=BOSCH_HALE_DT_C2,
        c3=BOSCH_HALE_DT_C3,
        c4=BOSCH_HALE_DT_C4,
        c5=BOSCH_HALE_DT_C5,
        c6=BOSCH_HALE_DT_C6,
        c7=BOSCH_HALE_DT_C7,
        t_min_kev=BOSCH_HALE_DT_T_MIN_KEV,
        t_max_kev=BOSCH_HALE_DT_T_MAX_KEV,
        reaction="DT",
    )


def dhe3_reactivity_m3_s(temperature_kev: float) -> float:
    """Thermal <sigma v> for D-He3 in m^3/s."""
    return bosch_hale_reactivity_m3_s(
        temperature_kev,
        bg_sqrt_kev=BOSCH_HALE_DHE3_BG_SQRT_KEV,
        reduced_mass_energy_kev=BOSCH_HALE_DHE3_REDUCED_MASS_ENERGY_KEV,
        c1=BOSCH_HALE_DHE3_C1,
        c2=BOSCH_HALE_DHE3_C2,
        c3=BOSCH_HALE_DHE3_C3,
        c4=BOSCH_HALE_DHE3_C4,
        c5=BOSCH_HALE_DHE3_C5,
        c6=BOSCH_HALE_DHE3_C6,
        c7=BOSCH_HALE_DHE3_C7,
        t_min_kev=BOSCH_HALE_DHE3_T_MIN_KEV,
        t_max_kev=BOSCH_HALE_DHE3_T_MAX_KEV,
        reaction="D-He3",
    )


def pb11_reactivity_m3_s(temperature_kev: float) -> float:
    """Low-temperature Nevins-Swain <sigma v> for p-11B, in m^3/s.

    The non-resonant piece is equation (10) and the 148 keV resonance is
    equation (15). Energies are converted to joules and the S-factor to
    joule square-metres before those equations are evaluated, so the
    printed prefactors return m^3/s. The high-temperature Padé is absent.
    """
    temperature = float(temperature_kev)
    if (
        temperature < NEVINS_SWAIN_PB11_T_MIN_KEV
        or temperature > NEVINS_SWAIN_PB11_T_MAX_KEV
    ):
        raise ValueError(
            "p-B11 reactivity is outside the Nevins-Swain low-temperature "
            f"interval [{NEVINS_SWAIN_PB11_T_MIN_KEV}, {NEVINS_SWAIN_PB11_T_MAX_KEV}] keV: "
            f"{temperature} keV"
        )
    joule_per_kev = EV_PER_KEV * JOULE_PER_EV
    joule_per_mev = 1.0e6 * JOULE_PER_EV
    s_per_mev_barn = joule_per_mev * BARN_M2
    thermal = temperature * joule_per_kev
    gamow = NEVINS_SWAIN_PB11_EG_MEV * joule_per_mev
    resonance_energy = NEVINS_SWAIN_PB11_EL_KEV * joule_per_kev
    resonance_width = NEVINS_SWAIN_PB11_DEL_KEV * joule_per_kev
    reduced_mass = (NEVINS_SWAIN_PB11_MUC2_MEV * joule_per_mev) / SPEED_OF_LIGHT_M_S ** 2
    c0 = NEVINS_SWAIN_PB11_C0_MEV_B * s_per_mev_barn
    c1 = NEVINS_SWAIN_PB11_C1_MEV_B * s_per_mev_barn
    c2 = NEVINS_SWAIN_PB11_C2_MEV_B * s_per_mev_barn
    amplitude = NEVINS_SWAIN_PB11_AL_MEV_B * s_per_mev_barn
    gamow_peak = (gamow / 4.0) ** (1.0 / 3.0) * thermal ** (2.0 / 3.0)
    peak_width = 4.0 * math.sqrt(thermal * gamow_peak / 3.0)
    tau = 3.0 * gamow_peak / thermal
    s_effective = (
        c0 * (1.0 + 5.0 / (12.0 * tau))
        + c1 * ((gamow_peak + (35.0 / 36.0) * thermal) / joule_per_kev)
        + c2 * ((gamow_peak ** 2 + (89.0 / 36.0) * gamow_peak * thermal) / joule_per_kev ** 2)
    )
    nonresonant = (
        math.sqrt(2.0 / reduced_mass)
        * peak_width
        * s_effective
        / thermal ** 1.5
        * math.exp(-tau)
    )
    resonant = (
        math.sqrt(8.0 * math.pi * thermal / reduced_mass)
        * (joule_per_kev / thermal) ** 2
        * (amplitude / resonance_width)
        * math.exp(-(math.sqrt(gamow / resonance_energy) + resonance_energy / thermal))
    )
    return nonresonant + resonant
