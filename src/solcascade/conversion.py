"""Energy ledger for one pulse.

Charged fusion products that escape can be collected by a direct converter.
The converter efficiency is an input, not a derived plasma result, and it
cannot exceed one.

Neutron energy deposited in the shell, and alphas that stop in the first
wall, are heat. They become electricity only when ``shell_channel`` is
``carnot_fraction``. That channel is capped by the Carnot limit between the
hottest shell temperature and the configured cold temperature. The default
channel is ``none``: the heat is reported, and it is not counted as
electricity.

Coil ohmic dissipation is always a cost. Magnetic energy that is returned
to the supply is credited by ``coil_energy_recovery_fraction``. Energy that
is still stored in the coils at the end of the schedule is reported on its
own and is not counted as destroyed.

Fission heat is a separate term. It is zero unless a fertile blanket is
configured. The coil supply is not changed by that term. ``thermal_efficiency``
is a hardware factor, the same kind of input as the coil recovery fraction.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SupplyLedger:
    ohmic_j: float
    magnetic_invested_j: float
    magnetic_returned_j: float
    magnetic_unrecovered_j: float
    magnetic_energy_change_j: float
    drawn_j: float


@dataclass(frozen=True)
class ConversionLedger:
    fusion_j: float
    alpha_j: float
    neutron_j: float
    alpha_escaped_j: float
    alpha_wall_j: float
    alpha_retained_j: float
    neutron_deposited_j: float
    neutron_not_deposited_j: float
    shell_heat_j: float
    electric_alpha_j: float
    electric_shell_j: float
    electric_gross_j: float
    thermal_ceiling_j: float
    carnot_efficiency: float
    shell_channel: str
    supply: SupplyLedger
    net_electrical_j: float
    fusion_to_electric_efficiency: float | None
    fission_j: float = 0.0
    electric_fission_j: float = 0.0
    thermal_efficiency: float = 0.0


def supply_energy(
    time_s: np.ndarray,
    magnetic_energy_j: np.ndarray,
    ohmic_power_w: np.ndarray,
    recovery_fraction: float,
) -> SupplyLedger:
    if not 0.0 <= recovery_fraction <= 1.0:
        raise ValueError("coil energy recovery fraction must be in [0, 1]")
    invested = 0.0
    returned = 0.0
    for index in range(time_s.size - 1):
        change = float(magnetic_energy_j[index + 1] - magnetic_energy_j[index])
        if change >= 0.0:
            invested += change
        else:
            returned += -change
    ohmic = float(np.trapezoid(ohmic_power_w, time_s))
    unrecovered = returned * (1.0 - recovery_fraction)
    change = float(magnetic_energy_j[-1] - magnetic_energy_j[0])
    drawn = invested - recovery_fraction * returned + ohmic
    return SupplyLedger(
        ohmic_j=ohmic,
        magnetic_invested_j=invested,
        magnetic_returned_j=returned,
        magnetic_unrecovered_j=unrecovered,
        magnetic_energy_change_j=change,
        drawn_j=drawn,
    )


def carnot_efficiency(hot_k: float, cold_k: float) -> float:
    if cold_k <= 0.0 or hot_k <= 0.0:
        raise ValueError("Carnot temperatures must be positive")
    if hot_k <= cold_k:
        return 0.0
    return 1.0 - cold_k / hot_k


def convert(
    fusion_j: float,
    alpha_fraction: float,
    neutron_deposited_fraction: float,
    alpha_escape_fraction: float,
    alpha_wall_fraction: float,
    alpha_direct_efficiency: float,
    shell_channel: str,
    carnot_fraction: float,
    hot_k: float,
    cold_k: float,
    supply: SupplyLedger,
    fission_j: float = 0.0,
    thermal_efficiency: float = 0.0,
) -> ConversionLedger:
    if alpha_escape_fraction < 0.0 or alpha_wall_fraction < 0.0:
        raise ValueError("alpha fractions must be >= 0")
    if alpha_escape_fraction + alpha_wall_fraction > 1.0 + 1.0e-12:
        raise ValueError("alpha escape and wall fractions sum to more than 1")
    if not 0.0 <= alpha_direct_efficiency <= 1.0:
        raise ValueError("alpha direct efficiency must be in [0, 1]")
    if not 0.0 <= carnot_fraction <= 1.0:
        raise ValueError("carnot fraction must be in [0, 1]")
    if shell_channel not in {"none", "carnot_fraction"}:
        raise ValueError(f"unknown shell channel {shell_channel!r}")
    if neutron_deposited_fraction < 0.0:
        raise ValueError("neutron deposited fraction must be >= 0")
    if fission_j < 0.0:
        raise ValueError("fission energy must be >= 0")
    if not 0.0 <= thermal_efficiency <= 1.0:
        raise ValueError("thermal efficiency must be in [0, 1]")
    alpha_j = fusion_j * alpha_fraction
    neutron_j = fusion_j - alpha_j
    escaped = alpha_j * alpha_escape_fraction
    wall = alpha_j * alpha_wall_fraction
    retained = alpha_j - escaped - wall
    # The fraction is heating divided by neutron kinetic energy. OpenMC
    # heating can exceed the kinetic energy when reaction Q is deposited.
    deposited = neutron_j * neutron_deposited_fraction
    not_deposited = max(neutron_j - deposited, 0.0)
    heat = deposited + wall
    carnot = carnot_efficiency(hot_k, cold_k)
    ceiling = carnot * heat
    if shell_channel == "none":
        electric_shell = 0.0
    else:
        electric_shell = carnot_fraction * ceiling
    electric_alpha = alpha_direct_efficiency * escaped
    electric_fission = thermal_efficiency * fission_j
    gross = electric_alpha + electric_shell + electric_fission
    net = gross - supply.drawn_j
    efficiency = None if fusion_j <= 0.0 else net / fusion_j
    return ConversionLedger(
        fusion_j=fusion_j,
        alpha_j=alpha_j,
        neutron_j=neutron_j,
        alpha_escaped_j=escaped,
        alpha_wall_j=wall,
        alpha_retained_j=retained,
        neutron_deposited_j=deposited,
        neutron_not_deposited_j=not_deposited,
        shell_heat_j=heat,
        electric_alpha_j=electric_alpha,
        electric_shell_j=electric_shell,
        electric_gross_j=gross,
        thermal_ceiling_j=ceiling,
        carnot_efficiency=carnot,
        shell_channel=shell_channel,
        supply=supply,
        net_electrical_j=net,
        fusion_to_electric_efficiency=efficiency,
        fission_j=fission_j,
        electric_fission_j=electric_fission,
        thermal_efficiency=thermal_efficiency,
    )
