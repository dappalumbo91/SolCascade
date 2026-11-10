# Production refinement

The coil draws are the four supplies locked in `scripts/apply_fsot.py`.
Fusion energy is the integral of Bosch-Hale power over the checked-in
radius schedule. Electricity is the checked-in direct converter:
escaped alphas only, `shell_channel` none. The warp-shell rest energy
is not this score.

The compressed sphere holds `2.18627e-10` kg of deuterium-tritium,
`1.31092e-10` kg of it tritium. Complete burn releases
`73778.1` J.

## Topologies

Checked-in schedule. A closed row has a coil draw above the complete-burn ceiling.

| Topology | Coil supply (J) | Fusion (J) | Burn fraction | Electric (J) | Net (J) | Ceiling / supply |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline | 1.42387e+06 | 0.497703 | 6.74595e-06 | 0.0137774 | -1.42387e+06 | 0.051815 |
| same_b_flat | 2106.22 | 0.497703 | 6.74595e-06 | 0.0137774 | -2106.2 | 35.0288 |
| same_b_optimal | 2020.44 | 0.497703 | 6.74595e-06 | 0.0137774 | -2020.42 | 36.5159 |
| same_b_helmholtz | 827.501 | 0.497703 | 6.74595e-06 | 0.0137774 | -827.487 | 89.1577 |

The baseline store stays closed on this fuel load. The flat, optimal,
and Helmholtz same-B draws sit under the ceiling. At the checked-in
yield every net is negative.

## Schedule on the Bosch-Hale curve

Initial temperature sets the adiabatic trajectory. Every temperature
on a row stays inside 0.2–100 keV.
The walk yield is peak power times `gyro_walk_time_s` at that peak.
It is a different integral from the schedule.

| T0 (keV) | Peak (keV) | Fusion (J) | Burn fraction | Peak beta | Walk time (s) | Walk yield (J) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.5 | 3.85802 | 0.497703 | 6.74595e-06 | 2.38876 | 8.81033e-07 | 0.000361779 |
| 1.296 | 10 | 11.3505 | 0.000153847 | 6.19168 | 3.39905e-07 | 0.00305344 |
| 2.592 | 20 | 45.9661 | 0.000623032 | 12.3834 | 1.69952e-07 | 0.00581869 |
| 6.48 | 50 | 104.49 | 0.00141628 | 30.9584 | 6.79809e-08 | 0.00464887 |
| 12.96 | 100 | 115.021 | 0.00155901 | 61.9168 | 3.39905e-08 | 0.0022703 |

The hottest schedule in the window releases `115.021` J
at a peak of `100` keV. Its burn fraction
is `0.00155901`. The Helmholtz draw is
`827.501` J, and the burn fraction that matches it
on this ceiling is `0.0112161`.
The schedule remains short of that draw by a factor of
`7.19435`.

Matching the Helmholtz draw from the hottest peak power would take
`0.0123892` s at that power. The gyro-walk time at that peak is
`3.39905e-08` s, and the walk yield is
`0.0022703` J. The walk yield at 100 keV is below
the walk yield at 20 keV (`0.00581869` J).
Heating shortens the walk. The confinement time does not close the
Helmholtz gap.

Peak beta on the checked-in schedule is `2.38876`.
It rises with temperature. Beta above 1 means the thermal pressure
exceeds the magnetic pressure of the frozen-flux field. The schedule
integral is the prescribed compression.

## Inventory

Mass below is more fuel at the same density, temperature, and burn
fraction. The joule target is the catalog energy of that reference unit.

- One JET pulse, `one JET deuterium-tritium pulse, 59 MJ fusion`: `5.9e+07` J.
  Checked-in burn fraction needs `0.025917` kg.
  The hottest burn fraction needs `0.000112145` kg.
- One Pu-239 kilogram, `1 kg Pu-239, complete fission`:
  `8.07233e+13` J of fission energy, before the catalog efficiency.
  Checked-in burn fraction needs `35459.5` kg of DT.
  The hottest burn fraction needs `153.435` kg of DT.

A denser sphere would raise the schedule burn fraction in proportion
to density, and it would raise beta by the same factor. The checked-in
peak is already above beta 1, so that direction moves farther from
magnetic balance.

## Tritium and DT neutrons

Each DT reaction consumes one tritium nucleus and emits one neutron.
The neutron energy in this ledger is `14048.2` keV,
`14.0482` MeV, from the two-body split of Q.
Loaded tritium scales with fuel mass. Consumed tritium scales with
the number of reactions.

| Account | Fuel (kg) | Tritium loaded (kg) | Tritium consumed (kg) | Neutrons | Neutron energy (J) |
| --- | ---: | ---: | ---: | ---: | ---: |
| checked-in schedule | 2.18627e-10 | 1.31092e-10 | 8.84342e-16 | 1.76609e+11 | 0.397504 |
| hottest schedule | 2.18627e-10 | 1.31092e-10 | 2.04375e-13 | 4.08148e+13 | 91.8646 |
| JET mass at the hottest burn fraction | 0.000112145 | 6.72438e-05 | 1.04834e-07 | 2.0936e+19 | 4.71219e+07 |
| Pu-239 kilogram at the hottest burn fraction | 153.435 | 92.0024 | 0.143433 | 2.86444e+25 | 6.44718e+13 |

The p-B11 burn fraction at the checked-in peak, one walk at
`3.85802` keV, is `5.82048e-17`.
That branch does not carry this tritium account.
