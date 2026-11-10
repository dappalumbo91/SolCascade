# Failure trace

Net electrical energy is electric output minus coil supply. The chain
below is that difference on the checked-in pulse: pressure balance,
the burn integral, the converter, and the stored field. The warp-shell
rest energy is outside this score.

## Checked-in product

The operating converter returns `0.027682` J of electricity
per fusion joule. That product is the alpha share of Q times the YAML
escape fraction times the YAML direct efficiency.
The alpha share is `0.201323`. The neutron share is
`0.798677`. `shell_channel` is `none`,
so the neutron share is heat on this ledger.

Checked-in fusion is `0.497703` J. Electricity is
`0.0137774` J. The baseline supply is `1.42387e+06` J,
and the net is `-1.42387e+06` J.
The Helmholtz ledger draw is `827.501` J, and the
net against that draw is `-827.487` J.

## Pressure

Peak beta on the checked-in schedule is `2.38876`.
Beta is thermal pressure divided by magnetic pressure. Above 1, the
thermal pressure exceeds the frozen-flux field that the schedule prescribes.

Cooling at fixed density divides the initial temperature by that beta.
The evolved initial temperature is `0.209313` keV
and the evolved peak beta is `1`. The evolved
schedule releases `0.00953724` J, electric
`0.00026401` J. The partner-ion e-folding time at that
peak is `9141.59` s. The gyro-walk there is
`2.10458e-06` s.

The density-limited rows keep each temperature trajectory and divide
the density by the raw peak beta. Fusion then scales as the square of
that factor because the integrator is undepleted, and the ceiling
scales as the density. Peak beta on each of those rows is 1.

| Peak (keV) | Raw beta | Density scale | Fusion (J) | Electric (J) | Walk yield (J) | E-fold (s) | Complete electric (J) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3.85802 | 2.38876 | 0.418627 | 0.0872216 | 0.00241447 | 6.34012e-05 | 429.19 | 854.971 |
| 10 | 6.19168 | 0.161507 | 0.296073 | 0.0081959 | 7.96476e-05 | 50.8515 | 329.85 |
| 20 | 12.3834 | 0.0807536 | 0.299752 | 0.00829772 | 3.79445e-05 | 26.685 | 164.925 |
| 50 | 30.9584 | 0.0323014 | 0.109023 | 0.00301798 | 4.85055e-06 | 33.3999 | 65.97 |
| 100 | 61.9168 | 0.0161507 | 0.0300027 | 0.000830535 | 5.92199e-07 | 68.3926 | 32.985 |

The strongest density-limited schedule releases `0.299752` J
at `20` keV, electric `0.00829772` J.
Its walk yield is `3.79445e-05` J. Its e-folding time is
`26.685` s beside a walk of `1.69952e-07` s and a
schedule duration of `0.0035` s.
The schedule burns a fraction of the ceiling.

Complete burn of that reduced density would deliver `164.925` J
at the checked-in converter. The schedule on the same row is short of
that ceiling by `19875.9`.
The Helmholtz practical supply is `8517.06` J, and that
ceiling is short of that supply by `51.642`.

## Shell heat

The built shell stores `279035` J per kelvin.
Checked-in neutron energy is `0.397504` J and
would raise that shell by `1.42457e-06` K if it all
stopped in the shell. Complete burn of the original sphere would raise
it by `0.211174` K. The configured shell temperature
is the heat boundary the converter reads. The rises above are the pulse.

## Converter bound

The installed `convert` function, with every alpha collected and every
neutron joule passed through the Carnot factor of the configured shell
and the configured cold sink, returns `0.600662` J per fusion
joule. The operating product is `0.027682` J. The ratio
of those two products is `21.6987`.

On the strongest density-limited schedule the bound times the integral is
`0.180049` J. On the widest density-limited
ceiling the bound times complete burn is `18551.7` J.
Widest converter-cap complete burn versus the practical supply: above.
The checked-in ledger keeps the operating product. The bound is the
parameter cap of `convert()`. The fusion column is the schedule integral.

## Stored field

The confining on-axis field at the first current knot is `0.21639` T.
The four-coil store is `62850.9` J at that knot
and `2.79337e+06` J at the peak scale. The baseline supply
on the locked ledger is the draw of that ramp.

The Helmholtz pair at the same on-axis field stores `15379.1` J.
Its locked ledger draw is `827.501` J. At recovery
`0.5` the practical supply is the ledger plus the
unrecovered store, `8517.06` J. A zero-resistivity winding
sets the ohmic ledger to zero and leaves the unrecovered store,
`7689.56` J.

Original-density complete burn at the checked-in converter is `2042.32` J.
Original-density complete burn versus the Helmholtz ledger: above.
Original-density complete burn versus the unrecovered store: below.
The widest density-limited ceiling is the `3.85802` keV row,
`854.971` J at the checked-in converter. Its schedule
electricity is `0.00241447` J and its e-folding time is
`429.19` s.
Widest complete burn versus the Helmholtz ledger: above.
Widest complete burn versus the practical supply: below.
Widest complete burn versus the unrecovered store: below.

## Field scaling

At beta 1 the density scales with the square of B. Undepleted schedule
fusion scales with the square of density, so with B to the fourth,
while magnetic energy and ohmic energy on a fixed winding and a fixed
pulse length scale with B squared. The undepleted ratio therefore
improves with B squared while the burn fraction stays below 1.
At complete burn the ceiling scales with density, so with B squared,
and electric energy over the store stays at the ratio for that temperature.
The highest complete-burn ratio in this table is the `3.85802` keV row, `0.100383`.
Hotter rows sit lower because beta is higher and the beta = 1 ceiling
is smaller. Raising the field while holding density and temperature
multiplies the store by the present beta and leaves the fusion integral
where it is.

## FSOT record

Pin `aeb2adad6e80f487` matches the loaded engine prefix
`aeb2adad6e80f487`. The scalars below are dimensionless.
The FSOT scalars are recorded and stay out of the joule product.

Engine valve: `2.620867e-17`.
Plasma beta on the checked-in peak: `2.38876`.
Engine valve and plasma beta are different quantities.

| Scalar | Value |
| --- | ---: |
| Particle_Physics | 0.950197 |
| High_Energy_Physics | 0.88733 |
| Electromagnetism | 0.955729 |
| Nuclear_Physics | 0.936388 |
| Thermodynamics | 0.936388 |
| Condensed_Matter | 0.941271 |
| Chemistry | 0.955289 |
| Quantum_Mechanics | 0.950197 |
| Optics | 0.953659 |
| Acoustics | 0.953659 |
| Fluid_Dynamics | -0.617907 |
| Atmospheric_Physics | -0.597972 |
| Oceanography | -0.597972 |
| Geophysics | -0.578458 |
| K | 0.420109 |
| phi^-4 | 0.145898 |
| gamma_FSOT | 2.24054 |

K is the scalar rest unit. gamma_FSOT is the orifice witness.
phi^-4 is the filament-fraction ceiling recorded on the coil scan.
The enthalpy identity in the FSOT physics notes uses a power that is
the logarithm of one chemistry leaf. This trace assigns no power to
the fusion yield.

## Standing

The confined schedule remains a net consumer on every beta = 1 row.
On the widest ceiling the burn factor is `354103`, the store
factor is `9.96181`, and the converter factor is
`21.6987`. The burn factor is complete-burn electricity
divided by the schedule electricity. The store factor is the practical
supply divided by that complete-burn electricity. The converter factor
is the Carnot-and-alpha cap divided by the operating product. Holding
that peak for one e-folding takes `4.87144e+08` walks. The walk
is the confinement time already installed. The FSOT record does not
move reactivity, density, field, or the converter product.
