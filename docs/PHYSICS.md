# Physics assumptions

SolCascade evaluates a pulse you specify. It does not search the pulse for a larger yield.

## Constants

CODATA 2022 fixes the inverse fine-structure constant used here. The vacuum permeability is computed from

    mu_0 = 2 h alpha / (e^2 c)

with the exact 2019 SI values of h, e, and c. It is not typed as 4 pi e-7.

The DT Q value is the AME2020 mass-excess difference

    Q = Delta(D) + Delta(T) - Delta(He-4) - Delta(n).

Atomic excesses already balance the two electrons. The kinetic split puts the reactants at rest and gives the neutron the alpha's share of Q, using the bare alpha mass (helium atomic mass minus two electron rest energies). Electron binding is not applied.

The Bosch-Hale DT fit is used only from 0.2 keV to 100 keV. Outside that interval the reactivity function raises.

## Coils

Each coil is a filamentary circular bundle. The field is the elliptic-integral field of a loop, multiplied by the turn count. On the axis it matches

    B_z = mu_0 I R^2 / (2 (R^2 + z^2)^(3/2)).

Self-inductance is the low-frequency round-wire formula, scaled by N^2. Mutual inductance is the Neumann coaxial-loop formula, scaled by N_i N_j. Stored energy is (1/2) I^T L I in the conductor-current basis.

Ohmic power uses the material resistivity at `conversion.conductor_temperature_K`, with the linear temperature coefficient from the material file.

## Plasma

The plasma is a uniform sphere inside the shell.

- Density scales as R^-3.
- Temperature scales as R^(3(1-gamma)). For a monatomic gas gamma = 1 + 2/3, so temperature scales as R^-2.
- `frozen_axial_flux` holds B R^2 at the vacuum on-axis value from the first schedule time. `vacuum_only` uses the coil field at each time.
- `coil_current_scale` multiplies every coil. `coil_scales` gives each coil its own series on that same time base. The first row is the frozen-flux current. A later row does not change the frozen internal field.
- Electrons and ions are both Z = 1 and, in this version, Te = Ti. Pressure is n_i k T_i + n_e k T_e.
- Fusion power is n_D n_T <sigma v> Q times the sphere volume. There is no beam-target term.
- The confinement time used in the triple product is the smaller of the schedule duration and the gyro-period random walk. One thermal gyroradius per cyclotron period gives D = kT / (4 π e B), so τ = 2 π a² e B / kT. The ion mass cancels. This time is π/4 of the Bohm time a² / (2 D_B) with D_B = kT / (16 e B). The 16 stays in that comparison. Because both times scale as B/T, n T τ still scales with density and field. The ratio to 3e21 keV s / m^3 is a benchmark quoted near 10–20 keV, not an ignition switch.

## Shell transport

The source is a point at the center, so every uncollided ray is radial. In a cell the collision probability is the difference of the exponentials of the optical depths at its faces. A kerma fraction of that collision is deposited in the cell. Uncollided leakage and the non-kerma share are not deposited. Scattered neutrons are not followed. The one-group macroscopic totals in `configs/materials.yaml` are illustrative inputs, not an evaluated cross-section library.

OpenMC, when it actually returns a heating tally, replaces the neutron part of the cell power. Layer heating is spread across that layer's cells with the analytical shape as the weight. OpenMC heating can exceed the neutron kinetic energy when reaction Q is deposited. The ledger allows that. The analytical path cannot, because its kerma fractions are at most one. Each OpenMC material is set to `heat.initial_temperature_K`, so the cross sections are the shell's starting temperature. The library must contain that temperature or a neighbor OpenMC will accept.

Alphas that stop in the wall are added as power in the innermost cell. They are not part of the neutron tally.

## Heat

The pulse is short compared with conduction through a thick shell, so the native solve is backward Euler on the shell cells over the schedule duration, with the pulse-average power held constant during that interval. Zero conductivity reproduces q dt / (rho c_p). The inner face is adiabatic. The outer face loses heat by convection to the sink. A long, well-resolved run of the same operator approaches the steady spherical solution.

The MOOSE input is that transient problem on a 1D spherical mesh.

## Stress

Internal pressure on the shell is the larger of the peak plasma pressure and the vacuum magnetic pressure at the inner wall on the midplane. The choice is recorded.

Radial and hoop stress from that pressure are the Lamé solution for a hollow sphere with zero external pressure. Tensile stress is positive, so the inner radial stress equals the negative of the internal pressure.

Thermal stress is the free hollow-sphere solution on each solid layer, using that layer's own temperature samples. Layers do not constrain each other in this estimate. A uniform temperature change produces no thermal stress. A material with zero Young's modulus is a liquid and is skipped.

Von Mises stress for equal hoop components is the absolute difference of hoop and radial stress. It is compared with the yield strength from the material file.

## Electricity

Alpha energy splits into an escaped fraction, a wall fraction, and the remainder, which stays in the plasma. The escaped part is multiplied by `alpha_direct_efficiency`. That efficiency is an assumption you set. It is rejected if it lies outside 0 to 1.

Neutron heating and wall-alpha heating are shell heat. With `shell_channel: none` they produce no electricity. With `carnot_fraction` the electrical share is that fraction of the Carnot limit between the hottest shell temperature and `cold_temperature_K`. The Carnot ceiling is always reported.

The coil supply pays ohmic dissipation and any magnetic energy that is not handed back. `coil_energy_recovery_fraction` is the share of a downward change in stored energy that returns to the supply. A closed current waveform with perfect recovery and zero resistance draws no magnetic energy. Energy still sitting in the coils at the end of the schedule is reported as a change in stored energy; it is part of the supply draw only to the extent the supply had to provide it.

Net electrical energy is gross electricity minus supply draw. It may be negative.
