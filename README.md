# SolCascade

SolCascade is a local research code for one question: given a prescribed multi-coil compression pulse and a spherical shell, how much of a deuterium-tritium fusion output can be turned into electricity without inventing a thermal cycle that was not asked for.

The minimal pipeline does four things.

1. It builds an axisymmetric stack of circular coils and evaluates the vacuum field.
2. It steps a uniform spherical plasma through a radius and current schedule and integrates fusion power with the Bosch-Hale reactivity.
3. It deposits neutron energy in the shell, by a first-collision model or by OpenMC when that code and a cross-section library are installed.
4. It solves transient radial heat conduction, estimates elastic stress in the solid layers, and writes an energy ledger.

Shell heat is reported as heat. It is counted as electricity only if the config sets `conversion.shell_channel` to `carnot_fraction`. The default is `none`. Escaped alpha particles use a direct-converter efficiency that is an input in the YAML file. Net electrical energy is gross electricity minus the energy drawn from the coil supply. On the checked-in lab-scale pulse that net is negative. That is a result, not a failure of the ledger.

SolCascade is not a qualified nuclear design tool, and it does not adjust the schedule to chase yield.

## Layout

| Path | Role |
| --- | --- |
| `src/solcascade/constants.py` | CODATA 2022, AME2020, the Bosch-Hale coefficients, and the scored warp-shell mass. The only place those numbers are typed. |
| `src/solcascade/interfaces.py` | Flat numeric slots. An FSOT provider replaces slots by key. |
| `src/solcascade/geometry/` | Coils and the spherical shell. |
| `src/solcascade/compression/` | Staged radius and current timing. |
| `src/solcascade/plasma/` | Zero-dimensional confinement metrics. TokaMaker is optional. |
| `src/solcascade/transport/` | First-collision deposition and the OpenMC model writer. |
| `src/solcascade/mechanics/` | Radial heat conduction and thick-shell stress. |
| `src/solcascade/conversion.py` | Direct-conversion and coil-supply ledger. |
| `src/solcascade/hybrid/` | Fertile blanket, fission ledger inputs, Bateman depletion, and the r-process comparison library. |
| `src/solcascade/energy_baseline/` | Mechanism catalog and the warp-shell ranking. Phase one is the table only. |
| `src/solcascade/energy_baseline/phase2/` | Observational FSOT comparison against that table. The catalog stays read-only. |
| `src/solcascade/energy_baseline/phase3/` | Derived confinement time and the re-ranked table. Overrides live outside the catalog. |
| `src/solcascade/coupling/` | MOOSE, Cardinal, TokaMaker, Griffin, and Bison files. |
| `configs/` | Example YAML. |
| `docs/PHYSICS.md` | What each model assumes. |
| `docs/ARCHITECTURE.md` | How a later scalar module plugs in. |

## Install

Python 3.12 or newer. The pipeline uses NumPy, SciPy, and PyYAML. It does not call a network service and it does not need an account.

pip, from the repository root:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
pytest
python -m solcascade run configs/minimal_pipeline.yaml
```

conda:

```bash
conda env create -f environment.yml
conda activate solcascade
python -m pip install -e ".[dev]"
pytest
python -m solcascade run configs/minimal_pipeline.yaml
```

On Windows the same pip commands work in PowerShell. OpenMC, MOOSE, and TokaMaker are built for Linux and macOS.

## What the run writes

The work directory defaults to `runs/minimal`.

- `summary.json` and `report.md` — ledger, peak plasma numbers, which backends ran
- `plasma_history.csv` — radius, field, fusion power, triple-product ratio
- `radial_map.csv` — deposition, volumetric source, temperature, von Mises stress
- `field_grid.npz` — vacuum Br and Bz at the strongest coil current
- `openmc_model.py` — spherical-shell OpenMC model
- `heat_conduction.i` — MOOSE heat-transfer input for the same source
- `cardinal_skeleton.i` — Cardinal coupling skeleton
- `tokamaker_equilibrium.py` — TokaMaker handoff
- `griffin_interface.json` and `bison_interface.json` — tables, not executions

List every number a provider is allowed to replace:

```bash
python -m solcascade keys configs/minimal_pipeline.yaml
```

## OpenMC

The analytical deposition always runs, so the ledger does not depend on nuclear data. Set `transport.backend` to `auto` (the example) or `openmc`.

On Linux or macOS, conda-forge is the straightforward install:

```bash
conda install -c conda-forge openmc
```

PyPI has no Windows wheel. On this machine OpenMC 0.16.0 is installed in Ubuntu WSL, conda environment `solcascade` under `/opt/miniforge3`. The nuclear data is the official ENDF/B-VIII.0 HDF5 library at `/opt/nuclear-data/endfb-viii.0-hdf5/cross_sections.xml`. `scripts/openmc-wsl.sh` activates that environment and exports `OPENMC_CROSS_SECTIONS`.

```bash
wsl -d Ubuntu -- bash "/mnt/c/Users/damia/Desktop/SolCascade fsot/scripts/openmc-wsl.sh" run configs/minimal_pipeline.yaml --workdir runs/openmc-minimal
```

Continuous-energy transport needs an HDF5 cross-section library. OpenMC publishes that data without an account. Point `OPENMC_CROSS_SECTIONS` at `cross_sections.xml`. The generated model is a fixed-source problem: a point isotropic neutron at the DT energy from `constants.py`, spherical shell cells at `heat.initial_temperature_K`, and a heating tally. If the run finishes, layer heating replaces the first-collision neutron source before the heat solve. Photon transport is off unless `transport.openmc.photon_transport` is true.

Stable OpenMC runs in parallel with OpenMP (`transport.openmc.threads`) and MPI. This repository does not pass a CUDA switch. If your OpenMC executable was built with a GPU backend, set `transport.openmc.executable` to that binary. The coil-field grid can use CuPy when you install the `gpu` extra and set `field.array_backend` to `auto` or `cupy`. The default field method is the elliptic-integral formula on the CPU.

## MOOSE and Cardinal

The reported temperature is the native backward-Euler solve on the shell cells. `heat_conduction.i` is the same problem for the open MOOSE heat-transfer module: one-dimensional mesh, `coord_type = RSPHERICAL`, piecewise source, transient conduction. The inner face is adiabatic because no boundary condition is applied there. The outer face is `ConvectiveHeatFluxBC`.

Build MOOSE from the public repository. No account is required for the framework. Then:

```bash
export SOLCASCADE_MOOSE_EXE=/path/to/moose-opt
```

or set `heat.moose_executable`. `heat.execute_moose` is true in the example. If the executable is missing, the file is still written and the native solve stands.

[Cardinal](https://cardinal.cels.anl.gov/) is the open MOOSE application that wraps OpenMC. `cardinal_skeleton.i` names `OpenMCCellAverageProblem` and a sub-app pointing at the heat input. It is a skeleton: the mesh-to-cell map has to be filled for a real Cardinal build. The Python pipeline already copies the heating into the heat source, which is the data Cardinal would move.

## Griffin, Bison, and TokaMaker

Griffin and Bison are NEAMS applications. Obtaining the full apps goes through the Idaho National Laboratory distribution, which asks for a registration. SolCascade does not contact that service and does not require it. The pipeline writes the power-density and temperature tables those codes would consume. Stress in the summary comes from the Python Lamé solution plus a per-layer thermal-stress estimate. Liquid layers, entered with zero Young's modulus, carry no elastic stress.

[TokaMaker](https://github.com/OpenFUSIONToolkit/OpenFUSIONToolkit) is the Open FUSION Toolkit Grad-Shafranov solver. `pip install openfusiontoolkit` is available for some platforms, and a source build uses `OFT_ROOTPATH`. The generated script imports TokaMaker, prints the coil set, and stops before a solve unless a mesh path is set. It does not invent a plasma boundary. The equilibrium in `summary.json` is the vacuum coil field plus the zero-dimensional pressure.

## Injecting FSOT scalars

Every physical number the solvers read is copied into a slot map first. A provider implements `value(key) -> float | None`. `None` keeps the YAML value. Unknown keys raise. The checked-in stub is `examples/fsot_scalars.example.json`. It injects scalars only. It does not import an FSOT formula module.

```bash
python -m solcascade run configs/minimal_pipeline.yaml --provider examples/fsot_scalars.example.json
```

Three slots are the useful first injection:

- `plasma.ion_density_m3` and `plasma.ion_temperature_keV` set the fuel state. The example raises density from the YAML value `1.0e20` to `1.0e22` and leaves temperature at `0.5` keV.
- `plasma.thermal_energy_density_J_m3`, when it is a positive number, rewrites ion temperature and keeps density. `null` in the example means that slot is not injected.
- `coils.<name>.current_a` replaces a coil current. The example repeats the YAML value, `+40000` A on each coil, so the coil supply does not move.

`conversion.shell_channel` is not a slot. It stays in the YAML. Keys that start with `_` are comments. Any other unknown key is rejected.

Injected slots are listed under `overrides` in `summary.json`. A value taken from the YAML is marked `config` and is omitted there. A provider value is marked `mapping`. If the thermal-energy slot was injected, the ion-temperature entry ends in `->temperature`.

## Hybrid Pivot

The pure-fusion path above is unchanged and remains the comparison baseline. On the checked-in pulse the fusion yield stays about 0.498 J and the coil supply stays on the order of a megajoule, so the net stays negative. A fertile blanket does not rewrite that fusion yield or that coil supply. It adds a second term: fission heat in the blanket, converted at `conversion.thermal_efficiency`.

The fusion plasma is a 14.1 MeV neutron driver. The lithium-lead absorber stays in place and still breeds tritium. A U-238 or Th-232 zone is wrapped outside the steel. `blanket.enrichment` is the weight fraction of that fertile nuclide, not a U-235 enrichment. `blanket.density_g_cm3` is converted to kg/m3 at the config boundary. `blanket.neutron_yield_scale` multiplies the number of source neutrons seen by the blanket. It does not multiply fusion joules. Per-source multiplication M is the analytical kerma fraction plus fissions per source neutron times the 200 MeV ledger quantum, divided by the DT neutron energy. That M does not move when the yield scale moves. Fission joules do. The 200 MeV figure is the ledger quantum in `constants.py`. It is not an ENDF Q-value.

The analytical transport model is a first-collision surrogate. One collision removes the neutron. It cannot build the slowed cascade that a hybrid blanket uses for extra fissions and for breeding Pu-239 or U-233. Published energy multiplications of about 10 to 90 are literature anchors for that cascade. They are not installed as a target. OpenMC, when it runs, replaces the fission-per-source and tritium-per-source numbers with tallies. Heating from OpenMC is recorded separately and is not substituted for M. If that heating was applied to the heat equation, the analytical fission power is not added again.

`shell_channel` stays `none` on these configs, so shell heat is still not electricity. Fission electricity is `thermal_efficiency` times fission thermal energy, with the Rankine starting value 0.33 stored in the YAML.

A neutron-star merger is the natural case of a neutron flux driving fission recycling. Neutron-rich ejecta captures neutrons faster than beta decay (the r-process). Past the actinides, fission fragments are captured again and the cycle repeats. Freiburghaus, Rosswog, and Thielemann, Astrophysical Journal 525, L121 (1999), calculated r-process nucleosynthesis in merger ejecta. Korobkin, Rosswog, Arcones, and Winteler, Monthly Notices of the Royal Astronomical Society 426, 1940 (2012), is the fission-cycling calculation that later heating papers cite; the 2015 citation Korobkin et al., Monthly Notices 448, 2544, is recorded as a citation, and no yield number from that page is stored. Barnes, Kasen, Wu, and Martínez-Pinedo, Astrophysical Journal 829, 110 (2016), give a representative fission share of about 10 percent of the radioactive energy before one day, on ejecta of 0.005 solar masses at 0.2c. GW170817 / GRB 170817A (Abbott et al., Physical Review Letters 119, 161101, 2017) and the kilonova AT2017gfo are the observational anchor. `src/solcascade/hybrid/rprocess.py` stores those read numbers and can score a caller-supplied fission-energy fraction against the Barnes 10 percent. It does not simulate a merger and it does not score this blanket's M against that 10 percent.

The hybrid proposal itself is Hans Bethe, "The Fusion Hybrid," Physics Today (1979). Later blanket studies quoted here as scale anchors, not as this pulse's result, include the Chinese Academy of Sciences FFHER concept near 13 times multiplication, FDS-EM near 90, thorium blankets near 13 to 16, the FDS-series papers of Wu et al. in Fusion Engineering and Design, and the IAEA technical reports on fusion-fission hybrid systems. One-group tallies in this repository use ENDF/B-VIII.0 point values at 14.1 MeV. A point value is not those reactor multiplications.

Fission backends sit beside `transport.backend`:

- `serpent` writes `blanket_serpent.txt`. Serpent 2 is issued by VTT for non-commercial research after registration. It is not open source, it is not redistributed here, and it is not downloaded. The pipeline runs `sss2` or `serpent` only when that binary is already on PATH and `fission.execute` is true. Otherwise `executed` is false.
- `mcnp` writes `blanket_mcnp.txt`. MCNP6 is distributed by RSICC and is export controlled. It is not open source and it is not downloaded. The same PATH rule applies to `mcnp6` or `mcnp`.
- `griffin` keeps the existing table handoff in `griffin_interface.json`. A Griffin binary is not launched.
- `none` writes no fission input.

The Bateman tracker still runs when `fission.pulses` is positive. It follows U-238 to Pu-239 and Th-232 to U-233 across the chosen interpulse gap. That gap is a schedule input. On `hybrid_u238.yaml` it is one day, long enough for the beta decays and not a measured repetition rate. Decay is not applied inside the 3.5 ms pulse. The inventory is `isotopic_inventory.json`. At 14.1 MeV the capture cross section is small, so the surrogate breeds almost no Pu-239 or U-233. Breeding in a real blanket is carried by slowed neutrons.

```bash
py -3.14 -m solcascade run configs/hybrid_u238.yaml --workdir runs/hybrid/depletion-u238
py -3.14 -m solcascade run configs/hybrid_th232.yaml --workdir runs/hybrid/th232-analytical
py -3.14 -m solcascade run configs/hybrid_u238.yaml --provider examples/hybrid_neutron_yield.json --workdir runs/hybrid/fsot-yield
```

OpenMC for the same blankets, from Windows, uses the WSL environment and the ENDF/B-VIII.0 library:

```bash
wsl -d Ubuntu -- bash "/mnt/c/Users/damia/Desktop/SolCascade fsot/scripts/openmc-wsl.sh" run configs/hybrid_u238.yaml --workdir runs/hybrid/openmc-u238
wsl -d Ubuntu -- bash "/mnt/c/Users/damia/Desktop/SolCascade fsot/scripts/openmc-wsl.sh" run configs/hybrid_th232.yaml --workdir runs/hybrid/openmc-th232
```

`examples/hybrid_neutron_yield.json` sets `blanket.neutron_yield_scale` to 2. That value is a slot check. It is not an FSOT yield.

## Energy Baseline

The catalog is additive. It does not change the fusion ledger, the hybrid blanket, or any solver. This phase writes the comparison table. It does not derive an FSOT law, it does not edit `vendor/fsot_compute.py`, and it does not mark any mechanism as mapped.

SolCascade did not previously store a warp-shell joule requirement. `constants.py` now holds `WARP_SHELL_MASS_KG` for the shell already scored in the local warp ledger: compactness 1/2, outer radius 110 m, mass `3.703e28` kg. `warp_shell_energy_j()` returns that mass times `c` squared. The ranker calls that function. It does not retype the joule product.

Each mechanism is one YAML file in `src/solcascade/energy_baseline/data/`. The file records the governing equation and its source, the efficiency ceiling, the best demonstrated efficiency with a year, the energy density and the reference unit it belongs to, the power density, the input energy required to initiate and sustain the same unit, whether the process is a chain or needs a continuous drive, a TRL from 1 to 9 with the main engineering barrier, and an FSOT status. Unknown numbers are null. A fringe claim also records the source, the claimed mechanism, the measured watts when a number was actually read, and the replication status. The claim is kept and flagged. `sources.md` lists every DOI, arXiv id, and institutional page, grouped by mechanism id.

Every file in this catalog has `fsot_status: open_gap`, with `fsot_scalar` and `fsot_mapping` null. The note may name the nearest domain in the Lean engine. `domain_scalar` on that domain is a dimensionless fold. It does not set an efficiency, an energy density, or an input energy. No mapping is emitted in this phase.

Delivered joules for one reference unit are `energy_density_j * demonstrated_efficiency` when both numbers are present. If the stored density is already the reported output, the file sets the efficiency to 1 and says so in the notes, so the product is not applied twice. Coverage is delivered joules divided by `warp_shell_energy_j()`. Output over input is delivered joules divided by `input_energy_j` when that input is positive and covers the same unit and the same time window. A zero or negative input leaves the ratio blank. Rows with a blank coverage sort last. Equal coverage sorts by mechanism id. A fringe flag does not move a row.

The reference units are not interchangeable. One kilogram of U-235, one square-meter-year of sunlight, and one fusion pulse are each scored on the joules stated for that unit. The sorted table is a comparison of those stated units against the shell. It is not a choice of which technology to build.

To add a mechanism, copy one YAML file, keep every key the ranker requires, put null in for an unknown number, and cite the source in the file and in `sources.md`. Leave `fsot_status` as `open_gap`. Do not invent a scalar or a mapping. Phase two is the later step that asks which mechanism has a free parameter the framework can collapse, and it starts from this table.

From the project directory:

```bash
py -3.14 -m solcascade.energy_baseline.rank
```

The command prints the table and writes `runs/energy_baseline/ranking.md`. The `runs/` directory is gitignored.

## Energy Baseline Phase Two

Phase two reads the catalog and scores FSOT against it. It does not edit a mechanism file, a solver, a hybrid config, or `vendor/fsot_compute.py`. A row that looks like a mapping lives only in the phase-two markdown. Phase one stays `open_gap`.

The engine load matches `scripts/apply_fsot.py`: the file at `FSOT-2.1-Lean/vendor/fsot_compute.py`, pin prefix `aeb2adad6e80f487`, imported by path. The joule target is `warp_shell_energy_j()`. It is not retyped.

For each mechanism the script names the parameter that moves the scored coverage most, then one candidate scalar. Fission uses the burnup fraction and the nuclear scalar. Fusion uses confinement time and the nuclear scalar. MHD and the electrostatic generator use the electromagnetic scalar. Other rows use the domain already named in the phase-one note. `K`, `phi^-4`, and `gamma_FSOT` (the condensed-matter scalar divided by `K`) are captured with the domain scalars. They are not substituted unless a reading names them. No product of seeds is built to chase a hit.

Coverage is the phase-one value, delivered joules over the shell. A dimensionless parameter that multiplies coverage linearly is in range only when the scalar is at least ten times the demonstrated value and still at or under the physical cap. Coverage if fixed is then catalog coverage times scalar divided by that demonstrated value. The note reports the orders that remain. A scalar outside the window, a cap that cannot admit a factor of ten, or a dimensionful parameter (seconds, meters, electron-volts, a speed) is an open gap. The engine has no zero-parameter conversion from a dimensionless fold into those units, and this pass does not invent one.

To add a scalar to the comparison, add its domain name to `CAPTURED_DOMAINS` in `map_fsot.py`, or extend `capture_scalars` for another stable attribute the engine already prints. Point `ParameterReading.candidate` at that name for the mechanism that should use it. Leave the engine file and the catalog YAML alone.

From the project directory:

```bash
py -3.14 -m solcascade.energy_baseline.phase2.map_fsot
py -3.14 -m solcascade.energy_baseline.phase2.gaps
```

The first command writes `runs/energy_baseline/phase2_mapping.md`. The second writes `runs/energy_baseline/gaps.md`. Computable upside is the coverage gained by taking a dimensionless parameter to its cap. That list comes first. Dimensionful gaps follow, ordered by the coverage already scored, because a missing unit is not given a made-up closure. Unscored rows sort last by id.

## Phase Three — Derived Values

Phase three installs dimensionful quantities this project already derives and substitutes them where the derivation matches the dominant parameter. Mechanism YAML files stay as phase one wrote them. `vendor/fsot_compute.py` stays unread-for-writing. The overrides live in `src/solcascade/energy_baseline/phase3/overrides.yaml`.

The confinement time is `gyro_walk_time_s` in `src/solcascade/plasma/metrics.py`: `D = kT / (4 π e B)`, `τ = a² / (2 D)`. The radius, temperature, and field are the peak state recorded in `runs/fsot-apply/report.json` and stored in the overrides file. The function runs when the ranker loads that file. `gamma_enters_the_time` is false. The confined DT yield is `E = n_D n_T <σv>(T) Q V τ` at that peak, with density and temperature taken from the checked-in pulse by the same compression `evolve` uses. The three DT rows take that yield as delivered joules. A scored catalog yield then moves by the factor yield/delivered.

The same peak supplies the remaining units. The confinement radius replaces hydraulic head, tidal range, and the Casimir plate gap. Coverage scales linearly with head and with the square of the range. Casimir cavity energy scales as the inverse cube of the gap, and the coverage stays 0 because the demonstrated efficiency is 0. The ion thermal speed inside the walk, `v = sqrt(kT / m)`, replaces the wind speed, and coverage scales as speed cubed. The photovoltaic bandgap is that `kT` in electron-volts. The factor is the fraction of 6000 K blackbody energy above it, the solar temperature in Shockley and Queisser, J. Appl. Phys. 32, 510 (1961). D-He3 uses `dhe3_reactivity_m3_s` on the Bosch-Hale interval 0.5–190 keV. The burn fraction `n_He3 <σv> τ` multiplies the catalog joules per reference kilogram. p-B11 uses `pb11_reactivity_m3_s`, the Nevins and Swain low-temperature branch, Nucl. Fusion 40, 865 (2000), on 1 keV to 70 keV. The burn fraction `n_B <σv> τ` multiplies the catalog joules per reference kilogram. The high-temperature branch starts at 50 keV and is not evaluated at this peak.

`FSOT-2.1-Lean/scripts/tier_gap_fill_lib.py` `_fsot_scaled` and the magnetic-confinement panel in `scripts/tier71_fusion_lab_lib.py` were read. They rescale a measured triple product by `(1 + |S| * factor)`. The walk clock is the formula this pass uses. `FSOT-2.1-Neural` is absent on this machine. `FSOT-Genetics/formulas/FSOT_PROTEIN_DERIVATIONS.md` and the engine chemistry bond table were read and left off the ranking.

From the project directory:

```bash
py -3.14 -m solcascade.energy_baseline.phase3.rank_derived
```

The command writes `runs/energy_baseline/phase3_ranking.md` and `runs/energy_baseline/phase3_derivations.md`, and refreshes `runs/energy_baseline/gaps.md` with the phase-two table plus the phase-three section.

## Production standing

The phase-one, phase-two, and phase-three tables still divide delivered joules by `warp_shell_energy_j()`. That figure is the warp-shell rest energy. Useful production is the net joules of a device set against the energy that safe, operating plants already deliver on the reference units in the catalog.

The checked-in schedule deposits `0.4977033387900156` J of fusion and draws `1423874.990828079` J from the baseline coil supply. Fusion over that supply is `3.495e-7`. Net electrical energy on that pulse is about `-1.423875e6` J.

The compressed sphere holds `2.186e-10` kg of deuterium-tritium. Burning every pair in it releases `73778` J. That ceiling is 19.3 times short of the baseline coil draw, so the baseline store remains a net consumer at every burn fraction of this sphere. The schedule integral burns `6.746e-6` of the ceiling. The walk-time yield, `3.618e-4` J, burns `4.904e-9` of it. Those two integrals are different.

The same fusion figure against the four supplies locked in `scripts/apply_fsot.py`:

| Topology | Coil supply (J) | Fusion / supply | Ceiling / supply |
| --- | ---: | ---: | ---: |
| baseline | 1423874.990828079 | 3.495e-7 | 0.05182 |
| same_b_flat | 2106.216 | 2.363e-4 | 35.03 |
| same_b_optimal | 2020.4381637452752 | 2.463e-4 | 36.52 |
| same_b_helmholtz | 827.5009597961755 | 6.015e-4 | 89.16 |

The Helmholtz, flat, and optimal same-B draws sit under the fuel ceiling. Matching the Helmholtz draw takes the schedule yield up by a factor of 1663, to a burn fraction of `0.01122`. A zero-resistivity winding removes the ohmic channel and leaves the stored field.

Bosch-Hale reactivity at the confinement-check peak, 3.858 keV, is `5.194e-24` m³/s. At 10 keV it is 21.88 times that value. At 20 keV it is 83.38 times that value. Those factors belong to the reactivity. The ledger is recomputed from the schedule. Temperature on this curve supplies part of the factor of 1663. Confinement time supplies the rest.

One catalog JET pulse is `5.9e7` J, 800 times this sphere's complete-burn ceiling. Complete fission of one kilogram of Pu-239 in the catalog is `8.072e13` J, and the delivered row at efficiency 0.40 is `3.229e13` J. Complete DT burn is `3.375e14` J per kilogram of fuel, 4.18 times that Pu-239 fission energy density. A JET-class pulse or a fission kilogram takes a larger reacting inventory than the `2.186e-10` kg in this sphere.

Safe, operating catalog rows, on the units the files state: fission-thermal delivers `2.709e13` J per kilogram of U-235 at efficiency 0.33 (TRL 9); wind delivers `1.578e11` J per square meter of rotor over 25 years (TRL 9); hydroelectric delivers `882.9` J per kilogram of water through a 100 m head (TRL 9); photovoltaic delivers `4.325e10` J per square meter over 30 years (TRL 9). The phase-three wind, hydro, tidal, photovoltaic, and Casimir factors place the confinement radius, the ion thermal speed, and `kT` into those catalog laws. Each substitution uses a quantity the pulse already carries. D-He3 at this peak yields `453.4` J per reference kilogram. p-B11 yields `0.004072` J per reference kilogram. Hybrid fission electricity on this pulse is `0.0347` J on the analytical U-238 blanket and `0.558` J on the OpenMC U-238 case. The pulse net stays on the fusion and coil figures above.

Refinement order:

1. Score net joules, burn fraction, and coil draw on the four topologies in the table. Keep `warp_shell_energy_j()` as the warp-shell score.
2. Treat the baseline 1.424 MJ store as closed for net production on this fuel load. The Helmholtz draw and the other same-B supplies are the draws under the ceiling.
3. Raise the schedule burn fraction from `6.746e-6` toward `0.0112` on the installed Bosch-Hale curve and the confinement time, and recompute the ledger.
4. Grow the reacting inventory when the goal is a JET-class pulse or a fission kilogram.
5. Write the tritium inventory and the 14 MeV neutron account before a larger DT pulse is treated as a safe plant. The p-B11 burn fraction at this peak is `5.820e-17`.

## Tests

```bash
pytest
```

The tests cover the on-axis coil field, a Biot-Savart cross-check, inductance, Bosch-Hale's published temperature window, compression scalings, neutron-fraction balance, the adiabatic heat limit, zero thermal stress at uniform temperature, the Lamé inner-wall balance, the electrical ledger, and one full run of the example config. They do not need OpenMC, MOOSE, a GPU, or a network.

## License

MIT. See `LICENSE`.
