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

Bosch-Hale reactivity at the confinement-check peak, 3.858 keV, is `5.194e-24` m³/s. At 10 keV it is 21.88 times that value. At 20 keV it is 83.38 times that value. Those factors belong to the reactivity. Production refinement recomputes the schedule integral on this curve and sets the gyro-walk time beside it.

One catalog JET pulse is `5.9e7` J, 800 times this sphere's complete-burn ceiling. Complete fission of one kilogram of Pu-239 in the catalog is `8.072e13` J, and the delivered row at efficiency 0.40 is `3.229e13` J. Complete DT burn is `3.375e14` J per kilogram of fuel, 4.18 times that Pu-239 fission energy density. A JET-class pulse or a fission kilogram takes a larger reacting inventory than the `2.186e-10` kg in this sphere.

Safe, operating catalog rows, on the units the files state: fission-thermal delivers `2.709e13` J per kilogram of U-235 at efficiency 0.33 (TRL 9); wind delivers `1.578e11` J per square meter of rotor over 25 years (TRL 9); hydroelectric delivers `882.9` J per kilogram of water through a 100 m head (TRL 9); photovoltaic delivers `4.325e10` J per square meter over 30 years (TRL 9). The phase-three wind, hydro, tidal, photovoltaic, and Casimir factors place the confinement radius, the ion thermal speed, and `kT` into those catalog laws. Each substitution uses a quantity the pulse already carries. D-He3 at this peak yields `453.4` J per reference kilogram. p-B11 yields `0.004072` J per reference kilogram. Hybrid fission electricity on this pulse is `0.0347` J on the analytical U-238 blanket and `0.558` J on the OpenMC U-238 case. The pulse net stays on the fusion and coil figures above.

Refinement order:

1. Score net joules, burn fraction, and coil draw on the four topologies in the table. Keep `warp_shell_energy_j()` as the warp-shell score.
2. Treat the baseline 1.424 MJ store as closed for net production on this fuel load. The Helmholtz draw and the other same-B supplies are the draws under the ceiling.
3. Raise the schedule burn fraction from `6.746e-6` toward `0.0112` on the installed Bosch-Hale curve and the confinement time, and recompute the ledger.
4. Grow the reacting inventory when the goal is a JET-class pulse or a fission kilogram.
5. Write the tritium inventory and the 14 MeV neutron account before a larger DT pulse is treated as a safe plant. The p-B11 burn fraction at this peak is `5.820e-17`.

## Production refinement

From the project directory:

```bash
py -3.14 -m solcascade.production.refine
```

The command writes `runs/production/refinement.md`. The four coil draws stay the supplies locked in `scripts/apply_fsot.py`. Fusion is the Bosch-Hale integral over the checked-in radius schedule. Electricity is escaped alphas at the checked-in converter settings, with `shell_channel` none. `runs/` stays gitignored.

On that schedule the baseline draw remains above the `73778.1` J ceiling, so that store stays closed. The same-B draws sit under the ceiling. Every checked-in net is negative. Direct electric output on the checked-in yield is `0.0137774` J.

Raising the initial temperature along the adiabatic trajectory, with every point kept inside 0.2–100 keV, moves the schedule burn fraction from `6.74595e-06` to `0.00155901` at a 100 keV peak. That yield is `115.021` J. The Helmholtz draw is `827.501` J, and the burn fraction that matches it on this ceiling is `0.0112161`. The schedule remains short of that draw by a factor of `7.19435`.

The gyro-walk yield at that 100 keV peak is `0.0022703` J, below the walk yield at 20 keV (`0.00581869` J). Heating shortens the walk. The confinement time does not close the Helmholtz gap. Peak beta on the checked-in schedule is `2.38876` and rises with temperature. The schedule integral is the prescribed compression.

At the 100 keV burn fraction, one catalog JET pulse takes `0.000112145` kg of DT and one kilogram of complete Pu-239 fission takes `153.435` kg of DT. The checked-in burn fraction would take `0.025917` kg and `35459.5` kg for those same targets. A denser sphere would raise beta, which is already above 1.

The DT neutron in this ledger carries `14.0482` MeV. The checked-in pulse loads `1.31092e-10` kg of tritium and consumes `8.84342e-16` kg, in `1.76609e+11` neutrons. The JET-mass account at the hottest burn fraction loads `6.72438e-05` kg of tritium. The p-B11 burn fraction at the checked-in peak is `5.82048e-17`.

## Failure trace

From the project directory:

```bash
py -3.14 -m solcascade.production.trace
```

The command writes `runs/production/failure_trace.md`. It loads pin `aeb2adad6e80f487` through the phase-two engine loader and records the domain scalars, K, phi^-4, gamma_FSOT, and the engine valve `2.620867e-17`. Those numbers stay out of the joule product. Plasma beta on the checked-in peak is `2.38876`. The valve and that pressure ratio are different quantities.

The operating converter returns `0.027682` J of electricity per fusion joule: the alpha share `0.201323` times the YAML escape fraction times the YAML direct efficiency. The neutron share is `0.798677`, and `shell_channel` is none. Checked-in electricity is `0.0137774` J against the baseline supply `1.42387e+06` J.

The built shell stores `279035` J/K. The checked-in neutron energy would raise it by `1.42457e-06` K. Complete burn of this sphere would raise it by `0.211174` K. The configured 600 K shell is the heat boundary. The pulse supplies the rises above.

Cooling at fixed density until peak beta is 1 sets the initial temperature to `0.209313` keV. The evolved schedule releases `0.00953724` J. The partner-ion e-folding time at that peak is `9141.59` s beside a gyro-walk of `2.10458e-06` s.

Dividing the density by the raw peak beta, so each hotter trajectory also sits at beta 1, scales the undepleted integral by beta squared. The strongest of those schedules releases `0.299752` J at 20 keV, electric `0.00829772` J. The widest ceiling in that table is the 3.858 keV row: complete-burn electricity `854.971` J, schedule electricity `0.00241447` J, e-folding time `429.19` s. That complete burn sits above the Helmholtz ledger draw `827.501` J and below the practical supply `8517.06` J. The practical supply is the locked ledger plus half of the `15379.1` J Helmholtz store. A zero-resistivity winding leaves the unrecovered store, `7689.56` J. Original-density complete burn at the checked-in converter is `2042.32` J, below that store.

`convert()` with every alpha collected and every neutron joule through the Carnot factor of the configured shell returns `0.600662` J per fusion joule, `21.6987` times the operating product. That cap times complete burn of the widest ceiling is `18551.7` J and sits above the practical supply. The schedule integral on the same row does not. The checked-in ledger keeps the operating product.

At beta 1, density and the complete-burn ceiling scale as B squared, and so does the coil store. The highest complete-burn ratio in the table, electricity over the practical supply, stays `0.100383` at every field that keeps this temperature and beta at 1. On that widest ceiling the burn factor is `354103`, the store factor is `9.96181`, and one e-folding takes `4.87144e+08` walks. The confined schedule remains a net consumer on every beta = 1 row.

## Dominant cost

From the project directory:

```bash
py -3.14 -m solcascade.production.cost
```

The command writes `runs/production/dominant_cost.md`. The baseline coil store at full current is `2.79337e6` J, above the `73778.1` J fuel ceiling. That is the largest energy in the ranking.

The Helmholtz coil store at the same on-axis field is `15379.1` J. The magnetic energy inside the checked-in plasma is `27.0976` J, so the coils hold `567.545` times the field energy in the fuel. At recovery 0.5 the unrecovered store is `7689.56` J. The ohmic ledger is `827.501` J. Beta = 1 complete-burn electricity is `854.971` J: above the ohmic ledger, below the unrecovered store. The schedule electricity on that load is `0.00241447` J.

The shell inner radius is 0.20 m. A beta = 1 sphere scaled by `2.15169`, initial radius `0.107584` m, meets the practical supply of `8517.06` J on a complete burn. The Helmholtz midplane field at that edge is `0.996268` of the on-axis field. The evolved schedule returns `0.0240525` J of electricity. The gyro-walk returns `8.09453e-05` J.

The bore-limited sphere uses 0.95 of the shell inner radius as its initial radius, the same compression ratio, and the same beta = 1 density. Evolved peak beta is 1. Complete fusion is `1.69475e6` J and complete-burn electricity is `46914` J, above the practical supply and below one catalog JET pulse. Schedule electricity is `0.132487` J. Walk electricity is `0.00139063` J in `1.27221e-05` s. Peak electric power is `109.308` W, and the e-folding time is `429.19` s. The checked-in schedule is not rewritten.

## Surface field

From the project directory:

```bash
py -3.14 -m solcascade.production.surface
```

The command writes `runs/production/surface_field.md`. A spherical surface current `(3/2) (B / μ0) sinθ` replaces the distant windings. The field inside the sphere is uniform and the field outside is the matching dipole. The stored energy is three halves of the internal magnetic energy. The quadrature of that exterior dipole is `0.500111`.

The shell-wall current sits on the shell inner radius `0.2` m and holds the vacuum field `0.21639` T. Its magnetic energy is `936.493` J, below the Helmholtz store `15379.1` J. At recovery 0.5 the pulsed supply is `468.247` J. The equator current is `258297` A/m. The checked-in compression still raises the field inside the plasma.

The plasma-surface current at radius `0.018` m holds the compressed field `1.66968` T in the plasma volume. Its magnetic energy is `40.6464` J and the pulsed supply is `20.3232` J. The equator current is `1.99303e+06` A/m. Beta = 1 complete-burn electricity is `854.971` J, above that pulsed supply. Schedule electricity is `0.00241447` J, below it. The gyro-walk on that load returns `1.75507e-06` J in `8.81033e-07` s.

The same surface current on the bore-limited sphere, radius `0.0684` m, stores `2230.35` J. The pulsed supply is `1115.17` J. Complete-burn electricity on that sphere is `46914` J, above the pulsed supply. Schedule electricity is `0.132487` J, below it. Walk electricity is `0.00139063` J. Peak electric power is `109.308` W.

A persistent surface current is already established, so the pulse draws no magnetic energy to recharge it. The interior field is uniform, so the parallel clock is the thermal transit. The schedule integral and the gyro-walk remain the cross-field figures. With the field removed, ions cross the `0.018` m sphere at `v = sqrt(kT / m)`. The transit is `4.67829e-08` s and the electricity in that transit is `9.31944e-08` J. The magnetic supply on that arrangement is zero. The checked-in schedule is not rewritten. The distant windings remain the ledger until a surface current is built.

## Nature

From the project directory:

```bash
py -3.14 -m solcascade.production.nature
```

The command writes `runs/production/nature.md`. The gyro-walk remains the cross-field clock. A particle mirrors along the field only when a throat is stronger than the field it occupies. The isotropic share that misses that throat is `1 - sqrt(1 - 1/R)`.

The beta = 1 sphere is `2.99726` thermal gyroradii across. The bore-limited sphere is `11.3896` thermal gyroradii across. The installed walk is `8.81033e-07` s, electricity `1.75507e-06` J. The thermal transit is `4.67829e-08` s, electricity `9.31944e-08` J.

The four-coil vacuum axis rises by `1.03738` between the midplane and the shell. An isotropic distribution at that vacuum field would lose `0.810182` of its particles. The throat over the compressed internal field is `0.134444`, so the compressed plasma sits above the throat and the internal loss fraction is 1. The Helmholtz axis is strongest at the midplane, vacuum rise `1`, throat over internal `0.1296`. The surface current is uniform inside the sphere, throat over internal `1`, and the dipole axis outside falls as `(a/r)³`. On all three fields an isotropic fuel leaves on the thermal transit. Axial electricity is below the gyro-walk electricity. On the bore-limited sphere the transit is `1.77775e-07` s and its electricity is `1.94323e-05` J, below that sphere's gyro-walk electricity.

Pin `aeb2adad6e80f487` is recorded with the domain scalars. A domain scalar is not a mirror ratio. The engine valve `2.62087e-17` stays out of the mirror ratio. gamma_FSOT does not multiply the walk. The checked-in schedule is not rewritten.

## Throat

From the project directory:

```bash
py -3.14 -m solcascade.production.throat
```

The command writes `runs/production/throat.md`. The gyro-walk remains the cross-field clock. Half of an isotropic distribution is trapped when the missing share is one half, and the mirror ratio is then 4/3.

The frozen internal field is `1.66968` T. The checked-in coils at peak scale reach `0.896294` of that field. The designed pair is a steady bore winding on the checked-in wire, centerline `0.008` m at `±0.191833` m. That is the smallest centerline the wire allows, against the inner wall. One turn at `9335.05` A holds a throat field of `2.22733` T. The internal field on that current is `1.67049` T. The stored energy is `0.288622` J. At recovery 0.5 the pulsed supply is `0.144311` J. A persistent throat current draws no magnetic energy from the supply. Trapped electricity is `4.38768e-07` J and stays on the gyro-walk.

The ratio-1 row on the same pair stores `0.0165795` J and has no trapped set. The one-tenth loss row stores `28.9011` J. The strongest vacuum rise outside the shell is `6.68754`. At large current that mirror ratio approaches `0.866705`, so a steady pair outside the shell stays below the internal field. A pulsed pair outside the shell stores `65037.3` J to reach the same half-trapped ratio. Checked-in fusion remains `0.497703` J. The checked-in schedule is not rewritten.

## Working plants

From the project directory:

```bash
py -3.14 -m solcascade.production.plants
```

The command writes `runs/production/plants.md`. The thermal-reactor row is an operating, self-sustaining chain, TRL 9, with continuous external drive set to false. Delivered energy on its reference unit, one kilogram of U-235 fissioned completely, is `2.70929e+13` J at the file efficiency `0.33`. The file ceiling is `0.45`, a headroom of `1.36364`. The fast-reactor row is the operating fast plant, TRL 8, efficiency `0.4` and ceiling `0.5`, delivered `3.22893e+13` J on one kilogram of Pu-239. Its demonstrated efficiency is `1.21212` times the thermal efficiency.

The thermal notes cite `45` GWd per tonne of heavy metal, which is `3.888e+12` J per kilogram of heavy metal. That discharge stands at `0.047357` of the ranked pure-U-235 density. A heavy-metal kilogram and a pure U-235 kilogram are different loads. The enrichment estimate scores against the ranked kilogram: thermal output over that estimate is `752.582`.

Pin `aeb2adad6e80f487` records Nuclear_Physics and Thermodynamics at `0.936388`. A domain scalar is not a plant efficiency. The nuclear scalar stays out of the fission joules. gamma_FSOT does not multiply the plant. The phase-two burnup on the ranked unit is already `1` with a cap of `1`, so the tenfold window is empty. The refinement inside these plants is the file efficiency ceiling, and the fast row already running on the actinides a thermal reactor leaves behind. The checked-in schedule is not rewritten.

## Fission observables

From the project directory:

```bash
py -3.14 -m solcascade.production.fission_observables
```

The command writes `runs/production/fission_observables.md`. The pinned engine emits eight dimensional nuclear readings. Every one of them is inside the 0.5 percent green gate and, where the row prints a percent tag, inside that tag. The comparison target on that table is the number stored beside the formula.

The published bars are the central values and absolute uncertainties written in the accuracy note. A leaf is installed where that note names one formula. A derived dressing is installed when that dressing sits inside the published bar. The helium-4 leaf residual is `-0.0267531`, the triton leaf residual is `0.129238`, and the proton-moment leaf residual is `-0.236892`. Those three sit inside one uncertainty. The deuteron named-seed windows return three dressings, all inside one uncertainty and all still low: `pi^{-4}` on alpha squared at `-0.838545`, phi on alpha cubed at `-0.382061`, and the square root of e on alpha cubed at `-0.315742`. The suction window holds no named seed. Those are different dressings. Short products of e and phi fill the same windows. The scored deuteron binding is `sqrt(e)+phi^-4` on alpha cubed, residual `-0.000436849`. That dressing also sits inside the mass-excess bar. The mass-excess nearest is a different product, `e^3/phi^5`, and the rounded match is wider than the helium standard `4.23431e-06`. The scored deuteron moment is `sqrt(K)/psi` on alpha squared, residual `0.111115`. It crosses the CODATA central, and its quotient match `5.22163e-06` is wider than the helium standard. The bare moment remains about 21,280 uncertainties low. Both scored readings sit inside one published uncertainty. The wave neutron-proton mass difference is `0.764616` and the wave proton radius is `0.871352`, both inside one uncertainty on the wave formula. The 200 MeV ledger quantum is the plant energy, and no wave result carries it. The channel search stays closed. The checked-in schedule is not rewritten.

## Fission paths

From the project directory:

```bash
py -3.14 -m solcascade.production.fission_paths
```

The command writes `runs/production/fission_paths.md`. The thermal plant and the thermionic row are the two catalog converters on one kilogram of U-235 fissioned completely. The thermal plant delivers `2.70929e+13` J at efficiency `0.33`. The thermionic row delivers `4.10499e+12` J at efficiency `0.05`, which is `0.151515` of the thermal delivery. Its file ceiling `0.3` still delivers `2.46299e+13` J, below the thermal plant. The fast row stays on its own Pu-239 kilogram. The refinement inside the operating plants is the file efficiency ceiling.

The thermal governing source states that about 85 percent of the released energy starts as fission-fragment kinetic energy, and that in solid fuel the fragments stop in a microscopic distance so the energy becomes heat. Neutron, gamma, and delayed releases add about 30 MeV on a total of about 200 MeV. Subtracting that addendum leaves the share `0.85`. The fragment kinetic inventory on the ranked kilogram is `6.97848e+13` J. That inventory is already inside the ranked density, and it is heat in solid fuel. The catalog files no collector efficiency for a fragment that leaves the fuel. The operating electricity on this kilogram is the Rankine product.

Pin `aeb2adad6e80f487` records Nuclear_Physics at `0.936388`, a relative gap of `0.101632` from the share. The comparison ranks products of pi, e, phi, gamma, and G. The nearest is `G^2` at `0.838993`, relative gap `0.0129494`. The next is `e/pi` at `0.865256`. The page figure is a round share. The nearest seed is not installed as the share. A domain scalar is not the fragment share. The nuclear scalar stays out of the fission joules. gamma_FSOT does not multiply the plant. The checked-in schedule is not rewritten.

## Fission balance

From the project directory:

```bash
py -3.14 -m solcascade.production.fission_balance
```

The command writes `runs/production/fission_balance.md`. The same physics page assigns the balance, after the fragment kinetic energy, to gamma rays and to neutron kinetic energy. That balance is the `30` MeV addendum, the share `0.15`, and `1.2315e+13` J on the ranked kilogram. Added to the fragment inventory it makes the ranked density. Its electricity at efficiency `0.33` is `4.06394e+12` J, already inside the Rankine product of `2.70929e+13` J. Scoring it again would count the same joules twice.

The catalog files no shield-recovery efficiency. Neutrons and gamma rays that stop outside the coolant stay a shield load. The page's delayed-neutron proportions are `0.66` percent for U-235, `0.27` percent for U-233, and `0.23` percent for Pu-239. Those figures count neutrons. They are the control margin, and the control margin stays in the core.

The file-ceiling headroom on the whole ledger is `9.85198e+12` J. An upper bound that converts the balance at one, and leaves the rest at the demonstrated efficiency, adds `8.25103e+12` J and is not a filed plant. The headroom is larger. The refinement inside the operating plant is the file efficiency ceiling.

The same foundational products put `gamma/pi` nearest the balance, at `0.183733`, relative gap `0.22489`. The complement of the fragment nearest, `1-G^2` at `0.161007`, sits closer, relative gap `0.0733802`, and is not installed. The nearest seed to the U-235 delayed proportion is the same `gamma/pi`, relative gap `26.8384`. The fine-structure leaf is `0.00729735`, relative gap `0.105659` from that proportion. No seed is installed as the balance or as the delayed proportion. The nuclear scalar stays out of the fission joules. The checked-in schedule is not rewritten.

## Fission ceiling

From the project directory:

```bash
py -3.14 -m solcascade.production.fission_ceiling
```

The command writes `runs/production/fission_ceiling.md`. The thermal ceiling source states a steam cycle near `560` K rejected near `310` K. The Carnot bound of those temperatures is `0.446429`. The catalog ceiling stays `0.45` and sits above that Carnot. The refined bound is the temperature Carnot. On the ranked U-235 kilogram the demonstrated plant delivers `2.70929e+13` J, the temperature Carnot delivers `3.66517e+13` J, and the catalog ceiling delivers `3.69449e+13` J. The step up to the temperature Carnot is `9.55876e+12` J. The catalog ceiling adds `2.93214e+11` J beyond it, and that extra is the round-file margin.

The fast ceiling source states no temperature pair, so the fast catalog ceiling stays `0.5` and that row is not tightened. The nearest seed to the demonstrated efficiency `0.33` is `gamma^2` at `0.333178`, relative gap `0.00963007`. The nearest seed to the temperature Carnot is `phi/pi` at `0.515036`, relative gap `0.153681`. The nearest seed to the fast ceiling is `phi/pi`, relative gap `0.0300724`. The closest of these seeds is still wider than one part in a thousand, so none is installed as an efficiency. The nuclear scalar stays out of the fission joules. The checked-in schedule is not rewritten.

## Fission process

From the project directory:

```bash
py -3.14 -m solcascade.production.fission_process
```

The command writes `runs/production/fission_process.md`. The operating core is light water, with the chain held at k_eff of 1. The filed power density, `1e+08` W per cubic metre, does not set the joules per kilogram. The catalog carries three loads: the ranked pure U-235 density, an enrichment estimate of `200` SWU at `50` kWh per SWU for one highly enriched kilogram, and a heavy-metal discharge of `45` GWd per tonne.

The association's fuel-cycle overview, at that same 45 GWd/t, 33 percent, and 50 kWh per SWU, spends `182000` SWU on `24.3` tonnes of 4.5 percent enriched uranium. That is `7.48971` SWU per kilogram, and the highly enriched estimate is `26.7033` times that separative work. Discharge electricity on the heavy-metal kilogram is `1.28304e+12` J. The temperature Carnot raises it to `1.73571e+12` J, a step of `1.35281`. The fabrication page carries burn-up past `60` GWd/t as enrichment rises toward a normal plant limit of 5 percent U-235, and the overview anchors 60 GWd/t at about 6.5 percent of the uranium atoms. That burn-up step is `1.33333`. The two steps together are `1.80375`, or `2.31429e+12` J per kilogram of heavy metal, and that stack is not a filed plant. The 60 GWd/t discharge is still `0.0631427` of the ranked pure U-235 density. The fast row remains the filed route for the actinides left behind. The nearest seed to each of these process fractions is `gamma/pi`, and the gaps are wider than one, so no seed is installed. The checked-in schedule is not rewritten.

## Fission safety

From the project directory:

```bash
py -3.14 -m solcascade.production.fission_safety
```

The command writes `runs/production/fission_safety.md`. Routine cooling is set by the steam efficiency. At `0.33` the plant dumps `7.30909` MJ of heat per kilowatt-hour, the page's round `7.3` MJ. Plants now being built near `0.36` dump `6.4` MJ. The temperature Carnot dumps `4.464` MJ. The heat per kilowatt-hour falls as the efficiency rises.

Fission stops at shutdown and decay heat continues. The cooling page states no fraction for that heat, so none is scored. The page's path is a turbine bypass into the condenser, then residual heat removal, with emergency core cooling if the primary circuit fails. Kept in circulation, that heat ceases to be a problem after a day or two. Large units with a loss of power are the case the page flags. Decay heat stays a cooling load. The delayed-neutron margin stays in the core.

The small class is `300` MWe or less, against the `1000` MWe nameplate. At the filed power density the large core is `30.303` cubic metres and the `300` MWe boundary is `9.09091` cubic metres. If the shape is kept, the linear scale is `0.669433`. A full-power year at that boundary is `9.4608e+15` J, three tenths of the large nameplate year. The joules per kilogram stay on the discharge, and the fuel stays inside the `5` percent fabrication cap. The fast row's sodium barrier is unchanged. The nearest seed to the size ratio `0.3` is `G/pi` at `0.291561`, relative gap `0.0281303`. The nearest seed to `0.36` is `gamma/phi` at `0.356739`, relative gap `0.00905861`. Both gaps are wider than one part in a thousand, so neither is installed. The checked-in schedule is not rewritten.

## Fission decay

From the project directory:

```bash
py -3.14 -m solcascade.production.fission_decay
```

The command writes `runs/production/fission_decay.md`. The cooling page states no decay-heat fraction. Equation (2) of Untermyer and Weills, AECD-3454 (25 February 1952), is the calorimetry fit for irradiated uranium. It includes the decay heat of U-239 and Np-239, and it leaves delayed-neutron fission out. The run uses the filed full-power year, `3.1536e+07` s. Operating heat is `3.0303e+09` W on the large nameplate and `9.09091e+08` W on the `300` MWe boundary.

At one second the fit gives `0.0582203` of operating heat, `1.76425e+08` W, inside the report's ±50 percent band. The long-operation limit is `0.0588892`. At one day the fit gives `0.00661589`, `2.00482e+07` W, inside the ±10 percent band. The limit is `0.00728386`. At two days the fit gives `0.0052866`, `1.602e+07` W, inside the same ±10 percent band. The limit is `0.00595368`. The small-class loads are `5.29276e+07` W, `6.01445e+06` W, and `4.806e+06` W.

From one second to one day the one-year load falls to `0.113635` of the one-second load. From one day to two days it falls to `0.799076` of the one-day load. The nearest seed to each scored fraction, and to the one-day ratio, is `gamma/pi`. The nearest seed to the two-day ratio is `G^2` at `0.838993`, relative gap `0.0499534`. That is the closest seed. The gap is wider than one part in a thousand, so no seed is installed. The one-day heat ratio and the U-235 delayed-neutron count of 0.66 percent stay apart, with relative gap `0.00240798`. Decay heat stays a cooling load. The checked-in schedule is not rewritten.

## Fission decay seed

From the project directory:

```bash
py -3.14 -m solcascade.production.fission_decay_seed
```

The command writes `runs/production/fission_decay_seed.md`. The watts stay the 1952 fit. The naming pass uses the e-and-phi depth-two ladder from the deuteron check, `1224` expressions with `180` inside (0, 1), then the same operations on pi, e, phi, gamma, and G, `7380` expressions with `1395` inside (0, 1). Median spacings are `0.032156` and `0.00349101`.

The one-second fraction tightens from `gamma/pi` at gap `2.15583` to `e^-6+phi^-6` at gap `0.000231991`, and the five-atom nearest is `phi^-6/sqrt(G)` at gap `0.000137647`. The one-day fraction's five-atom nearest is `pi^-4*G^5` at gap `0.000478739`. The two-day fraction's is `pi^-5*phi` at gap `0.000141649`. The one-day drop's is `gamma^3*G^6` at gap `0.000516464`. Those four meet the one-part-in-a-thousand review bar. The one-day limit stays on `e^-3*phi^-4` at gap `0.00274929`. The two-day drop's five-atom nearest is `pi^-3+G^3` at gap `0.00208217`. Every gap stays wider than the helium quotient match `4.23431e-06`. No expression is installed, and the cooling loads remain `1.76425e+08` W, `2.00482e+07` W, and `1.602e+07` W. The checked-in schedule is not rewritten.

## Fission decay outside

From the project directory:

```bash
py -3.14 -m solcascade.production.fission_decay_outside
```

The command writes `runs/production/fission_decay_outside.md`. The four figures still outside the five-atom review bar are the one-second limit, the one-day limit, the two-day limit, and the two-day drop. The added operation is the complement of a five-atom piece: `2274` expressions, `1395` distinct values inside (0, 1), median spacing `0.000333282`.

The two-day limit complement is `1-(gamma^6+sqrt(G))` at `0.00595383`, relative gap `2.46138e-05`. The two-day drop complement is `1-(gamma^3/sqrt(G))` at `0.799056`, relative gap `2.538e-05`. Both meet the review bar and stay wider than the helium match `4.23431e-06`, so neither is installed. The one-second limit stays with `pi^-5+phi^-6` at gap `0.00181089`. The one-day limit's closest leaf is the fine-structure value `0.00729735`, gap `0.00185228`. The cooling loads remain `1.76425e+08` W, `2.00482e+07` W, and `1.602e+07` W. The checked-in schedule is not rewritten.

## Fission decay dress

From the project directory:

```bash
py -3.14 -m solcascade.production.fission_decay_dress
```

The command writes `runs/production/fission_decay_dress.md`. The two limits the complement left outside are dressed from their five-atom anchors. The factor is the fine-structure leaf `0.00729735` times one base piece, with either sign. Each anchor has `116` such dressings, median spacing `0.000511364`.

The one-second limit's nearest is `(pi^-5+phi^-6)*(1-alpha*phi^-3)` at `0.0588942`, relative gap `8.50946e-05`. Seven dressings of that anchor meet the review bar, and one is inside one part in ten thousand. The one-day limit's nearest is `(e^-3*phi^-4)*(1+alpha*phi^-2)` at `0.00728408`, relative gap `3.03874e-05`. Two dressings of that anchor meet the review bar, and one is inside one part in ten thousand. Both gaps stay wider than the helium match `4.23431e-06`, so neither dressing is installed. The cooling loads remain `1.76425e+08` W, `2.00482e+07` W, and `1.602e+07` W. The checked-in schedule is not rewritten.

## Fission decay recorded

From the project directory:

```bash
py -3.14 -m solcascade.production.fission_decay_recorded
```

The command writes `runs/production/fission_decay_recorded.md`. The six figures the first dressing left untouched already had a recorded nearest inside the review bar. The bare value is that nearest. The factor is the same fine-structure leaf `0.00729735` times one base piece, with either sign. The two limits already dressed stay on that score.

The one-second fraction's nearest is `(phi^-6/sqrt(G))*(1-alpha*e^-4)`, relative gap `3.97292e-06`. The two-day limit's nearest is `(1-(gamma^6+sqrt(G)))*(1-alpha*pi^-5)`, relative gap `7.67233e-07`. The two-day drop's nearest is `(1-(gamma^3/sqrt(G)))*(1+alpha*pi^-5)`, relative gap `1.53457e-06`. Those three meet the helium match `4.23431e-06`. One dressing of each of those anchors meets that match.

The one-day fraction's nearest is `(pi^-4*G^5)*(1-alpha*gamma^5)`, relative gap `1.0935e-05`. The two-day fraction's nearest is `(pi^-5*phi)*(1-alpha*e^-4)`, relative gap `7.97431e-06`. The one-day drop's nearest is `(gamma^3*G^6)*(1+alpha*gamma^5)`, relative gap `4.91253e-05`. Those three stay wider than the helium match. No expression is installed as a decay-heat fraction. The cooling loads remain `1.76425e+08` W, `2.00482e+07` W, and `1.602e+07` W. The checked-in schedule is not rewritten.

## Fission decay open

From the project directory:

```bash
py -3.14 -m solcascade.production.fission_decay_open
```

The command writes `runs/production/fission_decay_open.md`. The five figures still wider than the helium match take the same base-piece factor again. Every added anchor has `116` dressings, median spacing `0.000511364`.

The one-day fraction closes as `((pi^-4*G^5)*(1-alpha*gamma^5))*(1-alpha*pi^-6)`, relative gap `3.34448e-06`. The two-day fraction closes as `((pi^-5*phi)*(1-alpha*e^-4))*(1-alpha*pi^-6)`, relative gap `3.83823e-07`. The one-day drop closes as `((gamma^3*G^6)*(1+alpha*gamma^5))*(1+alpha*e^-5)`, relative gap `4.14453e-08`. The one-second limit needs one factor that stays wider, `((pi^-5+phi^-6)*(1-alpha*phi^-3))*(1-alpha*pi^-4)` at gap `1.01737e-05`, then closes as `(((pi^-5+phi^-6)*(1-alpha*phi^-3))*(1-alpha*pi^-4))*(1-alpha*pi^-6)`, relative gap `2.58319e-06`. The one-day limit needs `((e^-3*phi^-4)*(1+alpha*phi^-2))*(1-alpha*pi^-5)` at gap `6.54067e-06`, then closes as `(((e^-3*phi^-4)*(1+alpha*phi^-2))*(1-alpha*pi^-5))*(1-alpha*pi^-6)`, relative gap `1.0498e-06`.

`pi^-6` closes four of the five with a minus sign. `e^-5` closes the one-day drop with a plus sign. One dressing of each closing anchor meets the helium match `4.23431e-06`. No expression is installed as a decay-heat fraction. The cooling loads remain `1.76425e+08` W, `2.00482e+07` W, and `1.602e+07` W. The checked-in schedule is not rewritten.

## Fission decay load

From the project directory:

```bash
py -3.14 -m solcascade.production.fission_decay_load
```

The command writes `runs/production/fission_decay_load.md`. The cooling load is the helium-tight FSOT fraction times the filed operating heat. That heat is the nameplate divided by the filed efficiency, `3.0303e+09` W on the large nameplate and `9.09091e+08` W on the small-class boundary. The 1952 fit remains the calorimetry witness.

The one-second load is `(phi^-6/sqrt(G))*(1-alpha*e^-4)` times that heat, `1.76426e+08` W, against the fit `1.76425e+08` W. The one-day load is `((pi^-4*G^5)*(1-alpha*gamma^5))*(1-alpha*pi^-6)`, `2.00482e+07` W at six figures. The two-day load is `((pi^-5*phi)*(1-alpha*e^-4))*(1-alpha*pi^-6)`, `1.602e+07` W at six figures. The small-class loads are `5.29278e+07` W, `6.01447e+06` W, and `4.80601e+06` W. The fraction gaps are `3.97292e-06`, `3.34448e-06`, and `3.83823e-07`. The ratio of the derived one-day load to the derived one-second load is `0.113635`, gap `6.28438e-07`. The ratio of the derived two-day load to the derived one-day load is `0.799074`, gap `2.96064e-06`. Every gap meets the helium match `4.23431e-06`. Decay heat stays a cooling load. The checked-in schedule is not rewritten.

## Tests

```bash
pytest
```

The tests cover the on-axis coil field, a Biot-Savart cross-check, inductance, Bosch-Hale's published temperature window, compression scalings, neutron-fraction balance, the adiabatic heat limit, zero thermal stress at uniform temperature, the Lamé inner-wall balance, the electrical ledger, and one full run of the example config. They do not need OpenMC, MOOSE, a GPU, or a network.

## License

MIT. See `LICENSE`.
