# Architecture

The pipeline in `solcascade.pipeline` is the only place the pieces are ordered. Each package accepts plain data and returns plain data.

```
YAML + scalar provider
        |
        v
   slot map  ---- coils, shell, materials, schedule
        |
        v
   plasma history (fusion power, fields, triple product)
        |
        v
   neutron deposition (first collision, or OpenMC heating)
        |
        +---> transient heat ----> stress
        |
        v
   conversion ledger, MOOSE input, Cardinal / TokaMaker / NEAMS files
```

## Slots

`collect_slots` flattens every physical number into a string key before a solver runs. `apply_provider` asks the provider for each key. A mapping provider that names an unknown key raises, so a typo does not vanish.

`plasma.thermal_energy_density_J_m3` is special. If the provider sets it, the ion temperature is rewritten so that a monatomic Te = Ti plasma has that energy density at the resolved ion density. If the provider does not set it, the slot is recomputed from density and temperature.

Non-numeric choices stay in YAML: backend names, the flux-conservation law, element symbols, and discretization counts.

A later FSOT module should implement `ScalarProvider` and return already-computed floats. It should not be imported by the coil, transport, or heat modules. Those modules never see where a number came from. Provenance is stored on the summary as `config` or the provider name.

## Backends

| Concern | Always available | Optional local code |
| --- | --- | --- |
| Magnetic field | Elliptic-integral loops | CuPy Biot-Savart grid |
| Neutron heating | First-collision shell | OpenMC fixed source |
| Heat | Backward-Euler radial solve | MOOSE heat-transfer input |
| OpenMC inside MOOSE | Cardinal skeleton | A built Cardinal app |
| Equilibrium | Vacuum field and 0-D pressure | TokaMaker, if imported, and only with a mesh |
| Neutronics coupling | JSON power table | Griffin, not executed |
| Structural response | Lamé and thermal stress | JSON table for Bison, not executed |

Optional codes are detected at run time. A missing executable writes the input and continues with the native result. `transport.backend: openmc` is the exception: if OpenMC cannot run, the pipeline raises, because that setting asked for OpenMC specifically.

## Adding a scalar without refactoring

1. Run `solcascade keys` and choose an existing key, or add one field to `collect_slots` and to the object builder that already reads slots.
2. Return the value from a provider.
3. Leave the solver equations alone.

Material target properties are already keys of the form `materials.<name>.<field>`. Field geometry keys are `coils.<name>.radius_m`, `z_m`, `turns`, `current_a`, and `wire_radius_m`.
