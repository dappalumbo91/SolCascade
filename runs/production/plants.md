# Working plants

The thermal-reactor row is an operating chain. It is self-sustaining.
The catalog sets continuous external drive to false.
The fast-reactor row is the operating fast plant on the same fission energy.
The chain sustains once critical.

Delivered joules are the catalog product of energy density and
demonstrated efficiency, on each row's own reference unit.

| Plant | TRL | Efficiency | Ceiling | Delivered (J) | At the ceiling (J) | Output / enrichment estimate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Thermal-reactor fission | 9 | 0.33 | 0.45 | 2.70929e+13 | 3.69449e+13 | 752.582 |
| Fast-reactor fission | 8 | 0.4 | 0.5 | 3.22893e+13 | 4.03617e+13 | 896.926 |

The thermal reference unit is 1 kg U-235, complete fission.
The fast reference unit is 1 kg Pu-239, complete fission.
Thermal headroom, ceiling over demonstrated efficiency, is `1.36364`.
Fast headroom is `1.25`.
The fast demonstrated efficiency over the thermal efficiency is `1.21212`.

The thermal notes cite `45` GWd per tonne of heavy metal.
That discharge is `3.888e+12` J per kilogram of heavy metal.
The ranked thermal density is `8.20998e+13` J per kilogram of pure U-235.
Per kilogram, the discharge citation stands at `0.047357` of that ranked density.
A heavy-metal kilogram and a pure U-235 kilogram are different loads.
The enrichment estimate belongs to one kilogram of highly enriched uranium
and is scored against the ranked kilogram. The fast row uses that same
enrichment estimate.

## FSOT record

Pin `aeb2adad6e80f487` matches the loaded engine prefix `aeb2adad6e80f487`.
The scalars below are dimensionless.
A domain scalar is not a plant efficiency.
The nuclear scalar stays out of the fission joules.
Nuclear_Physics: `0.936388`.
Thermodynamics: `0.936388`.
gamma_FSOT: `2.24054`.
gamma_FSOT does not multiply the plant.
The engine valve stays out of the fission joules.
Engine valve: `2.62087e-17`.

The phase-two burnup on the ranked unit is `1`
with a cap of `1`.
The tenfold window is empty.

The refinement inside these plants is the file efficiency ceiling.
The fast row is the operating plant on the actinides a thermal reactor
leaves behind. The discharge citation is the gap between a heavy-metal
tonne and the ranked pure kilogram. The 200 MeV ledger quantum, the
efficiencies, and that discharge stay on the catalog figures.

The checked-in schedule is not rewritten.
