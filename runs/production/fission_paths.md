# Fission paths

The thermal plant and the thermionic row are the two catalog converters
on one kilogram of U-235 fissioned completely.
Delivered joules are the catalog product of that shared density and each
row's demonstrated efficiency.

| Converter | TRL | Efficiency | Ceiling | Delivered (J) | At the ceiling (J) | Output / enrichment estimate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Thermal-reactor fission | 9 | 0.33 | 0.45 | 2.70929e+13 | 3.69449e+13 | 752.582 |
| Thermionic converter | 7 | 0.05 | 0.3 | 4.10499e+12 | 2.46299e+13 | 114.028 |

The thermal reference unit is 1 kg U-235, complete fission.
The thermionic reference unit is 1 kg U-235, complete fission, thermionic conversion.
The shared density is `8.20998e+13` J per kilogram.
The thermionic row is the other filed converter on this same U-235 kilogram.
The thermionic notes name TOPAZ as a flown thermionic converter.
The thermionic row delivers `4.10499e+12` J.
The thermionic row delivers less electricity than the thermal plant.
The delivered ratio, thermionic over thermal, is `0.151515`.
Even at its file ceiling the thermionic row delivers less than the thermal plant.
Both rows use the same enrichment estimate.
The fast reference unit is 1 kg Pu-239, complete fission.
The fast row stays on that unit.
The refinement inside the operating plants is the file efficiency ceiling.

## Fragment inventory

The thermal governing source is https://world-nuclear.org/information-library/nuclear-fuel-cycle/introduction/physics-of-nuclear-energy.
That page states that about 85 percent of the released energy starts as fission-fragment kinetic energy.
In solid fuel the fragments stop in a microscopic distance and the energy becomes heat.
Neutron, gamma, and delayed releases add about `30` MeV on a total of about `200` MeV.
Subtracting that addendum from the ledger leaves the share `0.85`.
The catalog density stays on the 200 MeV ledger.
The fragment kinetic inventory on the ranked kilogram is `6.97848e+13` J.
The fragment share is already inside that density.
The Rankine efficiency multiplies the whole ranked density.
The fragment inventory is heat in solid fuel.
The fragment inventory is not electricity.
The catalog files no collector efficiency for a fragment that leaves the fuel.
The operating electricity on this kilogram is the Rankine product.

## FSOT record

Pin `aeb2adad6e80f487` matches the loaded engine prefix `aeb2adad6e80f487`.
The scalars below are dimensionless.
A domain scalar is not a plant efficiency.
A domain scalar is not the fragment share.
The nuclear scalar stays out of the fission joules.
Nuclear_Physics: `0.936388`.
Thermodynamics: `0.936388`.
The domain scalar's relative gap to the share is `0.101632`.
The domain scalar sits farther from the share than the nearest seed.
gamma_FSOT: `2.24054`.
gamma_FSOT does not multiply the plant.
The engine valve stays out of the fission joules.
Engine valve: `2.62087e-17`.

The comparison uses each of pi, e, phi, gamma, and G, its square, its reciprocal,
and each pair's product and two quotients.
`45` products are named and `19` lie strictly between zero and one.
The nearest seed is `G^2` at `0.838993`.
Its relative gap to the share is `0.0129494`.
The next seed is `e/pi` at `0.865256`.
Its relative gap to the share is `0.0179482`.
The page figure is a round share.
The nearest seed is not installed as the share.

The checked-in schedule is not rewritten.
