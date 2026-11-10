# Fission balance

The thermal governing source is https://world-nuclear.org/information-library/nuclear-fuel-cycle/introduction/physics-of-nuclear-energy.
That page assigns the balance, after the fragment kinetic energy, to gamma rays and to the kinetic energy of the neutrons.
Neutron, gamma, and delayed releases add about `30` MeV on a total of about `200` MeV.
That balance is the share `0.15`.
The catalog density stays on the 200 MeV ledger.
The balance inventory on the ranked kilogram is `1.2315e+13` J.
The fragment inventory on the ranked kilogram is `6.97848e+13` J.
The fragment inventory and the balance inventory add to the ranked density.
The balance is already inside the ranked density.
The Rankine product already includes it.
Its electricity at the demonstrated efficiency is `4.06394e+12` J, already inside the Rankine product of `2.70929e+13` J.
Scoring that electricity again would count the same joules twice.

## Shield and control

The catalog files no shield-recovery efficiency.
Neutrons and gamma rays that stop outside the coolant are a shield load.
The shield stays a shield.
The page gives these delayed-neutron proportions.

| Nucleus | Delayed neutrons (percent) |
| --- | ---: |
| U-235 | 0.66 |
| U-233 | 0.27 |
| Pu-239 | 0.23 |

Those proportions count neutrons.
They are not an energy share of the ledger.
Delayed neutrons enable control of the chain.
The control margin stays in the core.

## Ceiling

The file-ceiling headroom on the whole ledger is `9.85198e+12` J.
An upper bound that converts the balance at one, and leaves the rest of the ledger at the demonstrated efficiency, adds `8.25103e+12` J.
That bound is not a filed plant.
The file-ceiling headroom is larger than that bound.
The refinement inside the operating plant is the file efficiency ceiling.

## FSOT record

Pin `aeb2adad6e80f487` matches the loaded engine prefix `aeb2adad6e80f487`.
The scalars below are dimensionless.
A domain scalar is not a plant efficiency.
A domain scalar is not the balance share.
The nuclear scalar stays out of the fission joules.
Nuclear_Physics: `0.936388`.
Thermodynamics: `0.936388`.
gamma_FSOT: `2.24054`.
gamma_FSOT does not multiply the plant.
The engine valve stays out of the fission joules.
Engine valve: `2.62087e-17`.

The balance uses the same foundational products as the fragment share.
The balance grammar's nearest seed is `gamma/pi` at `0.183733`.
Its relative gap to the balance is `0.22489`.
The complement of the fragment nearest is `1-G^2` at `0.161007`.
Its relative gap to the balance is `0.0733802`.
The complement sits closer than the balance grammar.
The complement is not installed as the balance.
The nearest seed to the U-235 delayed proportion is `gamma/pi` at `0.183733`.
Its relative gap to that proportion is `26.8384`.
The fine-structure leaf is `0.00729735`.
Its relative gap to the U-235 delayed proportion is `0.105659`.
The fine-structure leaf is not the delayed-neutron proportion.
The page figures stay round figures.
No seed is installed as the balance or as the delayed proportion.

The checked-in schedule is not rewritten.
