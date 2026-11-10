# Fission process

The operating implementation is a light-water core. Water is the moderator and the coolant.
The chain is held at k_eff of 1.
The core power density is `1e+08` W per cubic metre.
That power density does not set the joules per kilogram.
The fragments stop in the solid fuel, and the Rankine cycle takes the heat.

## Three loads

The ranked density is `8.20998e+13` J per kilogram of pure U-235.
The enrichment basis is `200` SWU at `50` kWh per SWU.
That product is `3.6e+10` J per kilogram of highly enriched uranium.
The notes cite `45` GWd per tonne of heavy metal.
That discharge is `3.888e+12` J per kilogram of heavy metal.
Per kilogram, the discharge stands at `0.047357` of the ranked density.
A heavy-metal kilogram, a highly enriched kilogram, and a pure U-235 kilogram are three loads.

## Normal reload

The material balance is https://world-nuclear.org/information-library/nuclear-fuel-cycle/introduction/nuclear-fuel-cycle-overview.
That page's typical year is a 1000 MWe light-water reactor at 4.5 percent enrichment, 45 GWd/t, and 33 percent.
The overview spends `182000` SWU on `24.3` tonnes of enriched uranium.
That is `7.48971` SWU per kilogram.
At the same `50` kWh per SWU, the enrichment electricity is `1.34815e+09` J per kilogram.
The highly enriched estimate is `26.7033` times that separative work per kilogram.
Discharge electricity over that reload enrichment is `951.705`.
The nameplate year at `1000` MWe for `8760` hours is `3.1536e+16` J.
The fuel-side year, reload mass times discharge electricity, is `3.11779e+16` J.
Those two annual figures are the page's round companions.

## Where the output can still move

Discharge electricity at the demonstrated efficiency is `1.28304e+12` J per kilogram of heavy metal.
The temperature Carnot on that same discharge is `1.73571e+12` J.
The Carnot step on the discharge is `1.35281`.
The fabrication page is https://world-nuclear.org/information-library/nuclear-fuel-cycle/conversion-enrichment-and-fabrication/fuel-fabrication.
That page carries burn-up from about 40 GWd/t past 60 GWd/t as enrichment rises toward 5 percent U-235.
The overview anchors 60 GWd/t at about 6.5 percent of the uranium atoms.
The normal fabrication limit is 5 percent U-235.
The higher burnup band stops at that limit.
The 60 GWd/t discharge is `5.184e+12` J per kilogram.
At the demonstrated efficiency that discharge delivers `1.71072e+12` J.
The burnup step, 60 over 45, is `1.33333`.
The temperature Carnot on the 60 GWd/t discharge is `2.31429e+12` J.
The two steps together are `1.80375`.
That stack is not a filed plant.
The 60 GWd/t discharge stands at `0.0631427` of the ranked density.
The stacked electricity stands at `0.0854203` of the ranked electricity at the demonstrated efficiency.
The fast row remains the filed route for the actinides a thermal discharge leaves behind.

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

The nearest seed to the discharge ratio is `gamma/pi` at `0.183733`.
Its relative gap is `2.87975`.
The nearest seed to the 4.5 percent enrichment is `gamma/pi` at `0.183733`.
Its relative gap is `3.08297`.
The nearest seed to the 6.5 percent atomic burn-up is `gamma/pi` at `0.183733`.
Its relative gap is `1.82667`.
No seed is installed as a process fraction.

The checked-in schedule is not rewritten.
