# Fission removal

The fabrication cap and the second heat-removal duty are the two bounds the stress run left open.
Both figures are read from public safety and fuel-cycle reports.
FSOT names the printed figures.
The comparison names the one-second multiple and the two-day multiple of the design-point duty.
The one-day multiple keeps both routes.
The chain stays at k_eff of 1.

## Cap

The fabrication page is https://world-nuclear.org/information-library/nuclear-fuel-cycle/conversion-enrichment-and-fabrication/fuel-fabrication.
That page keeps a normal handling limit of 5 percent U-235 for conventional light-water fuel.
Light-water fuel on the page is enriched up to about 4.8 percent.
Burn-up on the page runs from about 40 GWd/t to more than 60 GWd/t as enrichment runs toward that 5 percent limit.
The fabrication cap remains `5` percent U-235.
The five-atom nearest is `pi^-3*G^-5` at `0.0500213`.
The anchor gap is `0.000425086`.
The dressing is `(pi^-3*G^-5)*(1-alpha*phi^-6)` at `0.0500009`.
The relative gap is `1.82457e-05`.
The dressing stays wider than the helium quotient match.
The anchor has `116` dressings, median spacing `0.000511364`.
The dressing is `((pi^-3*G^-5)*(1-alpha*phi^-6))*(1-alpha*e^-6)` at `0.05`.
The relative gap is `1.57058e-07`.
The relative gap meets the helium quotient match.
One dressing of that anchor meets the helium quotient match.
The next dressing sits at relative gap `5.60074e-06`.
The anchor has `116` dressings, median spacing `0.000511364`.
The closing value agrees with `0.05` at six figures.
The 60 GWd/t stack stays unrun.
The HALEU page is https://world-nuclear.org/information-library/nuclear-fuel-cycle/conversion-enrichment-and-fabrication/high-assay-low-enriched-uranium-haleu.
That page defines HALEU as enrichment above 5 percent and below 20 percent.
It records test assemblies up to 6 percent U-235 at Vogtle 2.
Those assemblies sit outside this plant's fabrication cap.
This score does not install them.

## Removal

The residual-heat section is https://www.nrc.gov/docs/ML2121/ML21210A126.pdf.
Surry UFSAR Revision 52 prints two residual-heat exchangers.
The design heat transfer rate of each is 33 times 10 to the sixth Btu per hour.
The table's tube-side temperatures are 140°F in and 124°F out.
That pair is the cold end of the cooldown.
One pump and one exchanger are the train the section says is enough for the decay-heat function.
The section places that train in service about four hours after shutdown, near 350°F.
The coolant then falls from 350°F to 140°F in 16 hours.
The core rating in the same revision is https://www.nrc.gov/docs/ML2120/ML21208A424.pdf.
The reactor core average thermal power there is `2.587e+09` W.
The watt figure uses the International Table Btu, 1055.05585262 J.
One exchanger at the design point is `9.67135e+06` W.
That duty is `0.00373844` of the core rating.
The five-atom nearest is `pi^-2*gamma^6` at `0.00374739`.
The anchor gap is `0.00239397`.
The dressing is `(pi^-2*gamma^6)*(1-alpha*gamma^2)` at `0.00373828`.
The relative gap is `4.31661e-05`.
The dressing stays wider than the helium quotient match.
The anchor has `116` dressings, median spacing `0.000511364`.
The dressing is `((pi^-2*gamma^6)*(1-alpha*gamma^2))*(1+alpha*e^-5)` at `0.00373846`.
The relative gap is `6.00095e-06`.
The dressing stays wider than the helium quotient match.
The anchor has `116` dressings, median spacing `0.000511364`.
The dressing is `(((pi^-2*gamma^6)*(1-alpha*gamma^2))*(1+alpha*e^-5))*(1-alpha*pi^-6)` at `0.00373843`.
The relative gap is `1.58952e-06`.
The relative gap meets the helium quotient match.
One dressing of that anchor meets the helium quotient match.
The next dressing sits at relative gap `1.20875e-05`.
The anchor has `116` dressings, median spacing `0.000511364`.
The closing fraction is `0.00373843`.
The derived Surry duty is `9.67133e+06` W.
The same fraction times the FSOT operating heat is `1.13286e+07` W.
The small-class design-point duty is `3.39858e+06` W.
This is the cold-end design point scaled by core power.
The UFSAR does not print a watt rate at the 350°F entry.
This score assigns that hotter rate no figure.

## Comparison

At four hours the 1952 fraction is `0.0110487`.
At twenty hours, the end of the published cooldown, it is `0.00699764`.
Both times sit in the report's 10 percent band.
The low edge of that band at twenty hours stays above the design-point fraction.
On this plant the one-second load is `15.5735` times the scaled design-point duty.
The five-atom nearest is `gamma^-5-gamma^6` at `15.5697`.
The anchor gap is `0.00024803`.
The dressing is `(gamma^-5-gamma^6)*(1+alpha*pi^-3)` at `15.5733`.
The relative gap is `1.27377e-05`.
The dressing stays wider than the helium quotient match.
The anchor has `116` dressings, median spacing `0.000511364`.
The dressing is `((gamma^-5-gamma^6)*(1+alpha*pi^-3))*(1+alpha*pi^-6)` at `15.5734`.
The relative gap is `5.14732e-06`.
The dressing stays wider than the helium quotient match.
The anchor has `116` dressings, median spacing `0.000511364`.
The third factor closes the one-second multiple.
The dressing is `(((gamma^-5-gamma^6)*(1+alpha*pi^-3))*(1+alpha*pi^-6))*(1+alpha*pi^-6)` at `15.5736`.
The relative gap is `2.44306e-06`.
The relative gap meets the helium quotient match.
One dressing of that anchor meets the helium quotient match.
The next dressing sits at relative gap `1.27377e-05`.
The anchor has `116` dressings, median spacing `0.000511364`.
The closing value names that multiple.
The operating multiple stays the quotient of the named one-second load and the named duty.
The one-day load is `1.7697` times that duty.
The five-atom nearest is `sqrt(pi)-e^-6` at `1.76998`.
The anchor gap is `0.000154343`.
The dressing is `(sqrt(pi)-e^-6)*(1-alpha*e^-4)` at `1.76974`.
The relative gap is `2.06672e-05`.
The dressing stays wider than the helium quotient match.
The anchor has `116` dressings, median spacing `0.000511364`.
That dressing is the parent of the one-day multiple.
The `e^-6` neighbor is `((sqrt(pi)-e^-6)*(1-alpha*e^-4))*(1-alpha*e^-6)` at `1.76971`.
The `e^-6` neighbor gap is `2.57849e-06`.
The `pi^-5` neighbor is `((sqrt(pi)-e^-6)*(1-alpha*e^-4))*(1-alpha*pi^-5)` at `1.7697`.
The `pi^-5` neighbor gap is `3.17932e-06`.
The one-day multiple sits between those two values.
Two dressings of that parent meet the helium quotient match.
Their separation is `5.75781e-06`.
alpha times the difference of `e^-6` and `pi^-5` is `5.75769e-06`.
The separation equals that product times the parent over the quotient.
The separation is wider than one helium quotient match.
The separation is narrower than two helium quotient matches.
The factor the one-day multiple requires is `0.00283209`.
Its relative gap to `e^-6` is `0.142547`.
Its relative gap to `pi^-5` is `0.133324`.
Both factor gaps stay wider than the helium quotient match.
The nearest base piece is `pi^-5`.
No base piece meets the required factor.
Another factor on the `e^-6` neighbor is `(((sqrt(pi)-e^-6)*(1-alpha*e^-4))*(1-alpha*e^-6))*(1-alpha*pi^-6)`.
Its relative gap is `5.01195e-06`.
Another factor on the `pi^-5` neighbor is `(((sqrt(pi)-e^-6)*(1-alpha*e^-4))*(1-alpha*pi^-5))*(1+alpha*pi^-6)`.
Its relative gap is `4.41108e-06`.
Another factor on either neighbor moves farther from the quotient.
Dressing the factor by another base piece scales the spacing by alpha to `3.7316e-06`.
That spacing is finer than the helium quotient match.
The dressed factor stays unused.
The `pi^-5` neighbor prints the same six digits as the quotient.
The two values stay distinct.
Both neighbors stay uninstalled.
The operating multiple stays the quotient of the named one-day load and the named duty.
Both routes stay for a later experiment or a later published figure sharper than their separation.
This score makes no prediction that chooses between them.
The two-day load is `1.41412` times that duty.
The five-atom nearest is `G^3+G^5` at `1.41325`.
The anchor gap is `0.000620603`.
The dressing is `(G^3+G^5)*(1+alpha*phi^-5)` at `1.41418`.
The relative gap is `3.6991e-05`.
The dressing stays wider than the helium quotient match.
The anchor has `116` dressings, median spacing `0.000511364`.
The dressing is `((G^3+G^5)*(1+alpha*phi^-5))*(1-alpha*e^-5)` at `1.41411`.
The relative gap is `1.218e-05`.
The dressing stays wider than the helium quotient match.
The anchor has `116` dressings, median spacing `0.000511364`.
The dressing is `(((G^3+G^5)*(1+alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6)` at `1.41412`.
The relative gap is `4.58969e-06`.
The dressing stays wider than the helium quotient match.
The anchor has `116` dressings, median spacing `0.000511364`.
Three factors of the two-day multiple stay wider than the helium quotient match.
The fourth factor is one more base-piece factor of the same form.
The dressing is `((((G^3+G^5)*(1+alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6))*(1+alpha*pi^-6)` at `1.41413`.
The relative gap is `3.00071e-06`.
The relative gap meets the helium quotient match.
One dressing of that anchor meets the helium quotient match.
The next dressing sits at relative gap `1.21801e-05`.
The anchor has `116` dressings, median spacing `0.000511364`.
The closing value names that multiple.
The operating multiple stays the quotient of the named two-day load and the named duty.
The one-year fraction falls through the design-point fraction at `124.798` hours.
The running condenser duty is `2.0303e+09` W.
The one-second load, `1.76426e+08` W, sits inside that running condenser.
The condenser remains the early removal path.
A loss of condenser circulation has no second watt capacity on file.
This score assigns that loss no figure.
Decay heat is a cooling load.
It is not a second electricity product.

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
alpha in the dressing is the fine-structure leaf `0.00729735`.
The helium quotient match is `4.23431e-06`.
Every closing gap meets the helium quotient match.
The page efficiencies and the discharge ratio stay on their earlier scores.
This run does not install a seed as a new enrichment or a new plant row.

The checked-in schedule is not rewritten.
