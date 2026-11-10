# Fission decay

The decay-heat fit is equation (2) of Untermyer and Weills, AECD-3454, 25 February 1952, https://www.osti.gov/servlets/purl/4376281.
The cooling page states no decay-heat fraction.
The fraction below is that 1952 fit, evaluated on the filed plant.
It is a cooling load.
It is not a second electricity product.
The fit includes the decay heat of U-239 and Np-239.
The 1952 report leaves delayed-neutron fission out of the fit.
The delayed-neutron margin stays in the core.
The report states that equation (2) sits about 22 percent above its fission-product sum.
The fit is scored as published, with U-239 and Np-239 left in.
The filed operation time is `3.1536e+07` s, one full-power year.
Operating heat on the large nameplate is `3.0303e+09` W.
Operating heat on the small-class boundary is `9.09091e+08` W.

## Cooling times

At one second after the filed full-power year, the fit gives `0.0582203` of operating heat.
The report's stated band at that time is ±50 percent.
That cooling load is `1.76425e+08` W on the large nameplate.
On the small-class boundary the same fraction is `5.29276e+07` W.
The formula's limit as operation time grows is `0.0588892` of operating heat, `1.78452e+08` W.
The nearest seed to the one-year fraction is `gamma/pi` at `0.183733`.
Its relative gap is `2.15583`.
The nearest seed to the formula's limit is `gamma/pi` at `0.183733`.
Its relative gap is `2.11998`.

At one day after the filed full-power year, the fit gives `0.00661589` of operating heat.
The report's stated band at that time is ±10 percent.
That cooling load is `2.00482e+07` W on the large nameplate.
On the small-class boundary the same fraction is `6.01445e+06` W.
The formula's limit as operation time grows is `0.00728386` of operating heat, `2.20723e+07` W.
The nearest seed to the one-year fraction is `gamma/pi` at `0.183733`.
Its relative gap is `26.7715`.
The nearest seed to the formula's limit is `gamma/pi` at `0.183733`.
Its relative gap is `24.2247`.

At two days after the filed full-power year, the fit gives `0.0052866` of operating heat.
The report's stated band at that time is ±10 percent.
That cooling load is `1.602e+07` W on the large nameplate.
On the small-class boundary the same fraction is `4.806e+06` W.
The formula's limit as operation time grows is `0.00595368` of operating heat, `1.80415e+07` W.
The nearest seed to the one-year fraction is `gamma/pi` at `0.183733`.
Its relative gap is `33.7545`.
The nearest seed to the formula's limit is `gamma/pi` at `0.183733`.
Its relative gap is `29.8605`.

## How the load falls

From one second to one day the one-year load falls to `0.113635` of the one-second load.
The nearest seed to that ratio is `gamma/pi` at `0.183733`.
Its relative gap is `0.616868`.
From one day to two days the one-year load falls to `0.799076` of the one-day load.
The nearest seed to that ratio is `G^2` at `0.838993`.
Its relative gap is `0.0499534`.
The closest of these seeds is still wider than one part in a thousand.
No seed is installed as a decay-heat fraction.
The one-day figure counts heat.
The U-235 figure of 0.66 percent counts delayed neutrons.
Their relative gap is `0.00240798`.
That gap is wider than one part in a thousand.
The one-day heat ratio stays off that neutron count.

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

The checked-in schedule is not rewritten.
