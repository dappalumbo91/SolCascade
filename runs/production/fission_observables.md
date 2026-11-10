# Fission observables

The pinned engine emits the reading. The unit is the unit in the name.
The wave target is the number stored beside the formula.
The printed percent tag is the number stored on that row.
The green gate is a relative error of 0.5 percent.
The engine summary gate is a relative error under 5 percent.

| Reading | Formula | Computed | Wave target | Error % | Tag % | Green |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Deuteron_binding_MeV | √e/e + φ | 2.22456 | 2.22457 | 6.07551e-05 | 0.7 | inside |
| Triton_binding_MeV | e² + G⁻¹ | 8.4808 | 8.482 | 0.0141457 | 0.1 | inside |
| He4_binding_MeV | πγ/γ⁵ | 28.3007 | 28.3 | 0.00253334 | 0.5 | inside |
| m_n-m_p_MeV | γ⁴/|Chaos| + P_var | 1.29333 | 1.29333 | 0.000216538 | 0.3 | inside |
| mu_p_muN | π·(1 − γ⁴) | 2.79285 | 2.79285 | 7.93549e-05 | 0.2 | inside |
| r_p_fm | G⁷ + P_new | 0.841249 | 0.8414 | 0.0180008 | 0.1 | inside |
| Neutron_lifetime_s | π⁷·θ_S | 878.593 | 878.4 | 0.0219548 | 0.4 | inside |
| Deuteron_mu_muN | G⁴ + Poof | 0.857391 | 0.8574 | 0.00100092 | — | inside |

Every wave reading is inside the green gate.
Every tagged reading is inside its printed percent tag.

The published bar is the central value and the absolute uncertainty
written in the accuracy note. The residual counts those uncertainties.
Positive means the reading is above the central value.
A leaf is installed where the accuracy note names one formula.
A derived dressing is installed when that dressing sits inside the published bar.

| Reading | Central | Uncertainty | Bare residual | Leaf residual | Scored |
| --- | ---: | ---: | ---: | ---: | --- |
| Deuteron_binding_MeV | 2.2245662 | 4e-07 | -3.87884 | -0.000436848 | inside |
| Triton_binding_MeV | 8.4817962 | 9e-07 | -1106.71 | 0.129238 | inside |
| He4_binding_MeV | 28.2956624 | 8e-07 | 6318.17 | -0.0267531 | inside |
| m_n-m_p_MeV | 1.29333251 | 3.8e-07 | 0.764616 | — | inside |
| mu_p_muN | 2.79284734463 | 8.2e-10 | 5941.02 | -0.236892 | inside |
| r_p_fm | 0.8409 | 0.0004 | 0.871352 | — | inside |
| Deuteron_mu_muN | 0.8574382335 | 2.2e-09 | -21279.7 | 0.111115 | inside |

The helium-4 leaf multiplies the wave formula by one minus alpha squared times pi plus P_base.
The triton leaf multiplies the wave formula by one plus (Poof times Suction) squared times gamma times psi_con squared.
The proton moment leaf is the g-factor leaf divided by two.
The alpha-squared window holds pi to the minus four.
The alpha-cubed window holds phi and the square root of e.
The suction window holds no named seed.
Those are different dressings.
`pi^{-4}` on alpha squared reads `-0.838545`.
`phi` on alpha cubed reads `-0.382061`.
`sqrt(e)` on alpha cubed reads `-0.315742`.
The helium quotient matches pi plus P_base at `4.23431e-06`.
The rounded same-side nearest is `sqrt(e)+phi^-4` on alpha cubed at `-0.000436849`.
The mass-excess nearest is `e^3/phi^5` on alpha cubed at `-0.0338697`.
The two centrals return different products.
The rounded match is wider than the helium standard.
The scored deuteron binding is the rounded same-side nearest.
That dressing sits inside the mass-excess bar.
The helium-4 leaf sits inside one published uncertainty.
The triton leaf sits inside one published uncertainty.
The proton-moment leaf sits inside one published uncertainty.
The deuteron binding sits inside one published uncertainty.
`sqrt(K)/psi` on alpha squared reads `0.111115`.
The moment quotient match is `5.22163e-06`.
That dressing crosses the central.
Its quotient match is wider than the helium standard.
The scored deuteron moment is `sqrt(K)/psi` on alpha squared.
The deuteron moment sits inside one published uncertainty.

## Fission quantum

A domain scalar is not a megaelectronvolt.
The nuclear scalar stays out of the fission joules.
Nuclear_Physics: `0.936388`.
Thermodynamics: `0.936388`.
gamma_FSOT: `2.24054`.
gamma_FSOT does not multiply the plant.
The engine valve stays out of the fission joules.
Engine valve: `2.62087e-17`.
Pin `aeb2adad6e80f487` matches the loaded engine prefix `aeb2adad6e80f487`.

The 200 MeV ledger quantum is the plant energy. The file value is `200` MeV.
No wave result carries that quantum.
The chain, the efficiencies, and the discharge stay on the catalog figures.

The channel search stays closed.

The checked-in schedule is not rewritten.
