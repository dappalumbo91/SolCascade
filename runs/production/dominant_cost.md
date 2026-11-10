# Dominant cost

The baseline coil store is the largest energy in this ranking.
It is the magnetic energy of the four checked-in coils at full current.
The fuel ceiling of this sphere sits under that store, so every burn
fraction of this load stays under the baseline draw.

The Helmholtz winding is the smaller store at the same on-axis field.
Its practical supply is the ohmic ledger plus the unrecovered part of
the coil energy. The unrecovered Helmholtz store is larger than the
ohmic ledger and larger than the beta = 1 complete-burn electricity.

## Ranking

| Energy | Joules |
| --- | ---: |
| Baseline coil store at full current | 2.79337e+06 |
| Baseline supply draw | 1.42387e+06 |
| Original fuel ceiling | 73778.1 |
| Helmholtz coil store | 15379.1 |
| Helmholtz practical supply | 8517.06 |
| Unrecovered Helmholtz store | 7689.56 |
| Beta = 1 complete-burn electricity | 854.971 |
| Helmholtz ohmic ledger | 827.501 |
| Magnetic energy inside the checked-in plasma | 27.0976 |
| Beta = 1 schedule electricity | 0.00241447 |
| Beta = 1 walk electricity | 1.75507e-06 |

The Helmholtz store is `567.545` times the magnetic
energy inside the checked-in plasma. The coils carry the field that
the reacting volume does not.

Beta = 1 complete burn versus the ohmic ledger: above.
Beta = 1 complete burn versus the unrecovered store: below.

## Sphere inside the bore

The shell inner radius is `0.2` m. The store-matching
sphere is the beta = 1 density at the checked-in temperature, with the
radius schedule scaled so a complete burn returns the practical supply.
Its initial radius is `0.107584` m and its compressed
radius is `0.0387304` m, a factor `2.15169`
on the checked-in radii. The midplane field of the Helmholtz pair at
the initial edge is `0.996268` of the on-axis field.
Evolved peak beta is `1`. Schedule electricity is
`0.0240525` J. Walk electricity is `8.09453e-05` J.
Complete-burn electricity is `8517.06` J.
The store-matching sphere fits inside the shell.

The bore-limited sphere sets the initial radius to `0.95` of
the shell inner radius and keeps the same compression ratio and the
same beta = 1 density. Initial radius `0.19` m,
compressed radius `0.0684` m, factor `3.8`.
The midplane field at that initial edge is `0.957682`
of the on-axis field. Evolved peak beta is `1`.
Schedule fusion is `4.78603` J and schedule electricity is
`0.132487` J. Walk electricity is `0.00139063` J
over `1.27221e-05` s. Complete fusion is `1.69475e+06` J
and complete-burn electricity is `46914` J.
Peak fusion power is `3948.72` W. Peak electric power is
`109.308` W. The e-folding time at that peak is
`429.19` s.

The bore-limited complete burn sits above the practical supply.
The bore-limited schedule stays under the practical supply.
Bore-limited complete fusion versus one catalog JET pulse: below.

## Approach

The payment to bring under the fuel is the external coil store.
On the Helmholtz winding the unrecovered half of that store is the
term above the beta = 1 electricity. A complete burn of a larger
sphere inside this bore meets that term. The evolved schedule and the
gyro-walk on that sphere remain the small integrals. Peak electric
power is the rate at the compressed state for as long as the plasma
is held. This ranking does not rewrite the checked-in schedule.
