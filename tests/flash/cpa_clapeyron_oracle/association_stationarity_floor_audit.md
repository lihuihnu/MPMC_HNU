# State-4-liquid association-stationarity floor audit

## Scope

Diagnostic-only audit for the pinned unmodified external implementation:

`ClapeyronThermo/Clapeyron.jl@229b09452f36c2f812486150df0bb43b197bb4e5`

Frozen profile:

`CPA/SRK-physical/simplified-rdf-1.9eta/explicit-site-pairs/v1`

Only `state=4 phase=liquid` is exercised.

No MPMC_HNU production code, CPA parameters, frozen Gate-E thresholds, or pinned Clapeyron source are modified.

Evidence workflow run:

`35297338425`

Evidence code head:

`765448170e7e2c08b8ae73c60c4332bfc1daf9d1`

Association-stationarity artifact digest:

`sha256:c3203e3727c094fad845a74ed66f440c6840f0d64d2ad5fc7042947155650f26`

## Frozen protocol

BigFloat arithmetic precisions:

`256, 512 bits`.

Association solver matrix:

| label | rtol=atol | max_iters |
| --- | ---: | ---: |
| baseline | `1e-16` | 4096 |
| tighter_18 | `1e-18` | 8192 |
| tighter_20 | `1e-20` | 8192 |
| tighter_24 | `1e-24` | 16384 |
| tighter_28 | `1e-28` | 32768 |
| tighter_32 | `1e-32` | 65536 |

All rows retain:

- damping factor `0.5`;
- implicit AD enabled;
- explicit source-complete CR-1 pairs;
- runtime `:nocombining`;
- the exact same frozen state coordinates and parameter snapshot.

The derivative probe is fixed rather than rescanned:

- centered 7-point O(h^6);
- `h/V=1e-6`.

Every shifted scalar value comes only from
`Clapeyron.eos_res(model,V,T,z)`.

## Observed stationarity result

Tightening `rtol/atol` by sixteen orders of magnitude and increasing the iteration budget from 4096 to 65536 has **no effect whatsoever** on this state.

At both 256 and 512 bits, every solver row produces exactly the same:

- mass-action residual;
- four association site fractions `X`;
- reduced residual Helmholtz value `F_res`;
- 7-point `dA_res/dV`;
- reconstructed pressure.

The audit records zero deltas versus the baseline for every tighter row:

`max_abs_x_delta_vs_baseline = 0`

`abs_f_res_delta_vs_baseline = 0`

`abs_dAres_dV_delta_vs_baseline = 0 Pa`.

### Center-state values

At 256 bits:

`X = [`
`0.07217354945206787267277671697942690202820666434685915515112319503651912724011061,`
`0.07217354945206787267277671697942690202820666434685915515112319503651912724011061,`
`0.1830794007502537763200478479156393516863491585841312009147483152660297297547323,`
`0.1830794007502537763200478479156393516863491585841312009147483152660297297547323`
`]`.

Mass-action residual vector:

`[`
`6.2968808828794948067664854843108389520883833499709206097242670643186769527143e-14,`
`6.2968808828794948067664854843108389520883833499709206097242670643186769527143e-14,`
`1.147063547859161267258374540583576555366087759654751974479106634628823493209157e-13,`
`1.147063547859161267258374540583576555366087759654751974479106634628823493209157e-13`
`]`.

Therefore

`max mass-action residual = 1.147063547859161267...e-13`.

The reduced scalar is:

`F_res = -6.097343836105055511041671997355782010185586739670461339580846860328572673380125`.

The fixed high-precision derivative probe gives:

`dA_res/dV = 7.530185116875904996337762534851074566072338815562252447058034867097361549356312e7 Pa`.

and

`P_fd = 72832.00004892297580193620421566404314014153702784098619450056864682576599075169 Pa`.

Thus:

- `abs(P_fd - target) = 4.8922975801936206e-5 Pa`;
- `abs(P_fd - frozen Clapeyron internal P) = 8.396822224905284e-10 Pa`.

The same printed values are obtained for every solver setting at 512 bits.

## Why the solver settings do not matter for this state

The solver-path probe on the exact state reports:

`original_dimension=4`

`compression_selected=true`

`reduced_dimension=2`

`mapping=1,1,2,2`

`initializer_success=true`.

The pinned source performs the following sequence:

1. `X_and_Δ` constructs the four-site association matrix.
2. `__maybe_compress(K)` selects compression.
3. `compress_assoc_matrix(K)` reduces the 4x4 problem to 2x2.
4. `assoc_matrix_solve` first calls `assoc_matrix_x0!`.
5. For a 2x2 matrix, `assoc_matrix_x0!` calls `X_exact2!`.
6. `X_exact2!` returns `success=true`.
7. `assoc_matrix_solve` immediately returns that solution **before** entering the general iterative solver that uses `AssocOptions.rtol`, `AssocOptions.atol`, and `max_iters`.

Therefore the entire pre-frozen solver matrix is intentionally demonstrated to be inactive for state 4 liquid.

The pinned `X_exact2!` path has its own internal stopping scale:

`epsilon = 1e-12 * one(y)`

with at most 100 internal updates.

This hard-coded exact-2 initializer path is independent of the user-facing association solver tolerances.

The observed fixed `~1.147e-13` expanded four-site mass-action residual is consistent with that path. The present audit does not modify the pinned source to prove a causal coefficient-level relationship; it establishes that the user-configurable association stopping tolerances cannot reduce this residual because they are bypassed.

## Answer to the audit question

For state 4 liquid, the remaining derivative offset is **not controlled by the configurable association `rtol/atol/max_iters`**.

Tightening those settings from `1e-16/4096` through `1e-32/65536` changes none of:

- mass-action residual;
- `X`;
- `F_res`;
- `dA_res/dV`;
- pressure.

The reason is structural in the pinned implementation: the state is compressed from 4 association unknowns to 2 and solved by the `X_exact2!` initializer, which returns before the configurable iterative convergence criteria are consulted.

This further narrows the Gate-E discrepancy from a generic "association solver tolerance floor" to the pinned compressed exact-2 association path / its fixed internal numerical behavior for state 4 liquid.

It still does **not** establish a CPA formulation error. The scalar residual Helmholtz potential remains aligned, and the BigFloat scalar derivative remains tightly consistent with Clapeyron's own internal derivative.

## Gate consequence

The frozen Gate-E thresholds remain unchanged:

- pressure: `5e-6 Pa`;
- `ln(phi)`: `1e-10`.

State 4 liquid remains outside the frozen envelope.

Therefore:

- Gate E remains **BLOCKED**;
- Gate F remains **NOT READY**;
- no production source-of-truth switch is authorized;
- no threshold widening is authorized.
