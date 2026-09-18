# Ten-state independent compressed-association diagnostic

## Scope

This diagnostic extends the previously validated independent compressed-association correction to the complete frozen ten-state Clapeyron phase-kernel set.

Pinned external implementation:

`ClapeyronThermo/Clapeyron.jl@229b09452f36c2f812486150df0bb43b197bb4e5`

Frozen profile:

`CPA/SRK-physical/simplified-rdf-1.9eta/explicit-site-pairs/v1`

No raw pinned oracle value is replaced. No Gate-E threshold, MPMC_HNU production code, CPA parameter, or pinned Clapeyron source is modified.

Evidence workflow run:

`35301193333`

Evidence code head:

`2035c8253f5ab8767e5a2e9c96b2631f06730d9b`

Ten-state diagnostic artifact digest:

`sha256:d78d8c02615585e62f7d72c1cec2172fb05daa9693f24592a6a8c2374c66c438`

## Frozen protocol

The diagnostic uses the same ten frozen phase states as the Gate-E oracle.

BigFloat precisions:

`256, 512 bits`.

For every state:

1. build the pinned original association matrix;
2. probe the pinned compression / initializer path;
3. only compressed 2x2 exact-initializer states may be corrected;
4. solve the reduced mass-action system independently with the previously frozen analytic-Jacobian Newton solver;
5. require `max|R| <= 1e-30` at the center and all six shifted-volume points;
6. expand through the exact pinned compression mapping;
7. rebuild the association residual scalar and total residual scalar;
8. reconstruct pressure with the fixed centered 7-point O(h^6) stencil at `h/V=1e-6`.

Z-only fugacity is evaluated with the already-frozen Clapeyron residual chemical potentials:

`ln(phi_i)_diag = mu_i^res/(RT) - ln(Z_diag)`

with

`Z_diag = P_diag*V/(nRT)`.

The comparison target isolates only the pressure/Z contribution:

`ln(phi_i)_target-Z = mu_i^res/(RT) - ln(Z_target)`.

Frozen diagnostic envelopes remain:

- pressure: `5e-6 Pa`;
- Z-only `ln(phi)`: `1e-10`.

## Solver-path result

All ten frozen states enter the same corrected path:

`compressed_exact2_corrected`.

Therefore:

- corrected states: `10`;
- uncorrected raw-path states: `0`.

The independent correction is not selectively applied only to the previously failing liquid states; it is applicable to the complete ten-state frozen set under the same pinned model topology.

## Complete ten-state result

The 256- and 512-bit results are identical at the reported scale.

| State | abs pressure delta vs target [Pa] | max abs Z-only ln(phi) delta | Pressure gate | ln(phi) gate |
| --- | ---: | ---: | --- | --- |
| state 0 liquid | `5.323330668403745e-7` | `1.3571962033602398e-11` | PASS | PASS |
| state 0 vapor | `5.637441791602343e-12` | `1.437279604212412e-16` | PASS | PASS |
| state 1 liquid | `9.045130259360312e-8` | `1.851537349416336e-12` | PASS | PASS |
| state 1 vapor | `5.504969250096571e-12` | `1.126866709673416e-16` | PASS | PASS |
| state 2 liquid | `2.2093047546199933e-7` | `3.899782451853953e-12` | PASS | PASS |
| state 2 vapor | `1.1082274507083665e-11` | `1.9562018123073614e-16` | PASS | PASS |
| state 3 liquid | `1.958186021763269e-7` | `3.059761276549681e-12` | PASS | PASS |
| state 3 vapor | `9.325687883535563e-12` | `1.4571842688108322e-16` | PASS | PASS |
| state 4 liquid | `1.537448861669867e-8` | `2.1109524133211077e-13` | PASS | PASS |
| state 4 vapor | `4.119539308383041e-12` | `5.656221589937172e-17` | PASS | PASS |

### Global maxima

Across all ten states:

`max abs pressure delta = 5.323330668403745e-7 Pa`

versus the frozen `5e-6 Pa` envelope.

Margin:

`5e-6 / 5.323330668403745e-7 ~= 9.39`.

For the Z-only fugacity contribution:

`max abs ln(phi) delta = 1.3571962033602398e-11`

versus the frozen `1e-10` envelope.

Margin:

`1e-10 / 1.3571962033602398e-11 ~= 7.37`.

Both 256- and 512-bit summaries report:

`result=PASS`.

## Interpretation

The earlier raw pinned Gate-E pressure and ln(phi) failures disappear across the complete ten-state set when the same pinned association matrices are solved to high stationarity independently.

This is stronger than the previous state-4-only result:

- the correction is applicable to every frozen phase state;
- no state requires a different model, parameter, topology, solver class, or finite-difference step;
- no state remains outside the original frozen pressure / Z-only ln(phi) numerical envelope;
- 256- and 512-bit results agree at the reported scale.

The diagnostic therefore supports the following scientific conclusion:

The complete frozen ten-state pressure/Z discrepancy is attributable to the pinned compressed exact-initializer numerical path rather than to a formulation mismatch between MPMC_HNU and the audited Clapeyron CPA profile.

The dominant corrected residual is now state 0 liquid, at only `5.3233e-7 Pa`, already well inside the original `5e-6 Pa` envelope.

## Gate consequence

This diagnostic does **not** replace the raw pinned oracle.

The literal pre-existing Gate-E contract compares against the unmodified pinned Clapeyron output. Those raw output values remain unchanged and still exceed the frozen pressure / ln(phi) thresholds.

Therefore, under the literal frozen contract:

- Gate E remains **BLOCKED**;
- Gate F remains **NOT READY**.

At the same time, the independent diagnostic now demonstrates that the same formulation and same ten frozen state coordinates satisfy the original numerical envelopes once the identified pinned `X_exact2!` stationarity defect is removed without changing the model or thresholds.

No threshold widening and no production source-of-truth switch is authorized by this diagnostic.
