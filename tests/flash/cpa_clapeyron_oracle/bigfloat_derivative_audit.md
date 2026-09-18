# Pinned Clapeyron BigFloat scalar-derivative audit

## Scope

This is a diagnostic-only follow-up to the ordinary-double derivative-floor audit.

Frozen external implementation:

`ClapeyronThermo/Clapeyron.jl@229b09452f36c2f812486150df0bb43b197bb4e5`

Frozen MPMC_HNU profile:

`CPA/SRK-physical/simplified-rdf-1.9eta/explicit-site-pairs/v1`

No MPMC_HNU production formula, parameter, Gate-E threshold, or pinned Clapeyron source is modified.

Evidence workflow run:

`35295045764`

Evidence code head:

`c7a925b954603a40127d22ed1cd49d9fbc2be99c`

BigFloat audit artifact digest:

`sha256:157d4e677069d20a60ab9a29f90d0e317095c43cb99a50f2d3d39285b7508649`

## Capability result

The pinned CPA scalar path **does support BigFloat state evaluation without source modification**.

The capability probe reports:

- `Clapeyron.eos_res(model,V,T,z)` returns `BigFloat` when `V,T,z` are BigFloat;
- association fractions are backed by `BigFloat`;
- repeated scalar evaluation is bit-for-bit stable at the active BigFloat precision;
- the source checkout remains clean.

Observed marker:

`CLAPEYRON_BIGFLOAT_CAPABILITY_OK direct_state_bigfloat=true scalar_type=BigFloat association_eltype=BigFloat`

The generic call `Clapeyron.promote_model(BigFloat, model)` returns without throwing, but the resulting CPA model still reports `eltype == Float64` because the pinned `CPAParam` container is concretely Float64-typed. Therefore **this audit does not rely on full parameter-container promotion**. It uses the exact already-audited Float64 parameter snapshot with BigFloat state variables and BigFloat association/scalar arithmetic. This is the correct isolation for testing whether the previous derivative floor was caused by Float64 state-arithmetic cancellation.

## Frozen protocol

Problem states were frozen before seeing results:

- state 0 liquid
- state 2 liquid
- state 3 liquid
- state 4 liquid

These are the four dense-liquid states that exceeded the ordinary-double finite-difference pressure gate.

BigFloat precisions were frozen as:

`128, 256, 512 bits`.

Relative volume steps `r=h/V` were frozen as:

`2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 5e-5, 2e-5, 1e-5, 5e-6, 2e-6, 1e-6, 5e-7, 2e-7, 1e-7, 5e-8, 2e-8, 1e-8, 2e-9, 1e-10, 1e-11, 1e-12`.

Only two centered stencils were used:

- 5-point, O(h^4)
- 7-point, O(h^6)

Every shifted value comes only from the public scalar API:

`Clapeyron.eos_res(model,V,T,z)`.

Pressure is reconstructed as

`P_fd = nRT/V - dA_res/dV`.

The Z-only fugacity contribution is reconstructed using the frozen Clapeyron residual chemical potentials:

`ln(phi_i)_fdZ = mu_i^res/(RT) - ln(P_fd V/(nRT))`.

The association numerical options are unchanged from the frozen Clapeyron oracle:

- `rtol = 1e-16`
- `atol = 1e-16`
- `max_iters = 4096`
- damping `0.5`
- implicit AD enabled
- explicit source-complete CR-1 pairs
- runtime `:nocombining`.

## Precision-convergence result

Increasing arithmetic precision from 128 to 256 to 512 bits does **not** move the best derivative toward the frozen target pressure.

### Best common step across all four problem states

| Precision | Stencil | Best `r` vs target | max abs dP vs target | max abs dP vs Clapeyron internal | max Z-only dln(phi) vs target |
| ---: | --- | ---: | ---: | ---: | ---: |
| 128 | 5-point | `1e-4` | `4.326496737109466e-5 Pa` | `5.29363703457927e-7 Pa` | `6.059291075280396e-10` |
| 128 | 7-point | `1e-10` | `4.8922975801936185e-5 Pa` | `3.2222783163239513e-7 Pa` | `6.717236350162679e-10` |
| 256 | 5-point | `1e-4` | `4.326496737109466e-5 Pa` | `5.293637034579253e-7 Pa` | `6.059291075280396e-10` |
| 256 | 7-point | `2e-6` | `4.8922975801936206e-5 Pa` | `3.2222783163239513e-7 Pa` | `6.717236350162682e-10` |
| 512 | 5-point | `1e-4` | `4.326496737109466e-5 Pa` | `5.293637034579253e-7 Pa` | `6.059291075280396e-10` |
| 512 | 7-point | `2e-6` | `4.8922975801936206e-5 Pa` | `3.2222783163239513e-7 Pa` | `6.717236350162682e-10` |

The 128/256/512-bit results are numerically identical at the scale relevant to Gate E.

This rules out ordinary Float64 cancellation as the controlling source of the remaining `~1e-5 to 1e-4 Pa` discrepancy.

## Direct convergence to Clapeyron's internal derivative

The strongest evidence is state 4 liquid.

At 256 bits with the 7-point stencil:

| `r=h/V` | abs dP vs frozen target | abs dP vs Clapeyron internal |
| ---: | ---: | ---: |
| `1e-3` | `5.46700889937e-5 Pa` | `5.74627350957e-6 Pa` |
| `5e-4` | `4.90127630149e-5 Pa` | `8.8947530708e-8 Pa` |
| `2e-4` | `4.89233435573e-5 Pa` | `4.71926906276e-10 Pa` |
| `1e-4` | `4.89229815481e-5 Pa` | `8.33936074933e-10 Pa` |
| `2e-5` | `4.89229758023e-5 Pa` | `8.39681854738e-10 Pa` |
| `1e-6` | `4.89229758019e-5 Pa` | `8.39682222491e-10 Pa` |
| `1e-10` | `4.89229758019e-5 Pa` | `8.39682222491e-10 Pa` |
| `1e-12` | `4.89229758019e-5 Pa` | `8.39682222491e-10 Pa` |

The high-precision scalar derivative therefore converges tightly to Clapeyron's own internal derivative pressure, not to the frozen target pressure.

The residual difference to the internal derivative is sub-nanopascal while the target offset remains about `4.8923e-5 Pa`.

## Per-state result at 512 bits

Using the best pre-frozen 7-point step available for each state:

| State | best abs dP vs target | corresponding Z-only dln(phi) | Gate-E pressure/ln(phi) |
| --- | ---: | ---: | --- |
| state 0 liquid | `3.5391273005818973e-6 Pa` | `9.023091810268379e-11` | PASS / PASS |
| state 2 liquid | `4.079693024529503e-6 Pa` | `7.201322149937527e-11` | PASS / PASS |
| state 3 liquid | `1.6603970193107485e-6 Pa` | `2.5944514192462415e-11` | PASS / PASS |
| state 4 liquid | `4.8922975801936206e-5 Pa` | `6.717236350162682e-10` | **FAIL / FAIL** |

The BigFloat audit therefore localizes the unresolved Gate-E discrepancy to **state 4 liquid** under the frozen state set.

The 5-point stencil reaches a slightly smaller target-pressure error for state 4 at a coarser step (`4.160531232226014e-5 Pa` at `r=2e-4`), but this is a truncation-error crossing rather than convergence to the target: as the step is reduced, the 5-point result also converges to the same Clapeyron internal derivative plateau near `4.89229758e-5 Pa` from the target.

## Association residual versus precision

The center-state association-equation residuals are stable rather than improving with higher arithmetic precision.

At 512 bits:

- state 0 liquid: `2.0518670148010931e-13`
- state 2 liquid: `8.835754689685336e-14`
- state 3 liquid: `1.9277171022962471e-13`
- state 4 liquid: `1.1470635478591613e-13`.

The leading digits are unchanged from 128 to 256 to 512 bits.

This is evidence that the remaining floor is not set by Float64 arithmetic precision. It is consistent with the pinned association solve / stationarity / derivative path reaching its own numerical floor under the already-frozen association solver settings.

That interpretation does **not** establish a CPA formulation error. The scalar residual Helmholtz values remain aligned and the BigFloat derivative converges to Clapeyron's own internal derivative. The unresolved issue is specifically the difference between that pinned Clapeyron derivative path and the frozen target in the densest liquid state.

## Answer to the audit question

The remaining `10^-5 to 10^-4 Pa` derivative discrepancy is **not an ordinary Float64 cancellation floor**.

BigFloat state arithmetic is genuinely active, and increasing precision from 128 to 512 bits does not reduce the target-pressure discrepancy. High-precision 5/7-point derivatives instead converge to Clapeyron's internal derivative with sub-micropascal, and in state 4 sub-nanopascal, agreement.

Therefore the remaining Gate-E blocker belongs to the pinned Clapeyron association/scalar-derivative numerical path under the frozen solver settings, not to Float64 finite-difference cancellation.

## Gate consequence

The frozen Gate-E limits remain unchanged:

- pressure: `5e-6 Pa`
- `ln(phi)`: `1e-10`.

State 4 liquid still exceeds both.

Therefore:

- Gate E remains **BLOCKED**;
- Gate F remains **NOT READY**;
- no production Helmholtz source-of-truth switch is authorized;
- no threshold widening is authorized.
