# Pinned Clapeyron derivative-floor audit

## Scope

Diagnostic audit only for:

- MPMC_HNU profile:
  `CPA/SRK-physical/simplified-rdf-1.9eta/explicit-site-pairs/v1`
- pinned unmodified external implementation:
  `ClapeyronThermo/Clapeyron.jl@229b09452f36c2f812486150df0bb43b197bb4e5`
- the same ten frozen `(T,V,n)` phase states used by the Gate-E comparison.

No production MPMC formula, parameter, or frozen Gate-E threshold is changed.

Evidence workflow run:

`35292787546`

Evidence code head:

`4edf2f23c372b9953dd728296865a617e1e59d63`

The uploaded derivative-floor artifact has digest:

`sha256:1b75cad62c3a3256dab19e6e806513d4d2e73f4ace32020564af4cc64120c2c9`

## Frozen numerical protocol

Only Clapeyron's public scalar residual Helmholtz API is sampled:

`Clapeyron.eos_res(model,V,T,z)`.

The relative volume step is `r=h/V`, frozen before seeing results as:

`1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 5e-5, 2e-5, 1e-5, 5e-6, 2e-6, 1e-6, 5e-7, 2e-7, 1e-7, 5e-8, 2e-8, 1e-8`.

Three centered first-derivative stencils are evaluated:

- 3-point, O(h^2):
  `(A(V+h)-A(V-h))/(2h)`
- 5-point, O(h^4):
  `(A(V-2h)-8A(V-h)+8A(V+h)-A(V+2h))/(12h)`
- 7-point, O(h^6):
  `(-A(V-3h)+9A(V-2h)-45A(V-h)+45A(V+h)-9A(V+2h)+A(V+3h))/(60h)`

Pressure is reconstructed as

`P_fd = nRT/V - dA_res/dV`.

With the already-frozen Clapeyron residual chemical potential held fixed, the same reconstructed pressure gives

`Z_fd = P_fd V/(nRT)`

and

`ln(phi_i)_fdZ = mu_i^res/(RT) - ln(Z_fd)`.

This isolates the V-derivative / Z contribution from the composition derivative.

## Reproducibility checks

The center-point scalar residual Helmholtz values reproduced the frozen oracle exactly:

- maximum `F_res` delta versus the frozen oracle: `0`
- maximum repeat-evaluation delta in `A_res`: `0 J`

Therefore the observed derivative floor is not caused by nondeterministic center-point scalar evaluations.

## Global step-convergence result

The following values are the best common frozen step for all ten states for each stencil.

| Stencil | best `r` | max abs pressure delta vs target | max abs pressure delta vs Clapeyron internal derivative | max abs Z-only `ln(phi)` delta vs target-Z |
| --- | ---: | ---: | ---: | ---: |
| 3-point O(h^2) | `5e-6` | `6.826743483543396e-2 Pa` | `6.826703915430699e-2 Pa` | `1.397432799699061e-6` |
| 5-point O(h^4) | `2e-4` | `8.723437786102295e-4 Pa` | `8.72739459737204e-4 Pa` | `1.7856868872456744e-8` |
| 7-point O(h^6) | `1e-3` | `1.6723392764106393e-4 Pa` | `1.674343948252499e-4 Pa` | `3.423277099290356e-9` |

The frozen Gate-E limits remain:

- pressure: `5e-6 Pa`
- `ln(phi)`: `1e-10`.

None of the finite-difference stencils reaches those global limits.

## U-shaped error evidence

The higher-order stencils show the expected truncation-error / cancellation-noise tradeoff.

### 5-point pressure, max over ten states

| `r` | max abs pressure delta vs target |
| ---: | ---: |
| `1e-3` | `1.9825716316699982e-1 Pa` |
| `5e-4` | `1.2172073125839233e-2 Pa` |
| `2e-4` | **`8.723437786102295e-4 Pa`** |
| `1e-4` | `1.594025205122307e-3 Pa` |
| `5e-5` | `2.5708824396133423e-3 Pa` |
| `1e-5` | `1.0932639241218567e-2 Pa` |

### 7-point pressure, max over ten states

| `r` | max abs pressure delta vs target |
| ---: | ---: |
| `5e-3` | `4.22602117061615e-1 Pa` |
| `2e-3` | `1.6250759363174438e-3 Pa` |
| `1e-3` | **`1.6723392764106393e-4 Pa`** |
| `5e-4` | `4.8138201236724854e-4 Pa` |
| `2e-4` | `8.417604622081853e-4 Pa` |
| `1e-4` | `1.7198026180267334e-3 Pa` |

Reducing the step beyond the optimum worsens the result, so the Gate-E pressure gap cannot be removed by simply taking a smaller ordinary-double finite-difference step.

## State-specific 7-point minima

Even allowing each state to choose its own best step from the frozen sequence, several dense liquid states remain above the `5e-6 Pa` Gate-E pressure threshold.

| State | best abs pressure delta vs target | best `r` |
| --- | ---: | ---: |
| state 0 liquid | `1.1965632438659668e-4 Pa` | `1e-3` |
| state 0 vapor | `2.6082852855324745e-7 Pa` | `5e-4` |
| state 1 liquid | `7.301568984985352e-7 Pa` | `1e-3` |
| state 1 vapor | `3.857145202346146e-6 Pa` | `2e-5` |
| state 2 liquid | `1.0542571544647217e-4 Pa` | `1e-3` |
| state 2 vapor | `4.3637555791065097e-7 Pa` | `5e-4` |
| state 3 liquid | `2.9593706130981445e-5 Pa` | `1e-3` |
| state 3 vapor | `2.765227691270411e-7 Pa` | `1e-4` |
| state 4 liquid | `2.0757317543029785e-5 Pa` | `5e-4` |
| state 4 vapor | `1.2751843314617872e-7 Pa` | `1e-3` |

Four of the five dense liquid states still fail the frozen pressure threshold after per-state step optimization.

## Interpretation

The pinned Clapeyron internal derivative pressure itself differs from the frozen target pressure by at most

`4.8923815484158695e-5 Pa`

(state 4 liquid), and the corresponding internal Z-only `ln(phi)` difference is at most

`6.717355560681426e-10`.

The best global 7-point scalar finite difference is less accurate:

- pressure: `1.6723392764106393e-4 Pa`
- Z-only `ln(phi)`: `3.423277099290356e-9`.

Therefore this audit does **not** show that Clapeyron's internal AD pressure is an avoidable numerical mistake that can be corrected by ordinary-double finite differencing of the public scalar API. It shows instead that direct finite differencing of that scalar API has a coarser derivative floor than Clapeyron's current internal derivative path for these states.

The scalar potential itself remains strongly aligned: the frozen Gate-E scalar `F_res` row passes, and this audit reproduces the center-point scalar values exactly. The unresolved discrepancy is derivative-level and is dominated by the dense liquid states.

## Gate consequence

Gate E remains **BLOCKED** under the already-frozen pressure and `ln(phi)` envelopes.

This audit does not authorize:

- widening the `5e-6 Pa` pressure threshold;
- widening the `1e-10` `ln(phi)` threshold;
- replacing Clapeyron's internal derivative with finite differences;
- switching MPMC_HNU production pressure / residual chemical potential / `ln(phi)` to the Helmholtz single source of truth.

