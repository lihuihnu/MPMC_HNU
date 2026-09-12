# CPA PT stability adapter

## Scope

This gate connects the validated CPA fixed-composition PT phase-property kernel to the existing model-neutral tangent-plane-distance search. It does not add a second stability algorithm and does not implement phase split.

The CPA-specific convention is:

`CPA/PT/stability/minimum-Gibbs-mechanically-admissible-root/v1`.

`CpaStabilityEvaluator` owns one configured `CpaPtPhase`. For every composition requested by the generic TPD driver it performs a fresh bounded CPA density-root search and returns one `StabilityPhase`.

## Same-composition root policy

A root is eligible only when the CPA PT kernel has returned a resolved simple root and `dP/drho > 0`, which is equivalent to `dP/dv < 0` at fixed `T,n`. This is mechanical admissibility only; it is not compositional stability.

If multiple mechanically admissible simple roots exist at the same `p,T,x`, the adapter compares

```text
g_k/(RT) difference = sum_i x_i ln(phi_i,k)
```

because ideal/reference terms cancel at common `p,T,x`. The minimum-Gibbs root supplies the fugacity coefficients to TPD.

If two admissible roots are equal within the existing floating-point Gibbs comparison guard, the selected value is retained but `StabilityPhase::smooth=false`. The generic TPD search may use a robust negative TPD point as instability evidence, but it may not accept a nonsmooth root-envelope point as a stationary stability conclusion.

No root index is a persistent physical phase identity. `StabilityPhase::branch` is only the current density-root slot for diagnostics.

## Failure mapping

CPA density-root results remain conservative at the flash boundary:

- `near_multiple` or root-count quota -> `root_topology`;
- density-root iteration/evaluation quota -> `root_iteration_limit`;
- no representable root / phase-property failure -> root-range/property failure evidence;
- no mechanically admissible simple root -> `no_admissible_branch`;
- nonfinite Gibbs/fugacity arithmetic -> `nonfinite_properties`.

The generic stability driver converts these provider failures to `indeterminate` unless a separate robust negative-TPD witness has already been retained. A finite CPA root scan is not promoted to global root certification, and the TPD multistart search keeps `global_stability_proven=false`.

## Validation

Current fixtures remain explicitly synthetic.

1. A non-associating pure SRK-limit state has three simple roots. The adapter must select the lower Gibbs value among the two mechanically admissible roots and preserve its `ln(phi)` exactly.
2. A low-pressure binary feed exercises the accepted reference/stationary path and verifies ordered model/provenance publication.
3. A low-temperature/high-pressure binary fixture has an independently pre-audited robust negative TPD start; the generic search must retain it as `unstable` without requiring stationarity.
4. Reversing runtime component order must preserve the TPD value and remap the witness composition by component identity.
5. Exhausting the CPA density-root evaluation budget at the feed reference must produce `indeterminate` stability with explicit root-iteration/resource evidence.

## Boundary

After this gate CPA has a valid feed/trial phase provider for the existing generic TPD algorithm. It still does not provide:

- a two-phase equilibrium solve;
- phase-fraction/material-balance solution;
- final common-tangent review of an accepted pair;
- max-three-phase orchestration;
- phase-set publication/backend integration;
- physical CPA parameter validation.

The next gate should reuse the existing generic `iterate_pt_split(...)` / `solve` contract with a CPA evaluator, verify material balance and fugacity equality, then perform the same final common-tangent review before any two-phase result is published.
