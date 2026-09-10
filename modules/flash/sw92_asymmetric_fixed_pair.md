# SW92 Xu-style Gate 3B.1 family-aware fixed-pair primitive

## Scope

This increment implements only the **Gate 3B.1 fixed-family-pair numerical primitive** under the
already reserved parent equilibrium profile

```text
SW92-equilibrium/xu-asymmetric-gibbs/v1
```

with numerical primitive identity

```text
SW92-equilibrium/xu-asymmetric-gibbs/fixed-pair-logK-SSI-RR/v1
```

It does **not** perform Gate-3A-witness orchestration, compare/deduplicate multiple family-pair
candidates, select an overall joint phase set, perform final two-family phase-set stability, or
publish an accepted maximum-two-phase equilibrium. Those remain Gate 3B.2/3B.3 work.

The thermodynamic profile remains

```text
SW92/corrected-original/PR76-base/NaCl-molality
```

and NaCl molality remains one externally fixed coordinate. Only the ordered EOS-component vector is
material-balanced.

## Family-aware phase contract

A fixed-pair call receives an explicit ordered pair

```text
(F0, F1), F0,F1 in {AQ,NA}.
```

Repeated families are legal. `AQ+AQ`, `AQ+NA`, `NA+AQ`, and `NA+NA` are numerical assignments; slot
0/1 order is not a physical liquid/vapor label. Gate 3B.2 will later remove slot-swapped duplicate
representations during candidate orchestration.

Each returned phase stores:

- explicit `SwPhaseFamily`;
- mole phase fraction;
- EOS-component composition;
- selected same-family minimum-Gibbs `StabilityPhase` including cubic-root branch;
- compressibility factor Z as property/diagnostic metadata.

Z and root index never determine AQ/NA family identity.

## Behavior-preserving SW92 stability extension

`Sw92FamilyStabilityEvaluator` now exposes `evaluate_selected(...)`, returning
`Sw92FamilySelectedPhase { activity, compressibility_factor }` from the **same single
minimum-Gibbs-root selection** already used by `operator()`.

The original `operator()` still returns `StabilityPhase` and delegates to that selection. This avoids
re-solving the cubic merely to recover Z and avoids duplicating the tested same-family minimum-Gibbs
root/tie/conditioning logic. Existing SW92 stability, one-family VLE and Gate 3A behavior must remain
unchanged under regression.

## Reusable imposed-common-tangent asymmetric stability helper

Gate 3A now exposes

```text
test_sw92_pt_asymmetric_stability_against(...)
```

for two required AQ/NA finite TPD searches against exactly one caller-supplied reduced tangent. The
helper owns:

- preflight of both per-family search contracts before thermodynamic work;
- AQ/NA evaluator construction;
- exact propagation of the same `common_log_activity` to both searches;
- independent family resource/evaluation records;
- family-tagged robust negative witnesses;
- combined `unstable / no_instability_found / indeterminate` semantics;
- `global_stability_proven=false`.

The existing Gate-3A feed-reference entry reuses this helper after constructing its lower-family feed
tangent. The generic `pt_stability.hpp` is unchanged.

This helper is not itself a phase-set stability decision; Gate 3B.3 will later call it with a tangent
formed from a selected joint pair and with both candidate phase compositions supplied as starts to
both family searches.

## Fixed-pair material balance

For active feed support, the primitive defines

```text
logK_i = ln(x_i^1 / x_i^0)
```

and uses the existing Rachford-Rice kernel privately:

```text
x_i^0 = z_i / (1-beta + beta*K_i)
x_i^1 = K_i*x_i^0
sum_i z_i*(K_i-1)/(1-beta+beta*K_i) = 0.
```

The old `RachfordRiceResult` field names `liquid`, `vapor`, and `vapor_fraction` are immediately
mapped internally to phase0 composition, phase1 composition, and phase1 fraction. They are not
exposed as physical labels by the new public state.

The input feed is normalized only within the repository's existing roundoff-only composition
contract. No positive floor is added; zero-overall-feed components remain outside the active
support.

## Fixed-family chemical-potential equations

At one material-balanced state:

```text
m_i^0 = ln(x_i^0) + ln(phi_i^F0(x^0))
m_i^1 = ln(x_i^1) + ln(phi_i^F1(x^1))
r_i   = m_i^0 - m_i^1.
```

The fixed-pair equations reach their numerical gate only when

```text
max_i |r_i| <= chemical_potential_tolerance
```

for all active components.

Each family evaluation uses the mechanically admissible **same-family minimum-Gibbs cubic root**.
A same-family root-envelope tie is nonsmooth and cannot be converted to a low-Z/high-Z branch
choice.

## Log-K iteration and Gibbs line search

With the family assignment fixed during an attempt, the same component-inventory variation used by
the existing PT-VLE derivation gives

```text
d(G/RT) = -sum_i r_i d(v_i)
```

for phase-1 component amounts `v_i`. The Rachford-Rice kinematic relation is independent of which
SW92 family supplies each phase. Therefore the baseline direction is

```text
logK <- logK + alpha * r / scale
```

with the same principles as the existing solver:

- cap the log step;
- require resolved Gibbs decrease when the predicted decrease exceeds arithmetic uncertainty;
- otherwise require representable residual progress;
- reject line-search property/RR failures without silently clipping compositions;
- keep property and iteration budgets explicit.

The option type is new and family-aware. It intentionally has no relative-Z phase-distinction
criterion.

## Equation gate ordering

When the chemical-potential residual first satisfies its tolerance, the primitive checks in order:

1. EOS-component material balance;
2. phase disappearance threshold;
3. composition distinction using active-support `logK` contrast;
4. lower-envelope family assignment for both phases.

A phase at/below `minimum_phase_fraction` is reported as `phase_disappearance` only after the balance
and chemical-potential gates pass. It is not converted to accepted single-phase evidence.

## Mandatory lower-envelope family assignment check

A fixed-family equation root is not sufficient. At each converged phase composition the primitive
also evaluates the **opposite** SW92 family using its same-family minimum-Gibbs root and computes

```text
Delta_g/(RT)
  = sum_i x_i [ln(phi_i^AQ) - ln(phi_i^NA)].
```

The same arithmetic guard principle used by Gate 3A is applied.

For each phase:

- assigned family resolved lower -> `assigned_lower`;
- opposite family resolved lower -> `assigned_dominated`;
- AQ/NA family Gibbs tie -> `family_tie`;
- opposite-family minimum-root nonsmooth -> `opposite_family_nonsmooth`;
- property failure -> explicit failure.

The fixed-pair result is `converged` only if **both** phase assignments are resolved lower. A
dominated/nonsmooth pair retains its equation-converged point for diagnostics/reference but
`equations_converged()` returns false because the pair is not an admissible lower-envelope candidate.
The primitive never switches family in place after detecting dominance.

## Result status is not overall equilibrium acceptance

`Sw92AsymmetricFixedPairStatus::converged` means only that this one explicit family assignment passed:

- interior Rachford-Rice material balance;
- common reduced chemical-potential equations;
- phase-fraction and composition-distinction gates;
- same-family root smoothness/conditioning;
- per-phase lower-envelope family-dominance checks.

It does **not** mean the pair is the lowest-Gibbs candidate among all assignments, stable to all AQ/NA
trials, or accepted as the maximum-two-phase equilibrium. There is deliberately no
`accepted_phase_set()` in Gate 3B.1.

## Independent numerical references

The dedicated `tests/flash/sw92_asymmetric_fixed_pair/reference_decimal.py` uses Python standard-library
`Decimal(80)` only and imports no production code. It independently rebuilds:

- corrected SW92 water alpha and CO2/water AQ BIP;
- PR pure and mixture coefficients;
- cubic roots and mechanically admissible same-family minimum-Gibbs root selection;
- fixed-family-pair common-chemical-potential equations;
- material-balanced phase fraction for the repository CO2/H2O feed;
- total reduced pair Gibbs;
- AQ-minus-NA family Gibbs ordering at each converged phase.

At `p=3 MPa`, `T=340 K`, fresh water, `z_CO2=0.7`, three explicit assignments are anchored:

- `AQ+AQ`: equation-converged and both phases are lower-envelope AQ at their compositions;
- `AQ+NA`: fixed-pair equations have an interior root, but the NA phase is AQ-dominated;
- `NA+NA`: fixed-pair equations have an interior root, but both assigned NA phases are AQ-dominated.

This distinction deliberately verifies that solving fixed-pair equations is not treated as an
accepted asymmetric phase set.

Additional regressions cover:

- rough-seed SSI/Gibbs-line-search convergence for `AQ+AQ`;
- imposed-common-tangent helper equivalence with Gate 3A;
- synthetic AQ/NA family tie blocking candidate promotion;
- phase disappearance only after equation/balance convergence;
- component permutation;
- explicit property/dominance evaluation budgets and root failure;
- public-header self containment.

These are model numerical references, not experimental validation and not a global equilibrium
proof.

## Deferred Gate 3B work

### Gate 3B.2

Still deferred:

- Gate-3A family-tagged witness -> family-pair attempt planning;
- required same-family alternatives for opposite-family witnesses;
- resource-neutral ordering across AQ/NA witnesses;
- slot-swap/equivalent-candidate deduplication;
- guarded Gibbs ranking across distinct fixed-pair candidates;
- stable single-phase top-level flow.

### Gate 3B.3

Still deferred:

- selected-pair common tangent as the final phase-set reference;
- both candidate compositions as starts in both AQ and NA final searches;
- pair-vs-lower-feed Gibbs comparison;
- final two-family stability acceptance/rejection;
- authoritative family-aware `accepted_phase_set()`.

Only after Gate 3B.3 passes review may the repository claim a maximum-two-phase jointly
material-balanced SW92/Xu candidate phase set. Even then finite multistart stability remains
`global_stability_proven=false`, and fixed molality still does not conserve salt inventory.
