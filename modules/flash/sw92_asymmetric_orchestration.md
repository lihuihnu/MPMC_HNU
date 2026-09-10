# SW92 Xu-style Gate 3B.2 witness orchestration and candidate selection

## Scope

This increment implements only **Gate 3B.2** under

```text
SW92-equilibrium/xu-asymmetric-gibbs/v1
```

with orchestration identity

```text
SW92-equilibrium/xu-asymmetric-gibbs/max2-witness-orchestration-gibbs-selection/v1
```

It consumes the already validated Gate 3A asymmetric stability result and Gate 3B.1 family-aware
fixed-pair primitive. It builds family-aware pair attempts from robust negative TPD witnesses,
executes a family-neutral finite attempt plan, removes equivalent/slot-swapped fixed-pair candidates,
and selects the uniquely lowest resolved candidate Gibbs value when possible.

It deliberately **does not** implement Gate 3B.3 final two-family phase-set stability, pair-vs-feed
Gibbs acceptance, or `accepted_phase_set()`. A selected two-phase candidate remains explicitly
pending final stability.

The thermodynamic model remains

```text
SW92/corrected-original/PR76-base/NaCl-molality
```

with fixed NaCl molality as an external coordinate and material balance only over the ordered EOS
components.

## 1. Gate 3A witness -> fixed-pair seed contract

Gate 3A must first return `StabilityStatus::unstable` with a resolved lower-feed
`reference_family`. Each retained `Sw92AsymmetricNegativeWitness` already has an explicit family,
trial index, composition, TPD value and arithmetic guard.

For witness family `Fw`, the primary Gate 3B.2 assignment is

```text
(reference_family, Fw).
```

The witness composition is placed in phase slot 1 as a finite incipient phase. A complementary phase
slot 0 is constructed from the normalized feed using the same positive-support material-balance seed
principle as the established PT split baseline. The resulting primitive coordinate is

```text
logK_i = ln(w_i) - ln(x0_i).
```

No composition floor is introduced. A witness that cannot produce a finite positive complementary
phase remains in the plan with an explicit unrepresentable-seed diagnostic.

If `Fw != reference_family`, the same witness also generates the mandatory same-family alternative

```text
(Fw, Fw).
```

with the same material-balanced `logK` seed. This is essential: an opposite-family instability does
not prove the final two-phase minimum must retain the original lower-feed family. In particular, a
feed with AQ as the lower single-phase surface must still be able to reach an `NA+NA` candidate.

When `Fw == reference_family`, the primary assignment is already the same-family alternative and no
duplicate entry is created.

## 2. Family-neutral finite attempt order

Gate 3A stores diagnostic results by family, but that storage order must not become scientific
priority. `build_sw92_asymmetric_pair_plan(...)` therefore reconstructs an explicit attempt order.

Within each family, robust witnesses are ranked primarily by their common-reference TPD value: more
negative values represent stronger instability evidence. For each witness rank, the primary entry
from every family that has evidence is inserted **before** any same-family alternative for that
rank. Only after both family-primary entries have been represented are the required opposite-family
same-family alternatives appended.

Consequences:

- with negative evidence from both AQ and NA, the first two primary plan slots represent different
  witness families;
- one family cannot consume all early attempt slots merely because Gate 3A happened to store its
  diagnostics first;
- reversing the stored `negative_witnesses` vector does not change the semantic plan;
- deterministic tie fallbacks affect execution order only. If the configured finite budget truncates
  the plan, Gate 3B.2 withholds candidate selection entirely, so such fallback order cannot become an
  accepted scientific decision.

The public result retains the complete plan, the attempted prefix, and `attempt_limit_reached`.

## 3. Attempt budget semantics

`max_pair_attempts` limits the number of plan entries executed. Each consumed plan entry counts even
if its material-balanced seed is unrepresentable. Every actual fixed-pair call receives the same
configured `Sw92AsymmetricFixedPairOptions`; family identity does not change the per-attempt numerical
budget.

If

```text
attempts.size() < plan.size()
```

then the result is `attempt_limit_reached` and `selected_candidate_class` remains empty even if the
attempted prefix already contains an admissible pair. This prevents a partial plan from silently
promoting an order-dependent provisional minimum.

The result still retains any completed attempts/candidate classes as diagnostics.

## 4. Fixed-pair attempt semantics are inherited, not weakened

Each planned entry calls the reviewed Gate 3B.1 primitive

```text
iterate_sw92_asymmetric_fixed_pair(...)
```

with an explicit `(F0,F1)`. Gate 3B.2 does not bypass or recompute its acceptance gates. A pair can
enter candidate deduplication only when `candidate_admissible()` is true, which already requires:

- common reduced chemical-potential equations within tolerance;
- EOS-component material balance within tolerance;
- non-disappearing, compositionally distinct phases;
- same-family minimum-Gibbs root smoothness/conditioning;
- assigned family resolved lower than the opposite family at both phase compositions.

Equation-converged but dominated/tied/disappearing/resource-limited fixed-pair results remain in the
attempt record but do not compete in candidate Gibbs selection.

## 5. Equivalent-candidate and slot-swap deduplication

Different Gate 3A witnesses can converge to the same physical two-phase candidate. Phase-slot swap
also produces a numerically different representation of the same family-aware set.

Gate 3B.2 groups admissible attempts into `Sw92AsymmetricCandidateClass` objects. Two states are
considered equivalent only when either direct slot matching or complete slot swapping satisfies:

- matched phases have identical `SwPhaseFamily`;
- mole phase fractions agree within a derived numerical equivalence tolerance;
- active-support log compositions agree within the fixed-pair `log_k_separation` scale;
- both source attempts have already passed the fixed-pair equation/material-balance gates.

The phase-fraction equivalence tolerance is derived from the existing fixed-pair material-balance
thresholds and machine roundoff; it is not a new phase-acceptance threshold. The log-composition
threshold reuses the existing fixed-pair composition-distinction scale.

Within one equivalence class, the representative is chosen by numerical quality: lower
chemical-potential residual, then lower material-balance residual, then smaller arithmetic Gibbs
guard. Candidate Gibbs ranking is performed **between classes**, not between duplicate witness
solutions.

## 6. Guarded candidate Gibbs selection

Every Gate 3B.1 candidate already records

```text
G_pair/(RT)
gibbs_roundoff_guard
```

under the common AQ/NA component-reference gauge audited in Gate 3A.

After deduplication, Gate 3B.2 finds the numerically lowest candidate-class Gibbs value. If any
other distinct class differs from that minimum by no more than the sum of the two recorded arithmetic
guards, the result is

```text
candidate_gibbs_tie
```

with no selected candidate. Attempt order is not used to break a distinct near-tie.

If one class is uniquely lower beyond the combined guards and the full plan was executed, the result
is

```text
pair_candidate_selected_pending_final_stability
```

and `selected_pair_candidate_pending_final_stability()` exposes only that selected Gate 3B.1 state.
The name is deliberate: this is **not** an accepted equilibrium phase set.

Gate 3B.2 does not yet compare the selected pair Gibbs with the original lower-envelope feed Gibbs;
that remains a Gate 3B.3 acceptance check together with final asymmetric stability.

## 7. Stable and indeterminate initial Gate 3A paths

When Gate 3A returns `no_instability_found` with a resolved lower-feed family, Gate 3B.2 retains a
family-aware single-phase candidate at the normalized feed composition. No forced split is run, and
`global_stability_proven` remains false because the initial search is finite multistart.

When Gate 3A is indeterminate, including family-envelope tie, same-family root nonsmoothness or
required-search resource failure, Gate 3B.2 performs no pair planning and does not coerce the state to
single phase.

## 8. Result/publication boundary

`Sw92AsymmetricPairSelectionResult` records:

- `maximum_phase_count = 2` as an algorithm capability only;
- p/T/feed/fixed molality and SW92 dataset/component metadata;
- Gate 3A initial asymmetric stability;
- the complete family-neutral plan;
- attempted plan entries and per-attempt Gate 3B.1 diagnostics;
- total pair property-evaluation accounting;
- candidate equivalence tolerances;
- deduplicated candidate classes;
- optional uniquely selected class;
- the explicit Gate 3B.2 orchestration identity;
- `global_stability_proven=false`.

There is intentionally **no** `accepted_phase_set()` and no final asymmetric stability field in this
increment.

A pair selected here still has to pass Gate 3B.3:

1. selected-pair versus original lower-feed Gibbs check;
2. common tangent formed from the selected pair's two reduced chemical-potential vectors;
3. both phase compositions added as starts to both AQ and NA searches;
4. final two-family common-tangent stability with no robust negative witness and no indeterminate
   family search;
5. only then family-aware maximum-two-phase publication.

## 9. Validation scope

The Gate 3B.2 regression uses the already independently anchored CO2/H2O state at 3 MPa, 340 K,
fresh water and `z_CO2=0.7`.

The Gate 3B.1 Decimal(80) reference establishes that at this state:

- the `AQ+AQ` fixed-pair solution is family-admissible and has the lowest anchored pair Gibbs among
  the three explicit assignments examined;
- `AQ+NA` has an equation root but its NA phase is AQ-dominated;
- `NA+NA` has an equation root but both NA assignments are AQ-dominated.

Gate 3B.2 verifies that family-tagged Gate 3A witnesses generate the required primary and same-family
alternative attempts, that diagnostic witness-vector reversal leaves the plan invariant, that the
full finite plan selects the anchored AQ+AQ candidate, and that a truncated two-slot budget first
represents both witness families but withholds selection.

Additional structural tests verify:

- stable Gate-3A feed -> family-aware single-phase candidate with no forced split;
- exact synthetic family-envelope tie -> initial indeterminate with no pair plan;
- slot-swapped candidates -> one equivalence class;
- two structurally distinct candidate classes tied within their Gibbs guards -> no attempt-order
  selection;
- ordered-component permutation invariance;
- fixed-pair property-budget failure attribution;
- public-header self containment.

The dedicated workflow reruns the existing Gate 3A and Gate 3B.1 suites and their independent
Decimal(80) reference regeneration. No new thermodynamic parameter or experimental datum is added by
Gate 3B.2.

## Deferred work

Gate 3B.3 remains explicitly deferred:

- pair-vs-feed lower-envelope Gibbs acceptance;
- selected-pair common-reference allowance;
- candidate phase compositions as starts in both final family searches;
- final AQ/NA `unstable / indeterminate / no_instability_found` phase-set decision;
- authoritative family-aware `accepted_phase_set()`.

Three-phase, salt-inventory conservation, SW92 derivatives/implicit sensitivity, physics closure and
global interval stability also remain outside this increment.
