# CPA association continuation production safety contract v1

## Status

This document freezes the **minimum safety contract that must be satisfied before
association-state continuation may be enabled in production**.

Contract identifier:

```text
MPMC_HNU/CPA/association-continuation-safety/v1
```

This contract is design-only. It does **not** enable continuation in the current
production CPA association solver, density-root search, stability search, or PT
flash solver.

The contract is motivated by the completed test-only continuity audit on the
validated 333.15 K methanol-water five-state workload. That audit established a
credible performance opportunity, but also showed that tiny differences in
association site fractions can be amplified in the pressure observable at dense
states. Production reuse therefore must remain advisory and fail closed.

## Normative invariants

A production implementation conforming to v1 MUST preserve all of the following.

### 1. Cache ownership is exactly one density-root search

The continuation cache belongs to exactly one invocation of
`CpaPtPhase::roots(...)`.

It MUST be created empty on entry to one root search and destroyed or cleared on
exit from that same root search.

Cache state MUST NOT cross any of the following boundaries:

- another `CpaPtPhase::roots(...)` invocation;
- another composition, even at identical `T` and `P`;
- another pressure or temperature state;
- another model/parameter snapshot;
- another thread/task/root-search instance;
- stability and split calls that happen to request numerically similar states.

No global, static, model-owned, session-owned, or cross-flash association cache is
permitted by v1.

This ownership rule is both a correctness rule and a concurrency rule: two root
searches must never observe or mutate each other's continuation state.

### 2. Seed selection is nearest density within that cache only

For a new reduced density `u`, the only admissible warm seed under v1 is the
converged site-fraction state associated with the cached density nearest to `u`
inside the same root-search cache.

If the cache is empty, the association solve MUST use the current production cold
initialization.

An exact-density cache hit may reuse that exact state, subject to the same
acceptance rules as every other warm attempt.

No extrapolation, interpolation, cross-composition transport, cross-root-search
reuse, or learned predictor is part of v1.

Nearest-density selection is defined in density space rather than execution order
because bisection may jump back into a previously scanned bracket.

### 3. The fixed-point equations are unchanged

Continuation may change only the initial vector supplied to the existing
association fixed-point iteration.

The following MUST remain bit-for-bit identical in mathematical definition to the
current production cold solver:

- association-strength equations;
- site topology and multiplicities;
- mixing/cross-association parameters;
- simplified RDF definition;
- update equation;
- damping rule;
- iteration ordering;
- representability checks;
- exceptional/failure semantics.

Continuation MUST NOT introduce a second approximate association model.

### 4. The production convergence tolerance remains `1e-12`

The existing association site-fraction stopping tolerance remains unchanged.

```text
site_fraction_tolerance = 1e-12
```

Warm starts MUST NOT obtain acceptance through a looser stopping condition,
reduced iteration accuracy, pressure-tuned tolerance, or any tolerance that
varies with whether a seed was available.

The existing iteration limit also remains authoritative. A warm attempt that
hits the existing failure semantics is not converted into success.

### 5. Continuation is advisory, never authoritative

A warm association solve produces a **candidate state only**.

Convergence of the fixed-point iteration is necessary but not sufficient for that
candidate to replace the cold path.

The candidate MUST pass an observable-based acceptance certificate before its
site fractions can be used as the accepted density-state association result or
inserted into the cache as an accepted continuation state.

If the candidate cannot produce a valid certificate, the implementation MUST
fall back deterministically to the current cold association solve for that same
density state.

### 6. Acceptance must be observable-based

The continuity audit demonstrated that checking `X_A` alone is insufficient:
very small site-fraction differences can be magnified in the association-pressure
term at dense states.

Therefore v1 requires an acceptance certificate that bounds the thermodynamic
observables affected by association reuse.

At minimum, a future production design MUST certify that the warm candidate
cannot violate the already frozen CPA numerical-parity scales relevant to the
phase kernel, including:

```text
total-pressure error scale <= 5e-6 Pa
accepted-root ln(phi) error scale <= 1e-10
```

These values reference the existing ThermoPack parity v1 envelope. They are not
new physical tolerances and MUST NOT be relaxed merely to make continuation pass.

The certificate MAY be stricter.

### 7. Production acceptance may not require an unconditional cold reference solve

The test-only continuity audit could compare warm and cold results directly
because its purpose was characterization. A production fast path cannot require a
full cold association solve for every warm candidate and still claim the measured
continuation saving.

Therefore the production acceptance certificate MUST be computable without an
unconditional full cold reference solve at the same density.

Acceptable future mechanisms include only mechanisms that are independently
shown to be conservative, for example an a-posteriori bound derived from the warm
fixed-point residual/contraction and the sensitivity of the affected observable.
This document does **not** yet bless any particular formula.

Until such a no-cold-reference certificate is derived, reviewed, and validated,
production continuation MUST remain disabled.

A conditional cold solve remains mandatory whenever the certificate cannot be
established.

### 8. Cold fallback is deterministic and exact

Fallback MUST invoke the current production cold association path with the same:

- all-ones initialization;
- equations;
- options;
- `1e-12` site-fraction tolerance;
- iteration budget;
- error/status semantics.

Fallback is not allowed to use a second relaxed warm attempt.

If the cold solve fails, the caller MUST observe the same failure that the current
production implementation would have produced without continuation.

### 9. Only accepted states may advance the cache

After one density evaluation:

- an accepted warm candidate may be inserted into the local cache;
- a rejected warm candidate MUST NOT be inserted;
- after fallback, only the successfully converged cold result may be inserted;
- a failed cold result inserts nothing.

This prevents a rejected or numerically suspect warm state from contaminating
subsequent density evaluations.

### 10. Root topology remains controlled by the accepted thermodynamic result

Continuation MUST NOT change:

- density scan locations;
- scan interval count;
- pressure-sign classification;
- bracket construction;
- bisection ordering;
- root merging;
- root slope classification;
- pressure/root tolerances;
- admissible-root selection.

A density-state pressure used by root topology must come only from an association
state accepted under this contract (warm candidate with valid certificate, or the
cold fallback result).

### 11. No cross-composition cache reuse

Composition identity is a hard cache boundary, not a proximity criterion.

Two composition vectors that are numerically close are still different root
searches/states and MUST NOT share continuation state under v1.

The same rule applies to component ordering and parameter revision. Runtime
component reordering or replacement invalidates all prior association cache state.

### 12. Thread safety follows ownership, not locking a shared cache

The preferred v1 design is stack/local ownership by one root-search invocation.
No shared mutable cache is introduced, so no global lock is needed and no
cross-thread ordering becomes part of numerical behavior.

An implementation that requires shared mutable cache state is non-conforming to
v1 even if protected by a mutex.

## Required state carried by one cache entry

A conforming v1 cache entry must contain enough identity to make accidental reuse
across incompatible site layouts impossible. At minimum:

```text
reduced_density
converged site fractions in production site ordering
site-layout identity/size sufficient for defensive validation
```

The cache itself inherits `T`, composition, parameter set, and association options
from its owning root-search invocation; those values are not global lookup keys.

A future implementation may store additional diagnostics, but those diagnostics
must not alter the scientific equations.

## Required warm-attempt decision sequence

The normative sequence for one density evaluation is:

```text
1. locate nearest accepted density state in this root-search cache
2. if none exists -> run current cold association solve
3. otherwise run the unchanged association fixed-point equations from that seed
4. if warm solve does not converge -> cold fallback
5. compute the no-cold-reference observable acceptance certificate
6. if certificate fails or is unavailable -> cold fallback
7. if warm candidate is accepted -> use it and cache it
8. if fallback succeeds -> use/caches only the cold result
9. if fallback fails -> preserve current production failure semantics
```

There is no third success path.

## Required diagnostics before production merge

The first production implementation must expose test-visible counters sufficient
to prove, without changing scientific results:

```text
warm_attempts
warm_accepts
cold_fallbacks
warm_nonconvergence_fallbacks
certificate_fallbacks
exact_density_hits
association_sweeps_warm
association_sweeps_cold_fallback
```

These counters may be test instrumentation or returned audit diagnostics; they
must not become externally visible frontend/API requirements.

## Acceptance gates for a future production implementation

A production implementation is not acceptable merely because it is faster. It
must simultaneously demonstrate:

1. identical root-search status/topology on the fixed audit workload;
2. unchanged existing physical-validation outcomes;
3. unchanged material-balance, fugacity, TPD and common-tangent gates;
4. unchanged failure/indeterminate semantics;
5. all frozen ThermoPack parity v1 thresholds still pass;
6. deterministic fallback behavior under deliberately bad/unusable warm seeds;
7. cache isolation across composition, root search and concurrent calls;
8. a repeated Release performance improvement against the frozen baseline after
   charging fallback cost.

No gate may be weakened in the same change that introduces continuation.

## Evidence from the pre-production continuity audit

The test-only audit that motivated this contract covered the full fixed
five-state workload:

```text
root searches                     = 497
density evaluations               = 305,928
seeded evaluations                = 305,431
accepted warm continuations       = 297,415
fallbacks                         = 8,016
fallback rate                     = 2.6244880186%
nonconvergence fallbacks          = 0
site-fraction guard fallbacks     = 0
pressure guard fallbacks          = 8,016
accepted-root ln(phi) fallbacks   = 0
```

With every rejected warm attempt charged for both the attempted warm work and a
full cold fallback, fixed-point sweeps decreased from `12,548,901` to
`7,068,639`, an effective reduction of `43.671250574%` on that workload.

This evidence shows that nearest-density continuation is worth designing. It does
not waive any requirement above.

## Explicit non-goals of v1

The following are outside this contract and require a new design/version:

- cross-root-search or cross-flash cache reuse;
- interpolation/extrapolation of site fractions;
- adaptive density scan reduction;
- relaxed association/root tolerances;
- changing damping based on continuation state;
- learned or heuristic seed predictors;
- skipping cold fallback after certificate failure;
- three-phase/max3-specific shared continuation state.

## Current blocker before production implementation

The remaining design blocker is precise and intentional:

> derive and independently validate a conservative **no-cold-reference observable
> acceptance certificate** for warm association candidates.

The previous shadow audit used direct warm-vs-cold pressure comparison to decide
fallback. That is valid for auditing but cannot be the production fast path.

Until the certificate exists, this v1 contract is frozen and continuation remains
an audit-only capability.
