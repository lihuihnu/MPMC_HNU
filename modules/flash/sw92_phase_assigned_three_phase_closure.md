# SW92 Profile-C C2b.2 W-present topology closure

## Scope

`sw92_phase_assigned_three_phase_closure.hpp` implements the post-C2b.1
**W-present nested-topology review** for

```text
C1 W(AQ)+H(NA)
  -> C2a1 additional-NA witness
  -> C2b.1 W(AQ)+H0(NA)+H1(NA) candidate
  -> C2b.2 W-present H-multiplicity review
```

Algorithm identity:

```text
SW92-equilibrium/phase-assigned-aq-na-joint/
w-present-c2b1-nested-topology-na-review/v1
```

C2b.2 is deliberately morphology-neutral. It never calls H0/H1 liquid or
vapor, never publishes `accepted_phase_set()`, and keeps
`global_stability_proven=false`.

## Fresh source-chain revalidation

The public review consumes matching C1, C2a1 and C2b.1 results from one ordered
SW92 snapshot. It does not trust cached status enums alone.

It verifies:

- dataset/revision/component ordering, p, T and fixed NaCl molality;
- the exact C2a1 witness selected by the C2b.1 source adapter;
- C2b.1 initial H0/H1 logK seeds and source-fraction split provenance;
- C1 material balance, relative-water topology and fresh AQ/NA minimum-Gibbs
  root properties;
- C2b.1 material balance and fresh AQ/NA/NA minimum-Gibbs root properties;
- common reduced chemical potentials across W/H0/H1;
- positive/simplex phase fractions or explicit disappearance boundaries;
- W water-richer than both retained H compositions;
- pairwise composition distinction for an interior three-phase candidate.

Fresh root/property failures and same-family nonsmoothness remain explicit
`indeterminate` evidence.

## Nested Gibbs consistency

C1 and C2b.1 reduced Gibbs values are recomputed from the fresh selected roots.
Because the transition keeps the same physical model assignment

```text
W -> AQ
H, H0, H1 -> NA
```

C2b.2 requires

```text
G_C2b1 <= G_C1 + arithmetic_roundoff_guard
```

for the nested phase-addition path. This comparison is local to the retained
W-present topology chain. It is not a feed-level AQ/NA family competition and
must not be reused to choose W-only versus H-only.

## Final NA-only review

For an interior C2b.1 candidate, C2b.2 reconstructs the three-phase common
log-activity tangent and runs the existing generic finite TPD driver with a
fixed SW92 NA evaluator.

Mandatory starts include:

1. retained W composition evaluated under NA, as a role-guard diagnostic;
2. retained H0 composition;
3. retained H1 composition;
4. the selected original C2a1 witness composition.

Automatic generic starts and caller-supplied starts may be added within the
existing resource quotas.

The TPD tolerance is enlarged only by the explicitly recomputed common-reference
arithmetic/equation allowance. No SW92 thermodynamic tolerance is relaxed.

## Additional-H witness contract

A robustly negative final NA trial blocks the retained three-phase candidate
only when all three conditions hold:

```text
trial distinct from H0
trial distinct from H1
trial water-poorer than W by a roundoff-guarded margin
```

Negative trials that are equivalent to H0/H1 are retained as diagnostics but
do not invent an extra phase. Negative W-like trials are also retained but fail
the physical-role guard.

If a usable additional-H witness remains, C2b.2 returns
`higher_phase_count_witness_found`. Under the current maximum-three-phase scope
this means the candidate is not accepted: the evidence may indicate a higher
phase count or a wrong/local candidate.

If every finite trial terminates as stationary or robust-negative evidence and
none is a usable additional-H witness, the result is only

```text
w_present_h_multiplicity_locally_closed
```

This is finite-search evidence, not global phase-number proof.

## Disappearance routing

C2b.2 does not silently convert a disappearing phase into an accepted lower
phase count.

- one H disappears -> `route_to_w_h`, returning to the C1 topology path;
- W disappears while both H phases remain -> `route_to_no_w_h0_h1`;
- W and an H disappear together -> `single_phase_endpoint_unresolved`.

The no-W and autonomous single-phase topology paths are separate future gates.

## Scientific boundary

C2b.2 intentionally does **not** implement:

- PIP or any L/V morphology classifier;
- relative/absolute Z as LV-vs-LL evidence;
- no-W H/H0+H1/H0+H1+H2 orchestration;
- autonomous W-only versus H-only role selection;
- authoritative Profile-C phase-set publication;
- salt-inventory conservation, derivatives or physics coupling.

This keeps the shortest path toward a correct three-phase flash: first close
the thermodynamic/topology graph, then expose the resulting phase instances to
higher layers without forcing unsupported morphology metadata.

## Validation

`tests/flash/sw92_phase_assigned_three_phase_closure/reference_decimal.py` is a
stdlib-only Decimal(80) structural reference. It reuses the already-independent
C2b.1 equation transcription and verifies that the nested three-phase state has
lower reduced Gibbs than its metastable C1 source, and that H0/H1 and the old
C2a1 witness lie on the final three-phase common tangent.

The underlying CH4/CO2 non-water `kij=0` remains explicitly synthetic test data;
there is no experimental three-phase claim.

Focused C++ regression covers:

- W-present locally closed review;
- W-like/duplicate/additional-H witness classification guards;
- H disappearance, W disappearance and one-phase endpoint routing;
- runtime component permutation;
- finite-search resource exhaustion and fresh root/property failure;
- source-chain mismatch/provenance guards;
- public-header self containment.

The dedicated hosted workflow also reruns affected C2b.1, C2a1 and C1 tests and
all independent Decimal references in that dependency chain.
