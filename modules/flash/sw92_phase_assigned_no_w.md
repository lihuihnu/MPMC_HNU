# SW92 Profile-C no-W topology adapter

## Scope

`sw92_phase_assigned_no_w.hpp` supplies the missing Profile-C topology path in
which no AQ-family phase is retained:

```text
H
or
H0 + H1
```

`no-W` does **not** remove water from the thermodynamic system. Water remains an
ordinary EOS component in every retained phase. All retained H phases use the
SW92 NA family with the same ordered parameter snapshot and fixed NaCl molality.

Algorithm identity:

```text
SW92-equilibrium/phase-assigned-aq-na-joint/
no-w-na-flash-targeted-aq-appearance/v1
```

This adapter is morphology-neutral: low/high-Z numerical slots are only a
deterministic representation inherited from the existing one-family split. They
are exposed as unordered H instances, never as physical liquid/vapor labels.

## Stage 1: fixed-NA phase multiplicity

The adapter reuses the already validated `solve_sw92_pt_family_vle` primitive
with `SwPhaseFamily::nonaqueous` fixed for all property calls.

That primitive already performs:

1. minimum-Gibbs-root feed stability inside NA;
2. witnessed material-balanced two-phase solve when needed;
3. common reduced chemical potentials;
4. Gibbs comparison against the feed;
5. final same-NA common-tangent finite stability review.

Consequently C2 no-W does not add a second, competing two-phase solver.

The accepted finite-search outcomes are represented as:

```text
one NA phase -> one unresolved H instance
two NA phases -> unordered H0 + H1 instances
```

A final same-NA instability is retained as
`higher_h_multiplicity_or_wrong_candidate`; it is not truncated to two phases.
Phase-disappearance/numerical/property failures remain unresolved.

## Stage 2: targeted W appearance

A no-W H state is not authoritative until the topology question

```text
should a distinct water-rich AQ phase appear?
```

has been tested.

The adapter reconstructs the retained H common tangent and evaluates AQ trials
against that tangent. The default AQ search is intentionally **targeted**, not a
whole-simplex automatic AQ/NA competition.

The mandatory start set is:

- one dynamically constructed water-richer trial on the same active component
  support;
- every retained H composition evaluated under AQ only as a diagnostic;
- any explicit caller starts that pass the same support/resource contract.

The targeted water start uses a water fraction halfway between the maximum
retained-H water fraction and one, with the remaining amount distributed over
positive non-water feed components in proportion to their feed inventory. If a
numerically distinct water-richer start cannot be represented, the result stays
`single_phase_role_unresolved` rather than inventing a phase-domain cutoff.

If total water inventory is exactly zero, an AQ phase cannot appear by closed
nonreactive material balance, so no AQ search is required.

## W-witness acceptance

A robustly negative AQ TPD trial becomes a usable W seed only when it is:

1. compositionally distinct from **every** retained H phase on active support;
2. water-richer than the most water-rich retained H phase by a roundoff-scaled
   guard.

A negative AQ trial at an existing H composition is retained as model-extension
diagnostics but cannot become a physical W phase. This is the same role-vs-model
separation used by the earlier H-side NA witness gate.

Outcomes:

```text
no_w_single_h_locally_closed
no_w_two_h_locally_closed
aqueous_phase_witness_found
higher_h_multiplicity_or_wrong_candidate
phase_disappearance_unresolved
single_phase_role_unresolved
indeterminate
```

All finite-search closure states retain:

```text
global_stability_proven = false
accepted_phase_set_published = false
morphology_resolved = false
```

## Why this is not Profile B

The adapter never asks which of `g_AQ(x)` and `g_NA(x)` is globally lower at an
arbitrary composition. NA is fixed for retained H instances; AQ is evaluated
only as a candidate **new physical water-rich phase**. A role-inadmissible AQ
negative trial cannot reassign an H phase to AQ.

## Boundary semantics

The adapter is suitable both as an autonomous no-W candidate generator and as
the destination of C2b.2 `route_to_no_w_h0_h1`.

- zero water inventory closes W appearance by material balance;
- pure/near-pure water support where a distinct W-vs-H topology cannot be
  represented remains role-unresolved;
- an accepted H pair whose final same-NA review finds another instability is
  not published under the maximum-two-H adapter;
- phase disappearance remains unresolved rather than being silently converted
  to a lower phase count.

## Validation boundary

The implementation reuses the independent Decimal(80) one-family NA anchors and
adds a standalone Decimal cross-check of targeted AQ appearance relative to the
NA common tangent. Existing C1 reference states are used only as independent
water-rich trial compositions; no production result is imported.

Focused C++ tests cover:

- wet CO2/H2O no-W candidate producing a targeted W witness;
- dry H-only finite closure;
- zero-water material-balance closure;
- runtime component permutation;
- witness role/distinction guards;
- resource and root/property failures;
- public-header self containment.

This gate still does not implement PIP/LV/LL morphology, authoritative phase-set
publication, salt inventory, derivatives or physics coupling.
