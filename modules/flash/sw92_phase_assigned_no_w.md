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
ordinary EOS component in every NA phase evaluation. All retained numerical
candidates use the SW92 NA family with the same ordered parameter snapshot and
fixed NaCl molality.

Algorithm identity:

```text
SW92-equilibrium/phase-assigned-aq-na-joint/
no-w-na-flash-targeted-aq-appearance/v1
```

This adapter is morphology-neutral. More strongly, the fixed-NA phases are
initially **NA numerical candidates**, not automatically authoritative physical
H roles. This distinction is necessary because the NA model can mathematically
produce a very water-rich branch that a later physical-role topology solve must
replace rather than mislabel as hydrocarbon phase.

## Stage 1: fixed-NA phase multiplicity

The adapter reuses the already validated `solve_sw92_pt_family_vle` primitive
with `SwPhaseFamily::nonaqueous` fixed for all property calls.

That primitive already performs:

1. minimum-Gibbs-root feed stability inside NA;
2. witnessed material-balanced two-phase solve when needed;
3. common reduced chemical potentials;
4. Gibbs comparison against the feed;
5. final same-NA common-tangent finite stability review.

Consequently the no-W adapter does not add a second competing two-phase solver.
Low/high-Z split slots remain deterministic numerical representation only; they
are never liquid/vapor labels.

A final same-NA instability is retained as
`higher_h_multiplicity_or_wrong_candidate`; it is not truncated to two phases.
Phase-disappearance/numerical/property failures remain unresolved.

## Stage 2: targeted W appearance

A fixed-NA candidate state is not a final physical no-W topology until the
question

```text
can a distinct water-rich AQ phase be added?
```

has been tested.

The adapter reconstructs the retained NA common tangent and evaluates AQ trials
against that tangent. The default AQ search is intentionally **targeted**, not a
whole-simplex automatic AQ/NA competition.

The start set is ordered deliberately:

1. one dynamically constructed water-targeted composition on the same active
   component support — the **only** candidate-generating W start in v1;
2. each retained NA candidate composition under AQ — diagnostics only;
3. caller-supplied starts — diagnostics only in v1.

This start provenance is part of the physical-role guard. A negative AQ descent
originating from an H-like/NA-like seed is model-extension evidence and cannot
become a physical W phase merely by moving slightly toward water.

The targeted water start uses a water fraction halfway between the maximum
retained-NA-candidate water fraction and one, with the remaining amount
distributed over positive non-water feed components in proportion to feed
inventory. This is initialization only; it is not an absolute physical phase
cutoff. If such a distinct start cannot be represented, the result remains
`single_phase_role_unresolved`.

If total water inventory is exactly zero, an AQ phase cannot appear by closed
nonreactive material balance, so no AQ search is required.

## Why W is not required to be richer than every fixed-NA candidate

Independent Decimal evidence at the audited binary CO2/H2O state, `p=3 MPa`,
`T=340 K`, shows:

```text
fixed-NA low-Z candidate water  ~= 0.999716
Profile-C AQ-W water            ~= 0.994055
fixed-NA high-Z candidate water ~= 0.011274
```

Thus a fixed-NA mathematical pair can contain a W-like candidate even though
both property evaluations used `SwPhaseFamily::nonaqueous`. Requiring an AQ
trial to be water-richer than **every** fixed-NA candidate would reject the
known joint Profile-C W/H solution for the wrong reason.

Therefore W-candidate generation requires:

1. robust negative AQ TPD from the explicit water-targeted start;
2. composition distinct from every retained NA candidate;
3. water fraction greater than the **water-poorest** retained NA candidate by a
   roundoff-scaled guard.

This is deliberately weaker than final physical-role acceptance. Once a W seed
enters the joint `W+H` or `W+H0+H1` solver, the existing Profile-C contract again
requires W to be water-richer than **all final retained H phases**. Any W-like NA
mathematical candidate that cannot survive that joint solve is discarded.

This separation keeps candidate generation permissive enough to find the right
topology while keeping final role acceptance strict.

## Outcomes

```text
no_w_single_h_locally_closed
no_w_two_h_locally_closed
aqueous_phase_witness_found
higher_h_multiplicity_or_wrong_candidate
phase_disappearance_unresolved
single_phase_role_unresolved
indeterminate
```

`no_w_*_locally_closed` means the fixed-NA multiplicity and the **targeted** W
appearance start are locally closed under the declared finite searches. AQ
negative trials originating from retained NA diagnostics do not prevent that
closure because they are forbidden model-reassignment evidence, not physical W
candidate generation.

All finite-search states retain:

```text
global_stability_proven = false
accepted_phase_set_published = false
morphology_resolved = false
```

## Why this is not Profile B

The adapter never asks which of `g_AQ(x)` and `g_NA(x)` is globally lower at an
arbitrary composition. NA is fixed for the existing numerical state; AQ is
used only as a candidate **new water-rich physical phase** from explicit W-like
initialization. An AQ negative descended from an NA/H-like diagnostic start
cannot reassign that phase to AQ.

## Boundary semantics

The adapter is suitable both as an autonomous no-W candidate generator and as
the destination of C2b.2 `route_to_no_w_h0_h1`.

- zero water inventory closes W appearance by material balance;
- pure/near-pure support where a distinct W-vs-H topology cannot be represented
  remains role-unresolved;
- a fixed-NA pair whose final same-family review finds another instability is
  not published under the maximum-two-candidate adapter;
- phase disappearance remains unresolved rather than being silently converted
  to a lower phase count.

## Validation boundary

The implementation reuses the independent Decimal(80) one-family NA anchors and
adds a standalone Decimal cross-check. For the wet CO2/H2O state the independent
reference verifies both:

- negative AQ TPD relative to the fixed-NA pair common tangent;
- the fixed-NA pair itself contains a water-richer mathematical branch than the
  independently solved AQ-W state, proving why `NA family == physical H` is not
  an admissible shortcut.

For dry `z_CO2=0.995`, the same independent AQ-W composition has positive TPD
against the minimum-Gibbs NA feed tangent. These are model numerical checks,
not experimental validation or global-stability proofs.

Focused C++ tests cover:

- wet CO2/H2O targeted W-witness generation;
- dry one-NA/H candidate finite closure;
- zero-water material-balance closure;
- runtime component permutation;
- witness provenance/role/distinction guards;
- resource and root/property failures;
- public-header self containment.

This gate still does not implement PIP/LV/LL morphology, authoritative phase-set
publication, salt inventory, derivatives or physics coupling.
