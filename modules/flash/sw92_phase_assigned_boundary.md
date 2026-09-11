# SW92 Profile-C three-phase disappearance re-solve

## Purpose

`sw92_phase_assigned_boundary.hpp` closes a deliberate gap between C2b.1/C2b.2
and an authoritative maximum-three-phase flash.

A C2b.1 state can satisfy material balance and chemical-potential equality while
one phase fraction is at or below `minimum_phase_fraction`. C2b.2 correctly
classifies such states as a route to a neighboring topology, but the existing
`solve_sw92_phase_assigned_pt(...)` deliberately stops at `topology_unresolved`.
It must not simply delete the small phase.

This adapter performs the missing operation: **re-solve the neighboring
physical topology from the retained boundary compositions and re-run the
neighbor-specific stability/topology review**.

## Source-chain provenance

Boundary re-solve is allowed only for one consistent retained chain. The public
resolver verifies the C1/C2a1/C2b.1/C2b.2 p, T, feed, fixed NaCl molality,
dataset/revision and ordered component snapshot. It also reuses the existing
C2a1-to-C2b.1 provenance checks: the selected C2a1 witness index and the C2b.1
initial H0/H1 logK/fraction seeds must still identify the same phase-addition
path. A mismatched C2a1 result is therefore `source_chain_inconsistent` before
any neighboring solve is attempted.

## Rules

### H disappearance: `W+H0+H1 -> W+H`

The surviving W(AQ) and H(NA) compositions seed a fresh C1 joint solve using the
same p/T/feed/model snapshot and fixed NaCl molality. The disappearing H
composition is supplied as an explicit C2a1 diagnostic start.

The lower topology closes only when C1 equations/material balance/role guard
converge again and C2a1 finds no robust additional-NA witness. If the
disappearing H still produces a robust negative witness, the phase is **not**
dropped even if the original C2b.1 fraction was below the numerical threshold.

### W disappearance: `W+H0+H1 -> H0+H1`

Both retained H(NA) compositions are supplied to a fresh no-W solve as initial
and final stability starts. The lower topology closes only when the no-W solver
itself reports locally closed one/two-H topology. A renewed W witness,
higher-H evidence, root/property failure or unresolved disappearance prevents
phase deletion.

For the physical Mortezazadeh-Rasaei Sample-6 H0+H1 geometric edge, the fresh
no-W solve still finds a robust targeted AQ-W appearance witness. The edge is a
triple-line boundary with an incipient W phase, not evidence that W can be
removed and H0+H1 published as an independently closed no-W state.

### Simultaneous endpoint

When W and an H disappear together, the surviving H composition seeds the same
fresh no-W path. The endpoint is accepted only if that no-W path closes; there
is no direct `three phases -> one phase` clipping rule.

## Boundary-aware wrapper

`solve_sw92_phase_assigned_pt_boundary_aware(...)` preserves all interior
behavior of the existing top-level driver. If the old driver stops after C2a1
because C2b.1 reaches disappearance, the wrapper re-runs bounded C2b.1 witness
attempts, retains only equation-converged disappearance states satisfying the
nested Gibbs guard, obtains an explicit C2b.2 route, and re-solves the requested
neighboring topology. A lower phase count is returned only when that neighbor
review closes.

Specific higher-phase evidence from the fresh neighbor is not collapsed to a
generic unresolved status. If a no-W fixed-NA maximum-two-phase state still has
a further same-family instability, the wrapper propagates
`higher_phase_count_or_wrong_candidate`.

This adapter still does **not** publish an authoritative phase set:

```text
global_stability_proven = false
accepted_phase_set_published = false
morphology_resolved = false
```

## Validation

The physical Sample-6 independent Decimal(80) oracle is regenerated before the
C++ tests. Boundary regression verifies that the Sample-6 W+H edge closes after
a fresh C1/C2a1 re-solve, whereas the Sample-6 H0+H1 geometric edge retains an
AQ-W appearance witness and therefore refuses W removal. Tagged synthetic tests
verify that large numerical disappearance thresholds cannot delete physically
required phases, and a provenance regression rejects a mismatched C2a1 chain
before re-solving.

Numerical/test head `1e77cb988b8ca1d123a36abd1a8544bc451db99e`
was validated by GitHub Actions run `34567768458` on GCC Debug + ASan/UBSan,
Clang Release and MSVC Release. The same run passed the affected top-level
Profile-C PT, C2b.2, no-W, C2a1 and C1 suites. Commits after that validated head
only synchronize this document; no production code, test, workflow, tolerance,
parameter or reference anchor changed after validation.

Frontend work remains frozen until maximum-three-phase validation and
authoritative publication are complete.
