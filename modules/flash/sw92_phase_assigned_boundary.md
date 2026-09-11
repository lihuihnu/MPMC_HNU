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

## Rules

### H disappearance: `W+H0+H1 -> W+H`

The surviving W(AQ) and H(NA) compositions seed a fresh C1 joint solve using the
same p/T/feed/model snapshot and fixed NaCl molality. The disappearing H
composition is supplied as an explicit C2a1 diagnostic start.

The lower topology closes only when:

1. C1 equations/material balance/role guard converge again; and
2. C2a1 finds no robust additional-NA witness.

If the disappearing H still produces a robust negative witness, the phase is
**not** dropped even if the original C2b.1 fraction was below the numerical
threshold.

### W disappearance: `W+H0+H1 -> H0+H1`

Both retained H(NA) compositions are supplied to a fresh no-W solve as initial
and final stability starts. The existing no-W solver must independently recover
the same-family split and targeted AQ-W appearance review.

The lower topology closes only when the no-W solver itself reports locally
closed one/two-H topology. A renewed W witness, higher-H evidence, root/property
failure or unresolved disappearance prevents phase deletion.

### Simultaneous endpoint

When W and an H disappear together, the surviving H composition seeds the same
fresh no-W path. The endpoint is accepted only if that no-W path closes; there
is no direct `three phases -> one phase` clipping rule.

## Boundary-aware wrapper

`solve_sw92_phase_assigned_pt_boundary_aware(...)` preserves all interior
behavior of the existing top-level driver. If the old driver stops after C2a1
because every C2b.1 attempt reached a disappearance state, the wrapper:

1. re-runs the bounded C2b.1 witness attempts;
2. retains only equation-converged disappearance states that satisfy the nested
   C1/C2b.1 Gibbs guard;
3. runs C2b.2 to obtain an explicit route;
4. re-solves the requested neighboring topology;
5. returns a lower locally-closed phase candidate only after the neighbor review
   closes.

This adapter still does **not** publish an authoritative phase set:

```text
global_stability_proven = false
accepted_phase_set_published = false
morphology_resolved = false
```

It supplies the missing boundary evidence required by the later authoritative
publication layer.

## Validation strategy

Two types of regression are required.

1. **True coexistence-edge re-solves** using the physical Mortezazadeh-Rasaei
   Sample-6 phase vertices from the independent Decimal(80) regression. An edge
   feed constructed from W+H or H0+H1 lies on the same three-phase common tangent
   with the third phase at zero amount; the neighboring solver must recover the
   lower topology without changing formulas or tolerances.
2. **Artificial large disappearance thresholds** using the tagged synthetic
   C2b.2 fixture. These are negative tests: if a physically non-negligible phase
   is merely declared small by a large threshold, the neighbor stability review
   must find that it is still required and refuse to drop it.

Frontend work remains frozen until boundary re-solve, maximum-three-phase
coverage and authoritative publication are complete.
