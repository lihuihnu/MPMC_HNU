# SW92 Profile-C top-level PT topology driver

## Scope

`sw92_phase_assigned_pt.hpp` is the first single backend entry point that
orchestrates the already-audited Profile-C topology gates for a given
`p, T, z` state:

```text
no-W fixed-NA candidate generation + targeted W appearance
        |
        +-- locally closed 1 NA candidate
        +-- locally closed 2 NA candidates
        |
        +-- W witness
                |
                v
          C1 W(AQ)+H(NA)
                |
                v
          C2a1 H-side NA addition
             /        \
       no witness    witness
          |             |
          v             v
       W + H       C2b.1 W+H0+H1
                         |
                         v
                    C2b.2 closure
```

Algorithm identity:

```text
SW92-equilibrium/phase-assigned-aq-na-joint/
topology-orchestration-pt/v1
```

The driver adds **no new thermodynamic equation**. It is an orchestration layer
around the validated no-W, C1, C2a1, C2b.1 and C2b.2 primitives.

## Result contract for frontend/service use

The result exposes a stable topology status vocabulary:

```text
no_w_single_h_locally_closed
no_w_two_h_locally_closed
w_h_locally_closed
w_h0_h1_locally_closed
higher_phase_count_or_wrong_candidate
topology_unresolved
numerical_indeterminate
```

For a locally closed candidate, `phases` contains one, two or three explicit
phase instances. Each instance carries:

- physical role: `aqueous` or `nonaqueous_unclassified`;
- SW92 thermodynamic family: AQ or NA;
- molar phase fraction;
- ordered component mole fractions;
- optional compressibility factor when already available from the retained
  primitive.

This is deliberately enough for a frontend to display phase count candidates,
compositions, fractions and diagnostics without reimplementing flash logic.

It is **not** an authoritative physical phase-set publication:

```text
global_stability_proven = false
accepted_phase_set_published = false
morphology_resolved = false
```

No nonaqueous phase is labeled liquid or vapor.

## no-W routing

The driver starts from the morphology-neutral no-W adapter. If its fixed-NA
multiplicity and targeted W-appearance search are locally closed, the driver
returns the corresponding one- or two-NA candidate state directly.

If the no-W adapter finds a targeted AQ W seed, the driver continues into joint
Profile-C equations rather than accepting that seed as a phase.

## C1 candidate attempts

A W witness may coexist with multiple retained NA numerical candidates. The
driver therefore tries each NA candidate that is water-poorer than the W seed
and has representable common support.

For each attempt:

```text
logK_i = log(x_i^H / x_i^W)
```

initializes the existing C1 `W(AQ)+H(NA)` joint primitive.

A converged C1 candidate is retained only if its recomputed reduced Gibbs does
not exceed the retained no-W common-tangent Gibbs beyond the summed arithmetic
roundoff guards. This is a nested physical-topology comparison, not global
AQ/NA model competition. Among admissible attempts the lowest reduced-Gibbs C1
candidate is retained.

Attempt quotas remain explicit. Exhaustion never becomes success.

## C2a1 and C2b.1 routing

The retained C1 H composition is reviewed by C2a1. Distinct retained no-W NA
candidates are supplied as extra H-side starts, which improves candidate
recovery without changing C2a1's role guards.

If C2a1 finds no additional H witness, the top-level state is
`w_h_locally_closed`.

If C2a1 finds usable additional-NA witnesses, the driver tries the existing
C2b.1 source adapter for each witness within a bounded attempt quota. A
three-phase candidate must also have reduced Gibbs no higher than the retained
C1 candidate beyond arithmetic guards. The lowest admissible candidate is sent
to C2b.2.

Only a C2b.2 result of
`w_present_h_multiplicity_locally_closed` becomes the top-level
`w_h0_h1_locally_closed` status.

Higher-H evidence, disappearance routing, source-chain inconsistency and
numerical/property failures remain explicit unresolved/failure states.

## Why the driver does not implement PIP

The current platform goal is to finish thermodynamic/topology flash correctness
and then integrate the frontend. Hydrocarbon morphology is orthogonal to these
joint equations because all H instances use the same SW92 NA family.

Consequently v1 deliberately exposes

```text
nonaqueous_unclassified
```

rather than blocking the backend/frontend boundary on an unvalidated
`LV / LL` classifier.

A future morphology module can annotate these phase instances after PIP/source/
derivative validation without changing the top-level PT material-balance or
chemical-potential contract.

## Validation

The focused regression must cover at least:

- traceable wet CO2/H2O -> locally closed `W+H`;
- dry CO2/H2O -> locally closed no-W single NA/H candidate;
- zero-water inventory -> no-W single NA/H candidate;
- synthetic CH4/CO2/H2O C2b.1 structural fixture -> locally closed
  `W+H0+H1` through the **top-level** entry point;
- runtime component permutation;
- bounded candidate-attempt failure semantics;
- public-header self containment.

The dedicated hosted workflow reruns all directly affected no-W/C1/C2a1/C2b.1/
C2b.2 regressions and their independent Decimal references. No formula,
parameter, tolerance or reference anchor is changed merely to make orchestration
pass.
