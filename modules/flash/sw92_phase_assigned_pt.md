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
        |       |
        |       v
        | C1 W(AQ)+H(NA)
        |       |
        |       v
        | C2a1 H-side NA addition
        |    /        \
        | no witness  witness
        |    |           |
        |    v           v
        | W + H     C2b.1 W+H0+H1
        |                |
        |                v
        |           C2b.2 closure
        |
        +-- fixed-NA pair has further NA instability
                |
                v
          targeted W-rival review
                |
          if W seed exists, continue into C1/C2a/C2b
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

This is deliberately enough for a frontend/service boundary to display phase
count candidates, compositions, fractions and diagnostics without duplicating
EOS or flash logic.

It is **not** an authoritative physical phase-set publication:

```text
global_stability_proven = false
accepted_phase_set_published = false
morphology_resolved = false
```

No nonaqueous phase is labeled liquid or vapor.

## no-W routing

The driver starts from the morphology-neutral no-W adapter. If fixed-NA
multiplicity and targeted W-appearance are locally closed, the driver returns
the corresponding one- or two-NA candidate state directly.

If the no-W adapter finds a targeted AQ W seed, the driver continues into joint
Profile-C equations rather than accepting that seed as a phase.

### Unstable fixed-NA pair is not a terminal topology

A fixed-NA two-phase result may fail its same-family final stability review.
That evidence means only that the all-NA maximum-two-phase candidate is not
closed. It does **not** prove the physical Profile-C state has more than three H
phases: the correct solution may instead be

```text
W(AQ) + H0(NA) + H1(NA).
```

Therefore the top-level driver reconstructs the selected fixed-NA pair and runs
the same targeted AQ/W appearance review before terminating at
`higher_phase_count_or_wrong_candidate`.

If a usable targeted W seed exists, the graph proceeds through C1/C2a1/C2b.1/
C2b.2. This behavior is required by the synthetic structural three-phase
regression and does not modify the underlying fixed-family VLE semantics.

## C1 candidate attempts

A W witness may coexist with multiple retained NA numerical candidates. The
driver tries each candidate that is water-poorer than the W seed and has
representable common component support.

For each attempt:

```text
logK_i = log(x_i^H / x_i^W)
```

initializes the existing C1 `W(AQ)+H(NA)` joint primitive.

### No cross-topology AQ/NA Gibbs ranking

The driver deliberately does **not** reject C1 by comparing

```text
G(all phases evaluated with NA)
```

against

```text
G(W evaluated with AQ + H evaluated with NA).
```

Those are different physical-role/model assignments in the asymmetric SW92
Profile-C contract. Treating their feed-level Gibbs surfaces as globally
competing would recreate Profile-B lower-envelope semantics.

The targeted W instability creates the rival physical topology; C1 joint
equations, material balance, common reduced chemical potentials and relative
water-role guard decide whether an admissible `W+H` candidate exists. If
multiple C1 attempts converge, the lowest reduced-Gibbs candidate **within that
same W(AQ)+H(NA) topology** is retained.

Attempt quotas remain explicit. Exhaustion never becomes success.

## C2a1 and C2b.1 routing

The retained C1 H composition is reviewed by C2a1. Distinct retained no-W NA
candidates are supplied as extra H-side starts, improving candidate recovery
without changing C2a1's role guards.

If C2a1 finds no additional H witness, the top-level state is
`w_h_locally_closed`.

If C2a1 finds usable additional-NA witnesses, the driver tries the existing
C2b.1 source adapter for each witness within a bounded attempt quota.

Here a reduced-Gibbs comparison **is** valid as a nested topology guard because

```text
W(AQ)+H(NA)
    ->
W(AQ)+H0(NA)+H1(NA)
```

keeps the physical model assignment fixed while splitting one NA phase into two
NA phase instances. C2b.1 must not increase reduced Gibbs beyond arithmetic
roundoff guards. The lowest admissible three-phase candidate is sent to C2b.2.

Only C2b.2
`w_present_h_multiplicity_locally_closed` becomes the top-level
`w_h0_h1_locally_closed` status.

Higher-H evidence, disappearance routing, source-chain inconsistency and
numerical/property failures remain explicit unresolved/failure states.

## Why the driver does not implement PIP

The current platform priority is to finish thermodynamic/topology flash
correctness and integrate the frontend. Hydrocarbon morphology is orthogonal to
these joint equations because all H instances use the same SW92 NA family.

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

The focused regression covers:

- traceable wet CO2/H2O -> locally closed `W+H` and independent C1 Decimal
  golden values;
- dry CO2/H2O -> locally closed no-W single NA/H candidate;
- zero-water inventory -> no-W single NA/H candidate;
- synthetic CH4/CO2/H2O structural fixture -> locally closed
  `W+H0+H1` through the **top-level single entry point**, including recovery from
  an unstable fixed-NA pair into the W-rival topology;
- runtime component permutation;
- bounded candidate-attempt failure semantics;
- public-header self containment.

Numerical code head `372b84b59ac7cdd1710172e2108290021b2a382f`
was validated by GitHub Actions run `34561821136`:

- GCC Debug + ASan/UBSan: success;
- Clang Release: success;
- MSVC Release: success;
- top-level Profile-C PT focused regression: success on all three platforms;
- affected no-W, C2b.2, C2b.1, C2a1 and C1 regression suites: success on all
  three platforms;
- all independent Decimal references in that dependency chain: success.

No formula, parameter, tolerance or reference anchor was changed to make the
orchestration pass.
