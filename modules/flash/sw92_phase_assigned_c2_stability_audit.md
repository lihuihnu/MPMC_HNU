# SW92 Profile-C Gate C2 constrained-stability audit

## Status and controlling decision

This is the post-C1, write-before-code audit for Gate C2 under

```text
SW92-equilibrium/phase-assigned-aq-na-joint/v1
```

against repository baseline

```text
main@011d1a728680f6e354a81230e9d0da08243e1e7a
```

C1 is implemented as the fixed two-phase primitive

```text
SW92-equilibrium/phase-assigned-aq-na-joint/fixed-two-phase-logK-SSI-RR/v1
```

with one physically aqueous phase fixed to `SwPhaseFamily::aqueous`, one physically non-aqueous
phase fixed to `SwPhaseFamily::nonaqueous`, one material balance, cross-phase reduced chemical-
potential equality, same-family minimum-Gibbs roots, and a relative water-rich/water-poor topology
gate.

**Audit verdict: C2 authoritative publication is BLOCKED under the previously proposed unrestricted
AQ/NA TPD review.** Independent Decimal(80) evidence from the already accepted C1 binary references
shows that an unrestricted opposite-family TPD search necessarily reintroduces the Xu/Profile-B
lower-envelope model choice at an already assigned physical phase. That contradicts the Profile-C
algorithm identity established by the profile-separation audit.

No C2 production finalizer or `accepted_phase_set()` is added by this audit. C1 remains valid as a
fixed-role equation candidate primitive with `final_stability_checked=false` and
`global_stability_proven=false`.

## 1. Why the unrestricted common-tangent review fails

Let a converged C1 pair have aqueous composition `x`, non-aqueous composition `y`, and common
reduced chemical potentials

```text
d_i = ln(x_i) + ln(phi_i^AQ(x))
    = ln(y_i) + ln(phi_i^NA(y)).
```

The C1 role assignment is part of the thermodynamic problem:

```text
physical aqueous    -> AQ model
physical nonaqueous -> NA model.
```

If C2 now evaluates an **unrestricted AQ trial** at the already assigned non-aqueous composition
`y`, its tangent-plane distance is

```text
D_AQ(y)
 = sum_i y_i [ln(y_i) + ln(phi_i^AQ(y)) - d_i]
 = sum_i y_i [ln(phi_i^AQ(y)) - ln(phi_i^NA(y))]
 = g_AQ(y) - g_NA(y).
```

Thus the sign of this TPD is exactly the opposite-family Gibbs ordering at the physical NA phase.
A negative value says that the AQ model surface is lower at the same composition. Treating that as
an admissible destabilizing phase is mathematically the same discrete model choice used by the
Xu lower-envelope profile.

The analogous identity at the physical aqueous phase is

```text
D_NA(x) = g_NA(x) - g_AQ(x).
```

Therefore running the generic family searches over the whole composition simplex is **not merely a
neutral stability check** for Profile C. It changes the admissible model space and collapses the
final review back toward Profile-B family competition.

## 2. Independent binary evidence

The standalone audit reference

```text
tests/flash/sw92_phase_assigned_joint/reference_c2_domain_decimal.py
```

imports only the already independent stdlib Decimal C1 equation implementation. It does not import
C++ production code or a Profile-A/Profile-B result.

### CO2/H2O

State:

```text
p = 3 MPa
T = 340 K
m_NaCl = 0
z_CO2 = 0.7
```

The accepted C1 equation reference is approximately

```text
x_CO2^AQ = 0.00594506516618610144...
y_CO2^NA = 0.98878773430905531066...
```

At the physical NA composition,

```text
D_AQ(y) = g_AQ(y) - g_NA(y)
        = -0.0013073695742952295...
```

which is orders of magnitude below the existing finite-search negative-TPD threshold scale. An
unrestricted AQ final search can therefore reject the C1 state immediately at the already existing
NA composition, without discovering a new physical phase topology.

At the physical AQ composition,

```text
D_NA(x) = g_NA(x) - g_AQ(x)
        = +0.01756589035909635...
```

so the asymmetry is not a roundoff effect.

### CH4/H2O

State:

```text
p = 10 MPa
T = 350 K
m_NaCl = 1 mol/kg H2O
z_CH4 = 0.5
```

The accepted C1 equation reference is approximately

```text
x_CH4^AQ = 0.00091079954727104846...
y_CH4^NA = 0.99467854446703175471...
```

The independent cross-role values are

```text
D_AQ(y) = g_AQ(y) - g_NA(y)
        = -0.0031092923866916092...

D_NA(x) = g_NA(x) - g_AQ(x)
        = +0.0045262164216192922...
```

Again, an unrestricted AQ search would reject the phase-assigned NA state for the same reason that
Profile B marks that family assignment as lower-envelope dominated.

These values are model-internal numerical evidence, not experimental validation.

## 3. Consequence for Profile-C algorithm identity

Profile C was introduced precisely because original SW92 phase-specific empirical semantics and Xu
lower-envelope semantics answer different questions.

A C2 rule of the form

```text
search AQ over the full simplex
search NA over the full simplex
reject on any negative TPD
```

makes an AQ model at `y` an admissible competitor to the physical NA phase and an NA model at `x`
an admissible competitor to the physical AQ phase. That restores the discrete family reassignment
that Profile C explicitly forbids.

Therefore the following implementation is **not authorized**:

- build the C1 common tangent;
- add both C1 phase compositions to both unrestricted family start sets;
- call the existing Xu/Gate-3A two-family common-tangent helper unchanged;
- publish only if both unrestricted searches report no instability.

That procedure is suitable for the Xu-style model space, not for the phase-assigned Profile-C model
space.

## 4. What a scientifically valid C2 would require

A phase-assigned stability calculation needs explicit admissible **physical-role domains**, not only
family-specific property evaluators. Conceptually it needs

```text
Omega_AQ(candidate state)
Omega_NA(candidate state)
```

and must minimize the AQ and NA TPD only over the corresponding domains.

The current C1 topology rule provides one necessary pairwise fact:

```text
x_water(aqueous-role phase) > x_water(nonaqueous-role phase).
```

It does **not** yet define a complete optimization domain for an arbitrary trial phase. A tempting
extension would require an AQ trial to remain water-richer than the retained NA phase and an NA trial
to remain water-poorer than the retained AQ phase. That idea is not yet sufficient for production
because:

1. the resulting TPD problem is inequality-constrained and can have a minimum on the role boundary;
2. the current generic stability solver assumes an unconstrained composition-simplex stationary
   condition and has no KKT/boundary acceptance contract for an additional role constraint;
3. rejecting an out-of-domain property call is not equivalent to solving the constrained minimum
   and may turn a valid boundary minimum into a numerical failure;
4. the meaning of an additional AQ or NA trial phase already touches phase-role multiplicity
   (`AQ+AQ+NA` or `AQ+NA+NA`), which the earlier audit assigned to a later topology gate;
5. no current source evidence supplies an absolute water-fraction cutoff, and this project must not
   invent one merely to make the finalizer pass.

## 5. Revised gate status

The post-C1 evidence revises the earlier conditional C2 decision:

| Gate | Revised status | Reason |
| --- | --- | --- |
| C1 fixed AQ+NA joint equations | **PASS / implemented** | Direct Decimal(80), material balance, chemical-potential equality and relative role topology are available. |
| C2 unrestricted AQ/NA common-tangent TPD | **REJECTED** | Reintroduces Profile-B family dominance at an assigned physical phase. |
| C2 role-constrained stability | **BLOCKED pending C2a audit** | Physical-role domains and constrained/KKT numerical semantics are not yet defined. |
| Authoritative Profile-C `accepted_phase_set()` | **BLOCKED** | No scientifically valid final stability contract yet. |
| C3 autonomous single-phase/phase-number | **BLOCKED** | Single-phase physical role remains unresolved. |
| C4 physical three-component Profile-C regression | **BLOCKED until C2** | C1 equations alone are not authoritative equilibrium publication. |
| C5 three-phase / physics coupling | **BLOCKED** | Depends on role multiplicity and authoritative equilibrium semantics. |

## 6. Required C2a design audit

The next increment should be a design/scientific audit, not a solver patch. It must decide:

1. whether Profile-C stability is defined only within the fixed `AQ+NA` two-phase topology, or also
   against creation of additional same-role phases;
2. if trial roles are constrained by relative water ordering, the exact feasible sets for AQ and NA
   trials and their roundoff/tie boundary;
3. the KKT/stationarity condition on the role boundary;
4. how a finite search can distinguish `no_negative_found` from failure to resolve a constrained
   boundary minimum;
5. how phase disappearance and a trial crossing the physical-role boundary are reported;
6. whether a binary-only constrained stability primitive should be validated first before any
   runtime-multicomponent generalization;
7. which independent Decimal(80) constrained minima and negative-witness cases are required before
   authoritative publication.

Until this audit is closed, C1 candidates must remain non-authoritative.

## Verification policy

This increment changes no production C++ or thermodynamic parameter. The new standalone Decimal(80)
script regenerates the two cross-role TPD identities above from the independent C1 equations. The
existing Profile-C hosted workflow is extended only to run that audit reference in addition to the
unchanged C1 tests/reference. No tolerance or existing anchor is modified.