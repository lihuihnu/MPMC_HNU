# SW92 Profile-C post-C2b.1 rival-topology and hydrocarbon-morphology audit

## Status, scope, and controlling decision

This is a docs-only scientific/architecture audit for the state reached after merging
Profile-C Gate C2b.1:

```text
W(AQ) + H0(NA) + H1(NA)
```

where `H0/H1` are unordered hydrocarbon phase instances and neither is yet called liquid or
vapor.

Baseline:

```text
main@46cce1e720e3324cc0e51b6bef285b73e5b7c8db
```

The baseline already contains:

- C1 `W(AQ)+H(NA)` joint candidate equations;
- C2a0 physical-topology separation;
- C2a1 H-side additional-NA witness search;
- C2a2 separation of phase instance, AQ/NA family, cubic root/Z and physical L/V morphology;
- C2b.1 unordered `W(AQ)+H0(NA)+H1(NA)` joint three-phase candidate equations.

This audit changes **no** production C++, SW92 equation, parameter, tolerance, root selector, generic
TPD implementation, numerical reference, CI workflow, derivative or phase-set publication API.

### Controlling verdict

The post-C2b.1 problem must be split into two distinct gates:

```text
C2b.2  rival-topology closure while morphology is still unresolved
        |
        v
C2c    post-solve H-pair morphology evidence: LV / LL / unresolved
```

The gates must not be collapsed into one finalizer.

**GO** for a next production increment that closes the *W-present nested topology* around the
already-converged C2b.1 candidate without assigning L/V.

**CONDITIONAL GO** for a future derivative-based PIP morphology diagnostic/result contract.

**NO-GO** for authoritative `LV`/`LL` publication from Z, root index, density ordering,
pseudo-critical temperature, or an unvalidated PIP threshold.

**NO-GO** for Profile-C `accepted_phase_set()` immediately after C2b.1. Rival no-water and
single-phase topologies are still unresolved, and morphology evidence is not authoritative.

SW92 physics coupling therefore remains blocked for role-dependent constitutive laws.

---

## 1. Evidence update after C2a2

The C2a2 source hierarchy remains controlling. This audit additionally rechecked the two primary
source gaps that matter most to post-C2b.1 orchestration.

### 1.1 Venkatarathnam-Oellrich 2011

Primary publication metadata and publisher abstract are available for:

G. Venkatarathnam and L. R. Oellrich,
*Identification of the phase of a fluid using partial derivatives of pressure, volume, and
temperature without reference to saturation properties: Applications in phase equilibria
calculations*, Fluid Phase Equilibria 301 (2011) 225-233,
DOI `10.1016/j.fluid.2010.12.001`.

The abstract explicitly states that the derivative-based method identifies liquid/vapor without
saturation properties and is intended for phase-equilibrium applications including liquid-liquid
and vapor-liquid-liquid calculations.

The complete 2011 primary article is still not present in the project evidence bundle. Exact
production reproduction of its mixture handling, critical/supercritical branch logic and all
special cases therefore remains blocked.

### 1.2 Jayanti-Venkatarathnam 2016 strengthens, but does not remove, the PIP caution

P. C. Jayanti and G. Venkatarathnam,
*Identification of the phase of a substance from the derivatives of pressure, volume and
temperature, without prior knowledge of saturation properties: Extension to solid phase*,
Fluid Phase Equilibria 425 (2016) 269-277,
DOI `10.1016/j.fluid.2016.06.001`, reproduces the fluid PIP criterion from the earlier work.

For molar volume `v`, the parameter is given as

```text
        [ (d2p/dv dT)       (d2p/dv2)_T ]
Pi = v [ ------------  -  --------------- ]
        [ (dp/dT)_v          (dp/dv)_T    ]
```

and the paper states the fluid-side classification in the form

```text
Pi > 1   -> liquid / liquid-like vapor
Pi <= 1  -> vapor
```

This is materially stronger evidence than Z/root ordering because it is built from local
thermodynamic derivatives.

However, the same paper explicitly notes a critical limitation: at high temperatures a vapor may
be identified as liquid-like vapor, and additional procedure is required to distinguish those
states correctly. The threshold is therefore not an authorization to force every state into a
binary physical liquid/vapor enum.

For MPMC_HNU this reinforces, rather than removes, the need for

```text
near_critical_or_coalescing
insufficient_evidence
out_of_applicability
numerical_indeterminate
```

states.

### 1.3 Li-Nghiem 1986 remains an orchestration provenance gap

Y.-K. Li and L. X. Nghiem,
*Phase equilibria of oil, gas and water/brine mixtures from a cubic equation of state and Henry's
law*, Canadian Journal of Chemical Engineering 64 (1986) 486-496,
DOI `10.1002/cjce.5450640319`, is available through publisher metadata/abstract.

The abstract establishes cubic-EOS oil/gas plus Henry-law aqueous modeling and salinity treatment,
but the full primary algorithm remains unavailable to the project. Exact historical Li-Nghiem
phase-type branch logic must not be copied from later descriptions and called primary reproduction.

### 1.4 Existing P1 sources remain controlling for topology

The directly inspected Michelsen stability, Whitson-Michelsen negative flash,
Li-Firoozabadi staged stability/split, Sabet-Gahrooei, Mortezazadeh-Rasaei, Nazari et al. and
Soreide-Whitson sources continue to support the following architecture:

```text
phase-addition evidence
    -> solve candidate topology
    -> revalidate whole candidate
    -> resolve rival topology / disappearance
    -> only then resolve physical morphology metadata
```

---

## 2. C2b.1 closed the equations, not the topology graph

A converged C2b.1 state proves, within its numerical tolerances, that one particular imposed
three-phase candidate satisfies

```text
one common EOS-component inventory
beta_W + beta_H0 + beta_H1 = 1
mu_i^W(AQ) = mu_i^H0(NA) = mu_i^H1(NA)
```

plus the existing pair-distinction and relative-water topology guards.

It does **not** prove that the physical system must contain exactly these three phases.

In particular, C2b.1 deliberately does not decide whether the correct neighboring topology is

```text
W + H
H0 + H1
H
W
```

and does not decide whether the two retained H instances are

```text
L + V
L1 + L2
unresolved
```

Therefore

```text
C2b.1 candidate_admissible()
```

must continue to mean an equation/topology candidate, not an accepted physical phase set.

---

## 3. Rival topology must be represented as a graph, not as whole-simplex AQ/NA competition

Profile C cannot return to the rejected rule

```text
at arbitrary composition x:
choose min(g_AQ(x), g_NA(x))
```

because that is the Xu-style Profile-B lower-envelope model.

The Profile-C topology graph should instead be explicit.

For the project's maximum-three-phase scope, the relevant physical nodes are conceptually:

```text
W                         one aqueous phase
H                         one hydrocarbon phase
W + H                     C1
H0 + H1                   two hydrocarbon phases
W + H0 + H1               C2b.1
H0 + H1 + H2              no-water three-hydrocarbon topology needed for eventual V+L+L scope
```

The last node is not part of C2b.2 implementation scope, but it matters to the statement that a
general arbitrary-mixture Profile-C `accepted_phase_set()` is still blocked.

Phase-addition/disappearance edges, rather than global family competition, are the transferable
control structure:

```text
W        <-> W+H        by role-targeted H appearance/disappearance
H        <-> W+H        by role-targeted W appearance/disappearance
H        <-> H0+H1      by same-NA H splitting
W+H      <-> W+H0+H1    by same-NA H splitting (C2a1 -> C2b.1)
H0+H1    <-> W+H0+H1    by W appearance/disappearance
H0+H1    <-> H0+H1+H2   by another NA phase addition
```

This keeps the physical topology as the orchestration variable and AQ/NA as the model assignment
inside that topology.

---

## 4. What C2b.2 can close now without L/V information

The next numerical gate may be implemented now because it needs no morphology classifier.

Recommended identity:

```text
C2b.2 — W-present nested-topology closure for C2b.1
```

Its scope is intentionally narrower than an authoritative phase-number solver.

### 4.1 Revalidate all retained evidence

C2b.2 should consume an admissible chain

```text
C1 -> C2a1 -> C2b.1
```

from one ordered model snapshot and revalidate the acceptance-critical information instead of
trusting cached status enums alone:

- common ordered component IDs, dataset and revision;
- p, T and fixed NaCl molality;
- C1 material balance and W/H common chemical potentials;
- C2a1 selected robust negative NA witness and role admissibility;
- C2b.1 three-phase material balance and both W-H chemical-potential residual sets;
- all three positive fractions above the declared threshold;
- pairwise composition distinction;
- W water-richer than both H phases;
- same-family minimum-Gibbs root smoothness.

### 4.2 Require the added phase to resolve the source instability

At the converged C2b.1 common tangent, the old C2a1 witness direction should no longer remain a
robust negative direction merely because the added phase was ignored.

At minimum, the final H0/H1 compositions must both be explicit starts of the post-C2b.1 NA search,
and the selected C2a1 witness must be retained as a diagnostic start.

### 4.3 Compare the nested split against its source C1 state

Within the nested

```text
W(AQ)+H(NA)
    ->
W(AQ)+H0(NA)+H1(NA)
```

transition, model assignment is unchanged: one W remains AQ and the hydrocarbon material remains in
NA. A solved phase addition that arose from a robust negative NA TPD should not produce a clearly
higher reduced Gibbs candidate than the revalidated C1 source.

Therefore C2b.2 may use a guarded consistency check

```text
G_C2b1 <= G_C1 + arithmetic_guard
```

as a rejection/indeterminate gate.

This is **not** permission to rank arbitrary `W-only` AQ and `H-only` NA states by lower family Gibbs
at feed composition. That would recreate Profile B.

### 4.4 Final H-side phase-addition search must remain role-constrained

A post-C2b.1 NA-only finite TPD search is scientifically meaningful because every hydrocarbon phase
uses the same NA family.

However, just as in C2a1, a mathematically negative NA trial near the retained W composition is not
an admissible new hydrocarbon phase.

A usable additional-H witness must therefore be:

- robustly negative under the declared tolerance/roundoff semantics;
- distinct from **both** H0 and H1;
- water-poorer than retained W by a guarded margin;
- on the same active feed support;
- produced without changing the SW92 family/root contract.

If such a witness exists while W, H0 and H1 are already present, the candidate requires an
additional phase or indicates a wrong local candidate. Under the current maximum-three-phase scope,
C2b.2 must report an explicit higher-phase-count/out-of-scope or unresolved status. It must not
silently discard the witness and accept W+H0+H1.

### 4.5 No-negative-witness is still finite-search evidence

If all required NA searches complete without a usable additional-H witness, the correct statement is
only

```text
W-present H multiplicity locally closed under the declared finite search
```

not

```text
global phase number proven
```

`global_stability_proven=false` remains mandatory.

---

## 5. Phase disappearance is routing evidence, not final publication

C2b.1 already distinguishes aqueous and hydrocarbon phase disappearance. C2b.2 should convert
those statuses into topology-routing evidence, not accepted results.

### 5.1 One H disappears

If one H fraction reaches the disappearance threshold after equation convergence, the neighboring
candidate is

```text
W + H
```

and should be routed back to the C1 topology path for revalidation.

It must not be interpreted as proof that the third phase physically cannot exist beyond the local
boundary.

### 5.2 W disappears

If W reaches its disappearance boundary, the neighboring candidate is

```text
H0 + H1
```

not a failed three-phase solve and not an accepted no-water state.

A dedicated no-W NA/NA candidate/review path is still required.

### 5.3 Two phases disappear

A one-phase endpoint requires an autonomous role contract:

```text
W-only AQ
or
H-only NA
```

The current Profile-C architecture cannot choose those roles by simply comparing `g_AQ(z)` and
`g_NA(z)` at the feed. Single-phase autonomous role remains a separate blocker.

---

## 6. Why an accepted Profile-C phase set is still blocked after C2b.2

Even a successful W-present nested-topology closure cannot publish the general Profile-C phase set.
At least three scientific gaps remain:

1. **no-W topology orchestration** — `H`, `H0+H1`, and eventually `H0+H1+H2` for the project's
   `V+L+L` scope;
2. **single-phase physical-role selection** — W-only versus H-only without whole-simplex family
   competition;
3. **hydrocarbon morphology** — L/V labels required by role-dependent physics.

Accordingly, the strongest publication after C2b.2 should remain diagnostic/candidate evidence such
as

```text
w_present_h_multiplicity_locally_closed
```

not `accepted_phase_set()`.

---

## 7. Morphology is a post-topology metadata problem

Once a topology containing two distinct H phases has been selected for morphology review, the
thermodynamic equations no longer depend on whether the pair is called

```text
L + V
or
L1 + L2
```

because both phases already use SW92 NA.

This is the key reason morphology must remain post-solve.

A pair-level result should distinguish at least:

```text
resolved_liquid_vapor_pair
resolved_liquid_liquid_pair
near_critical_or_coalescing
insufficient_evidence
evidence_conflict
out_of_applicability
numerical_indeterminate
```

The last five are first-class physical/software states, not exceptional fallbacks.

---

## 8. Evidence hierarchy for a future pair-morphology resolver

The resolver should not average or vote across unrelated heuristics. Evidence has different
scientific weight.

Recommended hierarchy:

### Tier M0 — hard validity/coalescence gates

Before asking L/V:

- both H phases must be admissible, distinct and have positive fractions;
- root/property evaluation must be smooth and within model applicability;
- composition and molar-volume separation must not be below declared coalescence guards;
- derivative denominators required by a morphology metric must be representable and well
  conditioned.

Failure here returns `near_critical_or_coalescing`, `out_of_applicability` or
`numerical_indeterminate` before any L/V label is attempted.

### Tier M1 — validated thermodynamic phase-identification metric

A validated Venkatarathnam-Oellrich-style PIP or an equivalently audited derivative-based quantity
is the best current candidate for primary morphology evidence.

It is not yet production-authorized for arbitrary SW92 mixtures.

### Tier M2 — independent auxiliary phase-type evidence

After separate audits, supporting evidence may include:

- a Sabet-style projected hydrocarbon-only L/V calculation;
- trusted continuation from a nearby already-resolved state;
- independently available phase-envelope/saturation information.

These may corroborate M1 or identify a conflict. They may not override a hard M0 invalidity gate.

### Tier M3 — compatibility diagnostics only

Mortezazadeh-2017 pseudo-critical PSD may be retained as a named compatibility diagnostic after
sourced `Vc` and validation are added.

### Representation-only evidence

Relative Z/molar volume belongs below all morphology tiers. It remains useful for deterministic
slot ordering and for orienting dense/light members *after* an `LV` morphology is independently
resolved. It never decides `LV` versus `LL`.

---

## 9. PIP contract required before implementation

A future SW92-NA PIP implementation must define derivatives at fixed composition and fixed external
NaCl molality.

For the PR-form pressure function

```text
p = p(v, T, x, family=NA, molality)
```

PIP requires, at minimum,

```text
(dp/dT)_v
(dp/dv)_T
(d2p/dv dT)
(d2p/dv2)_T
```

with the full SW92 temperature dependence included:

- pure-component alpha(T);
- water alpha(T,m);
- temperature-dependent AQ correlations where applicable to other uses;
- family BIPs and all terms entering `a_mix(T,x,m)`;
- fixed composition and fixed molality semantics made explicit.

The current SW92 flash path is floating-point primal and does not expose this derivative contract.
No finite-difference approximation should be hidden behind an API named as authoritative PIP.

### Required singularity semantics

The PIP expression contains derivative denominators. If

```text
(dp/dT)_v
or
(dp/dv)_T
```

is zero, nonrepresentable or too ill-conditioned for a trustworthy ratio, the result must be
`numerical_indeterminate` or a named near-spinodal/critical state. It must not be regularized into a
forced label by adding an arbitrary epsilon.

---

## 10. How PIP evidence may eventually map to LV / LL / unresolved

The following is a **future validation target**, not an implemented rule.

After SW92-specific mixture validation and critical-region guards exist:

### 10.1 Strong LV evidence

If one H phase is confidently vapor-like under the validated PIP contract while the other is
confidently liquid/liquid-like, and M0 gates show the phases are well separated, the pair may be
resolved as

```text
LV
```

with the individual vapor/liquid instance determined by the phase-identification evidence, not slot
number.

### 10.2 Candidate LL evidence

If both distinct H phases are confidently liquid/liquid-like, an `LL` result is plausible and is
consistent with the project's required liquid-liquid capability.

But `Pi > 1` includes liquid-like vapor in the 2016 formulation. Therefore an arbitrary-mixture
`both Pi > 1 -> LL` rule is **not** authorized until dedicated LLE and dense-supercritical mixture
validation demonstrates the intended behavior.

### 10.3 Both vapor-like or conflicting evidence

If both phases appear vapor-like, or PIP and an independent audited diagnostic disagree, the
resolver must return `evidence_conflict` / `insufficient_evidence`. It must not invent `VV`, choose
by Z ordering, or silently discard one phase.

### 10.4 Threshold/critical neighborhood

If either PIP lies within its numerical/validation guard around the decision boundary, or the two
phase states approach coalescence, return

```text
near_critical_or_coalescing
```

or `insufficient_evidence` as appropriate.

A stable property evaluation does not imply that a discrete morphology label is authoritative.

---

## 11. Continuation is supporting evidence only

A resolved state at a neighboring pressure/temperature may help track phase instances, but the
history cannot create current identity.

Continuation must be invalidated or downgraded when:

- a phase appears/disappears;
- H0/H1 coalesce;
- a critical endpoint is crossed;
- branch/root topology changes;
- two liquid branches exchange density ordering;
- current derivative evidence conflicts with the previous label.

A stale `liquid` or `vapor` tag must never override current thermodynamic evidence.

---

## 12. Sabet-style water-removed L/V evidence remains auxiliary

Sabet and Gahrooei provide a useful architecture precedent: their water-removed L/V calculation is
used to determine phase number/type, not final full-system equilibrium properties.

For MPMC_HNU, adopting such a diagnostic still requires a dedicated projected-subsystem contract:

- remove water by explicit component projection;
- renormalize the hydrocarbon inventory;
- retain sourced non-water/non-water BIPs;
- define what happens when projected support is empty or nearly singular;
- use an NA-compatible hydrocarbon model without pretending it is the original full SW92 snapshot;
- map only the **diagnostic** back to the full-component phase instance.

It must not change the converged C2b.1 compositions/fractions or replace full-system material
balance.

---

## 13. Mortezazadeh PSD remains a named compatibility diagnostic

The pseudo-critical-temperature PSD remains useful for reproducing a published engineering path,
but its status is unchanged:

- it assumes a restricted oil/gas topology;
- it requires critical-volume `Vc` data absent from the current parameter contract;
- it is not a Gibbs/stability theorem;
- it is not adequate for arbitrary `LL` / `VLL` morphology.

If implemented later, it must have a distinct identity such as

```text
hydrocarbon-role/mortezazadeh2017-pseudocritical-psd/v1
```

and cannot silently become the Profile-C default.

---

## 14. Result/data-model implications

Phase instance identity, family and morphology remain distinct fields.

Conceptually:

```text
instance_id / canonical_slot
thermodynamic_family = AQ or NA
physical_role = aqueous / hydrocarbon_liquid / hydrocarbon_vapor /
                hydrocarbon_unclassified
morphology_evidence = ...
```

For a pair result:

```text
resolved_liquid_vapor_pair
    -> one H instance may be published liquid, the other vapor

resolved_liquid_liquid_pair
    -> both H instances may be published hydrocarbon_liquid

all unresolved states
    -> H instances remain hydrocarbon_unclassified
```

Role enums are not unique keys; two liquid instances are valid.

Every resolved morphology must retain the evidence actually used: PIP values/guards, derivative
status, auxiliary diagnostic result, continuation provenance and any conflict flags.

---

## 15. What is explicitly rejected after C2b.1

This audit rejects:

- `H0 = liquid`, `H1 = vapor` by slot convention;
- low-Z H = liquid and high-Z H = vapor before independent LV resolution;
- `second NA phase = vapor`;
- `both Pi > 1 = LL` before mixture/LLE/supercritical validation;
- forcing a PIP label when its derivative ratios are singular/ill-conditioned;
- using Mortezazadeh PSD as an undocumented default;
- using a water-removed auxiliary flash as final equilibrium;
- comparing feed-level AQ and NA Gibbs to choose W-only versus H-only;
- treating C2b.1 equation convergence as phase-number proof;
- accepting W+H0+H1 when a robust distinct role-admissible additional-H witness remains;
- silently converting W disappearance to an accepted H0+H1 state;
- publishing role-dependent SW92 physics from unresolved H morphology.

---

## 16. Gate decisions

| Gate / question | Decision | Required meaning |
| --- | --- | --- |
| C2b.1 unordered W+H0+H1 equations | **PASS / implemented** | Candidate equations only; H morphology unresolved. |
| C2b.2 revalidate C1->C2a1->C2b.1 chain | **GO** | Do not trust cached statuses alone. |
| C2b.2 guarded C2b.1-vs-C1 Gibbs consistency | **GO** | Only nested same-role-map comparison; not global AQ/NA competition. |
| C2b.2 final NA-only additional-H search | **GO** | Must guard W-like trials and distinguish from both H phases. |
| Additional admissible H witness from W+H0+H1 | **NO ACCEPT** | Higher phase count / wrong candidate / out-of-scope; remain unresolved. |
| One H disappears | **ROUTE** | Return to W+H topology review; do not auto-accept. |
| W disappears | **ROUTE** | Enter dedicated H0+H1 no-W path; do not auto-accept. |
| Whole-simplex AQ/NA lower-Gibbs finalizer | **REJECTED** | This is Profile B, not Profile C. |
| Single-phase W-only vs H-only by feed family Gibbs | **REJECTED** | Autonomous physical-role contract still missing. |
| Relative Z/molar-volume LV-vs-LL classifier | **REJECTED** | Canonicalization only. |
| PIP evidence/result API | **CONDITIONAL GO** | Exact derivative contract + guards + validation required. |
| PIP as immediate authoritative arbitrary-mixture resolver | **BLOCKED** | 2011 full-text/mix handling and SW92 derivative validation incomplete. |
| Sabet hydrocarbon-only diagnostic | **CONDITIONAL GO** | Auxiliary evidence only after projected-subsystem audit. |
| Mortezazadeh PSD diagnostic | **CONDITIONAL GO** | Named compatibility path only after sourced Vc + validation. |
| W+H0+H1 -> W+L+V / W+L1+L2 publication | **BLOCKED** | Requires validated morphology evidence and rival-topology closure. |
| Profile-C accepted_phase_set() | **BLOCKED** | no-W/single-phase topology and morphology not closed. |
| SW92 role-dependent physics | **BLOCKED** | End-to-end physical roles still not authoritative. |

---

## 17. Recommended next production increment

The next implementation should be **C2b.2 only**, not a morphology classifier:

```text
C2b.2 — W-present post-C2b.1 nested-topology closure
```

Minimum contract:

1. consume matching C1, C2a1 and C2b.1 evidence from one model snapshot;
2. independently revalidate the acceptance-critical equations/balance/topology guards;
3. recompute/verify the three-phase common tangent;
4. retain H0, H1 and the selected C2a1 witness as mandatory NA-search starts;
5. run a finite NA-only search against that tangent;
6. classify a negative trial as an additional-H witness only when it is distinct from both H phases
   and remains water-poorer than W;
7. compare C2b.1 and source C1 reduced Gibbs with an arithmetic guard as a nested consistency check;
8. route H disappearance back to W+H review;
9. route W disappearance to a future no-W H0+H1 path;
10. retain explicit resource/root/nonsmooth/indeterminate states;
11. publish no L/V labels;
12. publish no authoritative phase set;
13. keep `global_stability_proven=false`.

Independent validation should include:

- the existing synthetic C2b.1 structural anchor;
- a controlled extra-H negative witness case if one can be produced without falsifying thermodynamic
  provenance;
- no-extra-H locally closed case;
- source C1/C2a1/C2b.1 mismatch guards;
- H disappearance and W disappearance routing;
- component permutation;
- root/resource failure propagation;
- unchanged C1/C2a1/C2b.1 regressions.

No PIP implementation is required for C2b.2.

---

## 18. Morphology work that should follow C2b.2

After C2b.2, the morphology track should proceed separately:

```text
C2c.0  exact PIP/source/derivative contract audit
C2c.1  SW92-NA PIP derivative kernel + independent derivative validation
C2c.2  mixture morphology regression: V, L, LV, LL, near-critical/coalescing
C2c.3  pair morphology evidence resolver with explicit unresolved/conflict states
```

Before C2c.3 may publish authoritative `LV` or `LL`:

- obtain and inspect the complete Venkatarathnam-Oellrich 2011 primary paper or otherwise close the
  exact primary-method provenance gap;
- validate all required pressure derivatives independently;
- include real, traceable mixture states covering vapor, liquid, liquid-liquid and vapor-liquid
  behavior;
- include critical/supercritical and phase-coalescence regression;
- demonstrate that `Pi > 1` liquid-like-vapor behavior does not silently convert unresolved states
  into `LL`;
- define evidence-conflict behavior with any auxiliary diagnostic.

---

## 19. Final recommendation

Do not implement `LV / LL` classification next.

Implement C2b.2 first so that the existing C2b.1 solution is checked as a whole W-present topology
candidate and any additional hydrocarbon instability, H disappearance or W disappearance is routed
explicitly.

In parallel, close the 2011 PIP primary-source gap and design the SW92 pressure-derivative contract.
Only after those derivatives and mixture morphology states are independently validated should a
post-solve resolver publish `LV` or `LL`.

This ordering preserves the project's audit-first rule, prevents physical labels from changing the
thermodynamic equations, and keeps `unresolved` as a scientifically valid outcome rather than a
failure to force a label.
