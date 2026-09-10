# SW92 Xu-style Gate 3B maximum-two-phase joint-split design audit

## Status and decision

This document is the write-before-code design audit for Gate 3B under the already reserved parent
profile

```text
SW92-equilibrium/xu-asymmetric-gibbs/v1
```

against repository baseline

```text
main@ce547350e24f6f13d4f79629777a8c975706a3b9
```

Gate 3A is already implemented as

```text
SW92-equilibrium/xu-asymmetric-gibbs/stability-foundation/v1
```

and provides a lower-family feed reference, one common reduced feed tangent, two family-specific
finite TPD searches against that same tangent, family-tagged witnesses and explicit resource/tie
semantics.

**Design verdict:** conditional **GO** for a staged Gate 3B implementation. The repository now has the
scientific reference-gauge and stability foundation required to define a joint maximum-two-phase
problem. The first implementation must use a new family-aware result and a new fixed-family-pair
split contract. It must not reuse `solve_pt_vle`/`iterate_pt_split` as the scientific API, must not
force one AQ plus one NA phase, and must not publish an accepted state until final two-family
common-tangent stability is resolved.

This audit does **not** implement the solver, modify numerical tolerances, change SW92 formulas,
introduce salt conservation, add derivatives, or enter three-phase flash.

## Sources and existing repository contracts reviewed

Scientific sources:

- G. Xu, W. D. Haynes and M. A. Stadtherr, *Reliable Phase Stability Analysis for Asymmetric
  Models*, Fluid Phase Equilibria 235 (2005) 152-165, DOI `10.1016/j.fluid.2005.06.016`.
- M. L. Michelsen, *The isothermal flash problem. Part I. Stability*, Fluid Phase Equilibria 9
  (1982) 1-19, DOI `10.1016/0378-3812(82)85001-2`.
- M. L. Michelsen, *The isothermal flash problem. Part II. Phase-split calculation*, Fluid Phase
  Equilibria 9 (1982) 21-40, DOI `10.1016/0378-3812(82)85002-4`.
- I. Soreide and C. H. Whitson, *Peng-Robinson predictions for hydrocarbons, CO2, N2, and H2S with
  pure water and NaCl brine*, Fluid Phase Equilibria 77 (1992) 217-240,
  DOI `10.1016/0378-3812(92)85105-H`, with the project-recorded authors' errata PDF.

Repository contracts reviewed:

- `modules/flash/sw92_xu_asymmetric_audit.md`
- `modules/flash/sw92_asymmetric_stability.md`
- `modules/flash/include/mpmc/flash/sw92_asymmetric_stability.hpp`
- `modules/flash/include/mpmc/flash/sw92_stability.hpp`
- `modules/flash/include/mpmc/flash/pt_stability.hpp`
- `modules/flash/include/mpmc/flash/rachford_rice.hpp`
- `modules/flash/include/mpmc/flash/pt_split.hpp`
- `modules/flash/include/mpmc/flash/pt_phase_set.hpp`
- `modules/flash/include/mpmc/flash/sw92_split.hpp`
- `modules/thermodynamics/include/mpmc/thermodynamics/sw92_phase.hpp`

The design below preserves the already audited fixed-molality common reduced-chemical-potential
gauge. It does not reopen the Gate 3A source-model decision.

## 1. What “joint” means in Gate 3B

A Gate 3B result is joint only when all actually retained phases participate in **one** component
inventory and **one** equilibrium problem:

```text
sum_alpha beta_alpha = 1
z_i = sum_alpha beta_alpha * x_i_alpha
sum_i x_i_alpha = 1
mu_i_alpha/(RT) = common_i  for every active transferable component
```

under one common

```text
p, T, ordered EOS-component snapshot, fixed NaCl molality.
```

“Joint” does **not** mean that one phase must use `SwPhaseFamily::aqueous` and another must use
`SwPhaseFamily::nonaqueous`.

The Xu-style thermodynamic model is the lower envelope of two family surfaces. A maximum-two-phase
solver therefore has three family multisets available in principle:

```text
AQ + AQ
AQ + NA
NA + NA
```

Forcing `AQ + NA` would impose an extra physical heuristic that is absent from the lower-envelope
model and could discard a lower-Gibbs same-family pair. The names `aqueous` and `nonaqueous` remain
**thermodynamic family/model identities** in this algorithm. They are not a post-hoc guarantee that
a converged phase is physically water-rich or water-poor.

This distinction is especially important because the SW92 family-specific BIPs are empirical and
each family is mathematically continued beyond the composition region in which its experimental
support is strongest. Gate 3B is therefore a model-internal research equilibrium under the explicit
Xu-style algorithm identity, not a claim that lower-envelope family selection repairs or removes
the original SW92 consistency/validity caveats.

## 2. Result contract — a new family-aware SW92 phase set is mandatory

The current generic `PtCandidatePhase` stores phase fraction, composition, `StabilityPhase` and
optional Z, but no thermodynamic family/model identity. Gate 3B cannot publish through that type
without losing the discrete state variable that defines the asymmetric model.

The first implementation should introduce an SW92-specific result, conceptually:

```text
Sw92AsymmetricCandidatePhase
  family : SwPhaseFamily
  mole_phase_fraction
  composition
  activity : StabilityPhase   // ln(phi), same-family minimum-Gibbs root branch, smooth flag
  compressibility_factor      // diagnostic/property metadata, not a phase-family selector

Sw92AsymmetricCandidatePhaseSet
  phases                       // one or two entries in Gate 3B

Sw92AsymmetricPtResult
  capability.maximum_phase_count = 2
  status
  p, T, normalized feed
  fixed NaCl molality
  dataset/revision/component order
  thermodynamics profile
  equilibrium profile
  split convention
  initial asymmetric stability
  split attempts and selected attempt
  final asymmetric common-tangent stability
  optional family-aware candidate phase set
  diagnostics/resource accounting
  global_stability_proven = false
```

The model-independent `PtPhaseSetStatus`/maximum-phase-count idea may be reused, but the generic
`PtCandidatePhaseSet` itself must not be the authoritative Gate 3B result until it can preserve
phase-model identity without SW92-specific coupling.

`accepted_phase_set()` on the SW92 result must return only the family-aware set and only after all
Gate 3B acceptance checks succeed.

### 2.1 Phase vector order is not a physical label

Vector slot 0/1, initialization order, root index and compressibility-factor ordering are all
numerical diagnostics. Callers must use the explicit `family` stored on each phase. No `liquid`,
`vapor`, `water-rich`, or `hydrocarbon-rich` label is implied by phase-slot order.

Repeated family identity is legal: an accepted maximum-two-phase candidate may contain two AQ
entries or two NA entries if that is what the lower-envelope model and final stability checks
support.

## 3. Fixed-family-pair nonlinear problem

Gate 3B should solve one family assignment at a time. For a fixed ordered pair

```text
(F0, F1),  F0,F1 in {AQ,NA},
```

let phase 0 have fraction `1-beta` and composition `x`, and phase 1 have fraction `beta` and
composition `y`.

For active feed support, define

```text
logK_i = ln(y_i / x_i).
```

`K` here is only an algebraic composition ratio between the two candidate phases. It is not a
vapor/liquid label and it is unrelated to the Whitson Profile-A `cross_model_equilibrium_ratio`.

For fixed `logK`, the ordinary two-phase material balance remains

```text
x_i = z_i / (1-beta + beta*K_i)
y_i = K_i*x_i
sum_i z_i*(K_i-1)/(1-beta+beta*K_i) = 0.
```

This algebra does not depend on whether F0 and F1 are the same thermodynamic family. Therefore the
existing `solve_rachford_rice` numerical kernel is scientifically reusable **as a private
fixed-K material-balance engine**.

The legacy field names `liquid`, `vapor` and `vapor_fraction` in `RachfordRiceResult` must not leak
into the new public result. Gate 3B should immediately map them internally to

```text
phase0 composition
phase1 composition
phase1 fraction
```

and expose only family-aware terminology.

### 3.1 Common chemical-potential residual

At one material-balanced state, evaluate the same-family minimum-Gibbs root for each fixed family:

```text
m_i^0 = ln(x_i) + ln(phi_i^F0(x))
m_i^1 = ln(y_i) + ln(phi_i^F1(y))
r_i   = m_i^0 - m_i^1.
```

The pair equations are converged only when

```text
max_i |r_i| <= audited fugacity/chemical-potential tolerance
```

on every active feed component **and** the returned material balance satisfies its absolute and
relative tolerances.

Zero-overall-feed components remain exactly outside support under the existing nonreactive closed
feed contract. They are not activated by floors and are not included in logarithms/residual norms.

### 3.2 Log-K successive-substitution direction remains valid for a fixed family pair

The legacy `iterate_pt_split` API cannot be reused directly because it requests
`liquid_candidate/vapor_candidate` roots and applies a Z-based phase distinction. Its underlying
log-K direction, however, has a model-independent thermodynamic derivation for a **fixed family
assignment**.

With component amounts in phase 1 denoted by `v_i = beta*y_i`, fixed total inventory gives

```text
d(G/RT) = sum_i (m_i^1 - m_i^0) d(v_i)
        = -sum_i r_i d(v_i).
```

The Rachford-Rice kinematic identities relating `d(v_i)` to `d(logK_i)` are independent of the EOS.
Each fixed SW92 family is internally a differentiable PR mixture away from root-envelope
nonsmoothness. Consequently the existing Gibbs-descent argument remains valid when F0 and F1 are
held fixed during an attempt, and the baseline update may retain the form

```text
logK <- logK + alpha * r / scale
```

with a Gibbs/residual line search.

This proof does **not** justify differentiating through a same-family root switch or through a
change of `SwPhaseFamily`. A nonsmooth minimum-root evaluation terminates the attempt as
indeterminate/property-nonsmooth; family reassignment is handled by separate attempts, not by
silently switching models inside one iteration.

### 3.3 New iteration options, not `PtSplitIterationOptions` by alias

The first Gate 3B implementation should define a dedicated option type. The tested PR76/VLE numeric
values may be adopted as an initial baseline where dimensions/meaning are identical, but the type
must omit assumptions that are not valid here.

In particular:

- keep explicit fugacity/chemical-potential tolerance;
- keep absolute/relative EOS-component material-balance tolerances;
- keep minimum phase fraction and logK/composition-separation controls;
- keep max log step, line-search controls, iterations and property budgets;
- reuse `RachfordRiceOptions` internally;
- keep separate AQ/NA root options;
- **do not use relative Z separation as a universal phase distinction rule**.

Two phases with different family models can have similar Z, and family identity is not density
ordering. Z remains a diagnostic/property quantity.

## 4. Initialization from Gate 3A negative witnesses

Gate 3A already returns a resolved lower-feed `reference_family` and family-tagged robust negative
TPD witnesses whenever its combined status is `unstable`.

For each witness with family `Fw` and composition `w`, a material-balanced seed should be created by
placing `w` in one incipient phase and constructing the complementary phase from the overall feed,
using the same positive-support/no-floor rules as the existing VLE seed construction.

The first family-pair seed must be

```text
(reference_family, Fw).
```

If `Fw != reference_family`, Gate 3B must also consider

```text
(Fw, Fw).
```

using the same witness composition as a seed. This second attempt is necessary because an
opposite-family instability does not prove that the final maximum-two-phase minimum contains the
original feed family. Without it, a feed whose lower single-phase surface is AQ could never reach an
`NA+NA` two-phase candidate.

When a witness has the same family as the feed reference, the natural pair is already

```text
(Fw, Fw).
```

With negative witnesses available from both Gate-3A family searches, this construction can cover all
three possible family multisets without hard-coding one AQ plus one NA phase.

### 4.1 No duplicate orientation attempts

For a fixed family pair, swapping both phase slots and replacing `logK` by `-logK` describes the
same material-balanced physical candidate. The legacy VLE path tries both witness orientations
because low-Z/high-Z requested roles are different. Gate 3B uses explicit fixed families and
same-family minimum-Gibbs roots, so it should not duplicate an attempt merely to swap slot names.

### 4.2 Resource-order policy

`Sw92AsymmetricStabilityResult::negative_witnesses` is diagnostic data and its current collection
order must not become an AQ-first scientific priority.

Before applying a finite `max_pair_attempts` budget, Gate 3B should build a deterministic plan from
all robust witnesses. At minimum:

1. partition/identify witnesses by family;
2. use common-reference TPD value as the primary evidence ranking because both family searches use
   the same dimensionless tangent;
3. ensure the best witness from each family that has negative evidence is represented before one
   family can consume all attempt slots, when the configured budget permits;
4. record unattempted planned seeds and `attempt_limit_reached` explicitly.

Exhausting initialization attempts can yield `indeterminate`; it must never be converted to a stable
single phase.

## 5. Pair-state acceptance before candidate selection

A fixed-family-pair iteration is only an **equation candidate**. Before it may compete as the
selected Gate 3B phase set, all of the following must hold.

### 5.1 Chemical-potential equality and material balance

The active-component residual and EOS-component material balance must both satisfy their audited
thresholds. A small fugacity residual cannot hide a balance failure and vice versa.

### 5.2 Phase disappearance

Only after fugacity and material-balance checks pass may a phase fraction at or below
`minimum_phase_fraction` be classified as `phase_disappearance`.

As in the existing VLE contract, disappearance is a numerical boundary state. It is not proof that
the original feed is stable and must not automatically publish the remaining phase as an accepted
single-phase equilibrium when Gate 3A found the feed unstable.

### 5.3 Distinct-phase check

Gate 3B should use composition/logK distinction, not Z ordering. If the phase compositions become
numerically indistinguishable, the two-phase parameterization is degenerate and the attempt is not
an accepted pair.

Different `SwPhaseFamily` values alone do not justify two nonzero phase fractions at identical
composition. At a common composition, a lower family can replace a higher family without changing
component inventory; a family tie is an envelope nonsmooth point rather than evidence for two
separate macroscopic phases.

### 5.4 Same-family root smoothness

Every phase property evaluation uses that phase's mechanically admissible same-family minimum-Gibbs
root. If the root envelope is nonsmooth/ill-conditioned, the attempt remains unresolved. Gate 3B
must not fall back to lowest-Z/highest-Z to force a branch.

## 6. Mandatory discrete family-dominance check at every converged phase

Solving chemical-potential equations for a fixed family pair is not enough. At each converged phase
composition `x_alpha`, Gate 3B must evaluate **both** family minimum-Gibbs surfaces at that exact
composition and compare them with the same family-Gibbs arithmetic guard principle used at the
Gate-3A feed.

For each phase:

```text
Delta_g_alpha/(RT)
  = sum_i x_i_alpha * [ln(phi_i^AQ(x_alpha)) - ln(phi_i^NA(x_alpha))].
```

Then:

- if the phase's assigned family is the resolved lower family: the discrete assignment is locally
  admissible;
- if the opposite family is resolved lower: the candidate is `family_assignment_dominated` and
  cannot compete for acceptance;
- if the two family surfaces tie within the arithmetic guard: the lower envelope is nonsmooth at an
  actual phase composition, so family identity is unresolved and the candidate is
  `family_assignment_nonsmooth` / indeterminate.

This check is required even though final asymmetric TPD searches are also performed. An opposite
family that is lower at the exact same composition would normally appear as negative TPD when that
composition is supplied as a start, but an exact/near family tie may yield zero TPD while still
leaving the discrete phase-family identity nonsmooth. Gate 3B must not silently choose by execution
order.

The first implementation should **not** switch a phase's family in place after detecting dominance.
That would introduce a nonsmooth discrete update into the continuous log-K iteration. Instead, reject
that fixed-family attempt and let the separately planned family assignment solve its own equations.

## 7. Candidate reduced Gibbs and selection across attempts

For a converged, family-admissible pair,

```text
G_pair/(RT)
  = (1-beta) * sum_i x_i [ln(x_i) + ln(phi_i^F0(x))]
  + beta     * sum_i y_i [ln(y_i) + ln(phi_i^F1(y))].
```

The omitted component-specific reference terms are linear in the total component inventory and are
identical for every candidate at the same `p,T,z,molality`, so they cancel in candidate-to-candidate
and pair-to-feed Gibbs comparisons under the Gate-3A common gauge.

Among all equation-converged, family-admissible candidates, the solver should prefer the lowest
resolved reduced Gibbs value.

### 7.1 Candidate Gibbs ties require deduplication before declaring nonsmoothness

Multiple witnesses can converge to the same phase set. Slot-swapped states can also represent the
same set. Therefore a near-equal Gibbs value must not immediately be called a distinct phase-set
tie.

The selection procedure must first identify numerically equivalent candidates, including the
phase-slot-swapped representation. Equivalence must compare:

- family identity per matched phase;
- phase fractions up to complement under a swap;
- compositions on active support;
- converged chemical-potential/material-balance residual scale.

After equivalent duplicates are removed, two **distinct** candidate phase sets whose reduced Gibbs
values are tied within their combined arithmetic guards make candidate selection nonsmooth. Gate 3B
must then return `indeterminate`; it must not pick the first attempted state.

### 7.2 Pair Gibbs must not exceed the original lower-envelope feed state

Store the Gate-3A selected-feed reduced Gibbs and compare the chosen pair against it. If the pair is
resolved higher than the feasible original feed by more than the combined arithmetic guard, it
cannot be accepted as the equilibrium state.

This condition is distinct from final TPD stability evidence and should have its own diagnostic. It
must not be hidden by a relaxed final TPD tolerance.

## 8. Final common tangent for a converged pair

For a converged pair, define

```text
m_i^0 = ln(x_i) + ln(phi_i^F0)
m_i^1 = ln(y_i) + ln(phi_i^F1)
d_i   = midpoint(m_i^0, m_i^1).
```

The pair's finite residual means this common tangent is known only to the same numerical
approximation as the chemical-potential equality. As in the existing tested VLE path, a reasonable
baseline numerical reference disagreement is

```text
common_reference_allowance = 0.5 * max_i |m_i^0 - m_i^1|.
```

This is a numerical allowance, not a rigorous error enclosure and not a physical tolerance.

Final asymmetric stability must then run **both** family searches against exactly the same `d_i`.
The final TPD tolerance may include this explicit common-reference allowance, provided the base
stability tolerance and the added allowance are both retained in result metadata.

## 9. Gate 3A should expose/reuse one imposed-common-tangent two-family helper

The existing public Gate-3A entry first constructs a feed tangent and then directly runs two calls to
`test_pt_stability_against`. Gate 3B needs the same two-family search semantics for a tangent produced
by a converged phase pair.

The recommended minimal refactor is an SW92-specific helper/result such as conceptually

```text
test_sw92_pt_asymmetric_stability_against(
    p, T, support_feed, common_log_activity,
    model, fixed_molality,
    per-family options, per-family starts)
```

which owns:

- AQ/NA evaluator construction;
- preflight validation of both family search budgets/starts;
- exact propagation of one imposed tangent to both searches;
- separate family results/evaluation counts;
- family-tagged robust negative witnesses;
- combined `unstable / no_instability_found / indeterminate` status;
- `global_stability_proven=false`.

The existing Gate-3A feed-reference entry can retain its current public result and behavior while
reusing the helper internally after feed-tangent construction. This refactor must be behavior-
preserving and rerun the existing Gate-3A regression suite.

The generic `pt_stability.hpp` does not need SW92 family knowledge and should remain unchanged.

### 9.1 Candidate phases must be explicit starts in both final family searches

For final review, both converged phase compositions must be appended as extra starts to **both** the
AQ and NA searches, in addition to caller-supplied final starts.

This directly tests whether changing only the family model at an already converged composition can
lower Gibbs. It also strengthens same-family/two-family phase-set review without inventing a
composition cutoff.

Start quotas for both required family searches must be validated before final property calls.

## 10. Final top-level acceptance semantics

A future maximum-two-phase Gate 3B entry should have the following phase-count flow.

### 10.1 Initial Gate-3A result: no instability found

If Gate 3A returns `no_instability_found` with a resolved lower feed family, Gate 3B may publish one
family-aware single-phase candidate at the normalized feed composition, preserving
`global_stability_proven=false`.

The finite search result is accepted only in the same qualified sense already used by the project:
“no instability found in the declared finite search,” never “globally proven stable.”

### 10.2 Initial Gate-3A result: indeterminate

Do not start a joint split and do not coerce the result to single phase. Preserve the feed-reference,
tie/nonsmooth/property/resource diagnostics.

### 10.3 Initial Gate-3A result: unstable

Generate the family-aware attempt plan from robust witnesses and solve fixed-family pairs.

If no admissible equation-converged pair exists, return `indeterminate` with attempt/resource
summary. A phase-disappearance attempt is not single-phase acceptance evidence.

### 10.4 Selected pair: final two-family search finds a negative witness

Return `phase_set_unstable` and retain the candidate plus family-tagged final witness evidence. Do
not publish it through `accepted_phase_set()`.

The first Gate 3B baseline may stop here rather than automatically launching another phase-addition
cycle. Automatic alternation between split and stability can be a later independently audited
extension.

### 10.5 Selected pair: final search indeterminate

Return `indeterminate`, retaining the equation-converged candidate and the complete AQ/NA final
search diagnostics. Equation convergence never overrides unresolved final stability.

### 10.6 Selected pair: both final family searches find no instability

The pair may be accepted only if all of the following are also true:

- EOS-component material balance passed;
- common chemical-potential residual passed;
- both phase fractions exceed the disappearance threshold;
- phases are compositionally distinct;
- each phase's same-family minimum-root evaluation is smooth/conditioned;
- each assigned family passed the explicit lower-envelope family-dominance/tie check;
- selected candidate Gibbs is not resolved above the original lower-envelope feed Gibbs;
- no distinct candidate phase-set Gibbs tie remains after deduplication;
- both final family searches used the same recorded common tangent and completed with
  `no_instability_found`.

Even then:

```text
global_stability_proven = false
```

because both initial and final stability searches remain finite multistart searches.

## 11. Why `iterate_pt_split` / `solve_pt_vle` are not the Gate 3B solver

The current VLE code contains valuable tested mathematics but also specific semantics that do not
fit the asymmetric family problem:

- it requests `PtPhaseRole::liquid_candidate` / `vapor_candidate`;
- the SW92 one-family adapter maps those roles to low-Z/high-Z candidate branches;
- it uses relative Z separation as part of distinct-phase acceptance;
- its result names phases liquid/vapor and does not carry `SwPhaseFamily`;
- its full driver assumes one model/reference provider for both phases and one family-blind final
  stability provider.

Therefore Gate 3B must not call it and reinterpret its roles as AQ/NA.

The following pieces are reusable at the mathematical/implementation level when their semantics are
kept explicit:

- `solve_rachford_rice` as a private fixed-K two-phase material-balance solver;
- compensated summation utilities;
- no-floor active-support rules;
- the log-K successive-substitution/Gibbs-line-search derivation for a fixed family pair;
- convergence, failure and resource-accounting patterns;
- `test_pt_stability_against` through the SW92-specific two-family common-tangent helper.

The first implementation may duplicate a small amount of fixed-pair state machinery rather than
prematurely refactor stable PR76/VLE production code. Shared extraction should be considered only
when equivalence tests prove that a model-independent abstraction is genuinely cleaner.

## 12. Fixed-molality and component-balance boundary remains unchanged

Every phase property call, family comparison, split iteration and final stability search uses the
same externally supplied

```text
NaCl molality [mol NaCl / kg H2O].
```

Only the ordered EOS-component vector is conserved by Gate 3B.

No accepted result may claim:

- conservation of total NaCl inventory;
- salt partitioning between phases;
- automatic molality change with phase fraction;
- a closed electrolyte equilibrium.

The fixed molality and this limitation must remain in result metadata/documentation. A conserved-salt
state formulation is a separate thermodynamic problem.

## 13. Independent validation required for the implementation

The Gate 3B implementation must add a new independent reference path that does not import production
code. Existing `Decimal(80)` SW92 infrastructure may be used as a pattern, but the reference must
independently rebuild the required equations.

Minimum evidence set:

1. **Fixed-family-pair binary equations.** Independently solve material balance plus common reduced
   chemical potentials for representative `AQ+AQ`, `AQ+NA`, and `NA+NA` assignments where such
   interior roots are representable. Do not assume beforehand that `AQ+NA` is the lowest-Gibbs pair.
2. **Candidate Gibbs ranking.** Compare independently computed total reduced Gibbs across distinct
   converged family assignments and verify guarded lowest-candidate selection.
3. **Family-dominance at each phase.** At every reference phase composition independently evaluate
   AQ and NA minimum-root Gibbs ordering. A dominated or tied assigned family must not be accepted.
4. **Final common tangent.** Independently construct `d_i` from a converged pair and evaluate
   prescribed AQ/NA TPD values against that same tangent.
5. **Final-instability rejection.** Include a converged fixed-family assignment that satisfies its
   pair equations but is rejected when the opposite family/lower-envelope final search supplies a
   robust negative witness, if a traceable reference state exposes such a case.
6. **Stable single-phase path.** Reuse an already audited Gate-3A stable feed to verify family-aware
   one-phase publication without running a forced split.
7. **Phase disappearance.** Construct a traceable small-phase state from an independently solved
   coexistence anchor and verify that balance/chemical-potential tolerances pass before disappearance
   status is emitted; it must not become accepted single phase.
8. **Candidate equivalence/deduplication.** Different witness seeds and slot-swapped representations
   that converge to the same family-aware phase set must not create a false Gibbs-tie
   indeterminate state.
9. **True distinct candidate tie handling.** If a deterministic synthetic-test fixture can create a
   distinct near-tie without inventing physical claims, verify `indeterminate` rather than
   attempt-order selection.
10. **Family-envelope tie at a candidate phase.** A synthetic exact/near AQ/NA family tie at an
    actual candidate composition must block accepted family identity.
11. **Component permutation invariance.** Permuting the ordered snapshot/feed must permute
    compositions/metadata without changing the physical set or status.
12. **Resource/failure isolation.** Invalid options/starts are preflighted before thermodynamic work;
    family/root/property failures remain attributable; split and final AQ/NA evaluation counts are
    explicit.
13. **Same-family pair support.** At least one regression must demonstrate that the solver does not
    require one AQ and one NA phase.
14. **Public-header self containment and hosted cross-platform CI.** GCC Debug + ASan/UBSan, Clang
    Release and MSVC Release, plus independent Decimal regeneration.

Reference values are model numerical anchors, not experimental validation and not a global phase-
equilibrium certificate.

## 14. Proposed algorithm/result identities

Keep the already merged identities unchanged:

```text
thermodynamics profile:
  SW92/corrected-original/PR76-base/NaCl-molality

parent equilibrium profile:
  SW92-equilibrium/xu-asymmetric-gibbs/v1

Gate 3A finite stability convention:
  SW92-equilibrium/xu-asymmetric-gibbs/stability-foundation/v1
```

Reserve a distinct Gate 3B numerical convention such as

```text
SW92-equilibrium/xu-asymmetric-gibbs/max2-logK-SSI-RR/common-tangent/v1
```

The exact string should be fixed in the implementation PR before public code is merged and tested as
metadata. `max2` is an algorithm capability, not a thermodynamic assertion that SW92/Xu can never
admit three phases.

## 15. Recommended staged implementation order

Gate 3B should not appear as one large patch. The safest sequence is:

### Gate 3B.1 — family-aware fixed-pair primitive

Implement first:

- SW92 family-aware candidate/result types without generic phase-set projection;
- an SW92 two-family imposed-common-tangent stability helper, preserving Gate-3A behavior;
- fixed-family-pair state evaluation with `F0,F1` explicit and equality allowed;
- private Rachford-Rice mapping to phase0/phase1 terminology;
- log-K SSI + Gibbs/residual line search;
- balance/fugacity/disappearance/composition-distinction/root-smoothness semantics;
- post-convergence per-phase AQ/NA family-dominance/tie check;
- independent fixed-pair Decimal(80) anchors and affected Gate-3A regression.

Do **not** yet orchestrate all witnesses or publish an accepted overall equilibrium if this primitive
has not independently passed review.

### Gate 3B.2 — witness orchestration and candidate selection

Then add:

- Gate-3A witness -> family-pair seed planning, including same-family alternatives;
- family-neutral resource ordering;
- attempt accounting;
- equivalent-candidate deduplication;
- guarded reduced-Gibbs selection and distinct-candidate tie semantics;
- one-phase flow for initial `no_instability_found`.

### Gate 3B.3 — final acceptance

Finally add:

- common tangent from the selected equation-converged pair;
- both candidate phase compositions as starts in both final family searches;
- pair-vs-feed Gibbs check;
- final two-family stability acceptance/rejection;
- authoritative family-aware `accepted_phase_set()`.

Only after Gate 3B.3 passes independent binary references, boundary/failure tests and hosted CI may
the repository claim a **maximum-two-phase, jointly material-balanced SW92/Xu candidate phase set**.
It must still state that stability is finite-search only and salt inventory is not conserved.

## 16. Explicitly deferred work

Gate 3B does not include:

- three-phase split or phase-addition iteration;
- a rule forcing one water-rich AQ phase plus one water-poor NA phase;
- universal physical phase labeling;
- salt-inventory conservation/electrolyte components;
- SW92 derivatives or implicit flash sensitivity;
- physics closure;
- interval-Newton/global optimization;
- changing the SW92 corrected-original thermodynamic parameters;
- replacing Whitson Profile A.

Profile A remains the evidence-aligned engineering compatibility path. Profile B/Gate 3B is a
separately identified generalized lower-envelope research equilibrium.

## Final decision

Gate 3A has removed the scientific blocker on defining a joint two-phase nonlinear problem, but
Gate 3B must preserve the discrete family model explicitly rather than forcing it through legacy
liquid/vapor/root roles.

The design is approved to proceed in the staged order above, beginning with **Gate 3B.1
family-aware fixed-pair primitive**. The key invariants are:

- one common EOS-component material balance;
- one common reduced chemical-potential gauge;
- explicit phase family per candidate phase, with repeated family allowed;
- minimum-Gibbs root inside each fixed family;
- explicit lower-envelope family-dominance/tie checks at converged phases;
- no density/root-order family inference;
- guarded Gibbs candidate selection only after duplicate phase sets are removed;
- final AQ and NA searches against one common converged-pair tangent;
- no accepted result when final asymmetric stability is negative or indeterminate;
- `global_stability_proven=false`;
- fixed molality remains external and salt is not conserved.
