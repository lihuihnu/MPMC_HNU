# SW92 Xu-style asymmetric-Gibbs pre-implementation audit

## Status

This document is the write-before-code scientific audit for the reserved equilibrium profile

```text
SW92-equilibrium/xu-asymmetric-gibbs/v1
```

against repository baseline

```text
main@508ffc24bcfba4bee8a93878ce615ded14133af6
```

It does **not** implement a cross-family stability adapter, a joint phase-split solver, a salt
inventory balance, derivatives, or three-phase flash.

**Audit verdict:** the two existing SW92 AQ/NA family models can be placed on one common reduced
chemical-potential gauge and therefore support a scientifically defined **generalized asymmetric
model** in the sense of Xu, Haynes and Stadtherr (2005). The next numerical increment may implement
only the cross-family feed-reference and common-reference stability foundation. A production joint
AQ/NA split remains blocked until that foundation has independent numerical anchors and until the
result contract preserves explicit phase-family identity.

This profile remains distinct from the already implemented
`SW92-equilibrium/whitson-dual-model-observables/v1`. Xu-style asymmetric equilibrium is a new
thermodynamic construction over the two SW92 phase-family models; it is not claimed to be the
original 1992 SW92 flash algorithm.

## Sources reviewed

Primary asymmetric-model source:

- G. Xu, W. D. Haynes and M. A. Stadtherr, *Reliable Phase Stability Analysis for Asymmetric
  Models*, Fluid Phase Equilibria 235 (2005) 152-165, DOI `10.1016/j.fluid.2005.06.016`.
  Author preprint: `https://academicweb.nd.edu/~markst/xu_fpe2005.pdf`.

Supporting author review:

- M. A. Stadtherr, G. Xu, G. I. Burgos-Solorzano and W. D. Haynes, *Reliable computation of phase
  stability and equilibrium using interval methods*.
  Author copy: `https://www3.nd.edu/~markst/stadtherr-ijrs2007.pdf`.

SW92 model source:

- I. Soreide and C. H. Whitson, *Peng-Robinson predictions for hydrocarbons, CO2, N2, and H2S with
  pure water and NaCl brine*, Fluid Phase Equilibria 77 (1992) 217-240,
  DOI `10.1016/0378-3812(92)85105-H`.
- The repository implementation was audited against the project-supplied corrected PDF recorded in
  `modules/thermodynamics/sw92.md`, SHA-256
  `cb5b1d5034d78d934e887449ce0c89692d43431835d1371da5606d95b2c6bf58`.

Repository contracts and kernels reviewed:

- `modules/thermodynamics/include/mpmc/thermodynamics/sw92_pure.hpp`
- `modules/thermodynamics/include/mpmc/thermodynamics/sw92_mixture.hpp`
- `modules/thermodynamics/include/mpmc/thermodynamics/sw92_phase.hpp`
- `modules/flash/include/mpmc/flash/sw92_stability.hpp`
- `modules/flash/include/mpmc/flash/pt_stability.hpp`
- `modules/flash/include/mpmc/flash/pt_phase_set.hpp`
- `modules/flash/include/mpmc/flash/pt_split.hpp`
- `modules/flash/sw92_equilibrium_algorithm.md`

## 1. What Xu's asymmetric construction requires

For a symmetric model, tangent-plane stability uses one reduced Gibbs surface `g(x)`. Xu et al.
consider different Gibbs models for different phase types. The thermodynamic surface used for
stability is the lower envelope of the model-specific surfaces. For two models A and B,
conceptually,

```text
g(x) = min(g_A(x), g_B(x)).
```

At the feed `z`, the tangent/reference data must be taken from the model surface that is lower at
that same composition. Every model-specific tangent-plane-distance function must then use that
**same feed tangent plane**. For SW92 AQ/NA families the corresponding reduced stability functions
are

```text
D_AQ(w) = g_AQ(w) - common_feed_tangent(w)
D_NA(w) = g_NA(w) - common_feed_tangent(w)
D(w)    = min(D_AQ(w), D_NA(w)).
```

Xu et al. introduce a binary model variable in a pseudo-TPD formulation to represent model choice
inside one deterministic optimization problem. Their global guarantee comes from interval analysis,
not merely from the pseudo-TPD definition.

For MPMC_HNU, running one finite multistart search on each family against the same imposed tangent is
a valid finite approximation to the lower-envelope minimum because

```text
min_w min(D_AQ(w), D_NA(w))
  = min(min_w D_AQ(w), min_w D_NA(w)).
```

It is **not** equivalent to Xu et al.'s deterministic interval-Newton guarantee. The repository must
retain `global_stability_proven=false`.

## 2. Common chemical-potential reference audit — PASS

### 2.1 Family-independent pure-component basis

The current `Sw92Pure` kernel is independent of `SwPhaseFamily`:

- every selected component uses the same sourced `Tc`, `Pc`, and acentric factor;
- non-water components use the same PR76 pure alpha in AQ and NA;
- water uses the same SW92 Eq. (9) alpha at the same `T` and fixed NaCl molality in AQ and NA;
- every component uses the same printed PR `a_c`, `b`, and repository gas-constant convention.

`SwPhaseFamily` first appears in the **mixture interaction** path. The family changes selected
`k_ij`, not the component identity or pure-component reference basis.

A particularly strong invariant follows at every pure-component simplex vertex. All cross terms
vanish, so AQ and NA reduce to the same pure-component cubic EOS and therefore the same admissible
root set and fugacity coefficients at fixed `p,T,molality`. This must become an independent
regression before cross-family production code is accepted.

### 2.2 Same fugacity-coefficient convention

Both families call the same `Sw92Phase::evaluate` Peng-Robinson mixture fugacity expression. At
fixed `p,T` the component chemical potential may therefore be written, up to the same
component-specific ideal/reference term in both families, as

```text
mu_i^F / (R T) = common_i(T,p) + ln(x_i) + ln(phi_i^F),
F in {AQ, NA}.
```

The common component term is linear in composition and is identical for AQ and NA because the
ordered component identities and pure-component basis are the same. It cancels in family-to-family
Gibbs comparisons at the same composition and in a consistently formed common tangent.

Thus the reduced activity coordinate already used by the generic stability machinery,

```text
d_i = ln(x_i) + ln(phi_i),
```

is a valid common **gauge** for AQ/NA comparison under the exact same `p,T,molality`, ordered
snapshot, model profile, and component support.

This conclusion is stronger than observing that the formulas have the same shape: the two family
models share the component basis and coincide at the pure-component limits.

### 2.3 Required implementation guard

`test_pt_stability_against` is algebraically suitable for a provider evaluated against externally
supplied `d_i`: internally it evaluates

```text
q_i(w) = ln(w_i) + ln(phi_i^trial(w)) - d_i.
```

Its current public comment describes the reference-state cancellation for the "same p,T/model".
Before it becomes a documented cross-family building block, that contract should be clarified to
"same p,T and compatible component reference gauge" without changing its numerical algorithm.
The SW92 adapter, not generic TPD, must own the proof that the gauge is compatible.

## 3. Gibbs-surface domain audit — PASS as a mathematical model; empirical scope LIMITED

The current AQ and NA kernels share the same ordered component simplex and evaluate the same PR
mixing/fugacity algebra. For a validated snapshot and fixed finite `p>0`, `T>0`, molality `>=0`, each
family is mathematically evaluable for normalized compositions for which the existing cubic/property
kernel remains representable. There is no composition-based AQ/NA cutoff in production.

Within one family, `Sw92FamilyStabilityEvaluator` already:

1. enumerates all mechanically admissible cubic roots;
2. evaluates each admissible root with that family;
3. selects the minimum-Gibbs root at the requested composition;
4. marks a same-family root Gibbs near-tie nonsmooth;
5. rejects inadequate root conditioning.

This defines the family surface needed by a future asymmetric stability search as

```text
g_F(x) = min over mechanically admissible roots r of g_(F,r)(x).
```

A higher-Gibbs root at the same `x` cannot represent a stable phase of that family because replacing
it by the lower root changes no component inventory and lowers Gibbs energy.

The important limitation is scientific, not algebraic. SW92 AQ and NA BIPs are phase-specific
empirical parameterizations. Evaluating each family over the whole active composition simplex is a
model continuation outside regions where experimental support may be strong. Missing declared
applicability bounds mean **unknown validity**, not global validation. Therefore a future Xu-profile
result is a generalized **model-internal equilibrium** under SW92/corrected-original; it is not an
experimental validity certificate.

No `x_water` clipping or family cutoff should be invented to hide this limitation. Such a cutoff
would change the Gibbs model and break the lower-envelope formulation.

## 4. Fixed NaCl molality audit — PASS as an external conditional coordinate

The current thermodynamic profile is explicitly

```text
SW92/corrected-original/PR76-base/NaCl-molality.
```

NaCl molality enters the water alpha and selected BIP correlations, but Na+ and Cl- are not EOS
components and are not part of the flash material-balance vector. A Xu-style joint calculation can
therefore be well defined **conditional on one externally fixed molality**:

```text
(p, T, z_EOS, m_NaCl) -> conditional asymmetric equilibrium.
```

The same molality must be used in every AQ and NA feed/trial/phase property call. The accepted
component balances may cover only the ordered EOS components.

This profile must not claim:

- conservation of total NaCl inventory;
- redistribution of salt among phases;
- automatic change of molality as water moves between phases;
- equivalence to a closed H2O/NaCl/gas equilibrium problem.

A later conserved-salt formulation is a different state model and requires a separate thermodynamic
and variable-coordinate audit.

## 5. SW92 thermodynamic-consistency caveat — retained, not "fixed"

The project-supplied 1992 source audit records the authors' warning that using different
phase-specific attraction/mixing behavior is thermodynamically inconsistent in the ordinary
symmetric-EOS sense. The original paper intentionally uses distinct AQ and NA interaction
parameterizations for mutual-solubility engineering predictions.

The Xu construction does **not** make that warning disappear and must not be described as a
correction to original SW92. Instead, Profile B explicitly defines a new asymmetric model whose
system Gibbs surface is the lower envelope of two separately defined empirical family surfaces.

Required identity separation remains:

```text
thermodynamics: SW92/corrected-original/PR76-base/NaCl-molality
algorithm:      SW92-equilibrium/xu-asymmetric-gibbs/v1
```

Validation reports and public results must preserve both identities.

## 6. Feed-reference contract for the next numerical increment

The next implementation should stop at the asymmetric **stability foundation**.

For normalized feed `z` with active support `I`:

1. Evaluate `Sw92FamilyStabilityEvaluator` at `z` for AQ and NA independently. Each evaluation must
   resolve its same-family minimum-Gibbs cubic root.
2. If either family reference is unavailable, ill-conditioned, or same-family nonsmooth, the
   asymmetric feed reference is `indeterminate`.
3. Compare the two family reduced Gibbs values at the same feed. The common ideal-mixing term is
   identical, so family ordering may be computed as

   ```text
   Delta_g/(RT) = sum_i z_i [ln(phi_i^AQ(z)) - ln(phi_i^NA(z))].
   ```

4. Use a documented roundoff guard based on the magnitude of both family contributions. Do not use
   execution order, root index, or water fraction to break a near-tie.
5. If `Delta_g` is resolved, choose the lower family `F0` and construct the common feed tangent

   ```text
   d_i = ln(z_i) + ln(phi_i^F0(z)),  i in I.
   ```

6. If AQ and NA feed Gibbs values are tied within the numerical guard, report a **family-envelope
   nonsmooth feed** as `indeterminate`. Do not choose AQ-first or NA-first.

Zero-overall-feed components remain outside the active support under the current nonreactive closed
feed contract. No mole-fraction floor should be added.

## 7. Common-reference asymmetric stability contract

After a resolved feed tangent exists:

1. Run one AQ family TPD search against exactly that imposed `d_i`.
2. Run one NA family TPD search against exactly the same imposed `d_i`.
3. Use independent, explicit per-family resource options and report evaluations/trials separately.
   Defaults may be symmetric, but one family must not silently consume the other's budget.
4. Execute both required family searches; do not make AQ-first or NA-first execution order determine
   which evidence exists or which witness is retained.
5. Preserve family identity on every negative witness and lowest sampled point.

Combined status semantics should be:

- if either family supplies a robust negative TPD witness: `unstable`;
- otherwise, if either required family search is `indeterminate`: `indeterminate`;
- otherwise, when both finite searches report no instability: `no_instability_found`;
- always `global_stability_proven=false`.

A negative witness is sufficient evidence of instability even if the opposite-family search is
indeterminate. Absence of a negative witness is not a global proof.

The first MPMC_HNU implementation need not reproduce Xu's interval-Newton machinery or binary
variable solver. It must, however, name its search convention accurately as a finite two-family
multistart approximation to the Xu lower-envelope stability criterion.

## 8. Why the existing `solve_pt_vle` must NOT be reused as the joint solver

`solve_pt_vle` assumes two numerical candidate roles under one thermodynamic model/reference:

```text
liquid_candidate -> requested low-Z admissible branch
vapor_candidate  -> requested high-Z admissible branch.
```

That contract was deliberately designed so root order is not a universal physical phase identity.
A Xu SW92 joint split has a different discrete state variable: **phase family**. An actual phase
must carry `SwPhaseFamily::{aqueous,nonaqueous}` explicitly, while its cubic root remains a local
branch diagnostic inside that family.

Mapping

```text
liquid_candidate -> AQ
vapor_candidate  -> NA
```

would incorrectly equate density/root role with family identity and would reintroduce a heuristic
already rejected by the repository contracts.

Therefore the future joint split must have an asymmetric phase-family-specific nonlinear contract.
It may reuse generic material-balance utilities, convergence ideas, and diagnostics, but not the
legacy two-role provider API as the scientific definition of the joint problem.

At each phase composition the family property evaluation should use that family's mechanically
admissible minimum-Gibbs root. A same-family root switch/near-tie is nonsmooth and must remain
indeterminate unless a separately validated treatment is added.

## 9. Result/publication contract blocker

The current generic `PtCandidatePhase` stores composition, phase fraction, activity/root branch, and
optional Z, but it has no thermodynamic **phase-family/model identity** field.

Profile B requires every actual phase to retain explicit AQ/NA family identity. Consequently, the
future joint solver must not publish an accepted generic `PtPhaseSetResult` if doing so loses that
identity.

The smallest safe first implementation is a SW92-specific asymmetric result containing, for each
candidate phase:

```text
SwPhaseFamily family
mole phase fraction
composition
selected root/ln(phi)/Z diagnostics
```

and a phase-set status/acceptance guard equivalent in rigor to the generic contract. A later,
separate architecture audit may decide whether generic `PtCandidatePhase` should gain a
model-independent phase-model metadata field and whether projection is then appropriate.

This result-contract issue is a **blocking prerequisite for joint-solver publication**, not for the
next cross-family stability-foundation increment.

## 10. Independent validation required before joint split

The next stability-foundation PR should add independent calculations that do not import production
code. Existing traceable CO2/H2O and CH4/H2O data are sufficient for this first evidence set; no new
experimental constants need to be invented.

Required anchors/invariants:

1. **Pure-component gauge equality:** for every selected component at fixed `p,T,m`, AQ and NA give
   identical pure-component root/property results because all cross BIPs vanish.
2. **Cross-family feed Gibbs ordering:** independently compute AQ and NA minimum-root reduced Gibbs
   values at prescribed binary feed states and verify the selected lower family.
3. **Common tangent identity:** at a resolved feed, the selected feed family must give `D_F0(z)=0`
   to numerical roundoff, while the nonselected family gives the positive feed-family Gibbs gap.
4. **Common-reference TPD anchors:** prescribe trial compositions and independently compute
   `D_AQ(w)` and `D_NA(w)` against the same `d_i`, including at least one robust negative witness.
5. **Exact/near family tie structural regression:** construct a synthetic-test snapshot/state where
   AQ and NA pair interactions coincide, then require feed-reference `indeterminate` rather than an
   execution-order family choice.
6. **Component permutation invariance.**
7. **Family-isolated property/resource failures** and explicit separate evaluation counts.
8. **Same-family root nonsmoothness propagation** into asymmetric-feed/reference diagnostics.

The reference implementation should use high precision (the current repository pattern is
stdlib-only `Decimal(80)` for binary SW92 anchors) and should independently rebuild the corrected
SW92 correlations, PR mixing, cubic roots, minimum-root Gibbs ranking, feed-family comparison, and
common-reference TPD.

These are model numerical anchors, not experimental validation and not a global stability proof.

## 11. Gates after this audit

### Gate 3A — next recommended implementation

Implement only:

- a family-tagged asymmetric feed-reference result;
- robust AQ/NA feed Gibbs comparison and tie semantics;
- one common reduced feed tangent;
- two family-specific finite TPD searches against that tangent;
- combined asymmetric stability status and family-tagged witnesses;
- independent Decimal(80) anchors listed above;
- focused GitHub-hosted GCC/Clang/MSVC CI.

Do **not** solve a joint phase split in Gate 3A.

### Gate 3B — joint two-family, maximum-two-phase split

Only after Gate 3A review, implement a new asymmetric joint split that enforces one common material
balance and common component chemical potentials with explicit phase families. The first solver may
advertise maximum phase count 2, but that is an algorithm capability, not a claim that the
thermodynamic model cannot admit a third phase.

A converged pair is accepted only after final AQ **and** NA common-reference stability searches find
no robust negative witness within their declared finite search scope.

### Gate 4 — three-phase

Only after independent joint two-family binary/two-phase validation should the same thermodynamic
formulation be generalized to three phases. The likely physically relevant pattern may include one
AQ phase and two distinct NA phases, so three-phase design must not assume one phase per family.

## Final decision

The scientific preconditions for defining `SW92-equilibrium/xu-asymmetric-gibbs/v1` are satisfied
**conditionally** under the existing fixed-molality thermodynamic profile:

- common component reference gauge: **PASS**;
- mathematical AQ/NA Gibbs surfaces on the active simplex: **PASS**;
- empirical validity across the whole searched simplex: **LIMITED / model-extrapolation caveat**;
- fixed molality as an external state coordinate: **PASS, no salt-conservation claim**;
- original SW92 consistency warning: **retained by distinct asymmetric algorithm identity**;
- existing generic TPD algebra for an imposed common tangent: **reusable after SW92-specific
  reference-gauge contract and tests**;
- existing `solve_pt_vle` as a joint asymmetric solver: **NO**;
- current generic phase-set publication without family metadata: **NO**.

Therefore the repository is scientifically ready for **Gate 3A asymmetric stability foundation**,
but not yet for a production joint AQ/NA phase-split solver.
