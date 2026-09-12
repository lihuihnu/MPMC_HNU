# CPA PT VLE and maximum-three-phase flash baseline

## Scope

This increment connects the existing CPA PT phase-property and stability kernels to the model-neutral flash primitives. It implements:

- CPA vapor-liquid PT split through the existing generic `solve_pt_vle(...)` route;
- fixed three-phase CPA equilibrium through the existing generic `iterate_pt_three_phase(...)` / generalized-RR3 route;
- two-phase final-instability `2 -> 3` orchestration;
- three-phase disappearance followed by a **fresh** two-phase resolve for `3 -> 2`;
- generic phase-set publication;
- `CpaPtFlashBackend` integration with the unified PT flash capability and transition contracts.

It does **not** add morphology classification, mathematical global-stability certification, frontend/service transport, flash sensitivities, or physics behavior.

## CPA phase/root semantics

CPA pressure roots are obtained from the existing finite density-root scan. Stability always evaluates the minimum-Gibbs mechanically admissible root at the current composition.

A fixed split or three-phase phase slot cannot use that rule for every slot, because doing so would collapse distinct candidate slots onto the same lower Gibbs envelope. The CPA split/max3 adapters therefore use numerical density-root sides:

```text
lower_density_admissible
upper_density_admissible
```

For the two-phase baseline, the historical generic candidate roles map as:

```text
liquid_candidate -> upper-density admissible root
vapor_candidate  -> lower-density admissible root
```

These names are numerical solver coordinates only. They are **not** a general CPA vapor/liquid morphology classifier, and a root index is never persisted as a permanent phase identity across composition changes.

Near-multiple topology, root-search limits, missing mechanically admissible roots and nonrepresentable properties remain explicit property failures/indeterminate states. No fallback root is silently substituted.

## Two-phase baseline

`solve_cpa_pt_vle(...)` reuses the generic sequence unchanged:

```text
feed TPD stability using CPA minimum-Gibbs all-root evaluation
    -> negative-TPD witness
    -> material-balanced RR/logK seed
    -> explicit density-side two-phase SSI
    -> material balance + fugacity equality
    -> common-tangent final TPD review using minimum-Gibbs all-root CPA
```

A converged pair is not accepted until the final TPD review returns no sampled instability. A phase-fraction endpoint is not converted to a single-phase answer by clipping/deleting a phase.

The structural two-phase regression uses an explicitly `synthetic_test` non-associating CPA/SRK-limit binary at `1 MPa, 160 K`. Its offline coexistence anchors are retained only as numerical regression data. A 90/10 phase mixture is used so the existing generic 10% witness seed reconstructs the intended basin without changing generic split tolerances or seeding rules.

An additional associating synthetic fixture contains an explicit association site/pair and is exercised through the CPA split phase-property adapter. Thus the production adapter is association-aware. The current full accepted two-phase topology regression, however, is **not** a physical associating-mixture validation.

## Fixed three-phase primitive

`CpaThreePhaseEvaluator` adapts the same CPA phase-property kernel to the generic three-phase solver. The generic unknowns remain two independent log-K vectors plus two independent mole phase fractions, with the third fraction dependent.

Acceptance of the fixed candidate requires the unchanged generic gates:

- generalized RR3 balance;
- total component material balance;
- two independent chemical-potential equality sets;
- pairwise phase-composition separation;
- no phase-fraction disappearance.

The structural ternary regression is an explicitly synthetic non-associating CPA/SRK-limit fixture. At its reference point, the exact structural start closes with chemical-potential and material-balance residuals well below the production thresholds. This proves the CPA adapter and generic RR3 coupling; it is not experimental CPA three-phase validation.

## Maximum-three-phase orchestration

`solve_cpa_pt_max3(...)` first executes the full CPA two-phase path. It enters three-phase orchestration only when the two-phase final common-tangent review reports `phase_set_unstable`.

Caller-supplied `CpaPtThreePhaseStart` values are initialization hints only. They cannot create phase-count evidence, cannot bypass the base two-phase instability gate, and cannot bypass final three-phase stability. Malformed starts are rejected before the scientific solve.

Continuation-hint attempts and automatic negative-TPD-witness attempts have independent bounded budgets, so exhausting one initialization source cannot suppress the other.

The `2 -> 3` route is:

```text
accepted equations for a two-phase candidate
    -> final two-phase TPD finds additional-phase evidence
    -> fixed three-phase generalized RR / chemical-potential solve
    -> three-phase common-tangent final TPD review
    -> publish three phases only if that review closes
```

## Three-phase disappearance and fresh neighbor

A small third phase is never deleted directly. If the fixed three-phase solver reaches `phase_disappearance`, the two surviving compositions are used only as starts for a fresh complete `solve_cpa_pt_vle(...)` call.

Only when that fresh neighbor independently completes:

```text
feed stability
+ two-phase equilibrium
+ material balance
+ fugacity equality
+ final common-tangent review
```

may max3 publish a two-phase result and `3 -> 2 accepted_target` transition evidence.

The structural regression locks this behavior at the default three-phase minimum-fraction threshold: a third-phase fraction of `1e-10` reaches the disappearance route while the base pair is still found incomplete, and the fresh two-phase neighbor closes. At smaller structural fractions the base two-phase solve may already close directly, in which case no synthetic three-phase attempt is manufactured.

## Unified backend

`CpaPtFlashBackend` is a coarse runtime adapter only; inner CPA root, association, stability and equilibrium loops remain strongly typed.

Its current capability publishes:

```text
phase counts: 1, 2, 3
1 -> 2: fresh_target_resolve
2 -> 1: detection_only
2 -> 3: fresh_target_resolve
3 -> 2: fresh_target_resolve
initial stability search: yes
final phase-set review: yes
fresh boundary neighbor resolve: yes
global_stability_proven: false
morphology_resolved: false
```

Focused conformance compares direct CPA max3 publication against runtime `PtFlashBackend` publication using the same configured starts and requires identical phase payloads. Transition evidence is retained without inventing phase morphology.

## Validation boundary

Hosted validation uses GCC Debug + ASan/UBSan, Clang Release and MSVC Release.

Current regression evidence covers:

- explicit lower/upper density-root side selection;
- complete structural CPA VLE stability/split/final-review chain;
- component permutation for the two-phase result;
- an association-aware CPA phase-provider smoke;
- fixed CPA three-phase generalized-RR/equilibrium closure;
- accepted `2 -> 3` with final all-root stability review;
- default-threshold `3 -> 2` disappearance with fresh accepted two-phase neighbor;
- ternary component permutation;
- direct max3 publication versus runtime backend equivalence;
- continuation hints not acting as phase-count evidence;
- malformed three-phase-start rejection;
- public-header self containment.

The structural phase-split/max3 fixtures are `synthetic_test`. They validate software architecture and numerical state-machine behavior, not experimental CPA accuracy. A traceable **associating** two-/three-phase reference remains a separate validation gate; missing literature parameters or phase results must not be guessed.

`global_stability_proven` remains false throughout: all TPD decisions are finite searches, not mathematical global proofs.
