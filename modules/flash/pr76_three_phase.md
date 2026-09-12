# PR76 generic max-three-phase PT baseline

## Scope

This increment extends the existing PR76 PT flash route from a validated maximum of two phases to a **maximum-three-phase baseline**. It deliberately does not classify the three returned phase instances as vapor/liquid/L1/L2. The public generic phase slots are numerical equilibrium instances only.

The implementation has four layers:

1. model-neutral fixed-three-phase equilibrium primitive in `pt_three_phase.hpp`;
2. PR76 phase-property/root-side adapter and max3 orchestration in `pr76_three_phase.hpp`;
3. generic phase-set publication in `pr76_max3_phase_set.hpp`;
4. unified `PtFlashBackend` exposure through `Pr76PtFlashBackend` v2.

`global_stability_proven` remains false. All stability statements are finite-search results under the declared TPD search contract.

## Fixed-three-phase coordinates

Let phase 0 be only the numerical reference phase. Phases 1 and 2 are represented with two independent log-ratio vectors

```text
logK1_i = log(x1_i / x0_i)
logK2_i = log(x2_i / x0_i)
```

and two independent mole phase fractions `beta1`, `beta2`; `beta0=1-beta1-beta2`.

For fixed log-K values, the generalized three-phase Rachford-Rice solve reconstructs three compositions while enforcing the feed material balance. The nonlinear outer solve then drives two chemical-potential residual sets

```text
r1_i = mu0_i - mu1_i
r2_i = mu0_i - mu2_i
```

to the declared tolerance. Candidate acceptance additionally checks:

- generalized-RR residual;
- absolute and positive-feed relative material balance;
- positive finite phase properties;
- all three phase fractions above the disappearance threshold;
- pairwise log-composition separation.

The fixed-state primitive intentionally performs **no final phase-set stability search**. A converged fixed-three-phase state is only a candidate until the PR76 orchestration performs the common-tangent review.

## PR76 root handling

PR76 phase slots use only the numerical root-side labels

```text
lower_admissible
upper_admissible
```

where admissible roots have the existing mechanically stable PR76 root condition. These labels are not universal liquid/vapor identities and are not exposed as morphology.

Raw cubic root indices are not treated as persistent phase identities because the number/order of algebraic roots can change along an iteration. A caller-supplied composition with one admissible root may be attempted on both numerical sides so later root splitting does not silently impose a physical label.

## 2 -> 3 trigger and initialization

The max3 driver always begins by running the established `solve_pr76_pt_vle(...)` path. It does **not** start a three-phase solve merely because a caller supplied three phase guesses.

A three-phase route is eligible only when the existing two-phase final common-tangent review returns `phase_set_unstable`, i.e. it has retained negative-TPD evidence that the current pair is incomplete.

Two initialization sources are then allowed:

1. automatic: source two-phase compositions plus negative final-stability trial witnesses;
2. optional caller-supplied `Pr76PtThreePhaseStart` triples.

Caller-supplied triples are continuation/initialization hints only. Before any EOS call they are bounded and checked for:

- configured count/entry quotas;
- exact component dimension;
- finite normalized compositions under the existing roundoff contract;
- identical active-component support to the feed;
- a feasible phase-fraction seed.

They cannot create 2->3 evidence, cannot bypass fixed-three-phase equations, and cannot bypass final stability. This is useful for PT continuation or a previous nearby three-phase state while keeping phase-count decisions thermodynamic rather than history-driven.

The current automatic witness seed remains a finite local initialization strategy; failure to converge it is reported as unresolved/indeterminate rather than repaired by weakening tolerances.

## Three-phase final review

For a converged candidate, the driver constructs the common log-activity reference from the three computed phases and calls the existing all-root PR76 stability evaluator via `test_pt_stability_against(...)`.

The three-phase candidate is accepted only when that final finite search returns `no_instability_found`. Additional negative TPD evidence is retained as a higher-phase/wrong-candidate state. An indeterminate search remains indeterminate.

No finite multistart result is promoted to a proof of global stability.

## 3 -> 2 disappearance

A small phase fraction is not a deletion rule.

If the fixed-three-phase equations converge while one phase fraction is at or below the declared disappearance threshold, the primitive returns explicit `phase_disappearance` evidence. The max3 orchestration then:

1. keeps the two surviving compositions as initialization hints;
2. runs a **fresh complete** `solve_pr76_pt_vle(...)` at the same `p,T,z`;
3. requires that fresh neighbor to pass its own feed stability, phase split equations, material balance, fugacity equality and final phase-set review.

Only then may the topology be published as two phase and the generic transition report mark `3 -> 2 accepted_target`.

If that fresh neighbor does not close, the phase is not removed and the boundary remains unresolved.

## Generic publication and backend

`project_pr76_pt_max3_phase_set(...)` publishes the accepted 1/2/3-phase state through the existing `PtPhaseSetResult`. Maximum phase count is 3; no provider-specific morphology metadata are added.

`Pr76PtFlashBackend` v2 declares:

```text
supported phase counts: 1, 2, 3
1 -> 2: fresh_target_resolve
2 -> 1: detection_only
2 -> 3: fresh_target_resolve
3 -> 2: fresh_target_resolve
```

The legacy `project_pr76_pt_transition_report(PtSplitResult)` remains available for VLE-only callers; the backend itself uses the max3 transition projection.

## Validation

The tests intentionally separate algorithm verification from physical validation.

### Generic mathematical fixture

A manufactured thermodynamically consistent three-phase provider verifies exact generalized-RR material balance, chemical-potential equality and disappearance semantics. It is an algebra/software test only.

### PR76 structural fixture

A symmetric three-component PR76 dataset is tagged `SourceKind::synthetic_test` and requires `DataPolicy::allow_synthetic_tests`. It is used only to exercise a reproducible three-phase topology, component permutation, publication, final stability and backend conformance. It is **not experimental validation, parameter fitting, or evidence that the synthetic fluid represents a real system**.

An independently solved structural three-phase composition triple is stored as a numerical regression anchor. The exact fixed-three-phase candidate is tested directly before it is used as a continuation hint.

The same synthetic topology has a zero-share three-phase edge regression where the surviving pair must fresh-resolve through the complete PR76 VLE route before two-phase acceptance.

Hosted validation uses the project matrix:

- GCC Debug + ASan/UBSan;
- Clang Release;
- MSVC Release.

The existing unified backend conformance and affected downstream physics-inventory regression are also rerun. These tests establish implementation consistency, not experimental accuracy of a real PR76 three-phase system.

## Non-capabilities

This baseline does not provide:

- liquid/vapor/L1/L2 morphology classification;
- mathematical global phase-stability certification;
- guaranteed discovery of every three-phase basin from default automatic starts;
- CPA thermodynamics or flash;
- service/frontend transport;
- new physics/reservoir-flow behavior.

A subsequent gate should strengthen physical PR76 three-phase/reference validation and PT continuation/transition regressions before using this route as a general-purpose production three-phase benchmark for additional thermodynamic backends.
