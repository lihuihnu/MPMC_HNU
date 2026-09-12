# PT phase-transition / boundary contract

## Scope

`mpmc/flash/pt_phase_transition.hpp` defines the model-neutral contract used to describe phase-count transition capability and evidence around PT flash boundaries.

It does **not** add or remove phases, solve a new equilibrium problem, rerun stability, classify morphology, or certify global stability. The existing model-specific solver remains authoritative. The transition layer records what that solver can prove about moving between phase counts and whether the target topology was actually re-solved.

Public convention:

```text
PT/phase-transition-boundary/v1
```

The central rule is:

> A small/disappearing phase is evidence that a neighboring topology must be reconsidered. It is never, by itself, permission to delete a phase from the source solution.

Likewise, a negative final stability trial proves that the current phase set is incomplete, but it does not by itself prove a specific replacement phase count.

## Two-layer model

### Capability edge

`PtPhaseTransitionCapability` contains directed edges

```text
(source_phase_count -> target_phase_count)
```

with one of two support levels:

- `detection_only`: the backend can identify evidence that the target topology should be considered, but does not own a complete fresh target solve/review path for that edge;
- `fresh_target_resolve`: the backend owns a path that solves/reviews the target topology and may publish it only after that path closes.

Every current transition edge sets `requires_fresh_target_solve=true`. A capability with duplicate/self/out-of-range edges is structurally invalid.

### Per-solve evidence

`PtPhaseTransitionReport` contains zero or more `PtPhaseTransitionEvidence` records. Each record stores:

- source phase count;
- optional target phase count;
- trigger;
- resolution;
- whether a fresh target solve was attempted;
- whether the target topology actually closed;
- provider evidence profile and diagnostic.

The target count is optional deliberately. For example, a lower-Gibbs trial against an accepted-looking two-phase common tangent proves that the pair is incomplete, but a PR76 VLE backend does not thereby prove that the replacement must be a three-phase state.

## Triggers

Current generic triggers are:

- `initial_stability_witness`: a feed/reference phase is unstable and a higher-count candidate path is entered;
- `final_phase_set_instability`: final phase-set/common-tangent review finds evidence against the current set;
- `phase_disappearance`: an equilibrium iterate reaches a phase-fraction endpoint;
- `provider_topology_witness`: a model-specific topology witness requests another phase count;
- `provider_boundary_route`: a model-specific boundary resolver routes to a neighboring topology.

These triggers do not define aqueous/liquid/vapor morphology.

## Resolutions

`accepted_target` is the only resolution that authorizes the target topology. The structural contract requires:

- a known target phase count;
- a declared edge with `fresh_target_resolve` support;
- `fresh_target_solve_attempted=true`;
- `target_topology_closed=true`.

Other resolutions remain explicitly non-authoritative:

- `target_resolve_required`: evidence exists but no fresh target solve has been performed;
- `target_resolve_failed`: a fresh target solve was attempted but did not close;
- `candidate_not_accepted`: a candidate topology was solved but did not pass its acceptance/review chain;
- `broader_topology_required`: the current phase set is known to be incomplete but the replacement phase count is not proven;
- `indeterminate`: the transition decision did not complete reliably.

## PR76 mapping

`project_pr76_pt_transition_report(const PtSplitResult&)` is a pure projection of the existing VLE result.

Current PR76 capability:

| edge | support | meaning |
| --- | --- | --- |
| `1 -> 2` | `fresh_target_resolve` | feed instability can lead to a freshly solved/reviewed two-phase VLE state |
| `2 -> 1` | `detection_only` | a converged endpoint/phase-disappearance attempt is retained as evidence; no phase is silently deleted |

Specific semantics:

- `two_phase_no_instability_found` -> `1->2 accepted_target`;
- a `phase_disappearance` split attempt -> `2->1 target_resolve_required`;
- `phase_set_unstable` retains the two-phase candidate as not accepted and publishes `broader_topology_required` with unknown target phase count.

Therefore this contract does **not** advertise a PR76 three-phase solver.

## SW92 Profile-C mapping

`project_sw92_profile_c_transition_report(const Sw92PhaseAssignedBoundaryAwareResult&)` is a pure projection of the existing Profile-C topology/boundary chain.

Current SW92 Profile-C capability:

| edge | support |
| --- | --- |
| `1 -> 2` | `fresh_target_resolve` |
| `2 -> 1` | `detection_only` |
| `2 -> 3` | `fresh_target_resolve` |
| `3 -> 2` | `fresh_target_resolve` |
| `3 -> 1` | `fresh_target_resolve` |

Examples:

- accepted W(AQ)+H0(NA)+H1(NA) after an additional-H witness -> `2->3 accepted_target`;
- three-phase H disappearance followed by a fresh W+H re-solve and H-side review -> `3->2 accepted_target`;
- W disappearance followed by a fresh no-W re-solve may close at one or two NA/H phases -> `3->1` or `3->2 accepted_target`;
- a fresh neighboring-topology solve that still finds the removed phase necessary remains `target_resolve_failed` and does not publish the lower phase count.

AQ/NA family, H0/H1 representation and physical roles remain provider-specific metadata. The generic transition contract does not infer LV/LL morphology.

## Backend integration

`PtFlashBackendCapability` owns a `transition_capability`, and `PtFlashBackendResult` owns a `transition_report`. The backend structural guard validates the report against the configured capability.

Runtime dispatch remains coarse-grained. No virtual dispatch is added to EOS properties, TPD iterations, RR, or fixed topology equilibrium iterations.

## Validation

The PT backend conformance suite covers:

- existing PR76 and SW92 direct-vs-backend phase-set equality;
- Sample-6 authoritative three-phase state as `2->3 accepted_target`;
- PR76 disappearance projection as `2->1 target_resolve_required` rather than accepted single phase;
- PR76 final phase-set instability with unknown target count;
- the existing physical SW92 Sample-6 W+H boundary fresh re-solve projected as `3->2 accepted_target`;
- a physical no-W neighboring topology that fails to close projected as `target_resolve_failed`;
- duplicate transition edges and detection-only edges mislabelled as accepted being rejected by structural guards;
- public-header self containment;
- GCC Debug + ASan/UBSan, Clang Release and MSVC Release hosted builds/tests.

## Non-capabilities

This increment does not implement:

- PR76 three-phase equilibrium;
- a new generic three-phase solver;
- CPA;
- a continuation/history state machine across a sequence of PT points;
- global stability certification;
- a service/wire protocol or frontend visualization of transition evidence;
- any new physics/discretization behavior.

The next algorithmic gate can now use this contract to add PR76 three-phase support without inventing another phase-boundary vocabulary.
