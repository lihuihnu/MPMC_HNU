# Generic PT candidate-phase / accepted phase-set contract

## Scope

This increment introduces a model-independent **representation contract** for PT flash phase sets. It does not change the validated VLE iteration, does not add a three-phase solver, does not introduce SW/CPA, and does not modify or extend `physics`.

The existing `PtSplitResult`/`Pr76PtSplitResult` VLE APIs remain available unchanged. New code may consume the generic representation through:

```cpp
#include <mpmc/flash/pt_vle_phase_set.hpp>
const auto generic = mpmc::flash::project_pt_vle_phase_set(legacy_solution);
```

or, for the existing PR76 path:

```cpp
#include <mpmc/flash/pr76_phase_set.hpp>
const auto result = mpmc::flash::solve_pr76_pt_phase_set(
    p_pa, t_k, z, evaluator);
```

`solve_pr76_pt_phase_set` runs the established `solve_pr76_pt_vle` path once and then performs a structural projection. It is not a second flash implementation.

## Why the contract is separate from `PtSplitState`

The current VLE iteration is intentionally specialized around two requested candidate roles and Rachford–Rice variables. Replacing its internal `liquid`/`vapor` fields with a dynamic phase vector would couple this representation change to mature two-phase numerics, converged-solution sensitivity, and downstream compatibility.

Instead, the generic contract is a new boundary type. Existing VLE internals remain fixed; future LLE or joint three-phase algorithms can produce the same boundary type without inheriting VLE-specific storage.

## Status is independent of phase count

`PtPhaseSetStatus` has only:

- `accepted`: the candidate set is the phase set accepted by the current finite-search contract;
- `phase_set_unstable`: a candidate set may exist, but the final stability review found lower-Gibbs evidence;
- `indeterminate`: the current calculation cannot make an accepted phase-set decision.

A successful phase count is therefore **not** encoded in the status. Use:

```cpp
result.accepted_phase_count()
```

or `result.accepted_phase_set()->phases.size()`.

This allows the same result shape to represent accepted one-, two-, or future three-phase states.

As everywhere in the current flash module, `global_stability_proven` remains false. `accepted` means accepted under the stated finite-search/numerical contract; it is not a global thermodynamic proof.

## Candidate phase semantics

`PtCandidatePhase` contains:

- `mole_phase_fraction`: mole phase fraction, never pore-volume saturation;
- `composition`: ordered component mole fractions;
- `activity`: `ln(phi)`, provider-local branch diagnostic, and smoothness flag;
- optional `compressibility_factor`.

The phase-vector index and `activity.branch` are **diagnostics, not universal physical identities**. The generic contract does not define `phase[0] == liquid`, `phase[1] == vapor`, or any aqueous label. Multiple phase instances may in future come from the same provider family or branch at different compositions.

`compressibility_factor` is optional because the generic stability provider does not require Z. For the existing accepted two-phase VLE projection it is copied exactly from the converged legacy pair; for an accepted single phase projected only from `StabilityPhase`, it is absent rather than recomputed.

## Candidate versus accepted set

`candidate_phase_set` may be present even when the status is `phase_set_unstable` or `indeterminate`. This preserves a converged pair for diagnostics without presenting it as accepted.

`accepted_phase_set()` returns a pointer only when status is `accepted`, the candidate set is present and nonempty, and its phase count does not exceed `capability.maximum_phase_count`; otherwise it returns null. This is only a structural publication guard, not a thermodynamic validator. The rvalue overload is deleted to prevent dangling pointers.

For legacy VLE projection:

| Legacy `PtSplitStatus` | Generic status | Candidate set | Accepted set |
| --- | --- | --- | --- |
| `single_phase_no_instability_found` | `accepted` | one phase | one phase |
| `two_phase_no_instability_found` | `accepted` | two phases | two phases |
| `phase_set_unstable` | `phase_set_unstable` | converged pair if available | none |
| `indeterminate` | `indeterminate` | converged pair if available | none |

Malformed manually constructed legacy results are handled conservatively: an accepted legacy status without the required reference/candidate is projected as `indeterminate`, never as fabricated success.

## No thermodynamic recomputation in the compatibility projection

`project_pt_vle_phase_set` performs no EOS/provider call, no root solve, no Rachford–Rice solve, no normalization, no phase relabeling, no stability search, and no equilibrium iteration. The legacy compositions, activity data and Z values are copied directly. Legacy VLE stores only the vapor mole phase fraction explicitly, so the projected first-phase fraction is the contract-defined complement `1 - beta_V`; this is structural conversion, not a new flash calculation.

The established PR76 metadata (`dataset_id`, `revision`, ordered component IDs, model profile, PT convention and root options) remains outside the model-independent payload and is retained by `Pr76PtPhaseSetResult`.

## Capability semantics

`PtPhaseSetCapability::maximum_phase_count` describes the **algorithm/entry point being called**. The current VLE projection reports `2`. This must not be interpreted as a statement that PR76, SW, CPA, or a fluid system can never admit more than two phases.

A future joint three-phase solver can report `3` while using the same `PtPhaseSetResult` and `PtCandidatePhase` types.

## Explicit non-capabilities

This increment does not implement LLE equations, joint three-phase residuals, generalized Rachford–Rice, phase switching, phase appearance/disappearance derivatives, SW, CPA, water models, frontend/service types, or any physics coupling. Existing two-phase sensitivity remains on the established VLE result contract and is not generalized here.
