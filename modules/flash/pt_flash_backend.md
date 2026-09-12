# Unified PT flash backend and capability contract

## Scope

`mpmc/flash/pt_flash_backend.hpp` is the model-neutral top-level PT flash boundary used to select a configured thermodynamic backend without duplicating its scientific solver. Runtime polymorphism occurs once at the coarse `solve(p,T,z)` boundary; property, stability and equilibrium inner loops remain strongly typed.

Current adapters:

- `Pr76PtFlashBackend` v2: PR76 maximum-three-phase route;
- `Sw92ProfileCPtFlashBackend`: SW92 Profile-C boundary-aware 1/2/3-phase route.

CPA is not implemented yet. A future CPA adapter must satisfy the same backend, phase-set and transition contracts after its thermodynamic/stability/equilibrium kernels are independently validated.

## Capability snapshot

Each configured backend instance owns a `PtFlashBackendCapability` containing:

- backend/model/algorithm/publication/configuration identities;
- dataset/revision and ordered component IDs;
- supported phase counts;
- optional opaque scalar configuration provenance `(id,value,unit)`;
- initial-stability, final-review and fresh-boundary-resolve capabilities;
- model-neutral `PtPhaseTransitionCapability` edges;
- explicit `global_stability_proven`;
- optional provider-specific phase-metadata namespace.

These are implemented-route capabilities, not statements that a thermodynamic model can never admit another phase.

Current high-level capabilities are:

| backend | phase counts | initial stability | final review | fresh boundary neighbor | global proof |
| --- | --- | --- | --- | --- | --- |
| PR76 PT max3 | 1, 2, 3 | yes | yes | yes for 3→2; 2→1 remains detection-only | no |
| SW92 Profile-C | 1, 2, 3 | yes | yes, within declared finite topology/search contract | yes | no |

Current transition edges:

```text
PR76:
  1 -> 2  fresh_target_resolve
  2 -> 1  detection_only
  2 -> 3  fresh_target_resolve
  3 -> 2  fresh_target_resolve

SW92 Profile-C:
  1 -> 2  fresh_target_resolve
  2 -> 1  detection_only
  2 -> 3  fresh_target_resolve
  3 -> 2  fresh_target_resolve
  3 -> 1  fresh_target_resolve
```

Accepted lower-phase-count boundaries still require the target topology to be solved and reviewed fresh; a small phase fraction alone is not acceptance evidence.

## Configured-backend provenance

The generic request stays small, but a configured backend must not lose thermodynamic configuration that materially identifies the solve. The capability therefore carries `configuration_profile` plus optional scalar settings. Generic code validates these fields but does not interpret their scientific meaning.

Current configuration profiles:

```text
PR76/PT/max3/backend-configuration/v2
SW92/Profile-C/fixed-molality/backend-configuration/v1
```

PR76 currently has no generic scalar model setting. SW92 Profile-C records its prescribed NaCl molality as `mol/kg_H2O`.

A future wire/service API should still use typed PR/SW/CPA configuration rather than arbitrary string maps.

## Request

`PtFlashRequest` contains:

```text
pressure_pa
temperature_k
ordered feed
```

The ordered component identity comes from `backend.capability().component_ids`. The generic layer does not normalize, clip, reorder or repair the feed; existing backend validation remains authoritative.

Model-specific initialization/tuning belongs to the configured adapter. In particular, PR76 max3 may receive bounded caller-supplied three-phase continuation starts. They are initialization hints only and are ignored unless the existing two-phase final review has already established additional-phase instability.

## Result

`PtFlashBackendResult` owns:

- the capability/configuration snapshot;
- the generic `PtPhaseSetResult`;
- the model-neutral `PtPhaseTransitionReport`;
- provider result convention/provenance;
- optional backend-namespaced phase metadata;
- morphology-resolved flag.

The structural guard checks phase-count capability, feed dimension, provider metadata alignment and transition/report consistency before exposing an accepted phase set.

## Provider phase metadata

Provider metadata remain opaque to the generic layer.

SW92 Profile-C uses:

```text
SW92/Profile-C/phase-metadata/v1
role_id:   aqueous | nonaqueous_unclassified
family_id: aqueous | nonaqueous
```

H0/H1 remain `nonaqueous_unclassified`; AQ/NA family or cubic-root branch is not promoted to LV/LL morphology.

PR76 publishes no provider-specific phase morphology. Its numerical lower/upper admissible-root sides and the historical VLE candidate roles are solver coordinates, not universal liquid/vapor classification.

## PR76 max3 adapter

`Pr76PtFlashBackend` v2 calls `solve_pr76_pt_max3(...)`, which retains the existing PR76 VLE route as its first stage.

The max3 route is:

```text
feed stability
    -> existing two-phase solve
    -> existing two-phase final common-tangent review
    -> if stable: publish 1/2 phase
    -> if unstable: additional-phase evidence
         -> fixed three-phase material balance + chemical-potential solve
         -> three-phase common-tangent final review
         -> accept 3 phase OR retain unresolved/higher-topology evidence
```

A converged three-phase disappearance does not delete a phase. The two surviving compositions seed a fresh complete PR76 VLE solve; only a independently closed neighbor may publish `3 -> 2 accepted_target`.

Caller-supplied three-phase starts are bounded, validated continuation hints. They cannot trigger 2→3 and cannot bypass final stability. See [pr76_three_phase.md](pr76_three_phase.md) for the detailed numerical contract.

## SW92 Profile-C adapter

`Sw92ProfileCPtFlashBackend` continues to run the established boundary-aware Profile-C driver once and purely projects its authoritative publication and transition evidence. It does not rerun TPD, change topology, classify H morphology, or reinterpret AQ/NA semantics.

## Conformance requirement

An adapter is not accepted merely because it compiles. Focused regression requires direct model-specific and runtime-backend paths to reproduce the same generic publication.

Current coverage includes:

- PR76 max3 direct orchestration/publication versus runtime backend using the same explicit structural continuation hints;
- PR76 capability `{1,2,3}` and `2→3` accepted transition evidence;
- SW92 wet-binary direct publication versus runtime backend;
- physical SW92 Sample-6 authoritative three-phase publication;
- SW92 fresh-neighbor boundary projection;
- capability/configuration/transition structural guards;
- provider metadata preservation without morphology invention;
- public-header self containment;
- GCC Debug + ASan/UBSan, Clang Release and MSVC Release hosted builds/tests.

Any future backend, including CPA, should add equivalent direct-vs-adapter conformance before frontend/service discovery exposes it.

## Non-capabilities

This backend contract still does not provide:

- mathematical global-stability certification;
- guaranteed discovery of all PR76 three-phase basins from default finite automatic starts;
- morphology classification;
- CPA thermodynamics or flash;
- backend registry/plugin ABI;
- wire/service API or frontend model discovery;
- new flash sensitivities or physics behavior.

Those remain separate gates. The unified boundary now has a common 1/2/3-phase result and transition vocabulary for both PR76 and SW92, which is the prerequisite for later CPA and frontend integration.
