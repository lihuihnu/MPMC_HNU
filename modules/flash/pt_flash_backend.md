# Unified PT flash backend and capability contract

## Scope

`mpmc/flash/pt_flash_backend.hpp` is the model-neutral top-level PT flash boundary used to select a configured thermodynamic backend without duplicating its scientific solver. Runtime polymorphism occurs once at the coarse `solve(p,T,z)` boundary; property, stability and equilibrium inner loops remain strongly typed.

Current adapters:

- `Pr76PtFlashBackend` v2: PR76 maximum-three-phase route;
- `Sw92ProfileCPtFlashBackend`: SW92 Profile-C boundary-aware 1/2/3-phase route;
- `CpaPtFlashBackend` v1: CPA stability + explicit-density-side VLE + generic-RR3 maximum-three-phase route.

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

| backend | phase counts | initial stability | final review | fresh boundary neighbor | global proof |
| --- | --- | --- | --- | --- | --- |
| PR76 PT max3 | 1, 2, 3 | yes | yes | yes for 3→2; 2→1 detection-only | no |
| SW92 Profile-C | 1, 2, 3 | yes | yes, within declared finite topology/search contract | yes | no |
| CPA PT max3 | 1, 2, 3 | yes | yes | yes for 3→2; 2→1 detection-only | no |

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

CPA:
  1 -> 2  fresh_target_resolve
  2 -> 1  detection_only
  2 -> 3  fresh_target_resolve
  3 -> 2  fresh_target_resolve
```

Accepted lower-phase-count boundaries require the target topology to be solved and reviewed fresh; a small phase fraction alone is not acceptance evidence.

## Configured-backend provenance

The generic request stays small, but a configured backend must not lose thermodynamic configuration that materially identifies the solve. The capability therefore carries `configuration_profile` plus optional scalar settings. Generic code validates these fields but does not interpret their scientific meaning.

Current configuration profiles include:

```text
PR76/PT/max3/backend-configuration/v2
SW92/Profile-C/fixed-molality/backend-configuration/v1
CPA/PT/max3/backend-configuration/v1
```

PR76 and the current CPA adapter have no generic scalar model setting. SW92 Profile-C records prescribed NaCl molality as `mol/kg_H2O`.

A future wire/service API should still use typed PR/SW/CPA configuration rather than arbitrary string maps.

## Request and result

`PtFlashRequest` contains positive finite pressure, temperature and the ordered feed. Ordered component identity comes from `backend.capability().component_ids`. The generic layer does not normalize, clip, reorder or repair the feed; existing backend validation remains authoritative.

`PtFlashBackendResult` owns:

- the capability/configuration snapshot;
- generic `PtPhaseSetResult`;
- model-neutral `PtPhaseTransitionReport`;
- provider result convention/provenance;
- optional backend-namespaced phase metadata;
- morphology-resolved flag.

The structural guard checks phase-count capability, feed dimension, provider metadata alignment and transition/report consistency before exposing an accepted phase set.

Model-specific initialization belongs to the configured adapter. PR76 and CPA max3 may receive bounded caller-supplied three-phase continuation starts. Such starts are initialization hints only: they are ignored unless the current two-phase final review has already established additional-phase instability, and they never bypass final three-phase stability.

## Provider phase metadata and morphology

Provider metadata remain opaque to the generic layer.

SW92 Profile-C uses its AQ/NA/physical-role namespace. H0/H1 remain `nonaqueous_unclassified`; AQ/NA family or root branch is not promoted to LV/LL morphology.

PR76 and CPA publish no provider-specific phase morphology. PR lower/upper admissible-root sides, CPA lower/upper density-root sides and historical VLE candidate roles are numerical solver coordinates, not universal phase classification.

## PR76 max3 adapter

`Pr76PtFlashBackend` calls the established PR76 max3 route:

```text
feed stability
 -> two-phase solve
 -> two-phase final common-tangent review
 -> if unstable: fixed three-phase equilibrium
 -> three-phase final common-tangent review
 -> accept 3 phase, or fresh-resolve a disappearance neighbor
```

A converged three-phase disappearance never directly deletes a phase. See [pr76_three_phase.md](pr76_three_phase.md).

## SW92 Profile-C adapter

`Sw92ProfileCPtFlashBackend` runs the established boundary-aware Profile-C driver once and purely projects its authoritative phase-set publication and transition evidence. It does not rerun TPD, change topology, classify H morphology, or reinterpret AQ/NA semantics.

## CPA max3 adapter

`CpaPtFlashBackend` calls `solve_cpa_pt_max3(...)`.

CPA stability uses the minimum-Gibbs mechanically admissible density root at each composition. Fixed VLE and fixed three-phase slots use explicit lower-/upper-density admissible root sides so distinct candidate slots do not collapse onto the same Gibbs-envelope root. Those density sides are numerical branch coordinates only.

The route is:

```text
CPA all-root/minimum-Gibbs feed stability
 -> explicit-density-side two-phase RR/logK solve
 -> all-root final two-phase TPD review
 -> if unstable: generic generalized-RR3 + chemical-potential solve
 -> all-root final three-phase TPD review
 -> accept 3 phase OR fresh-resolve a disappearance neighbor
```

A `3 -> 2` result is published only after the surviving compositions seed a fresh complete CPA VLE solve whose own final TPD review closes. See [cpa_flash.md](cpa_flash.md).

## Conformance requirement

An adapter is not accepted merely because it compiles. Focused regression requires direct model-specific and runtime-backend paths to reproduce the same generic publication.

Current coverage includes:

- PR76 max3 direct publication versus runtime backend;
- SW92 wet-binary / Sample-6 publication and boundary projection;
- CPA max3 direct publication versus runtime backend using the same structural starts;
- PR76 and CPA `2→3` accepted transition evidence;
- CPA default-threshold `3→2` fresh-neighbor evidence;
- capability/configuration/transition structural guards;
- provider metadata preservation without morphology invention;
- public-header self containment;
- GCC Debug + ASan/UBSan, Clang Release and MSVC Release hosted builds/tests.

## Validation boundary and non-capabilities

This contract does not provide:

- mathematical global-stability certification;
- guaranteed discovery of every phase basin from finite automatic starts;
- morphology classification;
- a backend registry/plugin ABI;
- wire/service API or frontend model discovery;
- new flash sensitivities or physics behavior.

CPA production code is association-aware, but the current complete CPA VLE/max3 topology regressions use explicitly synthetic non-associating SRK-limit fixtures; association is independently exercised in CPA phase-property/split-provider tests. Traceable associating two-/three-phase physical validation is a separate gate and must not be fabricated from missing parameters or phase data.

The unified boundary now has a common 1/2/3-phase result and transition vocabulary across PR76, SW92 Profile-C and CPA, which is the required backend foundation for later service/frontend integration.
