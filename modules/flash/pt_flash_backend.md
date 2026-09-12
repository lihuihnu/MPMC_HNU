# Unified PT flash backend and capability contract

## Scope

`mpmc/flash/pt_flash_backend.hpp` is the model-neutral top-level PT flash boundary used to select a configured thermodynamic backend without duplicating its scientific solver. Runtime polymorphism occurs once at the coarse `solve(p,T,z)` boundary; property, stability and equilibrium inner loops remain strongly typed.

Current adapters:

- `Pr76PtFlashBackend` v2: PR76 maximum-three-phase route;
- `Sw92ProfileCPtFlashBackend`: SW92 Profile-C boundary-aware 1/2/3-phase route;
- `CpaPtFlashBackend` v1: CPA stability + explicit-density-side VLE + generic-RR3 maximum-three-phase route.

## Service/frontend handoff freeze

The following model-neutral public semantics are frozen as the **v1 handoff contract** for the next service/frontend layer:

```text
PT/flash-backend/capability-and-dispatch/v1
PT/flash-backend/result/v1
PT/phase-transition-boundary/v1
```

Conformance tests instantiate PR76, SW92 Profile-C and CPA together and lock the common shape rather than relying on documentation alone.

For all three configured backends the frozen common capability includes:

- supported phase counts `{1,2,3}`;
- initial feed stability search;
- final accepted-phase-set review;
- fresh boundary-neighbor resolution capability;
- `1 -> 2  fresh_target_resolve`;
- `2 -> 1  detection_only`;
- `2 -> 3  fresh_target_resolve`;
- `3 -> 2  fresh_target_resolve`;
- `global_stability_proven=false`.

Provider differences remain explicit rather than being forced into a false common taxonomy:

- SW92 Profile-C additionally supports fresh `3 -> 1` and retains its provider-specific AQ/NA/physical-role metadata namespace;
- PR76 and CPA do not advertise `3 -> 1` and publish no provider morphology metadata;
- SW92 retains prescribed NaCl molality as configured scalar provenance;
- PR76 and the current CPA adapter have no generic scalar setting.

A future change that alters the meaning or required shape of the generic request/result/capability/transition contract must introduce a new generic convention version rather than silently changing v1. Backend-specific solver algorithms, parameter datasets and provider result/profile revisions may continue to evolve behind this boundary, but every adapter must continue to pass the frozen conformance suite or explicitly adopt a new generic contract version.

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

The wire/service API should use typed PR/SW/CPA configuration rather than arbitrary string maps.

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
- a three-backend freeze executable that compares the common phase-count/search/transition capability surface;
- PR76 and CPA `2→3` accepted transition evidence;
- CPA default-threshold `3→2` fresh-neighbor evidence;
- capability/configuration/transition structural guards;
- provider metadata preservation without morphology invention;
- public-header self containment;
- GCC Debug + ASan/UBSan, Clang Release and MSVC Release hosted builds/tests.

## CPA physical-validation status

CPA now has a traceable associating **two-phase** physical validation for methanol(2B) + water(4C) at 333.15 K using published CPA pure/association parameters, explicit CR-1 cross-association records, the published CR-1 binary interaction parameter, and Kurihara et al. experimental P-x-y data.

The full production CPA VLE route closes all five retained interior literature states and keeps active association. The focused five-point regression is frozen at:

```text
mean |Delta x(MeOH)| <= 0.01
mean |Delta y(MeOH)| <= 0.01
max  |Delta x(MeOH)| <= 0.015
max  |Delta y(MeOH)| <= 0.015
```

See [cpa_physical_validation.md](cpa_physical_validation.md) for the provenance and validation semantics.

This does **not** upgrade the existing synthetic CPA max3 structural regression into a physical three-phase validation. A physical associating VLLE oracle remains separate until a complete compatible literature chain supplies all CPA pure/association/BIP data and three phase compositions without inference or parameter guessing.

## Validation boundary and non-capabilities

This frozen v1 contract does not provide:

- mathematical global-stability certification;
- guaranteed discovery of every phase basin from finite automatic starts;
- morphology classification;
- a backend registry/plugin ABI;
- wire/service transport itself or frontend model discovery UI;
- new flash sensitivities or physics behavior;
- a traceable physical associating CPA three-phase oracle.

The common 1/2/3-phase result and transition vocabulary across PR76, SW92 Profile-C and CPA is now the stable boundary on which the service/frontend integration layer may be built.

The first transport-neutral consumer is now
[`mpmc::runtime::PtService`](../runtime/README.md). It registers configured
backend instances, exposes runtime component/capability discovery, and projects
the frozen result without adding EOS, flash, or acceptance logic. Wire transport
and the frontend connection remain separate increments.
