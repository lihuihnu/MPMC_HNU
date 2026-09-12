# Unified PT flash backend and capability contract

## Scope

`mpmc/flash/pt_flash_backend.hpp` defines the model-neutral top-level PT flash backend boundary. It exists so application/service/frontend layers can select among configured thermodynamic backends without duplicating the scientific solver or hard-coding PR/SW-specific result types.

The contract deliberately does **not** replace the existing flash algorithms. Runtime polymorphism occurs once at the coarse `solve(p,T,z)` boundary; after dispatch, each concrete adapter calls its established strongly typed solver and all property/stability/equilibrium inner loops remain unchanged.

Current adapters:

- `Pr76PtFlashBackend`: existing PR76 VLE/stability path projected through `solve_pr76_pt_phase_set(...)`;
- `Sw92ProfileCPtFlashBackend`: existing SW92 Profile-C boundary-aware authoritative 1/2/3-phase path through `solve_sw92_profile_c_pt_phase_set(...)`.

CPA is not implemented yet. A future CPA adapter must satisfy the same backend contract after its own thermodynamic/stability/two-/three-phase kernels are independently validated.

## Capability snapshot

Each configured backend instance owns a `PtFlashBackendCapability` containing:

- backend ID;
- thermodynamic model profile;
- algorithm profile;
- publication profile;
- backend configuration profile;
- dataset/revision identity;
- ordered component IDs;
- explicitly supported phase counts;
- optional opaque scalar configuration provenance `(id,value,unit)`;
- whether the backend performs an initial stability search;
- whether it performs a final phase-set stability/review gate;
- whether it provides a fresh neighboring-topology re-solve for phase-boundary routing;
- whether mathematical global stability is proven;
- an optional backend-specific phase-metadata namespace.

The capability describes the implemented route. It is not permission to infer capabilities that the route does not publish. In particular, `maximum_phase_count` or support for `{1,2,3}` is an algorithm capability, not a proof that the thermodynamic model cannot admit another phase.

Current declared capabilities are:

| backend | phase counts | initial stability | final review | fresh boundary neighbor | global proof |
| --- | --- | --- | --- | --- | --- |
| PR76 PT VLE | 1, 2 | yes | yes | no | no |
| SW92 Profile-C | 1, 2, 3 | yes | yes, within declared finite topology/search contract | yes | no |

PR76 is therefore **not** advertised as a three-phase backend by this adapter. SW92 Profile-C keeps `global_stability_proven=false` even when it authoritatively publishes a phase set under its declared finite-search/topology contract.

### Configured-backend provenance

The common request intentionally stays small, but a configured backend must not silently lose thermodynamic configuration that materially identifies the solve.

The capability therefore carries a `configuration_profile` plus optional opaque scalar settings. The generic layer only validates finiteness, non-empty IDs/units and unique setting IDs; it does not interpret them.

Current profiles:

```text
PR76/PT-VLE/backend-configuration/v1
SW92/Profile-C/fixed-molality/backend-configuration/v1
```

PR76 currently publishes no scalar model setting through this layer. SW92 Profile-C publishes:

```text
id    = nacl_molality_mol_per_kg_water
value = configured prescribed molality
unit  = mol/kg_H2O
```

This preserves the physically relevant prescribed molality in the capability/result snapshot. It still does not turn arbitrary strings into the future service API: a versioned wire contract should expose PR/SW/CPA-specific configuration with explicit typed/`oneof` fields.

## Request

`PtFlashRequest` contains only:

```text
pressure_pa
temperature_k
ordered feed
```

The ordered component identity comes from `backend.capability().component_ids`.

This layer does not normalize, clip, reorder or repair compositions. Existing backend validation remains authoritative so wrapping a solver does not create a second input contract with different tolerances.

Model-specific configuration is owned by the configured adapter rather than placed into this generic request. For example, SW92 Profile-C adapter options own prescribed NaCl molality and the existing Profile-C solver options. The material scalar configuration is copied into the capability snapshot as described above; solver-tuning options are not flattened into a generic string map.

## Result

`PtFlashBackendResult` owns:

- the capability snapshot used for the solve, including configured-backend provenance;
- the existing generic `PtPhaseSetResult` without reinterpretation;
- provider result convention/provenance;
- optional backend-namespaced phase metadata;
- a morphology-resolved flag.

The structural guard checks that the capability itself is valid, including scalar configuration provenance, and that the generic result's advertised maximum phase count, feed dimension, supported accepted phase count and optional provider phase metadata are consistent with the capability snapshot.

### Phase metadata

Provider metadata is intentionally opaque to the generic layer.

SW92 Profile-C currently uses namespace:

```text
SW92/Profile-C/phase-metadata/v1
```

with the already-authorized role/family IDs:

```text
role_id:   aqueous | nonaqueous_unclassified
family_id: aqueous | nonaqueous
```

This does not create a universal aqueous/liquid/vapor taxonomy. H0/H1 remain `nonaqueous_unclassified`; AQ/NA family and cubic-root branch are not promoted to LV/LL morphology.

The PR76 adapter publishes no provider-specific phase metadata. Its legacy liquid/vapor **candidate roles** are numerical branch requests inside the existing VLE algorithm and are not promoted by the unified backend into universal physical morphology.

## PR76 adapter

`Pr76PtFlashBackend` owns a reference to the existing sequential `Pr76VleEvaluator` plus existing `PtSplitOptions`/optional start sets.

Its `solve()` calls exactly:

```cpp
solve_pr76_pt_phase_set(...)
```

which itself runs the established `solve_pr76_pt_vle(...)` once and only projects the owned result. The adapter introduces no second stability search, RR solve, EOS/root selection, convergence tolerance or publication decision.

The configured backend instance inherits the evaluator's sequential scratch semantics and must not be shared concurrently without external synchronization.

## SW92 Profile-C adapter

`Sw92ProfileCPtFlashBackend` owns a reference to the existing SW92 model plus prescribed NaCl molality and existing Profile-C options.

Its `solve()` calls exactly:

```cpp
solve_sw92_profile_c_pt_phase_set(...)
```

The established boundary-aware Profile-C driver remains the sole solver/orchestration source. The adapter does not rerun TPD, change topology, select roots, normalize compositions or classify H morphology.

## Conformance requirement

A backend adapter is not accepted merely because it compiles. The focused conformance regression requires the adapter path to reproduce the direct model-specific entry point exactly for the generic phase-set payload.

Current coverage includes:

- PR76 binary direct `solve_pr76_pt_phase_set(...)` versus runtime `PtFlashBackend::solve(...)`;
- SW92 wet binary direct Profile-C publication versus runtime backend solve;
- physical Mortezazadeh–Rasaei Sample-6 authoritative three-phase publication through the unified backend;
- capability identity/phase-count/stability/boundary/configuration flags;
- SW prescribed-molality configuration provenance;
- SW provider metadata preservation without H morphology invention;
- invalid capability, duplicate/non-finite scalar-setting and invalid PR feed guards;
- public-header self containment;
- GCC Debug + ASan/UBSan, Clang Release and MSVC Release hosted builds/tests.

Any future backend (including CPA) should add equivalent direct-vs-adapter conformance tests before it is exposed through service or frontend capability discovery.

## Non-capabilities

This contract does not yet provide:

- a PR76 generic three-phase solver;
- CPA thermodynamics, stability or flash;
- a generic 1↔2↔3 transition state machine beyond the capabilities already owned by each adapter;
- a backend registry/plugin ABI;
- a wire/service API;
- frontend model discovery;
- global-stability certification;
- new flash sensitivities or physics behavior.

Those are separate gates. The immediate purpose of this contract is to make the next phase-transition, PR three-phase, CPA and service/frontend work converge on one top-level backend boundary instead of adding more model-specific application APIs.
