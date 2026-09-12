# Model-neutral PT service boundary v1

## Scope

`mpmc::runtime::PtService` is the transport-neutral C++ application boundary in
front of the frozen model-neutral `mpmc::flash::PtFlashBackend` contract. It
provides configured-backend discovery, an explicit runtime component inventory,
ID-keyed PT requests, self-contained result provenance, variable-cardinality
accepted phase sets, and stable indeterminate/error semantics.

Public conventions are:

```text
PT/service-boundary/v1
PT/service-component-inventory/v1
PT/service-capability-discovery/v1
PT/service-request/v1
PT/service-result/v1
```

The module is a header-only C++20 target:

```cmake
add_subdirectory(modules/runtime)
target_link_libraries(my_worker PRIVATE mpmc::runtime)
```

It depends on `mpmc::flash`; no thermodynamic or flash module depends on
`runtime`. No Protobuf, gRPC, HTTP, browser, database, plugin ABI, or third-party
runtime dependency is introduced by v1.

## Scientific delegation boundary

The service chooses exactly one configured `PtFlashBackend` and calls its
existing `solve(p,T,z)` method once. It does not:

- evaluate or select an EOS root;
- run or repeat material-balance, fugacity, TPD, common-tangent, phase-count, or
  phase-transition logic;
- normalize, clip, fill, or fit a composition;
- infer liquid/vapor/aqueous morphology from vector position, density, branch,
  or compressibility factor;
- retry an indeterminate solve or silently fall back to another model.

The only numeric checks in `runtime` are transport-independent request-domain
and result-integrity checks: finite positive SI pressure/temperature, finite
mole fractions in `[0,1]`, a sum-to-one roundoff guard, array dimensions,
finite returned scalars, and consistency with the discovered backend snapshot.
They do not constitute a second scientific acceptance calculation.

## Backend registry and capability discovery

Construction takes non-owning pointers to already configured backends. The
service snapshots and validates each capability and rejects null backends,
duplicate backend IDs, invalid capabilities, or configured resource-limit
violations before accepting requests. The backends and all objects referenced by
them must outlive the service.

`discover_capabilities()` returns an immutable ordered snapshot. Each
`PtServiceBackendDescriptor` contains:

- a service-local `configured_backend_id` selecting this exact registered
  instance;
- the complete frozen `PtFlashBackendCapability`, including backend/model/
  algorithm/publication/configuration identities;
- dataset ID and revision;
- supported phase counts and transition capability;
- configured scalar settings and search/review limitations;
- provider phase-metadata namespace, if any;
- `PtRuntimeComponentInventory`, with one stable component ID and canonical
  `feed_index` per backend component.

The inventory deliberately does not invent display names, formulas, aliases, or
chemical identity records because the frozen backend capability currently
provides only stable IDs. A later inventory version may add sourced presentation
metadata without changing scientific model selection.

`configured_backend_id` is deliberately distinct from
`PtFlashBackendCapability::backend_id`: the latter identifies the adapter
implementation and can be shared by several dataset or scalar-configuration
instances. The worker assigns a unique configured ID to every registration, so
multiple PR76, SW92, or CPA parameter snapshots can coexist without losing
provenance.

The registry is immutable after construction and has no global registration
side effects. It adds no locking: concurrent calls are allowed only when the
selected configured backend and its referenced evaluator/model are safe for that
use, or when the worker supplies external serialization.

## Request contract

`PtServiceRequest` contains:

```text
configured_backend_id
pressure_pa
temperature_k
feed[] = (component_id, mole_fraction)
```

The feed may arrive in any order because identity is explicit. The service
rejects missing, duplicate, unknown, empty, over-limit, non-finite, negative, or
greater-than-one entries, then maps the unchanged values into the discovered
canonical backend order. It never trims an ID or changes a mole fraction.

The mole-fraction sum must differ from one by no more than
`64 * epsilon(double)`, matching the generic flash composition roundoff
contract. This is a dimensionless arithmetic guard, not a normalization
tolerance. A request outside it is rejected; a value inside it is passed through
unchanged. Any backend-owned roundoff-only representation in the returned feed
is reported as the evaluated feed.

`PtServiceLimits` bounds configured backend count, components per backend, and
identifier bytes. The optional [native gRPC process adapter](../runtime_grpc/README.md)
adds pre-deserialization and post-decode request limits, response size, gRPC
memory quota, solve concurrency admission, and deadline/cancellation
publication guards. These process controls remain outside the transport-neutral
C++ v1 boundary.

## Result, provenance, and variable phase count

`PtServiceResponse` is an envelope containing either a computation result or a
service error, never both. A computation result owns:

- the discovered backend descriptor used for this solve;
- generic backend-result, phase-set, and transition convention IDs;
- provider result convention;
- evaluated `p`, `T`, and feed in canonical component order;
- zero or more **accepted** phases;
- the backend's transition report, diagnostic, global-stability flag, and
  morphology-resolution flag.

Every accepted `PtServicePhase` contains its diagnostic vector index, molar
phase fraction, optional `Z`, provider branch/smoothness diagnostics, and one
entry per canonical component with mole fraction and copied `ln(phi)`. Optional
role/family strings remain opaque inside the backend-declared metadata namespace.
Neither the phase index nor provider branch is a universal physical phase label.

There is no fixed one-/two-/three-phase union. The accepted phase count is
`response.result->phases.size()` (or `accepted_phase_count()`), and must be one of
the selected backend's discovered supported counts. Backend candidates that were
not accepted are intentionally not placed in this phase vector.

The provenance snapshot makes each returned calculation self-describing. It
records the exact backend, model, algorithm, publication and configuration
profiles; dataset/revision; scalar configuration; ordered component inventory;
generic convention chain; and provider result convention. It is provenance and
configuration identity, not proof that a dataset or finite search is globally
valid.

## Outcome and error semantics

`PtServiceOutcome` separates scientific decision states from service failures:

| outcome | result | accepted phases | service error | meaning |
| --- | --- | ---: | --- | --- |
| `accepted` | present | 1..N | absent | Backend accepted the phase set under its declared finite-search contract. |
| `phase_set_unstable` | present | 0 | absent | Final review found evidence against the candidate; no candidate is promoted. |
| `indeterminate` | present | 0 | absent | The backend could not make an accepted decision; provenance, transition evidence and diagnostic remain available. |
| `error` | absent | 0 | present | Request, dispatch, backend contract, or execution failed before a publishable computation response existed. |

In particular, `indeterminate` is not converted to accepted, HTTP success/failure
is not decided in this C++ layer, and an unchanged request is not automatically
retried. A future wire adapter should preserve this envelope rather than map
`indeterminate` to an empty successful phase set or an opaque transport failure.

`PtServiceErrorCode` has stable categories:

| category | codes |
| --- | --- |
| selection | `invalid_configured_backend_id`, `configured_backend_not_found` |
| PT state | `invalid_pressure`, `invalid_temperature` |
| inventory/composition | `component_count_mismatch`, `invalid_component_id`, `duplicate_component`, `unknown_component`, `invalid_mole_fraction`, `composition_not_normalized` |
| delegated backend | `backend_rejected_request`, `backend_contract_violation`, `backend_execution_failure` |

Bad allocation propagates rather than being mislabeled as an ordinary request or
scientific outcome. Invalid service construction is also separate: it throws
`PtServiceConfigurationError` with a configuration error code during worker
startup.

## Validation

The focused regression covers:

- ordered runtime component inventories and capability lookup;
- actual PR76, SW92 Profile-C, and CPA backend discovery through one registry;
- ID-keyed request permutation into the backend's canonical order without value
  repair;
- complete result provenance and opaque provider metadata preservation;
- accepted one-, two-, and three-phase result vectors;
- suppression of non-accepted diagnostic candidates;
- distinct `phase_set_unstable`, `indeterminate`, request-error, contract-error,
  and backend-execution paths;
- registry/configuration limits and public-header self containment;
- GCC Debug with ASan/UBSan, Clang Release, and MSVC Release hosted CI.

The scripted backend used by boundary tests is explicitly a software fixture; it
is not experimental data or a thermodynamic validation. Existing backend suites
remain the authority for material balance, fugacity equality, final common
tangent, physical regression, and model-specific provenance.

## Protobuf/gRPC-Web and frontend handoff

The separately versioned [`mpmc.runtime.v1` wire contract](../../api/README.md)
and generated browser client now map discovery, requests, variable accepted phase
sets, provenance, transition evidence and the response envelope. The frontend:

1. completes capability discovery before enabling solve;
2. builds read-only component identities from the selected runtime inventory;
3. submits only configured backend ID, Pa, K and ID-keyed mole fractions;
4. renders the variable accepted phase vector and provider-opaque metadata;
5. shows `phase_set_unstable`, `indeterminate`, `PtServiceError`, and gRPC/wire
   failures as distinct states without synthesizing phase data.

The C++ `runtime` target remains transport-neutral and has no Protobuf/gRPC
dependency. [`runtime_grpc`](../runtime_grpc/README.md) now provides an optional
thin native adapter and the repository provides a development/CI Envoy edge, but
no model-configured production worker, TLS/auth policy, or deployed endpoint.
Any worker still constructs its backends and registry outside the adapter; only
`PtService` may dispatch a solve.
