# Model-neutral PT wire contract v1

## Scope and dependency direction

[`pt_service.proto`](proto/mpmc/runtime/v1/pt_service.proto) is the versioned
process/browser mapping for `mpmc::runtime::PtService`. Its Protobuf package is
`mpmc.runtime.v1`; its gRPC service is `PtFlashService` with two unary methods:

```text
DiscoverPtCapabilities
SolvePtFlash
```

The schema and browser adapter depend conceptually on the transport-neutral
runtime contract. Protobuf and gRPC types do not enter `runtime`, `flash`, or
`thermodynamics`, and the adapter does not evaluate EOS properties, choose roots,
run stability analysis, alter phase count, normalize feed, retry an indeterminate
solve, or select a fallback backend.

The repository now provides the schema, generated TypeScript descriptors, strict
browser mapping, gRPC-Web client, and an optional
[native C++ process adapter](../modules/runtime_grpc/README.md). It does not ship
a model-configured production worker or deployed endpoint. A worker constructs
its registered backends outside the adapter; the adapter projects to/from the
existing `PtService` and calls `PtService::solve` exactly once per admitted
`SolvePtFlash` RPC.

## Version and compatibility rules

The binary API version is carried by the package path and package name:

```text
api/proto/mpmc/runtime/v1/pt_service.proto
mpmc.runtime.v1
```

Messages also carry `mpmc.runtime.v1/PT-flash-service/v1` and the existing C++
convention strings. The browser rejects a response whose package-generated shape
is valid but whose convention chain does not match v1.

Within v1:

- existing field numbers and meanings are immutable;
- removed fields must be `reserved`, never reused;
- compatible additions use new field numbers and safe defaults/presence rules;
- an incompatible semantic or shape change requires `mpmc.runtime.v2`;
- enum zero values are `UNSPECIFIED` and are rejected at the application mapping;
- semantically required scalars use proto3 `optional` so absence is distinguishable
  from valid zero/false values;
- ordered inventories, feeds and phases use `repeated`, not `map`, because order is
  part of the runtime contract.

## C++ to Protobuf mapping

| Runtime C++ contract | Protobuf v1 mapping | Notes |
| --- | --- | --- |
| `PtService::discover_capabilities()` | `DiscoverPtCapabilitiesResponse.backends` | Immutable owning snapshot; configured IDs must be unique. |
| `PtServiceBackendDescriptor` | `PtBackendDescriptor` | Keeps the service selector distinct from adapter `backend_id`. |
| `PtFlashBackendCapability` | `PtBackendCapability` | Complete model/algorithm/publication/configuration/dataset/revision, phase-count, setting, transition and limitation snapshot. |
| `PtRuntimeComponentInventory` | `PtRuntimeComponentInventory` | Stable ID plus canonical `feed_index`; must exactly match capability `component_ids`. |
| `PtServiceRequest` | `SolvePtFlashRequest` | Pa, K, configured backend ID and ID-keyed feed; values are not repaired. |
| `PtServiceComputationResult` | `PtComputationResult` | Evaluated state, provenance, transition report, diagnostic and variable accepted phase vector. |
| `PtServicePhase::provider_branch` (`size_t`) | `optional uint64 provider_branch` | Browser maps the generated `bigint` to a decimal string; it is diagnostic, not physical identity. |
| component/phase indices and phase counts | `optional uint32` | Worker must checked-cast; current runtime limits are far below `uint32` maximum. |
| `PtServiceOutcome::{accepted, phase_set_unstable, indeterminate}` | `PtComputationResult.outcome` | Always the `result` oneof arm. Non-accepted states publish zero phases. |
| `PtServiceOutcome::error` + `PtServiceError` | `SolvePtFlashResponse.error` | Application/service failure arm, never a scientific outcome. |

Result provenance contains the complete backend descriptor plus backend-result,
phase-set, transition and provider-result conventions. Provider role/family IDs
remain opaque and are displayed only with the capability's metadata namespace.
The browser never converts branch index, vector order, `Z`, or provider metadata
into universal liquid/vapor/aqueous morphology.

## gRPC-Web outcome and error mapping

`SolvePtFlashResponse` uses a Protobuf `oneof`, making result and service error
mutually exclusive:

| gRPC status / payload | Browser meaning | Phases shown |
| --- | --- | ---: |
| OK + result / `ACCEPTED` | Backend published under its declared contract | 1..N |
| OK + result / `PHASE_SET_UNSTABLE` | Scientific candidate rejection | 0 |
| OK + result / `INDETERMINATE` | Scientific decision unavailable; provenance and diagnostics preserved | 0 |
| OK + `error` | `PtService` request, dispatch, contract, or execution error | 0 |
| non-OK gRPC status | Transport/process/auth/deadline/cancellation failure | 0 |
| OK + malformed/unknown v1 payload | Client wire-contract failure | 0 |

In particular, `INDETERMINATE` must not become a non-OK gRPC status or a service
error. Conversely, an empty or malformed result must not be presented as
indeterminate. The client performs no automatic retry; one UI submission produces
one RPC and the form permits only one in-flight solve. Discovery defaults to a
10-second deadline, solve to 120 seconds, and both accept `AbortSignal` cancellation.
A deadline/cancellation is a transport failure even if server-side work may have
started.

## Native adapter policy

`mpmc::runtime_grpc::PtGrpcServiceAdapter` implements the mapping obligations:

1. bound serialized request bytes before handler deserialization and again after
   decode, plus bound serialized responses and gRPC memory;
2. reject rather than queue work above the configured solve-concurrency limit;
3. require and cap RPC deadlines, observe cancellation before entry/publication,
   and suppress a response when the post-solve check observes either condition;
4. validate required field presence, v1 convention strings, and checked integer
   conversions;
5. forward the unchanged ID-keyed request to `PtService` and map its response
   envelope without recomputing scientific acceptance;
6. return `PtServiceError` as the Protobuf error arm with gRPC OK;
7. reserve non-OK gRPC statuses for process/transport/auth/resource/deadline or
   cancellation failures outside the completed `PtServiceResponse`;
8. expose ordinary native gRPC and leave gRPC-Web/CORS to the audited Envoy edge.

The adapter cannot preempt an already-entered backend because the frozen generic
backend contract has no cancellation token. Such work retains its solve permit
until it returns and its late response is discarded. A process host remains
responsible for constructing real registered backends, credentials, TLS,
authentication/authorization, deployed CORS origins, address binding, lifecycle,
and operational observability.

The browser additionally rejects more than 64 discovered backends, more than 256
components per backend, more than 64 phases/phase-count declarations, excessive
settings/transition items, invalid or over-limit IDs, absent required fields,
non-finite numbers, inconsistent inventory/provenance snapshots, unknown enum
values and illegal phase publication. These checks protect the client and
presentation contract; they are not a second thermodynamic acceptance calculation.

## Generation and verification

Generation is local and deterministic with exact npm versions:

```bash
cd frontend
npm ci --ignore-scripts --no-audit --no-fund
npm run proto:format-check
npm run proto:lint
npm run proto:generate
git diff --exit-code -- src/gen
```

Do not hand-edit `frontend/src/gen`. The checked-in output allows the frontend to
build without a network-time code-generation step while CI verifies that it still
matches the canonical schema.
