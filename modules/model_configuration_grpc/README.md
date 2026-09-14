# Versioned model configuration service

This optional C++ gRPC service connects the bounded PR76 registry to real unary
RPCs. Its package is `mpmc.model_configuration.v1`; every request/response uses
`mpmc.model_configuration.v1/model-service/v1`. It does not reinterpret or change
`mpmc.runtime.v1.PtFlashService`, configured backend IDs or existing clients.

```cmake
add_subdirectory(modules/model_configuration_grpc)
target_link_libraries(my_host PRIVATE mpmc::model_configuration_grpc)
```

| RPC | Request | Result |
| --- | --- | --- |
| `CreateModel` | Complete definition + explicit preset ID or complete versioned settings | Opaque handle + resolved model/settings/capability/host-policy snapshot |
| `DescribeModel` | Handle | Owning snapshot of the immutable model |
| `SolveModel` | Handle, P [Pa], T [K], ordered mole feed | Complete native result, including diagnostic candidate phases |
| `ReleaseModel` | Handle | Invalidate future use; already admitted registry work retains ownership |

Creation never infers a preset. Unknown versions/families/presets and incomplete
parameters/settings reject without fallback. PR76 custom models are supported;
SW92/CPA custom creation remains unsupported. Existing curated SW/CPA solves
remain available through the original service. All 57 numerical/boolean settings
have explicit protobuf presence. The factory still owns domain and quota checks;
zero/false cannot turn into an implicit native default.

The canonical new schema is bundled under this module's `proto/` import root,
with the existing `api/proto` root supplying frozen v1 capability/transition data
types. Code generation uses protoc/gRPC, never edited generated code. Keeping
this additive schema with its optional service avoids generating unused frontend
clients while UI integration is postponed. Later clients must compile this same
schema rather than invent another representation.

## Host composition and lifetime

The host owns a `Pr76ModelRegistry` and a `ModelGrpcServiceAdapter` constructed
with a **required** thread-safe authorization callback. The callback must authorize
the caller for that registry/session using the host's authenticated identity or
session policy; neither a model handle nor this module authenticates creation.
It is not supplied by an RPC. There is no unauthenticated production default.
A callback returning false produces `PERMISSION_DENIED`; unexpected callback
exceptions produce a sanitized `INTERNAL` error.

The host selects listening address and credentials, then calls
`configure_model_grpc_server(builder, adapter)` before `BuildAndStart`. Keep the
adapter, registry and callback dependencies alive until server `Shutdown` and
`Wait` finish. Close session registries on host session teardown. Configure a
separate listener or compatible shared message/quota policies when cohosting;
message size and gRPC resource quota are builder-wide. The coexistence test uses
identical 64 KiB receive limits for both services, so old pre-parse limits are not
weakened. This increment does not change production host flags, Envoy routes,
Electron/Android IPC, user/session routing or UI.

Default adapter limits: 1 MiB request, 4 MiB response, 64 MiB gRPC resource quota,
4 concurrent mapped requests, 1 concurrent solve, 120 s maximum client deadline.
The host can tighten them. Configure applies receive bounds before protobuf
request deserialization and caps outbound messages. The handler checks serialized
size again, rejects unknown fields recursively and checks model array shape before
DTO copies. Client-supplied host ceilings/hints are not ignored as unknown fields.
Unknown enum values reject. Native allocations remain subject to registry and
per-model limits; the gRPC quota is not a universal process-heap limit.

Every RPC requires a finite deadline within host policy. Lifecycle is checked
before work and before publication. Cancellation/deadline does not interrupt the
current numerical kernel; its lease remains held until the existing solve returns.
The adapter suppresses a late result and preserves native scientific outcomes.
Request/solve admission is released through RAII, including exceptions. No queue,
retry, hidden parameter rebuild or global solver workspace is added.

If creation succeeds locally but encoding, response size, cancellation or deadline
prevents publication, the adapter releases its newly created handle. If the server
returns OK and the network subsequently loses that response, creation has an
ambiguous outcome: there is no exactly-once/idempotency protocol in this increment.
Do not automatically retry Create; bounded session teardown reclaims unknown
handles. Release is committed even if its response is lost/cancelled; another
release returns NOT_FOUND. This is the registry's existing strict release contract.

## Results and errors

`FullPtResult` carries capability, all state/outcome fields, optional candidate
phase set (beta, compositions, ln(phi), provider branch/smoothness, optional Z),
transition evidence, provenance conventions, diagnostic, global-stability limit,
provider metadata and morphology flag. **Candidate presence is not acceptance.**
Only the outcome establishes acceptance. The original v1 `PtComputationResult`
continues to expose accepted phases only. An indeterminate numerical solve is a
successful RPC with an indeterminate scientific outcome, not an INTERNAL error.

Adapter-generated non-OK statuses contain a serialized `ModelServiceError` in
`grpc::Status.error_details`: exact wire version, stable domain-qualified code and
bounded field context. Internal `what()`, request bodies and bearer handles are
not echoed. Transport failures before dispatch (including receive limits and
client-observed cancellation/deadline) can lack these details.

| gRPC status | Representative detail codes |
| --- | --- |
| INVALID_ARGUMENT | `wire.missing_field`, `wire.unsupported_version`, `wire.unknown_field`, `configuration.missing_parameter`, `configuration.invalid_settings`, `registry.invalid_handle`, `request.rejected` |
| UNIMPLEMENTED | `configuration.unsupported_version`, `configuration.unsupported_family`, `configuration.unsupported_preset` |
| NOT_FOUND | `registry.model_not_found` (unknown, foreign, released; duplicate release) |
| RESOURCE_EXHAUSTED | `configuration.resource_limit`, `registry.capacity_exceeded`, `model.busy`, `rpc.concurrency_limit`, `wire.request_limit`, `wire.response_limit`, `rpc.memory_exhausted` |
| FAILED_PRECONDITION | `registry.closed` |
| UNAVAILABLE | `registry.entropy_unavailable` |
| PERMISSION_DENIED | `rpc.not_authorized` |
| CANCELLED / DEADLINE_EXCEEDED | `rpc.cancelled` / `rpc.deadline_exceeded` when emitted by the adapter |
| INTERNAL | Sanitized `rpc.internal_failure` or `registry.internal_failure` |

All existing `ModelConfigurationErrorCode` values have explicit code mappings.
Finite/deadline-policy errors use `rpc.deadline_required` / `rpc.deadline_limit`
with INVALID_ARGUMENT. Codes are service semantics, not exception-message parsing.

## Audit and verification scope

Baseline: PR #82 `da8859471c6eb9f1465de41fa4629e98a6a6aacf`. Before writing, reviewed
root AGENTS, the Gate, registry reservations/leases, factory snapshots, existing
v1 protobuf/result projection, gRPC admission/deadline/authentication/size policy,
process ownership and focused dependency CI. Extending v1 IDs/results would break
its frozen meanings. A separate optional service reuses existing native validation
and RPC dependencies, with its own full-candidate result and no UI changes.

Consulted protobuf [field-presence guidance](https://protobuf.dev/programming-guides/field_presence/)
and gRPC [deadline behavior](https://grpc.io/docs/guides/deadlines/): explicit presence
preserves zero-valued settings; application work must handle lifecycle itself.
No new external parameter data, numerical formula or tolerance is introduced.

`tests/model_configuration/grpc` has 16 real-loopback RPC cases and runs 6 unchanged
v1 adapter cases. Tests cover create/describe/solve/release, all 57 missing numeric
fields, preset/custom identity, pseudo/provenance/applicability roundtrips,
full-result single/binary/unseeded-ternary parity, errors/quotas, busy and retained
release, pre-dispatch receive limits, response rollback, authorization, internal
failure redaction, deadline/cancellation cleanup, and legacy cohosting. A test-only
wire decoder compares results against directly constructed native PR76 models;
expected results do not reuse production result encoding. Fixtures remain the
repository's attributed binary and explicitly synthetic structural examples.

The focused official-hosted workflow runs GCC ASan/UBSan, Clang and MSVC. Linux
reuses existing pinned apt packages; Windows restores the existing locked Conan
dependency graph without upgrading it. No product packaging or unrelated UI tests
are selected. Scientific kernels, existing v1 schema/adapter and configured
backends stay unchanged. TLS/Envoy deployment, cross-language clients, UI, public
hints, one-sided applicability and full product integration remain later gates.
