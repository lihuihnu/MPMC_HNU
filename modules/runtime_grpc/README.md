# PT native gRPC process adapter v1

## Boundary

`mpmc::runtime_grpc::PtGrpcServiceAdapter` is the optional native gRPC mapping
for the versioned [`mpmc.runtime.v1`](../../api/README.md) protocol. It depends on
the transport-neutral [`mpmc::runtime::PtService`](../runtime/README.md), generated
Protobuf messages, and gRPC C++. Neither `runtime`, `flash`, nor
`thermodynamics` depends on this module.

For each admitted `SolvePtFlash` RPC, the adapter:

1. checks process policy and required wire presence/version;
2. copies the unchanged ID-keyed request into `PtServiceRequest`;
3. calls `PtService::solve()` exactly once;
4. maps the returned result-or-error envelope without changing its decision.

It does not register models, select EOS roots, invoke a model-specific solver,
normalize composition, retry, infer morphology, or rerun material balance,
fugacity, TPD, common-tangent, phase-count, or phase-transition logic. Backend
construction and registration remain host-application choices. This module is a
library adapter, not a production worker executable.

## Process policy

`PtGrpcAdapterLimits` has explicit, validated defaults:

| Control | Default | Enforcement |
| --- | ---: | --- |
| serialized request | 64 KiB | gRPC pre-deserialization limit and adapter post-decode check |
| serialized response | 4 MiB | adapter pre-publication check and gRPC send limit |
| concurrent solves | 1 | non-blocking admission gate; excess calls receive `RESOURCE_EXHAUSTED` |
| gRPC resource quota | 64 MiB | server-builder gRPC memory quota |
| gRPC worker threads | 16 | resource-quota thread cap for synchronous handlers |
| discovery deadline | 10 s | finite client deadline required and capped |
| solve deadline | 120 s | finite client deadline required and capped |

The default single-solve gate is compatible with backends that have not declared
thread safety. A host may raise it only after every concurrently reachable
backend and referenced evaluator/cache has been audited for concurrent use. The
adapter queues no rejected solve and performs no automatic retry.

Request size is bounded by `ServerBuilder::SetMaxReceiveMessageSize()` before a
Protobuf message reaches the handler. `ByteSizeLong()` is checked again so direct
handler use cannot bypass policy. `configure_pt_grpc_server()` also applies the
response limit and `grpc::ResourceQuota` memory/thread bounds; it intentionally leaves address,
credentials, TLS, authentication, and process lifetime to the host.

## Deadline and cancellation semantics

Both RPCs require a finite client deadline. Admission is rejected if the
remaining deadline exceeds the method policy, has already elapsed, or the call
is already cancelled. State is checked again immediately before entering
`PtService` and immediately before response mapping.

The frozen `PtFlashBackend::solve(p,T,z)` interface has no cancellation token.
Consequently, deadline/cancellation can prevent admission and suppress a late
response, but cannot safely preempt a solve that has already entered a backend.
That solve keeps its concurrency permit until `PtService::solve()` returns. An
elapsed deadline, or cancellation observed at the post-solve checkpoint,
discards its result at the RPC boundary. The client and gRPC runtime still make
their own deadline/cancellation decision if cancellation races the final check.
This avoids both deliberately publishing stale work and starting overlapping
work against a backend assumed to be serialized.

`indeterminate` remains an OK gRPC response on the scientific result arm with
zero published phases. `PtServiceError` remains an OK gRPC response on the
service-error arm. Only wire, resource, process, deadline, cancellation, or
unexpected adapter failures use non-OK gRPC status.

## Build and dependencies

The module generates native C++ bindings at build time from the canonical schema;
generated C++ is not committed:

```bash
cmake -S tests/runtime/pt_grpc_adapter -B build/pt-grpc-adapter
cmake --build build/pt-grpc-adapter --parallel 2
ctest --test-dir build/pt-grpc-adapter --verbose
```

The focused official-runner workflow pins Ubuntu 24.04 packages for Protobuf
3.21.12 and gRPC 1.51.1, and validates the edge with the exact Envoy 1.36.9
container tag. These are isolated process-edge dependencies; the numerical core
keeps its existing dependency set. No upstream source was copied. Protocol
Buffers uses its upstream BSD-style license; gRPC and Envoy use Apache-2.0.

## Browser edge and verification

Browsers connect through the audited development/CI
[`Envoy gRPC-Web edge`](../../deploy/pt-grpc-web/README.md). CORS and gRPC-Web
framing stay at that edge rather than being reimplemented in C++.

Native regressions cover exact discovery/result mapping, service-error and
scientific-indeterminate separation, request/message limits, concurrency
admission, required/capped deadlines, explicit cancellation, permit draining,
and public-header containment. The cross-language golden starts a synthetic C++
backend behind the real adapter and Envoy, then exercises it with the generated
TypeScript gRPC-Web client, including CORS preflight. The fixture is explicitly
model-neutral and contains no EOS or physical reference data; existing backend
suites remain authoritative for all scientific gates.
