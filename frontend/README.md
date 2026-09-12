# MPMC_HNU model-neutral PT frontend

This React + TypeScript + Vite application consumes the versioned
[`mpmc.runtime.v1.PtFlashService`](../api/README.md) gRPC-Web contract. It first
discovers configured PT backends and their exact component inventories, then
builds the solve form from the selected immutable capability snapshot.

## Scientific and service boundary

The browser is not an EOS or flash implementation. It does not:

- evaluate thermodynamic properties or roots;
- decide phase count, stability, common tangency or acceptance;
- normalize, clip, fill or fit feed composition;
- infer morphology from phase order, branch, `Z`, role or family IDs;
- retry an indeterminate solve or fall back to another backend.

The selected backend determines the component IDs and canonical order. Users edit
only `p`, `T`, and one mole fraction per discovered component. Configured scalar
settings such as a backend-owned salinity value are discovery provenance and are
read-only; they are not model-neutral solve inputs.

Provider role/family metadata is shown verbatim with its declared namespace. It is
not relabeled as liquid, vapor or aqueous by generic frontend code.

## Connection

Set the gRPC-Web base URL when building or serving the application:

```bash
VITE_MPMC_GRPC_WEB_BASE_URL=https://example.test/pt-api npm run dev
```

The client sends binary gRPC-Web requests to:

```text
<base-url>/mpmc.runtime.v1.PtFlashService/DiscoverPtCapabilities
<base-url>/mpmc.runtime.v1.PtFlashService/SolvePtFlash
```

The endpoint or proxy must provide the required browser CORS response. If the
variable is absent, the frontend remains explicitly unconfigured and sends no
request. This repository does not ship a fake production backend or a deployed
C++ worker.

Startup discovery must complete before solve is enabled. Switching the configured
backend switches the inventory; component IDs cannot be added, removed or edited
in the browser. A result provenance snapshot must exactly match the selected
discovery descriptor or the response is rejected as a wire-contract failure.

## Result states

The presentation keeps four boundaries visible:

| State | Source | Presentation |
| --- | --- | --- |
| `accepted` | Protobuf result arm | Variable 1..N accepted phase cards and provenance |
| `phase_set_unstable` | Protobuf result arm | Scientific rejection, zero phase cards |
| `indeterminate` | Protobuf result arm | Scientific uncertainty with diagnostic/provenance, zero phase cards |
| `PtServiceError` | Protobuf error arm | Service error code, field and diagnostic |
| gRPC or wire failure | RPC / client adapter | Connection/deadline/cancel/contract failure |

An indeterminate computation is therefore never displayed as “service unavailable,”
and an RPC/service error is never presented as a thermodynamic decision.

The form allows one in-flight solve. Discovery has a 10-second client deadline,
solve has a 120-second deadline, both accept cancellation, and neither is retried
automatically.

## Development

Requires Node 24.x. Exact direct dependencies and the complete transitive graph are
pinned by `package.json` and `package-lock.json`.

```bash
cd frontend
npm ci --ignore-scripts --no-audit --no-fund
npm run proto:format-check
npm run proto:lint
npm run proto:generate
git diff --exit-code -- src/gen
npm run typecheck
npm test
npm run build
```

The generated `src/gen/mpmc/runtime/v1/pt_service_pb.ts` is committed and must not
be hand-edited. Tests cover binary Protobuf discovery round-trip, request value/order
preservation, capability/inventory integrity, variable accepted phases, distinct
indeterminate/service-error oneof arms, malformed-response rejection, gRPC timeout
options and no automatic retry.

The frontend remains independent of the C++ build. The browser dependencies are
confined here and do not become prerequisites of `runtime`, `flash` or
`thermodynamics`.
