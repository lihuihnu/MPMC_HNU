# MPMC_HNU model-neutral PT frontend

This React + TypeScript + Vite application consumes the versioned
[`mpmc.runtime.v1.PtFlashService`](../api/README.md) contract. A hosted browser
uses gRPC-Web; the Electron desktop shell reuses the same React renderer through
a versioned context-isolated preload bridge and native gRPC. Both first discover
configured PT backends and their exact component inventories, then build the
solve form from the selected immutable capability snapshot.

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

The endpoint or proxy must provide the required browser CORS response. The
repository includes an audited [development/CI Envoy edge](../deploy/pt-grpc-web/README.md),
a thin [native C++ adapter](../modules/runtime_grpc/README.md), and a
repository-curated three-backend process host, but no deployed endpoint. If the
variable is absent outside Electron, the frontend remains explicitly
unconfigured and sends no request.

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

## Electron desktop vertical slice v1

Electron embeds the built React files without giving the renderer Node access.
The sandboxed preload exposes only discovery, solve, and cancel under
`MPMC/PT/desktop-bridge/v1`. Main-process IPC enforces sender identity, 64 KiB
request shape, at most four active bridge calls, exact request IDs, deadlines,
and cancellation. The main process alone starts the staged
`mpmc_pt_service_host`, connects using native HTTP/2 gRPC, and maps the existing
Protobuf contract. It contains no EOS, flash, parameter, retry, fallback, or
scientific acceptance logic.

Each launch creates a 256-bit random bearer, sends it to the child only through
stdin, and accepts readiness only from an ephemeral `127.0.0.1` port reporting
exactly three configured backends. Parent shutdown or stdin EOF gracefully stops
the child. Production mTLS and deployed gRPC-Web remain separate modes; desktop
does not need Envoy, CORS, a deployment domain, or production certificates for
its same-device child session.

The hosted product gate loads the actual React page in Electron, waits for its
three-backend discovery state, and then traverses preload/IPC/native gRPC to ask
the real repository-curated PR76, SW92, and CPA backends for one solve each. The
smoke requires a scientific result arm and exact backend provenance but does not
require `accepted`; it is integration evidence, not a replacement for physical
regression or a license to reinterpret `indeterminate`.

The resulting Linux x64, Windows x64, and macOS arm64 directories are unsigned
engineering previews with `release_eligible=false`. They are not one-click public
installers yet: code signing, notarization, licensing, update policy, and
publisher identity remain fail-closed release gates.

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

An additional cross-language golden starts a synthetic C++ `PtService`, maps it
through the native gRPC adapter and Envoy, and exercises discovery, accepted,
indeterminate, service-error, provenance, variable-phase, and CORS behavior with
this generated TypeScript client. It has software-contract meaning only and is
not a physical regression.

The frontend remains independent of the C++ build. The browser dependencies are
confined here and do not become prerequisites of `runtime`, `flash` or
`thermodynamics`.
