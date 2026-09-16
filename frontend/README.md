# MPMC_HNU PT frontend

This React + TypeScript + Vite application consumes the versioned
[`mpmc.runtime.v1.PtFlashService`](../api/README.md) contract. A hosted browser
uses gRPC-Web; the Electron desktop shell reuses the same React renderer through
a versioned context-isolated preload bridge and native gRPC. The generic PT path
first discovers configured backends and their exact component inventories, then
builds the solve form from the selected immutable capability snapshot.

The [Android Product Shell](../products/pt_android/README.md) also reuses the
React PT UI through its app-local Capacitor/JNI `FlashClient`. Its build and
lifecycle contract is maintained in that product directory.

## Classic PR product workspace

Electron exposes the existing editable PR76 model path as the default **Classic
PR** workspace. It is local and requires no account or login. The user-facing
flow is deliberately narrower than the underlying service contracts:

1. define the ordered component set and each component's molar mass, critical
   temperature, critical pressure, and acentric factor;
2. provide every required binary interaction coefficient `kij` and, when needed,
   select the already-supported explicit PT solver settings;
3. apply the immutable PR76 model snapshot;
4. enter pressure, temperature, and the initial overall mole fractions;
5. run the native PR76 PT flash and view the accepted phase count, mole phase
   fractions, and component mole fractions in every accepted phase.

Accepted phase fractions and phase compositions are presented with responsive
bar charts plus an exact numerical table. Phase labels remain role-neutral
(`Phase 1`, `Phase 2`, ...): the frontend does not infer liquid, vapor, aqueous,
or other morphology from phase order, density, branch, or compressibility.
Native diagnostics, stability declarations, complete result JSON, and model
provenance remain available only under advanced disclosure sections.

The charts are presentation only. They do not normalize, clip, fit, or otherwise
change returned values. The PR76 EOS, stability search, phase splitting, final
phase-set review, model registry, and model-session lifetime remain native C++
backend responsibilities.

The shared CSS and React presentation are responsive so supported shells can
retain one visual language. Capability parity is still transport-specific: the
current Android shell and ordinary hosted browser expose the existing
model-neutral configured-PT path, while editable PR76 model ownership is provided
by the local Electron model-workbench bridge. This frontend PR does not invent a
new Web/Android model-configuration transport or authentication policy.

## Scientific and service boundary

The browser is not an EOS or flash implementation. It does not:

- evaluate thermodynamic properties or roots;
- decide phase count, stability, common tangency or acceptance;
- normalize, clip, fill or fit feed composition;
- infer morphology from phase order, branch, `Z`, role or family IDs;
- retry an indeterminate solve or fall back to another backend.

On the generic configured-PT path, the selected backend determines the component
IDs and canonical order. Users edit only `p`, `T`, and one mole fraction per
discovered component. Configured scalar settings such as a backend-owned salinity
value are discovery provenance and are read-only; they are not model-neutral
solve inputs.

Provider role/family metadata is shown verbatim with its declared namespace. It is
not relabeled as liquid, vapor or aqueous by generic frontend code.

## Connection

Set the gRPC-Web base URL when building or serving the generic hosted application:

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
variable is absent outside Electron, the generic frontend remains explicitly
unconfigured and sends no request.

Startup discovery must complete before generic PT solve is enabled. Switching the
configured backend switches the inventory; component IDs cannot be added, removed
or edited through that generic path. A result provenance snapshot must exactly
match the selected discovery descriptor or the response is rejected as a
wire-contract failure.

## Result states

The presentation keeps four boundaries visible:

| State | Source | Presentation |
| --- | --- | --- |
| `accepted` | Protobuf result arm | Accepted phase count/fractions/compositions and provenance |
| `phase_set_unstable` | Protobuf result arm | Scientific rejection, no accepted phase presentation |
| `indeterminate` | Protobuf result arm | Scientific uncertainty with diagnostic/provenance, no accepted phases |
| `PtServiceError` | Protobuf error arm | Service error code, field and diagnostic |
| gRPC or wire failure | RPC / client adapter | Connection/deadline/cancel/contract failure |

An indeterminate computation is therefore never displayed as “service unavailable,”
and an RPC/service error is never presented as a thermodynamic decision.

The generic form allows one in-flight solve. Discovery has a 10-second client
deadline, solve has a 120-second deadline, both accept cancellation, and neither
is retried automatically.

## Electron desktop vertical slice v1

Electron embeds the built React files without giving the renderer Node access.
The sandboxed preload exposes the model-neutral PT bridge plus the typed local
model-workbench capabilities used by the Classic PR workspace. Main-process IPC
enforces sender identity, bounded request shapes, exact request IDs, deadlines,
and lifecycle cleanup. The main process alone owns native process/session
transport. Renderer code contains no EOS, flash, parameter fitting, retry,
fallback, or scientific acceptance logic.

Each launch creates a 256-bit random bearer, sends it to the child only through
stdin, and accepts readiness only from an ephemeral `127.0.0.1` port. Parent
shutdown or stdin EOF gracefully stops the child. Production mTLS and deployed
gRPC-Web remain separate modes; desktop does not need Envoy, CORS, a deployment
domain, or production certificates for its same-device child session.

The existing model-neutral hosted product gate traverses preload/IPC/native gRPC
to exercise the repository-curated PR76, SW92, and CPA configured backends. The
editable Classic PR workspace is a separate PR76-only product surface layered on
the existing model-configuration/session client. This PR changes presentation and
entry-point behavior only; it does not alter those backend implementations.

The resulting Linux x64, Windows x64, and macOS arm64 directories are unsigned
engineering previews with `release_eligible=false`. Downstream packaging is
documented in the [product index](../products/pt/README.md): Windows has a
[desktop MSI and protected signed RC gate](../products/pt/desktop_installer/README.md).
Public-release signing/notarization, licensing, update and publisher requirements
remain governed by those product contracts.

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
be hand-edited. Model-configuration TypeScript is generated from the dedicated
versioned model-service proto before typecheck/test/build. Tests cover the generic
PT wire contract as well as Classic PR model editing, P/T/z input, no-login
product presentation, phase-count/fraction/composition output, chart values,
non-accepted candidate isolation, and structured error handling.

An additional cross-language golden starts a synthetic C++ `PtService`, maps it
through the native gRPC adapter and Envoy, and exercises discovery, accepted,
indeterminate, service-error, provenance, variable-phase, and CORS behavior with
this generated TypeScript client. It has software-contract meaning only and is
not a physical regression.

The frontend remains independent of the C++ build. Browser dependencies are
confined here and do not become prerequisites of `runtime`, `flash` or
`thermodynamics`.
