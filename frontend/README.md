# MPMC_HNU PT frontend

This React + TypeScript + Vite application consumes the versioned
[`mpmc.runtime.v1.PtFlashService`](../api/README.md) contract. A hosted browser
uses gRPC-Web; product shells reuse the same React product components through
platform-specific typed ownership adapters. The generic configured-PT path first
discovers configured backends and their exact component inventories, then builds
the solve form from the selected immutable capability snapshot.

The [Android Product Shell](../products/pt_android/README.md) reuses the same
Classic PR workspace and configured-PT compatibility surface through app-local
Capacitor/JNI adapters. Its packaging and native lifecycle contract remains in
that product directory.

## Classic PR product workspace

Classic PR is the shared editable PR76 product workflow. Electron, Android and an
explicitly authorized hosted-Web integration render the same React shell and
`ExpertPr76Workspace`; they differ only in the typed adapter that owns the model.
Electron and Android are local no-login products. Hosted Web instead requires an
authoritative browser identity provider and an authenticated model-session route:

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
change returned values. PR76 EOS evaluation, stability search, phase splitting,
final phase-set review, parameter validation and model lifetime remain native C++
backend responsibilities.

The shared CSS and React presentation are responsive so supported products retain
one visual language. Transport parity is platform-specific: Electron uses its
context-isolated model-workbench/IPC path; Android uses the app-local
Capacitor/JNI typed ownership adapter backed by the same native
`Pr76ExecutableModel`; hosted Web uses the existing same-origin authenticated
`webModelSession` connector, `ModelSessionClient` and `bindExpertModelOwner`.
Hosted ownership is fail-closed: without the exact trusted runtime capability,
the browser remains on its existing configured-PT product and does not open a
model session.

## Scientific and service boundary

The frontend is not an EOS or flash implementation. It does not:

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

The configured-PT client sends binary gRPC-Web requests to:

```text
<base-url>/mpmc.runtime.v1.PtFlashService/DiscoverPtCapabilities
<base-url>/mpmc.runtime.v1.PtFlashService/SolvePtFlash
```

The endpoint or proxy must provide the required browser CORS response. The
repository includes an audited [development/CI Envoy edge](../deploy/pt-grpc-web/README.md),
a thin [native C++ adapter](../modules/runtime_grpc/README.md), and a
repository-curated three-backend process host, but no deployed endpoint. If the
variable is absent outside a local product shell, the generic hosted frontend
remains explicitly unconfigured and sends no request.

Hosted Classic PR is deliberately **not** enabled by a build-time bearer or a
`VITE_*` token. A trusted hosting shell may provide this runtime capability only
after its deployment has authoritative end-user identity mapping and the required
same-origin model routes:

```ts
window.mpmcHostedWebModelOwnership = {
  convention: 'MPMC/model/hosted-web-ownership/v1',
  baseUrl: '/model-api',
  identity: authoritativeIdentityProvider,
};
```

`identity` is a `WebIdentityProvider` function object, not a token value. The
same-origin connector obtains the bearer only when the first model is applied,
opens one `ModelSessionClient` lease, and keeps that identity snapshot for the
session epoch. The token is not put in DOM state, browser storage, generated HTML
or user-facing errors. Page teardown disposes the session. There is no automatic
Create retry, reconnect or principal change.

The checked-in production edge currently exposes only `PtFlashService`; it does
**not** yet expose model-session/configuration routes or map a Web access token to
an authoritative native principal. Therefore production deployment must not
publish `mpmcHostedWebModelOwnership` merely because the frontend adapter exists.
That trusted-edge identity/routing gate remains a separate security increment.
See [`webModelSession.md`](src/api/webModelSession.md) and the
[edge deployment contract](../deploy/pt-grpc-web/README.md).

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
is retried automatically. The shared Classic PR owner similarly invalidates stale
model generations rather than publishing late results; the hosted owner opens its
authenticated session before dispatching the first Create.

## Product adapters

Electron embeds the built React files without giving the renderer Node access.
Its sandboxed preload exposes the configured-PT bridge plus the typed local
model-workbench capabilities used by Classic PR. Main-process IPC enforces sender
identity, bounded request shapes, exact request IDs, deadlines, cancellation and
lifecycle cleanup; native process/session transport remains outside React.

Android uses the same `DesktopProductShell`/`ExpertPr76Workspace` component tree
(the component name is historical, not a platform restriction). Its
`ModelWorkbenchBridge` implementation is app-local Capacitor/JNI transport. Java
only serializes work onto the native executor and transports bounded fields;
`Pr76ExecutableModel` remains authoritative for model preparation, settings,
solve and result semantics. No gRPC C++ runtime is added to the APK.

Hosted Web uses the same product component tree only when the trusted runtime
capability described above is present and valid. `hostedWebModelOwner.ts` lazily
connects the existing authenticated `ModelSessionClient` and adapts it with
`bindExpertModelOwner`; it does not expose handles, reproduce model lifecycle or
interpret scientific results. Missing or malformed capability falls back to the
existing configured-PT hosted application without attempting anonymous model
access.

The configured compatibility path remains available and continues to exercise
repository-curated PR76, SW92 and CPA backends. Editable Classic PR remains
PR76-only; this frontend work does not change those backend implementations or
add editable SW92/CPA semantics.

The resulting desktop directories and Android APK are engineering products with
separate release/signing gates. Windows downstream packaging is documented in the
[product index](../products/pt/README.md) and
[desktop installer README](../products/pt/desktop_installer/README.md). Android
packaging/lifecycle evidence is documented in the
[Android README](../products/pt_android/README.md).

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
PT wire contract as well as Classic PR model editing, P/T/z input, local/hosted
product presentation, phase-count/fraction/composition output, chart values,
non-accepted candidate isolation, hosted runtime-capability gating and structured
error handling.

An additional cross-language golden starts a synthetic C++ `PtService`, maps it
through the native gRPC adapter and Envoy, and exercises discovery, accepted,
indeterminate, service-error, provenance, variable-phase, and CORS behavior with
this generated TypeScript client. It has software-contract meaning only and is
not a physical regression.

Browser dependencies remain confined to the frontend/product shells and do not
become prerequisites of `runtime`, `flash` or `thermodynamics`.
