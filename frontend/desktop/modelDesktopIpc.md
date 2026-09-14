# Desktop model IPC v1

The additive sandboxed preload API is `window.mpmcModelDesktop`, convention
`MPMC/model/desktop-bridge/v1`. The existing `mpmcPtDesktop` v1 bridge and UI are
unchanged. Main registers fixed invoke/cancel channels once, explicitly attaches
trusted WebContents, and drains this registration before stopping the shared host.
Every attached window owns an independent `ModelSessionClient` and transport.

| Method | Arguments after unique request ID | Successful value |
| --- | --- | --- |
| connect | none | null |
| reconnect | none | null; all old local references invalid |
| create | `{definition, presetId}` or `{definition, settings}` | `{model, snapshot}` |
| describe | local model reference | complete ModelSnapshot JSON |
| solve | local reference, `{pressurePa, temperatureK, feed}` | complete FullPtResult JSON |
| release | local reference | null |
| cancel | none | no reply; closes this window session if ID is active |

Every invoke reply has `version` and either `{ok:true,value}` or
`{ok:false,error:{code,reason}}`. Codes use the existing gRPC numeric convention;
reasons are fixed IPC/shared-client identifiers. Exception text, host bearer,
server session IDs and raw model handles are never returned. JSON payloads use
the canonical generated Protobuf mapping: enum names, 64-bit integer strings,
and standard special floating-point strings. Do not convert these through an
unversioned UI DTO. Candidate phases are diagnostic unless the authoritative
outcome accepts them; no numerical result/status reinterpretation occurs here.

Main validates the exact attached sender and its current main frame at the exact
packaged entry URL (fragment allowed, query/other paths rejected). Navigation to
a different document blocks calls immediately and resets ownership; only a
trusted main-frame commit enables calls again. Prevented/failed navigations stay
fail-closed until the trusted page is loaded again. Subframes, detached senders
and unattached windows cannot invoke or cancel. Renderer-supplied window IDs,
wire handles or authentication metadata never select an owner.

Local references are random UUIDs mapped privately to in-memory shared-client
references. They are scoped to the window's document/session, not transferable
across windows or reloads. Create opens a session lazily. Explicit reconnect,
main-document reload, renderer crash, window destruction and app shutdown
invalidate references synchronously, abort calls and drain owned native streams.
Late create/solve replies cannot publish success after invalidation. Closing one
window does not invalidate another. Cancellation conservatively invalidates all
models and calls of that window, including ambiguous mutations, with no automatic
Create retry/recreation. Unknown cancellation IDs are harmless.

Bounds: 8 attached windows including closing/draining owners; 4 active requests
per window, retaining reservations until old work returns; 4 models including
pending creates via the shared client; 64 KiB request JSON and 4 MiB reply JSON.
Duplicate IDs and unknown version/operation/fields are rejected before dispatch.
The transport independently enforces protobuf message bounds and finite deadlines.
Electron structured cloning precedes main-process validation, so these are
application admission bounds, not a claim of pre-allocation limits inside Electron.
Dispose waits for local drain; server reclamation is asynchronous.

Verification: policy tests cover frame authorization, ownership, schemas, quota,
cancel and late replies; official Linux/Windows/macOS jobs additionally run the
actual sandboxed Electron preload against the native authenticated host. They
check full native snapshot/result equality, cross-window denial, capacity reuse,
independent reconnect, reload, crash, close with held actual native replies and
18 close/reopen cycles. An independent native client observes `session.not_found`
after teardown, and `registry.model_not_found` after release. Fixtures come from
the existing attributed C++ exporter; no new scientific validation is claimed.

The registration factory is an ownership seam, not a renderer-controlled API.
`PtDesktopGateway.createModelSession()` returns a caller-owned client; IPC disposes
all such clients before `gateway.stop()`. The legacy gateway `models` property
retains its original gateway-owned shutdown behavior.

The [typed renderer client](../src/api/rendererModelClient.md) now consumes this
bridge and preserves complete snapshots/results with explicit error states.
UI forms, Web routes and Android sessions remain future work. Sources: [Electron security](https://www.electronjs.org/docs/latest/tutorial/security)
and [WebContents lifecycle](https://www.electronjs.org/docs/latest/api/web-contents).


Desktop v2 adds optional versioned validation details while retaining the exact v1
preload/reply contract. Both fixed channel pairs share one window owner and its
resource/lifecycle policy. See [validation detail contract](../src/api/modelValidationDetail.md)
for negotiation, field path semantics and privacy checks.
