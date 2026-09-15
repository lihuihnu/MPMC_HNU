# Shared model session client

`ModelSessionClient` owns one server-streaming `OpenModelSession` lease and its
model references. It accepts a connection factory, so ownership is shared between
clients without importing Node, Electron, transport credentials or UI code into
the manager. The desktop gateway exposes `gateway.models` in the **main process**;
the additive model preload IPC now creates a separate owner for each window.
The host must opt into model sessions. The actual app does so on its authenticated private loopback child.

```ts
await gateway.models.connect();
const { model, snapshot } = await gateway.models.create(explicitDefinitionAndSolver);
const result = await gateway.models.solve(model, { pressurePa, temperatureK, feed });
await gateway.models.release(model);
await gateway.stop();
```

Inputs/results use generated types from the canonical C++ service proto and frozen
v1 capability types. No EOS, parameter, normalization or result-outcome algorithm
is duplicated. Creation still requires an explicit preset or complete settings;
an empty/invalid selection reaches the authoritative server validator. `solve`
returns the complete `FullPtResult`, including unresolved outcomes and candidate
phases. Candidate presence alone never implies acceptance.

The opaque frozen `ModelReference` contains no server handle. A private bounded
map binds object identity to a single live session. Foreign, released and previous
session references fail locally before RPC dispatch. Successful creation returns
an owning snapshot and a reference; snapshots contain no bearer capabilities.
Limits are four resident model references plus creation reservations, and four
active unary calls with no queue. Host quotas remain authoritative. Defaults:
10 s connection/handshake timeout, 29 min lease (below the host's 30 min maximum),
110 s unary timeout; per-call finite timeouts may only tighten that bound.

- `connect()` shares an existing opening. It does not automatically reconnect a
  lost session. A closing session must drain before a new opening is accepted.
- `reconnect()` invalidates current references synchronously, aborts old work,
  drains it, and opens a new lease. Simultaneous reconnects share one operation.
  A concurrent disconnect/dispose cancels that queued opening.
- `release()` invalidates the reference before awaiting the server reply. Its
  server-side capacity is reusable after a successful reply. Any release failure
  closes the session to reclaim a potentially unreleased model.
- `disconnect()` clears references, aborts the lease and all unary calls, closes
  the transport, and waits for the stream and admitted calls to finish. It is
  reusable; `dispose()` adds permanent closure and is idempotent.
- Stream loss and ambiguous RPC errors invalidate the session; every failed
  dispatched create also closes it conservatively. There is **no Create retry**
  or invisible model rebuild. Late results/handshakes cannot restore old refs.

The connection factory is trusted main-process plumbing. It must honor its abort
signal/settle promptly, and its idempotent `close()` must not throw. The manager
also closes any late-acquired transport. The desktop adapter owns a dedicated
`Http2SessionManager` and aborts it to release sockets, with 64 KiB write / 4 MiB
read limits. Authentication headers never enter client errors, returned model
objects, logs or renderer IPC. Errors carry a sanitized `Code` and local reason;
rich field-level domain-error presentation remains deferred.
The service uses standard `google.rpc.Status` rich errors with a typed
`ModelServiceError` detail, so native Connect error codes remain interoperable.

Client close completion means local work/transport drained, not a server cleanup
acknowledgment. The server observes stream cancellation/transport failure and
reclaims its registry; its finite lease bounds undetected network loss. Real
host regressions observe disappearance independently and reconnect more than
the host's 16-session capacity, so a leaked resident slot cannot hide behind a
locally cleared map. No Web endpoint is exposed before Web identity policy.

`npm run proto:model:generate` uses existing pinned Buf/Protobuf-ES packages and
`buf.model.yaml`'s two local import roots. Additive generated ES code is an ignored
build artifact, regenerated before typecheck and tests. Existing committed v1
code remains under its original regeneration/diff gate. No runtime download or
new package dependency is introduced.

Verification: deterministic fake-transport tests exercise ownership, abort,
late responses, admission and error handling. The official service workflow runs
real Node -> native host tests on Linux GCC/ASan, Windows and macOS, exporting the
existing attributed binary definition/PT state from the C++ fixture. These are
client/transport regressions, not independent physical validation. UI, Android
IPC and Web routes remain subsequent work.

Desktop IPC now owns one shared client per explicitly attached window. See
[desktop IPC v1](../../desktop/modelDesktopIpc.md) for the additive preload API,
local references, frame authorization, lifecycle and hosted regression scope.
