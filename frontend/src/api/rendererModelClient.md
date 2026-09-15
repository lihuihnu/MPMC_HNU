# Typed renderer model client

`rendererModelClient()` returns the single cached client for this window's
v2 `mpmcModelDesktopV2` preload bridge (v1 `mpmcModelDesktop` when v2 is absent), or null when no bridge exists. The adapter is
browser-only: it imports no Electron, Node HTTP/2, child-process or host credential
code. UI forms remain separate. Existing PT v1 clients/results are unchanged.

```ts
const client = rendererModelClient();
if (!client) throw new Error('Desktop model bridge is unavailable.');
const { model, snapshot } = await client.create({
  definition: explicitDefinition,
  solverSelection: { case: 'presetId', value: 'mpmc-balanced-default/v1' },
});
const computation = await client.solve(model, {
  pressurePa, temperatureK, feed, // Immutable model component order; SI units.
}, { signal });
if (computation.outcome === 'accepted') {
  // The complete result includes the authoritative flags and actual candidate phases.
}
await client.release(model);
```

Inputs use the canonical generated model initializer types, with private routing
fields omitted. Outputs retain the entire generated `ModelSnapshot` and
`FullPtResult`, not a reduced UI projection. Counts/settings/provider branches
remain bigint, optional presence stays explicit, and special diagnostic doubles
roundtrip via Protobuf JSON. No parameter defaults, inferred kij, normalization,
phase naming, EOS arithmetic or solver-status inference is added.

| Computation outcome | Meaning at this adapter |
| --- | --- |
| accepted | The backend accepted its phase set; global-stability/morphology flags remain independent. |
| phase_set_unstable | Return the full computation and any diagnostic candidates; do not publish it as accepted. |
| indeterminate | Return the full computation, diagnostics and optional candidates; no transport failure is invented. |

Every result retains capability, versions, evaluated PT/feed, maximum phase count,
candidate fractions/compositions/fugacity/branch data, transition evidence,
provider convention/metadata and validity flags. The wrapper adds only a typed
outcome. Unknown outcomes, missing protocol fields, incompatible conventions,
changed immutable snapshots and mismatched PT/component identity fail as data
loss. Accepted candidates require structural finite values; non-accepted
candidates retain non-finite diagnostics. This is protocol validation, not an
independent proof of thermodynamic acceptance. Feed identity uses the same
`2 * 64 * Number.EPSILON * max(1, |actual|, |requested|)` bound as existing
`ptWire.ts`, allowing native roundoff normalization without changing any solver
threshold or normalizing in the renderer.

Call failures throw `RendererModelError` with `code`, `category`, `reason` and
`source` (`client`, `ipc`, `transport`, `contract`). All 16 known gRPC codes have
separate categories; malformed error envelopes/codes become data_loss. Invoke
rejections become unavailable, local AbortSignal cancellation becomes cancelled,
and local deadlines become deadline_exceeded. Recognized fixed IPC/client reason
identifiers are retained; unknown reasons become `ipc.failed`. Exception text is
never forwarded. Desktop v2 additionally carries optional versioned `validation` details; v1 remains
status-only. See [validation errors](modelValidationDetail.md) for safe field paths,
version compatibility and fallback behavior. Missing field precision is never invented.

Returned model references are frozen opaque objects, not transferable strings.
A private map retains the main-process local token and an independent snapshot
clone, so caller edits cannot alter the expected snapshot. Foreign/released refs
fail before dispatch. `reconnect()` synchronously invalidates all old refs,
cancels current calls, drains actual IPC promises and opens a new session;
simultaneous reconnects share that operation. It never recreates models.

There are at most four model reservations and four unresolved IPC invocations.
Cancellation rejects promptly but retains admission until the actual IPC promise
settles. Local deadlines default to 125 seconds (10-second main handshake plus
110-second unary bound and transport headroom); callers may shorten them. A
reconnect's drain and opening each have that finite bound. Abort listeners/timers
are removed on settlement; late replies cannot publish references/results.
Repeated reconnect attempts cannot accumulate waits behind a broken preload.

Ambiguous mutation, transport and protocol errors set `requiresReconnect` and
reject further use until explicit reconnect. Cancellation follows IPC's
conservative whole-window session policy. A cancel arriving after main has
already completed a mutation may be too late to reclaim it immediately; explicit
reconnect or window teardown reclaims the main registry. No automatic Create
retry or exactly-once claim is made. If a broken preload never settles an invoke,
reconnect reports a drain timeout and subsequent attempts fail while it is still
draining; close the window for authoritative main-process cleanup.

The typed client has window lifetime and does not expose a fictitious disconnect
or disposal RPC: IPC v1 owns the native lease in main and reclaims it on reload,
crash, close or app shutdown. Callers can release models explicitly and cancel
active work. UI component unmount must not be described as window/session closure.

Verification uses synthetic schema-only fixtures for all outcomes/error statuses,
bigint/non-finite preservation, malformed replies, snapshot/input ownership,
quotas, cancel/deadline/drain and stale refs. The existing official three-platform
Electron test additionally runs a browser bundle of this actual adapter inside a
sandboxed renderer and compares complete snapshots/results against the native
attributed fixture, exercises native validation failures/reconnect, and observes
independent window cleanup. See [desktop IPC](../../desktop/modelDesktopIpc.md).

Reference: [Protocol Buffers JSON format](https://protobuf.dev/programming-guides/json/).
