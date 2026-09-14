import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { pathToFileURL } from 'node:url';
import { fromJson, toJson, type JsonObject, type JsonValue } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';
import { modelCleanupComplete } from './modelCleanupProbe';
import { app, session, type BrowserWindow } from 'electron';
import { MODEL_DESKTOP_CONVENTION, type ModelDesktopReply } from '../../src/api/modelDesktopContract';
import { MODEL_SESSION_HEADER, MODEL_WIRE_CONTRACT, ModelSessionClient } from '../../src/api/modelSessionClient';
import { CreateModelRequestSchema, FullPtResultSchema, ModelSnapshotSchema, SolveModelRequestSchema } from '../../src/gen/mpmc/model_configuration/v1/model_service_pb';
import { createDesktopModelConnection } from '../modelConnection';
import { registerModelDesktopIpc } from '../modelDesktopIpc';
import { PtHostSession } from '../hostSession';
import { PtDesktopGateway } from '../ptGateway';
import { createPtDesktopWindow, restrictPtDesktopSession } from '../window';

function requireThat(condition: unknown, message: string): asserts condition {
  if (!condition) throw new Error(message);
}
function deferred() {
  let resolve!: () => void;
  const promise = new Promise<void>(yes => { resolve = yes; });
  return { promise, resolve };
}
async function until(predicate: () => boolean, message: string) {
  const deadline = Date.now() + 10_000;
  while (!predicate()) {
    requireThat(Date.now() < deadline, message);
    await new Promise(yes => setTimeout(yes, 10));
  }
}
async function run() {
  const binary = process.env.MPMC_MODEL_CLIENT_HOST;
  const fixture = process.env.MPMC_MODEL_CLIENT_FIXTURE;
  requireThat(binary && fixture, 'Native model host and attributed fixture are required.');
  const data = JSON.parse(readFileSync(fixture, 'utf8')) as { create: JsonObject; solve: JsonObject };
  const definition = fromJson(CreateModelRequestSchema, data.create);
  const state = fromJson(SolveModelRequestSchema, data.solve);
  // Strip private wire routing fields; all parameter/settings values stay exact.
  const createInput = { ...data.create }; delete createInput.wireContract;
  const solveInput = { ...data.solve }; delete solveInput.wireContract; delete solveInput.modelHandle;
  const logs: string[] = [];
  const replies: ModelDesktopReply[] = [];
  const host = new PtHostSession(binary, line => logs.push(line), { enableModelSessions: true });
  const gateway = new PtDesktopGateway(host);
  const connection = await host.start();
  const observer = createDesktopModelConnection(connection);
  interface Observed { id: string; handle: string; closed: boolean }
  const observed: Observed[] = [];
  let hold: { operation: 'create' | 'solve'; entered: boolean; gate: ReturnType<typeof deferred> } | undefined;
  const entryPath = resolve('dist/index.html');
  const ipc = registerModelDesktopIpc(() => new ModelSessionClient(async signal => {
    const connected = await host.start(); signal.throwIfAborted();
    const transport = createDesktopModelConnection(connected);
    const entry: Observed = { id: '', handle: '', closed: false }; observed.push(entry);
    return {
      ...transport,
      sessions: {
        async *openModelSession(request, options) {
          for await (const opened of transport.sessions.openModelSession(request, options)) {
            entry.id = opened.sessionId!; yield opened;
          }
        },
      },
      models: {
        ...transport.models,
        async createModel(request, options) {
          const made = await transport.models.createModel(request, options); entry.handle = made.modelHandle!;
          const pending = hold;
          if (pending?.operation === 'create') { hold = undefined; pending.entered = true; await pending.gate.promise; }
          return made;
        },
        async solveModel(request, options) {
          const solved = await transport.models.solveModel(request, options);
          const pending = hold;
          if (pending?.operation === 'solve') { hold = undefined; pending.entered = true; await pending.gate.promise; }
          return solved;
        },
      },
      close() { entry.closed = true; transport.close(); },
    };
  }), pathToFileURL(entryPath).href);
  const windows = new Set<BrowserWindow>();
  async function window(attach = true) {
    const made = createPtDesktopWindow({ preloadPath: resolve('desktop/preload.cjs'), showWhenReady: false });
    windows.add(made);
    if (attach) ipc.attach(made.webContents);
    await made.loadFile(entryPath);
    return made;
  }
  let serial = 0;
  async function call(target: BrowserWindow, operation: string, ...args: JsonValue[]): Promise<ModelDesktopReply> {
    const reply = await target.webContents.executeJavaScript(
      `globalThis.mpmcModelDesktop[${JSON.stringify(operation)}](${[JSON.stringify(`ipc-${++serial}`), ...args.map(arg => JSON.stringify(arg))].join(',')})`,
    ) as ModelDesktopReply;
    requireThat(reply.version === MODEL_DESKTOP_CONVENTION, 'IPC reply version changed.'); replies.push(reply);
    return reply;
  }
  function value(reply: ModelDesktopReply): JsonValue {
    requireThat(reply.ok, !reply.ok ? `Real model IPC failed: ${reply.error.code} ${reply.error.reason}` : 'Unexpected reply.'); return reply.value;
  }
  async function make(target: BrowserWindow) {
    const made = value(await call(target, 'create', createInput)) as JsonObject;
    requireThat(typeof made.model === 'string', 'Missing local model reference.'); return made;
  }
  function options(id: string) {
    const headers = new Headers(observer.headers); headers.set(MODEL_SESSION_HEADER, id);
    return { headers, timeoutMs: 5_000 };
  }
  async function gone(entry: Observed, expected: 'session.not_found' | 'registry.model_not_found' = 'session.not_found') {
    const deadline = Date.now() + 10_000;
    while (true) {
      try {
        await observer.models.describeModel({ wireContract: MODEL_WIRE_CONTRACT, modelHandle: entry.handle }, options(entry.id));
      } catch (cause) {
        if (modelCleanupComplete(cause, expected)) return;
      }
      requireThat(Date.now() < deadline, 'Native host retained the closed model/session.');
      await new Promise(yes => setTimeout(yes, 10));
    }
  }
  function rejected(reply: ModelDesktopReply, code: Code) {
    requireThat(!reply.ok && reply.error.code === code, 'IPC rejection code mismatch.');
  }
  function equal(actual: JsonValue, expected: JsonValue) {
    requireThat(JSON.stringify(actual) === JSON.stringify(expected), 'Complete native snapshot/result changed across IPC.');
  }
  try {
    restrictPtDesktopSession(session.defaultSession);
    await gateway.models.connect();
    const baseline = await gateway.models.create(definition);
    const expectedSnapshot = toJson(ModelSnapshotSchema, baseline.snapshot);
    const expectedResult = toJson(FullPtResultSchema, await gateway.models.solve(baseline.model, state));
    const a = await window(); const b = await window(); const rogue = await window(false);
    rejected(await call(rogue, 'connect'), Code.PermissionDenied);
    const first = await make(a); const oldA = { ...observed.at(-1)! };
    const second = await make(b); const entryB = { ...observed.at(-1)! };
    equal(first.snapshot!, expectedSnapshot);
    equal(value(await call(a, 'describe', first.model!)), expectedSnapshot);
    equal(value(await call(a, 'solve', first.model!, solveInput)), expectedResult);
    for (const op of ['describe', 'solve', 'release']) {
      rejected(await call(b, op, first.model!, ...(op === 'solve' ? [solveInput] : [])), Code.NotFound);
    }
    rejected(await call(a, 'create', { ...createInput, wireContract: 'forbidden' }), Code.InvalidArgument);
    // Fill/release/reuse all four model reservations through the real preload.
    const extras = [];
    for (let i = 0; i < 3; ++i) extras.push(await make(a));
    rejected(await call(a, 'create', createInput), Code.ResourceExhausted);
    const released = { ...observed[0]! };
    value(await call(a, 'release', extras.at(-1)!.model!)); await gone(released, 'registry.model_not_found');
    await make(a);
    console.info('MODEL_IPC_ROUNDTRIP_OK');
    value(await call(a, 'reconnect')); await gone(oldA);
    rejected(await call(a, 'describe', first.model!), Code.NotFound);
    equal(value(await call(b, 'solve', second.model!, solveInput)), expectedResult);
    const reloaded = await make(a); const priorDocument = { ...observed.at(-1)! };
    const loaded = new Promise<void>(yes => a.webContents.once('did-finish-load', () => yes()));
    a.webContents.reload(); await loaded; await gone(priorDocument);
    rejected(await call(a, 'describe', reloaded.model!), Code.NotFound);
    await make(a); const beforeClose = { ...observed.at(-1)! };
    const closed = new Promise<void>(yes => a.once('closed', () => yes())); a.close(); await closed;
    await gone(beforeClose);
    equal(value(await call(b, 'solve', second.model!, solveInput)), expectedResult);

    console.info('MODEL_IPC_RECONNECT_RELOAD_CLOSE_OK');

    // Actual renderer loss must revoke the old document even if WebContents survives.
    const crashed = new Promise<void>(yes => b.webContents.once('render-process-gone', () => yes()));
    b.webContents.forcefullyCrashRenderer(); await crashed; await gone(entryB);
    await b.loadFile(entryPath);
    rejected(await call(b, 'describe', second.model!), Code.NotFound);
    const survivor = await make(b);

    // Hold an actual native reply at the transport boundary while its renderer dies.
    for (const operation of ['create', 'solve'] as const) {
      const transient = await window();
      const made = operation === 'solve' ? await make(transient) : undefined;
      const pending = { operation, entered: false, gate: deferred() }; hold = pending;
      void call(transient, operation, ...(made ? [made.model!, solveInput] : [createInput])).catch(() => {});
      try {
        await until(() => pending.entered, 'Native operation did not reach the held reply.');
        const live = observed.at(-1)!; const dying = { ...live };
        transient.destroy();
        await until(() => live.closed, 'Window destruction did not abort its native transport.');
        await gone(dying);
      } finally { pending.gate.resolve(); hold = undefined; }
    }
    console.info('MODEL_IPC_CRASH_INFLIGHT_CLOSE_OK');
    // More windows than the host's 16-session capacity; the surviving B stays usable.
    for (let i = 0; i < 18; ++i) {
      const transient = await window(); await make(transient);
      const dying = { ...observed.at(-1)! }; transient.destroy(); await gone(dying);
    }
    equal(value(await call(b, 'solve', survivor.model!, solveInput)), expectedResult);
    // Run the production typed adapter in the actual sandboxed renderer, not in main.
    const rendererCode = readFileSync(resolve('desktop-model-test-dist/model-renderer-smoke.js'), 'utf8');
    const typedA = await window(); const typedB = await window();
    await typedA.webContents.executeJavaScript(`${rendererCode}\n;void 0;`);
    await typedB.webContents.executeJavaScript(`${rendererCode}\n;void 0;`);
    async function renderer(target: BrowserWindow, operation: string, ...args: JsonValue[]): Promise<JsonObject> {
      return await target.webContents.executeJavaScript(
        `MpmcModelRendererSmoke[${JSON.stringify(operation)}](${args.map(arg => JSON.stringify(arg)).join(',')})`,
      ) as JsonObject;
    }
    const typedInitial = await renderer(typedA, 'initialize', createInput, solveInput);
    equal(typedInitial.snapshot!, expectedSnapshot); equal(typedInitial.described!, expectedSnapshot);
    equal(typedInitial.result!, expectedResult); requireThat(typedInitial.outcome === 'accepted', 'Typed accepted outcome changed.');
    await renderer(typedB, 'initialize', createInput, solveInput); const typedEntryB = { ...observed.at(-1)! };
    const invalid = await renderer(typedA, 'invalidSolve');
    requireThat(invalid.code === Code.InvalidArgument && invalid.category === 'invalid_argument' && invalid.source === 'ipc',
      'Native validation error lost its typed renderer status.');
    equal((await renderer(typedA, 'solve')).result!, expectedResult);
    const releasedTyped = await renderer(typedA, 'releaseAndRecreate');
    requireThat((releasedTyped.released as JsonObject).reason === 'renderer.stale_reference', 'Released renderer reference was reused.');
    equal(releasedTyped.result!, expectedResult);
    const recovery = await renderer(typedA, 'invalidCreateAndReconnect');
    requireThat((recovery.failure as JsonObject).code === Code.InvalidArgument &&
      (recovery.blocked as JsonObject).reason === 'renderer.reconnect_required' &&
      (recovery.stale as JsonObject).reason === 'renderer.stale_reference', 'Renderer mutation recovery/status mapping changed.');
    equal(recovery.result!, expectedResult);
    const typedEntryA = { ...observed.at(-1)! }; typedA.destroy(); await gone(typedEntryA);
    equal((await renderer(typedB, 'solve')).result!, expectedResult);
    typedB.destroy(); await gone(typedEntryB);
    console.info('MODEL_RENDERER_CLIENT_OK typed snapshots/full results, native validation errors, release, reconnect invalidation and independent window cleanup');
    const last = [...observed].reverse().find(entry => !entry.closed)!;
    await ipc.dispose(); await ipc.dispose(); await gone(last);
    requireThat(observed.every(entry => entry.closed), 'An IPC transport survived disposal.');
    const output = JSON.stringify(replies) + logs.join('\n');
    requireThat(!output.includes(connection.bearerToken), 'Private host credential leaked.');
    for (const entry of observed) {
      requireThat(!output.includes(entry.id) && !output.includes(entry.handle), 'Private native routing data leaked.');
    }
    console.info('MODEL_DESKTOP_IPC_OK real preload roundtrip, window isolation, reconnect, reload, crash, in-flight close, 18 close/reopen cycles and server reclamation');
  } finally {
    for (const target of windows) if (!target.isDestroyed()) target.destroy();
    await ipc.dispose(); observer.close(); await gateway.stop();
  }
  requireThat(logs.some(line => line.includes('"event":"pt_process_stopped"')), 'Native host did not shut down.');
}
// Test lifetime extends beyond the final window through asynchronous host teardown.
app.on('window-all-closed', () => {});
const watchdog = setTimeout(() => { console.error('Model IPC smoke timed out.'); app.exit(1); }, 180_000);
void app.whenReady().then(async () => {
  try { await run(); clearTimeout(watchdog); console.info('MODEL_DESKTOP_IPC_COMPLETE'); app.exit(0); }
  catch (cause) { console.error(cause instanceof Error ? cause.message : 'Model IPC smoke failed.'); app.exit(1); }
});
