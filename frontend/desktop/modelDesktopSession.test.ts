import { create, type JsonObject } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';
import { expect, it, vi } from 'vitest';
import { MODEL_DESKTOP_CONVENTION, MODEL_DESKTOP_V2_CONVENTION, type ModelDesktopOperation, type ModelDesktopReply } from '../src/api/modelDesktopContract';
import { ModelClientError, type ModelReference } from '../src/api/modelSessionClient';
import { FullPtResultSchema, ModelSnapshotSchema } from '../src/gen/mpmc/model_configuration/v1/model_service_pb';
import { MODEL_VALIDATION_DETAIL_VERSION } from '../src/api/modelValidationDetail';
import { ModelDesktopSession } from './modelDesktopSession';

function deferred<T>() {
  let resolve!: (value: T) => void;
  const promise = new Promise<T>(yes => { resolve = yes; });
  return { promise, resolve };
}
// Policy-only stand-in: synthetic schema messages, no numerical validation.
function setup() {
  const snapshot = create(ModelSnapshotSchema, { settings: { version: 'test-version',
    eosRoot: { maxIterations: 9007199254740993n } } });
  const result = create(FullPtResultSchema, { diagnostic: 'policy-fixture',
    candidatePhaseSet: { phases: [{ compressibilityFactor: Number.NaN }] } });
  const client = {
    connected: false,
    connect: vi.fn(async () => { client.connected = true; }),
    disconnect: vi.fn(async () => { client.connected = false; }),
    dispose: vi.fn(async () => { client.connected = false; }),
    create: vi.fn(async () => ({ model: Object.freeze({}) as ModelReference, snapshot })),
    describe: vi.fn(async () => snapshot),
    solve: vi.fn(async () => result),
    release: vi.fn(async () => {}),
  };
  return { client, owner: new ModelDesktopSession(client) };
}
function request(operation: ModelDesktopOperation, requestId = 'test', extra = {}) {
  return { version: MODEL_DESKTOP_CONVENTION, operation, requestId, ...extra };
}
function value(reply: ModelDesktopReply): JsonObject {
  expect(reply.ok).toBe(true);
  if (!reply.ok) throw new Error('Expected successful policy reply.');
  return reply.value as JsonObject;
}
async function make(owner: ModelDesktopSession) {
  return value(await owner.invoke(request('create', 'create', { input: { presetId: 'explicit-test' } }))).model as string;
}

it('preserves canonical full snapshots/results and keeps references local to one window', async () => {
  const a = setup(); const b = setup();
  try {
    const token = await make(a.owner);
    for (const op of ['describe', 'solve', 'release'] as const) {
      expect(await b.owner.invoke(request(op, op, { model: token, ...(op === 'solve' ? { input: {} } : {}) })))
        .toMatchObject({ ok: false, error: { code: Code.NotFound, reason: 'ipc.stale_reference' } });
    }
    expect(b.client.describe).not.toHaveBeenCalled(); expect(b.client.solve).not.toHaveBeenCalled();
    expect(b.client.release).not.toHaveBeenCalled();
    expect(value(await a.owner.invoke(request('describe', 'd', { model: token }))))
      .toMatchObject({ settings: { eosRoot: { maxIterations: '9007199254740993' } } });
    expect(value(await a.owner.invoke(request('solve', 's', { model: token, input: { feed: [1] } }))))
      .toMatchObject({ candidatePhaseSet: { phases: [{ compressibilityFactor: 'NaN' }] } });
    expect(await a.owner.invoke(request('release', 'r', { model: token }))).toMatchObject({ ok: true, value: null });
    expect(await a.owner.invoke(request('describe', 'd2', { model: token }))).toMatchObject({ ok: false, error: { code: Code.NotFound } });
    expect(a.client.release).toHaveBeenCalledTimes(1);
  } finally { await a.owner.dispose(); await b.owner.dispose(); }
});

it('rejects invalid versions, fields, operations, shapes and oversized requests before dispatch', async () => {
  const { owner, client } = setup();
  try {
    for (const raw of [null, request('connect', ''), { ...request('connect'), version: 'v9' },
      { ...request('connect'), operation: [] }, request('connect', 'x', { input: {} }),
      request('create', 'x', { input: { wireContract: 'injected' } }),
      request('create', 'x', { input: { definition: { invented: true } } }),
      request('solve', 'x', { model: 'raw-server-handle', input: {} }),
      request('connect', 'x', { ownerWindow: 1 })]) {
      expect(await owner.invoke(raw)).toMatchObject({ ok: false, error: { code: Code.InvalidArgument } });
    }
    expect(await owner.invoke(request('create', 'x', { input: { presetId: 'x'.repeat(65536) } })))
      .toMatchObject({ ok: false, error: { code: Code.ResourceExhausted } });
    expect(client.connect).not.toHaveBeenCalled(); expect(client.create).not.toHaveBeenCalled();
  } finally { await owner.dispose(); }
});

it('reconnect invalidates old tokens without affecting another window', async () => {
  const a = setup(); const b = setup();
  try {
    const old = await make(a.owner); const other = await make(b.owner);
    expect(await a.owner.invoke(request('reconnect'))).toMatchObject({ ok: true });
    expect(await a.owner.invoke(request('describe', 'old', { model: old }))).toMatchObject({ error: { code: Code.NotFound } });
    expect(await b.owner.invoke(request('describe', 'other', { model: other }))).toMatchObject({ ok: true });
    expect(await make(a.owner)).not.toBe(old);
    expect(b.client.disconnect).not.toHaveBeenCalled();
  } finally { await a.owner.dispose(); await b.owner.dispose(); }
});

it('discards a late create after reset and waits for the prior client drain before reopening', async () => {
  const { owner, client } = setup();
  const late = deferred<Awaited<ReturnType<typeof client.create>>>();
  client.create.mockImplementationOnce(() => late.promise);
  const creating = owner.invoke(request('create', 'pending', { input: {} }));
  await vi.waitFor(() => expect(client.create).toHaveBeenCalledTimes(1));
  const drain = deferred<void>();
  client.disconnect.mockImplementationOnce(() => { client.connected = false; return drain.promise; });
  const resetting = owner.reset();
  const opening = owner.invoke(request('connect', 'next'));
  expect(client.connect).toHaveBeenCalledTimes(1);
  late.resolve({ model: {} as ModelReference, snapshot: create(ModelSnapshotSchema) });
  expect(await creating).toMatchObject({ ok: false, error: { code: Code.Canceled } });
  drain.resolve(); await resetting;
  expect(await opening).toMatchObject({ ok: true });
  expect(client.connect).toHaveBeenCalledTimes(2);
  await owner.dispose();
});

it('bounds pending calls, rejects duplicate IDs and scopes cancellation to one owner', async () => {
  const a = setup(); const b = setup();
  const token = await make(a.owner); const other = await make(b.owner);
  const late = deferred<Awaited<ReturnType<typeof a.client.solve>>>();
  a.client.solve.mockImplementation(() => late.promise);
  const calls = Array.from({ length: 4 }, (_, i) => a.owner.invoke(request('solve', `s${i}`, { model: token, input: {} })));
  await vi.waitFor(() => expect(a.client.solve).toHaveBeenCalledTimes(4));
  expect(await a.owner.invoke(request('connect', 's0'))).toMatchObject({ error: { code: Code.AlreadyExists } });
  expect(await a.owner.invoke(request('connect', 'fifth'))).toMatchObject({ error: { code: Code.ResourceExhausted } });
  b.owner.cancel('s0'); a.owner.cancel('missing');
  expect(a.client.disconnect).not.toHaveBeenCalled();
  a.owner.cancel('s0');
  expect(a.client.disconnect).toHaveBeenCalledTimes(1);
  expect(await b.owner.invoke(request('describe', 'b', { model: other }))).toMatchObject({ ok: true });
  late.resolve(create(FullPtResultSchema));
  for (const reply of await Promise.all(calls)) expect(reply).toMatchObject({ ok: false, error: { code: Code.Canceled } });
  expect(await a.owner.invoke(request('describe', 'old', { model: token }))).toMatchObject({ error: { code: Code.NotFound } });
  await a.owner.dispose(); await b.owner.dispose();
});

it('window disposal rejects late solves and cannot be reopened', async () => {
  const { owner, client } = setup(); const token = await make(owner);
  const late = deferred<Awaited<ReturnType<typeof client.solve>>>();
  client.solve.mockImplementationOnce(() => late.promise);
  const pending = owner.invoke(request('solve', 'pending', { model: token, input: {} }));
  await vi.waitFor(() => expect(client.solve).toHaveBeenCalledTimes(1));
  await owner.dispose(); await owner.dispose();
  late.resolve(create(FullPtResultSchema));
  expect(await pending).toMatchObject({ ok: false, error: { code: Code.Canceled } });
  expect(await owner.invoke(request('reconnect'))).toMatchObject({ ok: false, error: { reason: 'ipc.window_closed' } });
});

it('sanitizes unexpected errors and closes ambiguous creation state', async () => {
  const { owner, client } = setup();
  client.create.mockRejectedValueOnce(new Error('private-bearer-session-handle'));
  const reply = await owner.invoke(request('create', 'bad', { input: {} }));
  expect(reply).toMatchObject({ ok: false, error: { code: Code.Internal, reason: 'ipc.failed' } });
  expect(JSON.stringify(reply)).not.toContain('private-');
  expect(client.disconnect).toHaveBeenCalledTimes(1);
  await owner.dispose();
});

it('closes the window session when a full reply exceeds its boundary', async () => {
  const { owner, client } = setup(); const token = await make(owner);
  client.solve.mockResolvedValueOnce(create(FullPtResultSchema, { diagnostic: 'x'.repeat(4 * 1024 * 1024) }));
  expect(await owner.invoke(request('solve', 'large', { model: token, input: {} })))
    .toMatchObject({ ok: false, error: { reason: 'ipc.reply_size' } });
  expect(client.disconnect).toHaveBeenCalledTimes(1);
  await owner.dispose();
});

it('shares v1/v2 ownership and emits validation details exclusively on v2', async () => {
  const { owner, client } = setup(); const model = await make(owner);
  const validation = { version: MODEL_VALIDATION_DETAIL_VERSION, code: 'configuration.invalid_value', field: 'parameters.pure[nitrogen].critical_temperature_k' };
  try {
    const raw = request('solve', 'v1', { model, input: {} });
    client.solve.mockRejectedValueOnce(new ModelClientError(Code.InvalidArgument, 'rpc.failed', validation));
    expect(await owner.invoke(raw)).toEqual({ version: MODEL_DESKTOP_CONVENTION, ok: false,
      error: { code: Code.InvalidArgument, reason: 'rpc.failed' } });
    client.solve.mockRejectedValueOnce(new ModelClientError(Code.InvalidArgument, 'rpc.failed', validation));
    expect(await owner.invoke({ ...raw, version: MODEL_DESKTOP_V2_CONVENTION }, MODEL_DESKTOP_V2_CONVENTION))
      .toEqual({ version: MODEL_DESKTOP_V2_CONVENTION, ok: false, error: { code: Code.InvalidArgument, reason: 'rpc.failed', validation } });
    expect(await owner.invoke({ ...raw, version: MODEL_DESKTOP_V2_CONVENTION })).toMatchObject({ version: MODEL_DESKTOP_CONVENTION, error: { reason: 'ipc.unsupported_version' } });
    expect(await owner.invoke(raw, MODEL_DESKTOP_V2_CONVENTION)).toMatchObject({ version: MODEL_DESKTOP_V2_CONVENTION, error: { reason: 'ipc.unsupported_version' } });
    expect(client.connect).toHaveBeenCalledTimes(1);
    const v2 = { ...request('release', 'release', { model }), version: MODEL_DESKTOP_V2_CONVENTION };
    expect(await owner.invoke(v2, MODEL_DESKTOP_V2_CONVENTION)).toMatchObject({ ok: true });
    expect(await owner.invoke(request('describe', 'old', { model }))).toMatchObject({ error: { reason: 'ipc.stale_reference' } });
  } finally { await owner.dispose(); }
});
it('rechecks optional detail shape at main and shares pending admission across versions', async () => {
  const { owner, client } = setup(); const model = await make(owner);
  try {
    client.solve.mockRejectedValueOnce(new ModelClientError(Code.InvalidArgument, 'rpc.failed', {
      version: MODEL_VALIDATION_DETAIL_VERSION, code: 'configuration.invalid_value', field: 'components[mh1_private]',
    }));
    expect(await owner.invoke({ ...request('solve', 'bad', { model, input: {} }), version: MODEL_DESKTOP_V2_CONVENTION }, MODEL_DESKTOP_V2_CONVENTION))
      .toMatchObject({ error: { validation: { version: MODEL_VALIDATION_DETAIL_VERSION, code: 'configuration.invalid_value' } } });
    const pending = deferred<Awaited<ReturnType<typeof client.solve>>>(); client.solve.mockImplementation(() => pending.promise);
    const calls = Array.from({ length: 4 }, (_, i) => {
      const version = i % 2 ? MODEL_DESKTOP_V2_CONVENTION : MODEL_DESKTOP_CONVENTION;
      return owner.invoke({ ...request('solve', `held-${i}`, { model, input: {} }), version }, version);
    });
    await vi.waitFor(() => expect(client.solve).toHaveBeenCalledTimes(5));
    expect(await owner.invoke({ ...request('connect', 'fifth'), version: MODEL_DESKTOP_V2_CONVENTION }, MODEL_DESKTOP_V2_CONVENTION))
      .toMatchObject({ error: { reason: 'ipc.request_limit' } });
    owner.cancel('held-0'); pending.resolve(create(FullPtResultSchema));
    for (const reply of await Promise.all(calls)) expect(reply).toMatchObject({ error: { code: Code.Canceled } });
  } finally { await owner.dispose(); }
});
