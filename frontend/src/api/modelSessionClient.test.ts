import { create } from '@bufbuild/protobuf';
import { Code, ConnectError, type CallOptions } from '@connectrpc/connect';
import { describe, expect, it, vi } from 'vitest';
import {
  CreateModelResponseSchema, DescribeModelResponseSchema, FullPtResultSchema,
  ModelSessionOpenedSchema, ModelSnapshotSchema, ReleaseModelResponseSchema, SolveModelResponseSchema,
  type CreateModelResponse, type SolveModelResponse,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { MODEL_SESSION_CONTRACT, MODEL_SESSION_HEADER, MODEL_WIRE_CONTRACT,
  ModelSessionClient, type ModelReference, type ModelSessionConnection } from './modelSessionClient';

function deferred<T>() {
  let resolve!: (value: T) => void;
  let reject!: (error: unknown) => void;
  const promise = new Promise<T>((yes, no) => { resolve = yes; reject = no; });
  return { promise, resolve, reject };
}
// Transport-only fixtures: no numerical or physical validation claim.
class FakeConnection implements ModelSessionConnection {
  headers = { authorization: 'Bearer private-test-token' };
  closed = false;
  close = vi.fn(() => { this.closed = true; });
  opens = 0;
  serial = 0;
  handshake: Promise<void> = Promise.resolve();
  contract = MODEL_SESSION_CONTRACT;
  end = deferred<void>();
  lastOptions: CallOptions | undefined;
  constructor(readonly id: string) { void this.end.promise.catch(() => {}); }
  sessions = {
    openModelSession: (_request: unknown, options?: CallOptions) => this.stream(options),
  };
  async *stream(options?: CallOptions) {
    ++this.opens;
    await this.handshake;
    if (options?.signal?.aborted) return;
    yield create(ModelSessionOpenedSchema, { wireContract: this.contract, sessionId: this.id });
    if (options?.signal?.aborted) return;
    const abort = () => this.end.resolve();
    options?.signal?.addEventListener('abort', abort, { once: true });
    try { await this.end.promise; }
    finally { options?.signal?.removeEventListener('abort', abort); }
  }
  models = {
    createModel: vi.fn(async (_request: unknown, options?: CallOptions) => {
      this.lastOptions = options;
      return create(CreateModelResponseSchema, { wireContract: MODEL_WIRE_CONTRACT,
        modelHandle: `${this.id}-handle-${++this.serial}`, snapshot: create(ModelSnapshotSchema) });
    }),
    describeModel: vi.fn(async (_request: unknown, options?: CallOptions) => {
      this.lastOptions = options;
      return create(DescribeModelResponseSchema, { wireContract: MODEL_WIRE_CONTRACT, snapshot: create(ModelSnapshotSchema) });
    }),
    solveModel: vi.fn(async (_request: unknown, options?: CallOptions) => {
      this.lastOptions = options;
      return create(SolveModelResponseSchema, { wireContract: MODEL_WIRE_CONTRACT,
        result: create(FullPtResultSchema, { diagnostic: 'transport-fixture' }) });
    }),
    releaseModel: vi.fn(async (_request: unknown, options?: CallOptions) => {
      this.lastOptions = options;
      return create(ReleaseModelResponseSchema, { wireContract: MODEL_WIRE_CONTRACT });
    }),
  };
}
function setup() {
  const connections: FakeConnection[] = [];
  const connector = vi.fn(async () => {
    const connection = new FakeConnection(`session-${connections.length}`);
    connections.push(connection); return connection;
  });
  return { connections, connector, client: new ModelSessionClient(connector) };
}

describe('shared model session ownership', () => {
  it('shares one opening and binds every call to private session metadata', async () => {
    const { client, connector, connections } = setup();
    try {
      const opening = client.connect(); expect(client.connect()).toBe(opening);
      await opening;
      const made = await client.create({ solverSelection: { case: 'presetId', value: 'explicit-test-preset' } });
      expect(connector).toHaveBeenCalledTimes(1);
      expect(Object.keys(made.model)).toEqual([]);
      expect(JSON.stringify(made)).not.toContain('handle-');
      const c = connections[0]!;
      expect(new Headers(c.lastOptions?.headers).get(MODEL_SESSION_HEADER)).toBe(c.id);
      expect(await client.solve(made.model, { feed: [1], pressurePa: 1, temperatureK: 1 }))
        .toMatchObject({ diagnostic: 'transport-fixture' });
      expect(c.models.createModel.mock.calls[0]?.[0]).toMatchObject({ wireContract: MODEL_WIRE_CONTRACT,
        solverSelection: { case: 'presetId', value: 'explicit-test-preset' } });
    } finally { await client.dispose(); }
  });

  it('rejects foreign references and restores model capacity after release', async () => {
    const { client, connections } = setup();
    try {
      await client.connect();
      const made = await Promise.all(Array.from({ length: 4 }, () => client.create({})));
      await expect(client.create({})).rejects.toMatchObject({ reason: 'client.model_limit' });
      await expect(client.describe({} as ModelReference)).rejects.toMatchObject({ reason: 'model.stale_reference' });
      const releasing = client.release(made[0]!.model);
      await expect(client.solve(made[0]!.model, {})).rejects.toMatchObject({ reason: 'model.stale_reference' });
      await releasing;
      await expect(client.release(made[0]!.model)).rejects.toMatchObject({ code: Code.NotFound });
      await client.create({});
      expect(connections[0]!.models.createModel).toHaveBeenCalledTimes(5);
    } finally { await client.dispose(); }
  });

  it('invalidates old references on reconnect and shares simultaneous reconnects', async () => {
    const { client, connections } = setup();
    try {
      await client.connect(); const old = await client.create({});
      const reconnecting = client.reconnect(); expect(client.reconnect()).toBe(reconnecting);
      await reconnecting;
      await expect(client.describe(old.model)).rejects.toMatchObject({ reason: 'model.stale_reference' });
      await client.create({});
      expect(connections).toHaveLength(2);
      expect(connections[0]!.closed).toBe(true);
      expect(connections[1]!.models.describeModel).not.toHaveBeenCalled();
    } finally { await client.dispose(); }
  });

  it('drains pending work before reconnect and discards late solve results', async () => {
    const { client, connections } = setup();
    await client.connect(); const old = await client.create({});
    const late = deferred<SolveModelResponse>();
    connections[0]!.models.solveModel.mockImplementationOnce(async () => late.promise);
    const solving = client.solve(old.model, {});
    const rejected = expect(solving).rejects.toMatchObject({ code: Code.Canceled });
    const reconnecting = client.reconnect();
    expect(connections).toHaveLength(1);
    late.resolve(create(SolveModelResponseSchema, { wireContract: MODEL_WIRE_CONTRACT, result: create(FullPtResultSchema) }));
    await rejected; await reconnecting;
    expect(connections).toHaveLength(2);
    await client.dispose();
  });

  it('closes a late connector result after disposal and never opens its stream', async () => {
    const late = deferred<ModelSessionConnection>(); const c = new FakeConnection('late');
    const client = new ModelSessionClient(async () => late.promise);
    const opening = client.connect(); const rejected = expect(opening).rejects.toMatchObject({ code: Code.Canceled });
    const disposed = client.dispose(); late.resolve(c);
    await rejected; await disposed; await client.dispose();
    expect(c.closed).toBe(true); expect(c.opens).toBe(0);
    await expect(client.connect()).rejects.toMatchObject({ reason: 'client.disposed' });
  });

  it('cancels a pending handshake and an explicitly queued reconnect', async () => {
    const c = new FakeConnection('opening'); const gate = deferred<void>(); c.handshake = gate.promise;
    const client = new ModelSessionClient(async () => c);
    const opening = client.connect(); const rejected = expect(opening).rejects.toMatchObject({ code: Code.Canceled });
    await Promise.resolve();
    const restart = client.reconnect(); const cancelled = expect(restart).rejects.toMatchObject({ code: Code.Canceled });
    const disconnected = client.disconnect(); gate.resolve();
    await rejected; await disconnected; await cancelled;
    expect(client.connected).toBe(false); await client.dispose();
  });

  it('fails closed on malformed handshakes and unexpected stream termination', async () => {
    const c = new FakeConnection('invalid'); c.contract = 'wrong';
    const client = new ModelSessionClient(async () => c);
    await expect(client.connect()).rejects.toMatchObject({ code: Code.DataLoss });
    await client.dispose(); expect(c.closed).toBe(true);
    const normal = setup(); await normal.client.connect(); await normal.client.create({});
    normal.connections[0]!.end.resolve();
    await vi.waitFor(() => expect(normal.client.connected).toBe(false));
    expect(normal.connector).toHaveBeenCalledTimes(1);
    await normal.client.dispose();
  });

  it('never retries ambiguous creation and sanitizes transport errors', async () => {
    const { client, connections } = setup(); await client.connect();
    const c = connections[0]!;
    c.models.createModel.mockRejectedValueOnce(new ConnectError('Bearer private-test-token mh1_secret', Code.Unavailable));
    await expect(client.create({})).rejects.toMatchObject({ code: Code.Unavailable, message: 'Model client: rpc.failed' });
    expect(c.models.createModel).toHaveBeenCalledTimes(1); expect(client.connected).toBe(false);
    await client.dispose(); expect(c.closed).toBe(true);
  });

  it('bounds active calls and aborts them on disposal without publishing late data', async () => {
    const { client, connections } = setup(); await client.connect(); const made = await client.create({});
    const late = deferred<SolveModelResponse>();
    connections[0]!.models.solveModel.mockImplementation(async (_request, options) => {
      connections[0]!.lastOptions = options; return late.promise;
    });
    const calls = Array.from({ length: 4 }, () => client.solve(made.model, {}));
    const outcomes = Promise.allSettled(calls);
    await expect(client.solve(made.model, {})).rejects.toMatchObject({ reason: 'client.request_limit' });
    const disposed = client.dispose();
    expect(connections[0]!.lastOptions?.signal?.aborted).toBe(true);
    late.resolve(create(SolveModelResponseSchema, { wireContract: MODEL_WIRE_CONTRACT, result: create(FullPtResultSchema) }));
    expect((await outcomes).every(outcome => outcome.status === 'rejected')).toBe(true);
    await disposed;
  });

  it('closes the session after a lost release response', async () => {
    const { client, connections } = setup(); await client.connect(); const made = await client.create({});
    connections[0]!.models.releaseModel.mockRejectedValueOnce(new ConnectError('lost', Code.Unavailable));
    await expect(client.release(made.model)).rejects.toMatchObject({ code: Code.Unavailable });
    expect(client.connected).toBe(false); await client.dispose();
  });

  it('discards a late created handle after disconnect and drains its reservation', async () => {
    const { client, connections } = setup(); await client.connect();
    const late = deferred<CreateModelResponse>();
    connections[0]!.models.createModel.mockImplementationOnce(async () => late.promise);
    const creating = client.create({}); const rejected = expect(creating).rejects.toMatchObject({ code: Code.Canceled });
    const disconnected = client.disconnect();
    late.resolve(create(CreateModelResponseSchema, { wireContract: MODEL_WIRE_CONTRACT,
      modelHandle: 'late-secret-handle', snapshot: create(ModelSnapshotSchema) }));
    await rejected; await disconnected; await client.connect();
    await Promise.all(Array.from({ length: 4 }, () => client.create({})));
    await client.dispose();
  });
});
