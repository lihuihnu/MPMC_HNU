import { create, type JsonObject } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';
import { describe, expect, it, vi } from 'vitest';

import {
  MODEL_WORKBENCH_CONVENTION,
  type ModelWorkbenchOperation,
  type ModelWorkbenchReply,
} from '../src/api/modelWorkbenchContract';
import { ModelClientError, type ModelReference, type ModelSessionClient } from '../src/api/modelSessionClient';
import {
  FullPtResultSchema,
  ModelSnapshotSchema,
} from '../src/gen/mpmc/model_configuration/v1/model_service_pb';
import { ModelWorkbenchSession } from './modelWorkbenchSession';

function setup() {
  const snapshot = create(ModelSnapshotSchema, {
    definition: { datasetId: 'private-workbench-test', revision: 'r1' },
  });
  const result = create(FullPtResultSchema, { diagnostic: 'workbench-policy-fixture' });
  let serial = 0;
  const client = {
    connected: false,
    connect: vi.fn(async () => { client.connected = true; }),
    disconnect: vi.fn(async () => { client.connected = false; }),
    dispose: vi.fn(async () => { client.connected = false; }),
    create: vi.fn(async (_input?: unknown) => ({
      model: Object.freeze({ serial: ++serial }) as unknown as ModelReference,
      snapshot,
    })),
    describe: vi.fn(async (_model?: ModelReference) => snapshot),
    solve: vi.fn(async (_model?: ModelReference, _input?: unknown) => result),
    release: vi.fn(async (_model?: ModelReference) => {}),
  };
  return {
    client,
    owner: new ModelWorkbenchSession(client as unknown as ModelSessionClient),
  };
}

function request(operation: ModelWorkbenchOperation, requestId = 'test', input?: JsonObject) {
  return {
    version: MODEL_WORKBENCH_CONVENTION,
    operation,
    requestId,
    ...(input === undefined ? {} : { input }),
  };
}

function value(reply: ModelWorkbenchReply) {
  expect(reply.ok).toBe(true);
  if (!reply.ok) throw new Error('Expected successful workbench reply.');
  return reply.value;
}

describe('main-owned model workbench session', () => {
  it('returns snapshots/results without exposing a model reference', async () => {
    const { owner, client } = setup();
    try {
      const applied = value(await owner.invoke(request('apply', 'a', { presetId: 'test' })));
      expect(applied).toMatchObject({ definition: { datasetId: 'private-workbench-test' } });
      expect(JSON.stringify(applied)).not.toContain('modelHandle');
      expect(JSON.stringify(applied)).not.toContain('model-token');

      const solved = value(await owner.invoke(request('solve', 's', {
        pressurePa: 1_000_000,
        temperatureK: 300,
        feed: [1],
      })));
      expect(solved).toMatchObject({ diagnostic: 'workbench-policy-fixture' });
      expect(client.solve).toHaveBeenCalledTimes(1);
      expect(await owner.invoke(request('release', 'r'))).toMatchObject({ ok: true, value: null });
      expect(client.release).toHaveBeenCalledTimes(1);
    } finally {
      await owner.dispose();
    }
  });

  it('atomically retires the prior model when a replacement is applied', async () => {
    const { owner, client } = setup();
    try {
      expect(await owner.invoke(request('apply', 'a1', {}))).toMatchObject({ ok: true });
      expect(await owner.invoke(request('apply', 'a2', {}))).toMatchObject({ ok: true });
      expect(client.create).toHaveBeenCalledTimes(2);
      expect(client.release).toHaveBeenCalledTimes(1);

      expect(await owner.invoke(request('solve', 's', {}))).toMatchObject({ ok: true });
      const solvedReference = client.solve.mock.calls.at(-1)?.[0];
      const releasedReference = client.release.mock.calls.at(-1)?.[0];
      expect(solvedReference).not.toBe(releasedReference);
    } finally {
      await owner.dispose();
    }
  });

  it('drops prior ownership after a failed apply instead of reviving stale state', async () => {
    const { owner, client } = setup();
    try {
      expect(await owner.invoke(request('apply', 'first', {}))).toMatchObject({ ok: true });
      client.create.mockRejectedValueOnce(new ModelClientError(Code.InvalidArgument, 'rpc.failed'));
      expect(await owner.invoke(request('apply', 'bad', {}))).toMatchObject({
        ok: false,
        error: { code: Code.InvalidArgument, reason: 'rpc.failed' },
      });
      expect(client.disconnect).toHaveBeenCalledTimes(1);
      expect(await owner.invoke(request('solve', 'after', {}))).toMatchObject({
        ok: false,
        error: { code: Code.NotFound, reason: 'workbench.no_model' },
      });
    } finally {
      await owner.dispose();
    }
  });

  it('rejects renderer attempts to inject model/session ownership fields', async () => {
    const { owner, client } = setup();
    try {
      for (const raw of [
        null,
        { ...request('apply', 'x', {}), model: 'renderer-token' },
        { ...request('solve', 'x', {}), sessionId: 'renderer-session' },
        { ...request('release', 'x'), input: {} },
        { ...request('apply', 'x', {}), version: 'v0' },
      ]) {
        expect(await owner.invoke(raw)).toMatchObject({
          ok: false,
          error: { code: Code.InvalidArgument },
        });
      }
      expect(client.create).not.toHaveBeenCalled();
      expect(client.solve).not.toHaveBeenCalled();
      expect(client.release).not.toHaveBeenCalled();
    } finally {
      await owner.dispose();
    }
  });
});
