import { Code } from '@connectrpc/connect';
import { describe, expect, it, vi } from 'vitest';

import type { ModelWorkbenchBridge, ModelWorkbenchReply } from './modelWorkbenchContract';
import { MODEL_WORKBENCH_CONVENTION } from './modelWorkbenchContract';
import { ModelWorkbenchOwner } from './modelWorkbenchOwner';
import { RendererModelError } from './rendererModelWire';
import { expertResultJson, expertSnapshotJson } from '../test/modelInspectorFixtures';

function ok(value: unknown): ModelWorkbenchReply {
  return { version: MODEL_WORKBENCH_CONVENTION, ok: true, value: value as never };
}
function failure(code: Code, reason: string): ModelWorkbenchReply {
  return { version: MODEL_WORKBENCH_CONVENTION, ok: false, error: { code, reason } };
}

function bridge() {
  const apply = vi.fn(async () => ok(expertSnapshotJson()));
  const solve = vi.fn(async () => ok(expertResultJson()));
  const release = vi.fn(async () => ok(null));
  const cancel = vi.fn();
  const value: ModelWorkbenchBridge = {
    convention: MODEL_WORKBENCH_CONVENTION,
    apply,
    solve,
    release,
    cancel,
  };
  return { value, apply, solve, release, cancel };
}

const state = { pressurePa: 10_000_000, temperatureK: 330, feed: [0.7, 0.3] };

describe('renderer-safe model workbench owner', () => {
  it('applies and solves without any renderer-visible model token', async () => {
    const fake = bridge();
    const owner = new ModelWorkbenchOwner(fake.value);
    const model = await owner.create({});
    const result = await model.solve(state);

    expect(model.snapshot.definition?.datasetId).toBe('expert-inspector-fixture');
    expect(result.pressurePa).toBe(state.pressurePa);
    expect(result.feed).toEqual(state.feed);
    expect(fake.apply).toHaveBeenCalledTimes(1);
    expect(fake.solve).toHaveBeenCalledTimes(1);
    expect(fake.apply.mock.calls[0]).toHaveLength(2);
    expect(fake.solve.mock.calls[0]).toHaveLength(2);
    expect(JSON.stringify(fake.apply.mock.calls)).not.toContain('modelHandle');
    expect(JSON.stringify(fake.solve.mock.calls)).not.toContain('modelHandle');
  });

  it('makes replaced renderer objects stale without releasing the replacement', async () => {
    const fake = bridge();
    const owner = new ModelWorkbenchOwner(fake.value);
    const first = await owner.create({});
    const second = await owner.create({});

    await first.release();
    expect(fake.release).not.toHaveBeenCalled();
    await expect(first.solve(state)).rejects.toMatchObject({
      code: Code.NotFound,
      reason: 'workbench.stale_model',
    });
    await expect(first.source.describe()).rejects.toMatchObject({
      code: Code.NotFound,
      reason: 'workbench.stale_model',
    });

    await second.solve(state);
    await second.release();
    expect(fake.release).toHaveBeenCalledTimes(1);
  });

  it('invalidates prior ownership after a failed apply', async () => {
    const fake = bridge();
    const owner = new ModelWorkbenchOwner(fake.value);
    const previous = await owner.create({});
    fake.apply.mockImplementationOnce(async () => failure(Code.InvalidArgument, 'ipc.invalid_input'));

    await expect(owner.create({})).rejects.toBeInstanceOf(RendererModelError);
    await expect(previous.solve(state)).rejects.toMatchObject({
      code: Code.NotFound,
      reason: 'workbench.stale_model',
    });
    await previous.release();
    expect(fake.release).not.toHaveBeenCalled();
  });

  it('retains the current model after ordinary solve validation errors', async () => {
    const fake = bridge();
    const owner = new ModelWorkbenchOwner(fake.value);
    const model = await owner.create({});
    fake.solve.mockImplementationOnce(async () => failure(Code.InvalidArgument, 'ipc.invalid_input'));

    await expect(model.solve(state)).rejects.toMatchObject({ code: Code.InvalidArgument });
    expect((await model.solve(state)).feed).toEqual(state.feed);
  });
});
