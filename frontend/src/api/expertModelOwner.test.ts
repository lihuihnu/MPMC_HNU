import { equals, fromJson } from '@bufbuild/protobuf';
import { describe, expect, expectTypeOf, it, vi } from 'vitest';

import type { ModelReference, ModelSessionClient } from './modelSessionClient';
import type { RendererModelClient, RendererModelReference } from './rendererModelClient';
import {
  bindExpertModelOwner,
  type ExpertModelOwner,
} from './expertModelOwner';
import {
  ModelSnapshotSchema,
  type ModelSnapshot,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { expertSnapshotJson } from '../test/modelInspectorFixtures';

function bindWeb(client: ModelSessionClient): ExpertModelOwner {
  return bindExpertModelOwner(client);
}
function bindElectron(client: RendererModelClient): ExpertModelOwner {
  return bindExpertModelOwner(client);
}

describe('Expert immutable model ownership adapter', () => {
  it('accepts both typed client ownership surfaces', () => {
    expectTypeOf(bindWeb).returns.toMatchTypeOf<ExpertModelOwner>();
    expectTypeOf(bindElectron).returns.toMatchTypeOf<ExpertModelOwner>();
    expectTypeOf<ModelReference>().not.toEqualTypeOf<RendererModelReference>();
  });

  it('hides the reference, exposes live describe and releases exactly once', async () => {
    const snapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
    const reference = Object.freeze({ session: 'opaque-reference' });
    const create = vi.fn(async () => ({ model: reference, snapshot }));
    const describe = vi.fn(async (model: typeof reference) => {
      expect(model).toBe(reference);
      return snapshot;
    });
    const release = vi.fn(async (model: typeof reference) => {
      expect(model).toBe(reference);
    });
    const owner = bindExpertModelOwner({ create, describe, release });
    const made = await owner.create({});

    expect(Object.keys(made)).toEqual(['snapshot', 'source', 'released', 'release']);
    expect(equals(ModelSnapshotSchema, made.snapshot, snapshot)).toBe(true);
    snapshot.definition!.displayName = 'mutated caller object';
    expect(made.snapshot.definition?.displayName).toBe('Read-only PR76 fixture');
    expect((await made.source.describe()).snapshot.definition?.displayName).toBe('mutated caller object');
    expect(describe).toHaveBeenCalledTimes(1);

    expect(made.released).toBe(false);
    await made.release();
    expect(made.released).toBe(true);
    await made.release();
    expect(release).toHaveBeenCalledTimes(1);
  });

  it('does not retry a failed remote release after local ownership is consumed', async () => {
    const snapshot: ModelSnapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
    const release = vi.fn(async () => { throw new Error('transport lost after release admission'); });
    const owner = bindExpertModelOwner({
      create: async () => ({ model: Symbol('opaque'), snapshot }),
      describe: async () => snapshot,
      release,
    });
    const made = await owner.create({});
    await expect(made.release()).rejects.toThrow('transport lost');
    expect(made.released).toBe(true);
    await made.release();
    expect(release).toHaveBeenCalledTimes(1);
  });
});
