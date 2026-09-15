import { equals, fromJson } from '@bufbuild/protobuf';
import { describe, expect, expectTypeOf, it, vi } from 'vitest';

import type { ModelReference, ModelSessionClient } from './modelSessionClient';
import type { RendererModelClient, RendererModelReference } from './rendererModelClient';
import { bindExpertModelSource, type ExpertModelSource } from './expertModelSource';
import { ModelSnapshotSchema } from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { expertSnapshotJson } from '../test/modelInspectorFixtures';

function bindWeb(client: ModelSessionClient, model: ModelReference): ExpertModelSource {
  return bindExpertModelSource(client, model);
}
function bindElectron(client: RendererModelClient, model: RendererModelReference): ExpertModelSource {
  return bindExpertModelSource(client, model);
}

describe('Expert model owner handoff', () => {
  it('accepts both nominal typed reference/client pairs without exposing the reference', () => {
    expectTypeOf(bindWeb).returns.toMatchTypeOf<ExpertModelSource>();
    expectTypeOf(bindElectron).returns.toMatchTypeOf<ExpertModelSource>();
  });

  it('calls the owning describe path and retains the complete isolated snapshot', async () => {
    const snapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
    const reference = Object.freeze({ owner: 'session-a' });
    const describeModel = vi.fn(async (model: typeof reference) => {
      expect(model).toBe(reference);
      return snapshot;
    });
    const source = bindExpertModelSource({ describe: describeModel }, reference);

    expect(Object.keys(source)).toEqual(['describe']);
    const inspection = await source.describe();
    expect(describeModel).toHaveBeenCalledTimes(1);
    expect(equals(ModelSnapshotSchema, inspection.snapshot, snapshot)).toBe(true);

    snapshot.definition!.displayName = 'caller mutation';
    expect(inspection.snapshot.definition?.displayName).toBe('Read-only PR76 fixture');
  });

  it('forwards cancellation/timeout options without reconnecting or retrying', async () => {
    const snapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
    const controller = new AbortController();
    const describeModel = vi.fn(async (_model: symbol, options?: { signal?: AbortSignal; timeoutMs?: number }) => {
      expect(options?.signal).toBe(controller.signal);
      expect(options?.timeoutMs).toBe(4321);
      return snapshot;
    });
    const source = bindExpertModelSource({ describe: describeModel }, Symbol('owned-model'));

    await source.describe({ signal: controller.signal, timeoutMs: 4321 });
    expect(describeModel).toHaveBeenCalledTimes(1);
  });
});
