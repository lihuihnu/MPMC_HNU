import { create, equals, fromJson, toBinary } from '@bufbuild/protobuf';
import { describe, expect, expectTypeOf, it, vi } from 'vitest';

import type { ModelCallOptions, ModelReference, ModelSessionClient, ModelSolveInput } from './modelSessionClient';
import type { RendererModelClient, RendererModelReference } from './rendererModelClient';
import { bindExpertModelOwner, type ExpertModelOwner } from './expertModelOwner';
import { FullPtResultSchema, ModelSnapshotSchema, type ModelSnapshot } from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { expertSnapshotJson } from '../test/modelInspectorFixtures';

function bindWeb(client: ModelSessionClient): ExpertModelOwner { return bindExpertModelOwner(client); }
function bindElectron(client: RendererModelClient): ExpertModelOwner { return bindExpertModelOwner(client); }

// Serialization/ownership fixtures only; physical references live in native tests.
describe('Expert immutable model ownership adapter', () => {
  it('accepts both typed client ownership surfaces', () => {
    expectTypeOf(bindWeb).returns.toMatchTypeOf<ExpertModelOwner>();
    expectTypeOf(bindElectron).returns.toMatchTypeOf<ExpertModelOwner>();
    expectTypeOf<ModelReference>().not.toEqualTypeOf<RendererModelReference>();
  });

  it('hides the reference, exposes live describe and releases exactly once', async () => {
    const snapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
    const reference = Object.freeze({ session: 'opaque-reference' });
    const make = vi.fn(async () => ({ model: reference, snapshot }));
    const describe = vi.fn(async (model: typeof reference) => { expect(model).toBe(reference); return snapshot; });
    const release = vi.fn(async (model: typeof reference) => { expect(model).toBe(reference); });
    const solve = vi.fn(async () => create(FullPtResultSchema));
    const owner = bindExpertModelOwner({ create: make, describe, release, solve });
    const made = await owner.create({});

    expect(Object.keys(made)).toEqual(['snapshot', 'source', 'released', 'solve', 'release']);
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
    await expect(made.solve({})).rejects.toMatchObject({ reason: 'model.stale_reference' });
    expect(solve).not.toHaveBeenCalled();
  });

  it('does not retry a failed remote release after local ownership is consumed', async () => {
    const snapshot: ModelSnapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
    const release = vi.fn(async () => { throw new Error('transport lost after release admission'); });
    const owner = bindExpertModelOwner({
      create: async () => ({ model: Symbol('opaque'), snapshot }),
      describe: async () => snapshot,
      solve: async () => create(FullPtResultSchema),
      release,
    });
    const made = await owner.create({});
    await expect(made.release()).rejects.toThrow('transport lost');
    expect(made.released).toBe(true);
    await made.release();
    expect(release).toHaveBeenCalledTimes(1);
  });

  it.each(['shared', 'renderer'] as const)('preserves the complete %s result, hints and call options without exposing the reference', async (kind) => {
    const snapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
    const reference = Symbol('private-model');
    const result = create(FullPtResultSchema, {
      diagnostic: 'serialization-only unresolved diagnostic', globalStabilityProven: false,
      feed: [0.25, 0.75], pressurePa: 2e6, temperatureK: 300,
      candidatePhaseSet: { phases: [{ molePhaseFraction: 1, composition: [0.25, 0.75], lnFugacityCoefficient: [0, 0], providerBranch: 123n }] },
    });
    const solve = vi.fn(async (_ref: typeof reference, _input: ModelSolveInput, _options?: ModelCallOptions) =>
      kind === 'renderer' ? { result } : result);
    const made = await bindExpertModelOwner({
      create: async () => ({ model: reference, snapshot }), describe: async () => snapshot,
      solve, release: async () => {},
    }).create({});
    const signal = new AbortController().signal;
    const input: ModelSolveInput = { pressurePa: 2e6, temperatureK: 300, feed: [0.25, 0.75],
      hints: { version: 'pt-solve-hints/v1', initialStabilityStarts: [{ composition: [0.25, 0.75] }] } };
    const options = { signal, timeoutMs: 1000 };
    const solved = await made.solve(input, options);
    expect(solve).toHaveBeenCalledExactlyOnceWith(reference, input, options);
    expect(toBinary(FullPtResultSchema, solved)).toEqual(toBinary(FullPtResultSchema, result));
    expect(solved).not.toBe(result);
    expect(Object.values(made)).not.toContain(reference);
    result.feed[0] = 0.99;
    expect(solved.feed).toEqual([0.25, 0.75]);
  });
});
