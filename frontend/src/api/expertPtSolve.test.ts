import { clone, create, fromJson } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';
import { describe, expect, it, vi } from 'vitest';
import { ModelClientError } from './modelSessionClient';
import { buildExpertPtRequest, emptyExpertPtInput, ExpertPtSolveGate, type ExpertPtInput } from './expertPtSolve';
import { FullPtResultSchema, ModelSnapshotSchema, type FullPtResult } from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { expertSnapshotJson } from '../test/modelInspectorFixtures';

function snapshot() { return fromJson(ModelSnapshotSchema, expertSnapshotJson()); }
function input(): ExpertPtInput {
  return { pressurePa: '2e6', temperatureK: '300', feed: [
    { componentId: 'methane', fraction: '0.25' }, { componentId: 'ethane', fraction: '0.75' },
  ] };
}
function deferred<T>() {
  let resolve!: (value: T) => void;
  let reject!: (cause: unknown) => void;
  const promise = new Promise<T>((yes, no) => { resolve = yes; reject = no; });
  return { promise, resolve, reject };
}

describe('custom-model PT input and publication', () => {
  it('uses the applied component order and never invents missing PT/feed values', () => {
    const model = snapshot();
    expect(emptyExpertPtInput(model)).toEqual({ pressurePa: '', temperatureK: '', feed: [
      { componentId: 'methane', fraction: '' }, { componentId: 'ethane', fraction: '' },
    ] });
    expect(buildExpertPtRequest(model, input())).toEqual({ pressurePa: 2e6, temperatureK: 300, feed: [0.25, 0.75] });
    expect(() => buildExpertPtRequest(model, emptyExpertPtInput(model))).toThrow();
  });

  it('rejects stale component alignment after reorder/add/remove and resets new revisions', () => {
    const model = snapshot();
    const reordered = clone(ModelSnapshotSchema, model);
    reordered.definition!.components.reverse();
    expect(emptyExpertPtInput(reordered).feed.map(item => item.componentId)).toEqual(['ethane', 'methane']);
    expect(() => buildExpertPtRequest(reordered, input())).toThrow();
    expect(buildExpertPtRequest(reordered, { ...input(), feed: [...input().feed].reverse() })).toMatchObject({ feed: [0.75, 0.25] });
    reordered.definition!.components.pop();
    expect(emptyExpertPtInput(reordered).feed).toEqual([{ componentId: 'ethane', fraction: '' }]);
    expect(() => buildExpertPtRequest(reordered, input())).toThrow();
    reordered.definition!.components.push({ ...model.definition!.components[0]!, componentId: 'new-component' });
    expect(() => buildExpertPtRequest(reordered, input())).toThrow();
  });

  it('rejects missing/nonfinite/nondecimal/invalid scalars with structured field locations', () => {
    for (const value of ['', ' ', 'NaN', 'Infinity', '0x10', '1e999', '0', '-1']) {
      try { buildExpertPtRequest(snapshot(), { ...input(), pressurePa: value }); throw new Error('unexpected acceptance'); }
      catch (cause) { expect(cause).toMatchObject({ reason: 'expert.invalid_pt_input', validation: { field: 'pressure_pa' } }); }
    }
    for (const value of ['', 'NaN', '-0.1', '1.1']) {
      expect(() => buildExpertPtRequest(snapshot(), { ...input(), feed: [input().feed[0]!, { componentId: 'ethane', fraction: value }] })).toThrow();
    }
  });

  it('does not copy a normalization tolerance or silently normalize the supplied feed', () => {
    const state = { ...input(), feed: input().feed.map(entry => ({ ...entry, fraction: '0.4' })) };
    expect(buildExpertPtRequest(snapshot(), state)).toMatchObject({ feed: [0.4, 0.4] });
    expect(state.feed.map(entry => entry.fraction)).toEqual(['0.4', '0.4']);
    // This invalid sum is intentionally left for the existing native validator.
  });

  it('keeps one admitted solve while edits suppress stale results without cancelling', async () => {
    const pending = deferred<FullPtResult>();
    const solve = vi.fn(() => pending.promise);
    const gate = new ExpertPtSolveGate();
    const signal = new AbortController().signal;
    const request = buildExpertPtRequest(snapshot(), input());
    const running = gate.solve({ solve }, request, { signal });
    expect(gate.busy).toBe(true);
    gate.invalidate(); // P/T/feed/model draft changed during the native solve.
    await expect(gate.solve({ solve }, request)).rejects.toMatchObject({ reason: 'expert.solve_busy' });
    expect(solve).toHaveBeenCalledTimes(1);
    expect(signal.aborted).toBe(false);
    pending.resolve(create(FullPtResultSchema, { diagnostic: 'old-input' }));
    expect(await running).toBeNull();
    expect(gate.busy).toBe(false);
    const fresh = create(FullPtResultSchema, { diagnostic: 'new-input' });
    expect(await gate.solve({ solve: async () => fresh }, request)).toBe(fresh);
  });

  it('ignores late failures after unmount and preserves current typed errors', async () => {
    const pending = deferred<FullPtResult>();
    const gate = new ExpertPtSolveGate();
    const running = gate.solve({ solve: () => pending.promise }, {});
    gate.invalidate();
    pending.reject(new Error('late internal failure'));
    expect(await running).toBeNull();
    const failure = new ModelClientError(Code.InvalidArgument, 'rpc.failed');
    await expect(gate.solve({ solve: async () => { throw failure; } }, {})).rejects.toBe(failure);
    expect(gate.busy).toBe(false);
  });
});
