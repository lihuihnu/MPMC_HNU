import { create, fromJson } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';
import { renderToStaticMarkup } from 'react-dom/server';
import { describe, expect, it, vi } from 'vitest';

import type { ExpertModelOwner, ExpertOwnedModel } from '../api/expertModelOwner';
import { ModelClientError } from '../api/modelSessionClient';
import { MODEL_VALIDATION_DETAIL_VERSION } from '../api/modelValidationDetail';
import { FullPtResultSchema, ModelSnapshotSchema } from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { expertSnapshotJson } from '../test/modelInspectorFixtures';
import {
  ExpertPr76WorkspaceView,
  deferExpertWorkspaceRelease,
  expertWorkspaceFailure,
  retirePreviousExpertModel,
} from './ExpertPr76Workspace';

function owned(
  releaseWork: () => Promise<void> = async () => {},
): ExpertOwnedModel {
  const snapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
  let released = false;
  return Object.freeze({
    snapshot,
    source: { describe: async () => ({ snapshot }) },
    get released() { return released; },
    solve: async () => create(FullPtResultSchema),
    async release() {
      if (released) return;
      released = true;
      await releaseWork();
    },
  });
}

function owner(): ExpertModelOwner {
  return { create: vi.fn(async () => owned()) };
}

function render(
  current: ExpertOwnedModel | null,
  retirementFailure: ReturnType<typeof expertWorkspaceFailure> | null = null,
  unapplied = false,
) {
  return renderToStaticMarkup(
    <ExpertPr76WorkspaceView
      owner={owner()}
      current={current}
      generation={current ? 1 : 0}
      retirementFailure={retirementFailure}
      onCreated={() => {}}
      unapplied={unapplied}
    />,
  );
}

describe('PR76 workspace ownership and product presentation', () => {
  it('starts with editable PR fluid data and no invented live result', () => {
    const html = render(null);
    expect(html).toContain('data-expert-workspace="pr76"');
    expect(html).toContain('Peng–Robinson flash');
    expect(html).toContain('Define PR fluid');
    expect(html).toContain('Create a PR fluid model');
    expect(html).toContain('Add at least two components to define a binary interaction coefficient');
    expect(html).not.toContain('Read-only PR76 fixture');
    expect(html).not.toContain('data-expert-pt-solve');
    expect(html).toContain('Local calculation · no login');
  });

  it('derives the next immutable revision and exposes product-level flash inputs', () => {
    const html = render(owned());
    expect(html).toContain('Edit PR fluid');
    expect(html).toContain('<summary>Methane</summary>');
    expect(html).toContain('kij methane ↔ ethane');
    expect(html).toContain('Reading live model snapshot');
    expect(html).not.toContain('Create a PR fluid model');
    expect(html).toContain('Run PR flash');
    expect(html).toContain('PR feed methane');
    expect(html).toContain('PR feed ethane');
    expect(html).toContain('Advanced: applied model parameters and provenance');
  });

  it('blocks the applied solve while edits are not applied and does not expose native handles', () => {
    const html = render(owned(), null, true);
    expect(html).toContain('Fluid data changed');
    expect(html).toContain('type="submit" disabled=""');
    expect(html).not.toContain('modelHandle');
    expect(html).not.toContain('data-expert-result');
  });

  it('retires the previous model only after a distinct next model exists', async () => {
    const release = vi.fn(async () => {});
    const previous = owned(release);
    const next = owned();
    expect(await retirePreviousExpertModel(null, next)).toBeNull();
    expect(release).not.toHaveBeenCalled();
    expect(await retirePreviousExpertModel(previous, previous)).toBeNull();
    expect(release).not.toHaveBeenCalled();
    expect(await retirePreviousExpertModel(previous, next)).toBeNull();
    expect(release).toHaveBeenCalledTimes(1);
    expect(previous.released).toBe(true);
    expect(await retirePreviousExpertModel(previous, next)).toBeNull();
    expect(release).toHaveBeenCalledTimes(1);
  });

  it('cancels StrictMode replay cleanup and releases only the current model on real teardown', () => {
    const firstRelease = vi.fn(async () => {});
    const secondRelease = vi.fn(async () => {});
    const first = owned(firstRelease);
    const second = owned(secondRelease);
    let current: ExpertOwnedModel | null = first;
    let stillUnmounted = false;
    const queued: Array<() => void> = [];

    deferExpertWorkspaceRelease(
      () => current,
      () => stillUnmounted,
      (task) => queued.push(task),
    );
    current = second;
    queued.shift()?.();
    expect(firstRelease).not.toHaveBeenCalled();
    expect(secondRelease).not.toHaveBeenCalled();

    stillUnmounted = true;
    deferExpertWorkspaceRelease(
      () => current,
      () => stillUnmounted,
      (task) => queued.push(task),
    );
    queued.shift()?.();
    expect(firstRelease).not.toHaveBeenCalled();
    expect(secondRelease).toHaveBeenCalledTimes(1);
    expect(second.released).toBe(true);
  });

  it('surfaces typed ambiguous cleanup failure without leaking arbitrary error text', async () => {
    const validation = Object.freeze({
      version: MODEL_VALIDATION_DETAIL_VERSION,
      code: 'request.rejected',
      field: 'request',
    });
    const previous = owned(async () => {
      throw new ModelClientError(Code.Unavailable, 'rpc.failed', validation);
    });
    const failure = await retirePreviousExpertModel(previous, owned());
    expect(failure).toEqual({ reason: 'rpc.failed', code: Code.Unavailable, validation });
    const html = render(owned(), failure);
    expect(html).toContain('The previous model could not be fully released');
    expect(html).toContain('Technical details');
    expect(html).toContain('rpc.failed');
    expect(html).toContain('request');

    const unknown = expertWorkspaceFailure(new Error('private-session-token'));
    expect(unknown).toEqual({ reason: 'expert.release_failed' });
    expect(JSON.stringify(unknown)).not.toContain('private-session-token');
  });
});
