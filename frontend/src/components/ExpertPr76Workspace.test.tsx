import { fromJson } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';
import { renderToStaticMarkup } from 'react-dom/server';
import { describe, expect, it, vi } from 'vitest';

import type { ExpertModelOwner, ExpertOwnedModel } from '../api/expertModelOwner';
import { ModelClientError } from '../api/modelSessionClient';
import { MODEL_VALIDATION_DETAIL_VERSION } from '../api/modelValidationDetail';
import { ModelSnapshotSchema } from '../gen/mpmc/model_configuration/v1/model_service_pb';
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
) {
  return renderToStaticMarkup(
    <ExpertPr76WorkspaceView
      owner={owner()}
      current={current}
      generation={current ? 1 : 0}
      retirementFailure={retirementFailure}
      onCreated={() => {}}
    />,
  );
}

describe('PR76 Expert workspace ownership and presentation', () => {
  it('starts with the explicit editor and no invented live model', () => {
    const html = render(null);
    expect(html).toContain('data-expert-workspace="pr76"');
    expect(html).toContain('Create custom PR76 model');
    expect(html).toContain('No Expert model created yet');
    expect(html).toContain('No unordered component pair exists yet');
    expect(html).not.toContain('Read-only PR76 fixture');
  });

  it('derives the next immutable revision from the current owned snapshot and exposes live describe loading', () => {
    const html = render(owned());
    expect(html).toContain('Create immutable revision');
    expect(html).toContain('0: methane');
    expect(html).toContain('kij methane ↔ ethane');
    expect(html).toContain('Reading live model snapshot');
    expect(html).not.toContain('No Expert model created yet');
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
    expect(html).toContain('Previous model cleanup was not confirmed');
    expect(html).toContain('rpc.failed');
    expect(html).toContain('request');

    const unknown = expertWorkspaceFailure(new Error('private-session-token'));
    expect(unknown).toEqual({ reason: 'expert.release_failed' });
    expect(JSON.stringify(unknown)).not.toContain('private-session-token');
  });
});
