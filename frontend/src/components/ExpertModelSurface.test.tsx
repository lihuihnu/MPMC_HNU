import { fromJson } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';
import { renderToStaticMarkup } from 'react-dom/server';
import { describe, expect, it, vi } from 'vitest';

import { inspectionFromSnapshot, type ExpertModelInspection } from '../api/expertModelInspector';
import { ModelClientError } from '../api/modelSessionClient';
import { RendererModelError } from '../api/rendererModelClient';
import { MODEL_VALIDATION_DETAIL_VERSION } from '../api/modelValidationDetail';
import { ModelSnapshotSchema } from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { expertSnapshotJson } from '../test/modelInspectorFixtures';
import {
  ExpertModelSurfaceView,
  expertModelSurfaceFailure,
  readExpertModelSurface,
  type ExpertModelSurfaceState,
} from './ExpertModelSurface';

function html(state: ExpertModelSurfaceState): string {
  return renderToStaticMarkup(<ExpertModelSurfaceView state={state} onRefresh={() => {}} />);
}

function deferred<T>() {
  let resolve!: (value: T) => void;
  let reject!: (cause: unknown) => void;
  const promise = new Promise<T>((yes, no) => { resolve = yes; reject = no; });
  return { promise, resolve, reject };
}

describe('live read-only Expert surface presentation', () => {
  it('renders the complete live inspection only in ready state', () => {
    const snapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
    const markup = html({ status: 'ready', inspection: inspectionFromSnapshot(snapshot) });

    expect(markup).toContain('data-expert-state="ready"');
    expect(markup).toContain('Live <code>describe()</code> snapshot');
    expect(markup).toContain('Read-only PR76 fixture');
    expect(markup).toContain('0: methane');
    expect(markup).toContain('1: ethane');
    expect(markup).toContain('kij methane ↔ ethane');
    expect(markup).toContain('Refresh snapshot');
  });

  it('renders loading without inventing a snapshot', () => {
    const markup = html({ status: 'loading' });
    expect(markup).toContain('data-expert-state="loading"');
    expect(markup).toContain('Reading live model snapshot');
    expect(markup).not.toContain('Read-only PR76 fixture');
  });

  it('keeps typed validation detail and never exposes arbitrary exception text', () => {
    const validation = Object.freeze({
      version: MODEL_VALIDATION_DETAIL_VERSION,
      code: 'request.rejected',
      field: 'pressure_pa',
    });
    const rendererFailure = expertModelSurfaceFailure(
      new RendererModelError(Code.InvalidArgument, 'rpc.failed', 'ipc', validation),
    );
    const markup = html({ status: 'failed', failure: rendererFailure });

    expect(markup).toContain('data-expert-state="failed"');
    expect(markup).toContain('rpc.failed');
    expect(markup).toContain('request.rejected');
    expect(markup).toContain('pressure_pa');
    expect(markup).toContain('Retry describe');

    const webFailure = expertModelSurfaceFailure(
      new ModelClientError(Code.NotFound, 'model.stale_reference'),
    );
    expect(webFailure).toMatchObject({ code: Code.NotFound, reason: 'model.stale_reference' });

    const unknown = expertModelSurfaceFailure(new Error('private-session-token'));
    expect(unknown).toEqual({ reason: 'expert.describe_failed' });
    expect(JSON.stringify(unknown)).not.toContain('private-session-token');
  });

  it('skips a retired StrictMode setup and issues only the active describe without a cancellation signal', async () => {
    const inspection = inspectionFromSnapshot(fromJson(ModelSnapshotSchema, expertSnapshotJson()));
    const read = vi.fn(async () => inspection);
    const oldPublish = vi.fn();
    const publish = vi.fn();
    const source = { describe: read };
    readExpertModelSurface(source, oldPublish)();
    const stop = readExpertModelSurface(source, publish);
    await Promise.resolve();
    await Promise.resolve();
    expect(read).toHaveBeenCalledExactlyOnceWith();
    expect(oldPublish).not.toHaveBeenCalled();
    expect(publish).toHaveBeenCalledExactlyOnceWith({ status: 'ready', inspection });
    stop();
  });

  it.each(['success', 'failure'] as const)('discards late %s after source change/unmount without cancelling the owning model', async (outcome) => {
    const pending = deferred<ExpertModelInspection>();
    const describe = vi.fn(() => pending.promise);
    const publish = vi.fn();
    const stop = readExpertModelSurface({ describe }, publish);
    await Promise.resolve();
    expect(describe).toHaveBeenCalledExactlyOnceWith();
    stop();
    if (outcome === 'success') pending.resolve(inspectionFromSnapshot(fromJson(ModelSnapshotSchema, expertSnapshotJson())));
    else pending.reject(new Error('private late failure'));
    await Promise.resolve();
    await Promise.resolve();
    expect(publish).not.toHaveBeenCalled();
  });

  it('publishes current typed read errors through the existing safe error contract', async () => {
    const publish = vi.fn();
    readExpertModelSurface({ describe: async () => {
      throw new ModelClientError(Code.NotFound, 'model.stale_reference');
    } }, publish);
    await Promise.resolve();
    await Promise.resolve();
    expect(publish).toHaveBeenCalledExactlyOnceWith({ status: 'failed',
      failure: { code: Code.NotFound, reason: 'model.stale_reference' } });
  });
});
