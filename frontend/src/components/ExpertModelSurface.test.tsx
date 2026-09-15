import { fromJson } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';
import { renderToStaticMarkup } from 'react-dom/server';
import { describe, expect, it } from 'vitest';

import { inspectionFromSnapshot } from '../api/expertModelInspector';
import { ModelClientError } from '../api/modelSessionClient';
import { RendererModelError } from '../api/rendererModelClient';
import { MODEL_VALIDATION_DETAIL_VERSION } from '../api/modelValidationDetail';
import { ModelSnapshotSchema } from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { expertSnapshotJson } from '../test/modelInspectorFixtures';
import {
  ExpertModelSurfaceView,
  expertModelSurfaceFailure,
  type ExpertModelSurfaceState,
} from './ExpertModelSurface';

function html(state: ExpertModelSurfaceState): string {
  return renderToStaticMarkup(<ExpertModelSurfaceView state={state} onRefresh={() => {}} />);
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
});
