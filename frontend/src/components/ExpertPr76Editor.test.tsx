import { create, fromJson } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';
import { renderToStaticMarkup } from 'react-dom/server';
import { describe, expect, it, vi } from 'vitest';

import type { ExpertModelOwner, ExpertOwnedModel } from '../api/expertModelOwner';
import { ModelClientError, type ModelCallOptions, type ModelCreateInput } from '../api/modelSessionClient';
import { MODEL_VALIDATION_DETAIL_VERSION } from '../api/modelValidationDetail';
import { pr76ExpertDraftFromSnapshot, setPr76SolverMode } from '../api/pr76ExpertDraft';
import {
  CreateModelRequestSchema,
  ModelSnapshotSchema,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { expertSnapshotJson } from '../test/modelInspectorFixtures';
import {
  ExpertPr76Editor,
  createPr76ExpertModel,
  expertPr76CreateFailure,
} from './ExpertPr76Editor';

function userIdentity<T extends ReturnType<typeof pr76ExpertDraftFromSnapshot>>(draft: T): T {
  return {
    ...draft,
    displayName: 'Edited PR76 model',
    datasetId: 'edited-pr76',
    revision: 'r2',
    userRecord: {
      reference: 'expert-user-record',
      revision: 'record-r2',
      locator: 'ui-test',
      note: 'test entry',
      acquisition: 'interactive editor',
      usageTerms: 'test-only user supplied data',
    },
  };
}

function fakeOwned(snapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson())): ExpertOwnedModel {
  return {
    snapshot,
    source: { describe: async () => ({ snapshot }) },
    released: false,
    release: async () => {},
  };
}

function fakeOwner(): ExpertModelOwner {
  return {
    create: vi.fn(async (_input: ModelCreateInput, _options?: ModelCallOptions) => fakeOwned()),
  };
}

describe('PR76 Expert editor presentation and create boundary', () => {
  it('renders a blank explicit editor without inventing components or kij', () => {
    const html = renderToStaticMarkup(
      <ExpertPr76Editor owner={fakeOwner()} onCreated={() => {}} />,
    );
    expect(html).toContain('Create custom PR76 model');
    expect(html).toContain('mpmc-balanced-default/v1');
    expect(html).toContain('Add component');
    expect(html).toContain('No unordered component pair exists yet');
    expect(html).toContain('Create once with the preset to obtain a complete resolved snapshot');
    expect(html).not.toContain('kij methane');
  });

  it('renders a derived snapshot with frozen component identities, explicit kij and custom mode available', () => {
    const snapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
    const html = renderToStaticMarkup(
      <ExpertPr76Editor owner={fakeOwner()} seedSnapshot={snapshot} onCreated={() => {}} />,
    );
    expect(html).toContain('Create immutable revision');
    expect(html).toContain('0: methane');
    expect(html).toContain('1: ethane');
    expect(html).toContain('kij methane ↔ ethane');
    expect(html).toContain('<option value="custom">Complete custom snapshot</option>');
    expect(html).toContain('<option value="preset" selected="">Frozen preset</option>');
    expect(html).not.toContain('EOS root');
    expect(html.match(/readOnly=""/gu)?.length ?? 0).toBeGreaterThanOrEqual(2);
  });

  it('submits only the domain-built immutable request through ExpertModelOwner', async () => {
    const snapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
    const draft = setPr76SolverMode(userIdentity(pr76ExpertDraftFromSnapshot(snapshot)), 'custom');
    const createModel = vi.fn(async (
      _input: ModelCreateInput,
      _options?: ModelCallOptions,
    ) => fakeOwned(snapshot));
    const owner: ExpertModelOwner = { create: createModel };

    const made = await createPr76ExpertModel(owner, draft, { timeoutMs: 4321 });
    expect(made.snapshot).toBe(snapshot);
    expect(createModel).toHaveBeenCalledTimes(1);
    const [initializer, options] = createModel.mock.calls[0]!;
    const request = create(CreateModelRequestSchema, initializer);
    expect(options).toEqual({ timeoutMs: 4321 });
    expect(request.definition?.datasetId).toBe('edited-pr76');
    expect(request.definition?.components.map(component => component.componentId)).toEqual(['methane', 'ethane']);
    expect(request.solverSelection.case).toBe('settings');
    expect(JSON.stringify(initializer)).not.toContain('modelHandle');
  });

  it('preserves typed creation errors and sanitizes arbitrary exception text', () => {
    const validation = Object.freeze({
      version: MODEL_VALIDATION_DETAIL_VERSION,
      code: 'configuration.invalid_value',
      field: 'parameters.binary[methane,ethane].kij',
    });
    expect(expertPr76CreateFailure(
      new ModelClientError(Code.InvalidArgument, 'rpc.failed', validation),
    )).toEqual({ reason: 'rpc.failed', code: Code.InvalidArgument, validation });

    const unknown = expertPr76CreateFailure(new Error('private-host-token'));
    expect(unknown).toEqual({ reason: 'expert.create_failed' });
    expect(JSON.stringify(unknown)).not.toContain('private-host-token');
  });
});
