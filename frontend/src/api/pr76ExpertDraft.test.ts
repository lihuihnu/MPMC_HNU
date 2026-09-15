import { create, fromJson } from '@bufbuild/protobuf';
import { describe, expect, it } from 'vitest';

import type { ModelCreateInput } from './modelSessionClient';
import {
  addPr76Component,
  buildPr76CreateInput,
  editPr76ComponentScalar,
  editPr76ComponentText,
  editPr76Kij,
  editPr76Setting,
  movePr76Component,
  newPr76ExpertDraft,
  pr76ExpertDraftFromSnapshot,
  removePr76Component,
  setPr76SolverMode,
  validatePr76ExpertDraft,
  Pr76DraftValidationError,
} from './pr76ExpertDraft';
import {
  ComponentKind,
  CreateModelRequestSchema,
  ModelSnapshotSchema,
  SolverSettingsKind,
  SourceKind,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { expertSnapshotJson } from '../test/modelInspectorFixtures';

function identity<T extends ReturnType<typeof newPr76ExpertDraft>>(draft: T): T {
  return {
    ...draft,
    displayName: 'User PR76 model',
    datasetId: 'user-pr76-dataset',
    revision: 'user-r1',
    userRecord: {
      reference: 'user-record-2026-09-15',
      revision: 'record-r1',
      locator: 'expert-ui-session-1',
      note: 'Entered explicitly in the Expert editor.',
      acquisition: 'interactive Expert UI entry',
      usageTerms: 'user-supplied for this project',
    },
  };
}

function request(input: ModelCreateInput) {
  return create(CreateModelRequestSchema, input);
}

function pr76Request(input: ModelCreateInput) {
  const value = request(input);
  const definition = value.definition;
  if (!definition || definition.parameters.case !== 'pr76') {
    throw new Error('materialized request is missing PR76 definition');
  }
  return { request: value, definition, pr76: definition.parameters.value };
}

function addComponent(
  draft: ReturnType<typeof newPr76ExpertDraft>,
  id: string,
  name: string,
  tc: string,
  pc: string,
  omega: string,
) {
  let next = addPr76Component(draft);
  const component = next.components.at(-1)!;
  next = editPr76ComponentText(next, component.key, 'componentId', id);
  next = editPr76ComponentText(next, component.key, 'displayName', name);
  next = editPr76ComponentScalar(next, component.key, 'criticalTemperatureK', tc);
  next = editPr76ComponentScalar(next, component.key, 'criticalPressurePa', pc);
  next = editPr76ComponentScalar(next, component.key, 'acentricFactor', omega);
  return next;
}

describe('PR76 Expert draft/create contract', () => {
  it('fails closed for an empty draft and never invents missing kij', () => {
    const empty = newPr76ExpertDraft();
    expect(validatePr76ExpertDraft(empty).map(issue => issue.field)).toContain('components');
    expect(() => buildPr76CreateInput(empty)).toThrow(Pr76DraftValidationError);

    let draft = identity(newPr76ExpertDraft());
    draft = addComponent(draft, 'methane', 'Methane', '190.564', '4599200', '0.01142');
    draft = addComponent(draft, 'ethane', 'Ethane', '305.322', '4872200', '0.0995');
    expect(draft.pairs).toHaveLength(1);
    expect(draft.pairs[0]!.kij.text).toBe('');
    expect(validatePr76ExpertDraft(draft)).toContainEqual(expect.objectContaining({
      field: 'parameters.binary[0].kij',
    }));
    expect(() => buildPr76CreateInput(draft)).toThrow(Pr76DraftValidationError);

    draft = editPr76Kij(draft, draft.pairs[0]!.key, '0');
    const { request: value, definition, pr76 } = pr76Request(buildPr76CreateInput(draft));
    expect(value.solverSelection).toEqual({ case: 'presetId', value: 'mpmc-balanced-default/v1' });
    expect(definition.components.map(component => component.componentId)).toEqual(['methane', 'ethane']);
    expect(pr76.binary).toHaveLength(1);
    expect(pr76.binary[0]?.kij?.value).toBe(0);
    expect(pr76.binary[0]?.kij?.provenance?.kind).toBe(SourceKind.USER_SUPPLIED);
    expect(definition.provenance?.reference).toBe('user-record-2026-09-15');
    expect(definition.applicability?.temperatureLowerK).toBeUndefined();
    expect(definition.applicability?.provenance?.kind).toBe(SourceKind.USER_SUPPLIED);
  });

  it('keeps pair records explicit through reorder/remove/add operations', () => {
    let draft = identity(newPr76ExpertDraft());
    draft = addComponent(draft, 'a', 'A', '200', '4000000', '0.1');
    draft = addComponent(draft, 'b', 'B', '300', '5000000', '0.2');
    draft = editPr76Kij(draft, draft.pairs[0]!.key, '0.04');
    const a = draft.components[0]!;
    const b = draft.components[1]!;

    draft = movePr76Component(draft, b.key, -1);
    expect(draft.components.map(component => component.componentId)).toEqual(['b', 'a']);
    expect(draft.pairs).toHaveLength(1);
    expect(draft.pairs[0]!.kij.text).toBe('0.04');
    const { pr76 } = pr76Request(buildPr76CreateInput(draft));
    expect(pr76.binary[0]).toMatchObject({
      firstComponentId: 'a', secondComponentId: 'b',
    });

    draft = removePr76Component(draft, a.key);
    expect(draft.pairs).toHaveLength(0);
    draft = addComponent(draft, 'c', 'C', '350', '5500000', '0.3');
    expect(draft.pairs).toHaveLength(1);
    expect(draft.pairs[0]!.kij.text).toBe('');
    expect(() => editPr76ComponentText(draft, b.key, 'componentId', 'renamed-existing')).toThrow(
      'Existing component identity is immutable',
    );
  });

  it('preserves unchanged provenance and marks only edited values as user supplied', () => {
    const snapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
    let draft = identity(pr76ExpertDraftFromSnapshot(snapshot));
    let materialized = pr76Request(buildPr76CreateInput(draft));
    expect(materialized.pr76.pure[0]?.criticalTemperatureK?.provenance?.kind).toBe(SourceKind.LITERATURE);
    expect(materialized.pr76.binary[0]?.kij?.provenance?.kind).toBe(SourceKind.LITERATURE);
    expect(materialized.definition.provenance?.kind).toBe(SourceKind.USER_SUPPLIED);
    expect(materialized.definition.applicability?.provenance?.kind).toBe(SourceKind.LITERATURE);

    const methane = draft.components[0]!;
    draft = editPr76ComponentScalar(draft, methane.key, 'criticalTemperatureK', '191');
    draft = editPr76Kij(draft, draft.pairs[0]!.key, '0.02');
    materialized = pr76Request(buildPr76CreateInput(draft));
    expect(materialized.pr76.pure[0]?.criticalTemperatureK).toMatchObject({
      value: 191,
      originalUnit: 'K',
      conversion: 'identity: Expert UI canonical SI input',
      provenance: { kind: SourceKind.USER_SUPPLIED, reference: 'user-record-2026-09-15' },
    });
    expect(materialized.pr76.pure[0]?.criticalPressurePa?.provenance?.kind).toBe(SourceKind.LITERATURE);
    expect(materialized.pr76.binary[0]?.kij?.provenance?.kind).toBe(SourceKind.USER_SUPPLIED);
  });

  it('supports complete custom settings only after a resolved snapshot exists', () => {
    expect(() => setPr76SolverMode(newPr76ExpertDraft(), 'custom')).toThrow(
      'Custom settings require a resolved settings snapshot first',
    );
    const snapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
    let draft = identity(pr76ExpertDraftFromSnapshot(snapshot));
    draft = setPr76SolverMode(draft, 'custom');
    draft = editPr76Setting(draft, 'eosRoot', 'maxIterations', '17');
    draft = editPr76Setting(draft, 'initialStability', 'automaticMultistart', true);
    draft = editPr76Setting(draft, 'initialStability', 'tpdTolerance', '1e-8');
    const value = request(buildPr76CreateInput(draft));
    expect(value.solverSelection.case).toBe('settings');
    if (value.solverSelection.case !== 'settings') throw new Error('missing custom settings');
    const settings = value.solverSelection.value;
    expect(settings.kind).toBe(SolverSettingsKind.CUSTOM);
    expect(settings.presetId).toBe('');
    expect(settings.eosRoot?.maxIterations).toBe(17n);
    expect(settings.initialStability?.automaticMultistart).toBe(true);
    expect(settings.initialStability?.tpdTolerance).toBe(1e-8);
    expect(snapshot.settings?.eosRoot?.maxIterations).toBe(9007199254740993n);
  });

  it('validates finite/positive domains and component identity locally before mutation', () => {
    let draft = identity(newPr76ExpertDraft());
    draft = addComponent(draft, 'bad id', 'Bad', '0', 'NaN', 'Infinity');
    const issues = validatePr76ExpertDraft(draft);
    expect(issues).toEqual(expect.arrayContaining([
      expect.objectContaining({ field: 'components[0].component_id' }),
      expect.objectContaining({ field: 'components[0].critical_temperature_k' }),
      expect.objectContaining({ field: 'components[0].critical_pressure_pa' }),
      expect.objectContaining({ field: 'components[0].acentric_factor' }),
    ]));
    expect(draft.components[0]?.kind).toBe(ComponentKind.PURE);
  });
});
