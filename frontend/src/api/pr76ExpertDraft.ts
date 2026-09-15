import { clone, create } from '@bufbuild/protobuf';

import type { ModelCreateInput } from './modelSessionClient';
import {
  ComponentDefinitionSchema,
  ComponentKind,
  ModelApplicabilitySchema,
  ModelFamily,
  ModelProvenanceSchema,
  ModelScalarSchema,
  Pr76BinaryInteractionSchema,
  Pr76ParameterDefinitionSchema,
  Pr76PureParametersSchema,
  PtSolverSettingsSchema,
  SolverSettingsKind,
  SourceKind,
  ThermodynamicModelDefinitionSchema,
  type ComponentDefinition,
  type ModelApplicability,
  type ModelProvenance,
  type ModelScalar,
  type ModelSnapshot,
  type Pr76BinaryInteraction,
  type PtSolverSettings,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';

export const PR76_MODEL_DEFINITION_VERSION = 'thermodynamic-model/parameter-definition/v1' as const;
export const PR76_BALANCED_PRESET_ID = 'mpmc-balanced-default/v1' as const;

export interface ExpertUserRecordDraft {
  reference: string;
  revision: string;
  locator: string;
  note: string;
  acquisition: string;
  usageTerms: string;
}

export interface Pr76ScalarDraft {
  text: string;
  dirty: boolean;
  readonly canonicalUnit: 'K' | 'Pa' | 'dimensionless' | 'kg/mol';
  readonly positive: boolean;
  readonly optional: boolean;
  readonly original?: ModelScalar;
}

export interface Pr76ComponentDraft {
  readonly key: string;
  componentId: string;
  displayName: string;
  kind: ComponentKind;
  definitionDirty: boolean;
  readonly originalDefinition?: ComponentDefinition;
  molarMassKgPerMol: Pr76ScalarDraft;
  criticalTemperatureK: Pr76ScalarDraft;
  criticalPressurePa: Pr76ScalarDraft;
  acentricFactor: Pr76ScalarDraft;
}

export interface Pr76PairDraft {
  readonly key: string;
  readonly firstComponentKey: string;
  readonly secondComponentKey: string;
  kij: Pr76ScalarDraft;
  readonly original?: Pr76BinaryInteraction;
}

export type Pr76SettingsGroup =
  | 'eosRoot'
  | 'initialStability'
  | 'twoPhase'
  | 'finalTwoPhaseStability'
  | 'threePhase'
  | 'finalThreePhaseStability';

export interface Pr76ExpertDraft {
  displayName: string;
  datasetId: string;
  revision: string;
  userRecord: ExpertUserRecordDraft;
  components: Pr76ComponentDraft[];
  pairs: Pr76PairDraft[];
  readonly originalApplicability?: ModelApplicability;
  solverMode: 'preset' | 'custom';
  presetId: string;
  /** Complete resolved settings are available after deriving from a live snapshot. */
  customSettings?: PtSolverSettings;
}

export interface Pr76DraftIssue {
  readonly field: string;
  readonly message: string;
}

export class Pr76DraftValidationError extends Error {
  constructor(readonly issues: readonly Pr76DraftIssue[]) {
    super(`PR76 Expert draft has ${issues.length} validation issue${issues.length === 1 ? '' : 's'}`);
    this.name = 'Pr76DraftValidationError';
  }
}

function key(): string {
  return globalThis.crypto.randomUUID();
}

function scalarDraft(
  canonicalUnit: Pr76ScalarDraft['canonicalUnit'],
  positive: boolean,
  optional: boolean,
  original?: ModelScalar,
): Pr76ScalarDraft {
  return {
    text: original?.value === undefined ? '' : String(original.value),
    dirty: original === undefined,
    canonicalUnit,
    positive,
    optional,
    ...(original === undefined ? {} : { original: clone(ModelScalarSchema, original) }),
  };
}

function pairKey(a: string, b: string): string {
  return a < b ? `${a}\u001f${b}` : `${b}\u001f${a}`;
}

function syncPairs(
  components: readonly Pr76ComponentDraft[],
  prior: readonly Pr76PairDraft[],
): Pr76PairDraft[] {
  const existing = new Map(prior.map((pair) => [pair.key, pair]));
  const next: Pr76PairDraft[] = [];
  for (let i = 0; i < components.length; ++i) {
    for (let j = i + 1; j < components.length; ++j) {
      const first = components[i]!;
      const second = components[j]!;
      const id = pairKey(first.key, second.key);
      const old = existing.get(id);
      next.push(old ?? {
        key: id,
        firstComponentKey: first.key,
        secondComponentKey: second.key,
        kij: scalarDraft('dimensionless', false, false),
      });
    }
  }
  return next;
}

function blankUserRecord(): ExpertUserRecordDraft {
  return { reference: '', revision: '', locator: '', note: '', acquisition: '', usageTerms: '' };
}

export function newPr76ExpertDraft(): Pr76ExpertDraft {
  return {
    displayName: '',
    datasetId: '',
    revision: '',
    userRecord: blankUserRecord(),
    components: [],
    pairs: [],
    solverMode: 'preset',
    presetId: PR76_BALANCED_PRESET_ID,
  };
}

function requirePr76Snapshot(snapshot: ModelSnapshot) {
  const definition = snapshot.definition;
  if (!definition || definition.family !== ModelFamily.PR76 || definition.parameters.case !== 'pr76' ||
      !snapshot.settings) {
    throw new TypeError('Expert PR76 draft requires a complete PR76 model snapshot');
  }
  return { definition, pr76: definition.parameters.value, settings: snapshot.settings };
}

export function pr76ExpertDraftFromSnapshot(snapshot: ModelSnapshot): Pr76ExpertDraft {
  const { definition, pr76, settings } = requirePr76Snapshot(snapshot);
  const pure = new Map(pr76.pure.map((record) => [record.componentId, record]));
  const components: Pr76ComponentDraft[] = definition.components.map((component, index) => {
    const record = pure.get(component.componentId);
    if (!record) throw new TypeError(`PR76 snapshot missing pure record for component ${component.componentId}`);
    return {
      key: `snapshot-${index}`,
      componentId: component.componentId,
      displayName: component.displayName,
      kind: component.kind,
      definitionDirty: false,
      originalDefinition: clone(ComponentDefinitionSchema, component),
      molarMassKgPerMol: scalarDraft('kg/mol', true, true, component.molarMassKgPerMol),
      criticalTemperatureK: scalarDraft('K', true, false, record.criticalTemperatureK),
      criticalPressurePa: scalarDraft('Pa', true, false, record.criticalPressurePa),
      acentricFactor: scalarDraft('dimensionless', false, false, record.acentricFactor),
    };
  });
  const index = new Map(components.map((component) => [component.componentId, component.key]));
  const pairs: Pr76PairDraft[] = pr76.binary.map((record) => {
    const first = index.get(record.firstComponentId);
    const second = index.get(record.secondComponentId);
    if (!first || !second) throw new TypeError('PR76 snapshot binary pair references an unknown component');
    return {
      key: pairKey(first, second),
      firstComponentKey: first,
      secondComponentKey: second,
      kij: scalarDraft('dimensionless', false, false, record.kij),
      original: clone(Pr76BinaryInteractionSchema, record),
    };
  });
  return {
    displayName: definition.displayName,
    // A derived immutable model must receive explicit new identity before creation.
    datasetId: '',
    revision: '',
    userRecord: blankUserRecord(),
    components,
    pairs: syncPairs(components, pairs),
    ...(definition.applicability === undefined
      ? {}
      : { originalApplicability: clone(ModelApplicabilitySchema, definition.applicability) }),
    solverMode: settings.kind === SolverSettingsKind.CUSTOM ? 'custom' : 'preset',
    presetId: settings.presetId || PR76_BALANCED_PRESET_ID,
    customSettings: clone(PtSolverSettingsSchema, settings),
  };
}

export function addPr76Component(draft: Pr76ExpertDraft): Pr76ExpertDraft {
  const component: Pr76ComponentDraft = {
    key: key(), componentId: '', displayName: '', kind: ComponentKind.PURE,
    definitionDirty: true,
    molarMassKgPerMol: scalarDraft('kg/mol', true, true),
    criticalTemperatureK: scalarDraft('K', true, false),
    criticalPressurePa: scalarDraft('Pa', true, false),
    acentricFactor: scalarDraft('dimensionless', false, false),
  };
  const components = [...draft.components, component];
  return { ...draft, components, pairs: syncPairs(components, draft.pairs) };
}

export function removePr76Component(draft: Pr76ExpertDraft, componentKey: string): Pr76ExpertDraft {
  const components = draft.components.filter((component) => component.key !== componentKey);
  if (components.length === draft.components.length) return draft;
  return { ...draft, components, pairs: syncPairs(components, draft.pairs) };
}

export function movePr76Component(
  draft: Pr76ExpertDraft,
  componentKey: string,
  offset: -1 | 1,
): Pr76ExpertDraft {
  const index = draft.components.findIndex((component) => component.key === componentKey);
  const target = index + offset;
  if (index < 0 || target < 0 || target >= draft.components.length) return draft;
  const components = [...draft.components];
  [components[index], components[target]] = [components[target]!, components[index]!];
  return { ...draft, components, pairs: syncPairs(components, draft.pairs) };
}

export type Pr76ComponentTextField = 'componentId' | 'displayName';
export function editPr76ComponentText(
  draft: Pr76ExpertDraft,
  componentKey: string,
  field: Pr76ComponentTextField,
  value: string,
): Pr76ExpertDraft {
  const components = draft.components.map((component) => {
    if (component.key !== componentKey) return component;
    if (field === 'componentId' && component.originalDefinition) {
      throw new TypeError('Existing component identity is immutable; remove/add to replace it');
    }
    return { ...component, [field]: value, definitionDirty: true };
  });
  return { ...draft, components, pairs: syncPairs(components, draft.pairs) };
}

export function editPr76ComponentKind(
  draft: Pr76ExpertDraft,
  componentKey: string,
  kind: ComponentKind,
): Pr76ExpertDraft {
  return {
    ...draft,
    components: draft.components.map((component) => component.key === componentKey
      ? { ...component, kind, definitionDirty: true }
      : component),
  };
}

export type Pr76ComponentScalarField =
  | 'molarMassKgPerMol'
  | 'criticalTemperatureK'
  | 'criticalPressurePa'
  | 'acentricFactor';
export function editPr76ComponentScalar(
  draft: Pr76ExpertDraft,
  componentKey: string,
  field: Pr76ComponentScalarField,
  text: string,
): Pr76ExpertDraft {
  return {
    ...draft,
    components: draft.components.map((component) => component.key === componentKey
      ? { ...component, [field]: { ...component[field], text, dirty: true } }
      : component),
  };
}

export function editPr76Kij(
  draft: Pr76ExpertDraft,
  pairId: string,
  text: string,
): Pr76ExpertDraft {
  return {
    ...draft,
    pairs: draft.pairs.map((pair) => pair.key === pairId
      ? { ...pair, kij: { ...pair.kij, text, dirty: true } }
      : pair),
  };
}

export function setPr76SolverMode(
  draft: Pr76ExpertDraft,
  mode: 'preset' | 'custom',
): Pr76ExpertDraft {
  if (mode === 'custom' && !draft.customSettings) {
    throw new TypeError('Custom settings require a resolved settings snapshot first');
  }
  return { ...draft, solverMode: mode };
}

export function editPr76Setting(
  draft: Pr76ExpertDraft,
  group: Pr76SettingsGroup,
  field: string,
  value: string | boolean,
): Pr76ExpertDraft {
  if (!draft.customSettings) {
    throw new TypeError('Custom settings require a resolved settings snapshot first');
  }
  const settings = clone(PtSolverSettingsSchema, draft.customSettings);
  const block = settings[group];
  if (!block) throw new TypeError(`Missing resolved settings group ${group}`);
  const record = block as unknown as Record<string, unknown>;
  const current = record[field];
  if (typeof current === 'boolean') {
    if (typeof value !== 'boolean') throw new TypeError(`Setting ${group}.${field} requires boolean input`);
    record[field] = value;
  } else if (typeof current === 'bigint') {
    if (typeof value !== 'string' || !/^-?(0|[1-9][0-9]*)$/u.test(value.trim())) {
      throw new TypeError(`Setting ${group}.${field} requires integer input`);
    }
    record[field] = BigInt(value.trim());
  } else if (typeof current === 'number') {
    if (typeof value !== 'string' || value.trim() === '' || !Number.isFinite(Number(value))) {
      throw new TypeError(`Setting ${group}.${field} requires finite numeric input`);
    }
    record[field] = Number(value);
  } else {
    throw new TypeError(`Setting ${group}.${field} is not an editable scalar`);
  }
  settings.kind = SolverSettingsKind.CUSTOM;
  settings.presetId = '';
  return { ...draft, solverMode: 'custom', customSettings: settings };
}

function nonblank(value: string): boolean {
  return /\S/u.test(value) && !value.includes('\0');
}
function visibleAscii(value: string): boolean {
  return value.length > 0 && [...value].every((character) => {
    const code = character.charCodeAt(0);
    return code >= 33 && code <= 126;
  });
}
function numeric(scalar: Pr76ScalarDraft): number | undefined {
  if (scalar.optional && scalar.text.trim() === '') return undefined;
  if (scalar.text.trim() === '') return Number.NaN;
  const value = Number(scalar.text);
  if (!Number.isFinite(value) || (scalar.positive && value <= 0)) return Number.NaN;
  return value;
}

export function validatePr76ExpertDraft(draft: Pr76ExpertDraft): Pr76DraftIssue[] {
  const issues: Pr76DraftIssue[] = [];
  const requireText = (value: string, field: string) => {
    if (!nonblank(value)) issues.push({ field, message: 'nonblank text is required' });
  };
  requireText(draft.displayName, 'display_name');
  requireText(draft.datasetId, 'dataset_id');
  requireText(draft.revision, 'revision');
  for (const field of ['reference', 'revision', 'locator', 'acquisition', 'usageTerms'] as const) {
    requireText(draft.userRecord[field], `provenance.${field}`);
  }
  if (draft.components.length === 0) issues.push({ field: 'components', message: 'at least one component is required' });
  const seen = new Set<string>();
  for (let index = 0; index < draft.components.length; ++index) {
    const component = draft.components[index]!;
    const prefix = `components[${index}]`;
    if (!visibleAscii(component.componentId)) {
      issues.push({ field: `${prefix}.component_id`, message: 'visible ASCII component ID is required' });
    } else if (seen.has(component.componentId)) {
      issues.push({ field: `${prefix}.component_id`, message: 'component ID must be unique' });
    } else seen.add(component.componentId);
    requireText(component.displayName, `${prefix}.display_name`);
    if (component.kind !== ComponentKind.PURE && component.kind !== ComponentKind.PSEUDO) {
      issues.push({ field: `${prefix}.kind`, message: 'component kind must be pure or pseudo' });
    }
    for (const [name, scalar] of [
      ['molar_mass_kg_per_mol', component.molarMassKgPerMol],
      ['critical_temperature_k', component.criticalTemperatureK],
      ['critical_pressure_pa', component.criticalPressurePa],
      ['acentric_factor', component.acentricFactor],
    ] as const) {
      const value = numeric(scalar);
      if (value !== undefined && !Number.isFinite(value)) {
        issues.push({ field: `${prefix}.${name}`, message: scalar.positive ? 'finite positive value is required' : 'finite value is required' });
      }
    }
  }
  const expectedPairs = draft.components.length * (draft.components.length - 1) / 2;
  if (draft.pairs.length !== expectedPairs) {
    issues.push({ field: 'parameters.binary', message: 'every unordered component pair requires an explicit kij record' });
  }
  const known = new Set(draft.components.map((component) => component.key));
  for (let index = 0; index < draft.pairs.length; ++index) {
    const pair = draft.pairs[index]!;
    if (!known.has(pair.firstComponentKey) || !known.has(pair.secondComponentKey)) {
      issues.push({ field: `parameters.binary[${index}]`, message: 'pair references a removed component' });
    }
    const value = numeric(pair.kij);
    if (value === undefined || !Number.isFinite(value)) {
      issues.push({ field: `parameters.binary[${index}].kij`, message: 'explicit finite kij is required; missing is never zero' });
    }
  }
  if (draft.solverMode === 'preset') {
    requireText(draft.presetId, 'solver_selection.preset_id');
  } else if (!draft.customSettings) {
    issues.push({ field: 'solver_selection.settings', message: 'complete resolved settings are required for custom mode' });
  }
  return issues;
}

function userProvenance(draft: ExpertUserRecordDraft): ModelProvenance {
  return create(ModelProvenanceSchema, {
    kind: SourceKind.USER_SUPPLIED,
    reference: draft.reference,
    revision: draft.revision,
    locator: draft.locator,
    note: draft.note,
    acquisition: draft.acquisition,
    usageTerms: draft.usageTerms,
  });
}

function scalar(
  draft: Pr76ScalarDraft,
  provenance: ModelProvenance,
): ModelScalar | undefined {
  const value = numeric(draft);
  if (value === undefined) return undefined;
  if (!draft.dirty && draft.original) return clone(ModelScalarSchema, draft.original);
  return create(ModelScalarSchema, {
    value,
    provenance: clone(ModelProvenanceSchema, provenance),
    originalUnit: draft.canonicalUnit,
    conversion: 'identity: Expert UI canonical SI input',
  });
}

export function buildPr76CreateInput(draft: Pr76ExpertDraft): ModelCreateInput {
  const issues = validatePr76ExpertDraft(draft);
  if (issues.length > 0) throw new Pr76DraftValidationError(issues);
  const provenance = userProvenance(draft.userRecord);
  const componentByKey = new Map(draft.components.map((component) => [component.key, component]));
  const components = draft.components.map((component) => {
    const prior = component.originalDefinition;
    return create(ComponentDefinitionSchema, {
      componentId: component.componentId,
      displayName: component.displayName,
      kind: component.kind,
      provenance: !component.definitionDirty && prior?.provenance
        ? clone(ModelProvenanceSchema, prior.provenance)
        : clone(ModelProvenanceSchema, provenance),
      ...(scalar(component.molarMassKgPerMol, provenance) === undefined
        ? {}
        : { molarMassKgPerMol: scalar(component.molarMassKgPerMol, provenance)! }),
    });
  });
  const pure = draft.components.map((component) => create(Pr76PureParametersSchema, {
    componentId: component.componentId,
    criticalTemperatureK: scalar(component.criticalTemperatureK, provenance),
    criticalPressurePa: scalar(component.criticalPressurePa, provenance),
    acentricFactor: scalar(component.acentricFactor, provenance),
  }));
  const binary = draft.pairs.map((pair) => {
    const first = componentByKey.get(pair.firstComponentKey)!;
    const second = componentByKey.get(pair.secondComponentKey)!;
    return create(Pr76BinaryInteractionSchema, {
      firstComponentId: first.componentId,
      secondComponentId: second.componentId,
      kij: scalar(pair.kij, provenance),
    });
  });
  const applicability = draft.originalApplicability
    ? clone(ModelApplicabilitySchema, draft.originalApplicability)
    : create(ModelApplicabilitySchema, { provenance: clone(ModelProvenanceSchema, provenance) });
  const definition = create(ThermodynamicModelDefinitionSchema, {
    version: PR76_MODEL_DEFINITION_VERSION,
    family: ModelFamily.PR76,
    displayName: draft.displayName,
    datasetId: draft.datasetId,
    revision: draft.revision,
    provenance: clone(ModelProvenanceSchema, provenance),
    components,
    applicability,
    parameters: {
      case: 'pr76',
      value: create(Pr76ParameterDefinitionSchema, { pure, binary }),
    },
  });
  if (draft.solverMode === 'preset') {
    return { definition, solverSelection: { case: 'presetId', value: draft.presetId } };
  }
  const settings = clone(PtSolverSettingsSchema, draft.customSettings!);
  settings.kind = SolverSettingsKind.CUSTOM;
  settings.presetId = '';
  return { definition, solverSelection: { case: 'settings', value: settings } };
}
