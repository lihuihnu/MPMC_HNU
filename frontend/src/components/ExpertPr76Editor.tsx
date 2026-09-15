import {
  useEffect,
  useMemo,
  useRef,
  useState,
  type FormEvent,
} from 'react';
import { Code } from '@connectrpc/connect';

import type { ExpertModelOwner, ExpertOwnedModel } from '../api/expertModelOwner';
import type { ModelCallOptions } from '../api/modelSessionClient';
import { ModelClientError } from '../api/modelSessionClient';
import { RendererModelError } from '../api/rendererModelClient';
import type { ModelValidationDetail } from '../api/modelValidationDetail';
import {
  addPr76Component,
  buildPr76CreateInput,
  editPr76ComponentKind,
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
  type Pr76ComponentScalarField,
  type Pr76ExpertDraft,
  type Pr76SettingsGroup,
} from '../api/pr76ExpertDraft';
import {
  ComponentKind,
  type ModelSnapshot,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';

export interface ExpertPr76CreateFailure {
  readonly reason: string;
  readonly code?: Code;
  readonly validation?: ModelValidationDetail;
}

export function expertPr76CreateFailure(cause: unknown): ExpertPr76CreateFailure {
  if (cause instanceof ModelClientError || cause instanceof RendererModelError) {
    return Object.freeze({
      reason: cause.reason,
      code: cause.code,
      ...(cause.validation === undefined ? {} : { validation: cause.validation }),
    });
  }
  return Object.freeze({ reason: 'expert.create_failed' });
}

export function createPr76ExpertModel(
  owner: ExpertModelOwner,
  draft: Pr76ExpertDraft,
  options: ModelCallOptions = {},
): Promise<ExpertOwnedModel> {
  return owner.create(buildPr76CreateInput(draft), options);
}

export interface ExpertPr76EditorProps {
  owner: ExpertModelOwner;
  seedSnapshot?: ModelSnapshot;
  onCreated(model: ExpertOwnedModel): void;
}

const settingsGroups: ReadonlyArray<[Pr76SettingsGroup, string]> = [
  ['eosRoot', 'EOS root'],
  ['initialStability', 'Initial stability'],
  ['twoPhase', 'Two phase'],
  ['finalTwoPhaseStability', 'Final two-phase stability'],
  ['threePhase', 'Three phase'],
  ['finalThreePhaseStability', 'Final three-phase stability'],
];

function initialDraft(seedSnapshot?: ModelSnapshot): Pr76ExpertDraft {
  return seedSnapshot ? pr76ExpertDraftFromSnapshot(seedSnapshot) : newPr76ExpertDraft();
}

function scalarLabel(field: Pr76ComponentScalarField): [string, string] {
  switch (field) {
    case 'molarMassKgPerMol': return ['Molar mass', 'kg/mol'];
    case 'criticalTemperatureK': return ['Critical temperature', 'K'];
    case 'criticalPressurePa': return ['Critical pressure', 'Pa'];
    case 'acentricFactor': return ['Acentric factor', ''];
  }
}

function SettingsEditor({
  draft,
  onEdit,
}: {
  draft: Pr76ExpertDraft;
  onEdit(group: Pr76SettingsGroup, field: string, value: string | boolean): void;
}) {
  const settings = draft.customSettings;
  if (draft.solverMode !== 'custom' || !settings) return null;
  return (
    <div aria-label="Custom solver settings">
      {settingsGroups.map(([group, title]) => {
        const block = settings[group];
        const entries = block
          ? Object.entries(block as unknown as Record<string, unknown>).filter(
              ([name, value]) => name !== '$typeName' && value !== undefined,
            )
          : [];
        return (
          <details className="contract-details" open key={group}>
            <summary>{title}</summary>
            <div className="field-grid field-grid-two">
              {entries.map(([field, value]) => (
                <label key={`${group}.${field}`}>
                  <span>{field}</span>
                  {typeof value === 'boolean' ? (
                    <input
                      type="checkbox"
                      checked={value}
                      onChange={(event) => onEdit(group, field, event.currentTarget.checked)}
                    />
                  ) : (
                    <input
                      defaultValue={typeof value === 'bigint' ? value.toString() : String(value)}
                      inputMode="decimal"
                      onBlur={(event) => {
                        try {
                          onEdit(group, field, event.currentTarget.value);
                        } catch {
                          event.currentTarget.value = typeof value === 'bigint'
                            ? value.toString()
                            : String(value);
                        }
                      }}
                    />
                  )}
                </label>
              ))}
            </div>
          </details>
        );
      })}
    </div>
  );
}

/**
 * PR76 Expert create/edit UI. It edits only a local draft and every submit creates
 * a new immutable model through ExpertModelOwner. Raw model handles never enter
 * component state or the DOM.
 */
export function ExpertPr76Editor({ owner, seedSnapshot, onCreated }: ExpertPr76EditorProps) {
  const [draft, setDraft] = useState(() => initialDraft(seedSnapshot));
  const [submitted, setSubmitted] = useState(false);
  const [busy, setBusy] = useState(false);
  const [failure, setFailure] = useState<ExpertPr76CreateFailure | null>(null);
  const alive = useRef(true);

  useEffect(() => () => { alive.current = false; }, []);

  const issues = useMemo(() => validatePr76ExpertDraft(draft), [draft]);
  const componentsByKey = useMemo(
    () => new Map(draft.components.map((component) => [component.key, component])),
    [draft.components],
  );

  function setIdentity(field: 'displayName' | 'datasetId' | 'revision', value: string) {
    setDraft((current) => ({ ...current, [field]: value }));
  }

  function setUserRecord(
    field: keyof Pr76ExpertDraft['userRecord'],
    value: string,
  ) {
    setDraft((current) => ({
      ...current,
      userRecord: { ...current.userRecord, [field]: value },
    }));
  }

  async function submit(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    setSubmitted(true);
    setFailure(null);
    if (busy || issues.length > 0) return;
    setBusy(true);
    try {
      const made = await createPr76ExpertModel(owner, draft);
      if (!alive.current) return;
      let nextDraft: Pr76ExpertDraft | undefined;
      try {
        nextDraft = pr76ExpertDraftFromSnapshot(made.snapshot);
      } catch {
        setFailure({ reason: 'expert.created_snapshot_invalid' });
      }
      onCreated(made);
      if (nextDraft) {
        setDraft(nextDraft);
        setSubmitted(false);
      }
    } catch (cause) {
      if (alive.current) setFailure(expertPr76CreateFailure(cause));
    } finally {
      if (alive.current) setBusy(false);
    }
  }

  return (
    <form className="flash-form expert-pr76-editor" onSubmit={(event) => void submit(event)} noValidate>
      <div className="section-heading">
        <div>
          <p className="eyebrow">PR76 Expert model</p>
          <h2>{seedSnapshot ? 'Create immutable revision' : 'Create custom PR76 model'}</h2>
        </div>
        <span className="model-chip">Explicit parameters</span>
      </div>

      <div className="field-grid">
        {([
          ['displayName', 'Display name'],
          ['datasetId', 'New dataset ID'],
          ['revision', 'New revision'],
        ] as const).map(([field, label]) => (
          <label key={field}>
            <span>{label}</span>
            <input
              value={draft[field]}
              onChange={(event) => setIdentity(field, event.currentTarget.value)}
            />
          </label>
        ))}
      </div>

      <details className="contract-details" open>
        <summary>User-supplied provenance for new/changed records</summary>
        <div className="field-grid field-grid-two">
          {([
            ['reference', 'Reference'],
            ['revision', 'Source revision'],
            ['locator', 'Locator'],
            ['acquisition', 'Acquisition'],
            ['usageTerms', 'Usage terms'],
            ['note', 'Note (optional)'],
          ] as const).map(([field, label]) => (
            <label key={field}>
              <span>{label}</span>
              <input
                value={draft.userRecord[field]}
                onChange={(event) => setUserRecord(field, event.currentTarget.value)}
              />
            </label>
          ))}
        </div>
      </details>

      <div className="composition-header">
        <div>
          <h3>Ordered components</h3>
          <p>Replacing an existing physical identity requires remove + add.</p>
        </div>
        <button
          className="secondary-button cancel-button"
          type="button"
          disabled={busy}
          onClick={() => setDraft((current) => addPr76Component(current))}
        >
          Add component
        </button>
      </div>

      {draft.components.map((component, index) => (
        <details className="contract-details" open key={component.key} data-component-index={index}>
          <summary>{index}: {component.componentId || 'new component'}</summary>
          <div className="field-grid field-grid-two">
            <label>
              <span>Component ID</span>
              <input
                value={component.componentId}
                readOnly={component.originalDefinition !== undefined}
                onChange={(event) => setDraft((current) =>
                  editPr76ComponentText(current, component.key, 'componentId', event.currentTarget.value))}
              />
            </label>
            <label>
              <span>Display name</span>
              <input
                value={component.displayName}
                onChange={(event) => setDraft((current) =>
                  editPr76ComponentText(current, component.key, 'displayName', event.currentTarget.value))}
              />
            </label>
            <label>
              <span>Kind</span>
              <select
                value={component.kind}
                onChange={(event) => setDraft((current) =>
                  editPr76ComponentKind(current, component.key, Number(event.currentTarget.value) as ComponentKind))}
              >
                <option value={ComponentKind.PURE}>Pure</option>
                <option value={ComponentKind.PSEUDO}>Pseudo</option>
              </select>
            </label>
            {([
              'molarMassKgPerMol',
              'criticalTemperatureK',
              'criticalPressurePa',
              'acentricFactor',
            ] as const).map((field) => {
              const [label, unit] = scalarLabel(field);
              return (
                <label key={field}>
                  <span>{label}</span>
                  <div className="input-with-unit">
                    <input
                      inputMode="decimal"
                      value={component[field].text}
                      onChange={(event) => setDraft((current) =>
                        editPr76ComponentScalar(current, component.key, field, event.currentTarget.value))}
                    />
                    {unit ? <span>{unit}</span> : null}
                  </div>
                </label>
              );
            })}
          </div>
          <div className="submit-actions">
            <button type="button" className="secondary-button" disabled={busy || index === 0}
              onClick={() => setDraft((current) => movePr76Component(current, component.key, -1))}>
              Move up
            </button>
            <button type="button" className="secondary-button" disabled={busy || index === draft.components.length - 1}
              onClick={() => setDraft((current) => movePr76Component(current, component.key, 1))}>
              Move down
            </button>
            <button type="button" className="secondary-button" disabled={busy}
              onClick={() => setDraft((current) => removePr76Component(current, component.key))}>
              Remove
            </button>
          </div>
        </details>
      ))}

      <div className="composition-header">
        <div>
          <h3>Binary interactions</h3>
          <p>Every unordered selected pair requires an explicit finite kij; blank never means zero.</p>
        </div>
      </div>
      {draft.pairs.length === 0 ? (
        <div className="model-contract model-contract-empty">
          <span>No unordered component pair exists yet.</span>
        </div>
      ) : draft.pairs.map((pair) => {
        const first = componentsByKey.get(pair.firstComponentKey);
        const second = componentsByKey.get(pair.secondComponentKey);
        return (
          <label className="backend-selector" key={pair.key}>
            <span>kij {first?.componentId || '?'} ↔ {second?.componentId || '?'}</span>
            <input
              inputMode="decimal"
              value={pair.kij.text}
              onChange={(event) => setDraft((current) =>
                editPr76Kij(current, pair.key, event.currentTarget.value))}
            />
          </label>
        );
      })}

      <details className="contract-details" open>
        <summary>PT solver settings</summary>
        <label className="backend-selector">
          <span>Settings mode</span>
          <select
            value={draft.solverMode}
            disabled={busy}
            onChange={(event) => setDraft((current) =>
              setPr76SolverMode(current, event.currentTarget.value as 'preset' | 'custom'))}
          >
            <option value="preset">Frozen preset</option>
            <option value="custom" disabled={!draft.customSettings}>Complete custom snapshot</option>
          </select>
        </label>
        {draft.solverMode === 'preset' ? (
          <div className="model-contract">
            <span>Preset ID</span>
            <code>{draft.presetId}</code>
            {!draft.customSettings ? (
              <span>Create once with the preset to obtain a complete resolved snapshot before custom editing.</span>
            ) : null}
          </div>
        ) : (
          <SettingsEditor
            draft={draft}
            onEdit={(group, field, value) => setDraft((current) => editPr76Setting(current, group, field, value))}
          />
        )}
      </details>

      {submitted && issues.length > 0 ? (
        <div className="validation-panel" role="alert">
          <strong>Draft is incomplete</strong>
          <ul>
            {issues.map((issue, index) => (
              <li key={`${issue.field}-${index}`}>
                <code>{issue.field}</code>: {issue.message}
              </li>
            ))}
          </ul>
        </div>
      ) : null}

      {failure ? (
        <div className="validation-panel" role="alert" data-create-reason={failure.reason}>
          <strong>Model creation failed</strong>
          <p>
            Reason <code>{failure.reason}</code>
            {failure.code === undefined ? null : <> · gRPC code <code>{failure.code}</code></>}
          </p>
          {failure.validation ? (
            <p>
              <code>{failure.validation.code}</code> · field{' '}
              <code>{failure.validation.field ?? 'not supplied'}</code>
            </p>
          ) : null}
        </div>
      ) : null}

      <div className="submit-row">
        <span className="submit-note">
          Submit creates a new immutable PR76 model. No existing model handle is mutated.
        </span>
        <button className="primary-button" type="submit" disabled={busy}>
          {busy ? 'Creating…' : 'Create immutable model'}
        </button>
      </div>
    </form>
  );
}
