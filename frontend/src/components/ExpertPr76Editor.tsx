import {
  useEffect,
  useMemo,
  useRef,
  useState,
  type FormEvent,
  type SetStateAction,
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

/** A late create belongs to neither an unmounted editor nor a replacement view. */
export async function handoffExpertModel(
  model: ExpertOwnedModel,
  isCurrent: () => boolean,
  accept: (model: ExpertOwnedModel) => void,
): Promise<boolean> {
  if (isCurrent()) { accept(model); return true; }
  try { await model.release(); } catch { /* Existing typed session owns ambiguous cleanup. */ }
  return false;
}

export interface ExpertPr76EditorProps {
  owner: ExpertModelOwner;
  seedSnapshot?: ModelSnapshot;
  onCreated(model: ExpertOwnedModel): void;
  onDraftChanged?(): void;
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
          <details className="contract-details" key={group}>
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
 * PR76 create/edit UI. It edits only a local draft and every submit creates a
 * new immutable model through ExpertModelOwner. Raw model handles never enter
 * component state or the DOM.
 */
export function ExpertPr76Editor({ owner, seedSnapshot, onCreated, onDraftChanged }: ExpertPr76EditorProps) {
  const [draft, setDraftValue] = useState(() => initialDraft(seedSnapshot));
  const [submitted, setSubmitted] = useState(false);
  const [busy, setBusy] = useState(false);
  const [failure, setFailure] = useState<ExpertPr76CreateFailure | null>(null);
  const alive = useRef(true);
  const draftRef = useRef(draft);
  const creating = useRef(false);

  useEffect(() => {
    alive.current = true;
    return () => { alive.current = false; };
  }, []);

  function setDraft(update: SetStateAction<Pr76ExpertDraft>) {
    if (creating.current) return;
    const next = typeof update === 'function' ? update(draftRef.current) : update;
    draftRef.current = next;
    setDraftValue(next);
    onDraftChanged?.();
  }

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
    if (creating.current || validatePr76ExpertDraft(draftRef.current).length > 0) return;
    creating.current = true;
    onDraftChanged?.();
    setBusy(true);
    try {
      const made = await createPr76ExpertModel(owner, draftRef.current);
      if (!await handoffExpertModel(made, () => alive.current, onCreated)) return;
      if (!alive.current) return;
      try {
        const nextDraft = pr76ExpertDraftFromSnapshot(made.snapshot);
        draftRef.current = nextDraft;
        setDraftValue(nextDraft);
        setSubmitted(false);
      } catch {
        setFailure({ reason: 'expert.created_snapshot_invalid' });
      }
    } catch (cause) {
      if (alive.current) setFailure(expertPr76CreateFailure(cause));
    } finally {
      creating.current = false;
      if (alive.current) setBusy(false);
    }
  }

  return (
    <form
      className="flash-form expert-pr76-editor pr-fluid-editor"
      onSubmit={(event) => void submit(event)}
      onInput={() => onDraftChanged?.()}
      noValidate
    >
      <fieldset disabled={busy} style={{ border: 0, padding: 0, margin: 0, minWidth: 0 }}>
        <div className="section-heading">
          <div>
            <p className="eyebrow">Fluid definition</p>
            <h2>{seedSnapshot ? 'Edit PR fluid' : 'Define PR fluid'}</h2>
          </div>
          <span className="model-chip">Classic PR</span>
        </div>

        <label className="backend-selector pr-model-name">
          <span>Fluid / model name</span>
          <input
            value={draft.displayName}
            onChange={(event) => setIdentity('displayName', event.currentTarget.value)}
          />
        </label>

        <details className="contract-details advanced-details pr-record-details">
          <summary>Data source and model record (required)</summary>
          <p className="submit-note">
            Keep a traceable identifier and source for user-supplied component and interaction data.
          </p>
          <div className="field-grid field-grid-two">
            <label>
              <span>Dataset ID</span>
              <input
                value={draft.datasetId}
                onChange={(event) => setIdentity('datasetId', event.currentTarget.value)}
              />
            </label>
            <label>
              <span>Model revision</span>
              <input
                value={draft.revision}
                onChange={(event) => setIdentity('revision', event.currentTarget.value)}
              />
            </label>
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

        <div className="composition-header pr-editor-section">
          <div>
            <h3>Components</h3>
            <p>Add, remove or reorder the components used by this PR fluid.</p>
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
          <details className="contract-details pr-component-card" open key={component.key} data-component-index={index}>
            <summary>{component.displayName || component.componentId || `Component ${index + 1}`}</summary>
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
                <span>Component type</span>
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
            <div className="submit-actions pr-component-actions">
              <button
                type="button"
                className="secondary-button"
                disabled={busy || index === 0}
                onClick={() => setDraft((current) => movePr76Component(current, component.key, -1))}
              >
                Move up
              </button>
              <button
                type="button"
                className="secondary-button"
                disabled={busy || index === draft.components.length - 1}
                onClick={() => setDraft((current) => movePr76Component(current, component.key, 1))}
              >
                Move down
              </button>
              <button
                type="button"
                className="secondary-button"
                disabled={busy}
                onClick={() => setDraft((current) => removePr76Component(current, component.key))}
              >
                Remove
              </button>
            </div>
          </details>
        ))}

        <div className="composition-header pr-editor-section">
          <div>
            <h3>Binary interaction coefficients</h3>
            <p>Enter an explicit finite kij for every component pair.</p>
          </div>
        </div>
        {draft.pairs.length === 0 ? (
          <div className="model-contract model-contract-empty">
            <span>Add at least two components to define a binary interaction coefficient.</span>
          </div>
        ) : draft.pairs.map((pair) => {
          const first = componentsByKey.get(pair.firstComponentKey);
          const second = componentsByKey.get(pair.secondComponentKey);
          return (
            <label className="backend-selector pr-kij-row" key={pair.key}>
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

        <details className="contract-details advanced-details pr-solver-settings">
          <summary>Advanced: numerical solver settings</summary>
          <label className="backend-selector">
            <span>Settings mode</span>
            <select
              value={draft.solverMode}
              disabled={busy}
              onChange={(event) => setDraft((current) =>
                setPr76SolverMode(current, event.currentTarget.value as 'preset' | 'custom'))}
            >
              <option value="preset">Validated preset</option>
              <option value="custom" disabled={!draft.customSettings}>Custom settings</option>
            </select>
          </label>
          {draft.solverMode === 'preset' ? (
            <div className="model-contract">
              <span>Preset</span>
              <code>{draft.presetId}</code>
              {!draft.customSettings ? (
                <span>Apply once with the preset to obtain the complete editable settings snapshot.</span>
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
            <strong>Complete the required fluid data before applying.</strong>
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
            <strong>The PR fluid model could not be applied.</strong>
            <p>Check the entered model data and try again.</p>
            <details className="inline-details">
              <summary>Technical details</summary>
              <code>{failure.reason}</code>
              {failure.code === undefined ? null : <code> · {failure.code}</code>}
              {failure.validation?.field === undefined ? null : (
                <code> · {failure.validation.field}</code>
              )}
            </details>
          </div>
        ) : null}

        <div className="submit-row">
          <span className="submit-note">
            Applying creates a new immutable PR76 model from the entered data.
          </span>
          <button className="primary-button" type="submit" disabled={busy}>
            {busy ? 'Applying…' : 'Apply fluid model'}
          </button>
        </div>
      </fieldset>
    </form>
  );
}
