import { useEffect, useMemo, useState, type FormEvent } from 'react';

import {
  PT_COMPOSITION_ROUNDOFF,
  validatePtFlashRequest,
  type PtBackendDescriptor,
  type PtFlashRequest,
} from '../domain/flash';

interface FlashFormProps {
  backends: readonly PtBackendDescriptor[];
  discoveryLoading: boolean;
  busy: boolean;
  onBackendChange: () => void;
  onCancel: () => void;
  onSubmit: (request: PtFlashRequest, backend: PtBackendDescriptor) => void;
}

type FeedDrafts = ReadonlyMap<string, ReadonlyMap<string, string>>;

function parseFiniteInput(value: string): number {
  return value.trim() === '' ? Number.NaN : Number(value);
}

export function FlashForm({
  backends,
  discoveryLoading,
  busy,
  onBackendChange,
  onCancel,
  onSubmit,
}: FlashFormProps) {
  const [selectedBackendId, setSelectedBackendId] = useState(
    () => backends[0]?.configuredBackendId ?? '',
  );
  const [pressureMpa, setPressureMpa] = useState('10');
  const [temperatureK, setTemperatureK] = useState('350');
  const [feedDrafts, setFeedDrafts] = useState<FeedDrafts>(() => new Map());
  const [submitted, setSubmitted] = useState(false);

  useEffect(() => {
    setSelectedBackendId((current) => {
      if (backends.some((backend) => backend.configuredBackendId === current)) {
        return current;
      }
      return backends[0]?.configuredBackendId ?? '';
    });
  }, [backends]);

  const selectedBackend = useMemo(
    () =>
      backends.find(
        (backend) => backend.configuredBackendId === selectedBackendId,
      ) ?? null,
    [backends, selectedBackendId],
  );

  const request = useMemo<PtFlashRequest | null>(() => {
    if (selectedBackend === null) {
      return null;
    }
    const drafts = feedDrafts.get(selectedBackend.configuredBackendId);
    return {
      configuredBackendId: selectedBackend.configuredBackendId,
      pressurePa: parseFiniteInput(pressureMpa) * 1e6,
      temperatureK: parseFiniteInput(temperatureK),
      feed: selectedBackend.componentInventory.components.map((component) => ({
        componentId: component.componentId,
        moleFraction: parseFiniteInput(drafts?.get(component.componentId) ?? ''),
      })),
    };
  }, [feedDrafts, pressureMpa, selectedBackend, temperatureK]);

  const validationErrors = useMemo(
    () =>
      request === null || selectedBackend === null
        ? []
        : validatePtFlashRequest(request, selectedBackend),
    [request, selectedBackend],
  );

  const compositionSum = useMemo(() => {
    if (selectedBackend === null) {
      return 0;
    }
    const drafts = feedDrafts.get(selectedBackend.configuredBackendId);
    return selectedBackend.componentInventory.components.reduce((sum, component) => {
      const value = parseFiniteInput(drafts?.get(component.componentId) ?? '');
      return Number.isFinite(value) ? sum + value : sum;
    }, 0);
  }, [feedDrafts, selectedBackend]);

  function selectBackend(configuredBackendId: string) {
    setSelectedBackendId(configuredBackendId);
    setSubmitted(false);
    onBackendChange();
  }

  function updateMoleFraction(componentId: string, value: string) {
    if (selectedBackend === null) {
      return;
    }
    const backendId = selectedBackend.configuredBackendId;
    setFeedDrafts((current) => {
      const next = new Map(current);
      const backendDrafts = new Map(current.get(backendId));
      backendDrafts.set(componentId, value);
      next.set(backendId, backendDrafts);
      return next;
    });
  }

  function handleSubmit(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    setSubmitted(true);
    if (
      request === null ||
      selectedBackend === null ||
      validationErrors.length > 0 ||
      busy
    ) {
      return;
    }
    onSubmit(request, selectedBackend);
  }

  const capability = selectedBackend?.capability;

  return (
    <form className="flash-form" onSubmit={handleSubmit} noValidate>
      <div className="section-heading">
        <div>
          <p className="eyebrow">PT service request</p>
          <h2>Calculation input</h2>
        </div>
        <span className="model-chip">Runtime discovery</span>
      </div>

      <label className="backend-selector">
        <span>Configured backend</span>
        <select
          value={selectedBackendId}
          onChange={(event) => selectBackend(event.target.value)}
          disabled={busy || backends.length === 0}
          aria-label="Configured PT backend"
        >
          {backends.length === 0 ? (
            <option value="">
              {discoveryLoading ? 'Discovering backends…' : 'No backend available'}
            </option>
          ) : null}
          {backends.map((backend) => (
            <option
              value={backend.configuredBackendId}
              key={backend.configuredBackendId}
            >
              {backend.configuredBackendId} · {backend.capability.modelProfile}
            </option>
          ))}
        </select>
      </label>

      {capability === undefined ? (
        <div className="model-contract model-contract-empty">
          <span>
            Capability discovery must succeed before a component inventory or solve
            request can be created.
          </span>
        </div>
      ) : (
        <div className="capability-card">
          <div className="capability-primary">
            <span>Model profile</span>
            <code>{capability.modelProfile}</code>
          </div>
          <div className="capability-facts">
            <span>
              Phases <strong>{capability.supportedPhaseCounts.join(', ')}</strong>
            </span>
            <span>
              Dataset <strong>{capability.datasetId}</strong>
            </span>
            <span>
              Revision <strong>{capability.revision}</strong>
            </span>
          </div>
          {capability.scalarSettings.length > 0 ? (
            <div className="scalar-settings" aria-label="Configured backend settings">
              {capability.scalarSettings.map((setting) => (
                <span key={setting.id}>
                  {setting.id} = <strong>{setting.value}</strong> {setting.unit}
                </span>
              ))}
            </div>
          ) : null}
        </div>
      )}

      <div className="field-grid field-grid-two">
        <label>
          <span>Pressure</span>
          <div className="input-with-unit">
            <input
              inputMode="decimal"
              value={pressureMpa}
              onChange={(event) => setPressureMpa(event.target.value)}
              aria-label="Pressure in MPa"
            />
            <span>MPa</span>
          </div>
        </label>
        <label>
          <span>Temperature</span>
          <div className="input-with-unit">
            <input
              inputMode="decimal"
              value={temperatureK}
              onChange={(event) => setTemperatureK(event.target.value)}
              aria-label="Temperature in kelvin"
            />
            <span>K</span>
          </div>
        </label>
      </div>

      <div className="composition-header">
        <div>
          <h3>Backend component inventory</h3>
          <p>IDs and order come from discovery; only mole fractions are editable.</p>
        </div>
        <div
          className="composition-sum"
          data-valid={Math.abs(compositionSum - 1) <= PT_COMPOSITION_ROUNDOFF}
        >
          Σz = {compositionSum.toPrecision(8)}
        </div>
      </div>

      <div
        className="component-table inventory-component-table"
        role="group"
        aria-label="Overall composition"
      >
        <div className="component-table-head" aria-hidden="true">
          <span>Component ID</span>
          <span>Mole fraction z</span>
        </div>
        {selectedBackend?.componentInventory.components.map((component) => {
          const value =
            feedDrafts
              .get(selectedBackend.configuredBackendId)
              ?.get(component.componentId) ?? '';
          return (
            <div className="component-row" key={component.componentId}>
              <code>{component.componentId}</code>
              <input
                inputMode="decimal"
                value={value}
                onChange={(event) =>
                  updateMoleFraction(component.componentId, event.target.value)
                }
                placeholder="0.0"
                aria-label={`${component.componentId} mole fraction`}
              />
            </div>
          );
        })}
      </div>

      {submitted && validationErrors.length > 0 ? (
        <div className="validation-panel" role="alert">
          <strong>Input contract not satisfied</strong>
          <ul>
            {validationErrors.map((error) => (
              <li key={error}>{error}</li>
            ))}
          </ul>
        </div>
      ) : null}

      <div className="submit-row">
        <div className="submit-note">
          Browser checks are advisory. The service validates the unchanged request and
          delegates exactly one solve to the selected backend.
        </div>
        <div className="submit-actions">
          {busy ? (
            <button className="secondary-button cancel-button" type="button" onClick={onCancel}>
              Cancel
            </button>
          ) : null}
          <button
            className="primary-button"
            type="submit"
            disabled={selectedBackend === null || busy}
          >
            {busy ? 'Computing…' : 'Run PT flash'}
          </button>
        </div>
      </div>
    </form>
  );
}
