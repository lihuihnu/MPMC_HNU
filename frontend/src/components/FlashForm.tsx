import { useMemo, useRef, useState, type FormEvent } from 'react';

import {
  PROFILE_C_MODEL,
  type FlashComponentInput,
  type ProfileCPtRequest,
  validateProfileCPtRequest,
} from '../domain/flash';

interface FlashFormProps {
  backendConfigured: boolean;
  busy: boolean;
  onSubmit: (request: ProfileCPtRequest) => void;
}

interface ComponentDraft {
  key: number;
  id: string;
  moleFraction: string;
}

function parseFiniteInput(value: string): number {
  if (value.trim() === '') {
    return Number.NaN;
  }
  return Number(value);
}

export function FlashForm({ backendConfigured, busy, onSubmit }: FlashFormProps) {
  const nextKey = useRef(3);
  const [pressureMpa, setPressureMpa] = useState('10');
  const [temperatureK, setTemperatureK] = useState('350');
  const [naclMolality, setNaclMolality] = useState('0');
  const [components, setComponents] = useState<ComponentDraft[]>([
    { key: 1, id: '', moleFraction: '' },
    { key: 2, id: '', moleFraction: '' },
  ]);
  const [submitted, setSubmitted] = useState(false);

  const request = useMemo<ProfileCPtRequest>(() => {
    const parsedComponents: FlashComponentInput[] = components.map((component) => ({
      id: component.id,
      moleFraction: parseFiniteInput(component.moleFraction),
    }));
    return {
      pressurePa: parseFiniteInput(pressureMpa) * 1e6,
      temperatureK: parseFiniteInput(temperatureK),
      naclMolalityMolPerKgWater: parseFiniteInput(naclMolality),
      components: parsedComponents,
    };
  }, [components, naclMolality, pressureMpa, temperatureK]);

  const validationErrors = useMemo(
    () => validateProfileCPtRequest(request),
    [request],
  );

  const compositionSum = useMemo(
    () =>
      components.reduce((sum, component) => {
        const value = parseFiniteInput(component.moleFraction);
        return Number.isFinite(value) ? sum + value : sum;
      }, 0),
    [components],
  );

  function updateComponent(
    key: number,
    field: 'id' | 'moleFraction',
    value: string,
  ) {
    setComponents((current) =>
      current.map((component) =>
        component.key === key ? { ...component, [field]: value } : component,
      ),
    );
  }

  function addComponent() {
    const key = nextKey.current;
    nextKey.current += 1;
    setComponents((current) => [
      ...current,
      { key, id: '', moleFraction: '' },
    ]);
  }

  function removeComponent(key: number) {
    setComponents((current) => current.filter((component) => component.key !== key));
  }

  function handleSubmit(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    setSubmitted(true);
    if (!backendConfigured || validationErrors.length > 0 || busy) {
      return;
    }
    onSubmit(request);
  }

  const disabled = !backendConfigured || validationErrors.length > 0 || busy;

  return (
    <form className="flash-form" onSubmit={handleSubmit} noValidate>
      <div className="section-heading">
        <div>
          <p className="eyebrow">Profile-C PT flash</p>
          <h2>Calculation input</h2>
        </div>
        <span className="model-chip">SW92 corrected-original</span>
      </div>

      <div className="model-contract">
        <span>Model profile</span>
        <code>{PROFILE_C_MODEL}</code>
      </div>

      <div className="field-grid">
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
        <label>
          <span>NaCl molality</span>
          <div className="input-with-unit">
            <input
              inputMode="decimal"
              value={naclMolality}
              onChange={(event) => setNaclMolality(event.target.value)}
              aria-label="NaCl molality"
            />
            <span>mol/kg H₂O</span>
          </div>
        </label>
      </div>

      <div className="composition-header">
        <div>
          <h3>Ordered overall composition</h3>
          <p>Stable component IDs must match the backend parameter snapshot.</p>
        </div>
        <div className="composition-sum" data-valid={Math.abs(compositionSum - 1) <= 64 * Number.EPSILON}>
          Σz = {compositionSum.toPrecision(8)}
        </div>
      </div>

      <div className="component-table" role="group" aria-label="Overall composition">
        <div className="component-table-head" aria-hidden="true">
          <span>Component ID</span>
          <span>Mole fraction z</span>
          <span />
        </div>
        {components.map((component, index) => (
          <div className="component-row" key={component.key}>
            <input
              value={component.id}
              onChange={(event) => updateComponent(component.key, 'id', event.target.value)}
              placeholder={index === 0 ? 'e.g. carbon-dioxide' : 'stable ID'}
              aria-label={`Component ${index + 1} ID`}
            />
            <input
              inputMode="decimal"
              value={component.moleFraction}
              onChange={(event) =>
                updateComponent(component.key, 'moleFraction', event.target.value)
              }
              placeholder="0.0"
              aria-label={`Component ${index + 1} mole fraction`}
            />
            <button
              className="icon-button"
              type="button"
              onClick={() => removeComponent(component.key)}
              disabled={components.length === 1}
              aria-label={`Remove component ${index + 1}`}
              title="Remove component"
            >
              ×
            </button>
          </div>
        ))}
      </div>

      <button className="secondary-button" type="button" onClick={addComponent}>
        + Add component
      </button>

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
          The browser never normalizes z or evaluates EOS/flash equations.
        </div>
        <button className="primary-button" type="submit" disabled={disabled}>
          {busy ? 'Computing…' : 'Run PT flash'}
        </button>
      </div>
    </form>
  );
}
