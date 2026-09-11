import {
  physicalRoleLabel,
  topologyStatusLabel,
  topologyStatusSeverity,
  type ProfileCPtResponse,
} from '../domain/flash';

interface ResultPanelProps {
  result: ProfileCPtResponse | null;
  error: string | null;
}

function formatNumber(value: number, digits = 8): string {
  if (!Number.isFinite(value)) {
    return '—';
  }
  if (value === 0) {
    return '0';
  }
  const magnitude = Math.abs(value);
  if (magnitude < 1e-4 || magnitude >= 1e5) {
    return value.toExponential(6);
  }
  return value.toPrecision(digits);
}

export function ResultPanel({ result, error }: ResultPanelProps) {
  if (error) {
    return (
      <section className="result-panel result-empty" aria-live="polite">
        <p className="eyebrow">Compute service</p>
        <h2>Calculation unavailable</h2>
        <p>{error}</p>
      </section>
    );
  }

  if (!result) {
    return (
      <section className="result-panel result-empty" aria-live="polite">
        <p className="eyebrow">Result</p>
        <h2>No calculation yet</h2>
        <p>
          A validated backend response will appear here. The UI does not generate
          placeholder phase data.
        </p>
      </section>
    );
  }

  const severity = topologyStatusSeverity(result.status);

  return (
    <section className="result-panel" aria-live="polite">
      <div className="section-heading">
        <div>
          <p className="eyebrow">Result</p>
          <h2>PT topology candidate</h2>
        </div>
        <span className="status-chip" data-severity={severity}>
          {topologyStatusLabel(result.status)}
        </span>
      </div>

      <div className="result-contract-grid">
        <div>
          <span>Pressure</span>
          <strong>{formatNumber(result.pressurePa / 1e6)} MPa</strong>
        </div>
        <div>
          <span>Temperature</span>
          <strong>{formatNumber(result.temperatureK)} K</strong>
        </div>
        <div>
          <span>NaCl molality</span>
          <strong>{formatNumber(result.naclMolalityMolPerKgWater)} mol/kg H₂O</strong>
        </div>
        <div>
          <span>Returned phases</span>
          <strong>{result.phases.length}</strong>
        </div>
      </div>

      <div className="science-boundary">
        <strong>Finite-search candidate boundary</strong>
        <span>
          Global stability: not proven · Authoritative phase set: not published ·
          Nonaqueous L/V morphology: unresolved
        </span>
      </div>

      <div className="phase-stack">
        {result.phases.map((phase, phaseIndex) => (
          <article className="phase-card" key={`${phase.physicalRole}-${phaseIndex}`}>
            <div className="phase-card-heading">
              <div>
                <span className="phase-index">Phase {phaseIndex + 1}</span>
                <h3>{physicalRoleLabel(phase.physicalRole)}</h3>
              </div>
              <div className="phase-fraction">
                <span>β</span>
                <strong>{formatNumber(phase.molePhaseFraction)}</strong>
              </div>
            </div>

            <div className="phase-meta">
              <span>Family: {phase.thermodynamicFamily === 'aqueous' ? 'AQ' : 'NA'}</span>
              <span>
                Z: {phase.compressibilityFactor === undefined
                  ? 'not supplied'
                  : formatNumber(phase.compressibilityFactor)}
              </span>
            </div>

            <div className="composition-result-table">
              <div className="composition-result-head">
                <span>Component</span>
                <span>x</span>
              </div>
              {result.componentIds.map((componentId, componentIndex) => (
                <div className="composition-result-row" key={componentId}>
                  <span>{componentId}</span>
                  <code>
                    {formatNumber(phase.composition[componentIndex] ?? Number.NaN, 10)}
                  </code>
                </div>
              ))}
            </div>
          </article>
        ))}
      </div>

      <div className="diagnostic-block">
        <span>Backend diagnostic</span>
        <p>{result.diagnostic || 'No diagnostic text supplied.'}</p>
      </div>

      <details className="contract-details">
        <summary>Algorithm contract</summary>
        <dl>
          <div>
            <dt>Model</dt>
            <dd>{result.modelProfile}</dd>
          </div>
          <div>
            <dt>Algorithm</dt>
            <dd>{result.algorithmProfile}</dd>
          </div>
        </dl>
      </details>
    </section>
  );
}
