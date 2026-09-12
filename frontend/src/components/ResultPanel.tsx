import {
  ptOutcomeLabel,
  ptOutcomeSeverity,
  ptServiceErrorLabel,
  type PtFlashResponse,
} from '../domain/flash';

export interface PtClientFailure {
  kind: 'configuration' | 'transport' | 'contract';
  title: string;
  message: string;
}

interface ResultPanelProps {
  response: PtFlashResponse | null;
  clientFailure: PtClientFailure | null;
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

function identifierLabel(value: string): string {
  return value.replaceAll('_', ' ');
}

export function ResultPanel({ response, clientFailure }: ResultPanelProps) {
  if (clientFailure !== null) {
    return (
      <section
        className="result-panel result-empty client-failure"
        data-failure-kind={clientFailure.kind}
        aria-live="polite"
      >
        <p className="eyebrow">Client / RPC failure</p>
        <h2>{clientFailure.title}</h2>
        <p>{clientFailure.message}</p>
        <p className="failure-semantics">
          No scientific outcome was received. This state is separate from a backend
          computation reported as indeterminate.
        </p>
      </section>
    );
  }

  if (response === null) {
    return (
      <section className="result-panel result-empty" aria-live="polite">
        <p className="eyebrow">Result</p>
        <h2>No calculation yet</h2>
        <p>
          A service response will appear here. The UI does not generate placeholder
          phase data or infer a scientific outcome from an RPC status.
        </p>
      </section>
    );
  }

  if (response.kind === 'service_error') {
    return (
      <section className="result-panel result-empty service-error" aria-live="polite">
        <p className="eyebrow">PT service error</p>
        <h2>{ptServiceErrorLabel(response.error.code)}</h2>
        <div className="service-error-facts">
          <span>
            Code <code>{response.error.code}</code>
          </span>
          <span>
            Field <code>{response.error.field}</code>
          </span>
        </div>
        <p>{response.error.diagnostic}</p>
        <p className="failure-semantics">
          The service returned an error envelope, not an indeterminate thermodynamic
          result. No phases are presented.
        </p>
      </section>
    );
  }

  const result = response.result;
  const capability = result.provenance.backend.capability;
  const severity = ptOutcomeSeverity(result.outcome);
  const accepted = result.outcome === 'accepted';

  return (
    <section className="result-panel" aria-live="polite">
      <div className="section-heading">
        <div>
          <p className="eyebrow">Scientific computation result</p>
          <h2>{accepted ? 'Published PT phase set' : 'No accepted phase set'}</h2>
        </div>
        <span className="status-chip" data-severity={severity}>
          {ptOutcomeLabel(result.outcome)}
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
          <span>Accepted phases</span>
          <strong>{result.phases.length}</strong>
        </div>
        <div>
          <span>Supported counts</span>
          <strong>{capability.supportedPhaseCounts.join(', ')}</strong>
        </div>
      </div>

      <div className="science-boundary">
        <strong>Declared publication boundary</strong>
        <span>
          Global stability proven: {result.globalStabilityProven ? 'yes' : 'no'} ·
          Morphology resolved: {result.morphologyResolved ? 'yes' : 'no'} · Final
          phase-set review: {capability.performsFinalPhaseSetReview ? 'yes' : 'no'}
        </span>
      </div>

      {!accepted ? (
        <div className="scientific-outcome" data-outcome={result.outcome}>
          <strong>
            {result.outcome === 'indeterminate'
              ? 'Scientific decision remains indeterminate'
              : 'Candidate phase set was not accepted'}
          </strong>
          <p>
            {result.outcome === 'indeterminate'
              ? 'The RPC and service mapping succeeded, but the delegated backend could not publish an accepted phase set under its declared review contract.'
              : 'The delegated backend reported instability evidence during its phase-set review. Diagnostic candidates are intentionally withheld.'}
          </p>
        </div>
      ) : null}

      {accepted ? (
        <div className="phase-stack">
          {result.phases.map((phase) => (
            <article className="phase-card" key={phase.phaseIndex}>
              <div className="phase-card-heading">
                <div>
                  <span className="phase-index">Phase {phase.phaseIndex + 1}</span>
                  <h3>
                    {phase.providerMetadata === undefined
                      ? 'Provider phase instance'
                      : `Provider role: ${phase.providerMetadata.roleId}`}
                  </h3>
                </div>
                <div className="phase-fraction">
                  <span>β</span>
                  <strong>{formatNumber(phase.molePhaseFraction)}</strong>
                </div>
              </div>

              <div className="phase-meta">
                {phase.providerMetadata === undefined ? (
                  <span>Provider metadata: not supplied</span>
                ) : (
                  <>
                    <span>Family: {phase.providerMetadata.familyId}</span>
                    <span>Namespace: {capability.phaseMetadataNamespace}</span>
                  </>
                )}
                <span>Branch: {phase.providerBranch}</span>
                <span>Branch smooth: {phase.providerBranchSmooth ? 'yes' : 'no'}</span>
                <span>
                  Z:{' '}
                  {phase.compressibilityFactor === undefined
                    ? 'not supplied'
                    : formatNumber(phase.compressibilityFactor)}
                </span>
              </div>

              <div className="composition-result-table">
                <div className="composition-result-head">
                  <span>Component</span>
                  <span>x</span>
                  <span>ln φ</span>
                </div>
                {phase.components.map((component) => (
                  <div className="composition-result-row" key={component.componentId}>
                    <span>{component.componentId}</span>
                    <code>{formatNumber(component.moleFraction, 10)}</code>
                    <code>{formatNumber(component.lnFugacityCoefficient, 10)}</code>
                  </div>
                ))}
              </div>
            </article>
          ))}
        </div>
      ) : null}

      <div className="diagnostic-block">
        <span>Backend diagnostic</span>
        <p>{result.diagnostic || 'No diagnostic text supplied.'}</p>
      </div>

      {result.transitionReport.evidence.length > 0 ? (
        <details className="contract-details transition-details">
          <summary>Phase-transition evidence</summary>
          <ol>
            {result.transitionReport.evidence.map((evidence, index) => (
              <li key={`${evidence.providerEvidenceProfile}-${index}`}>
                <strong>{identifierLabel(evidence.trigger)}</strong>
                <span>
                  {evidence.sourcePhaseCount} → {evidence.targetPhaseCount ?? '?'} ·{' '}
                  {identifierLabel(evidence.resolution)}
                </span>
                <code>{evidence.providerEvidenceProfile}</code>
                {evidence.diagnostic ? <p>{evidence.diagnostic}</p> : null}
              </li>
            ))}
          </ol>
        </details>
      ) : null}

      <details className="contract-details">
        <summary>Backend and result provenance</summary>
        <dl>
          <div>
            <dt>Configured backend</dt>
            <dd>{result.provenance.backend.configuredBackendId}</dd>
          </div>
          <div>
            <dt>Adapter backend</dt>
            <dd>{capability.backendId}</dd>
          </div>
          <div>
            <dt>Model</dt>
            <dd>{capability.modelProfile}</dd>
          </div>
          <div>
            <dt>Algorithm</dt>
            <dd>{capability.algorithmProfile}</dd>
          </div>
          <div>
            <dt>Publication</dt>
            <dd>{capability.publicationProfile}</dd>
          </div>
          <div>
            <dt>Configuration</dt>
            <dd>{capability.configurationProfile}</dd>
          </div>
          <div>
            <dt>Dataset / revision</dt>
            <dd>
              {capability.datasetId} · {capability.revision}
            </dd>
          </div>
          <div>
            <dt>Provider result</dt>
            <dd>{result.provenance.providerResultConvention}</dd>
          </div>
        </dl>
      </details>
    </section>
  );
}
