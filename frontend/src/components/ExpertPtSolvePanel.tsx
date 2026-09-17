import { useEffect, useRef, useState, type FormEvent } from 'react';
import { toJson } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';

import type { ExpertOwnedModel } from '../api/expertModelOwner';
import {
  buildExpertPtRequest,
  emptyExpertPtInput,
  ExpertPtSolveGate,
  type ExpertPtInput,
} from '../api/expertPtSolve';
import { ModelClientError } from '../api/modelSessionClient';
import { RendererModelError } from '../api/rendererModelClient';
import type { ModelValidationDetail } from '../api/modelValidationDetail';
import {
  FullPtResultSchema,
  type FullPtResult,
  type ModelSnapshot,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { PtComputationOutcome } from '../gen/mpmc/runtime/v1/pt_service_pb';

export interface ExpertPtFailure {
  readonly reason: string;
  readonly code?: Code;
  readonly validation?: ModelValidationDetail;
}

export function expertPtFailure(cause: unknown): ExpertPtFailure {
  if (cause instanceof ModelClientError || cause instanceof RendererModelError) {
    return {
      reason: cause.reason,
      code: cause.code,
      ...(cause.validation === undefined ? {} : { validation: cause.validation }),
    };
  }
  return { reason: 'expert.solve_failed' };
}

function formatNumber(value: number | undefined, digits = 8): string {
  if (value === undefined || !Number.isFinite(value)) return 'not supplied';
  if (value === 0) return '0';
  const magnitude = Math.abs(value);
  if (magnitude < 1e-4 || magnitude >= 1e5) return value.toExponential(6);
  return value.toPrecision(digits);
}

function percentLabel(value: number | undefined): string {
  return value === undefined || !Number.isFinite(value)
    ? '—'
    : `${(value * 100).toFixed(3)}%`;
}

function percentWidth(value: number | undefined): string {
  if (value === undefined || !Number.isFinite(value)) return '0%';
  return `${Math.min(100, Math.max(0, value * 100))}%`;
}

function pressureMpaToPaInput(text: string): string {
  const value = text.trim();
  if (!/^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/iu.test(value)) return text;
  const parsed = Number(value);
  return Number.isFinite(parsed) ? String(parsed * 1e6) : text;
}

function feedSum(input: ExpertPtInput): number {
  return input.feed.reduce((sum, entry) => {
    const value = Number(entry.fraction);
    return Number.isFinite(value) ? sum + value : sum;
  }, 0);
}

/** Display native outcomes, never infer phase identity from density, Z or order. */
export function ExpertPtResultView({ result, snapshot }: { result: FullPtResult; snapshot: ModelSnapshot }) {
  const accepted = result.outcome === PtComputationOutcome.ACCEPTED;
  const phases = result.candidatePhaseSet?.phases ?? [];
  const components = snapshot.definition?.components ?? [];
  const label = accepted
    ? 'Equilibrium result'
    : result.outcome === PtComputationOutcome.PHASE_SET_UNSTABLE
      ? 'No accepted equilibrium phase set'
      : 'Calculation is indeterminate';

  const exactTable = (
    <div className="table-scroll">
      <table className="result-table">
        <thead>
          <tr>
            <th>Phase</th>
            <th>Phase fraction</th>
            {components.map((component, index) => (
              <th key={index}>x({component.componentId})</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {phases.map((phase, index) => (
            <tr key={index}>
              <th>Phase {index + 1}</th>
              <td>{formatNumber(phase.molePhaseFraction, 10)}</td>
              {components.map((_component, i) => (
                <td key={i}>{formatNumber(phase.composition[i], 10)}</td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );

  return (
    <section
      className="result-panel pr-result-panel"
      data-expert-result={accepted ? 'accepted' : 'not-accepted'}
      aria-live="polite"
    >
      <div className="section-heading">
        <div>
          <p className="eyebrow">PR flash result</p>
          <h2>{label}</h2>
        </div>
        {accepted ? <span className="status-chip" data-severity="accepted">Accepted</span> : null}
      </div>

      <div className="result-contract-grid pr-result-summary">
        <div>
          <span>Pressure</span>
          <strong>{result.pressurePa === undefined ? '—' : `${formatNumber(result.pressurePa / 1e6)} MPa`}</strong>
        </div>
        <div>
          <span>Temperature</span>
          <strong>{result.temperatureK === undefined ? '—' : `${formatNumber(result.temperatureK)} K`}</strong>
        </div>
        <div>
          <span>Phase count</span>
          <strong>{accepted ? phases.length : '—'}</strong>
        </div>
        <div>
          <span>Model</span>
          <strong>Peng–Robinson (PR76)</strong>
        </div>
      </div>

      {accepted ? (
        phases.length > 0 ? (
          <>
            <section className="result-chart-section" aria-label="Phase fraction chart">
              <div className="composition-header">
                <div>
                  <h3>Phase fractions</h3>
                  <p>Fraction of total moles assigned to each accepted phase.</p>
                </div>
              </div>
              <div className="fraction-chart">
                {phases.map((phase, index) => (
                  <div className="fraction-chart-row" key={index}>
                    <span>Phase {index + 1}</span>
                    <div className="fraction-chart-track" aria-hidden="true">
                      <span style={{ width: percentWidth(phase.molePhaseFraction) }} />
                    </div>
                    <strong>{percentLabel(phase.molePhaseFraction)}</strong>
                  </div>
                ))}
              </div>
            </section>

            <section className="result-chart-section" aria-label="Phase composition charts">
              <div className="composition-header">
                <div>
                  <h3>Composition in each phase</h3>
                  <p>Component mole fractions are shown without relabeling phase identity.</p>
                </div>
              </div>
              <div className="phase-stack">
                {phases.map((phase, phaseIndex) => (
                  <article className="phase-card pr-phase-card" key={phaseIndex}>
                    <div className="phase-card-heading">
                      <div>
                        <span className="phase-index">Phase {phaseIndex + 1}</span>
                        <h3>Composition</h3>
                      </div>
                      <div className="phase-fraction">
                        <span>phase fraction</span>
                        <strong>{percentLabel(phase.molePhaseFraction)}</strong>
                      </div>
                    </div>
                    <div className="composition-bars">
                      {components.map((component, componentIndex) => {
                        const fraction = phase.composition[componentIndex];
                        return (
                          <div className="composition-bar-row" key={`${phaseIndex}-${componentIndex}`}>
                            <span title={component.displayName || component.componentId}>
                              {component.displayName || component.componentId}
                            </span>
                            <div className="composition-bar-track" aria-hidden="true">
                              <span style={{ width: percentWidth(fraction) }} />
                            </div>
                            <strong>{percentLabel(fraction)}</strong>
                          </div>
                        );
                      })}
                    </div>
                  </article>
                ))}
              </div>
            </section>

            <section className="exact-results" aria-label="Exact phase mole fractions">
              <div className="composition-header">
                <div>
                  <h3>Exact mole fractions</h3>
                  <p>Numerical values corresponding to the charts above.</p>
                </div>
              </div>
              {exactTable}
            </section>
          </>
        ) : (
          <p role="alert">Invalid result: an accepted outcome contains no phase data.</p>
        )
      ) : (
        <div className="scientific-outcome" data-outcome={result.outcome === PtComputationOutcome.INDETERMINATE ? 'indeterminate' : 'phase_set_unstable'}>
          <strong>No equilibrium phase count or composition is published.</strong>
          <p>
            The calculation completed without an accepted final phase set. Adjust the
            inputs or model data rather than interpreting diagnostic candidates as results.
          </p>
          {phases.length > 0 ? (
            <details data-candidate-only="true" className="advanced-details">
              <summary>Diagnostic candidate phase set ({phases.length})</summary>
              {exactTable}
            </details>
          ) : null}
        </div>
      )}

      <details className="contract-details advanced-details">
        <summary>Advanced calculation details</summary>
        <div className="diagnostic-block">
          <span>Native diagnostic</span>
          <p>{result.diagnostic || 'No additional native diagnostic was supplied.'}</p>
        </div>
        <p>
          Initial stability search: {result.capability?.performsInitialStabilitySearch ? 'enabled' : 'not declared'};
          final phase-set review: {result.capability?.performsFinalPhaseSetReview ? 'enabled' : 'not declared'}.
        </p>
        <p>
          {result.globalStabilityProven === true
            ? 'Global stability is reported proven by the backend.'
            : 'Global stability is not proven; acceptance under configured searches is not a global proof.'}
        </p>
        <details className="contract-details">
          <summary>Complete native result and transition evidence</summary>
          <pre>{JSON.stringify(toJson(FullPtResultSchema, result), null, 2)}</pre>
        </details>
      </details>
    </section>
  );
}

export interface ExpertPtSolvePanelProps {
  model: ExpertOwnedModel;
  /** Unapplied parameter/component edits must not be solved using the old model. */
  blocked?: boolean;
}

/** The workspace keys this component by applied model revision. No login path. */
export function ExpertPtSolvePanel({ model, blocked = false }: ExpertPtSolvePanelProps) {
  const [input, setInput] = useState(() => emptyExpertPtInput(model.snapshot));
  const [pressureMpa, setPressureMpa] = useState('');
  const [result, setResult] = useState<FullPtResult | null>(null);
  const [failure, setFailure] = useState<ExpertPtFailure | null>(null);
  const [busy, setBusy] = useState(false);
  const [gate] = useState(() => new ExpertPtSolveGate());
  const alive = useRef(true);

  useEffect(() => {
    alive.current = true;
    return () => { alive.current = false; gate.invalidate(); };
  }, [gate]);
  useEffect(() => {
    gate.invalidate();
    setInput(emptyExpertPtInput(model.snapshot));
    setPressureMpa('');
    setResult(null);
    setFailure(null);
  }, [model, gate]);
  useEffect(() => {
    if (blocked) { gate.invalidate(); setResult(null); setFailure(null); }
  }, [blocked, gate]);

  function invalidateResult() {
    gate.invalidate();
    setResult(null);
    setFailure(null);
  }

  function edit(next: ExpertPtInput) {
    invalidateResult();
    setInput(next);
  }

  function editPressure(value: string) {
    invalidateResult();
    setPressureMpa(value);
  }

  async function submit(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    if (blocked || gate.busy || model.released) return;
    setResult(null);
    setFailure(null);
    try {
      const request = buildExpertPtRequest(model.snapshot, {
        ...input,
        pressurePa: pressureMpaToPaInput(pressureMpa),
      });
      setBusy(true);
      const solved = await gate.solve(model, request);
      if (alive.current && solved !== null) setResult(solved);
    } catch (cause) {
      if (alive.current) setFailure(expertPtFailure(cause));
    } finally {
      if (alive.current) setBusy(false);
    }
  }

  const compositionTotal = feedSum(input);

  return (
    <section data-expert-pt-solve="true">
      <form className="flash-form pr-solve-form" onSubmit={event => void submit(event)} noValidate>
        <div className="section-heading">
          <div>
            <p className="eyebrow">Flash conditions</p>
            <h2>Set P, T and initial composition</h2>
          </div>
          <span className="model-chip">Peng–Robinson (PR76)</span>
        </div>
        <p className="form-intro">
          Values are sent exactly as entered except for the displayed MPa-to-Pa unit conversion.
          Mole fractions are never filled, clipped or normalized by the UI.
        </p>
        {blocked ? (
          <div className="connection-banner" role="status">
            <strong>Fluid data changed.</strong>
            <span>Apply the updated model before calculating again.</span>
          </div>
        ) : null}
        <div className="field-grid field-grid-two">
          <label>
            <span>Pressure</span>
            <div className="input-with-unit">
              <input
                aria-label="PR pressure MPa"
                inputMode="decimal"
                value={pressureMpa}
                onChange={event => editPressure(event.currentTarget.value)}
              />
              <span>MPa</span>
            </div>
          </label>
          <label>
            <span>Temperature</span>
            <div className="input-with-unit">
              <input
                aria-label="PR temperature K"
                inputMode="decimal"
                value={input.temperatureK}
                onChange={event => edit({ ...input, temperatureK: event.currentTarget.value })}
              />
              <span>K</span>
            </div>
          </label>
        </div>

        <div className="composition-header">
          <div>
            <h3>Initial overall composition</h3>
            <p>Enter mole fraction z for every component in the applied model.</p>
          </div>
          <div className="composition-sum">Σz = {compositionTotal.toPrecision(8)}</div>
        </div>
        <div className="component-table inventory-component-table" role="group" aria-label="PR overall composition">
          <div className="component-table-head" aria-hidden="true">
            <span>Component</span>
            <span>Mole fraction z</span>
          </div>
          {input.feed.map((entry, index) => (
            <div className="component-row" key={entry.componentId}>
              <code>{entry.componentId}</code>
              <input
                aria-label={`PR feed ${entry.componentId}`}
                inputMode="decimal"
                value={entry.fraction}
                onChange={event => {
                  const fraction = event.currentTarget.value;
                  edit({
                    ...input,
                    feed: input.feed.map((item, i) => i === index ? { ...item, fraction } : item),
                  });
                }}
              />
            </div>
          ))}
        </div>

        <div className="submit-row">
          <span className="submit-note">The native solver determines phase count and equilibrium compositions.</span>
          <button className="primary-button" type="submit" disabled={blocked || busy || model.released}>
            {busy ? 'Calculating…' : 'Run PR flash'}
          </button>
        </div>
      </form>

      {!blocked && failure ? (
        <div className="validation-panel pr-solve-failure" role="alert">
          <strong>Calculation could not be completed.</strong>
          <p>
            Check pressure, temperature, composition and applied component data, then try again.
            {failure.validation?.field ? <> Problem field: <code>{failure.validation.field}</code>.</> : null}
          </p>
          <details className="inline-details">
            <summary>Technical details</summary>
            <code>{failure.reason}</code>
            {failure.code === undefined ? null : <code> · {failure.code}</code>}
          </details>
        </div>
      ) : null}
      {!blocked && result ? <ExpertPtResultView result={result} snapshot={model.snapshot} /> : null}
    </section>
  );
}
