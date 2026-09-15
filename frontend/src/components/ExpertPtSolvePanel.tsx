import { useEffect, useRef, useState, type FormEvent } from 'react';
import { toJson } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';

import type { ExpertOwnedModel } from '../api/expertModelOwner';
import { buildExpertPtRequest, emptyExpertPtInput, ExpertPtSolveGate, type ExpertPtInput } from '../api/expertPtSolve';
import { ModelClientError } from '../api/modelSessionClient';
import { RendererModelError } from '../api/rendererModelClient';
import type { ModelValidationDetail } from '../api/modelValidationDetail';
import { FullPtResultSchema, type FullPtResult, type ModelSnapshot } from '../gen/mpmc/model_configuration/v1/model_service_pb';
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

/** Display native outcomes, never infer stability or phase labels from density/Z. */
export function ExpertPtResultView({ result, snapshot }: { result: FullPtResult; snapshot: ModelSnapshot }) {
  const accepted = result.outcome === PtComputationOutcome.ACCEPTED;
  const phases = result.candidatePhaseSet?.phases ?? [];
  const components = snapshot.definition?.components ?? [];
  const label = accepted ? 'Accepted phase set'
    : result.outcome === PtComputationOutcome.PHASE_SET_UNSTABLE ? 'Phase set is unstable'
    : 'Indeterminate — no accepted phase set';
  const table = (
    <div className="table-scroll">
      <table>
        <thead>
          <tr><th>Phase</th><th>Mole phase fraction</th>{components.map((component, index) => (
            <th key={index}>x({component.componentId})</th>
          ))}<th>Z</th></tr>
        </thead>
        <tbody>{phases.map((phase, index) => (
          <tr key={index}>
            <th>Phase {index + 1}</th>
            <td>{phase.molePhaseFraction ?? 'not supplied'}</td>
            {components.map((_component, i) => <td key={i}>{phase.composition[i] ?? 'not supplied'}</td>)}
            <td>{phase.compressibilityFactor ?? 'not supplied'}</td>
          </tr>
        ))}</tbody>
      </table>
    </div>
  );
  return (
    <section className="result-panel" data-expert-result={accepted ? 'accepted' : 'not-accepted'} aria-live="polite">
      <h2>{label}</h2>
      <p>Model <code>{snapshot.definition?.datasetId}</code> · revision <code>{snapshot.definition?.revision}</code></p>
      <p>P = {result.pressurePa ?? 'not supplied'} Pa · T = {result.temperatureK ?? 'not supplied'} K</p>
      {accepted ? (
        phases.length > 0 ? <><h3>Accepted phase count: {phases.length}</h3>{table}</>
          : <p role="alert">Invalid result: accepted outcome has no phase data.</p>
      ) : (
        <>
          <p>No equilibrium phase count or composition is accepted for this calculation.</p>
          {phases.length > 0 ? <details data-candidate-only="true"><summary>Diagnostic candidates — not accepted results ({phases.length})</summary>{table}</details> : null}
        </>
      )}
      <h3>Stability and convergence diagnostics</h3>
      <p>{result.diagnostic || 'No additional native diagnostic was supplied.'}</p>
      <p>Initial stability search: {result.capability?.performsInitialStabilitySearch ? 'enabled' : 'not declared'}; final phase-set review: {result.capability?.performsFinalPhaseSetReview ? 'enabled' : 'not declared'}.</p>
      <p>{result.globalStabilityProven === true
        ? 'Global stability is reported proven by the backend.'
        : 'Global stability is not proven. Acceptance under the configured searches is not a global proof.'}</p>
      <details className="contract-details"><summary>Complete native result and transition evidence</summary>
        <pre>{JSON.stringify(toJson(FullPtResultSchema, result), null, 2)}</pre>
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
    setResult(null);
    setFailure(null);
  }, [model, gate]);
  useEffect(() => {
    if (blocked) { gate.invalidate(); setResult(null); setFailure(null); }
  }, [blocked, gate]);

  function edit(next: ExpertPtInput) {
    gate.invalidate();
    setResult(null);
    setFailure(null);
    setInput(next);
  }

  async function submit(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    if (blocked || gate.busy || model.released) return;
    setResult(null);
    setFailure(null);
    try {
      const request = buildExpertPtRequest(model.snapshot, input);
      setBusy(true);
      const solved = await gate.solve(model, request);
      if (alive.current && solved !== null) setResult(solved);
    } catch (cause) {
      if (alive.current) setFailure(expertPtFailure(cause));
    } finally {
      if (alive.current) setBusy(false);
    }
  }

  return (
    <section data-expert-pt-solve="true">
      <form className="flash-form" onSubmit={event => void submit(event)} noValidate>
        <h2>Calculate with the applied PR76 model</h2>
        <p>Enter pressure in Pa, temperature in K and mole fractions in the displayed component order. No fractions are filled or normalized automatically.</p>
        {blocked ? <p role="status">Parameters have changed. Apply a new immutable model before calculating; old results are hidden.</p> : null}
        <div className="field-grid field-grid-two">
          <label><span>Pressure [Pa]</span><input aria-label="Expert pressure Pa" inputMode="decimal" value={input.pressurePa}
            onChange={event => edit({ ...input, pressurePa: event.currentTarget.value })} /></label>
          <label><span>Temperature [K]</span><input aria-label="Expert temperature K" inputMode="decimal" value={input.temperatureK}
            onChange={event => edit({ ...input, temperatureK: event.currentTarget.value })} /></label>
        </div>
        <div className="field-grid field-grid-two">{input.feed.map((entry, index) => (
          <label key={entry.componentId}><span>z({entry.componentId})</span>
            <input aria-label={`Expert feed ${entry.componentId}`} inputMode="decimal" value={entry.fraction}
              onChange={event => {
                const fraction = event.currentTarget.value;
                edit({ ...input, feed: input.feed.map((item, i) => i === index ? { ...item, fraction } : item) });
              }} />
          </label>
        ))}</div>
        <p className="submit-note">The C++ solver performs initial stability search, phase splitting up to three phases and final phase-set review. Editing inputs discards an in-flight result without cancelling the model.</p>
        <button className="primary-button" type="submit" disabled={blocked || busy || model.released}>
          {busy ? 'Computing…' : 'Compute PR76 flash'}
        </button>
      </form>
      {!blocked && failure ? <div className="validation-panel" role="alert">
        <strong>Calculation unavailable</strong><p>Reason <code>{failure.reason}</code>{failure.code === undefined ? null : <> · code {failure.code}</>}</p>
        {failure.validation ? <p><code>{failure.validation.code}</code> · field <code>{failure.validation.field ?? 'not supplied'}</code></p> : null}
        <p>No result is published. Correct the input; if the internal model has expired, re-open the local Expert workspace and apply the model again. No account is required.</p>
      </div> : null}
      {!blocked && result ? <ExpertPtResultView result={result} snapshot={model.snapshot} /> : null}
    </section>
  );
}
