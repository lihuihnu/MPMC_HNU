import type { ReactNode } from 'react';

import type { ExpertModelInspection } from '../api/expertModelInspector';
import type { ModelValidationDetail } from '../api/modelValidationDetail';
import {
  ComponentKind,
  ModelFamily,
  SolverSettingsKind,
  SourceKind,
  type FullPtResult,
  type ModelProvenance,
  type ModelScalar,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { PtComputationOutcome } from '../gen/mpmc/runtime/v1/pt_service_pb';

export interface ExpertModelInspectorProps {
  inspection: ExpertModelInspection;
  /** Structured service/client detail only; never inferred from an error string. */
  validation?: ModelValidationDetail;
  /** Optional complete solve envelope for read-only diagnostic presentation. */
  result?: FullPtResult;
}

function enumLabel(value: number | undefined, labels: Readonly<Record<number, string>>): string {
  return value === undefined ? 'not supplied' : labels[value] ?? `unknown (${value})`;
}

const familyLabels: Readonly<Record<number, string>> = {
  [ModelFamily.PR76]: 'PR76',
  [ModelFamily.SW92]: 'Søreide–Whitson (SW92)',
  [ModelFamily.CPA]: 'CPA',
};
const componentKindLabels: Readonly<Record<number, string>> = {
  [ComponentKind.PURE]: 'pure',
  [ComponentKind.PSEUDO]: 'pseudo',
};
const sourceKindLabels: Readonly<Record<number, string>> = {
  [SourceKind.LITERATURE]: 'literature',
  [SourceKind.DATABASE]: 'database',
  [SourceKind.USER_SUPPLIED]: 'user supplied',
  [SourceKind.ASSUMPTION]: 'assumption',
  [SourceKind.SYNTHETIC_TEST]: 'synthetic test',
};
const settingsKindLabels: Readonly<Record<number, string>> = {
  [SolverSettingsKind.PRESET]: 'preset',
  [SolverSettingsKind.CUSTOM]: 'custom',
};
const outcomeLabels: Readonly<Record<number, string>> = {
  [PtComputationOutcome.ACCEPTED]: 'accepted',
  [PtComputationOutcome.PHASE_SET_UNSTABLE]: 'phase set unstable',
  [PtComputationOutcome.INDETERMINATE]: 'indeterminate',
};

function valueText(value: unknown): string {
  if (typeof value === 'bigint') return value.toString();
  if (typeof value === 'boolean') return value ? 'true' : 'false';
  if (typeof value === 'number' || typeof value === 'string') return String(value);
  return 'not supplied';
}

function humanize(name: string): string {
  const spaced = name.replace(/([a-z0-9])([A-Z])/gu, '$1 $2').replaceAll('_', ' ');
  return `${spaced.slice(0, 1).toUpperCase()}${spaced.slice(1)}`;
}

function Provenance({ value }: { value: ModelProvenance | undefined }) {
  if (!value) return <span className="phase-meta">Provenance: not declared</span>;
  const facts: Array<[string, string]> = [
    ['Kind', enumLabel(value.kind, sourceKindLabels)],
    ['Reference', value.reference ?? 'not supplied'],
    ['Revision', value.revision ?? 'not supplied'],
    ['Locator', value.locator ?? 'not supplied'],
    ['Note', value.note ?? 'not supplied'],
    ['Acquisition', value.acquisition ?? 'not supplied'],
    ['Usage terms', value.usageTerms ?? 'not supplied'],
  ];
  return (
    <dl className="contract-details">
      {facts.map(([label, fact]) => (
        <div key={label}>
          <dt>{label}</dt>
          <dd>{fact}</dd>
        </div>
      ))}
    </dl>
  );
}

function Scalar({
  label,
  scalar,
  unit,
}: {
  label: string;
  scalar: ModelScalar | undefined;
  unit?: string;
}) {
  return (
    <div className="model-contract">
      <span>{label}</span>
      <strong>
        {scalar?.value === undefined ? 'not supplied' : `${scalar.value}${unit ? ` ${unit}` : ''}`}
      </strong>
      {scalar?.originalUnit ? <code>original unit: {scalar.originalUnit}</code> : null}
      {scalar?.conversion ? <code>conversion: {scalar.conversion}</code> : null}
      <Provenance value={scalar?.provenance} />
    </div>
  );
}

function SettingsBlock({ title, value }: { title: string; value: object | undefined }) {
  const entries = value
    ? Object.entries(value as Record<string, unknown>).filter(
        ([name, item]) => name !== '$typeName' && item !== undefined,
      )
    : [];
  return (
    <details className="contract-details" open>
      <summary>{title}</summary>
      {entries.length === 0 ? (
        <p>Not supplied.</p>
      ) : (
        <dl>
          {entries.map(([name, item]) => (
            <div key={name}>
              <dt>{humanize(name)}</dt>
              <dd>
                <code>{valueText(item)}</code>
              </dd>
            </div>
          ))}
        </dl>
      )}
    </details>
  );
}

function Endpoint({
  label,
  value,
  exclusive,
  unit,
}: {
  label: string;
  value: number | undefined;
  exclusive: boolean;
  unit: string;
}) {
  return (
    <div data-known={value !== undefined}>
      <span>{label}</span>
      <strong>
        {value === undefined ? 'unknown' : `${value} ${unit} · ${exclusive ? 'open' : 'closed'}`}
      </strong>
    </div>
  );
}

function DiagnosticSection({
  validation,
  result,
}: {
  validation: ModelValidationDetail | undefined;
  result: FullPtResult | undefined;
}) {
  if (!validation && !result) return null;
  return (
    <section aria-label="Structured model diagnostics">
      <div className="section-heading">
        <div>
          <p className="eyebrow">Diagnostics</p>
          <h2>Structured service state</h2>
        </div>
      </div>
      {validation ? (
        <div className="validation-panel" data-validation-code={validation.code}>
          <strong>Validation detail</strong>
          <ul>
            <li>
              Code <code>{validation.code}</code>
            </li>
            <li>
              Field <code>{validation.field ?? 'not supplied'}</code>
            </li>
          </ul>
        </div>
      ) : null}
      {result ? (
        <div className="diagnostic-block" data-outcome={enumLabel(result.outcome, outcomeLabels)}>
          <span>Complete result diagnostic</span>
          <p>
            Outcome: <strong>{enumLabel(result.outcome, outcomeLabels)}</strong> · Candidate phases:{' '}
            <strong>{result.candidatePhaseSet?.phases.length ?? 0}</strong> · Global stability proven:{' '}
            <strong>{result.globalStabilityProven ? 'yes' : 'no'}</strong> · Morphology resolved:{' '}
            <strong>{result.morphologyResolved ? 'yes' : 'no'}</strong>
          </p>
          <p>{result.diagnostic ?? 'No diagnostic text supplied.'}</p>
          <p>
            Transition evidence: <strong>{result.transitionReport?.evidence.length ?? 0}</strong>
          </p>
        </div>
      ) : null}
    </section>
  );
}

function Section({ title, eyebrow, children }: { title: string; eyebrow: string; children: ReactNode }) {
  return (
    <section>
      <div className="section-heading">
        <div>
          <p className="eyebrow">{eyebrow}</p>
          <h2>{title}</h2>
        </div>
      </div>
      {children}
    </section>
  );
}

/** Read-only presentation of the exact frozen model-service snapshot. */
export function ExpertModelInspector({ inspection, validation, result }: ExpertModelInspectorProps) {
  const { snapshot } = inspection;
  const definition = snapshot.definition;
  const settings = snapshot.settings;
  const capability = snapshot.capability;
  const applicability = definition?.applicability;
  const pr76 = definition?.parameters.case === 'pr76' ? definition.parameters.value : undefined;

  return (
    <article className="result-panel expert-model-inspector" aria-label="Expert thermodynamic model inspector">
      <div className="section-heading">
        <div>
          <p className="eyebrow">Expert model inspector</p>
          <h2>{definition?.displayName || definition?.datasetId || 'Unnamed thermodynamic model'}</h2>
        </div>
        <span className="model-chip">{enumLabel(definition?.family, familyLabels)}</span>
      </div>

      <div className="model-contract">
        <span>Snapshot identity</span>
        <code>definition: {definition?.version ?? 'not supplied'}</code>
        <code>dataset: {definition?.datasetId ?? 'not supplied'}</code>
        <code>revision: {definition?.revision ?? 'not supplied'}</code>
        <code>publication: {capability?.publicationProfile ?? 'not supplied'}</code>
      </div>

      <Section eyebrow="Ordered inventory" title="Components">
        <div className="component-table inventory-component-table">
          <div className="component-table-head">
            <span>Ordered component</span>
            <span>Kind</span>
          </div>
          {(definition?.components ?? []).map((component, index) => (
            <div className="component-row" key={`${component.componentId}-${index}`} data-component-index={index}>
              <div>
                <code>{index}: {component.componentId || 'missing id'}</code>
                <div className="phase-meta">
                  <span>{component.displayName || 'no display name'}</span>
                  {component.molarMassKgPerMol?.value === undefined ? null : (
                    <span>molar mass: {component.molarMassKgPerMol.value} kg/mol</span>
                  )}
                </div>
                <Provenance value={component.provenance} />
              </div>
              <strong>{enumLabel(component.kind, componentKindLabels)}</strong>
            </div>
          ))}
        </div>
      </Section>

      <Section eyebrow="Declared validity" title="Applicability endpoints">
        <div className="result-contract-grid">
          <Endpoint label="Temperature lower" value={applicability?.temperatureLowerK} exclusive={applicability?.temperatureLowerExclusive ?? false} unit="K" />
          <Endpoint label="Temperature upper" value={applicability?.temperatureUpperK} exclusive={applicability?.temperatureUpperExclusive ?? false} unit="K" />
          <Endpoint label="Pressure lower" value={applicability?.pressureLowerPa} exclusive={applicability?.pressureLowerExclusive ?? false} unit="Pa" />
          <Endpoint label="Pressure upper" value={applicability?.pressureUpperPa} exclusive={applicability?.pressureUpperExclusive ?? false} unit="Pa" />
        </div>
        <Provenance value={applicability?.provenance} />
      </Section>

      <Section eyebrow="PR76 parameter snapshot" title="Pure parameters and kij">
        {pr76 ? (
          <>
            {pr76.pure.map((pure, index) => (
              <details className="contract-details" open key={`${pure.componentId}-${index}`}>
                <summary>{index}: {pure.componentId || 'missing component id'}</summary>
                <Scalar label="Critical temperature" scalar={pure.criticalTemperatureK} unit="K" />
                <Scalar label="Critical pressure" scalar={pure.criticalPressurePa} unit="Pa" />
                <Scalar label="Acentric factor" scalar={pure.acentricFactor} />
              </details>
            ))}
            {pr76.binary.length === 0 ? <p>No binary interaction records declared.</p> : null}
            {pr76.binary.map((pair, index) => (
              <details className="contract-details" open key={`${pair.firstComponentId}-${pair.secondComponentId}-${index}`}>
                <summary>
                  kij {pair.firstComponentId || '?'} ↔ {pair.secondComponentId || '?'}
                </summary>
                <Scalar label="Binary interaction kij" scalar={pair.kij} />
              </details>
            ))}
          </>
        ) : (
          <p>No PR76 parameter payload is present in this snapshot.</p>
        )}
      </Section>

      <Section eyebrow="Frozen solver contract" title="Solver settings">
        <div className="model-contract">
          <span>Settings identity</span>
          <code>version: {settings?.version ?? 'not supplied'}</code>
          <code>kind: {enumLabel(settings?.kind, settingsKindLabels)}</code>
          <code>preset: {settings?.presetId || 'none (custom snapshot)'}</code>
        </div>
        <SettingsBlock title="EOS root" value={settings?.eosRoot} />
        <SettingsBlock title="Initial stability" value={settings?.initialStability} />
        <SettingsBlock title="Two phase" value={settings?.twoPhase} />
        <SettingsBlock title="Final two-phase stability" value={settings?.finalTwoPhaseStability} />
        <SettingsBlock title="Three phase" value={settings?.threePhase} />
        <SettingsBlock title="Final three-phase stability" value={settings?.finalThreePhaseStability} />
      </Section>

      <DiagnosticSection validation={validation} result={result} />
    </article>
  );
}
