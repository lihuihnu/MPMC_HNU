import { create, fromJson } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';
import { renderToStaticMarkup } from 'react-dom/server';
import { describe, expect, it, vi } from 'vitest';
import type { ExpertOwnedModel } from '../api/expertModelOwner';
import { ModelClientError } from '../api/modelSessionClient';
import { MODEL_VALIDATION_DETAIL_VERSION } from '../api/modelValidationDetail';
import { FullPtResultSchema, ModelSnapshotSchema } from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { PtComputationOutcome } from '../gen/mpmc/runtime/v1/pt_service_pb';
import { expertSnapshotJson } from '../test/modelInspectorFixtures';
import { ExpertPtResultView, ExpertPtSolvePanel, expertPtFailure } from './ExpertPtSolvePanel';

// Deliberately synthetic presentation data, not a three-phase physical reference.
function fixture(outcome: PtComputationOutcome) {
  return create(FullPtResultSchema, {
    outcome, pressurePa: 2e6, temperatureK: 300, feed: [0.25, 0.75],
    maximumPhaseCount: 3, globalStabilityProven: false, diagnostic: 'native-test-diagnostic',
    candidatePhaseSet: { phases: [0.2, 0.3, 0.5].map(molePhaseFraction => ({
      molePhaseFraction, composition: [0.25, 0.75], lnFugacityCoefficient: [0, 0],
    })) },
  });
}
function model(): ExpertOwnedModel {
  const snapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
  return { snapshot, source: { describe: async () => ({ snapshot }) }, released: false,
    solve: vi.fn(async () => fixture(PtComputationOutcome.INDETERMINATE)), release: async () => {} };
}

describe('custom PR76 calculation presentation', () => {
  it('renders product-level P, T and ordered feed inputs without invented results', () => {
    const owned = model();
    const html = renderToStaticMarkup(<ExpertPtSolvePanel model={owned} />);
    expect(html).toContain('PR pressure MPa');
    expect(html).toContain('PR temperature K');
    expect(html).toContain('PR feed methane');
    expect(html).toContain('PR feed ethane');
    expect(html).toContain('Run PR flash');
    expect(html).toContain('Mole fractions are never filled, clipped or normalized');
    expect(html).not.toContain('data-expert-result');
    expect(owned.solve).not.toHaveBeenCalled();
  });

  it('blocks solving the old model while parameter edits are unapplied', () => {
    const html = renderToStaticMarkup(<ExpertPtSolvePanel model={model()} blocked />);
    expect(html).toContain('Fluid data changed');
    expect(html).toContain('type="submit" disabled=""');
    expect(html).not.toContain('data-expert-result');
  });

  it('shows accepted three-phase count, fractions, compositions and exact mole fractions', () => {
    const html = renderToStaticMarkup(
      <ExpertPtResultView result={fixture(PtComputationOutcome.ACCEPTED)} snapshot={model().snapshot} />,
    );
    expect(html).toContain('data-expert-result="accepted"');
    expect(html).toContain('Equilibrium result');
    expect(html).toContain('Phase count');
    expect(html).toContain('Phase fractions');
    expect(html).toContain('20.000%');
    expect(html).toContain('30.000%');
    expect(html).toContain('50.000%');
    expect(html).toContain('Composition in each phase');
    expect(html).toContain('x(methane)');
    expect(html).toContain('x(ethane)');
    expect(html).toContain('Advanced calculation details');
    expect(html).toContain('Global stability is not proven');
    expect(html).toContain('native-test-diagnostic');
    expect(html).not.toContain('data-candidate-only');
  });

  it.each([PtComputationOutcome.INDETERMINATE, PtComputationOutcome.PHASE_SET_UNSTABLE])(
    'keeps diagnostic candidate phases separate from results for outcome %s',
    (outcome) => {
      const html = renderToStaticMarkup(
        <ExpertPtResultView result={fixture(outcome)} snapshot={model().snapshot} />,
      );
      expect(html).toContain('data-expert-result="not-accepted"');
      expect(html).toContain('data-candidate-only="true"');
      expect(html).not.toContain('Equilibrium result');
      expect(html).toContain('Advanced calculation details');
      expect(html).toContain('Complete native result and transition evidence');
    },
  );

  it('preserves structured errors and never prints arbitrary backend exception text', () => {
    const validation = { version: MODEL_VALIDATION_DETAIL_VERSION, code: 'request.rejected', field: 'feed[1]' };
    expect(expertPtFailure(new ModelClientError(Code.InvalidArgument, 'rpc.failed', validation)))
      .toEqual({ reason: 'rpc.failed', code: Code.InvalidArgument, validation });
    expect(expertPtFailure(new Error('private native pointer and credential'))).toEqual({ reason: 'expert.solve_failed' });
  });
});
