import { beforeAll, describe, expect, it } from 'vitest';

import type { PtBackendDescriptor, PtFlashRequest } from '../domain/flash';
import { createGrpcWebFlashClient } from './flashClient';

const endpoint = process.env.MPMC_PT_GOLDEN_ENDPOINT;
const golden = endpoint === undefined ? describe.skip : describe;
const goldenEndpoint = endpoint ?? 'http://127.0.0.1:0';

function request(
  configuredBackendId: string,
  pressurePa = 5e6,
): PtFlashRequest {
  return {
    configuredBackendId,
    pressurePa,
    temperatureK: 325,
    // Intentionally reverse canonical inventory order. PtService, not the
    // browser or wire adapter, owns the stable-ID-to-backend-order mapping.
    feed: [
      { componentId: 'water', moleFraction: 0.6 },
      { componentId: 'methane', moleFraction: 0.4 },
    ],
  };
}

golden('TypeScript ↔ C++ PT process golden', () => {
  const client = createGrpcWebFlashClient(goldenEndpoint);
  let acceptedBackend: PtBackendDescriptor | undefined;
  let indeterminateBackend: PtBackendDescriptor | undefined;
  let discoveredBackendIds: string[] = [];

  beforeAll(async () => {
    const discovery = await client.discoverPtCapabilities({ timeoutMs: 2_000 });
    discoveredBackendIds = discovery.backends.map(
      (backend) => backend.configuredBackendId,
    );
    acceptedBackend = discovery.backends[0]!;
    indeterminateBackend = discovery.backends[1]!;
  });

  it('discovers the exact C++ backend and component inventory snapshots', () => {
    expect(discoveredBackendIds).toEqual([
      'golden.accepted',
      'golden.indeterminate',
    ]);
    expect(acceptedBackend).toMatchObject({
      convention: 'PT/service-capability-discovery/v1',
      configuredBackendId: 'golden.accepted',
      capability: {
        convention: 'PT/flash-backend/capability-and-dispatch/v1',
        backendId: 'golden/backend/v1',
        modelProfile: 'golden/model-neutral-fixture/v1',
        algorithmProfile: 'golden/scripted-no-eos/v1',
        publicationProfile: 'golden/wire-structure-only/v1',
        configurationProfile: 'golden/configuration/v1',
        datasetId: 'synthetic-wire-fixture/no-physical-data',
        revision: 'golden-r1',
        componentIds: ['methane', 'water'],
        supportedPhaseCounts: [1, 2, 3],
        scalarSettings: [
          { id: 'fixture-setting', value: 2.5, unit: 'dimensionless' },
        ],
        transitionCapability: {
          convention: 'PT/phase-transition-boundary/v1',
          edges: [
            {
              sourcePhaseCount: 1,
              targetPhaseCount: 2,
              support: 'fresh_target_resolve',
              requiresFreshTargetSolve: true,
            },
            {
              sourcePhaseCount: 2,
              targetPhaseCount: 1,
              support: 'detection_only',
              requiresFreshTargetSolve: true,
            },
            {
              sourcePhaseCount: 2,
              targetPhaseCount: 3,
              support: 'fresh_target_resolve',
              requiresFreshTargetSolve: true,
            },
            {
              sourcePhaseCount: 3,
              targetPhaseCount: 2,
              support: 'fresh_target_resolve',
              requiresFreshTargetSolve: true,
            },
          ],
        },
        performsInitialStabilitySearch: true,
        performsFinalPhaseSetReview: true,
        performsBoundaryNeighborResolve: true,
        globalStabilityProven: false,
        phaseMetadataNamespace: 'golden/provider-phase-metadata/v1',
      },
      componentInventory: {
        convention: 'PT/service-component-inventory/v1',
        components: [
          { componentId: 'methane', feedIndex: 0 },
          { componentId: 'water', feedIndex: 1 },
        ],
      },
    });
  });

  it('round-trips one accepted solve with canonical feed and provenance', async () => {
    const response = await client.solvePtFlash(
      request('golden.accepted'),
      acceptedBackend!,
      { timeoutMs: 2_000 },
    );
    expect(response.kind).toBe('result');
    if (response.kind !== 'result') {
      throw new Error('golden accepted call returned service error');
    }
    expect(response.result).toMatchObject({
      outcome: 'accepted',
      pressurePa: 5e6,
      temperatureK: 325,
      feed: [
        { componentId: 'methane', moleFraction: 0.4 },
        { componentId: 'water', moleFraction: 0.6 },
      ],
      globalStabilityProven: false,
      morphologyResolved: false,
      diagnostic: 'synthetic accepted result for cross-language wire golden',
      phases: [
        {
          phaseIndex: 0,
          molePhaseFraction: 0.25,
          components: [
            {
              componentId: 'methane',
              moleFraction: 0.4,
              lnFugacityCoefficient: -0.125,
            },
            {
              componentId: 'water',
              moleFraction: 0.6,
              lnFugacityCoefficient: -0.25,
            },
          ],
          providerBranch: '7',
          providerBranchSmooth: true,
          compressibilityFactor: 0.8,
          providerMetadata: {
            roleId: 'opaque-role-0',
            familyId: 'opaque-family',
          },
        },
        {
          phaseIndex: 1,
          molePhaseFraction: 0.75,
          components: [
            {
              componentId: 'methane',
              moleFraction: 0.4,
              lnFugacityCoefficient: -0.125,
            },
            {
              componentId: 'water',
              moleFraction: 0.6,
              lnFugacityCoefficient: -0.25,
            },
          ],
          providerBranch: '9',
          providerBranchSmooth: true,
          compressibilityFactor: 0.9,
          providerMetadata: {
            roleId: 'opaque-role-1',
            familyId: 'opaque-family',
          },
        },
      ],
    });
    expect(response.result.provenance.backend).toEqual(acceptedBackend!);
    expect(response.result.provenance).toMatchObject({
      backendResultConvention: 'PT/flash-backend/result/v1',
      phaseSetConvention: 'PT/phase-set-v1',
      phaseTransitionConvention: 'PT/phase-transition-boundary/v1',
      providerResultConvention: 'golden/provider-result/v1',
    });
    expect(response.result.transitionReport.evidence).toEqual([
      {
        sourcePhaseCount: 3,
        targetPhaseCount: 2,
        trigger: 'phase_disappearance',
        resolution: 'accepted_target',
        freshTargetSolveAttempted: true,
        targetTopologyClosed: true,
        providerEvidenceProfile: 'golden/transition-evidence/v1',
        diagnostic: 'synthetic accepted-target evidence for wire mapping only',
      },
    ]);
  });

  it('keeps indeterminate as a phase-free result arm', async () => {
    const response = await client.solvePtFlash(
      request('golden.indeterminate'),
      indeterminateBackend!,
      { timeoutMs: 2_000 },
    );
    expect(response).toMatchObject({
      kind: 'result',
      result: {
        outcome: 'indeterminate',
        phases: [],
        transitionReport: {
          evidence: [
            {
              sourcePhaseCount: 2,
              targetPhaseCount: 3,
              trigger: 'final_phase_set_instability',
              resolution: 'target_resolve_failed',
              freshTargetSolveAttempted: true,
              targetTopologyClosed: false,
              providerEvidenceProfile: 'golden/transition-evidence/v1',
              diagnostic:
                'synthetic failed-target evidence for wire mapping only',
            },
          ],
        },
        diagnostic:
          'synthetic indeterminate result for cross-language wire golden',
      },
    });
  });

  it('keeps a PtService validation error on the protobuf error arm', async () => {
    const response = await client.solvePtFlash(
      request('golden.accepted', -1),
      acceptedBackend!,
      { timeoutMs: 2_000 },
    );
    expect(response).toEqual({
      kind: 'service_error',
      error: {
        code: 'invalid_pressure',
        field: 'pressure_pa',
        diagnostic: 'pressure must be finite and greater than 0 Pa',
      },
    });
  });

  it('answers an allowed-origin CORS preflight without reaching gRPC', async () => {
    const rpcUrl = new URL(
      '/mpmc.runtime.v1.PtFlashService/DiscoverPtCapabilities',
      goldenEndpoint,
    );
    const response = await fetch(rpcUrl, {
      method: 'OPTIONS',
      headers: {
        Origin: 'http://localhost:5173',
        'Access-Control-Request-Method': 'POST',
        'Access-Control-Request-Headers':
          'content-type,x-grpc-web,grpc-timeout,x-user-agent',
      },
    });
    expect(response.ok).toBe(true);
    expect(response.headers.get('access-control-allow-origin')).toBe(
      'http://localhost:5173',
    );
    expect(response.headers.get('access-control-allow-methods')).toContain('POST');
    expect(response.headers.get('access-control-allow-headers')).toContain(
      'x-grpc-web',
    );
  });

  it('does not authorize an origin outside the CORS allowlist', async () => {
    const rpcUrl = new URL(
      '/mpmc.runtime.v1.PtFlashService/DiscoverPtCapabilities',
      goldenEndpoint,
    );
    const response = await fetch(rpcUrl, {
      method: 'OPTIONS',
      headers: {
        Origin: 'https://not-allowed.example',
        'Access-Control-Request-Method': 'POST',
        'Access-Control-Request-Headers': 'content-type,x-grpc-web',
      },
    });
    expect(response.headers.has('access-control-allow-origin')).toBe(false);
  });
});
