import { describe, expect, it } from 'vitest';

import {
  PT_BACKEND_CAPABILITY_CONVENTION,
  PT_COMPONENT_INVENTORY_CONVENTION,
  PT_PHASE_TRANSITION_CONVENTION,
  ptOutcomeSeverity,
  ptServiceErrorLabel,
  validatePtFlashRequest,
  type PtBackendDescriptor,
  type PtFlashRequest,
} from './flash';

function backend(): PtBackendDescriptor {
  return {
    convention: 'PT/service-capability-discovery/v1',
    configuredBackendId: 'pr76-default',
    capability: {
      convention: PT_BACKEND_CAPABILITY_CONVENTION,
      backendId: 'PR76/pt-flash-backend/v1',
      modelProfile: 'PR76/default',
      algorithmProfile: 'PR76/max3/v1',
      publicationProfile: 'PT/phase-set-v1',
      configurationProfile: 'PR76/default/v1',
      datasetId: 'fixture/dataset',
      revision: 'fixture-r1',
      componentIds: ['methane', 'carbon-dioxide'],
      supportedPhaseCounts: [1, 2, 3],
      scalarSettings: [],
      transitionCapability: {
        convention: PT_PHASE_TRANSITION_CONVENTION,
        edges: [],
      },
      performsInitialStabilitySearch: true,
      performsFinalPhaseSetReview: true,
      performsBoundaryNeighborResolve: true,
      globalStabilityProven: false,
      phaseMetadataNamespace: '',
    },
    componentInventory: {
      convention: PT_COMPONENT_INVENTORY_CONVENTION,
      components: [
        { componentId: 'methane', feedIndex: 0 },
        { componentId: 'carbon-dioxide', feedIndex: 1 },
      ],
    },
  };
}

function request(overrides: Partial<PtFlashRequest> = {}): PtFlashRequest {
  return {
    configuredBackendId: 'pr76-default',
    pressurePa: 10e6,
    temperatureK: 350,
    feed: [
      { componentId: 'methane', moleFraction: 0.7 },
      { componentId: 'carbon-dioxide', moleFraction: 0.3 },
    ],
    ...overrides,
  };
}

describe('model-neutral PT request validation', () => {
  it('accepts the exact discovered inventory without changing values', () => {
    expect(validatePtFlashRequest(request(), backend())).toEqual([]);
  });

  it('accepts ID-keyed feed permutation because order is not identity', () => {
    expect(
      validatePtFlashRequest(
        request({
          feed: [
            { componentId: 'carbon-dioxide', moleFraction: 0.3 },
            { componentId: 'methane', moleFraction: 0.7 },
          ],
        }),
        backend(),
      ),
    ).toEqual([]);
  });

  it('does not normalize or silently repair an invalid composition sum', () => {
    const errors = validatePtFlashRequest(
      request({
        feed: [
          { componentId: 'methane', moleFraction: 0.69 },
          { componentId: 'carbon-dioxide', moleFraction: 0.3 },
        ],
      }),
      backend(),
    );
    expect(errors.some((error: string) => error.includes('must sum to 1'))).toBe(true);
  });

  it('rejects a stale backend selector and inventory drift', () => {
    const errors = validatePtFlashRequest(
      request({
        configuredBackendId: 'cpa-default',
        feed: [
          { componentId: 'methane', moleFraction: 0.5 },
          { componentId: 'water', moleFraction: 0.5 },
        ],
      }),
      backend(),
    );
    expect(errors).toEqual(
      expect.arrayContaining([
        'Selected backend does not match the discovered capability snapshot.',
        'Unknown component ID: water.',
        'Feed component IDs do not match the discovered inventory.',
      ]),
    );
  });

  it('rejects nonphysical pressure and temperature', () => {
    expect(
      validatePtFlashRequest(
        request({ pressurePa: 0, temperatureK: Number.NaN }),
        backend(),
      ),
    ).toEqual(
      expect.arrayContaining([
        'Pressure must be finite and greater than 0 Pa.',
        'Temperature must be finite and greater than 0 K.',
      ]),
    );
  });
});

describe('outcome presentation semantics', () => {
  it('keeps indeterminate in the scientific result severity', () => {
    expect(ptOutcomeSeverity('accepted')).toBe('accepted');
    expect(ptOutcomeSeverity('phase_set_unstable')).toBe('warning');
    expect(ptOutcomeSeverity('indeterminate')).toBe('indeterminate');
  });

  it('labels service failures independently of scientific outcomes', () => {
    expect(ptServiceErrorLabel('backend_execution_failure')).toBe(
      'Backend execution failure',
    );
  });
});
