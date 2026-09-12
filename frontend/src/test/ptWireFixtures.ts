import { create } from '@bufbuild/protobuf';

import {
  DiscoverPtCapabilitiesResponseSchema,
  PtBackendDescriptorSchema,
  PtComputationOutcome,
  PtServiceErrorCode,
  SolvePtFlashResponseSchema,
} from '../gen/mpmc/runtime/v1/pt_service_pb';
import {
  PT_BACKEND_CAPABILITY_CONVENTION,
  PT_BACKEND_RESULT_CONVENTION,
  PT_CAPABILITY_CONVENTION,
  PT_COMPONENT_INVENTORY_CONVENTION,
  PT_PHASE_SET_CONVENTION,
  PT_PHASE_TRANSITION_CONVENTION,
  PT_SERVICE_BOUNDARY_CONVENTION,
  PT_SERVICE_RESULT_CONVENTION,
  PT_WIRE_CONTRACT,
  type PtFlashRequest,
} from '../domain/flash';

export function domainRequest(): PtFlashRequest {
  return {
    configuredBackendId: 'pr76-default',
    pressurePa: 10e6,
    temperatureK: 350,
    feed: [
      { componentId: 'methane', moleFraction: 0.7 },
      { componentId: 'carbon-dioxide', moleFraction: 0.3 },
    ],
  };
}

export function wireBackend() {
  return create(PtBackendDescriptorSchema, {
    convention: PT_CAPABILITY_CONVENTION,
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
  });
}

export function wireDiscovery() {
  return create(DiscoverPtCapabilitiesResponseSchema, {
    wireContract: PT_WIRE_CONTRACT,
    serviceBoundaryConvention: PT_SERVICE_BOUNDARY_CONVENTION,
    capabilityConvention: PT_CAPABILITY_CONVENTION,
    backends: [wireBackend()],
  });
}

export function wireResult(
  outcome: PtComputationOutcome = PtComputationOutcome.ACCEPTED,
) {
  const phases =
    outcome === PtComputationOutcome.ACCEPTED
      ? [
          {
            phaseIndex: 0,
            molePhaseFraction: 1,
            components: [
              {
                componentId: 'methane',
                moleFraction: 0.7,
                lnFugacityCoefficient: -0.1,
              },
              {
                componentId: 'carbon-dioxide',
                moleFraction: 0.3,
                lnFugacityCoefficient: -0.2,
              },
            ],
            providerBranch: 0n,
            providerBranchSmooth: true,
            compressibilityFactor: 0.8,
          },
        ]
      : [];
  return create(SolvePtFlashResponseSchema, {
    wireContract: PT_WIRE_CONTRACT,
    serviceResultConvention: PT_SERVICE_RESULT_CONVENTION,
    payload: {
      case: 'result',
      value: {
        outcome,
        provenance: {
          backend: wireBackend(),
          backendResultConvention: PT_BACKEND_RESULT_CONVENTION,
          phaseSetConvention: PT_PHASE_SET_CONVENTION,
          phaseTransitionConvention: PT_PHASE_TRANSITION_CONVENTION,
          providerResultConvention: 'PR76/provider-result/v1',
        },
        pressurePa: 10e6,
        temperatureK: 350,
        feed: [
          { componentId: 'methane', moleFraction: 0.7 },
          { componentId: 'carbon-dioxide', moleFraction: 0.3 },
        ],
        phases,
        transitionReport: {
          convention: PT_PHASE_TRANSITION_CONVENTION,
          evidence: [],
        },
        globalStabilityProven: false,
        morphologyResolved: false,
        diagnostic: outcome === PtComputationOutcome.ACCEPTED ? 'accepted' : 'not accepted',
      },
    },
  });
}

export function wireServiceError() {
  return create(SolvePtFlashResponseSchema, {
    wireContract: PT_WIRE_CONTRACT,
    serviceResultConvention: PT_SERVICE_RESULT_CONVENTION,
    payload: {
      case: 'error',
      value: {
        code: PtServiceErrorCode.INVALID_PRESSURE,
        field: 'pressure_pa',
        diagnostic: 'pressure must be positive',
      },
    },
  });
}
