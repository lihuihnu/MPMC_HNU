import type { DescMessage, JsonObject } from '@bufbuild/protobuf';

import {
  ModelConfigurationLimitsSchema,
  PtEosRootSettingsSchema,
  PtSolverSafetyLimitsSchema,
  PtStabilitySettingsSchema,
  PtThreePhaseSettingsSchema,
  PtTwoPhaseSettingsSchema,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';
import {
  PT_BACKEND_CAPABILITY_CONVENTION,
  PT_BACKEND_RESULT_CONVENTION,
  PT_PHASE_SET_CONVENTION,
  PT_PHASE_TRANSITION_CONVENTION,
} from '../domain/flash';

function block(schema: DescMessage): JsonObject {
  return Object.fromEntries(
    schema.fields.map((field) => [field.jsonName, field.localName === 'automaticMultistart' ? false : 1]),
  );
}

export function expertSnapshotJson(): JsonObject {
  const provenance = {
    kind: 'SOURCE_KIND_LITERATURE',
    reference: 'Doe 2026',
    revision: 'source-r1',
    locator: 'table-2',
    note: 'inspection fixture only',
    acquisition: 'repository fixture',
    usageTerms: 'test-only',
  };
  const scalar = (value: number, originalUnit = '') => ({
    value,
    provenance,
    ...(originalUnit ? { originalUnit, conversion: `converted from ${originalUnit}` } : {}),
  });
  const capability = {
    convention: PT_BACKEND_CAPABILITY_CONVENTION,
    datasetId: 'expert-inspector-fixture',
    revision: 'r1',
    publicationProfile: 'expert-inspector-publication/v1',
    componentIds: ['methane', 'ethane'],
    supportedPhaseCounts: [1, 2, 3],
    transitionCapability: { convention: PT_PHASE_TRANSITION_CONVENTION },
    globalStabilityProven: false,
  };
  return {
    definition: {
      version: 'thermodynamic-model/parameter-definition/v1',
      family: 'MODEL_FAMILY_PR76',
      displayName: 'Read-only PR76 fixture',
      datasetId: 'expert-inspector-fixture',
      revision: 'r1',
      provenance,
      components: [
        {
          componentId: 'methane', displayName: 'Methane', kind: 'COMPONENT_KIND_PURE', provenance,
          molarMassKgPerMol: scalar(0.016043, 'kg/mol'),
        },
        {
          componentId: 'ethane', displayName: 'Ethane', kind: 'COMPONENT_KIND_PURE', provenance,
          molarMassKgPerMol: scalar(0.03007, 'kg/mol'),
        },
      ],
      applicability: {
        temperatureLowerK: 250,
        temperatureLowerExclusive: true,
        pressureUpperPa: 30_000_000,
        pressureUpperExclusive: false,
        provenance,
      },
      pr76: {
        pure: [
          {
            componentId: 'methane', criticalTemperatureK: scalar(190.564),
            criticalPressurePa: scalar(4_599_200), acentricFactor: scalar(0.01142),
          },
          {
            componentId: 'ethane', criticalTemperatureK: scalar(305.322),
            criticalPressurePa: scalar(4_872_200), acentricFactor: scalar(0.0995),
          },
        ],
        binary: [
          {
            firstComponentId: 'methane', secondComponentId: 'ethane',
            kij: scalar(0.0123),
          },
        ],
      },
    },
    capability,
    settings: {
      version: 'pt-solver-settings/v1',
      kind: 'SOLVER_SETTINGS_KIND_PRESET',
      presetId: 'mpmc-balanced-default/v1',
      eosRoot: { ...block(PtEosRootSettingsSchema), maxIterations: '9007199254740993' },
      initialStability: block(PtStabilitySettingsSchema),
      twoPhase: block(PtTwoPhaseSettingsSchema),
      finalTwoPhaseStability: block(PtStabilitySettingsSchema),
      threePhase: block(PtThreePhaseSettingsSchema),
      finalThreePhaseStability: block(PtStabilitySettingsSchema),
    },
    parameterLimits: block(ModelConfigurationLimitsSchema),
    solverLimits: block(PtSolverSafetyLimitsSchema),
  };
}

export function expertResultJson(): JsonObject {
  const snapshot = expertSnapshotJson();
  return {
    backendResultConvention: PT_BACKEND_RESULT_CONVENTION,
    phaseSetConvention: PT_PHASE_SET_CONVENTION,
    capability: snapshot.capability as JsonObject,
    outcome: 'PT_COMPUTATION_OUTCOME_INDETERMINATE',
    maximumPhaseCount: 3,
    pressurePa: 10_000_000,
    temperatureK: 330,
    feed: [0.7, 0.3],
    candidatePhaseSet: {
      phases: [
        {
          molePhaseFraction: 1,
          composition: [0.7, 0.3],
          lnFugacityCoefficient: [0.1, 0.2],
          providerBranch: '7',
          providerBranchSmooth: false,
          compressibilityFactor: 0.8,
        },
      ],
    },
    globalStabilityProven: false,
    diagnostic: 'Opaque backend diagnostic preserved verbatim.',
    transitionReport: {
      convention: PT_PHASE_TRANSITION_CONVENTION,
      evidence: [
        {
          sourcePhaseCount: 2,
          targetPhaseCount: 1,
          trigger: 'PT_PHASE_TRANSITION_TRIGGER_PHASE_DISAPPEARANCE',
          resolution: 'PT_PHASE_TRANSITION_RESOLUTION_INDETERMINATE',
          freshTargetSolveAttempted: true,
          targetTopologyClosed: false,
          providerEvidenceProfile: 'expert-fixture',
          diagnostic: 'transition evidence',
        },
      ],
    },
    providerResultConvention: 'expert-inspector-publication/v1',
    phaseMetadata: [],
    morphologyResolved: false,
  };
}
