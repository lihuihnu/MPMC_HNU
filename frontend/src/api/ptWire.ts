import type { MessageInitShape } from '@bufbuild/protobuf';

import {
  PtComputationOutcome as WireComputationOutcome,
  PtPhaseTransitionResolution as WireTransitionResolution,
  PtPhaseTransitionSupport as WireTransitionSupport,
  PtPhaseTransitionTrigger as WireTransitionTrigger,
  PtServiceErrorCode as WireServiceErrorCode,
  SolvePtFlashRequestSchema,
  type DiscoverPtCapabilitiesResponse as WireDiscoveryResponse,
  type PtBackendDescriptor as WireBackendDescriptor,
  type PtComputationResult as WireComputationResult,
  type PtPhaseTransitionReport as WireTransitionReport,
  type PtServiceError as WireServiceError,
  type SolvePtFlashResponse as WireSolveResponse,
} from '../gen/mpmc/runtime/v1/pt_service_pb';
import {
  PT_BACKEND_CAPABILITY_CONVENTION,
  PT_BACKEND_RESULT_CONVENTION,
  PT_CAPABILITY_CONVENTION,
  PT_COMPONENT_INVENTORY_CONVENTION,
  PT_PHASE_SET_CONVENTION,
  PT_PHASE_TRANSITION_CONVENTION,
  PT_SERVICE_BOUNDARY_CONVENTION,
  PT_SERVICE_REQUEST_CONVENTION,
  PT_SERVICE_RESULT_CONVENTION,
  PT_WIRE_CONTRACT,
  type PtBackendDescriptor,
  type PtCapabilityDiscovery,
  type PtComputationOutcome,
  type PtComputationResult,
  type PtFlashRequest,
  type PtFlashResponse,
  type PtPhaseTransitionResolution,
  type PtPhaseTransitionSupport,
  type PtPhaseTransitionTrigger,
  type PtServiceError,
  type PtServiceErrorCode,
} from '../domain/flash';

const MAX_BACKENDS = 64;
const MAX_COMPONENTS_PER_BACKEND = 256;
const MAX_PHASES = 64;
const MAX_SUPPORTED_PHASE_COUNTS = 64;
const MAX_SCALAR_SETTINGS = 256;
const MAX_TRANSITION_ITEMS = 4096;
const MAX_IDENTIFIER_BYTES = 256;
const MAX_DISPLAY_STRING_BYTES = 64 * 1024;

const utf8 = new TextEncoder();

export class PtWireContractError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'PtWireContractError';
  }
}

function fail(path: string, message: string): never {
  throw new PtWireContractError(`${path}: ${message}`);
}

function required<T>(value: T | undefined, path: string): T {
  return value === undefined ? fail(path, 'required field is absent') : value;
}

function boundedString(
  value: string | undefined,
  path: string,
  allowEmpty = false,
): string {
  const result = required(value, path);
  if (
    (!allowEmpty && result.length === 0) ||
    utf8.encode(result).length > MAX_DISPLAY_STRING_BYTES
  ) {
    fail(path, allowEmpty ? 'field exceeds the client size limit' : 'field is empty or over limit');
  }
  return result;
}

function identifier(value: string | undefined, path: string): string {
  const result = boundedString(value, path);
  if (
    utf8.encode(result).length > MAX_IDENTIFIER_BYTES ||
    !/^[\x21-\x7e]+$/.test(result)
  ) {
    fail(path, 'identifier must be visible ASCII within 256 UTF-8 bytes');
  }
  return result;
}

function exactConvention<T extends string>(
  value: string | undefined,
  expected: T,
  path: string,
): T {
  if (value !== expected) {
    fail(path, `expected ${expected}`);
  }
  return expected;
}

function finite(value: number | undefined, path: string): number {
  const result = required(value, path);
  return Number.isFinite(result) ? result : fail(path, 'value must be finite');
}

function positive(value: number | undefined, path: string): number {
  const result = finite(value, path);
  return result > 0 ? result : fail(path, 'value must be greater than zero');
}

function fraction(value: number | undefined, path: string): number {
  const result = finite(value, path);
  return result >= 0 && result <= 1
    ? result
    : fail(path, 'fraction must be in [0,1]');
}

function uint32(value: number | undefined, path: string): number {
  const result = required(value, path);
  return Number.isSafeInteger(result) && result >= 0 && result <= 0xffff_ffff
    ? result
    : fail(path, 'value must be a uint32');
}

function unique<T extends string | number>(values: readonly T[], path: string): void {
  if (new Set(values).size !== values.length) {
    fail(path, 'values must be unique');
  }
}

function normalized(values: readonly number[], path: string): void {
  let sum = 0;
  let correction = 0;
  for (const value of values) {
    const increment = value - correction;
    const next = sum + increment;
    correction = (next - sum) - increment;
    sum = next;
  }
  if (Math.abs(sum - 1) > 64 * Number.EPSILON) {
    fail(path, 'fractions do not sum to one within the service roundoff contract');
  }
}

function mapTransitionSupport(
  value: WireTransitionSupport | undefined,
  path: string,
): PtPhaseTransitionSupport {
  switch (required(value, path)) {
    case WireTransitionSupport.DETECTION_ONLY:
      return 'detection_only';
    case WireTransitionSupport.FRESH_TARGET_RESOLVE:
      return 'fresh_target_resolve';
    case WireTransitionSupport.UNSPECIFIED:
      return fail(path, 'unspecified transition support');
    default:
      return fail(path, 'unknown transition support');
  }
}

function mapTransitionTrigger(
  value: WireTransitionTrigger | undefined,
  path: string,
): PtPhaseTransitionTrigger {
  switch (required(value, path)) {
    case WireTransitionTrigger.INITIAL_STABILITY_WITNESS:
      return 'initial_stability_witness';
    case WireTransitionTrigger.FINAL_PHASE_SET_INSTABILITY:
      return 'final_phase_set_instability';
    case WireTransitionTrigger.PHASE_DISAPPEARANCE:
      return 'phase_disappearance';
    case WireTransitionTrigger.PROVIDER_TOPOLOGY_WITNESS:
      return 'provider_topology_witness';
    case WireTransitionTrigger.PROVIDER_BOUNDARY_ROUTE:
      return 'provider_boundary_route';
    case WireTransitionTrigger.UNSPECIFIED:
      return fail(path, 'unspecified transition trigger');
    default:
      return fail(path, 'unknown transition trigger');
  }
}

function mapTransitionResolution(
  value: WireTransitionResolution | undefined,
  path: string,
): PtPhaseTransitionResolution {
  switch (required(value, path)) {
    case WireTransitionResolution.ACCEPTED_TARGET:
      return 'accepted_target';
    case WireTransitionResolution.TARGET_RESOLVE_REQUIRED:
      return 'target_resolve_required';
    case WireTransitionResolution.TARGET_RESOLVE_FAILED:
      return 'target_resolve_failed';
    case WireTransitionResolution.CANDIDATE_NOT_ACCEPTED:
      return 'candidate_not_accepted';
    case WireTransitionResolution.BROADER_TOPOLOGY_REQUIRED:
      return 'broader_topology_required';
    case WireTransitionResolution.INDETERMINATE:
      return 'indeterminate';
    case WireTransitionResolution.UNSPECIFIED:
      return fail(path, 'unspecified transition resolution');
    default:
      return fail(path, 'unknown transition resolution');
  }
}

function mapOutcome(
  value: WireComputationOutcome | undefined,
  path: string,
): PtComputationOutcome {
  switch (required(value, path)) {
    case WireComputationOutcome.ACCEPTED:
      return 'accepted';
    case WireComputationOutcome.PHASE_SET_UNSTABLE:
      return 'phase_set_unstable';
    case WireComputationOutcome.INDETERMINATE:
      return 'indeterminate';
    case WireComputationOutcome.UNSPECIFIED:
      return fail(path, 'unspecified computation outcome');
    default:
      return fail(path, 'unknown computation outcome');
  }
}

function mapServiceErrorCode(
  value: WireServiceErrorCode | undefined,
  path: string,
): PtServiceErrorCode {
  switch (required(value, path)) {
    case WireServiceErrorCode.INVALID_CONFIGURED_BACKEND_ID:
      return 'invalid_configured_backend_id';
    case WireServiceErrorCode.CONFIGURED_BACKEND_NOT_FOUND:
      return 'configured_backend_not_found';
    case WireServiceErrorCode.INVALID_PRESSURE:
      return 'invalid_pressure';
    case WireServiceErrorCode.INVALID_TEMPERATURE:
      return 'invalid_temperature';
    case WireServiceErrorCode.COMPONENT_COUNT_MISMATCH:
      return 'component_count_mismatch';
    case WireServiceErrorCode.INVALID_COMPONENT_ID:
      return 'invalid_component_id';
    case WireServiceErrorCode.DUPLICATE_COMPONENT:
      return 'duplicate_component';
    case WireServiceErrorCode.UNKNOWN_COMPONENT:
      return 'unknown_component';
    case WireServiceErrorCode.INVALID_MOLE_FRACTION:
      return 'invalid_mole_fraction';
    case WireServiceErrorCode.COMPOSITION_NOT_NORMALIZED:
      return 'composition_not_normalized';
    case WireServiceErrorCode.BACKEND_REJECTED_REQUEST:
      return 'backend_rejected_request';
    case WireServiceErrorCode.BACKEND_CONTRACT_VIOLATION:
      return 'backend_contract_violation';
    case WireServiceErrorCode.BACKEND_EXECUTION_FAILURE:
      return 'backend_execution_failure';
    case WireServiceErrorCode.UNSPECIFIED:
      return fail(path, 'unspecified service error code');
    default:
      return fail(path, 'unknown service error code');
  }
}

function mapBackendDescriptor(
  source: WireBackendDescriptor,
  path: string,
): PtBackendDescriptor {
  const configuredBackendId = identifier(
    source.configuredBackendId,
    `${path}.configured_backend_id`,
  );
  const capability = required(source.capability, `${path}.capability`);
  const componentInventory = required(
    source.componentInventory,
    `${path}.component_inventory`,
  );

  const componentIds = capability.componentIds.map((value, index) =>
    identifier(value, `${path}.capability.component_ids[${index}]`),
  );
  if (componentIds.length === 0 || componentIds.length > MAX_COMPONENTS_PER_BACKEND) {
    fail(`${path}.capability.component_ids`, 'component count is outside client limits');
  }
  unique(componentIds, `${path}.capability.component_ids`);

  const supportedPhaseCounts = capability.supportedPhaseCounts.map((value, index) => {
    const count = uint32(value, `${path}.capability.supported_phase_counts[${index}]`);
    return count > 0 ? count : fail(path, 'supported phase counts must be positive');
  });
  if (supportedPhaseCounts.length === 0) {
    fail(`${path}.capability.supported_phase_counts`, 'at least one phase count is required');
  }
  if (
    supportedPhaseCounts.length > MAX_SUPPORTED_PHASE_COUNTS ||
    supportedPhaseCounts.some((count) => count > MAX_PHASES)
  ) {
    fail(`${path}.capability.supported_phase_counts`, 'phase counts exceed client limits');
  }
  unique(supportedPhaseCounts, `${path}.capability.supported_phase_counts`);
  const maximumPhaseCount = Math.max(...supportedPhaseCounts);

  const transition = required(
    capability.transitionCapability,
    `${path}.capability.transition_capability`,
  );
  exactConvention(
    transition.convention,
    PT_PHASE_TRANSITION_CONVENTION,
    `${path}.capability.transition_capability.convention`,
  );
  if (transition.edges.length > MAX_TRANSITION_ITEMS) {
    fail(`${path}.capability.transition_capability.edges`, 'edge count exceeds client limits');
  }
  const edgeKeys = new Set<string>();
  const edges = transition.edges.map((edge, index) => {
    const edgePath = `${path}.capability.transition_capability.edges[${index}]`;
    const sourcePhaseCount = uint32(edge.sourcePhaseCount, `${edgePath}.source_phase_count`);
    const targetPhaseCount = uint32(edge.targetPhaseCount, `${edgePath}.target_phase_count`);
    if (
      sourcePhaseCount === 0 ||
      targetPhaseCount === 0 ||
      sourcePhaseCount === targetPhaseCount ||
      sourcePhaseCount > maximumPhaseCount ||
      targetPhaseCount > maximumPhaseCount
    ) {
      fail(edgePath, 'invalid phase-count transition edge');
    }
    const support = mapTransitionSupport(edge.support, `${edgePath}.support`);
    const requiresFreshTargetSolve = required(
      edge.requiresFreshTargetSolve,
      `${edgePath}.requires_fresh_target_solve`,
    );
    if (support === 'fresh_target_resolve' && !requiresFreshTargetSolve) {
      fail(edgePath, 'fresh-target support must require a fresh target solve');
    }
    const key = `${sourcePhaseCount}:${targetPhaseCount}`;
    if (edgeKeys.has(key)) {
      fail(edgePath, 'duplicate transition edge');
    }
    edgeKeys.add(key);
    return { sourcePhaseCount, targetPhaseCount, support, requiresFreshTargetSolve };
  });

  const scalarIds = new Set<string>();
  if (capability.scalarSettings.length > MAX_SCALAR_SETTINGS) {
    fail(`${path}.capability.scalar_settings`, 'setting count exceeds client limits');
  }
  const scalarSettings = capability.scalarSettings.map((setting, index) => {
    const settingPath = `${path}.capability.scalar_settings[${index}]`;
    const id = boundedString(setting.id, `${settingPath}.id`);
    if (scalarIds.has(id)) {
      fail(`${settingPath}.id`, 'duplicate scalar setting ID');
    }
    scalarIds.add(id);
    return {
      id,
      value: finite(setting.value, `${settingPath}.value`),
      unit: boundedString(setting.unit, `${settingPath}.unit`),
    };
  });

  exactConvention(
    componentInventory.convention,
    PT_COMPONENT_INVENTORY_CONVENTION,
    `${path}.component_inventory.convention`,
  );
  if (componentInventory.components.length !== componentIds.length) {
    fail(`${path}.component_inventory`, 'inventory does not match capability component count');
  }
  const components = componentInventory.components.map((component, index) => {
    const componentPath = `${path}.component_inventory.components[${index}]`;
    const componentId = identifier(component.componentId, `${componentPath}.component_id`);
    const feedIndex = uint32(component.feedIndex, `${componentPath}.feed_index`);
    if (componentId !== componentIds[index] || feedIndex !== index) {
      fail(componentPath, 'inventory ID/order does not match the backend capability');
    }
    return { componentId, feedIndex };
  });

  return {
    convention: exactConvention(
      source.convention,
      PT_CAPABILITY_CONVENTION,
      `${path}.convention`,
    ),
    configuredBackendId,
    capability: {
      convention: exactConvention(
        capability.convention,
        PT_BACKEND_CAPABILITY_CONVENTION,
        `${path}.capability.convention`,
      ),
      backendId: identifier(capability.backendId, `${path}.capability.backend_id`),
      modelProfile: boundedString(capability.modelProfile, `${path}.capability.model_profile`),
      algorithmProfile: boundedString(
        capability.algorithmProfile,
        `${path}.capability.algorithm_profile`,
      ),
      publicationProfile: boundedString(
        capability.publicationProfile,
        `${path}.capability.publication_profile`,
      ),
      configurationProfile: boundedString(
        capability.configurationProfile,
        `${path}.capability.configuration_profile`,
      ),
      datasetId: boundedString(capability.datasetId, `${path}.capability.dataset_id`),
      revision: boundedString(capability.revision, `${path}.capability.revision`),
      componentIds,
      supportedPhaseCounts,
      scalarSettings,
      transitionCapability: {
        convention: PT_PHASE_TRANSITION_CONVENTION,
        edges,
      },
      performsInitialStabilitySearch: required(
        capability.performsInitialStabilitySearch,
        `${path}.capability.performs_initial_stability_search`,
      ),
      performsFinalPhaseSetReview: required(
        capability.performsFinalPhaseSetReview,
        `${path}.capability.performs_final_phase_set_review`,
      ),
      performsBoundaryNeighborResolve: required(
        capability.performsBoundaryNeighborResolve,
        `${path}.capability.performs_boundary_neighbor_resolve`,
      ),
      globalStabilityProven: required(
        capability.globalStabilityProven,
        `${path}.capability.global_stability_proven`,
      ),
      phaseMetadataNamespace: boundedString(
        capability.phaseMetadataNamespace,
        `${path}.capability.phase_metadata_namespace`,
        true,
      ),
    },
    componentInventory: {
      convention: PT_COMPONENT_INVENTORY_CONVENTION,
      components,
    },
  };
}

export function mapDiscoveryResponse(source: WireDiscoveryResponse): PtCapabilityDiscovery {
  exactConvention(source.wireContract, PT_WIRE_CONTRACT, 'discovery.wire_contract');
  exactConvention(
    source.serviceBoundaryConvention,
    PT_SERVICE_BOUNDARY_CONVENTION,
    'discovery.service_boundary_convention',
  );
  exactConvention(
    source.capabilityConvention,
    PT_CAPABILITY_CONVENTION,
    'discovery.capability_convention',
  );
  if (source.backends.length > MAX_BACKENDS) {
    fail('discovery.backends', 'backend count exceeds the client limit');
  }
  const backends = source.backends.map((backend, index) =>
    mapBackendDescriptor(backend, `discovery.backends[${index}]`),
  );
  unique(
    backends.map((backend) => backend.configuredBackendId),
    'discovery.backends.configured_backend_id',
  );
  return {
    wireContract: PT_WIRE_CONTRACT,
    serviceBoundaryConvention: PT_SERVICE_BOUNDARY_CONVENTION,
    capabilityConvention: PT_CAPABILITY_CONVENTION,
    backends,
  };
}

export function toWireSolveRequest(
  request: PtFlashRequest,
): MessageInitShape<typeof SolvePtFlashRequestSchema> {
  return {
    wireContract: PT_WIRE_CONTRACT,
    serviceRequestConvention: PT_SERVICE_REQUEST_CONVENTION,
    configuredBackendId: request.configuredBackendId,
    pressurePa: request.pressurePa,
    temperatureK: request.temperatureK,
    feed: request.feed.map((component) => ({
      componentId: component.componentId,
      moleFraction: component.moleFraction,
    })),
  };
}

function mapTransitionReport(
  source: WireTransitionReport,
  path: string,
  backend: PtBackendDescriptor,
): PtComputationResult['transitionReport'] {
  exactConvention(source.convention, PT_PHASE_TRANSITION_CONVENTION, `${path}.convention`);
  const maximumPhaseCount = Math.max(...backend.capability.supportedPhaseCounts);
  if (source.evidence.length > MAX_TRANSITION_ITEMS) {
    fail(`${path}.evidence`, 'evidence count exceeds client limits');
  }
  return {
    convention: PT_PHASE_TRANSITION_CONVENTION,
    evidence: source.evidence.map((item, index) => {
      const evidencePath = `${path}.evidence[${index}]`;
      const sourcePhaseCount = uint32(
        item.sourcePhaseCount,
        `${evidencePath}.source_phase_count`,
      );
      if (sourcePhaseCount === 0 || sourcePhaseCount > maximumPhaseCount) {
        fail(`${evidencePath}.source_phase_count`, 'phase count is outside capability');
      }
      const targetPhaseCount =
        item.targetPhaseCount === undefined
          ? undefined
          : uint32(item.targetPhaseCount, `${evidencePath}.target_phase_count`);
      const declaredEdge =
        targetPhaseCount === undefined
          ? undefined
          : backend.capability.transitionCapability.edges.find(
              (edge) =>
                edge.sourcePhaseCount === sourcePhaseCount &&
                edge.targetPhaseCount === targetPhaseCount,
            );
      if (
        targetPhaseCount !== undefined &&
        (targetPhaseCount === 0 ||
          targetPhaseCount > maximumPhaseCount ||
          targetPhaseCount === sourcePhaseCount ||
          declaredEdge === undefined)
      ) {
        fail(`${evidencePath}.target_phase_count`, 'target has no declared capability edge');
      }
      const resolution = mapTransitionResolution(
        item.resolution,
        `${evidencePath}.resolution`,
      );
      const freshTargetSolveAttempted = required(
        item.freshTargetSolveAttempted,
        `${evidencePath}.fresh_target_solve_attempted`,
      );
      const targetTopologyClosed = required(
        item.targetTopologyClosed,
        `${evidencePath}.target_topology_closed`,
      );
      switch (resolution) {
        case 'accepted_target':
          if (
            targetPhaseCount === undefined ||
            !freshTargetSolveAttempted ||
            !targetTopologyClosed ||
            declaredEdge?.support !== 'fresh_target_resolve' ||
            !declaredEdge.requiresFreshTargetSolve
          ) {
            fail(evidencePath, 'accepted target lacks a closed fresh-target resolve');
          }
          break;
        case 'target_resolve_required':
          if (
            targetPhaseCount === undefined ||
            freshTargetSolveAttempted ||
            targetTopologyClosed ||
            declaredEdge === undefined ||
            !declaredEdge.requiresFreshTargetSolve
          ) {
            fail(evidencePath, 'target-resolve-required evidence is inconsistent');
          }
          break;
        case 'target_resolve_failed':
          if (!freshTargetSolveAttempted || targetTopologyClosed) {
            fail(evidencePath, 'failed target resolve must be attempted and unclosed');
          }
          break;
        case 'candidate_not_accepted':
        case 'broader_topology_required':
        case 'indeterminate':
          if (targetTopologyClosed) {
            fail(evidencePath, 'non-accepted transition evidence cannot close a target');
          }
          break;
      }
      const base = {
        sourcePhaseCount,
        trigger: mapTransitionTrigger(item.trigger, `${evidencePath}.trigger`),
        resolution,
        freshTargetSolveAttempted,
        targetTopologyClosed,
        providerEvidenceProfile: boundedString(
          item.providerEvidenceProfile,
          `${evidencePath}.provider_evidence_profile`,
        ),
        diagnostic: boundedString(item.diagnostic, `${evidencePath}.diagnostic`, true),
      };
      return targetPhaseCount === undefined ? base : { ...base, targetPhaseCount };
    }),
  };
}

function sameBackendSnapshot(left: PtBackendDescriptor, right: PtBackendDescriptor): boolean {
  return JSON.stringify(left) === JSON.stringify(right);
}

function mapComputationResult(
  source: WireComputationResult,
  expectedBackend: PtBackendDescriptor,
  expectedRequest: PtFlashRequest,
): PtComputationResult {
  const path = 'solve.result';
  const provenance = required(source.provenance, `${path}.provenance`);
  const backend = mapBackendDescriptor(
    required(provenance.backend, `${path}.provenance.backend`),
    `${path}.provenance.backend`,
  );
  if (!sameBackendSnapshot(backend, expectedBackend)) {
    fail(`${path}.provenance.backend`, 'snapshot differs from the selected discovery capability');
  }
  if (expectedRequest.configuredBackendId !== backend.configuredBackendId) {
    fail(`${path}.provenance.backend`, 'snapshot does not match the requested backend');
  }

  const outcome = mapOutcome(source.outcome, `${path}.outcome`);
  if (source.feed.length > MAX_COMPONENTS_PER_BACKEND) {
    fail(`${path}.feed`, 'component count exceeds client limits');
  }
  const feed = source.feed.map((component, index) => ({
    componentId: identifier(component.componentId, `${path}.feed[${index}].component_id`),
    moleFraction: fraction(component.moleFraction, `${path}.feed[${index}].mole_fraction`),
  }));
  if (feed.length !== backend.componentInventory.components.length) {
    fail(`${path}.feed`, 'feed count differs from the provenance inventory');
  }
  feed.forEach((component, index) => {
    if (component.componentId !== backend.componentInventory.components[index]?.componentId) {
      fail(`${path}.feed[${index}]`, 'feed order differs from the provenance inventory');
    }
  });
  normalized(
    feed.map((component) => component.moleFraction),
    `${path}.feed`,
  );
  if (expectedRequest.feed.length !== feed.length) {
    fail(`${path}.feed`, 'evaluated feed count differs from the submitted request');
  }
  for (const component of feed) {
    const matches = expectedRequest.feed.filter(
      (expected) => expected.componentId === component.componentId,
    );
    if (matches.length !== 1) {
      fail(`${path}.feed`, 'evaluated feed cannot be matched to the submitted request');
    }
    const requested = matches[0]?.moleFraction;
    if (requested === undefined || !Number.isFinite(requested)) {
      fail(`${path}.feed`, 'submitted request contains an invalid feed value');
    }
    const scale = Math.max(1, Math.abs(requested), Math.abs(component.moleFraction));
    if (
      Math.abs(requested - component.moleFraction) >
      2 * 64 * Number.EPSILON * scale
    ) {
      fail(`${path}.feed`, 'evaluated feed differs from the submitted request');
    }
  }

  if (source.phases.length > MAX_PHASES) {
    fail(`${path}.phases`, 'phase count exceeds client limits');
  }
  const phases = source.phases.map((phase, phaseIndex) => {
    const phasePath = `${path}.phases[${phaseIndex}]`;
    const wirePhaseIndex = uint32(phase.phaseIndex, `${phasePath}.phase_index`);
    if (wirePhaseIndex !== phaseIndex) {
      fail(`${phasePath}.phase_index`, 'phase indices must be contiguous vector indices');
    }
    if (phase.components.length !== feed.length) {
      fail(`${phasePath}.components`, 'component count differs from the feed');
    }
    const components = phase.components.map((component, componentIndex) => {
      const componentPath = `${phasePath}.components[${componentIndex}]`;
      const componentId = identifier(component.componentId, `${componentPath}.component_id`);
      if (componentId !== feed[componentIndex]?.componentId) {
        fail(componentPath, 'phase component order differs from the evaluated feed');
      }
      return {
        componentId,
        moleFraction: fraction(component.moleFraction, `${componentPath}.mole_fraction`),
        lnFugacityCoefficient: finite(
          component.lnFugacityCoefficient,
          `${componentPath}.ln_fugacity_coefficient`,
        ),
      };
    });
    normalized(
      components.map((component) => component.moleFraction),
      `${phasePath}.components`,
    );

    const metadataNamespace = backend.capability.phaseMetadataNamespace;
    if ((metadataNamespace.length === 0) === (phase.providerMetadata !== undefined)) {
      fail(`${phasePath}.provider_metadata`, 'presence disagrees with metadata namespace');
    }
    const providerMetadata =
      phase.providerMetadata === undefined
        ? undefined
        : {
            roleId: boundedString(
              phase.providerMetadata.roleId,
              `${phasePath}.provider_metadata.role_id`,
            ),
            familyId: boundedString(
              phase.providerMetadata.familyId,
              `${phasePath}.provider_metadata.family_id`,
            ),
          };
    const compressibilityFactor =
      phase.compressibilityFactor === undefined
        ? undefined
        : positive(phase.compressibilityFactor, `${phasePath}.compressibility_factor`);
    const base = {
      phaseIndex: wirePhaseIndex,
      molePhaseFraction: fraction(
        phase.molePhaseFraction,
        `${phasePath}.mole_phase_fraction`,
      ),
      components,
      providerBranch: required(
        phase.providerBranch,
        `${phasePath}.provider_branch`,
      ).toString(),
      providerBranchSmooth: required(
        phase.providerBranchSmooth,
        `${phasePath}.provider_branch_smooth`,
      ),
    };
    return {
      ...base,
      ...(compressibilityFactor === undefined ? {} : { compressibilityFactor }),
      ...(providerMetadata === undefined ? {} : { providerMetadata }),
    };
  });

  if (outcome === 'accepted') {
    if (
      phases.length === 0 ||
      !backend.capability.supportedPhaseCounts.includes(phases.length)
    ) {
      fail(`${path}.phases`, 'accepted result has an unsupported phase count');
    }
    normalized(
      phases.map((phase) => phase.molePhaseFraction),
      `${path}.phases`,
    );
  } else if (phases.length !== 0) {
    fail(`${path}.phases`, 'non-accepted result must not publish phases');
  }

  const globalStabilityProven = required(
    source.globalStabilityProven,
    `${path}.global_stability_proven`,
  );
  const morphologyResolved = required(
    source.morphologyResolved,
    `${path}.morphology_resolved`,
  );
  if (globalStabilityProven && !backend.capability.globalStabilityProven) {
    fail(`${path}.global_stability_proven`, 'result exceeds discovered capability');
  }
  if (morphologyResolved && backend.capability.phaseMetadataNamespace.length === 0) {
    fail(`${path}.morphology_resolved`, 'resolved morphology has no metadata namespace');
  }

  const pressurePa = positive(source.pressurePa, `${path}.pressure_pa`);
  const temperatureK = positive(source.temperatureK, `${path}.temperature_k`);
  if (
    pressurePa !== expectedRequest.pressurePa ||
    temperatureK !== expectedRequest.temperatureK
  ) {
    fail(path, 'evaluated PT state differs from the submitted request');
  }

  return {
    outcome,
    provenance: {
      backend,
      backendResultConvention: exactConvention(
        provenance.backendResultConvention,
        PT_BACKEND_RESULT_CONVENTION,
        `${path}.provenance.backend_result_convention`,
      ),
      phaseSetConvention: exactConvention(
        provenance.phaseSetConvention,
        PT_PHASE_SET_CONVENTION,
        `${path}.provenance.phase_set_convention`,
      ),
      phaseTransitionConvention: exactConvention(
        provenance.phaseTransitionConvention,
        PT_PHASE_TRANSITION_CONVENTION,
        `${path}.provenance.phase_transition_convention`,
      ),
      providerResultConvention: boundedString(
        provenance.providerResultConvention,
        `${path}.provenance.provider_result_convention`,
      ),
    },
    pressurePa,
    temperatureK,
    feed,
    phases,
    transitionReport: mapTransitionReport(
      required(source.transitionReport, `${path}.transition_report`),
      `${path}.transition_report`,
      backend,
    ),
    globalStabilityProven,
    morphologyResolved,
    diagnostic: boundedString(source.diagnostic, `${path}.diagnostic`, true),
  };
}

function mapServiceError(source: WireServiceError): PtServiceError {
  return {
    code: mapServiceErrorCode(source.code, 'solve.error.code'),
    field: boundedString(source.field, 'solve.error.field'),
    diagnostic: boundedString(source.diagnostic, 'solve.error.diagnostic'),
  };
}

export function mapSolveResponse(
  source: WireSolveResponse,
  expectedBackend: PtBackendDescriptor,
  expectedRequest: PtFlashRequest,
): PtFlashResponse {
  exactConvention(source.wireContract, PT_WIRE_CONTRACT, 'solve.wire_contract');
  exactConvention(
    source.serviceResultConvention,
    PT_SERVICE_RESULT_CONVENTION,
    'solve.service_result_convention',
  );
  switch (source.payload.case) {
    case 'result':
      return {
        kind: 'result',
        result: mapComputationResult(
          source.payload.value,
          expectedBackend,
          expectedRequest,
        ),
      };
    case 'error':
      return { kind: 'service_error', error: mapServiceError(source.payload.value) };
    case undefined:
      return fail('solve.payload', 'response contains neither result nor service error');
  }
}
