export const PT_WIRE_CONTRACT = 'mpmc.runtime.v1/PT-flash-service/v1';
export const PT_SERVICE_BOUNDARY_CONVENTION = 'PT/service-boundary/v1';
export const PT_COMPONENT_INVENTORY_CONVENTION =
  'PT/service-component-inventory/v1';
export const PT_CAPABILITY_CONVENTION = 'PT/service-capability-discovery/v1';
export const PT_SERVICE_REQUEST_CONVENTION = 'PT/service-request/v1';
export const PT_SERVICE_RESULT_CONVENTION = 'PT/service-result/v1';
export const PT_BACKEND_CAPABILITY_CONVENTION =
  'PT/flash-backend/capability-and-dispatch/v1';
export const PT_BACKEND_RESULT_CONVENTION = 'PT/flash-backend/result/v1';
export const PT_PHASE_SET_CONVENTION = 'PT/phase-set-v1';
export const PT_PHASE_TRANSITION_CONVENTION =
  'PT/phase-transition-boundary/v1';

export type PtPhaseTransitionSupport =
  | 'detection_only'
  | 'fresh_target_resolve';

export interface PtPhaseTransitionEdgeCapability {
  sourcePhaseCount: number;
  targetPhaseCount: number;
  support: PtPhaseTransitionSupport;
  requiresFreshTargetSolve: boolean;
}

export interface PtPhaseTransitionCapability {
  convention: typeof PT_PHASE_TRANSITION_CONVENTION;
  edges: readonly PtPhaseTransitionEdgeCapability[];
}

export interface PtBackendScalarSetting {
  id: string;
  value: number;
  unit: string;
}

export interface PtBackendCapability {
  convention: typeof PT_BACKEND_CAPABILITY_CONVENTION;
  backendId: string;
  modelProfile: string;
  algorithmProfile: string;
  publicationProfile: string;
  configurationProfile: string;
  datasetId: string;
  revision: string;
  componentIds: readonly string[];
  supportedPhaseCounts: readonly number[];
  scalarSettings: readonly PtBackendScalarSetting[];
  transitionCapability: PtPhaseTransitionCapability;
  performsInitialStabilitySearch: boolean;
  performsFinalPhaseSetReview: boolean;
  performsBoundaryNeighborResolve: boolean;
  globalStabilityProven: boolean;
  phaseMetadataNamespace: string;
}

export interface PtRuntimeComponent {
  componentId: string;
  feedIndex: number;
}

export interface PtRuntimeComponentInventory {
  convention: typeof PT_COMPONENT_INVENTORY_CONVENTION;
  components: readonly PtRuntimeComponent[];
}

export interface PtBackendDescriptor {
  convention: typeof PT_CAPABILITY_CONVENTION;
  configuredBackendId: string;
  capability: PtBackendCapability;
  componentInventory: PtRuntimeComponentInventory;
}

export interface PtCapabilityDiscovery {
  wireContract: typeof PT_WIRE_CONTRACT;
  serviceBoundaryConvention: typeof PT_SERVICE_BOUNDARY_CONVENTION;
  capabilityConvention: typeof PT_CAPABILITY_CONVENTION;
  backends: readonly PtBackendDescriptor[];
}

export interface PtCompositionEntry {
  componentId: string;
  moleFraction: number;
}

export interface PtFlashRequest {
  configuredBackendId: string;
  pressurePa: number;
  temperatureK: number;
  feed: readonly PtCompositionEntry[];
}

export interface PtPhaseComponent {
  componentId: string;
  moleFraction: number;
  lnFugacityCoefficient: number;
}

export interface PtPhaseMetadata {
  roleId: string;
  familyId: string;
}

export interface PtPhase {
  phaseIndex: number;
  molePhaseFraction: number;
  components: readonly PtPhaseComponent[];
  providerBranch: string;
  providerBranchSmooth: boolean;
  compressibilityFactor?: number;
  providerMetadata?: PtPhaseMetadata;
}

export type PtPhaseTransitionTrigger =
  | 'initial_stability_witness'
  | 'final_phase_set_instability'
  | 'phase_disappearance'
  | 'provider_topology_witness'
  | 'provider_boundary_route';

export type PtPhaseTransitionResolution =
  | 'accepted_target'
  | 'target_resolve_required'
  | 'target_resolve_failed'
  | 'candidate_not_accepted'
  | 'broader_topology_required'
  | 'indeterminate';

export interface PtPhaseTransitionEvidence {
  sourcePhaseCount: number;
  targetPhaseCount?: number;
  trigger: PtPhaseTransitionTrigger;
  resolution: PtPhaseTransitionResolution;
  freshTargetSolveAttempted: boolean;
  targetTopologyClosed: boolean;
  providerEvidenceProfile: string;
  diagnostic: string;
}

export interface PtPhaseTransitionReport {
  convention: typeof PT_PHASE_TRANSITION_CONVENTION;
  evidence: readonly PtPhaseTransitionEvidence[];
}

export interface PtResultProvenance {
  backend: PtBackendDescriptor;
  backendResultConvention: typeof PT_BACKEND_RESULT_CONVENTION;
  phaseSetConvention: typeof PT_PHASE_SET_CONVENTION;
  phaseTransitionConvention: typeof PT_PHASE_TRANSITION_CONVENTION;
  providerResultConvention: string;
}

export type PtComputationOutcome =
  | 'accepted'
  | 'phase_set_unstable'
  | 'indeterminate';

export interface PtComputationResult {
  outcome: PtComputationOutcome;
  provenance: PtResultProvenance;
  pressurePa: number;
  temperatureK: number;
  feed: readonly PtCompositionEntry[];
  phases: readonly PtPhase[];
  transitionReport: PtPhaseTransitionReport;
  globalStabilityProven: boolean;
  morphologyResolved: boolean;
  diagnostic: string;
}

export type PtServiceErrorCode =
  | 'invalid_configured_backend_id'
  | 'configured_backend_not_found'
  | 'invalid_pressure'
  | 'invalid_temperature'
  | 'component_count_mismatch'
  | 'invalid_component_id'
  | 'duplicate_component'
  | 'unknown_component'
  | 'invalid_mole_fraction'
  | 'composition_not_normalized'
  | 'backend_rejected_request'
  | 'backend_contract_violation'
  | 'backend_execution_failure';

export interface PtServiceError {
  code: PtServiceErrorCode;
  field: string;
  diagnostic: string;
}

export type PtFlashResponse =
  | { kind: 'result'; result: PtComputationResult }
  | { kind: 'service_error'; error: PtServiceError };

export type PtResultSeverity = 'accepted' | 'warning' | 'indeterminate';

export const PT_COMPOSITION_ROUNDOFF = 64 * Number.EPSILON;

export function validatePtFlashRequest(
  request: PtFlashRequest,
  backend: PtBackendDescriptor,
): string[] {
  const errors: string[] = [];

  if (request.configuredBackendId !== backend.configuredBackendId) {
    errors.push('Selected backend does not match the discovered capability snapshot.');
  }
  if (!Number.isFinite(request.pressurePa) || request.pressurePa <= 0) {
    errors.push('Pressure must be finite and greater than 0 Pa.');
  }
  if (!Number.isFinite(request.temperatureK) || request.temperatureK <= 0) {
    errors.push('Temperature must be finite and greater than 0 K.');
  }

  const expectedIds = new Set(
    backend.componentInventory.components.map((component) => component.componentId),
  );
  if (request.feed.length !== expectedIds.size) {
    errors.push('Feed must contain every discovered component exactly once.');
  }

  const seen = new Set<string>();
  let sum = 0;
  let correction = 0;
  for (const component of request.feed) {
    if (!expectedIds.has(component.componentId)) {
      errors.push(`Unknown component ID: ${component.componentId || '(blank)'}.`);
    } else if (seen.has(component.componentId)) {
      errors.push(`Duplicate component ID: ${component.componentId}.`);
    } else {
      seen.add(component.componentId);
    }

    if (
      !Number.isFinite(component.moleFraction) ||
      component.moleFraction < 0 ||
      component.moleFraction > 1
    ) {
      errors.push(`Component ${component.componentId} must have z in [0, 1].`);
      continue;
    }

    const increment = component.moleFraction - correction;
    const next = sum + increment;
    correction = (next - sum) - increment;
    sum = next;
  }

  if (seen.size !== expectedIds.size) {
    errors.push('Feed component IDs do not match the discovered inventory.');
  }
  if (Math.abs(sum - 1) > PT_COMPOSITION_ROUNDOFF) {
    errors.push(
      `Overall mole fractions must sum to 1 within roundoff; current sum is ${sum.toPrecision(16)}.`,
    );
  }
  return errors;
}

export function ptOutcomeLabel(outcome: PtComputationOutcome): string {
  switch (outcome) {
    case 'accepted':
      return 'Accepted phase set';
    case 'phase_set_unstable':
      return 'Phase set unstable';
    case 'indeterminate':
      return 'Scientific result indeterminate';
  }
}

export function ptOutcomeSeverity(outcome: PtComputationOutcome): PtResultSeverity {
  switch (outcome) {
    case 'accepted':
      return 'accepted';
    case 'phase_set_unstable':
      return 'warning';
    case 'indeterminate':
      return 'indeterminate';
  }
}

export function ptServiceErrorLabel(code: PtServiceErrorCode): string {
  switch (code) {
    case 'invalid_configured_backend_id':
    case 'configured_backend_not_found':
      return 'Backend selection error';
    case 'invalid_pressure':
    case 'invalid_temperature':
    case 'component_count_mismatch':
    case 'invalid_component_id':
    case 'duplicate_component':
    case 'unknown_component':
    case 'invalid_mole_fraction':
    case 'composition_not_normalized':
      return 'Request rejected by service';
    case 'backend_rejected_request':
      return 'Backend rejected request';
    case 'backend_contract_violation':
      return 'Backend contract violation';
    case 'backend_execution_failure':
      return 'Backend execution failure';
  }
}
