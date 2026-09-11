export const PROFILE_C_MODEL = 'SW92/corrected-original/PR76-base/NaCl-molality';
export const PROFILE_C_PT_ALGORITHM =
  'SW92-equilibrium/phase-assigned-aq-na-joint/topology-orchestration-pt/v1';

export const topologyStatuses = [
  'no_w_single_h_locally_closed',
  'no_w_two_h_locally_closed',
  'w_h_locally_closed',
  'w_h0_h1_locally_closed',
  'higher_phase_count_or_wrong_candidate',
  'topology_unresolved',
  'numerical_indeterminate',
] as const;

export type TopologyStatus = (typeof topologyStatuses)[number];
export type PhysicalRole = 'aqueous' | 'nonaqueous_unclassified';
export type ThermodynamicFamily = 'aqueous' | 'nonaqueous';

export interface FlashComponentInput {
  id: string;
  moleFraction: number;
}

export interface ProfileCPtRequest {
  pressurePa: number;
  temperatureK: number;
  naclMolalityMolPerKgWater: number;
  components: readonly FlashComponentInput[];
}

export interface ProfileCPtPhase {
  physicalRole: PhysicalRole;
  thermodynamicFamily: ThermodynamicFamily;
  molePhaseFraction: number;
  composition: readonly number[];
  compressibilityFactor?: number;
}

export interface ProfileCPtResponse {
  status: TopologyStatus;
  pressurePa: number;
  temperatureK: number;
  naclMolalityMolPerKgWater: number;
  componentIds: readonly string[];
  phases: readonly ProfileCPtPhase[];
  diagnostic: string;
  modelProfile: typeof PROFILE_C_MODEL;
  algorithmProfile: typeof PROFILE_C_PT_ALGORITHM;
  globalStabilityProven: false;
  acceptedPhaseSetPublished: false;
  morphologyResolved: false;
}

export type ResultSeverity = 'candidate' | 'warning' | 'indeterminate';

const compositionRoundoff = 64 * Number.EPSILON;

export function validateProfileCPtRequest(request: ProfileCPtRequest): string[] {
  const errors: string[] = [];

  if (!Number.isFinite(request.pressurePa) || request.pressurePa <= 0) {
    errors.push('Pressure must be finite and greater than 0 Pa.');
  }
  if (!Number.isFinite(request.temperatureK) || request.temperatureK <= 0) {
    errors.push('Temperature must be finite and greater than 0 K.');
  }
  if (
    !Number.isFinite(request.naclMolalityMolPerKgWater) ||
    request.naclMolalityMolPerKgWater < 0
  ) {
    errors.push('NaCl molality must be finite and non-negative.');
  }
  if (request.components.length === 0) {
    errors.push('At least one ordered component is required.');
    return errors;
  }

  const ids = new Set<string>();
  let sum = 0;
  let correction = 0;
  let hasPositive = false;

  for (const component of request.components) {
    const id = component.id.trim();
    if (id.length === 0) {
      errors.push('Every component requires a stable component ID.');
    } else if (ids.has(id)) {
      errors.push(`Duplicate component ID: ${id}.`);
    } else {
      ids.add(id);
    }

    if (
      !Number.isFinite(component.moleFraction) ||
      component.moleFraction < 0 ||
      component.moleFraction > 1
    ) {
      errors.push(`Component ${id || '(blank)'} must have z in [0, 1].`);
      continue;
    }
    hasPositive ||= component.moleFraction > 0;

    // Kahan summation mirrors the backend's intent: do not hide invalid input
    // behind browser-side normalization or clipping.
    const increment = component.moleFraction - correction;
    const next = sum + increment;
    correction = (next - sum) - increment;
    sum = next;
  }

  if (!hasPositive) {
    errors.push('At least one component mole fraction must be positive.');
  }
  if (Math.abs(sum - 1) > compositionRoundoff) {
    errors.push(
      `Overall mole fractions must sum to 1 within roundoff; current sum is ${sum.toPrecision(16)}.`,
    );
  }

  return errors;
}

export function topologyStatusLabel(status: TopologyStatus): string {
  switch (status) {
    case 'no_w_single_h_locally_closed':
      return '1 phase · no aqueous phase';
    case 'no_w_two_h_locally_closed':
      return '2 phases · no aqueous phase';
    case 'w_h_locally_closed':
      return '2 phases · aqueous + nonaqueous';
    case 'w_h0_h1_locally_closed':
      return '3 phases · aqueous + 2 nonaqueous';
    case 'higher_phase_count_or_wrong_candidate':
      return 'Higher phase count / rival candidate detected';
    case 'topology_unresolved':
      return 'Topology unresolved';
    case 'numerical_indeterminate':
      return 'Numerically indeterminate';
  }
}

export function topologyStatusSeverity(status: TopologyStatus): ResultSeverity {
  switch (status) {
    case 'no_w_single_h_locally_closed':
    case 'no_w_two_h_locally_closed':
    case 'w_h_locally_closed':
    case 'w_h0_h1_locally_closed':
      return 'candidate';
    case 'higher_phase_count_or_wrong_candidate':
    case 'topology_unresolved':
      return 'warning';
    case 'numerical_indeterminate':
      return 'indeterminate';
  }
}

export function physicalRoleLabel(role: PhysicalRole): string {
  return role === 'aqueous'
    ? 'Aqueous'
    : 'Nonaqueous · morphology unresolved';
}
