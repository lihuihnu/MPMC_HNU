import { describe, expect, it } from 'vitest';

import {
  physicalRoleLabel,
  topologyStatusLabel,
  topologyStatusSeverity,
  validateProfileCPtRequest,
  type ProfileCPtRequest,
} from './flash';

function request(overrides: Partial<ProfileCPtRequest> = {}): ProfileCPtRequest {
  return {
    pressurePa: 10e6,
    temperatureK: 350,
    naclMolalityMolPerKgWater: 0,
    components: [
      { id: 'carbon-dioxide', moleFraction: 0.7 },
      { id: 'water', moleFraction: 0.3 },
    ],
    ...overrides,
  };
}

describe('Profile-C PT request validation', () => {
  it('accepts a finite ordered composition that already sums to one', () => {
    expect(validateProfileCPtRequest(request())).toEqual([]);
  });

  it('does not normalize or silently repair an invalid composition sum', () => {
    const errors = validateProfileCPtRequest(
      request({
        components: [
          { id: 'carbon-dioxide', moleFraction: 0.69 },
          { id: 'water', moleFraction: 0.3 },
        ],
      }),
    );

    expect(errors.some((error) => error.includes('must sum to 1'))).toBe(true);
  });

  it('rejects duplicate stable component IDs', () => {
    const errors = validateProfileCPtRequest(
      request({
        components: [
          { id: 'water', moleFraction: 0.5 },
          { id: 'water', moleFraction: 0.5 },
        ],
      }),
    );

    expect(errors).toContain('Duplicate component ID: water.');
  });

  it('rejects nonphysical pressure, temperature and salinity inputs', () => {
    const errors = validateProfileCPtRequest(
      request({
        pressurePa: 0,
        temperatureK: Number.NaN,
        naclMolalityMolPerKgWater: -1,
      }),
    );

    expect(errors).toEqual(
      expect.arrayContaining([
        'Pressure must be finite and greater than 0 Pa.',
        'Temperature must be finite and greater than 0 K.',
        'NaCl molality must be finite and non-negative.',
      ]),
    );
  });
});

describe('Profile-C PT presentation semantics', () => {
  it('keeps nonaqueous morphology explicitly unresolved', () => {
    expect(physicalRoleLabel('nonaqueous_unclassified')).toContain(
      'morphology unresolved',
    );
  });

  it('distinguishes locally closed candidates from unresolved states', () => {
    expect(topologyStatusSeverity('w_h0_h1_locally_closed')).toBe('candidate');
    expect(topologyStatusSeverity('topology_unresolved')).toBe('warning');
    expect(topologyStatusSeverity('numerical_indeterminate')).toBe('indeterminate');
    expect(topologyStatusLabel('w_h0_h1_locally_closed')).toContain('3 phases');
  });
});
