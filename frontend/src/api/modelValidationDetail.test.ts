import { Code } from '@connectrpc/connect';
import { expect, it } from 'vitest';
import { MODEL_VALIDATION_DETAIL_VERSION as version, readModelValidationDetail } from './modelValidationDetail';

it('preserves native key, pair and index selectors and returns an immutable copy', () => {
  for (const field of ['parameters.pure[nitrogen.r1].critical_temperature_k',
    'parameters.binary[nitrogen,ethane].kij.provenance.reference',
    'definition.components[1].molar_mass_kg_per_mol.value',
    'definition.pr76.pure[1].critical_temperature_k.provenance.kind',
    'definition.pr76.binary[0].kij.provenance', 'definition.applicability.provenance.kind',
    'definition.components[1].kind', 'definition.family', 'settings.kind',
    'components[0].molar_mass_kg_per_mol', 'eos_root.max_iterations', 'scalar.value', 'request']) {
    const raw = { version, code: 'configuration.missing_field', field };
    const detail = readModelValidationDetail(raw, Code.InvalidArgument)!;
    expect(detail).toEqual(raw); expect(Object.isFrozen(detail)).toBe(true);
    raw.field = 'changed'; expect(detail.field).toBe(field);
  }
});
it('omits unsafe fields and private values without losing the public domain code', () => {
  for (const field of ['parameters.pure[mh1_secret].critical_pressure_pa', 'components[ms1_private]',
    'parameters.pure[private-value].acentric_factor', 'authorization', 'model_handle',
    'components[secret-token]', 'components[a/b]', 'components[a]\n', 'components[a].invented',
    'components[a].', 'components[a]garbage', 'x'.repeat(257), { password: 'private-value' }]) {
    expect(readModelValidationDetail({ version, code: 'configuration.invalid_value', field }, Code.InvalidArgument, ['private-value']))
      .toEqual({ version, code: 'configuration.invalid_value' });
  }
});
it('requires known versions and exact validation status/domain pairs; never forwards arbitrary details', () => {
  for (const raw of [null, [], {}, { version: 'future/v2', code: 'configuration.invalid_value' },
    { version, code: 'configuration.invented' }, { version, code: 'rpc.internal_failure' },
    { version, code: 'session.not_found' }, { version, code: 'configuration.unsupported_preset' },
    { version, code: 'configuration.invalid_value', message: 'private exception' }]) {
    expect(readModelValidationDetail(raw, Code.InvalidArgument)).toBeUndefined();
  }
  expect(readModelValidationDetail({ version, code: 'configuration.unsupported_preset', field: 'preset_id' }, Code.Unimplemented)?.field).toBe('preset_id');
  expect(readModelValidationDetail({ version, code: 'configuration.resource_limit' }, Code.ResourceExhausted)?.code).toBe('configuration.resource_limit');
});

it('preserves solve scalar/aggregate fields and only canonical numeric feed indices', () => {
  for (const field of ['pressure_pa', 'temperature_k', 'feed', 'feed[0]', 'feed[1]', 'feed[127]']) {
    const raw = { version, code: 'request.rejected', field };
    const detail = readModelValidationDetail(raw, Code.InvalidArgument)!;
    expect(detail).toEqual(raw); expect(Object.isFrozen(detail)).toBe(true);
  }
  for (const field of ['feed[nitrogen]', 'feed[0,1]', 'feed[-1]', 'feed[01]', 'feed[1.0]', 'feed[1:2]', 'feed[]']) {
    expect(readModelValidationDetail({ version, code: 'request.rejected', field }, Code.InvalidArgument))
      .toEqual({ version, code: 'request.rejected' });
  }
});
