import { create, toBinary } from '@bufbuild/protobuf';
import { Code, ConnectError } from '@connectrpc/connect';
import { expect, it } from 'vitest';
import { MODEL_WIRE_CONTRACT } from '../../src/api/modelSessionClient';
import { ModelServiceErrorSchema } from '../../src/gen/mpmc/model_configuration/v1/model_service_pb';
import { modelCleanupComplete } from './modelCleanupProbe';

function response(status: Code, code: string, wireContract = MODEL_WIRE_CONTRACT) {
  return new ConnectError('private remote text', status, undefined, [{
    type: ModelServiceErrorSchema.typeName,
    value: toBinary(ModelServiceErrorSchema, create(ModelServiceErrorSchema, { wireContract, code })),
  }]);
}
it('requires terminal session absence after an admitted registry closes', () => {
  const sequence = [response(Code.FailedPrecondition, 'registry.closed'), response(Code.NotFound, 'session.not_found')];
  expect(sequence.map(error => modelCleanupComplete(error, 'session.not_found'))).toEqual([false, true]);
  // No number of registry.closed responses proves that the session has disappeared.
  expect(modelCleanupComplete(sequence[0], 'session.not_found')).toBe(false);
});
it('requires exact released-model absence and never accepts registry closure for release', () => {
  expect(modelCleanupComplete(response(Code.NotFound, 'registry.model_not_found'), 'registry.model_not_found')).toBe(true);
  expect(() => modelCleanupComplete(response(Code.FailedPrecondition, 'registry.closed'), 'registry.model_not_found')).toThrow('Unexpected cleanup');
});
it('rejects wrong domains, wire versions, status codes and absent details', () => {
  for (const error of [response(Code.FailedPrecondition, 'session.host_closed'),
    response(Code.NotFound, 'registry.closed'), response(Code.NotFound, 'registry.model_not_found'),
    response(Code.NotFound, 'session.not_found', 'unknown/v2'),
    response(Code.Unavailable, 'session.not_found'), new ConnectError('private remote text', Code.NotFound)]) {
    expect(() => modelCleanupComplete(error, 'session.not_found')).toThrow('Unexpected cleanup');
  }
});
