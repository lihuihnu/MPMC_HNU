import { Code, ConnectError } from '@connectrpc/connect';
import { MODEL_WIRE_CONTRACT } from '../../src/api/modelSessionClient';
import { ModelServiceErrorSchema } from '../../src/gen/mpmc/model_configuration/v1/model_service_pb';

/** Test observer only: an admitted Describe can retain an Entry while its registry closes. */
export function modelCleanupComplete(cause: unknown, expected: 'session.not_found' | 'registry.model_not_found'): boolean {
  const error = ConnectError.from(cause);
  const details = error.findDetails(ModelServiceErrorSchema);
  const detail = details.length === 1 ? details[0] : undefined;
  if (detail?.wireContract === MODEL_WIRE_CONTRACT) {
    if (error.code === Code.NotFound && detail.code === expected) return true;
    if (expected === 'session.not_found' && error.code === Code.FailedPrecondition && detail.code === 'registry.closed') return false;
  }
  // Never print arbitrary remote text, credentials or routing identifiers.
  throw new Error(`Unexpected cleanup response (gRPC ${error.code}, expected ${expected}).`);
}
