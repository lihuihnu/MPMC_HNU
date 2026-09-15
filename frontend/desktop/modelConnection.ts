import { createClient } from '@connectrpc/connect';
import { createGrpcTransport, Http2SessionManager } from '@connectrpc/connect-node';
import { ModelConfigurationService, ModelSessionService } from '../src/gen/mpmc/model_configuration/v1/model_service_pb';
import type { ModelSessionConnection } from '../src/api/modelSessionClient';
import type { PtHostConnection } from './hostSession';

/** Main-process only: the host's private bearer never enters the renderer. */
export function createDesktopModelConnection(host: PtHostConnection): ModelSessionConnection {
  const manager = new Http2SessionManager(host.baseUrl);
  const transport = createGrpcTransport({ baseUrl: host.baseUrl, sessionManager: manager,
    readMaxBytes: 4 * 1024 * 1024, writeMaxBytes: 64 * 1024 });
  let closed = false;
  return {
    models: createClient(ModelConfigurationService, transport),
    sessions: createClient(ModelSessionService, transport),
    headers: { authorization: `Bearer ${host.bearerToken}` },
    close() { if (!closed) { closed = true; manager.abort(); } },
  };
}
