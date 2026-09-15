import { Code, createClient } from '@connectrpc/connect';
import { createGrpcWebTransport } from '@connectrpc/connect-web';

import {
  ModelConfigurationService,
  ModelSessionService,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';
import {
  ModelClientError,
  ModelSessionClient,
  type ModelSessionConnection,
  type ModelSessionConnector,
} from './modelSessionClient';

const TOKEN_MAX_LENGTH = 8192;
const BEARER_TOKEN = /^[A-Za-z0-9\-._~+/]+=*$/u;
const ENDPOINT_ANCHOR = 'https://mpmc.invalid';

export interface WebIdentityCredential {
  /** Opaque access token issued by the authoritative Web identity provider. */
  readonly accessToken: string;
  /** Optional absolute expiry in Unix epoch milliseconds. */
  readonly expiresAtEpochMs?: number;
}

export interface WebIdentityProvider {
  /**
   * Return the current authenticated browser credential. Implementations must not
   * persist the token in this module; reconnect explicitly calls this again.
   */
  getAccessToken(signal: AbortSignal): Promise<WebIdentityCredential | null>;
}

export interface AuthenticatedWebModelSessionOptions {
  /** Same-origin root-relative gRPC-Web base path, for example `/model-api`. */
  readonly baseUrl: string;
  readonly identity: WebIdentityProvider;
}

function failure(code: Code, reason: string): ModelClientError {
  return new ModelClientError(code, reason);
}

/**
 * Bearer-bearing model routes are deliberately same-origin in v1. Cross-origin
 * token forwarding belongs to the later reviewed edge-identity increment.
 */
export function normalizeAuthenticatedWebModelBaseUrl(baseUrl: string): string {
  const configured = baseUrl.trim();
  if (configured.length === 0 || !configured.startsWith('/') || configured.startsWith('//')) {
    throw failure(Code.InvalidArgument, 'web.endpoint_same_origin_required');
  }
  if (/\\|[\u0000-\u001f\u007f]/u.test(configured)) {
    throw failure(Code.InvalidArgument, 'web.endpoint_invalid');
  }

  let parsed: URL;
  try {
    parsed = new URL(configured, ENDPOINT_ANCHOR);
  } catch {
    throw failure(Code.InvalidArgument, 'web.endpoint_invalid');
  }
  if (
    parsed.origin !== ENDPOINT_ANCHOR ||
    parsed.username !== '' ||
    parsed.password !== '' ||
    parsed.search !== '' ||
    parsed.hash !== ''
  ) {
    throw failure(Code.InvalidArgument, 'web.endpoint_same_origin_required');
  }

  const normalized = parsed.pathname.replace(/\/+$/u, '');
  return normalized.length === 0 ? '/' : normalized;
}

function validateCredential(value: unknown): WebIdentityCredential {
  if (typeof value !== 'object' || value === null || !('accessToken' in value)) {
    throw failure(Code.Unauthenticated, 'web.identity_missing');
  }
  const accessToken = (value as { accessToken?: unknown }).accessToken;
  if (
    typeof accessToken !== 'string' ||
    accessToken.length === 0 ||
    accessToken.length > TOKEN_MAX_LENGTH ||
    !BEARER_TOKEN.test(accessToken)
  ) {
    throw failure(Code.Unauthenticated, 'web.identity_invalid');
  }

  const expiresAtEpochMs = (value as { expiresAtEpochMs?: unknown }).expiresAtEpochMs;
  if (expiresAtEpochMs !== undefined) {
    if (typeof expiresAtEpochMs !== 'number' || !Number.isFinite(expiresAtEpochMs)) {
      throw failure(Code.Unauthenticated, 'web.identity_invalid');
    }
    if (expiresAtEpochMs <= Date.now()) {
      throw failure(Code.Unauthenticated, 'web.identity_expired');
    }
    return Object.freeze({ accessToken, expiresAtEpochMs });
  }
  return Object.freeze({ accessToken });
}

async function acquireCredential(
  identity: WebIdentityProvider,
  signal: AbortSignal,
): Promise<WebIdentityCredential> {
  if (signal.aborted) throw failure(Code.Canceled, 'web.identity_cancelled');

  return new Promise<WebIdentityCredential>((resolve, reject) => {
    let settled = false;
    const finish = (work: () => void) => {
      if (settled) return;
      settled = true;
      signal.removeEventListener('abort', onAbort);
      work();
    };
    const onAbort = () => finish(() => reject(failure(Code.Canceled, 'web.identity_cancelled')));
    signal.addEventListener('abort', onAbort, { once: true });

    void Promise.resolve()
      .then(() => identity.getAccessToken(signal))
      .then(
        (credential) => {
          if (signal.aborted) {
            finish(() => reject(failure(Code.Canceled, 'web.identity_cancelled')));
            return;
          }
          try {
            const validated = validateCredential(credential);
            finish(() => resolve(validated));
          } catch (cause) {
            const error = cause instanceof ModelClientError
              ? cause
              : failure(Code.Unauthenticated, 'web.identity_invalid');
            finish(() => reject(error));
          }
        },
        () => finish(() => reject(
          signal.aborted
            ? failure(Code.Canceled, 'web.identity_cancelled')
            : failure(Code.Unauthenticated, 'web.identity_unavailable'),
        )),
      );
  });
}

/**
 * Build a connector for ModelSessionClient using one immutable identity snapshot
 * per session epoch. The same Authorization header is used for the lease and all
 * model RPCs in that epoch. Explicit reconnect obtains a fresh credential.
 */
export function createAuthenticatedWebModelSessionConnector(
  options: AuthenticatedWebModelSessionOptions,
): ModelSessionConnector {
  const endpoint = normalizeAuthenticatedWebModelBaseUrl(options.baseUrl);
  if (!options.identity || typeof options.identity.getAccessToken !== 'function') {
    throw failure(Code.InvalidArgument, 'web.identity_provider_required');
  }

  return async (signal: AbortSignal): Promise<ModelSessionConnection> => {
    const credential = await acquireCredential(options.identity, signal);
    if (signal.aborted) throw failure(Code.Canceled, 'web.identity_cancelled');

    const transport = createGrpcWebTransport({
      baseUrl: endpoint,
      useBinaryFormat: true,
    });
    const headers = new Headers();
    headers.set('authorization', `Bearer ${credential.accessToken}`);

    return {
      models: createClient(ModelConfigurationService, transport),
      sessions: createClient(ModelSessionService, transport),
      headers,
      // connect-web owns requests rather than a persistent browser socket. The
      // ModelSessionClient AbortSignal cancels the stream and admitted RPCs.
      close() {},
    };
  };
}

export function createAuthenticatedWebModelSessionClient(
  options: AuthenticatedWebModelSessionOptions,
): ModelSessionClient {
  return new ModelSessionClient(createAuthenticatedWebModelSessionConnector(options));
}
