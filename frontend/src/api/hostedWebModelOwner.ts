import { Code } from '@connectrpc/connect';

import {
  bindExpertModelOwner,
  type ExpertModelOwner,
  type ExpertOwnedModel,
} from './expertModelOwner';
import {
  ModelClientError,
  type ModelCallOptions,
  type ModelCreateInput,
  type ModelSessionClient,
} from './modelSessionClient';
import {
  createAuthenticatedWebModelSessionClient,
  type AuthenticatedWebModelSessionOptions,
  type WebIdentityProvider,
} from './webModelSession';

export const HOSTED_WEB_MODEL_OWNERSHIP_CONVENTION =
  'MPMC/model/hosted-web-ownership/v1' as const;

/**
 * Runtime-only host capability. It carries a token provider function, never a
 * bearer token value, and keeps the model endpoint same-origin by delegating to
 * webModelSession.ts validation.
 */
export interface HostedWebModelOwnershipCapability {
  readonly convention: typeof HOSTED_WEB_MODEL_OWNERSHIP_CONVENTION;
  readonly baseUrl: string;
  readonly identity: WebIdentityProvider;
}

export interface HostedWebExpertModelOwner extends ExpertModelOwner {
  /** Close the model-session lease and invalidate all session-local references. */
  dispose(): Promise<void>;
}

function validCapability(value: unknown): value is HostedWebModelOwnershipCapability {
  if (typeof value !== 'object' || value === null || Array.isArray(value)) return false;
  const capability = value as Record<string, unknown>;
  const keys = Object.keys(capability);
  if (
    keys.length !== 3 ||
    !keys.includes('convention') ||
    !keys.includes('baseUrl') ||
    !keys.includes('identity') ||
    capability.convention !== HOSTED_WEB_MODEL_OWNERSHIP_CONVENTION ||
    typeof capability.baseUrl !== 'string' ||
    typeof capability.identity !== 'object' ||
    capability.identity === null ||
    typeof (capability.identity as { getAccessToken?: unknown }).getAccessToken !== 'function'
  ) {
    return false;
  }
  return true;
}

class HostedWebExpertModelOwnerImpl implements HostedWebExpertModelOwner {
  readonly #delegate: ExpertModelOwner;
  #opening: Promise<void> | undefined;
  #disposed = false;

  constructor(private readonly session: ModelSessionClient) {
    this.#delegate = bindExpertModelOwner(session);
  }

  #startOpening(): Promise<void> {
    if (this.#opening !== undefined) return this.#opening;
    const opening = this.session.connect().finally(() => {
      if (this.#opening === opening) this.#opening = undefined;
    });
    this.#opening = opening;
    return opening;
  }

  async #ensureConnected(options: ModelCallOptions): Promise<void> {
    if (this.#disposed) {
      throw new ModelClientError(Code.FailedPrecondition, 'client.disposed');
    }
    if (this.session.connected) return;
    if (options.signal?.aborted) {
      throw new ModelClientError(Code.Canceled, 'web.owner_cancelled');
    }

    const opening = this.#startOpening();
    const signal = options.signal;
    if (signal === undefined) {
      await opening;
      return;
    }

    await new Promise<void>((resolve, reject) => {
      let settled = false;
      const finish = (work: () => void) => {
        if (settled) return;
        settled = true;
        signal.removeEventListener('abort', onAbort);
        work();
      };
      const onAbort = () => {
        void this.session.disconnect().catch(() => {});
        finish(() => reject(new ModelClientError(Code.Canceled, 'web.owner_cancelled')));
      };
      signal.addEventListener('abort', onAbort, { once: true });
      void opening.then(
        () => finish(resolve),
        (cause) => finish(() => reject(cause)),
      );
    });
  }

  async create(
    input: ModelCreateInput,
    options: ModelCallOptions = {},
  ): Promise<ExpertOwnedModel> {
    await this.#ensureConnected(options);
    if (this.#disposed) {
      throw new ModelClientError(Code.FailedPrecondition, 'client.disposed');
    }
    return this.#delegate.create(input, options);
  }

  async dispose(): Promise<void> {
    if (this.#disposed) return;
    this.#disposed = true;
    await this.session.dispose();
  }
}

export type HostedWebModelSessionFactory = (
  options: AuthenticatedWebModelSessionOptions,
) => ModelSessionClient;

/**
 * Return a product owner only for an explicit, versioned runtime capability.
 * Malformed or absent host configuration fails closed to the existing hosted
 * configured-PT product. No identity acquisition or network call occurs here;
 * the first model create opens the authenticated session lazily.
 */
export function hostedWebExpertModelOwner(
  capability: unknown = typeof window === 'undefined'
    ? undefined
    : window.mpmcHostedWebModelOwnership,
  createSession: HostedWebModelSessionFactory = createAuthenticatedWebModelSessionClient,
): HostedWebExpertModelOwner | null {
  if (!validCapability(capability)) return null;
  try {
    const session = createSession({
      baseUrl: capability.baseUrl,
      identity: capability.identity,
    });
    return new HostedWebExpertModelOwnerImpl(session);
  } catch {
    return null;
  }
}

declare global {
  interface Window {
    /** Supplied by the trusted hosted application shell, never by build-time secrets. */
    mpmcHostedWebModelOwnership?: HostedWebModelOwnershipCapability;
  }
}
