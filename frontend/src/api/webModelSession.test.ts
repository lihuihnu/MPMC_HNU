import { Code } from '@connectrpc/connect';
import { describe, expect, it, vi } from 'vitest';

import { ModelClientError, ModelSessionClient } from './modelSessionClient';
import {
  createAuthenticatedWebModelSessionClient,
  createAuthenticatedWebModelSessionConnector,
  normalizeAuthenticatedWebModelBaseUrl,
  type WebIdentityProvider,
} from './webModelSession';

function provider(accessToken: () => string): WebIdentityProvider {
  return {
    getAccessToken: vi.fn(async (_signal: AbortSignal) => ({ accessToken: accessToken() })),
  };
}

describe('authenticated Web model session identity binding', () => {
  it('accepts only same-origin root-relative model endpoints', () => {
    expect(normalizeAuthenticatedWebModelBaseUrl('/model-api/')).toBe('/model-api');
    expect(normalizeAuthenticatedWebModelBaseUrl('/')).toBe('/');

    for (const endpoint of [
      '',
      'model-api',
      'https://example.test/model-api',
      '//example.test/model-api',
      '/model-api?token=bad',
      '/model-api#fragment',
      '/model\\api',
    ]) {
      expect(() => normalizeAuthenticatedWebModelBaseUrl(endpoint)).toThrowError(
        expect.objectContaining({ code: Code.InvalidArgument }),
      );
    }
  });

  it('takes one opaque identity snapshot per connector invocation', async () => {
    let token = 'token-a.ABC_123';
    const identity = provider(() => token);
    const connector = createAuthenticatedWebModelSessionConnector({
      baseUrl: '/model-api',
      identity,
    });

    const first = await connector(new AbortController().signal);
    expect(new Headers(first.headers).get('authorization')).toBe('Bearer token-a.ABC_123');
    expect(identity.getAccessToken).toHaveBeenCalledTimes(1);

    token = 'token-b.DEF_456';
    expect(new Headers(first.headers).get('authorization')).toBe('Bearer token-a.ABC_123');
    const second = await connector(new AbortController().signal);
    expect(new Headers(second.headers).get('authorization')).toBe('Bearer token-b.DEF_456');
    expect(identity.getAccessToken).toHaveBeenCalledTimes(2);

    first.close();
    second.close();
  });

  it('settles cancellation even when the identity provider ignores its signal', async () => {
    const identity: WebIdentityProvider = {
      getAccessToken: vi.fn(async () => new Promise<never>(() => {})),
    };
    const connector = createAuthenticatedWebModelSessionConnector({
      baseUrl: '/model-api',
      identity,
    });
    const controller = new AbortController();
    const opening = connector(controller.signal);
    controller.abort();
    await expect(opening).rejects.toMatchObject({
      code: Code.Canceled,
      reason: 'web.identity_cancelled',
    });
  });

  it('fails closed for missing, malformed, expired and provider-failed credentials', async () => {
    const cases: Array<[WebIdentityProvider, string]> = [
      [{ getAccessToken: async () => null }, 'web.identity_missing'],
      [{ getAccessToken: async () => ({ accessToken: 'contains whitespace' }) }, 'web.identity_invalid'],
      [{ getAccessToken: async () => ({ accessToken: 'valid-token', expiresAtEpochMs: Date.now() - 1 }) }, 'web.identity_expired'],
    ];

    for (const [identity, reason] of cases) {
      const connector = createAuthenticatedWebModelSessionConnector({ baseUrl: '/model-api', identity });
      await expect(connector(new AbortController().signal)).rejects.toMatchObject({
        code: Code.Unauthenticated,
        reason,
      });
    }

    const secret = 'provider-private-secret';
    const failing: WebIdentityProvider = {
      getAccessToken: async () => { throw new Error(secret); },
    };
    const connector = createAuthenticatedWebModelSessionConnector({ baseUrl: '/model-api', identity: failing });
    let observed: unknown;
    try {
      await connector(new AbortController().signal);
    } catch (cause) {
      observed = cause;
    }
    expect(observed).toBeInstanceOf(ModelClientError);
    expect(observed).toMatchObject({ code: Code.Unauthenticated, reason: 'web.identity_unavailable' });
    expect(String(observed)).not.toContain(secret);
  });

  it('composes the authenticated connector with the shared ModelSessionClient', async () => {
    const client = createAuthenticatedWebModelSessionClient({
      baseUrl: '/model-api',
      identity: provider(() => 'opaque-web-token'),
    });
    expect(client).toBeInstanceOf(ModelSessionClient);
    expect(client.connected).toBe(false);
    await client.dispose();
  });
});
