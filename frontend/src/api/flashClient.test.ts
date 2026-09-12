import { Code, ConnectError } from '@connectrpc/connect';
import { describe, expect, it } from 'vitest';

import { type PtFlashRequest } from '../domain/flash';
import { wireDiscovery, wireResult } from '../test/ptWireFixtures';
import {
  GrpcWebFlashClient,
  PtTransportError,
  createGrpcWebFlashClient,
  type PtFlashRpc,
} from './flashClient';

describe('gRPC-Web PT client', () => {
  it('supports a same-origin root endpoint without treating it as empty', () => {
    expect(createGrpcWebFlashClient('/')).toMatchObject({
      configured: true,
      endpoint: '/',
    });
  });

  it('discovers before solving and makes exactly one solve call', async () => {
    let solveCalls = 0;
    const controller = new AbortController();
    const rpc: PtFlashRpc = {
      async discoverPtCapabilities(_request, options) {
        expect(options?.signal).toBe(controller.signal);
        expect(options?.timeoutMs).toBe(25);
        return wireDiscovery();
      },
      async solvePtFlash(request, options) {
        solveCalls += 1;
        expect(request.configuredBackendId).toBe('pr76-default');
        expect(request.feed).toHaveLength(2);
        expect(options?.signal).toBe(controller.signal);
        expect(options?.timeoutMs).toBe(50);
        return wireResult();
      },
    };
    const client = new GrpcWebFlashClient(rpc, 'https://example.test');
    const discovery = await client.discoverPtCapabilities({
      signal: controller.signal,
      timeoutMs: 25,
    });
    const backend = discovery.backends[0];
    expect(backend).toBeDefined();
    const request: PtFlashRequest = {
      configuredBackendId: 'pr76-default',
      pressurePa: 10e6,
      temperatureK: 350,
      feed: [
        { componentId: 'methane', moleFraction: 0.7 },
        { componentId: 'carbon-dioxide', moleFraction: 0.3 },
      ],
    };
    const response = await client.solvePtFlash(request, backend!, {
      signal: controller.signal,
      timeoutMs: 50,
    });
    expect(response.kind).toBe('result');
    expect(solveCalls).toBe(1);
  });

  it('surfaces a transport status without automatic retry', async () => {
    let calls = 0;
    const rpc: PtFlashRpc = {
      async discoverPtCapabilities() {
        return wireDiscovery();
      },
      async solvePtFlash() {
        calls += 1;
        throw new ConnectError('worker unavailable', Code.Unavailable);
      },
    };
    const client = new GrpcWebFlashClient(rpc, 'https://example.test');
    const backend = (await client.discoverPtCapabilities()).backends[0];
    expect(backend).toBeDefined();
    await expect(
      client.solvePtFlash(
        {
          configuredBackendId: 'pr76-default',
          pressurePa: 10e6,
          temperatureK: 350,
          feed: [
            { componentId: 'methane', moleFraction: 0.7 },
            { componentId: 'carbon-dioxide', moleFraction: 0.3 },
          ],
        },
        backend!,
      ),
    ).rejects.toMatchObject<Partial<PtTransportError>>({
      name: 'PtTransportError',
      code: Code.Unavailable,
    });
    expect(calls).toBe(1);
  });
});
