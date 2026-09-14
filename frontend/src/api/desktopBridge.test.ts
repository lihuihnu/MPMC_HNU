import { Code } from '@connectrpc/connect';
import { describe, expect, it } from 'vitest';

import type { PtFlashRequest } from '../domain/flash';
import { wireDiscovery, wireResult } from '../test/ptWireFixtures';
import {
  DesktopFlashClient,
  PT_DESKTOP_BRIDGE_CONVENTION,
  type DesktopBridgeReply,
  type PtDesktopBridge,
} from './desktopBridge';
import { PtTransportError } from './flashClient';
import { mapDiscoveryResponse, mapSolveResponse } from './ptWire';

function request(): PtFlashRequest {
  return {
    configuredBackendId: 'pr76-default',
    pressurePa: 10e6,
    temperatureK: 350,
    feed: [
      { componentId: 'methane', moleFraction: 0.7 },
      { componentId: 'carbon-dioxide', moleFraction: 0.3 },
    ],
  };
}

describe('Electron desktop PT bridge client', () => {
  it('preserves discovery and solve result envelopes', async () => {
    const discovery = mapDiscoveryResponse(wireDiscovery());
    const backend = discovery.backends[0]!;
    const solveResponse = mapSolveResponse(wireResult(), backend, request());
    let solveCalls = 0;
    const bridge: PtDesktopBridge = {
      convention: PT_DESKTOP_BRIDGE_CONVENTION,
      async discoverPtCapabilities() {
        return { ok: true, value: discovery };
      },
      async solvePtFlash(_id, actual) {
        solveCalls += 1;
        expect(actual).toEqual(request());
        return { ok: true, value: solveResponse };
      },
      cancel() {},
    };
    const client = new DesktopFlashClient(bridge);

    await expect(client.discoverPtCapabilities()).resolves.toEqual(discovery);
    await expect(client.solvePtFlash(request(), backend)).resolves.toEqual(
      solveResponse,
    );
    expect(solveCalls).toBe(1);
    expect(client.endpoint).toBe('desktop://local-pt-service');
  });

  it('maps desktop transport failure without creating a service result', async () => {
    const discovery = mapDiscoveryResponse(wireDiscovery());
    const bridge: PtDesktopBridge = {
      convention: PT_DESKTOP_BRIDGE_CONVENTION,
      async discoverPtCapabilities() {
        return { ok: true, value: discovery };
      },
      async solvePtFlash() {
        return {
          ok: false,
          error: { code: Code.Unavailable, message: 'desktop child exited' },
        };
      },
      cancel() {},
    };
    const client = new DesktopFlashClient(bridge);

    await expect(
      client.solvePtFlash(request(), discovery.backends[0]!),
    ).rejects.toMatchObject<Partial<PtTransportError>>({
      name: 'PtTransportError',
      code: Code.Unavailable,
      message: 'desktop child exited',
    });
  });

  it('forwards AbortSignal cancellation to the exact desktop call', async () => {
    const discovery = mapDiscoveryResponse(wireDiscovery());
    let pending:
      | ((reply: DesktopBridgeReply<typeof discovery>) => void)
      | undefined;
    let pendingId = '';
    const bridge: PtDesktopBridge = {
      convention: PT_DESKTOP_BRIDGE_CONVENTION,
      discoverPtCapabilities(id) {
        pendingId = id;
        return new Promise((resolve) => {
          pending = resolve;
        });
      },
      async solvePtFlash() {
        throw new Error('unexpected solve');
      },
      cancel(id) {
        expect(id).toBe(pendingId);
        pending?.({
          ok: false,
          error: { code: Code.Canceled, message: 'desktop call canceled' },
        });
      },
    };
    const client = new DesktopFlashClient(bridge);
    const controller = new AbortController();
    const call = client.discoverPtCapabilities({ signal: controller.signal });
    controller.abort();

    await expect(call).rejects.toMatchObject({ code: Code.Canceled });
  });
});
