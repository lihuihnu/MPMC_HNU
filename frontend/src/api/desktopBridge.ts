import { Code } from '@connectrpc/connect';

import type {
  PtBackendDescriptor,
  PtCapabilityDiscovery,
  PtFlashRequest,
  PtFlashResponse,
} from '../domain/flash';
import {
  PT_DESKTOP_BRIDGE_CONVENTION,
  type DesktopBridgeFailure,
  type DesktopBridgeReply,
  type PtDesktopBridge,
} from './desktopBridgeContract';
import { PtTransportError, type FlashClient, type PtClientCallOptions } from './flashClient';
import { PtWireContractError } from './ptWire';

export {
  PT_DESKTOP_BRIDGE_CONVENTION,
  type DesktopBridgeFailure,
  type DesktopBridgeReply,
  type PtDesktopBridge,
} from './desktopBridgeContract';

function requestId(): string {
  return globalThis.crypto.randomUUID();
}

function transportFailure(error: DesktopBridgeFailure): PtTransportError {
  const code =
    Number.isInteger(error.code) && error.code >= Code.Canceled && error.code <= Code.Unauthenticated
      ? (error.code as Code)
      : Code.Unknown;
  return new PtTransportError(code, error.message);
}

async function invokeDesktop<T>(
  bridge: PtDesktopBridge,
  operation: (id: string) => Promise<DesktopBridgeReply<T>>,
  options?: PtClientCallOptions,
): Promise<T> {
  const id = requestId();
  if (options?.signal?.aborted === true) {
    throw new PtTransportError(Code.Canceled, 'The desktop PT RPC was canceled.');
  }
  const cancel = () => bridge.cancel(id);
  options?.signal?.addEventListener('abort', cancel, { once: true });
  try {
    const reply = await operation(id);
    if (!reply.ok) {
      throw transportFailure(reply.error);
    }
    return reply.value;
  } finally {
    options?.signal?.removeEventListener('abort', cancel);
  }
}

export class DesktopFlashClient implements FlashClient {
  readonly configured = true;
  readonly endpoint = 'desktop://local-pt-service';

  constructor(private readonly bridge: PtDesktopBridge) {
    if (bridge.convention !== PT_DESKTOP_BRIDGE_CONVENTION) {
      throw new PtWireContractError('Unsupported PT desktop bridge convention.');
    }
  }

  discoverPtCapabilities(options?: PtClientCallOptions): Promise<PtCapabilityDiscovery> {
    return invokeDesktop(
      this.bridge,
      (id) => this.bridge.discoverPtCapabilities(id),
      options,
    );
  }

  solvePtFlash(
    request: PtFlashRequest,
    expectedBackend: PtBackendDescriptor,
    options?: PtClientCallOptions,
  ): Promise<PtFlashResponse> {
    if (request.configuredBackendId !== expectedBackend.configuredBackendId) {
      throw new PtWireContractError(
        'Selected backend does not match the desktop capability snapshot.',
      );
    }
    return invokeDesktop(
      this.bridge,
      (id) => this.bridge.solvePtFlash(id, request),
      options,
    );
  }
}

export function desktopFlashClient(): FlashClient | null {
  return window.mpmcPtDesktop === undefined
    ? null
    : new DesktopFlashClient(window.mpmcPtDesktop);
}
