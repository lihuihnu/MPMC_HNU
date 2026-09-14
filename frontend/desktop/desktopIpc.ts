import { Code, ConnectError } from '@connectrpc/connect';
import { ipcMain, type WebContents } from 'electron';

import type {
  DesktopBridgeReply,
  DesktopBridgeFailure,
} from '../src/api/desktopBridgeContract';
import type {
  PtCapabilityDiscovery,
  PtFlashResponse,
} from '../src/domain/flash';
import { PtWireContractError } from '../src/api/ptWire';
import {
  parseDesktopSolveRequest,
  PtDesktopCallGate,
} from './desktopIpcPolicy';
import { PtDesktopGateway } from './ptGateway';

export const DESKTOP_DISCOVER_CHANNEL = 'mpmc:pt:discover:v1';
export const DESKTOP_SOLVE_CHANNEL = 'mpmc:pt:solve:v1';
export const DESKTOP_CANCEL_CHANNEL = 'mpmc:pt:cancel:v1';

function failure(code: Code, message: string): DesktopBridgeFailure {
  return { code, message };
}

function failureFrom(cause: unknown): DesktopBridgeFailure {
  if (cause instanceof ConnectError) {
    return failure(cause.code, cause.rawMessage || 'The desktop PT RPC failed.');
  }
  if (cause instanceof PtWireContractError) {
    return failure(Code.DataLoss, cause.message);
  }
  return failure(
    Code.Unavailable,
    cause instanceof Error ? cause.message : 'The desktop PT host is unavailable.',
  );
}

export function registerPtDesktopIpc(
  gateway: PtDesktopGateway,
  allowedWebContents: () => WebContents | null,
): () => void {
  const gate = new PtDesktopCallGate();

  const authorized = (sender: WebContents): boolean => sender === allowedWebContents();

  ipcMain.handle(
    DESKTOP_DISCOVER_CHANNEL,
    async (event, id: unknown): Promise<DesktopBridgeReply<PtCapabilityDiscovery>> => {
      if (!authorized(event.sender)) {
        return { ok: false, error: failure(Code.PermissionDenied, 'Desktop IPC sender rejected.') };
      }
      let controller: AbortController | undefined;
      try {
        controller = gate.begin(id);
        return { ok: true, value: await gateway.discover(controller.signal) };
      } catch (cause) {
        return { ok: false, error: failureFrom(cause) };
      } finally {
        gate.finish(id, controller);
      }
    },
  );

  ipcMain.handle(
    DESKTOP_SOLVE_CHANNEL,
    async (
      event,
      id: unknown,
      rawRequest: unknown,
    ): Promise<DesktopBridgeReply<PtFlashResponse>> => {
      if (!authorized(event.sender)) {
        return { ok: false, error: failure(Code.PermissionDenied, 'Desktop IPC sender rejected.') };
      }
      let controller: AbortController | undefined;
      try {
        controller = gate.begin(id);
        const request = parseDesktopSolveRequest(rawRequest);
        return { ok: true, value: await gateway.solve(request, controller.signal) };
      } catch (cause) {
        return { ok: false, error: failureFrom(cause) };
      } finally {
        gate.finish(id, controller);
      }
    },
  );

  ipcMain.on(DESKTOP_CANCEL_CHANNEL, (event, id: unknown) => {
    if (authorized(event.sender)) {
      gate.cancel(id);
    }
  });

  return () => {
    gate.abortAll();
    ipcMain.removeHandler(DESKTOP_DISCOVER_CHANNEL);
    ipcMain.removeHandler(DESKTOP_SOLVE_CHANNEL);
    ipcMain.removeAllListeners(DESKTOP_CANCEL_CHANNEL);
  };
}
