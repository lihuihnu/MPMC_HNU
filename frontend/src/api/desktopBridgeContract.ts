import type {
  PtCapabilityDiscovery,
  PtFlashRequest,
  PtFlashResponse,
} from '../domain/flash';

export const PT_DESKTOP_BRIDGE_CONVENTION = 'MPMC/PT/desktop-bridge/v1' as const;

export interface DesktopBridgeFailure {
  code: number;
  message: string;
}

export type DesktopBridgeReply<T> =
  | { ok: true; value: T }
  | { ok: false; error: DesktopBridgeFailure };

export interface PtDesktopBridge {
  readonly convention: typeof PT_DESKTOP_BRIDGE_CONVENTION;
  discoverPtCapabilities(
    requestId: string,
  ): Promise<DesktopBridgeReply<PtCapabilityDiscovery>>;
  solvePtFlash(
    requestId: string,
    request: PtFlashRequest,
  ): Promise<DesktopBridgeReply<PtFlashResponse>>;
  cancel(requestId: string): void;
}

declare global {
  interface Window {
    mpmcPtDesktop?: PtDesktopBridge;
  }
}
