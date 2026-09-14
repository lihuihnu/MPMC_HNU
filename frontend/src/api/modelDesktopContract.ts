import type { JsonObject, JsonValue } from '@bufbuild/protobuf';

/** Additive to the frozen PT bridge. Payloads use canonical Protobuf JSON. */
export const MODEL_DESKTOP_CONVENTION = 'MPMC/model/desktop-bridge/v1' as const;
export const MODEL_DESKTOP_CHANNEL = 'mpmc:model:invoke:v1';
export const MODEL_DESKTOP_CANCEL_CHANNEL = 'mpmc:model:cancel:v1';
export type ModelDesktopOperation = 'connect' | 'reconnect' | 'create' | 'describe' | 'solve' | 'release';
export interface ModelDesktopRequest {
  version: typeof MODEL_DESKTOP_CONVENTION;
  requestId: string;
  operation: ModelDesktopOperation;
  model?: string;
  input?: JsonObject;
}
export type ModelDesktopReply =
  | { version: typeof MODEL_DESKTOP_CONVENTION; ok: true; value: JsonValue }
  | { version: typeof MODEL_DESKTOP_CONVENTION; ok: false; error: { code: number; reason: string } };
export interface ModelDesktopBridge {
  readonly convention: typeof MODEL_DESKTOP_CONVENTION;
  connect(requestId: string): Promise<ModelDesktopReply>;
  reconnect(requestId: string): Promise<ModelDesktopReply>;
  /** definition + exactly one explicit presetId/settings; no wireContract. */
  create(requestId: string, input: JsonObject): Promise<ModelDesktopReply>;
  describe(requestId: string, model: string): Promise<ModelDesktopReply>;
  /** pressurePa, temperatureK, feed in immutable component order. */
  solve(requestId: string, model: string, input: JsonObject): Promise<ModelDesktopReply>;
  release(requestId: string, model: string): Promise<ModelDesktopReply>;
  /** Cancelling an active call closes this window's session conservatively. */
  cancel(requestId: string): void;
}
declare global {
  interface Window { mpmcModelDesktop?: ModelDesktopBridge }
}
