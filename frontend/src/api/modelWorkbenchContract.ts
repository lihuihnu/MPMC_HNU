import type { JsonObject, JsonValue } from '@bufbuild/protobuf';
import type { ModelValidationDetail } from './modelValidationDetail';

/** Renderer-safe desktop contract: the main process owns sessions and model references. */
export const MODEL_WORKBENCH_CONVENTION = 'MPMC/model/workbench-bridge/v1' as const;
export const MODEL_WORKBENCH_CHANNEL = 'mpmc:model-workbench:invoke:v1';
export const MODEL_WORKBENCH_CANCEL_CHANNEL = 'mpmc:model-workbench:cancel:v1';

export type ModelWorkbenchOperation = 'apply' | 'solve' | 'release';

export interface ModelWorkbenchRequest {
  version: typeof MODEL_WORKBENCH_CONVENTION;
  requestId: string;
  operation: ModelWorkbenchOperation;
  input?: JsonObject;
}

export type ModelWorkbenchReply =
  | { version: typeof MODEL_WORKBENCH_CONVENTION; ok: true; value: JsonValue }
  | {
      version: typeof MODEL_WORKBENCH_CONVENTION;
      ok: false;
      error: { code: number; reason: string; validation?: ModelValidationDetail };
    };

/**
 * No session ID, native handle, model token, reconnect primitive, host credential,
 * or transport object is exposed to the renderer. Applying a model atomically
 * replaces the prior model owned by this window.
 */
export interface ModelWorkbenchBridge {
  readonly convention: typeof MODEL_WORKBENCH_CONVENTION;
  apply(requestId: string, input: JsonObject): Promise<ModelWorkbenchReply>;
  solve(requestId: string, input: JsonObject): Promise<ModelWorkbenchReply>;
  release(requestId: string): Promise<ModelWorkbenchReply>;
  cancel(requestId: string): void;
}

declare global {
  interface Window {
    mpmcModelWorkbench?: ModelWorkbenchBridge;
  }
}
