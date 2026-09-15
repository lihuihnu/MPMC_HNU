import {
  clone,
  create,
  toJson,
  type JsonObject,
  type JsonValue,
} from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';

import type { ExpertModelOwner, ExpertOwnedModel } from './expertModelOwner';
import { inspectionFromSnapshot } from './expertModelInspector';
import type { ExpertModelSource } from './expertModelSource';
import type { ModelCallOptions, ModelCreateInput, ModelSolveInput } from './modelSessionClient';
import {
  MODEL_WORKBENCH_CONVENTION,
  type ModelWorkbenchBridge,
  type ModelWorkbenchReply,
} from './modelWorkbenchContract';
import { readModelValidationDetail } from './modelValidationDetail';
import {
  RendererModelError,
  invalidReply,
  readModelResult,
  readModelSnapshot,
  record,
} from './rendererModelWire';
import {
  CreateModelRequestSchema,
  FullPtResultSchema,
  ModelSnapshotSchema,
  SolveModelRequestSchema,
  type ModelSnapshot,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';

const maximumCallMs = 125_000;
const knownReasons = new Set([
  'ipc.invalid_request', 'ipc.request_size', 'ipc.unsupported_version', 'ipc.unknown_operation',
  'ipc.invalid_reference', 'ipc.invalid_input', 'ipc.window_closed', 'ipc.duplicate_request',
  'ipc.request_limit', 'ipc.document_changed', 'ipc.stale_reference', 'ipc.reply_size',
  'ipc.failed', 'ipc.sender_rejected', 'client.disposed', 'client.request_limit',
  'client.invalid_timeout', 'client.model_limit', 'session.closing', 'session.open_timeout',
  'session.invalid_handshake', 'session.open_failed', 'session.ended', 'session.closed',
  'session.reconnect_cancelled', 'session.not_connected', 'session.changed',
  'model.stale_reference', 'model.invalid_response', 'rpc.failed',
  'workbench.no_model', 'workbench.busy', 'workbench.failed', 'workbench.sender_rejected',
  'workbench.invalid_request', 'workbench.request_size', 'workbench.unsupported_version',
  'workbench.unknown_operation', 'workbench.invalid_input', 'workbench.window_closed',
  'workbench.document_changed', 'workbench.reply_size', 'workbench.stale_model',
]);
const ownershipErrors = new Set([
  Code.Canceled, Code.Unknown, Code.DeadlineExceeded, Code.NotFound,
  Code.PermissionDenied, Code.FailedPrecondition, Code.Aborted, Code.Internal,
  Code.Unavailable, Code.DataLoss, Code.Unauthenticated,
]);

function inputObject(value: JsonValue): JsonObject {
  if (value === null || typeof value !== 'object' || Array.isArray(value)) {
    throw new RendererModelError(Code.InvalidArgument, 'workbench.invalid_input', 'client');
  }
  return value;
}

/** Treat the isolated preload reply as untrusted data. */
function readWorkbenchReply(raw: unknown): JsonValue {
  try {
    const text = JSON.stringify(raw);
    if (!text || new TextEncoder().encode(text).length > 4 * 1024 * 1024) invalidReply();
  } catch {
    invalidReply();
  }
  if (!record(raw) || raw.version !== MODEL_WORKBENCH_CONVENTION) invalidReply();
  if (raw.ok === false && Object.keys(raw).length === 3 && record(raw.error)) {
    const keys = Object.keys(raw.error);
    if (!(keys.length === 2 || (keys.length === 3 && keys.includes('validation'))) ||
        !keys.includes('code') || !keys.includes('reason')) invalidReply();
    const { code, reason } = raw.error;
    if (typeof code !== 'number' || !Number.isInteger(code) || code < 1 || code > 16 ||
        typeof reason !== 'string' || !reason || reason.length > 128) invalidReply();
    const validation = readModelValidationDetail(raw.error.validation, code as Code);
    throw new RendererModelError(
      code as Code,
      knownReasons.has(reason) ? reason : 'workbench.failed',
      'ipc',
      validation,
    );
  }
  if (raw.ok !== true || Object.keys(raw).length !== 3 || !Object.hasOwn(raw, 'value') || raw.value === undefined) {
    invalidReply();
  }
  return raw.value as JsonValue;
}

function createJson(input: ModelCreateInput): JsonObject {
  try {
    const json = inputObject(toJson(CreateModelRequestSchema, create(CreateModelRequestSchema, input)));
    delete json.wireContract;
    return json;
  } catch (cause) {
    if (cause instanceof RendererModelError) throw cause;
    throw new RendererModelError(Code.InvalidArgument, 'workbench.invalid_input', 'client');
  }
}

function solveMessage(input: ModelSolveInput) {
  try {
    return create(SolveModelRequestSchema, input);
  } catch {
    throw new RendererModelError(Code.InvalidArgument, 'workbench.invalid_input', 'client');
  }
}

function solveJson(message: ReturnType<typeof solveMessage>): JsonObject {
  try {
    const json = inputObject(toJson(SolveModelRequestSchema, message));
    delete json.wireContract;
    delete json.modelHandle;
    return json;
  } catch (cause) {
    if (cause instanceof RendererModelError) throw cause;
    throw new RendererModelError(Code.InvalidArgument, 'workbench.invalid_input', 'client');
  }
}

/**
 * Renderer ownership facade for the local desktop workbench. The renderer keeps
 * only immutable snapshots plus a local generation number; the main process owns
 * the actual model/session lifetime and atomically retires replaced models.
 */
export class ModelWorkbenchOwner implements ExpertModelOwner {
  #generation = 0;
  #activeRequest: string | undefined;

  constructor(private readonly bridge: ModelWorkbenchBridge) {
    if (bridge.convention !== MODEL_WORKBENCH_CONVENTION ||
        !['apply', 'solve', 'release', 'cancel'].every(
          key => typeof (bridge as unknown as Record<string, unknown>)[key] === 'function',
        )) {
      throw new RendererModelError(Code.DataLoss, 'renderer.unsupported_bridge', 'contract');
    }
  }

  #invalidate(): void {
    ++this.#generation;
  }

  #current(generation: number): boolean {
    return generation === this.#generation;
  }

  #timeout(options: ModelCallOptions): number {
    const timeout = options.timeoutMs ?? maximumCallMs;
    if (!Number.isFinite(timeout) || timeout <= 0 || timeout > maximumCallMs) {
      throw new RendererModelError(Code.InvalidArgument, 'client.invalid_timeout', 'client');
    }
    if (options.signal?.aborted) {
      throw new RendererModelError(Code.Canceled, 'workbench.document_changed', 'client');
    }
    return timeout;
  }

  #invoke(
    operation: (requestId: string) => Promise<ModelWorkbenchReply>,
    options: ModelCallOptions,
  ): Promise<JsonValue> {
    const timeout = this.#timeout(options);
    if (this.#activeRequest !== undefined) {
      throw new RendererModelError(Code.ResourceExhausted, 'workbench.busy', 'client');
    }
    const requestId = globalThis.crypto.randomUUID();
    this.#activeRequest = requestId;
    let settled = false;
    let timer: ReturnType<typeof setTimeout>;
    let resolve!: (value: JsonValue) => void;
    let reject!: (error: RendererModelError) => void;
    const clear = () => {
      clearTimeout(timer);
      options.signal?.removeEventListener('abort', abort);
    };
    const finishError = (error: RendererModelError) => {
      if (!settled) {
        settled = true;
        clear();
        reject(error);
      }
    };
    const abort = () => {
      try { this.bridge.cancel(requestId); } catch { /* sanitized below */ }
      finishError(new RendererModelError(Code.Canceled, 'workbench.document_changed', 'client'));
    };
    const result = new Promise<JsonValue>((yes, no) => {
      resolve = yes;
      reject = no;
    });
    timer = setTimeout(() => {
      try { this.bridge.cancel(requestId); } catch { /* sanitized below */ }
      finishError(new RendererModelError(Code.DeadlineExceeded, 'renderer.timeout', 'client'));
    }, timeout);
    options.signal?.addEventListener('abort', abort, { once: true });

    let wire: Promise<ModelWorkbenchReply>;
    try {
      wire = Promise.resolve(operation(requestId));
    } catch {
      wire = Promise.reject(new Error('invoke failed'));
    }
    void wire.then(raw => {
      if (settled) return;
      try {
        const value = readWorkbenchReply(raw);
        settled = true;
        clear();
        resolve(value);
      } catch (cause) {
        finishError(cause instanceof RendererModelError
          ? cause
          : new RendererModelError(Code.DataLoss, 'renderer.invalid_reply', 'contract'));
      }
    }, () => finishError(
      new RendererModelError(Code.Unavailable, 'renderer.invoke_failed', 'transport'),
    )).finally(() => {
      if (this.#activeRequest === requestId) this.#activeRequest = undefined;
    });
    return result;
  }

  async create(input: ModelCreateInput, options: ModelCallOptions = {}): Promise<ExpertOwnedModel> {
    const json = createJson(input);
    let snapshot: ModelSnapshot;
    try {
      snapshot = readModelSnapshot(await this.#invoke(
        id => this.bridge.apply(id, json), options,
      ));
    } catch (cause) {
      // Main-process workbench contract tears down ownership after any failed apply.
      this.#invalidate();
      throw cause;
    }

    const generation = ++this.#generation;
    const storedSnapshot = clone(ModelSnapshotSchema, snapshot);
    let released = false;
    const authority = this;
    const stale = () => released || !authority.#current(generation);
    const source: ExpertModelSource = Object.freeze({
      describe(describeOptions = {}) {
        authority.#timeout(describeOptions);
        if (stale()) {
          return Promise.reject(new RendererModelError(Code.NotFound, 'workbench.stale_model', 'client'));
        }
        return Promise.resolve(inspectionFromSnapshot(storedSnapshot));
      },
    });

    const owned: ExpertOwnedModel = {
      snapshot: clone(ModelSnapshotSchema, storedSnapshot),
      source,
      get released() {
        return released;
      },
      async solve(solveInput: ModelSolveInput, solveOptions: ModelCallOptions = {}) {
        if (stale()) {
          throw new RendererModelError(Code.NotFound, 'workbench.stale_model', 'client');
        }
        const message = solveMessage(solveInput);
        const jsonInput = solveJson(message);
        try {
          const value = await authority.#invoke(
            id => authority.bridge.solve(id, jsonInput), solveOptions,
          );
          return clone(
            FullPtResultSchema,
            readModelResult(value, storedSnapshot, message).result,
          );
        } catch (cause) {
          if (cause instanceof RendererModelError && ownershipErrors.has(cause.code)) {
            released = true;
            if (authority.#current(generation)) authority.#invalidate();
          }
          throw cause;
        }
      },
      async release(releaseOptions: ModelCallOptions = {}) {
        if (released) return;
        released = true;
        // A successful later apply already retired this model in main. Never let
        // stale React cleanup release the replacement model.
        if (!authority.#current(generation)) return;
        try {
          await authority.#invoke(id => authority.bridge.release(id), releaseOptions);
        } finally {
          if (authority.#current(generation)) authority.#invalidate();
        }
      },
    };
    return Object.freeze(owned);
  }
}

export function modelWorkbenchOwner(
  bridge = typeof window === 'undefined' ? undefined : window.mpmcModelWorkbench,
): ExpertModelOwner | null {
  return bridge ? new ModelWorkbenchOwner(bridge) : null;
}
