import { Buffer } from 'node:buffer';
import { randomUUID } from 'node:crypto';
import { Code } from '@connectrpc/connect';
import type { JsonValue } from '@bufbuild/protobuf';

import type { ModelSessionClient } from '../src/api/modelSessionClient';
import {
  MODEL_WORKBENCH_CONVENTION,
  type ModelWorkbenchReply,
  type ModelWorkbenchRequest,
} from '../src/api/modelWorkbenchContract';
import type { ModelValidationDetail } from '../src/api/modelValidationDetail';
import { MODEL_DESKTOP_V2_CONVENTION, type ModelDesktopReply } from '../src/api/modelDesktopContract';
import { validDesktopRequestId } from './desktopIpcPolicy';
import { ModelDesktopSession } from './modelDesktopSession';

function record(value: unknown): value is Record<string, unknown> {
  return value !== null && typeof value === 'object' && !Array.isArray(value);
}

function only(value: Record<string, unknown>, keys: readonly string[]): boolean {
  return Object.keys(value).every(key => keys.includes(key));
}

function fail(code: Code, reason: string, validation?: ModelValidationDetail): ModelWorkbenchReply {
  return {
    version: MODEL_WORKBENCH_CONVENTION,
    ok: false,
    error: { code, reason, ...(validation === undefined ? {} : { validation }) },
  };
}

function success(value: JsonValue): ModelWorkbenchReply {
  return { version: MODEL_WORKBENCH_CONVENTION, ok: true, value };
}

function parse(raw: unknown): ModelWorkbenchRequest | ModelWorkbenchReply {
  let value: unknown;
  try {
    const json = JSON.stringify(raw);
    if (json === undefined) return fail(Code.InvalidArgument, 'workbench.invalid_request');
    if (Buffer.byteLength(json, 'utf8') > 64 * 1024) {
      return fail(Code.ResourceExhausted, 'workbench.request_size');
    }
    value = JSON.parse(json);
  } catch {
    return fail(Code.InvalidArgument, 'workbench.invalid_request');
  }
  if (!record(value) || !only(value, ['version', 'requestId', 'operation', 'input']) ||
      !validDesktopRequestId(value.requestId)) {
    return fail(Code.InvalidArgument, 'workbench.invalid_request');
  }
  if (value.version !== MODEL_WORKBENCH_CONVENTION) {
    return fail(Code.InvalidArgument, 'workbench.unsupported_version');
  }
  if (typeof value.operation !== 'string' || !['apply', 'solve', 'release'].includes(value.operation)) {
    return fail(Code.InvalidArgument, 'workbench.unknown_operation');
  }
  const needsInput = value.operation === 'apply' || value.operation === 'solve';
  if (needsInput ? !record(value.input) : Object.hasOwn(value, 'input')) {
    return fail(Code.InvalidArgument, 'workbench.invalid_input');
  }
  return value as unknown as ModelWorkbenchRequest;
}

function desktopFailure(reply: ModelDesktopReply): ModelWorkbenchReply {
  if (reply.ok) return fail(Code.Internal, 'workbench.failed');
  return fail(reply.error.code as Code, reply.error.reason,
    'validation' in reply.error ? reply.error.validation : undefined);
}

const ownershipFailures = new Set([
  Code.Canceled, Code.Unknown, Code.DeadlineExceeded, Code.NotFound,
  Code.PermissionDenied, Code.FailedPrecondition, Code.Aborted, Code.Internal,
  Code.Unavailable, Code.DataLoss, Code.Unauthenticated,
]);

/**
 * One renderer document owns at most one applied model, but the renderer never
 * receives its session-local token. The low-level desktop session remains a
 * main-process implementation detail used for validation, transport and cleanup.
 */
export class ModelWorkbenchSession {
  #desktop: ModelDesktopSession;
  #currentModel: string | undefined;
  #activeRequest: string | undefined;
  #disposed = false;

  constructor(client: ModelSessionClient) {
    this.#desktop = new ModelDesktopSession(client);
  }

  async reset(): Promise<void> {
    this.#currentModel = undefined;
    await this.#desktop.reset();
  }

  async dispose(): Promise<void> {
    if (this.#disposed) return;
    this.#disposed = true;
    this.#currentModel = undefined;
    await this.#desktop.dispose();
  }

  cancel(id: unknown): void {
    if (!validDesktopRequestId(id) || id !== this.#activeRequest) return;
    this.#currentModel = undefined;
    this.#desktop.cancel(id);
    // Apply can be between the create reply and retirement of the previous model,
    // whose internal request ID is intentionally private. Reset cancels that whole
    // ownership epoch so a cancelled apply can never publish its replacement.
    void this.#desktop.reset().catch(() => {});
  }

  async #releaseToken(model: string): Promise<ModelDesktopReply> {
    return this.#desktop.invoke({
      version: MODEL_DESKTOP_V2_CONVENTION,
      requestId: randomUUID(),
      operation: 'release',
      model,
    }, MODEL_DESKTOP_V2_CONVENTION);
  }

  async invoke(raw: unknown): Promise<ModelWorkbenchReply> {
    if (this.#disposed) return fail(Code.FailedPrecondition, 'workbench.window_closed');
    const parsed = parse(raw);
    if ('ok' in parsed) return parsed;
    if (this.#activeRequest !== undefined) {
      return fail(
        this.#activeRequest === parsed.requestId ? Code.AlreadyExists : Code.ResourceExhausted,
        'workbench.busy',
      );
    }
    this.#activeRequest = parsed.requestId;
    try {
      if (parsed.operation === 'apply') {
        const created = await this.#desktop.invoke({
          version: MODEL_DESKTOP_V2_CONVENTION,
          requestId: parsed.requestId,
          operation: 'create',
          input: parsed.input,
        }, MODEL_DESKTOP_V2_CONVENTION);
        if (!created.ok) {
          await this.reset().catch(() => {});
          return desktopFailure(created);
        }
        if (!record(created.value) || Object.keys(created.value).length !== 2 ||
            typeof created.value.model !== 'string' || !Object.hasOwn(created.value, 'snapshot')) {
          await this.reset().catch(() => {});
          return fail(Code.DataLoss, 'workbench.failed');
        }
        const nextModel = created.value.model;
        const previous = this.#currentModel;
        if (previous !== undefined) {
          const retired = await this.#releaseToken(previous);
          if (!retired.ok) {
            // The new model must never become visible after ambiguous retirement.
            await this.#releaseToken(nextModel).catch(() => {});
            await this.reset().catch(() => {});
            return desktopFailure(retired);
          }
        }
        this.#currentModel = nextModel;
        return success(created.value.snapshot as JsonValue);
      }

      if (parsed.operation === 'solve') {
        const model = this.#currentModel;
        if (model === undefined) return fail(Code.NotFound, 'workbench.no_model');
        const solved = await this.#desktop.invoke({
          version: MODEL_DESKTOP_V2_CONVENTION,
          requestId: parsed.requestId,
          operation: 'solve',
          model,
          input: parsed.input,
        }, MODEL_DESKTOP_V2_CONVENTION);
        if (!solved.ok) {
          if (ownershipFailures.has(solved.error.code as Code)) {
            await this.reset().catch(() => {});
          }
          return desktopFailure(solved);
        }
        return success(solved.value);
      }

      const model = this.#currentModel;
      this.#currentModel = undefined;
      if (model === undefined) return success(null);
      const released = await this.#desktop.invoke({
        version: MODEL_DESKTOP_V2_CONVENTION,
        requestId: parsed.requestId,
        operation: 'release',
        model,
      }, MODEL_DESKTOP_V2_CONVENTION);
      if (!released.ok) {
        await this.reset().catch(() => {});
        return desktopFailure(released);
      }
      return success(null);
    } catch {
      await this.reset().catch(() => {});
      return fail(Code.Internal, 'workbench.failed');
    } finally {
      if (this.#activeRequest === parsed.requestId) this.#activeRequest = undefined;
    }
  }
}
