import { Buffer } from 'node:buffer';
import { randomUUID } from 'node:crypto';
import { fromJson, toJson, type JsonValue } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';
import { ModelClientError, type ModelReference, type ModelSessionClient } from '../src/api/modelSessionClient';
import { MODEL_DESKTOP_CONVENTION, MODEL_DESKTOP_V2_CONVENTION, type ModelDesktopVersion, type ModelDesktopReply, type ModelDesktopRequest } from '../src/api/modelDesktopContract';
import { CreateModelRequestSchema, FullPtResultSchema, ModelSnapshotSchema, SolveModelRequestSchema } from '../src/gen/mpmc/model_configuration/v1/model_service_pb';
import { readModelValidationDetail, type ModelValidationDetail } from '../src/api/modelValidationDetail';
import { validDesktopRequestId } from './desktopIpcPolicy';

type Client = Pick<ModelSessionClient, 'connected' | 'connect' | 'disconnect' | 'dispose' | 'create' | 'describe' | 'solve' | 'release'>;
function fail(code: Code, reason: string): never { throw new ModelClientError(code, reason); }
export function modelDesktopFailure(code: Code, reason: string, version: ModelDesktopVersion = MODEL_DESKTOP_CONVENTION,
  detail?: ModelValidationDetail): ModelDesktopReply {
  const validation = readModelValidationDetail(detail, code);
  if (version === MODEL_DESKTOP_V2_CONVENTION) return { version, ok: false, error: { code, reason, ...(validation ? { validation } : {}) } };
  return { version: MODEL_DESKTOP_CONVENTION, ok: false, error: { code, reason } };
}
function record(value: unknown): value is Record<string, unknown> {
  return value !== null && typeof value === 'object' && !Array.isArray(value);
}
function only(value: Record<string, unknown>, keys: readonly string[]): boolean {
  return Object.keys(value).every(key => keys.includes(key));
}
function parse(raw: unknown, version: ModelDesktopVersion): ModelDesktopRequest {
  let value: unknown;
  try {
    const json = JSON.stringify(raw);
    if (json === undefined) fail(Code.InvalidArgument, 'ipc.invalid_request');
    if (Buffer.byteLength(json, 'utf8') > 64 * 1024) fail(Code.ResourceExhausted, 'ipc.request_size');
    value = JSON.parse(json);
  } catch (cause) {
    if (cause instanceof ModelClientError) throw cause;
    fail(Code.InvalidArgument, 'ipc.invalid_request');
  }
  if (!record(value) || !only(value, ['version', 'requestId', 'operation', 'model', 'input']) ||
      !validDesktopRequestId(value.requestId)) fail(Code.InvalidArgument, 'ipc.invalid_request');
  if (value.version !== version) fail(Code.InvalidArgument, 'ipc.unsupported_version');
  const op = value.operation;
  if (typeof op !== 'string' || !['connect', 'reconnect', 'create', 'describe', 'solve', 'release'].includes(op)) {
    fail(Code.InvalidArgument, 'ipc.unknown_operation');
  }
  const needsModel = op === 'describe' || op === 'solve' || op === 'release';
  const needsInput = op === 'create' || op === 'solve';
  if (needsModel ? typeof value.model !== 'string' || !/^[a-f0-9-]{36}$/u.test(value.model) : 'model' in value) {
    fail(Code.InvalidArgument, 'ipc.invalid_reference');
  }
  if (needsInput ? !record(value.input) : 'input' in value) fail(Code.InvalidArgument, 'ipc.invalid_input');
  if (record(value.input) && !only(value.input, op === 'create'
    ? ['definition', 'settings', 'presetId'] : ['pressurePa', 'temperatureK', 'feed'])) {
    fail(Code.InvalidArgument, 'ipc.invalid_input');
  }
  return value as unknown as ModelDesktopRequest;
}

/** One attached window/document. No server handles or credentials cross IPC. */
export class ModelDesktopSession {
  #references = new Map<string, ModelReference>();
  #calls = new Map<string, AbortController>();
  #generation = 0;
  #disposed = false;
  #draining = false;
  #drain: Promise<void> = Promise.resolve();
  constructor(private readonly client: Client) {}

  /** Invalidates synchronously; remote cleanup and in-flight drain are asynchronous. */
  reset(): Promise<void> {
    ++this.#generation;
    this.#references.clear();
    for (const controller of this.#calls.values()) controller.abort();
    // Keep request reservations until each old operation actually returns.
    if (!this.#draining) {
      this.#draining = true;
      this.#drain = this.client.disconnect().finally(() => { this.#draining = false; });
    }
    return this.#drain;
  }
  async dispose(): Promise<void> {
    this.#disposed = true;
    await this.reset();
    await this.client.dispose();
  }
  cancel(id: unknown): void {
    if (validDesktopRequestId(id) && this.#calls.has(id)) void this.reset();
  }
  async invoke(raw: unknown, version: ModelDesktopVersion = MODEL_DESKTOP_CONVENTION): Promise<ModelDesktopReply> {
    let id: string | undefined;
    let controller: AbortController | undefined;
    try {
      if (this.#disposed) fail(Code.FailedPrecondition, 'ipc.window_closed');
      const request = parse(raw, version);
      if (this.#calls.has(request.requestId)) fail(Code.AlreadyExists, 'ipc.duplicate_request');
      if (this.#calls.size >= 4) fail(Code.ResourceExhausted, 'ipc.request_limit');
      // A reconnect cancels other work; it never recreates existing models.
      if (request.operation === 'reconnect') void this.reset();
      id = request.requestId; controller = new AbortController();
      this.#calls.set(id, controller);
      const generation = this.#generation;
      const current = () => {
        if (this.#disposed || generation !== this.#generation || controller!.signal.aborted) {
          fail(Code.Canceled, 'ipc.document_changed');
        }
      };
      await this.#drain; current();
      let value: JsonValue = null;
      const options = { signal: controller.signal };
      const reference = () => {
        const found = this.#references.get(request.model!);
        if (!found) fail(Code.NotFound, 'ipc.stale_reference');
        return found;
      };
      switch (request.operation) {
        case 'connect': case 'reconnect':
          if (!this.client.connected) this.#references.clear();
          await this.client.connect();
          break;
        case 'create': {
          // Decode before opening/dispatch. Unknown nested fields are rejected.
          let input;
          try { input = fromJson(CreateModelRequestSchema, request.input!); }
          catch { fail(Code.InvalidArgument, 'ipc.invalid_input'); }
          if (!this.client.connected) this.#references.clear();
          await this.client.connect(); current();
          const made = await this.client.create(input, options); current();
          const token = randomUUID();
          this.#references.set(token, made.model);
          value = { model: token, snapshot: toJson(ModelSnapshotSchema, made.snapshot) };
          break;
        }
        case 'describe':
          value = toJson(ModelSnapshotSchema, await this.client.describe(reference(), options));
          break;
        case 'solve': {
          const model = reference();
          let input;
          try { input = fromJson(SolveModelRequestSchema, request.input!); }
          catch { fail(Code.InvalidArgument, 'ipc.invalid_input'); }
          value = toJson(FullPtResultSchema, await this.client.solve(model, input, options));
          break;
        }
        case 'release': {
          const model = reference();
          this.#references.delete(request.model!);
          await this.client.release(model, options);
          break;
        }
      }
      current();
      const reply = { version, ok: true as const, value };
      if (Buffer.byteLength(JSON.stringify(reply), 'utf8') > 4 * 1024 * 1024) {
        void this.reset(); fail(Code.ResourceExhausted, 'ipc.reply_size');
      }
      return reply;
    } catch (cause) {
      if (cause instanceof ModelClientError) return modelDesktopFailure(cause.code, cause.reason, version, cause.validation);
      // Includes serialization/entropy failures after an ambiguous Create.
      void this.reset();
      return modelDesktopFailure(Code.Internal, 'ipc.failed', version);
    } finally {
      if (id !== undefined && this.#calls.get(id) === controller) this.#calls.delete(id);
      if (!this.client.connected) this.#references.clear();
    }
  }
}
