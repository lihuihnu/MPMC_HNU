import type { MessageInitShape } from '@bufbuild/protobuf';
import { Code, ConnectError, type CallOptions, type Client } from '@connectrpc/connect';
import { ModelServiceErrorSchema } from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { MODEL_VALIDATION_DETAIL_VERSION, readModelValidationDetail, type ModelValidationDetail } from './modelValidationDetail';
import type {
  CreateModelRequestSchema,
  FullPtResult,
  ModelConfigurationService,
  ModelSessionService,
  ModelSnapshot,
  SolveModelRequestSchema,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';

export const MODEL_WIRE_CONTRACT = 'mpmc.model_configuration.v1/model-service/v1';
export const MODEL_SESSION_CONTRACT = 'mpmc.model_configuration.v1/model-session/v1';
export const MODEL_SESSION_HEADER = 'x-mpmc-model-session';

declare const referenceBrand: unique symbol;
/** In-memory capability owned by exactly one client/session; never a wire handle. */
export interface ModelReference { readonly [referenceBrand]: true }
// Preserve Protobuf-ES's full-message / initializer union under strict optionals.
type OmitEach<T, K extends PropertyKey> = T extends unknown ? Omit<T, K> : never;
export type ModelCreateInput = OmitEach<MessageInitShape<typeof CreateModelRequestSchema>, 'wireContract'>;
export type ModelSolveInput = OmitEach<MessageInitShape<typeof SolveModelRequestSchema>, 'wireContract' | 'modelHandle'>;
export interface ModelCallOptions { signal?: AbortSignal; timeoutMs?: number }
export interface ModelSessionConnection {
  models: Client<typeof ModelConfigurationService>;
  sessions: Client<typeof ModelSessionService>;
  headers?: HeadersInit;
  /** Synchronous, idempotent transport abort. Must close owned sockets/streams. */
  close(): void;
}
/** Must settle on cancellation, or close any transport acquired after cancellation. */
export type ModelSessionConnector = (signal: AbortSignal) => Promise<ModelSessionConnection>;

export class ModelClientError extends Error {
  constructor(readonly code: Code, readonly reason: string, readonly validation?: ModelValidationDetail) {
    super(`Model client: ${reason}`);
    this.name = 'ModelClientError';
  }
}
function failure(code: Code, reason: string, validation?: ModelValidationDetail): ModelClientError {
  return new ModelClientError(code, reason, validation);
}
function validationDetail(error: ConnectError, privateValues: readonly string[]): ModelValidationDetail | undefined {
  // Only incoming binary details are accepted here; bound work before decoding metadata.
  if (error.details.length > 8) return undefined;
  let bytes = 0;
  for (const detail of error.details) {
    if (!('type' in detail) || !(detail.value instanceof Uint8Array)) return undefined;
    bytes += detail.value.byteLength;
    if (detail.value.byteLength > 4096 || bytes > 8192) return undefined;
  }
  const matching = error.details.filter(detail => 'type' in detail && detail.type === ModelServiceErrorSchema.typeName);
  if (matching.length !== 1) return undefined;
  const detail = error.findDetails(ModelServiceErrorSchema)[0];
  if (!detail || detail.wireContract !== MODEL_WIRE_CONTRACT) return undefined;
  return readModelValidationDetail({ version: MODEL_VALIDATION_DETAIL_VERSION, code: detail.code,
    ...(detail.field === undefined ? {} : { field: detail.field }) }, error.code, privateValues);
}
function deferred() {
  let resolve!: () => void;
  let reject!: (reason: unknown) => void;
  const promise = new Promise<void>((yes, no) => { resolve = yes; reject = no; });
  // The pump owns this promise even when its caller abandons an opening.
  void promise.catch(() => {});
  return { promise, resolve, reject };
}
interface Epoch {
  valid: boolean;
  finished: boolean;
  controller: AbortController;
  ready: ReturnType<typeof deferred>;
  done: Promise<void>;
  connection?: ModelSessionConnection;
  id?: string;
  models: Map<ModelReference, string>;
  reservations: number;
  calls: Set<Promise<unknown>>;
}
const ambiguous = new Set([Code.Canceled, Code.DeadlineExceeded, Code.Unavailable,
  Code.Unknown, Code.Internal, Code.DataLoss, Code.Unauthenticated]);

/** Owns one streaming lease. Reconnection is explicit and never retries Create. */
export class ModelSessionClient {
  #epoch: Epoch | undefined;
  #disposed = false;
  #restarting: Promise<void> | undefined;
  #restartToken: object | undefined;

  constructor(private readonly connector: ModelSessionConnector) {}

  get connected(): boolean { return this.#epoch?.valid === true && this.#epoch.id !== undefined; }

  connect(): Promise<void> { return this.#restarting ?? this.#open(); }

  #open(): Promise<void> {
    if (this.#disposed) return Promise.reject(failure(Code.FailedPrecondition, 'client.disposed'));
    if (this.#epoch?.valid) return this.#epoch.ready.promise;
    if (this.#epoch && !this.#epoch.finished) {
      return Promise.reject(failure(Code.FailedPrecondition, 'session.closing'));
    }
    const epoch: Epoch = { valid: true, finished: false, controller: new AbortController(),
      ready: deferred(), done: Promise.resolve(), models: new Map(), reservations: 0, calls: new Set() };
    this.#epoch = epoch;
    epoch.done = this.#pump(epoch);
    return epoch.ready.promise;
  }

  async #pump(epoch: Epoch): Promise<void> {
    const timer = setTimeout(() => this.#invalidate(epoch,
      failure(Code.DeadlineExceeded, 'session.open_timeout')), 10_000);
    try {
      const connection = await this.connector(epoch.controller.signal);
      epoch.connection = connection;
      if (!epoch.valid) return;
      const headers = new Headers(connection.headers);
      headers.delete(MODEL_SESSION_HEADER);
      for await (const opened of connection.sessions.openModelSession(
        { wireContract: MODEL_SESSION_CONTRACT },
        { headers, signal: epoch.controller.signal, timeoutMs: 1_740_000 },
      )) {
        if (!epoch.valid) break;
        if (epoch.id !== undefined || opened.wireContract !== MODEL_SESSION_CONTRACT ||
            !opened.sessionId || opened.sessionId.length > 256) {
          throw failure(Code.DataLoss, 'session.invalid_handshake');
        }
        epoch.id = opened.sessionId;
        clearTimeout(timer);
        epoch.ready.resolve();
      }
    } catch (cause) {
      const code = cause instanceof ModelClientError ? cause.code : ConnectError.from(cause).code;
      epoch.ready.reject(failure(code, 'session.open_failed'));
    } finally {
      clearTimeout(timer);
      this.#invalidate(epoch, failure(Code.Unavailable, 'session.ended'));
      // An aborted connector may resolve late. Its transport still belongs here.
      epoch.connection?.close();
      await Promise.allSettled([...epoch.calls]);
      epoch.finished = true;
    }
  }

  #invalidate(epoch: Epoch, error = failure(Code.Canceled, 'session.closed')): void {
    epoch.valid = false;
    epoch.models.clear();
    epoch.ready.reject(error);
    epoch.controller.abort();
    epoch.connection?.close();
  }

  disconnect(): Promise<void> {
    this.#restartToken = undefined;
    const epoch = this.#epoch;
    if (!epoch) return Promise.resolve();
    this.#invalidate(epoch);
    return epoch.done;
  }

  reconnect(): Promise<void> {
    if (this.#disposed) return Promise.reject(failure(Code.FailedPrecondition, 'client.disposed'));
    if (this.#restarting) return this.#restarting;
    const closed = this.disconnect();
    const token = {};
    this.#restartToken = token;
    const pending = closed.then(() => {
      if (this.#restartToken !== token) throw failure(Code.Canceled, 'session.reconnect_cancelled');
      return this.#open();
    }).finally(() => {
      if (this.#restarting === pending) this.#restarting = undefined;
    });
    this.#restarting = pending;
    return pending;
  }

  dispose(): Promise<void> {
    this.#disposed = true;
    return this.disconnect();
  }

  #active(): Epoch {
    if (this.#disposed) throw failure(Code.FailedPrecondition, 'client.disposed');
    const epoch = this.#epoch;
    if (!epoch?.valid || !epoch.id || !epoch.connection) {
      throw failure(Code.FailedPrecondition, 'session.not_connected');
    }
    return epoch;
  }
  #handle(epoch: Epoch, ref: ModelReference): string {
    const handle = epoch.models.get(ref);
    if (!handle) throw failure(Code.NotFound, 'model.stale_reference');
    return handle;
  }

  #call<T>(epoch: Epoch, options: ModelCallOptions, work: (options: CallOptions) => Promise<T>,
    mutation = false): Promise<T> {
    if (epoch.calls.size >= 4) return Promise.reject(failure(Code.ResourceExhausted, 'client.request_limit'));
    const timeoutMs = options.timeoutMs ?? 110_000;
    if (!Number.isFinite(timeoutMs) || timeoutMs <= 0 || timeoutMs > 110_000) {
      return Promise.reject(failure(Code.InvalidArgument, 'client.invalid_timeout'));
    }
    const headers = new Headers(epoch.connection?.headers);
    headers.set(MODEL_SESSION_HEADER, epoch.id!);
    const signal = options.signal
      ? AbortSignal.any([epoch.controller.signal, options.signal]) : epoch.controller.signal;
    const pending = (async () => {
      try {
        const result = await work({ headers, signal, timeoutMs });
        if (!epoch.valid || this.#epoch !== epoch) throw failure(Code.Canceled, 'session.changed');
        return result;
      } catch (cause) {
        const error = ConnectError.from(cause);
        const code = cause instanceof ModelClientError ? cause.code : error.code;
        // Capture private routing values before mutation failure clears the epoch's model map.
        const privateValues = [epoch.id ?? '', ...epoch.models.values(), ...headers.values()];
        const authorization = headers.get('authorization');
        if (authorization?.startsWith('Bearer ')) privateValues.push(authorization.slice(7));
        const validation = cause instanceof ModelClientError ? undefined : validationDetail(error, privateValues);
        if (mutation || ambiguous.has(code)) this.#invalidate(epoch);
        throw failure(code, 'rpc.failed', validation);
      }
    })();
    epoch.calls.add(pending);
    void pending.finally(() => epoch.calls.delete(pending)).catch(() => {});
    return pending;
  }

  async create(input: ModelCreateInput, options: ModelCallOptions = {}): Promise<{ model: ModelReference; snapshot: ModelSnapshot }> {
    const epoch = this.#active();
    if (epoch.models.size + epoch.reservations >= 4) throw failure(Code.ResourceExhausted, 'client.model_limit');
    ++epoch.reservations;
    try {
      const response = await this.#call(epoch, options, call => epoch.connection!.models.createModel(
        { ...input, wireContract: MODEL_WIRE_CONTRACT }, call), true);
      if (response.wireContract !== MODEL_WIRE_CONTRACT || !response.modelHandle ||
          response.modelHandle.length > 256 || !response.snapshot) {
        this.#invalidate(epoch);
        throw failure(Code.DataLoss, 'model.invalid_response');
      }
      if (!epoch.valid) throw failure(Code.Canceled, 'session.changed');
      const model = Object.freeze({}) as ModelReference;
      epoch.models.set(model, response.modelHandle);
      return { model, snapshot: response.snapshot };
    } finally { --epoch.reservations; }
  }

  async describe(model: ModelReference, options: ModelCallOptions = {}): Promise<ModelSnapshot> {
    const epoch = this.#active(); const modelHandle = this.#handle(epoch, model);
    const response = await this.#call(epoch, options, call => epoch.connection!.models.describeModel(
      { wireContract: MODEL_WIRE_CONTRACT, modelHandle }, call));
    if (!epoch.valid) throw failure(Code.Canceled, 'session.changed');
    if (response.wireContract !== MODEL_WIRE_CONTRACT || !response.snapshot) {
      this.#invalidate(epoch); throw failure(Code.DataLoss, 'model.invalid_response');
    }
    return response.snapshot;
  }

  async solve(model: ModelReference, input: ModelSolveInput, options: ModelCallOptions = {}): Promise<FullPtResult> {
    const epoch = this.#active(); const modelHandle = this.#handle(epoch, model);
    const response = await this.#call(epoch, options, call => epoch.connection!.models.solveModel(
      { ...input, wireContract: MODEL_WIRE_CONTRACT, modelHandle }, call));
    if (!epoch.valid) throw failure(Code.Canceled, 'session.changed');
    if (response.wireContract !== MODEL_WIRE_CONTRACT || !response.result) {
      this.#invalidate(epoch); throw failure(Code.DataLoss, 'model.invalid_response');
    }
    return response.result;
  }

  async release(model: ModelReference, options: ModelCallOptions = {}): Promise<void> {
    const epoch = this.#active(); const modelHandle = this.#handle(epoch, model);
    // Invalidate before yielding: concurrent callers cannot reuse a released ref.
    epoch.models.delete(model);
    try {
      const response = await this.#call(epoch, options, call => epoch.connection!.models.releaseModel(
        { wireContract: MODEL_WIRE_CONTRACT, modelHandle }, call), true);
      if (response.wireContract !== MODEL_WIRE_CONTRACT) throw failure(Code.DataLoss, 'model.invalid_response');
    } catch (cause) {
      this.#invalidate(epoch);
      throw cause;
    }
  }
}
