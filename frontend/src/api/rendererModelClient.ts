import { clone, create, equals, toJson, type JsonValue, type JsonObject } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';
import { MODEL_DESKTOP_CONVENTION, type ModelDesktopBridge } from './modelDesktopContract';
import type { ModelCreateInput, ModelSolveInput } from './modelSessionClient';
import { CreateModelRequestSchema, ModelSnapshotSchema, SolveModelRequestSchema, type ModelSnapshot } from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { RendererModelError, invalidReply, readModelReply, readModelResult, readModelSnapshot, record, type RendererModelResult } from './rendererModelWire';

export { RendererModelError, type RendererModelErrorCategory, type RendererModelResult } from './rendererModelWire';
export type { ModelCreateInput, ModelSolveInput } from './modelSessionClient';
declare const brand: unique symbol;
export interface RendererModelReference { readonly [brand]: true }
export interface RendererModelCallOptions { signal?: AbortSignal; timeoutMs?: number }
interface OwnedModel { token: string; snapshot: ModelSnapshot }
interface Pending { id: string; generation: number; reject(error: RendererModelError): void; done: Promise<void> }
const instances = new WeakMap<ModelDesktopBridge, RendererModelClient>();
const sessionErrors = new Set([Code.Canceled, Code.Unknown, Code.DeadlineExceeded, Code.NotFound,
  Code.PermissionDenied, Code.FailedPrecondition, Code.Aborted, Code.Internal, Code.Unavailable,
  Code.DataLoss, Code.Unauthenticated]);
function inputObject(value: JsonValue): JsonObject {
  if (value === null || typeof value !== 'object' || Array.isArray(value)) {
    throw new RendererModelError(Code.InvalidArgument, 'renderer.invalid_input');
  }
  return value;
}
const limitMs = 125_000; // Existing main handshake (10 s) + unary (110 s), with transport headroom.

/** Window-scoped adapter; main owns the native lease and reclaims it on window teardown. */
export class RendererModelClient {
  #models = new Map<RendererModelReference, OwnedModel>();
  #pending = new Set<Pending>();
  #reservations = 0;
  #generation = 0;
  #needsReconnect = false;
  #drainBlocked = false;
  #reconnecting: Promise<void> | undefined;
  private constructor(private readonly bridge: ModelDesktopBridge) {}

  static fromBridge(bridge: ModelDesktopBridge): RendererModelClient {
    if (bridge.convention !== MODEL_DESKTOP_CONVENTION ||
        !['connect', 'reconnect', 'create', 'describe', 'solve', 'release', 'cancel'].every(key =>
          typeof (bridge as unknown as Record<string, unknown>)[key] === 'function')) {
      throw new RendererModelError(Code.DataLoss, 'renderer.unsupported_bridge', 'contract');
    }
    let client = instances.get(bridge);
    if (!client) { client = new RendererModelClient(bridge); instances.set(bridge, client); }
    return client;
  }
  get requiresReconnect(): boolean { return this.#needsReconnect; }
  #ready(): void {
    if (this.#reconnecting || this.#needsReconnect) {
      throw new RendererModelError(Code.FailedPrecondition, 'renderer.reconnect_required');
    }
  }
  #owned(model: RendererModelReference): OwnedModel {
    const entry = this.#models.get(model);
    if (!entry) throw new RendererModelError(Code.NotFound, 'renderer.stale_reference');
    return entry;
  }
  #invalidate(except?: Pending): void {
    ++this.#generation; this.#models.clear(); this.#needsReconnect = true;
    for (const pending of this.#pending) {
      if (pending !== except) pending.reject(new RendererModelError(Code.Canceled, 'renderer.session_changed'));
      try { this.bridge.cancel(pending.id); } catch { /* Invoke failure is sanitized separately. */ }
    }
  }
  #timeout(options: RendererModelCallOptions): number {
    const timeout = options.timeoutMs ?? limitMs;
    if (!Number.isFinite(timeout) || timeout <= 0 || timeout > limitMs) {
      throw new RendererModelError(Code.InvalidArgument, 'renderer.invalid_timeout');
    }
    if (options.signal?.aborted) throw new RendererModelError(Code.Canceled, 'renderer.cancelled');
    return timeout;
  }
  #invoke<T>(operation: (id: string) => Promise<unknown>, decode: (value: JsonValue) => T,
    options: RendererModelCallOptions, mutation = false): Promise<T> {
    const timeout = this.#timeout(options);
    if (this.#pending.size >= 4) throw new RendererModelError(Code.ResourceExhausted, 'renderer.request_limit');
    const id = globalThis.crypto.randomUUID();
    let finish!: () => void;
    const done = new Promise<void>(yes => { finish = yes; });
    let reject!: (error: RendererModelError) => void;
    let resolve!: (value: T) => void;
    let settled = false;
    let timer: ReturnType<typeof setTimeout>;
    const clear = () => { clearTimeout(timer); options.signal?.removeEventListener('abort', abort); };
    const promise = new Promise<T>((yes, no) => {
      resolve = value => { if (!settled) { settled = true; clear(); yes(value); } };
      reject = error => { if (!settled) { settled = true; clear(); no(error); } };
    });
    const pending: Pending = { id, generation: this.#generation, reject, done };
    const abort = () => {
      this.#invalidate(pending); reject(new RendererModelError(Code.Canceled, 'renderer.cancelled'));
    };
    timer = setTimeout(() => {
      this.#invalidate(pending); reject(new RendererModelError(Code.DeadlineExceeded, 'renderer.timeout'));
    }, timeout);
    options.signal?.addEventListener('abort', abort, { once: true });
    this.#pending.add(pending);
    // Keep admission until the actual IPC promise settles, even after local cancellation.
    let wire: Promise<unknown>;
    try { wire = Promise.resolve(operation(id)); }
    catch { wire = Promise.reject(new RendererModelError(Code.Unavailable, 'renderer.invoke_failed', 'transport')); }
    void wire.then(raw => {
      if (settled || pending.generation !== this.#generation) return;
      try { resolve(decode(readModelReply(raw))); }
      catch (cause) {
        const error = cause instanceof RendererModelError ? cause
          : new RendererModelError(Code.DataLoss, 'renderer.invalid_reply', 'contract');
        if (mutation || sessionErrors.has(error.code)) this.#invalidate(pending);
        reject(error);
      }
    }, () => {
      if (!settled) {
        this.#invalidate(pending);
        reject(new RendererModelError(Code.Unavailable, 'renderer.invoke_failed', 'transport'));
      }
    }).finally(() => { this.#pending.delete(pending); finish(); });
    return promise;
  }
  #empty(value: JsonValue): void { if (value !== null) invalidReply(); }

  async connect(options: RendererModelCallOptions = {}): Promise<void> {
    this.#ready(); await this.#invoke(id => this.bridge.connect(id), value => this.#empty(value), options);
  }
  reconnect(options: RendererModelCallOptions = {}): Promise<void> {
    if (this.#reconnecting) return this.#reconnecting;
    if (this.#drainBlocked && this.#pending.size) {
      return Promise.reject(new RendererModelError(Code.FailedPrecondition, 'renderer.calls_draining'));
    }
    this.#drainBlocked = false;
    let timeout: number;
    try { timeout = this.#timeout(options); } catch (cause) { return Promise.reject(cause); }
    this.#invalidate();
    const drains = [...this.#pending].map(pending => pending.done);
    const work = (async () => {
      // A broken preload cannot grow a queue or keep reconnect pending indefinitely.
      await new Promise<void>((resolve, reject) => {
        let settled = false;
        const clear = () => { clearTimeout(timer); options.signal?.removeEventListener('abort', abort); };
        const fail = (error: RendererModelError) => {
          if (!settled) { settled = true; this.#drainBlocked = this.#pending.size > 0; clear(); reject(error); }
        };
        const abort = () => fail(new RendererModelError(Code.Canceled, 'renderer.cancelled'));
        const timer = setTimeout(() => fail(new RendererModelError(Code.DeadlineExceeded, 'renderer.drain_timeout')), timeout);
        options.signal?.addEventListener('abort', abort, { once: true });
        void Promise.all(drains).then(() => {
          if (options.signal?.aborted) abort();
          else if (!settled) { settled = true; clear(); resolve(); }
        });
      });
      await this.#invoke(id => this.bridge.reconnect(id), value => this.#empty(value), options, true);
      this.#needsReconnect = false;
    })().finally(() => { if (this.#reconnecting === work) this.#reconnecting = undefined; });
    this.#reconnecting = work; return work;
  }
  async create(input: ModelCreateInput, options: RendererModelCallOptions = {}): Promise<{ model: RendererModelReference; snapshot: ModelSnapshot }> {
    this.#ready(); this.#timeout(options);
    if (this.#models.size + this.#reservations >= 4) throw new RendererModelError(Code.ResourceExhausted, 'renderer.model_limit');
    let json: JsonObject;
    try { json = inputObject(toJson(CreateModelRequestSchema, create(CreateModelRequestSchema, input))); delete json.wireContract; }
    catch { throw new RendererModelError(Code.InvalidArgument, 'renderer.invalid_input'); }
    ++this.#reservations;
    try {
      return await this.#invoke(id => this.bridge.create(id, json), value => {
        if (!record(value) || Object.keys(value).length !== 2 || typeof value.model !== 'string' ||
            !/^[a-f0-9-]{36}$/u.test(value.model) || !Object.hasOwn(value, 'snapshot')) invalidReply();
        const snapshot = readModelSnapshot(value.snapshot as JsonValue);
        const model = Object.freeze({}) as RendererModelReference;
        this.#models.set(model, { token: value.model, snapshot: clone(ModelSnapshotSchema, snapshot) });
        return { model, snapshot };
      }, options, true);
    } finally { --this.#reservations; }
  }
  async describe(model: RendererModelReference, options: RendererModelCallOptions = {}): Promise<ModelSnapshot> {
    this.#ready(); const owned = this.#owned(model);
    return this.#invoke(id => this.bridge.describe(id, owned.token), value => {
      const snapshot = readModelSnapshot(value);
      if (!equals(ModelSnapshotSchema, snapshot, owned.snapshot)) invalidReply();
      return snapshot;
    }, options);
  }
  async solve(model: RendererModelReference, input: ModelSolveInput, options: RendererModelCallOptions = {}): Promise<RendererModelResult> {
    this.#ready(); const owned = this.#owned(model);
    let state; let json: JsonObject;
    try {
      state = clone(SolveModelRequestSchema, create(SolveModelRequestSchema, input)); json = inputObject(toJson(SolveModelRequestSchema, state));
      delete json.wireContract; delete json.modelHandle;
    } catch { throw new RendererModelError(Code.InvalidArgument, 'renderer.invalid_input'); }
    return this.#invoke(id => this.bridge.solve(id, owned.token, json), value => readModelResult(value, owned.snapshot, state), options);
  }
  async release(model: RendererModelReference, options: RendererModelCallOptions = {}): Promise<void> {
    this.#ready(); this.#timeout(options); const owned = this.#owned(model);
    // Only invalidate after local admission; a local limit error must retain a releasable ref.
    const releasing = this.#invoke(id => this.bridge.release(id, owned.token), value => this.#empty(value), options, true);
    this.#models.delete(model); await releasing;
  }
}
export function rendererModelClient(bridge = typeof window === 'undefined' ? undefined : window.mpmcModelDesktop): RendererModelClient | null {
  return bridge ? RendererModelClient.fromBridge(bridge) : null;
}
