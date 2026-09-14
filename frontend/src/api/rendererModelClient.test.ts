import { clone, fromJson, toJson, type DescMessage, type JsonObject, type JsonValue } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';
import { expect, it, vi } from 'vitest';
import { MODEL_DESKTOP_CONVENTION, MODEL_DESKTOP_V2_CONVENTION, type ModelDesktopBridge, type ModelDesktopReply } from './modelDesktopContract';
import { rendererModelClient, RendererModelError, type RendererModelReference } from './rendererModelClient';
import { MODEL_VALIDATION_DETAIL_VERSION } from './modelValidationDetail';
import { readModelReply, readModelResult, readModelSnapshot } from './rendererModelWire';
import {
  FullPtResultSchema, ModelSnapshotSchema, PtEosRootSettingsSchema, PtStabilitySettingsSchema,
  PtTwoPhaseSettingsSchema, PtThreePhaseSettingsSchema, ModelConfigurationLimitsSchema, PtSolverSafetyLimitsSchema,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { PT_BACKEND_CAPABILITY_CONVENTION, PT_BACKEND_RESULT_CONVENTION, PT_PHASE_SET_CONVENTION, PT_PHASE_TRANSITION_CONVENTION } from '../domain/flash';

function block(schema: DescMessage): JsonObject {
  return Object.fromEntries(schema.fields.map(field => [field.jsonName, field.localName === 'automaticMultistart' ? false : 1]));
}
// Serialization-only synthetic messages. Actual scientific fixtures are supplied by the native smoke.
function fixtures() {
  const capability: JsonObject = { convention: PT_BACKEND_CAPABILITY_CONVENTION, datasetId: 'mapping-test', revision: 'r1',
    publicationProfile: 'mapping-publication/v1', componentIds: ['synthetic-a'], supportedPhaseCounts: [1, 2, 3],
    transitionCapability: { convention: PT_PHASE_TRANSITION_CONVENTION }, globalStabilityProven: false };
  const snapshot: JsonObject = {
    definition: { version: 'thermodynamic-model/parameter-definition/v1', family: 'MODEL_FAMILY_PR76',
      datasetId: 'mapping-test', revision: 'r1', components: [{ componentId: 'synthetic-a' }] },
    capability,
    settings: { version: 'pt-solver-settings/v1', kind: 'SOLVER_SETTINGS_KIND_CUSTOM', presetId: '',
      eosRoot: { ...block(PtEosRootSettingsSchema), maxIterations: '9007199254740993' },
      initialStability: block(PtStabilitySettingsSchema), twoPhase: block(PtTwoPhaseSettingsSchema),
      finalTwoPhaseStability: block(PtStabilitySettingsSchema), threePhase: block(PtThreePhaseSettingsSchema),
      finalThreePhaseStability: block(PtStabilitySettingsSchema) },
    parameterLimits: block(ModelConfigurationLimitsSchema), solverLimits: block(PtSolverSafetyLimitsSchema),
  };
  const result: JsonObject = {
    backendResultConvention: PT_BACKEND_RESULT_CONVENTION, phaseSetConvention: PT_PHASE_SET_CONVENTION,
    capability, outcome: 'PT_COMPUTATION_OUTCOME_ACCEPTED', maximumPhaseCount: 3,
    pressurePa: 100000, temperatureK: 300, feed: [1],
    candidatePhaseSet: { phases: [{ molePhaseFraction: 1, composition: [1], lnFugacityCoefficient: [0.25],
      providerBranch: '18446744073709551615', providerBranchSmooth: false, compressibilityFactor: 0.5 }] },
    globalStabilityProven: false, diagnostic: 'Synthetic mapping only.',
    transitionReport: { convention: PT_PHASE_TRANSITION_CONVENTION, evidence: [{ sourcePhaseCount: 2,
      targetPhaseCount: 1, trigger: 'PT_PHASE_TRANSITION_TRIGGER_PHASE_DISAPPEARANCE',
      resolution: 'PT_PHASE_TRANSITION_RESOLUTION_ACCEPTED_TARGET', freshTargetSolveAttempted: true,
      targetTopologyClosed: true, providerEvidenceProfile: 'test', diagnostic: 'preserved evidence' }] },
    providerResultConvention: 'mapping-publication/v1', phaseMetadata: [{ roleId: 'test-role', familyId: 'test-family' }], morphologyResolved: false,
  };
  return { snapshot, result, state: { pressurePa: 100000, temperatureK: 300, feed: [1] } };
}
function success(value: JsonValue): ModelDesktopReply { return { version: MODEL_DESKTOP_CONVENTION, ok: true, value }; }
function failed(code: Code, reason = 'rpc.failed'): ModelDesktopReply { return { version: MODEL_DESKTOP_CONVENTION, ok: false, error: { code, reason } }; }
function setup() {
  const data = fixtures(); let serial = 0;
  const bridge = {
    convention: MODEL_DESKTOP_CONVENTION,
    connect: vi.fn(async (_id: string) => success(null)), reconnect: vi.fn(async (_id: string) => success(null)),
    create: vi.fn(async (_id: string, _input: JsonObject) => success({ model: `00000000-0000-4000-8000-${String(++serial).padStart(12, '0')}`, snapshot: structuredClone(data.snapshot) })),
    describe: vi.fn(async (_id: string, _model: string) => success(structuredClone(data.snapshot))),
    solve: vi.fn(async (_id: string, _model: string, _input: JsonObject) => success(structuredClone(data.result))),
    release: vi.fn(async (_id: string, _model: string) => success(null)), cancel: vi.fn((_id: string) => {}),
  } satisfies ModelDesktopBridge;
  return { ...data, bridge, client: rendererModelClient(bridge)! };
}
function deferred<T>() {
  let resolve!: (value: T) => void;
  const promise = new Promise<T>(yes => { resolve = yes; }); return { promise, resolve };
}

it('caches one client per bridge and rejects absent/unsupported preload contracts', () => {
  const { bridge, client } = setup(); expect(rendererModelClient(bridge)).toBe(client);
  expect(rendererModelClient()).toBeNull();
  expect(() => rendererModelClient({ ...bridge, convention: 'v2' } as unknown as ModelDesktopBridge)).toThrow(RendererModelError);
  expect(() => rendererModelClient({ convention: MODEL_DESKTOP_CONVENTION } as ModelDesktopBridge)).toThrow(RendererModelError);
});

it('roundtrips full snapshots/results, exact bigints, metadata and private immutable reference identity', async () => {
  const { bridge, client, snapshot, result, state } = setup();
  const made = await client.create({ solverSelection: { case: 'presetId', value: 'explicit-test' } });
  expect(bridge.create.mock.calls[0]?.[1]).toEqual({ presetId: 'explicit-test' });
  expect(Object.keys(made.model)).toEqual([]); expect(Object.isFrozen(made.model)).toBe(true);
  expect(made.snapshot.settings?.eosRoot?.maxIterations).toBe(9007199254740993n);
  expect(made.snapshot).toEqual(fromJson(ModelSnapshotSchema, snapshot));
  made.snapshot.settings!.eosRoot!.maxIterations = 2n; // Private expected snapshot is an independent clone.
  expect(await client.describe(made.model)).toEqual(fromJson(ModelSnapshotSchema, snapshot));
  const solved = await client.solve(made.model, state);
  expect(solved.outcome).toBe('accepted'); expect(solved.result).toEqual(fromJson(FullPtResultSchema, result));
  expect(solved.result.candidatePhaseSet?.phases[0]?.providerBranch).toBe(18446744073709551615n);
  expect(bridge.solve.mock.calls[0]?.[2]).toEqual(state);
  await client.release(made.model);
  await expect(client.describe(made.model)).rejects.toMatchObject({ code: Code.NotFound, category: 'not_found' });
  expect(bridge.describe).toHaveBeenCalledTimes(1);
});

it('maps all computation outcomes without promoting or dropping diagnostic candidates', () => {
  const { snapshot, result, state } = fixtures(); const expected = readModelSnapshot(snapshot);
  for (const [wire, outcome] of [['PT_COMPUTATION_OUTCOME_ACCEPTED', 'accepted'], ['PT_COMPUTATION_OUTCOME_PHASE_SET_UNSTABLE', 'phase_set_unstable'], ['PT_COMPUTATION_OUTCOME_INDETERMINATE', 'indeterminate']]) {
    const json = structuredClone(result); json.outcome = wire!;
    if (outcome !== 'accepted') {
      ((json.candidatePhaseSet as JsonObject).phases as JsonObject[])[0]!.lnFugacityCoefficient = ['NaN'];
      ((json.candidatePhaseSet as JsonObject).phases as JsonObject[])[0]!.compressibilityFactor = 'Infinity';
    }
    const mapped = readModelResult(json, expected, state);
    expect(mapped.outcome).toBe(outcome); expect(mapped.result).toEqual(fromJson(FullPtResultSchema, json));
    expect(toJson(FullPtResultSchema, mapped.result)).toEqual(toJson(FullPtResultSchema, fromJson(FullPtResultSchema, json)));
    if (outcome !== 'accepted') expect(mapped.result.candidatePhaseSet!.phases[0]!.lnFugacityCoefficient[0]).toBeNaN();
  }
  const noCandidate: JsonObject = { ...result, outcome: 'PT_COMPUTATION_OUTCOME_INDETERMINATE' }; delete noCandidate.candidatePhaseSet;
  expect(readModelResult(noCandidate, expected, state).result.candidatePhaseSet).toBeUndefined();
});

it('maps all 16 gRPC failures distinctly and sanitizes unknown reasons/exception text', async () => {
  const categories = ['cancelled', 'unknown', 'invalid_argument', 'deadline_exceeded', 'not_found', 'already_exists',
    'permission_denied', 'resource_exhausted', 'failed_precondition', 'aborted', 'out_of_range', 'unimplemented',
    'internal', 'unavailable', 'data_loss', 'unauthenticated'];
  for (let code = 1; code <= 16; ++code) {
    const { bridge, client, state } = setup(); const made = await client.create({});
    bridge.solve.mockResolvedValueOnce(failed(code as Code));
    await expect(client.solve(made.model, state)).rejects.toMatchObject({ code, category: categories[code - 1], source: 'ipc', reason: 'rpc.failed' });
  }
  expect(() => readModelReply(failed(Code.InvalidArgument, 'private-model-data'))).toThrow('ipc.failed');
  const { bridge, client } = setup(); bridge.connect.mockRejectedValueOnce(new Error('private-host-token'));
  await expect(client.connect()).rejects.toMatchObject({ code: Code.Unavailable, reason: 'renderer.invoke_failed', source: 'transport' });
  expect(client.requiresReconnect).toBe(true);
});

it('rejects malformed reply versions/shapes/error codes and bounds payload size', () => {
  for (const raw of [null, {}, { ...success(null), version: 'v2' }, { ...success(null), value: undefined },
    { ...success(null), extra: true }, { ...success(null), ok: 1 }, failed(0 as Code), failed(17 as Code),
    failed(1.5 as Code), { ...failed(Code.NotFound), error: { code: 5, reason: null } },
    success('x'.repeat(4 * 1024 * 1024))]) {
    expect(() => readModelReply(raw)).toThrow(RendererModelError);
  }
});

it('rejects incomplete/version-mismatched snapshots and misaligned result identity', () => {
  const { snapshot, result, state } = fixtures();
  const badSettings = structuredClone(snapshot);
  delete ((badSettings.settings as JsonObject).initialStability as JsonObject).automaticMultistart;
  for (const json of [{}, { ...snapshot, definition: {} }, { ...snapshot, capability: {} }, badSettings,
    { ...snapshot, settings: { ...(snapshot.settings as JsonObject), version: 'v2' } }, { ...snapshot, unexpected: true }]) {
    expect(() => readModelSnapshot(json)).toThrow(RendererModelError);
  }
  const expected = readModelSnapshot(snapshot);
  for (const compressibilityFactor of [0, -1, 'Infinity', 'NaN']) {
    const malformed = structuredClone(result);
    ((malformed.candidatePhaseSet as JsonObject).phases as JsonObject[])[0]!.compressibilityFactor = compressibilityFactor;
    expect(() => readModelResult(malformed, expected, state)).toThrow(RendererModelError);
  }
  expect(() => readModelResult({ ...result, pressurePa: -1 }, expected, { ...state, pressurePa: -1 })).toThrow(RendererModelError);
  expect(() => readModelResult({ ...result, feed: [-1] }, expected, { ...state, feed: [-1] })).toThrow(RendererModelError);
  for (const json of [{}, { ...result, outcome: 999 }, { ...result, outcome: 'PT_COMPUTATION_OUTCOME_UNSPECIFIED' },
    { ...result, backendResultConvention: 'v2' }, { ...result, pressurePa: 101000 }, { ...result, feed: [0.5] },
    { ...result, capability: { ...(result.capability as JsonObject), revision: 'r2' } },
    { ...result, candidatePhaseSet: { phases: [] } }, { ...result, unexpected: true }]) {
    expect(() => readModelResult(json, expected, state)).toThrow(RendererModelError);
  }
});

it('rejects changed describe snapshots and isolates caller mutations of solve inputs', async () => {
  const { bridge, client, state, snapshot, result } = setup(); const made = await client.create({});
  const late = deferred<ModelDesktopReply>(); bridge.solve.mockReturnValueOnce(late.promise);
  const solving = client.solve(made.model, state); state.feed[0] = 0.5;
  late.resolve(success(result)); expect((await solving).outcome).toBe('accepted');
  const changed = clone(ModelSnapshotSchema, fromJson(ModelSnapshotSchema, snapshot)); changed.settings!.eosRoot!.maxIterations = 7n;
  bridge.describe.mockResolvedValueOnce(success(toJson(ModelSnapshotSchema, changed)));
  await expect(client.describe(made.model)).rejects.toMatchObject({ source: 'contract', code: Code.DataLoss });
  expect(client.requiresReconnect).toBe(true);
});

it('rejects foreign, released and prior-session references before IPC dispatch', async () => {
  const a = setup(); const b = setup(); const made = await a.client.create({});
  await expect(b.client.describe(made.model)).rejects.toMatchObject({ reason: 'renderer.stale_reference' });
  await expect(a.client.describe({} as RendererModelReference)).rejects.toMatchObject({ code: Code.NotFound });
  const reconnecting = a.client.reconnect(); expect(a.client.reconnect()).toBe(reconnecting); await reconnecting;
  await expect(a.client.solve(made.model, a.state)).rejects.toMatchObject({ code: Code.NotFound });
  expect(a.bridge.solve).not.toHaveBeenCalled(); expect(b.bridge.describe).not.toHaveBeenCalled();
});

it('enforces model reservations and restores capacity after explicit release', async () => {
  const { client, bridge } = setup();
  const made = await Promise.all(Array.from({ length: 4 }, () => client.create({})));
  await expect(client.create({})).rejects.toMatchObject({ reason: 'renderer.model_limit' });
  await client.release(made[0]!.model); await client.create({});
  expect(bridge.create).toHaveBeenCalledTimes(5);
  expect(new Set(bridge.create.mock.calls.map(([id]) => id)).size).toBe(5);
});

it('keeps pending admission bounded after cancel, drains before reconnect and discards late success', async () => {
  const { client, bridge, state, result } = setup(); const made = await client.create({});
  const late = deferred<ModelDesktopReply>(); bridge.solve.mockReturnValue(late.promise);
  const calls = Array.from({ length: 4 }, () => client.solve(made.model, state).catch(error => error));
  await expect(client.describe(made.model)).rejects.toMatchObject({ reason: 'renderer.request_limit' });
  const reconnecting = client.reconnect(); expect(bridge.reconnect).not.toHaveBeenCalled();
  expect(bridge.cancel).toHaveBeenCalledTimes(4);
  for (const error of await Promise.all(calls)) expect(error).toMatchObject({ code: Code.Canceled });
  late.resolve(success(result)); await reconnecting;
  expect(bridge.reconnect).toHaveBeenCalledTimes(1); expect(client.requiresReconnect).toBe(false);
});

it('suppresses late creates after AbortSignal cancellation and never retries ambiguous mutations', async () => {
  const { client, bridge, snapshot } = setup(); const late = deferred<ModelDesktopReply>();
  bridge.create.mockReturnValueOnce(late.promise); const signal = new AbortController();
  const making = client.create({}, { signal: signal.signal }).catch(error => error); signal.abort();
  expect(await making).toMatchObject({ code: Code.Canceled, reason: 'renderer.cancelled' });
  late.resolve(success({ model: '00000000-0000-4000-8000-000000000001', snapshot }));
  await expect(client.create({})).rejects.toMatchObject({ reason: 'renderer.reconnect_required' });
  expect(bridge.create).toHaveBeenCalledTimes(1);
  await client.reconnect(); await client.create({}); expect(bridge.create).toHaveBeenCalledTimes(2);
});

it('handles pre-aborted calls and local deadlines without listener or reconnect queues', async () => {
  vi.useFakeTimers();
  try {
    const { client, bridge } = setup(); const aborted = new AbortController(); aborted.abort();
    await expect(client.connect({ signal: aborted.signal })).rejects.toMatchObject({ code: Code.Canceled });
    await expect(client.connect({ timeoutMs: 0 })).rejects.toMatchObject({ code: Code.InvalidArgument });
    expect(bridge.connect).not.toHaveBeenCalled();
    const late = deferred<ModelDesktopReply>(); bridge.connect.mockReturnValueOnce(late.promise);
    const signal = new AbortController(); const remove = vi.spyOn(signal.signal, 'removeEventListener');
    const pending = client.connect({ timeoutMs: 10, signal: signal.signal }).catch(error => error);
    await vi.advanceTimersByTimeAsync(10);
    expect(await pending).toMatchObject({ reason: 'renderer.timeout', category: 'deadline_exceeded' });
    expect(remove).toHaveBeenCalled();
    const reconnecting = client.reconnect({ timeoutMs: 10, signal: signal.signal }).catch(error => error);
    await vi.advanceTimersByTimeAsync(10);
    expect(await reconnecting).toMatchObject({ reason: 'renderer.drain_timeout' });
    await expect(client.reconnect()).rejects.toMatchObject({ reason: 'renderer.calls_draining' });
    expect(vi.getTimerCount()).toBe(0);
    late.resolve(success(null)); await Promise.resolve(); await Promise.resolve(); await Promise.resolve();
    await client.reconnect(); expect(client.requiresReconnect).toBe(false);
  } finally { vi.useRealTimers(); }
});

it('prefers v2, preserves immutable validation details and never downgrades failed calls', async () => {
  const { bridge } = setup();
  const validation = { version: MODEL_VALIDATION_DETAIL_VERSION, code: 'configuration.invalid_value', field: 'parameters.pure[nitrogen].critical_temperature_k' };
  const reply: ModelDesktopReply = { version: MODEL_DESKTOP_V2_CONVENTION, ok: false, error: { code: Code.InvalidArgument, reason: 'rpc.failed', validation } };
  const v2 = { ...bridge, convention: MODEL_DESKTOP_V2_CONVENTION, connect: vi.fn(async () => reply) };
  vi.stubGlobal('window', { mpmcModelDesktop: bridge, mpmcModelDesktopV2: v2 });
  try {
    const client = rendererModelClient()!;
    expect(client).toBe(rendererModelClient(v2));
    const error = await client.connect().catch(cause => cause);
    expect(error).toBeInstanceOf(RendererModelError); expect(error.validation).toEqual(validation);
    expect(Object.isFrozen(error.validation)).toBe(true); validation.field = 'changed';
    expect(error.validation.field).toBe('parameters.pure[nitrogen].critical_temperature_k');
    expect(bridge.connect).not.toHaveBeenCalled(); expect(v2.connect).toHaveBeenCalledTimes(1);
    vi.stubGlobal('window', { mpmcModelDesktop: bridge });
    expect(rendererModelClient()).toBe(rendererModelClient(bridge));
  } finally { vi.unstubAllGlobals(); }
});
it('enforces v1/v2 envelopes and safely ignores incompatible optional validation details', () => {
  const detail = { version: MODEL_VALIDATION_DETAIL_VERSION, code: 'configuration.invalid_value', field: 'eos_root.max_iterations' };
  const reply = { version: MODEL_DESKTOP_V2_CONVENTION, ok: false, error: { code: Code.InvalidArgument, reason: 'rpc.failed', validation: detail } };
  expect(() => readModelReply(reply)).toThrow('renderer.invalid_reply');
  expect(() => readModelReply({ ...reply, version: MODEL_DESKTOP_CONVENTION })).toThrow('renderer.invalid_reply');
  expect(() => readModelReply(failed(Code.InvalidArgument), MODEL_DESKTOP_V2_CONVENTION)).toThrow('renderer.invalid_reply');
  for (const validation of [{ ...detail, version: 'future/v2' }, { ...detail, code: 'rpc.internal_failure' }, { ...detail, message: 'private' }, null]) {
    try { readModelReply({ ...reply, error: { ...reply.error, validation } }, MODEL_DESKTOP_V2_CONVENTION); throw new Error('Expected error'); }
    catch (error) { expect(error).toMatchObject({ code: Code.InvalidArgument, reason: 'rpc.failed', validation: undefined }); }
  }
});
