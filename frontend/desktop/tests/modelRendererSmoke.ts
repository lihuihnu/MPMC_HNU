// Browser-only test entry: executed inside the same sandboxed renderer as production.
import { clone, fromJson, toJson, type JsonObject } from '@bufbuild/protobuf';
import { rendererModelClient, RendererModelError, type RendererModelReference } from '../../src/api/rendererModelClient';
import { CreateModelRequestSchema, FullPtResultSchema, ModelSnapshotSchema, SolveModelRequestSchema, PtSolverSettingsSchema, SolverSettingsKind, type PtSolverSettings, type CreateModelRequest, type SolveModelRequest } from '../../src/gen/mpmc/model_configuration/v1/model_service_pb';
import { MODEL_DESKTOP_V2_CONVENTION } from '../../src/api/modelDesktopContract';
const client = rendererModelClient();
let current: RendererModelReference | undefined;
let definition: CreateModelRequest;
let state: SolveModelRequest;
let settings: PtSolverSettings;
function error(cause: unknown) {
  if (!(cause instanceof RendererModelError)) throw new Error('Untyped renderer model error.');
  return { code: cause.code, category: cause.category, reason: cause.reason, source: cause.source, ...(cause.validation ? { validation: cause.validation } : {}) };
}
export async function initialize(createInput: JsonObject, solveInput: JsonObject) {
  if (!client || window.mpmcModelDesktopV2?.convention !== MODEL_DESKTOP_V2_CONVENTION) throw new Error('Model v2 preload unavailable in renderer.');
  definition = fromJson(CreateModelRequestSchema, createInput); state = fromJson(SolveModelRequestSchema, solveInput);
  const made = await client.create(definition); current = made.model;
  settings = clone(PtSolverSettingsSchema, made.snapshot.settings!);
  const described = await client.describe(current);
  const solved = await client.solve(current, state);
  return { snapshot: toJson(ModelSnapshotSchema, made.snapshot), described: toJson(ModelSnapshotSchema, described),
    outcome: solved.outcome, result: toJson(FullPtResultSchema, solved.result) };
}
export async function solve() {
  if (!client || !current) throw new Error('Renderer model missing.');
  const solved = await client.solve(current, state);
  return { outcome: solved.outcome, result: toJson(FullPtResultSchema, solved.result) };
}
export async function invalidSolve() {
  if (!client || !current) throw new Error('Renderer model missing.');
  try { await client.solve(current, { ...state, pressurePa: -1 }); }
  catch (cause) { return error(cause); }
  throw new Error('Invalid solve unexpectedly succeeded.');
}
export async function releaseAndRecreate() {
  if (!client || !current) throw new Error('Renderer model missing.');
  const old = current; await client.release(old);
  let released;
  try { await client.describe(old); } catch (cause) { released = error(cause); }
  current = (await client.create(definition)).model;
  return { released, ...(await solve()) };
}
export async function invalidCreateAndReconnect() {
  if (!client || !current) throw new Error('Renderer model missing.');
  const old = current; let failure;
  try { await client.create({ ...definition, solverSelection: { case: 'presetId', value: 'unknown-renderer-test-preset' } }); }
  catch (cause) { failure = error(cause); }
  if (!client.requiresReconnect) throw new Error('Ambiguous mutation did not invalidate renderer ownership.');
  let blocked;
  try { await client.solve(old, state); } catch (cause) { blocked = error(cause); }
  await client.reconnect(); let stale;
  try { await client.describe(old); } catch (cause) { stale = error(cause); }
  current = (await client.create(definition)).model;
  return { failure, blocked, stale, ...(await solve()) };
}

export async function invalidFieldsAndRecover() {
  if (!client || !current) throw new Error('Renderer model missing.');
  const bad = clone(CreateModelRequestSchema, definition);
  if (bad.definition?.parameters.case !== 'pr76') throw new Error('Expected attributed PR76 fixture.');
  const pure = bad.definition.parameters.value.pure[0]!;
  const expectedField = `parameters.pure[${pure.componentId}].critical_temperature_k`;
  pure.criticalTemperatureK = undefined;
  let parameter;
  try { await client.create(bad); } catch (cause) { parameter = error(cause); }
  await client.reconnect(); current = (await client.create(definition)).model;
  const custom = clone(PtSolverSettingsSchema, settings);
  custom.kind = SolverSettingsKind.CUSTOM; custom.presetId = ''; custom.eosRoot!.maxIterations = undefined;
  let setting;
  try { await client.create({ ...definition, solverSelection: { case: 'settings', value: custom } }); }
  catch (cause) { setting = error(cause); }
  await client.reconnect(); current = (await client.create(definition)).model;
  return { parameter, setting, expectedField, ...(await solve()) };
}

export async function invalidNestedFieldsAndRecover() {
  if (!client || !current) throw new Error('Renderer model missing.');
  const failures = [];
  for (const location of ['pureValue', 'binarySource', 'componentSource']) {
    const bad = clone(CreateModelRequestSchema, definition);
    if (bad.definition?.parameters.case !== 'pr76') throw new Error('Expected attributed PR76 fixture.');
    const parameters = bad.definition.parameters.value;
    if (location === 'pureValue') parameters.pure[1]!.criticalTemperatureK!.value = undefined;
    else if (location === 'binarySource') parameters.binary[0]!.kij!.provenance!.kind = undefined;
    else bad.definition.components[1]!.provenance!.kind = undefined;
    let failure;
    try { await client.create(bad); } catch (cause) { failure = error(cause); }
    if (!failure || !client.requiresReconnect) throw new Error('Nested invalid create was not rejected.');
    failures.push(failure);
    await client.reconnect(); current = (await client.create(definition)).model;
  }
  return { failures, ...(await solve()) };
}
