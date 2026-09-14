// Browser-only test entry: executed inside the same sandboxed renderer as production.
import { fromJson, toJson, type JsonObject } from '@bufbuild/protobuf';
import { rendererModelClient, RendererModelError, type RendererModelReference } from '../../src/api/rendererModelClient';
import { CreateModelRequestSchema, FullPtResultSchema, ModelSnapshotSchema, SolveModelRequestSchema, type CreateModelRequest, type SolveModelRequest } from '../../src/gen/mpmc/model_configuration/v1/model_service_pb';
const client = rendererModelClient();
let current: RendererModelReference | undefined;
let definition: CreateModelRequest;
let state: SolveModelRequest;
function error(cause: unknown) {
  if (!(cause instanceof RendererModelError)) throw new Error('Untyped renderer model error.');
  return { code: cause.code, category: cause.category, reason: cause.reason, source: cause.source };
}
export async function initialize(createInput: JsonObject, solveInput: JsonObject) {
  if (!client) throw new Error('Model preload unavailable in renderer.');
  definition = fromJson(CreateModelRequestSchema, createInput); state = fromJson(SolveModelRequestSchema, solveInput);
  const made = await client.create(definition); current = made.model;
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
