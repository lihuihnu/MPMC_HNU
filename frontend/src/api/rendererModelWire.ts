import { equals, fromJson, type DescMessage, type JsonValue, type Message } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';
import { MODEL_DESKTOP_CONVENTION, MODEL_DESKTOP_V2_CONVENTION, type ModelDesktopVersion } from './modelDesktopContract';
import { readModelValidationDetail, type ModelValidationDetail } from './modelValidationDetail';
import {
  FullPtResultSchema, ModelSnapshotSchema, SolverSettingsKind,
  PtEosRootSettingsSchema, PtStabilitySettingsSchema, PtTwoPhaseSettingsSchema,
  PtThreePhaseSettingsSchema, ModelConfigurationLimitsSchema, PtSolverSafetyLimitsSchema,
  type FullPtResult, type ModelSnapshot,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { PtBackendCapabilitySchema, PtComputationOutcome as WireOutcome } from '../gen/mpmc/runtime/v1/pt_service_pb';
import { PT_BACKEND_CAPABILITY_CONVENTION, PT_BACKEND_RESULT_CONVENTION,
  PT_PHASE_SET_CONVENTION, PT_PHASE_TRANSITION_CONVENTION, type PtComputationOutcome } from '../domain/flash';

const categories = {
  [Code.Canceled]: 'cancelled', [Code.Unknown]: 'unknown', [Code.InvalidArgument]: 'invalid_argument',
  [Code.DeadlineExceeded]: 'deadline_exceeded', [Code.NotFound]: 'not_found',
  [Code.AlreadyExists]: 'already_exists', [Code.PermissionDenied]: 'permission_denied',
  [Code.ResourceExhausted]: 'resource_exhausted', [Code.FailedPrecondition]: 'failed_precondition',
  [Code.Aborted]: 'aborted', [Code.OutOfRange]: 'out_of_range', [Code.Unimplemented]: 'unimplemented',
  [Code.Internal]: 'internal', [Code.Unavailable]: 'unavailable', [Code.DataLoss]: 'data_loss',
  [Code.Unauthenticated]: 'unauthenticated',
} as const;
export type RendererModelErrorCategory = typeof categories[Code];
export class RendererModelError extends Error {
  readonly category: RendererModelErrorCategory;
  constructor(readonly code: Code, readonly reason: string,
    readonly source: 'client' | 'ipc' | 'transport' | 'contract' = 'client', readonly validation?: ModelValidationDetail) {
    super(`Model request: ${reason}`); this.name = 'RendererModelError'; this.category = categories[code];
  }
}
export function invalidReply(): never {
  throw new RendererModelError(Code.DataLoss, 'renderer.invalid_reply', 'contract');
}
export function record(value: unknown): value is Record<string, unknown> {
  return value !== null && typeof value === 'object' && !Array.isArray(value);
}
function fields(value: Record<string, unknown>, keys: readonly string[]): boolean {
  return Object.keys(value).length === keys.length && keys.every(key => Object.hasOwn(value, key));
}
const reasons = new Set([
  'ipc.invalid_request', 'ipc.request_size', 'ipc.unsupported_version', 'ipc.unknown_operation',
  'ipc.invalid_reference', 'ipc.invalid_input', 'ipc.window_closed', 'ipc.duplicate_request',
  'ipc.request_limit', 'ipc.document_changed', 'ipc.stale_reference', 'ipc.reply_size',
  'ipc.failed', 'ipc.sender_rejected', 'client.disposed', 'client.request_limit',
  'client.invalid_timeout', 'client.model_limit', 'session.closing', 'session.open_timeout',
  'session.invalid_handshake', 'session.open_failed', 'session.ended', 'session.closed',
  'session.reconnect_cancelled', 'session.not_connected', 'session.changed',
  'model.stale_reference', 'model.invalid_response', 'rpc.failed',
]);

/** Treat even a typed preload promise as untrusted data at runtime. */
export function readModelReply(raw: unknown, version: ModelDesktopVersion = MODEL_DESKTOP_CONVENTION): JsonValue {
  try {
    const text = JSON.stringify(raw);
    if (!text || new TextEncoder().encode(text).length > 4 * 1024 * 1024) invalidReply();
  } catch { invalidReply(); }
  if (!record(raw) || raw.version !== version) invalidReply();
  if (raw.ok === false && fields(raw, ['version', 'ok', 'error']) && record(raw.error) &&
      (fields(raw.error, ['code', 'reason']) || (version === MODEL_DESKTOP_V2_CONVENTION && fields(raw.error, ['code', 'reason', 'validation'])))) {
    const { code, reason } = raw.error;
    if (typeof code !== 'number' || !Number.isInteger(code) || code < 1 || code > 16 ||
        typeof reason !== 'string' || !reason || reason.length > 128) invalidReply();
    const validation = version === MODEL_DESKTOP_V2_CONVENTION ? readModelValidationDetail(raw.error.validation, code as Code) : undefined;
    throw new RendererModelError(code as Code, reasons.has(reason) ? reason : 'ipc.failed', 'ipc', validation);
  }
  if (raw.ok !== true || !fields(raw, ['version', 'ok', 'value']) || raw.value === undefined) invalidReply();
  return raw.value as JsonValue;
}
function complete(schema: DescMessage, value: Message | undefined): void {
  if (!value || !schema.fields.every(field => Object.hasOwn(value, field.localName))) invalidReply();
}
export function readModelSnapshot(value: JsonValue): ModelSnapshot {
  let snapshot: ModelSnapshot;
  try { snapshot = fromJson(ModelSnapshotSchema, value); } catch { invalidReply(); }
  const { definition, settings, capability } = snapshot;
  if (!definition || definition.version !== 'thermodynamic-model/parameter-definition/v1' ||
      !settings || settings.version !== 'pt-solver-settings/v1' || !capability ||
      capability.convention !== PT_BACKEND_CAPABILITY_CONVENTION ||
      capability.transitionCapability?.convention !== PT_PHASE_TRANSITION_CONVENTION ||
      !definition.datasetId || !definition.revision || capability.datasetId !== definition.datasetId ||
      capability.revision !== definition.revision || !capability.publicationProfile) invalidReply();
  if (settings.kind !== SolverSettingsKind.PRESET && settings.kind !== SolverSettingsKind.CUSTOM) invalidReply();
  if (settings.kind === SolverSettingsKind.PRESET ? !settings.presetId : !!settings.presetId) invalidReply();
  const ids = definition.components.map(component => component.componentId);
  if (!ids.length || ids.length > 256 || ids.some(id => !id) || new Set(ids).size !== ids.length ||
      ids.length !== capability.componentIds.length || ids.some((id, index) => id !== capability.componentIds[index]) ||
      !capability.supportedPhaseCounts.length || capability.supportedPhaseCounts.some(n => n < 1 || n > 64)) invalidReply();
  complete(PtEosRootSettingsSchema, settings.eosRoot);
  complete(PtStabilitySettingsSchema, settings.initialStability);
  complete(PtTwoPhaseSettingsSchema, settings.twoPhase);
  complete(PtStabilitySettingsSchema, settings.finalTwoPhaseStability);
  complete(PtThreePhaseSettingsSchema, settings.threePhase);
  complete(PtStabilitySettingsSchema, settings.finalThreePhaseStability);
  complete(ModelConfigurationLimitsSchema, snapshot.parameterLimits);
  complete(PtSolverSafetyLimitsSchema, snapshot.solverLimits);
  return snapshot;
}

export interface RendererModelResult {
  readonly outcome: PtComputationOutcome;
  /** Complete native message, including diagnostic candidates for non-accepted outcomes. */
  readonly result: FullPtResult;
}
export interface EvaluatedModelState { pressurePa?: number | undefined; temperatureK?: number | undefined; feed: readonly number[] }
export function readModelResult(value: JsonValue, snapshot: ModelSnapshot, state: EvaluatedModelState): RendererModelResult {
  let result: FullPtResult;
  try { result = fromJson(FullPtResultSchema, value); } catch { invalidReply(); }
  let outcome: PtComputationOutcome;
  switch (result.outcome) {
    case WireOutcome.ACCEPTED: outcome = 'accepted'; break;
    case WireOutcome.PHASE_SET_UNSTABLE: outcome = 'phase_set_unstable'; break;
    case WireOutcome.INDETERMINATE: outcome = 'indeterminate'; break;
    default: invalidReply();
  }
  if (result.backendResultConvention !== PT_BACKEND_RESULT_CONVENTION || result.phaseSetConvention !== PT_PHASE_SET_CONVENTION ||
      !result.capability || !snapshot.capability || !equals(PtBackendCapabilitySchema, result.capability, snapshot.capability) ||
      result.transitionReport?.convention !== PT_PHASE_TRANSITION_CONVENTION ||
      result.providerResultConvention !== snapshot.capability.publicationProfile ||
      result.maximumPhaseCount !== Math.max(...snapshot.capability.supportedPhaseCounts) ||
      result.pressurePa !== state.pressurePa || result.temperatureK !== state.temperatureK ||
      !Number.isFinite(result.pressurePa) || !Number.isFinite(result.temperatureK) ||
      result.pressurePa! <= 0 || result.temperatureK! <= 0 ||
      typeof result.globalStabilityProven !== 'boolean' || typeof result.morphologyResolved !== 'boolean' ||
      typeof result.diagnostic !== 'string' || result.feed.length !== state.feed.length ||
      result.feed.length !== snapshot.capability.componentIds.length ||
      (result.globalStabilityProven && !snapshot.capability.globalStabilityProven) ||
      (result.morphologyResolved && !snapshot.capability.phaseMetadataNamespace)) invalidReply();
  for (const [i, actual] of result.feed.entries()) {
    const requested = state.feed[i]!;
    // Same transport identity allowance as ptWire.ts; no normalization or solver tolerance change.
    if (!Number.isFinite(actual) || !Number.isFinite(requested) || actual < 0 || actual > 1 ||
        Math.abs(actual - requested) > 2 * 64 * Number.EPSILON * Math.max(1, Math.abs(actual), Math.abs(requested))) invalidReply();
  }
  const phases = result.candidatePhaseSet?.phases ?? [];
  if (phases.length > result.maximumPhaseCount! ||
      (outcome === 'accepted' && (!phases.length || !snapshot.capability.supportedPhaseCounts.includes(phases.length)))) invalidReply();
  for (const phase of phases) {
    if (phase.composition.length !== result.feed.length || phase.lnFugacityCoefficient.length !== result.feed.length ||
        phase.molePhaseFraction === undefined || phase.providerBranch === undefined || phase.providerBranchSmooth === undefined) invalidReply();
    // Diagnostic candidates may contain non-finite data. Preserve them without promoting acceptance.
    if (outcome === 'accepted' && (![phase.molePhaseFraction, ...phase.composition, ...phase.lnFugacityCoefficient].every(Number.isFinite) ||
        phase.molePhaseFraction < 0 || phase.molePhaseFraction > 1 || phase.composition.some(x => x < 0 || x > 1) ||
        (phase.compressibilityFactor !== undefined && (!Number.isFinite(phase.compressibilityFactor) || phase.compressibilityFactor <= 0)))) invalidReply();
  }
  return { outcome, result };
}
