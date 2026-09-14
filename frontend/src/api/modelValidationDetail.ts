import { Code } from '@connectrpc/connect';
import {
  ComponentDefinitionSchema, ModelApplicabilitySchema, ModelProvenanceSchema, ModelScalarSchema,
  Pr76PureParametersSchema, Pr76BinaryInteractionSchema, Pr76ParameterDefinitionSchema,
  ThermodynamicModelDefinitionSchema, PtEosRootSettingsSchema, PtStabilitySettingsSchema,
  PtTwoPhaseSettingsSchema, PtThreePhaseSettingsSchema, PtSolverSettingsSchema,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';

export const MODEL_VALIDATION_DETAIL_VERSION = 'MPMC/model/validation-detail/v1' as const;
export interface ModelValidationDetail {
  readonly version: typeof MODEL_VALIDATION_DETAIL_VERSION;
  readonly code: string;
  /** Native advisory path, not a JSON Pointer. Missing means no safe location was supplied. */
  readonly field?: string;
}
const configurationCodes = {
  unsupported_version: Code.Unimplemented, unsupported_family: Code.Unimplemented,
  unsupported_preset: Code.Unimplemented, resource_limit: Code.ResourceExhausted,
  missing_field: Code.InvalidArgument, invalid_identifier: Code.InvalidArgument,
  duplicate_identifier: Code.InvalidArgument, unknown_component: Code.InvalidArgument,
  invalid_value: Code.InvalidArgument, invalid_unit: Code.InvalidArgument,
  invalid_source: Code.InvalidArgument, invalid_range: Code.InvalidArgument,
  duplicate_parameter: Code.InvalidArgument, missing_parameter: Code.InvalidArgument,
  invalid_pair: Code.InvalidArgument, invalid_settings: Code.InvalidArgument,
};
const statuses = new Map<string, Code>(Object.entries(configurationCodes).map(([code, status]) => [`configuration.${code}`, status]));
statuses.set('wire.missing_field', Code.InvalidArgument);
statuses.set('wire.unknown_field', Code.InvalidArgument);
statuses.set('wire.feed_limit', Code.ResourceExhausted);
statuses.set('request.rejected', Code.InvalidArgument);
const names = new Set([
  'configuration', 'definition', 'parameters', 'settings', 'solver_selection', 'scalar',
  'component_kind', 'source_kind', 'settings_kind', 'temperature_k', 'pressure_pa', 'PT',
  'feed', 'request', 'state', 'id',
  ...[ComponentDefinitionSchema, ModelApplicabilitySchema, ModelProvenanceSchema, ModelScalarSchema,
    Pr76PureParametersSchema, Pr76BinaryInteractionSchema, Pr76ParameterDefinitionSchema,
    ThermodynamicModelDefinitionSchema, PtEosRootSettingsSchema, PtStabilitySettingsSchema,
    PtTwoPhaseSettingsSchema, PtThreePhaseSettingsSchema, PtSolverSettingsSchema]
    .flatMap(schema => schema.fields.map(field => field.name)),
]);
function safeField(field: unknown, privateValues: readonly string[]): field is string {
  if (typeof field !== 'string' || !field || field.length > 256 ||
      /mh1_|ms1_|bearer|authorization|cookie|token|session|handle/iu.test(field) ||
      privateValues.some(value => value.length > 0 && field.includes(value))) return false;
  // Tokenize without splitting dots inside component IDs. Never reinterpret key selectors as indices.
  const segment = /([A-Za-z_][A-Za-z_0-9]*)(?:\[([A-Za-z0-9_.:-]{1,128}(?:,[A-Za-z0-9_.:-]{1,128})?)\])?(?:\.|$)/gyu;
  let offset = 0;
  while (offset < field.length) {
    segment.lastIndex = offset;
    const match = segment.exec(field);
    if (!match || !names.has(match[1]!)) return false;
    if (match[2] !== undefined) {
      if (match[1] === 'feed') {
        if (!/^(0|[1-9][0-9]*)$/u.test(match[2])) return false;
      } else if (!['components', 'pure', 'binary'].includes(match[1]!)) return false;
    }
    offset = segment.lastIndex;
  }
  return !field.endsWith('.');
}
/** Fail closed for optional details while retaining the authoritative RPC status. */
export function readModelValidationDetail(raw: unknown, status: Code, privateValues: readonly string[] = []): ModelValidationDetail | undefined {
  if (raw === null || typeof raw !== 'object' || Array.isArray(raw)) return undefined;
  const value = raw as Record<string, unknown>;
  if (Object.keys(value).some(key => !['version', 'code', 'field'].includes(key)) ||
      value.version !== MODEL_VALIDATION_DETAIL_VERSION || typeof value.code !== 'string' ||
      statuses.get(value.code) !== status) return undefined;
  const field = safeField(value.field, privateValues) ? value.field : undefined;
  return Object.freeze({ version: MODEL_VALIDATION_DETAIL_VERSION, code: value.code, ...(field === undefined ? {} : { field }) });
}
