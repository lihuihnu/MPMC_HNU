import { registerPlugin } from '@capacitor/core';

import type {
  FlashClient,
  PtClientCallOptions,
} from '../../../../frontend/src/api/flashClient';
import {
  PT_BACKEND_CAPABILITY_CONVENTION,
  PT_CAPABILITY_CONVENTION,
  PT_COMPONENT_INVENTORY_CONVENTION,
  PT_PHASE_TRANSITION_CONVENTION,
  PT_SERVICE_BOUNDARY_CONVENTION,
  PT_WIRE_CONTRACT,
  type PtBackendDescriptor,
  type PtCapabilityDiscovery,
  type PtFlashRequest,
  type PtFlashResponse,
} from '../../../../frontend/src/domain/flash';

export const PT_ANDROID_BRIDGE_CONVENTION =
  'MPMC/PT/android-product-bridge/v1' as const;

const DEFAULT_DISCOVERY_TIMEOUT_MS = 10_000;
const DEFAULT_SOLVE_TIMEOUT_MS = 120_000;

interface NativeJsonReply {
  json: string;
}

export interface MpmcPtNativePlugin {
  discover(): Promise<NativeJsonReply>;
  solve(options: { requestJson: string }): Promise<NativeJsonReply>;
}

class AndroidBridgeContractError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'AndroidBridgeContractError';
  }
}

type JsonObject = Record<string, unknown>;

function record(value: unknown, label: string): JsonObject {
  if (typeof value !== 'object' || value === null || Array.isArray(value)) {
    throw new AndroidBridgeContractError(`${label} must be a JSON object.`);
  }
  return value as JsonObject;
}

function text(value: unknown, label: string): string {
  if (typeof value !== 'string') {
    throw new AndroidBridgeContractError(`${label} must be a string.`);
  }
  return value;
}

function finiteNumber(value: unknown, label: string): number {
  if (typeof value !== 'number' || !Number.isFinite(value)) {
    throw new AndroidBridgeContractError(`${label} must be a finite number.`);
  }
  return value;
}

function integer(value: unknown, label: string): number {
  const result = finiteNumber(value, label);
  if (!Number.isInteger(result) || result < 0) {
    throw new AndroidBridgeContractError(`${label} must be a non-negative integer.`);
  }
  return result;
}

function booleanValue(value: unknown, label: string): boolean {
  if (typeof value !== 'boolean') {
    throw new AndroidBridgeContractError(`${label} must be boolean.`);
  }
  return value;
}

function array(value: unknown, label: string): readonly unknown[] {
  if (!Array.isArray(value)) {
    throw new AndroidBridgeContractError(`${label} must be an array.`);
  }
  return value;
}

function expectConvention(value: unknown, expected: string, label: string): void {
  if (value !== expected) {
    throw new AndroidBridgeContractError(
      `${label} does not match the supported Android PT contract.`,
    );
  }
}

function validateBackendDescriptor(value: unknown, label: string): PtBackendDescriptor {
  const descriptor = record(value, label);
  expectConvention(descriptor.convention, PT_CAPABILITY_CONVENTION, `${label}.convention`);
  const configuredBackendId = text(
    descriptor.configuredBackendId,
    `${label}.configuredBackendId`,
  );
  if (configuredBackendId.length === 0) {
    throw new AndroidBridgeContractError(`${label}.configuredBackendId is empty.`);
  }

  const capability = record(descriptor.capability, `${label}.capability`);
  expectConvention(
    capability.convention,
    PT_BACKEND_CAPABILITY_CONVENTION,
    `${label}.capability.convention`,
  );
  for (const key of [
    'backendId',
    'modelProfile',
    'algorithmProfile',
    'publicationProfile',
    'configurationProfile',
    'datasetId',
    'revision',
    'phaseMetadataNamespace',
  ] as const) {
    text(capability[key], `${label}.capability.${key}`);
  }

  const componentIds = array(
    capability.componentIds,
    `${label}.capability.componentIds`,
  ).map((item, index) => text(item, `${label}.capability.componentIds[${index}]`));
  if (componentIds.length === 0 || new Set(componentIds).size !== componentIds.length) {
    throw new AndroidBridgeContractError(
      `${label}.capability.componentIds must be nonempty and unique.`,
    );
  }
  array(capability.supportedPhaseCounts, `${label}.capability.supportedPhaseCounts`).forEach(
    (item, index) => integer(item, `${label}.capability.supportedPhaseCounts[${index}]`),
  );
  array(capability.scalarSettings, `${label}.capability.scalarSettings`).forEach(
    (item, index) => {
      const setting = record(item, `${label}.capability.scalarSettings[${index}]`);
      text(setting.id, `${label}.capability.scalarSettings[${index}].id`);
      finiteNumber(setting.value, `${label}.capability.scalarSettings[${index}].value`);
      text(setting.unit, `${label}.capability.scalarSettings[${index}].unit`);
    },
  );

  const transition = record(
    capability.transitionCapability,
    `${label}.capability.transitionCapability`,
  );
  expectConvention(
    transition.convention,
    PT_PHASE_TRANSITION_CONVENTION,
    `${label}.capability.transitionCapability.convention`,
  );
  array(transition.edges, `${label}.capability.transitionCapability.edges`).forEach(
    (item, index) => {
      const edge = record(item, `${label}.capability.transitionCapability.edges[${index}]`);
      integer(edge.sourcePhaseCount, `${label}.edge.sourcePhaseCount`);
      integer(edge.targetPhaseCount, `${label}.edge.targetPhaseCount`);
      const support = text(edge.support, `${label}.edge.support`);
      if (support !== 'detection_only' && support !== 'fresh_target_resolve') {
        throw new AndroidBridgeContractError(`${label}.edge.support is unsupported.`);
      }
      booleanValue(edge.requiresFreshTargetSolve, `${label}.edge.requiresFreshTargetSolve`);
    },
  );
  for (const key of [
    'performsInitialStabilitySearch',
    'performsFinalPhaseSetReview',
    'performsBoundaryNeighborResolve',
    'globalStabilityProven',
  ] as const) {
    booleanValue(capability[key], `${label}.capability.${key}`);
  }

  const inventory = record(descriptor.componentInventory, `${label}.componentInventory`);
  expectConvention(
    inventory.convention,
    PT_COMPONENT_INVENTORY_CONVENTION,
    `${label}.componentInventory.convention`,
  );
  const inventoryIds = array(
    inventory.components,
    `${label}.componentInventory.components`,
  ).map((item, index) => {
    const component = record(item, `${label}.componentInventory.components[${index}]`);
    const id = text(component.componentId, `${label}.componentInventory.componentId`);
    const feedIndex = integer(component.feedIndex, `${label}.componentInventory.feedIndex`);
    if (feedIndex !== index) {
      throw new AndroidBridgeContractError(`${label}.componentInventory feed order is invalid.`);
    }
    return id;
  });
  if (
    inventoryIds.length !== componentIds.length ||
    inventoryIds.some((id, index) => id !== componentIds[index])
  ) {
    throw new AndroidBridgeContractError(
      `${label} capability and component inventory orders disagree.`,
    );
  }

  return descriptor as unknown as PtBackendDescriptor;
}

function parseJson(json: string, label: string): JsonObject {
  try {
    return record(JSON.parse(json) as unknown, label);
  } catch (cause) {
    if (cause instanceof AndroidBridgeContractError) {
      throw cause;
    }
    throw new AndroidBridgeContractError(`${label} is not valid JSON.`);
  }
}

function parseDiscovery(json: string): PtCapabilityDiscovery {
  const payload = parseJson(json, 'Android discovery reply');
  expectConvention(
    payload.bridgeConvention,
    PT_ANDROID_BRIDGE_CONVENTION,
    'Android discovery bridgeConvention',
  );
  if (typeof payload.bridgeError === 'string') {
    throw new Error(payload.bridgeError);
  }
  expectConvention(
    payload.serviceBoundaryConvention,
    PT_SERVICE_BOUNDARY_CONVENTION,
    'Android discovery serviceBoundaryConvention',
  );
  expectConvention(
    payload.capabilityConvention,
    PT_CAPABILITY_CONVENTION,
    'Android discovery capabilityConvention',
  );
  const backends = array(payload.backends, 'Android discovery backends').map(
    (item, index) => validateBackendDescriptor(item, `Android discovery backends[${index}]`),
  );
  const ids = backends.map((backend) => backend.configuredBackendId);
  if (new Set(ids).size !== ids.length) {
    throw new AndroidBridgeContractError('Android discovery returned duplicate backend IDs.');
  }
  return {
    wireContract: PT_WIRE_CONTRACT,
    serviceBoundaryConvention: PT_SERVICE_BOUNDARY_CONVENTION,
    capabilityConvention: PT_CAPABILITY_CONVENTION,
    backends,
  };
}

function validateFlashResponse(
  json: string,
  expectedBackend: PtBackendDescriptor,
): PtFlashResponse {
  const payload = parseJson(json, 'Android solve reply');
  expectConvention(
    payload.bridgeConvention,
    PT_ANDROID_BRIDGE_CONVENTION,
    'Android solve bridgeConvention',
  );
  if (typeof payload.bridgeError === 'string') {
    throw new Error(payload.bridgeError);
  }
  const response = record(payload.response, 'Android solve response');
  const kind = text(response.kind, 'Android solve response.kind');
  if (kind === 'service_error') {
    const error = record(response.error, 'Android solve service error');
    text(error.code, 'Android solve service error.code');
    text(error.field, 'Android solve service error.field');
    text(error.diagnostic, 'Android solve service error.diagnostic');
    return response as unknown as PtFlashResponse;
  }
  if (kind !== 'result') {
    throw new AndroidBridgeContractError('Android solve response kind is unsupported.');
  }

  const result = record(response.result, 'Android solve result');
  const outcome = text(result.outcome, 'Android solve result.outcome');
  if (
    outcome !== 'accepted' &&
    outcome !== 'phase_set_unstable' &&
    outcome !== 'indeterminate'
  ) {
    throw new AndroidBridgeContractError('Android solve result outcome is unsupported.');
  }
  finiteNumber(result.pressurePa, 'Android solve result.pressurePa');
  finiteNumber(result.temperatureK, 'Android solve result.temperatureK');
  booleanValue(result.globalStabilityProven, 'Android solve result.globalStabilityProven');
  booleanValue(result.morphologyResolved, 'Android solve result.morphologyResolved');
  text(result.diagnostic, 'Android solve result.diagnostic');
  array(result.feed, 'Android solve result.feed');
  array(result.phases, 'Android solve result.phases');

  const provenance = record(result.provenance, 'Android solve result.provenance');
  const backend = validateBackendDescriptor(
    provenance.backend,
    'Android solve result.provenance.backend',
  );
  if (backend.configuredBackendId !== expectedBackend.configuredBackendId) {
    throw new AndroidBridgeContractError(
      'Android solve provenance does not match the selected capability snapshot.',
    );
  }
  return response as unknown as PtFlashResponse;
}

function timeoutMs(options: PtClientCallOptions | undefined, fallback: number): number {
  const configured = options?.timeoutMs ?? fallback;
  return Number.isFinite(configured) && configured > 0 ? configured : fallback;
}

async function boundedNativeCall<T>(
  operation: Promise<T>,
  options: PtClientCallOptions | undefined,
  fallbackTimeoutMs: number,
): Promise<T> {
  if (options?.signal?.aborted === true) {
    throw new Error('Android PT request was canceled before native execution.');
  }
  return await new Promise<T>((resolve, reject) => {
    let settled = false;
    const settle = (callback: () => void) => {
      if (!settled) {
        settled = true;
        clearTimeout(timer);
        options?.signal?.removeEventListener('abort', onAbort);
        callback();
      }
    };
    const onAbort = () => settle(() => reject(new Error('Android PT request was canceled.')));
    const timer = globalThis.setTimeout(
      () => settle(() => reject(new Error('Android PT native call exceeded its UI deadline.'))),
      timeoutMs(options, fallbackTimeoutMs),
    );
    options?.signal?.addEventListener('abort', onAbort, { once: true });
    operation.then(
      (value) => settle(() => resolve(value)),
      (cause) => settle(() => reject(cause)),
    );
  });
}

export class AndroidFlashClient implements FlashClient {
  readonly configured = true;
  readonly endpoint = 'android://in-process-pt-service';

  constructor(private readonly plugin: MpmcPtNativePlugin) {}

  async discoverPtCapabilities(
    options?: PtClientCallOptions,
  ): Promise<PtCapabilityDiscovery> {
    const reply = await boundedNativeCall(
      this.plugin.discover(),
      options,
      DEFAULT_DISCOVERY_TIMEOUT_MS,
    );
    return parseDiscovery(reply.json);
  }

  async solvePtFlash(
    request: PtFlashRequest,
    expectedBackend: PtBackendDescriptor,
    options?: PtClientCallOptions,
  ): Promise<PtFlashResponse> {
    if (request.configuredBackendId !== expectedBackend.configuredBackendId) {
      throw new AndroidBridgeContractError(
        'Selected backend does not match the Android capability snapshot.',
      );
    }
    const reply = await boundedNativeCall(
      this.plugin.solve({ requestJson: JSON.stringify(request) }),
      options,
      DEFAULT_SOLVE_TIMEOUT_MS,
    );
    return validateFlashResponse(reply.json, expectedBackend);
  }
}

export function createAndroidFlashClient(): FlashClient {
  return new AndroidFlashClient(registerPlugin<MpmcPtNativePlugin>('MpmcPt'));
}
