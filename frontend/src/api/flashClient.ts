import {
  Code,
  ConnectError,
  createClient,
  type CallOptions,
  type Client,
} from '@connectrpc/connect';
import { createGrpcWebTransport } from '@connectrpc/connect-web';

import { PtFlashService } from '../gen/mpmc/runtime/v1/pt_service_pb';
import type {
  PtBackendDescriptor,
  PtCapabilityDiscovery,
  PtFlashRequest,
  PtFlashResponse,
} from '../domain/flash';
import {
  mapDiscoveryResponse,
  mapSolveResponse,
  toWireSolveRequest,
} from './ptWire';

const DEFAULT_DISCOVERY_TIMEOUT_MS = 10_000;
const DEFAULT_SOLVE_TIMEOUT_MS = 120_000;

export interface PtClientCallOptions {
  signal?: AbortSignal;
  timeoutMs?: number;
}

export interface FlashClient {
  readonly configured: boolean;
  readonly endpoint: string | null;
  discoverPtCapabilities(options?: PtClientCallOptions): Promise<PtCapabilityDiscovery>;
  solvePtFlash(
    request: PtFlashRequest,
    expectedBackend: PtBackendDescriptor,
    options?: PtClientCallOptions,
  ): Promise<PtFlashResponse>;
}

export class BackendNotConfiguredError extends Error {
  constructor() {
    super('The frontend build has no configured MPMC_HNU gRPC-Web endpoint.');
    this.name = 'BackendNotConfiguredError';
  }
}

export class PtTransportError extends Error {
  readonly code: Code;

  constructor(code: Code, message: string) {
    super(message || 'The PT RPC failed.');
    this.name = 'PtTransportError';
    this.code = code;
  }
}

export type PtFlashRpc = Pick<
  Client<typeof PtFlashService>,
  'discoverPtCapabilities' | 'solvePtFlash'
>;

function callOptions(
  options: PtClientCallOptions | undefined,
  defaultTimeoutMs: number,
): CallOptions {
  const result: CallOptions = {
    timeoutMs: options?.timeoutMs ?? defaultTimeoutMs,
  };
  if (options?.signal !== undefined) {
    result.signal = options.signal;
  }
  return result;
}

async function invoke<T>(operation: () => Promise<T>): Promise<T> {
  try {
    return await operation();
  } catch (cause) {
    const error = ConnectError.from(cause);
    throw new PtTransportError(
      error.code,
      error.rawMessage || 'The gRPC-Web request failed.',
    );
  }
}

export class GrpcWebFlashClient implements FlashClient {
  readonly configured = true;
  readonly endpoint: string;

  constructor(
    private readonly rpc: PtFlashRpc,
    endpoint: string,
  ) {
    this.endpoint = endpoint;
  }

  async discoverPtCapabilities(
    options?: PtClientCallOptions,
  ): Promise<PtCapabilityDiscovery> {
    const response = await invoke(() =>
      this.rpc.discoverPtCapabilities(
        {},
        callOptions(options, DEFAULT_DISCOVERY_TIMEOUT_MS),
      ),
    );
    return mapDiscoveryResponse(response);
  }

  async solvePtFlash(
    request: PtFlashRequest,
    expectedBackend: PtBackendDescriptor,
    options?: PtClientCallOptions,
  ): Promise<PtFlashResponse> {
    const response = await invoke(() =>
      this.rpc.solvePtFlash(
        toWireSolveRequest(request),
        callOptions(options, DEFAULT_SOLVE_TIMEOUT_MS),
      ),
    );
    return mapSolveResponse(response, expectedBackend, request);
  }
}

export function createGrpcWebFlashClient(baseUrl: string): FlashClient {
  const configured = baseUrl.trim();
  if (configured.length === 0) {
    return unconfiguredFlashClient;
  }
  const endpoint = /^\/+$/u.test(configured)
    ? '/'
    : configured.replace(/\/+$/, '');
  const transport = createGrpcWebTransport({
    baseUrl: endpoint,
    useBinaryFormat: true,
  });
  return new GrpcWebFlashClient(createClient(PtFlashService, transport), endpoint);
}

export const unconfiguredFlashClient: FlashClient = {
  configured: false,
  endpoint: null,
  async discoverPtCapabilities(): Promise<PtCapabilityDiscovery> {
    throw new BackendNotConfiguredError();
  },
  async solvePtFlash(): Promise<PtFlashResponse> {
    throw new BackendNotConfiguredError();
  },
};

export function transportCodeLabel(code: Code): string {
  switch (code) {
    case Code.Canceled:
      return 'canceled';
    case Code.Unknown:
      return 'unknown';
    case Code.InvalidArgument:
      return 'invalid argument';
    case Code.DeadlineExceeded:
      return 'deadline exceeded';
    case Code.NotFound:
      return 'not found';
    case Code.AlreadyExists:
      return 'already exists';
    case Code.PermissionDenied:
      return 'permission denied';
    case Code.ResourceExhausted:
      return 'resource exhausted';
    case Code.FailedPrecondition:
      return 'failed precondition';
    case Code.Aborted:
      return 'aborted';
    case Code.OutOfRange:
      return 'out of range';
    case Code.Unimplemented:
      return 'unimplemented';
    case Code.Internal:
      return 'internal';
    case Code.Unavailable:
      return 'unavailable';
    case Code.DataLoss:
      return 'data loss';
    case Code.Unauthenticated:
      return 'unauthenticated';
  }
}
