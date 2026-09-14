import {
  Code,
  ConnectError,
  createClient,
  type CallOptions,
  type Client,
} from '@connectrpc/connect';
import { createGrpcTransport } from '@connectrpc/connect-node';

import type {
  PtBackendDescriptor,
  PtCapabilityDiscovery,
  PtFlashRequest,
  PtFlashResponse,
} from '../src/domain/flash';
import { PtFlashService } from '../src/gen/mpmc/runtime/v1/pt_service_pb';
import {
  mapDiscoveryResponse,
  mapSolveResponse,
  toWireSolveRequest,
} from '../src/api/ptWire';
import { PtHostSession } from './hostSession';

type PtClient = Client<typeof PtFlashService>;

interface ConnectedClient {
  readonly client: PtClient;
  readonly authorization: Headers;
}

export class PtDesktopGateway {
  private connected: Promise<ConnectedClient> | null = null;
  private discovery: PtCapabilityDiscovery | null = null;

  constructor(private readonly host: PtHostSession) {}

  private connect(): Promise<ConnectedClient> {
    if (this.connected !== null) {
      return this.connected;
    }
    this.connected = this.host.start().then((connection) => {
      const transport = createGrpcTransport({
        baseUrl: connection.baseUrl,
        interceptors: [],
      });
      return {
        client: createClient(PtFlashService, transport),
        authorization: new Headers({
          authorization: `Bearer ${connection.bearerToken}`,
        }),
      };
    });
    return this.connected;
  }

  private callOptions(
    connected: ConnectedClient,
    timeoutMs: number,
    signal?: AbortSignal,
  ): CallOptions {
    const options: CallOptions = {
      headers: connected.authorization,
      timeoutMs,
    };
    if (signal !== undefined) {
      options.signal = signal;
    }
    return options;
  }

  async discover(signal?: AbortSignal): Promise<PtCapabilityDiscovery> {
    const connected = await this.connect();
    const wire = await connected.client.discoverPtCapabilities(
      {},
      this.callOptions(connected, 10_000, signal),
    );
    const discovery = mapDiscoveryResponse(wire);
    this.discovery = discovery;
    return discovery;
  }

  async solve(request: PtFlashRequest, signal?: AbortSignal): Promise<PtFlashResponse> {
    const connected = await this.connect();
    const discovery = this.discovery ?? (await this.discover(signal));
    const expected = discovery.backends.find(
      (backend) => backend.configuredBackendId === request.configuredBackendId,
    );
    if (expected === undefined) {
      throw new ConnectError(
        'The selected backend is absent from the desktop discovery snapshot.',
        Code.NotFound,
      );
    }
    const wire = await connected.client.solvePtFlash(
      toWireSolveRequest(request),
      this.callOptions(connected, 120_000, signal),
    );
    return mapSolveResponse(wire, expected, request);
  }

  backend(configuredBackendId: string): PtBackendDescriptor | undefined {
    return this.discovery?.backends.find(
      (backend) => backend.configuredBackendId === configuredBackendId,
    );
  }

  stop(): Promise<void> {
    this.connected = null;
    this.discovery = null;
    return this.host.stop();
  }
}
