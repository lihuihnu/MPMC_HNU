import { create, fromBinary, toBinary } from '@bufbuild/protobuf';
import { describe, expect, it } from 'vitest';

import {
  DiscoverPtCapabilitiesResponseSchema,
  PtComputationOutcome,
  PtFlashService,
  SolvePtFlashResponseSchema,
} from '../gen/mpmc/runtime/v1/pt_service_pb';
import {
  PT_SERVICE_REQUEST_CONVENTION,
  PT_WIRE_CONTRACT,
  type PtFlashRequest,
} from '../domain/flash';
import {
  domainRequest,
  wireDiscovery,
  wireResult,
  wireServiceError,
} from '../test/ptWireFixtures';
import {
  PtWireContractError,
  mapDiscoveryResponse,
  mapSolveResponse,
  toWireSolveRequest,
} from './ptWire';

describe('mpmc.runtime.v1 protobuf mapping', () => {
  it('freezes the versioned service and unary method descriptors', () => {
    expect(PtFlashService.typeName).toBe('mpmc.runtime.v1.PtFlashService');
    expect(PtFlashService.method.discoverPtCapabilities.methodKind).toBe('unary');
    expect(PtFlashService.method.solvePtFlash.methodKind).toBe('unary');
  });

  it('round-trips capability discovery through protobuf binary encoding', () => {
    const message = wireDiscovery();
    const decoded = fromBinary(
      DiscoverPtCapabilitiesResponseSchema,
      toBinary(DiscoverPtCapabilitiesResponseSchema, message),
    );
    const discovery = mapDiscoveryResponse(decoded);
    expect(discovery.backends).toHaveLength(1);
    expect(discovery.backends[0]?.configuredBackendId).toBe('pr76-default');
    expect(
      discovery.backends[0]?.componentInventory.components.map(
        (component) => component.componentId,
      ),
    ).toEqual(['methane', 'carbon-dioxide']);
  });

  it('rejects a capability whose inventory snapshot drifts from component IDs', () => {
    const invalid = wireDiscovery();
    const inventory = invalid.backends[0]?.componentInventory;
    if (inventory === undefined || inventory.components[1] === undefined) {
      throw new Error('fixture inventory missing');
    }
    inventory.components[1].componentId = 'water';
    expect(() => mapDiscoveryResponse(invalid)).toThrow(/inventory ID\/order/);
  });

  it('maps request IDs and values without normalization or order changes', () => {
    const request: PtFlashRequest = {
      configuredBackendId: 'pr76-default',
      pressurePa: 10e6,
      temperatureK: 350,
      feed: [
        { componentId: 'carbon-dioxide', moleFraction: 0.3 },
        { componentId: 'methane', moleFraction: 0.7 },
      ],
    };
    const wire = toWireSolveRequest(request);
    expect(wire.wireContract).toBe(PT_WIRE_CONTRACT);
    expect(wire.serviceRequestConvention).toBe(PT_SERVICE_REQUEST_CONVENTION);
    expect(wire.feed).toEqual(request.feed);
  });

  it('maps a variable accepted phase vector with provenance', () => {
    const backend = mapDiscoveryResponse(wireDiscovery()).backends[0];
    expect(backend).toBeDefined();
    const wire = wireResult();
    const decoded = fromBinary(
      SolvePtFlashResponseSchema,
      toBinary(SolvePtFlashResponseSchema, wire),
    );
    const response = mapSolveResponse(decoded, backend!, domainRequest());
    expect(response.kind).toBe('result');
    if (response.kind === 'result') {
      expect(response.result.outcome).toBe('accepted');
      expect(response.result.phases).toHaveLength(1);
      expect(response.result.phases[0]?.providerBranch).toBe('0');
      expect(response.result.provenance.backend).toEqual(backend);
    }
  });

  it('keeps indeterminate on the computation-result branch with no phases', () => {
    const backend = mapDiscoveryResponse(wireDiscovery()).backends[0];
    expect(backend).toBeDefined();
    const response = mapSolveResponse(
      wireResult(PtComputationOutcome.INDETERMINATE),
      backend!,
      domainRequest(),
    );
    expect(response).toMatchObject({
      kind: 'result',
      result: { outcome: 'indeterminate', phases: [] },
    });
  });

  it('keeps a PtService error on the error branch', () => {
    const backend = mapDiscoveryResponse(wireDiscovery()).backends[0];
    expect(backend).toBeDefined();
    expect(mapSolveResponse(wireServiceError(), backend!, domainRequest())).toEqual({
      kind: 'service_error',
      error: {
        code: 'invalid_pressure',
        field: 'pressure_pa',
        diagnostic: 'pressure must be positive',
      },
    });
  });

  it('rejects a non-accepted response that tries to publish phases', () => {
    const backend = mapDiscoveryResponse(wireDiscovery()).backends[0];
    expect(backend).toBeDefined();
    const invalid = wireResult(PtComputationOutcome.INDETERMINATE);
    const accepted = wireResult();
    if (invalid.payload.case === 'result' && accepted.payload.case === 'result') {
      invalid.payload.value.phases.push(...accepted.payload.value.phases);
    }
    expect(() => mapSolveResponse(invalid, backend!, domainRequest())).toThrow(
      PtWireContractError,
    );
  });

  it('rejects an unknown computation enum rather than guessing a state', () => {
    const backend = mapDiscoveryResponse(wireDiscovery()).backends[0];
    expect(backend).toBeDefined();
    const invalid = wireResult();
    if (invalid.payload.case === 'result') {
      invalid.payload.value.outcome = 99 as PtComputationOutcome;
    }
    expect(() => mapSolveResponse(invalid, backend!, domainRequest())).toThrow(
      /unknown computation outcome/,
    );
  });

  it('rejects a response with no result/error oneof arm', () => {
    const backend = mapDiscoveryResponse(wireDiscovery()).backends[0];
    expect(backend).toBeDefined();
    const invalid = create(SolvePtFlashResponseSchema, {
      wireContract: PT_WIRE_CONTRACT,
      serviceResultConvention: 'PT/service-result/v1',
    });
    expect(() => mapSolveResponse(invalid, backend!, domainRequest())).toThrow(
      /neither result nor service error/,
    );
  });
});
