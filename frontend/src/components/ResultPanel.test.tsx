import { renderToStaticMarkup } from 'react-dom/server';
import { describe, expect, it } from 'vitest';

import { PtComputationOutcome } from '../gen/mpmc/runtime/v1/pt_service_pb';
import {
  domainRequest,
  wireDiscovery,
  wireResult,
  wireServiceError,
} from '../test/ptWireFixtures';
import { mapDiscoveryResponse, mapSolveResponse } from '../api/ptWire';
import { ResultPanel } from './ResultPanel';

function discoveredBackend() {
  const backend = mapDiscoveryResponse(wireDiscovery()).backends[0];
  if (backend === undefined) {
    throw new Error('fixture backend missing');
  }
  return backend;
}

describe('PT result presentation boundaries', () => {
  it('renders indeterminate as a scientific result, not a service failure', () => {
    const response = mapSolveResponse(
      wireResult(PtComputationOutcome.INDETERMINATE),
      discoveredBackend(),
      domainRequest(),
    );
    const html = renderToStaticMarkup(
      <ResultPanel response={response} clientFailure={null} />,
    );
    expect(html).toContain('Scientific decision remains indeterminate');
    expect(html).toContain('data-outcome="indeterminate"');
    expect(html).not.toContain('class="phase-card"');
    expect(html).not.toContain('PT service error');
  });

  it('renders the protobuf error arm as a service error', () => {
    const response = mapSolveResponse(
      wireServiceError(),
      discoveredBackend(),
      domainRequest(),
    );
    const html = renderToStaticMarkup(
      <ResultPanel response={response} clientFailure={null} />,
    );
    expect(html).toContain('PT service error');
    expect(html).toContain('Request rejected by service');
    expect(html).toContain('pressure_pa');
    expect(html).not.toContain('class="scientific-outcome"');
  });

  it('renders gRPC failure without inventing a scientific result', () => {
    const html = renderToStaticMarkup(
      <ResultPanel
        response={null}
        clientFailure={{
          kind: 'transport',
          title: 'Solve RPC deadline exceeded',
          message: 'request timed out',
        }}
      />,
    );
    expect(html).toContain('Client / RPC failure');
    expect(html).toContain('No scientific outcome was received');
    expect(html).not.toContain('Scientific computation result');
  });
});
