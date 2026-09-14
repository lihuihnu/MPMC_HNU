import { describe, expect, it } from 'vitest';

import { parsePtReadyEvent } from './hostSession';

describe('PT desktop host readiness parser', () => {
  it('accepts only the exact three-backend loopback child readiness shape', () => {
    expect(
      parsePtReadyEvent(
        JSON.stringify({
          event: 'pt_process_ready',
          selected_port: 49152,
          configured_backends: 3,
          snapshot_bundle: 'MPMC/PT/repository-curated-literature-snapshots/v1',
        }),
      ),
    ).toEqual({
      event: 'pt_process_ready',
      selected_port: 49152,
      configured_backends: 3,
    });
  });

  it('rejects malformed, fixed-zero-port, or incomplete readiness lines', () => {
    expect(parsePtReadyEvent('not-json')).toBeNull();
    expect(
      parsePtReadyEvent(
        JSON.stringify({
          event: 'pt_process_ready',
          selected_port: 0,
          configured_backends: 3,
        }),
      ),
    ).toBeNull();
    expect(
      parsePtReadyEvent(
        JSON.stringify({
          event: 'pt_process_ready',
          selected_port: 50051,
          configured_backends: 2,
        }),
      ),
    ).toBeNull();
  });
});
