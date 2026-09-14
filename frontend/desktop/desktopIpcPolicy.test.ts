import { Code, ConnectError } from '@connectrpc/connect';
import { describe, expect, it } from 'vitest';

import {
  parseDesktopSolveRequest,
  PtDesktopCallGate,
} from './desktopIpcPolicy';

function validRequest() {
  return {
    configuredBackendId: 'backend.fixture',
    pressurePa: 1e6,
    temperatureK: 300,
    feed: [{ componentId: 'component', moleFraction: 1 }],
  };
}

describe('PT desktop IPC process policy', () => {
  it('does not release an original call when a duplicate ID is rejected', () => {
    const gate = new PtDesktopCallGate();
    const original = gate.begin('same-id');
    expect(() => gate.begin('same-id')).toThrowError(
      expect.objectContaining<Partial<ConnectError>>({ code: Code.AlreadyExists }),
    );
    gate.finish('same-id', undefined);
    expect(gate.size).toBe(1);
    gate.cancel('same-id');
    expect(original.signal.aborted).toBe(true);
    gate.finish('same-id', original);
    expect(gate.size).toBe(0);
  });

  it('bounds active calls and aborts every owned call during teardown', () => {
    const gate = new PtDesktopCallGate();
    const controllers = ['a', 'b', 'c', 'd'].map((id) => gate.begin(id));
    expect(() => gate.begin('e')).toThrowError(
      expect.objectContaining<Partial<ConnectError>>({ code: Code.ResourceExhausted }),
    );
    gate.abortAll();
    expect(controllers.every((controller) => controller.signal.aborted)).toBe(true);
    expect(gate.size).toBe(0);
  });

  it('preserves valid request values and rejects malformed or oversized IPC input', () => {
    const request = validRequest();
    expect(parseDesktopSolveRequest(request)).toBe(request);
    expect(() => parseDesktopSolveRequest(undefined)).toThrowError(
      expect.objectContaining<Partial<ConnectError>>({ code: Code.InvalidArgument }),
    );
    expect(() => parseDesktopSolveRequest({ ...request, feed: [1] })).toThrowError(
      expect.objectContaining<Partial<ConnectError>>({ code: Code.InvalidArgument }),
    );
    expect(() =>
      parseDesktopSolveRequest({
        ...request,
        configuredBackendId: 'x'.repeat(70 * 1024),
      }),
    ).toThrowError(
      expect.objectContaining<Partial<ConnectError>>({ code: Code.ResourceExhausted }),
    );
  });
});
