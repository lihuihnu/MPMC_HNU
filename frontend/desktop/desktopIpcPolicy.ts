import { Buffer } from 'node:buffer';

import { Code, ConnectError } from '@connectrpc/connect';

import type { PtFlashRequest } from '../src/domain/flash';

const MAX_ACTIVE_DESKTOP_CALLS = 4;
const MAX_DESKTOP_REQUEST_BYTES = 64 * 1024;
const MAX_COMPONENTS = 256;
const MAX_IDENTIFIER_BYTES = 256;
const requestIdPattern = /^[A-Za-z0-9-]{1,128}$/u;

function isRecord(value: unknown): value is Record<string, unknown> {
  return typeof value === 'object' && value !== null && !Array.isArray(value);
}

export function validDesktopRequestId(value: unknown): value is string {
  return typeof value === 'string' && requestIdPattern.test(value);
}

export function parseDesktopSolveRequest(value: unknown): PtFlashRequest {
  let serialized: string | undefined;
  try {
    serialized = JSON.stringify(value);
  } catch {
    throw new ConnectError('The desktop solve request is not serializable.', Code.InvalidArgument);
  }
  if (typeof serialized !== 'string' || !isRecord(value)) {
    throw new ConnectError('The desktop solve request has an invalid wire shape.', Code.InvalidArgument);
  }
  if (Buffer.byteLength(serialized, 'utf8') > MAX_DESKTOP_REQUEST_BYTES) {
    throw new ConnectError(
      'The desktop solve request exceeds the process boundary.',
      Code.ResourceExhausted,
    );
  }
  if (
    typeof value.configuredBackendId !== 'string' ||
    Buffer.byteLength(value.configuredBackendId, 'utf8') > MAX_IDENTIFIER_BYTES ||
    typeof value.pressurePa !== 'number' ||
    typeof value.temperatureK !== 'number' ||
    !Array.isArray(value.feed) ||
    value.feed.length > MAX_COMPONENTS
  ) {
    throw new ConnectError('The desktop solve request has an invalid wire shape.', Code.InvalidArgument);
  }
  for (const entry of value.feed) {
    if (
      !isRecord(entry) ||
      typeof entry.componentId !== 'string' ||
      Buffer.byteLength(entry.componentId, 'utf8') > MAX_IDENTIFIER_BYTES ||
      typeof entry.moleFraction !== 'number'
    ) {
      throw new ConnectError('The desktop feed has an invalid wire shape.', Code.InvalidArgument);
    }
  }
  return value as unknown as PtFlashRequest;
}

export class PtDesktopCallGate {
  private readonly active = new Map<string, AbortController>();

  begin(id: unknown): AbortController {
    if (!validDesktopRequestId(id)) {
      throw new ConnectError('Invalid desktop PT request identifier.', Code.InvalidArgument);
    }
    if (this.active.has(id)) {
      throw new ConnectError('Duplicate desktop PT request identifier.', Code.AlreadyExists);
    }
    if (this.active.size >= MAX_ACTIVE_DESKTOP_CALLS) {
      throw new ConnectError('Too many desktop PT calls are active.', Code.ResourceExhausted);
    }
    const controller = new AbortController();
    this.active.set(id, controller);
    return controller;
  }

  finish(id: unknown, controller: AbortController | undefined): void {
    if (
      controller !== undefined &&
      validDesktopRequestId(id) &&
      this.active.get(id) === controller
    ) {
      this.active.delete(id);
    }
  }

  cancel(id: unknown): void {
    if (validDesktopRequestId(id)) {
      this.active.get(id)?.abort();
    }
  }

  abortAll(): void {
    for (const controller of this.active.values()) {
      controller.abort();
    }
    this.active.clear();
  }

  get size(): number {
    return this.active.size;
  }
}
