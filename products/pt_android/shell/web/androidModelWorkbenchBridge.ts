import { type JsonObject, type JsonValue } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';
import { registerPlugin } from '@capacitor/core';

import {
  MODEL_WORKBENCH_CONVENTION,
  type ModelWorkbenchBridge,
  type ModelWorkbenchReply,
} from '../../../../frontend/src/api/modelWorkbenchContract';

const MAX_REQUEST_BYTES = 64 * 1024;
const MAX_REPLY_BYTES = 4 * 1024 * 1024;

type JsonRecord = JsonObject;
interface NativeRecordsReply { records: string[] }
interface AndroidModelNativePlugin {
  modelApply(options: { requestId: string; recordsJson: string }): Promise<NativeRecordsReply>;
  modelSolve(options: { requestId: string; recordsJson: string }): Promise<NativeRecordsReply>;
  modelRelease(options: { requestId: string }): Promise<NativeRecordsReply>;
  modelCancel(options: { requestId: string }): Promise<void>;
}

const plugin = registerPlugin<AndroidModelNativePlugin>('MpmcPt');

// Protobuf JSON flattens oneof `definition.parameters` to the direct `pr76`
// field. The native value object keeps that oneof as a C++ variant named
// parameters, so only this mechanical path is translated by the adapter.
function nativePath(path: string): string {
  return path.replace(/^definition\.pr76(?=\.|$)/u, 'definition.parameters.pr76');
}
function rendererPath(path: string): string {
  return path.replace(/^definition\.parameters\.pr76(?=\.|$)/u, 'definition.pr76');
}

function flatten(value: JsonValue): string[] {
  const records: string[] = [];
  const push = (path: string, kind: string, raw: string) => {
    records.push(nativePath(path), kind, raw);
  };
  const visit = (current: JsonValue, path: string): void => {
    if (current === null) throw new Error('Android model bridge does not accept null values.');
    if (Array.isArray(current)) {
      push(path, 'array', String(current.length));
      current.forEach((item, index) => visit(item, `${path}[${index}]`));
      return;
    }
    if (typeof current === 'object') {
      for (const [key, item] of Object.entries(current)) {
        if (!key || /[.[\]]/u.test(key) || item === undefined) {
          throw new Error('Android model bridge received an invalid object field.');
        }
        visit(item as JsonValue, path ? `${path}.${key}` : key);
      }
      return;
    }
    if (typeof current === 'string') push(path, 'string', current);
    else if (typeof current === 'number' && Number.isFinite(current)) {
      push(path, 'number', String(current));
    } else if (typeof current === 'boolean') push(path, 'boolean', current ? 'true' : 'false');
    else throw new Error('Android model bridge received an unsupported scalar value.');
  };
  visit(value, '');
  return records;
}

interface Segment { key?: string; index?: number }
function segments(path: string): Segment[] {
  if (!path) return [];
  const result: Segment[] = [];
  for (const part of path.split('.')) {
    const match = /^([^\[]+)((?:\[\d+\])*)$/u.exec(part);
    if (!match) throw new Error('Android model bridge returned an invalid field path.');
    result.push({ key: match[1]! });
    for (const item of match[2]!.matchAll(/\[(\d+)\]/gu)) {
      result.push({ index: Number(item[1]) });
    }
  }
  return result;
}

function assign(root: JsonRecord, rawPath: string, value: JsonValue): void {
  const parts = segments(rendererPath(rawPath));
  if (!parts.length) throw new Error('Android model bridge returned an invalid root scalar.');
  let current: JsonRecord | JsonValue[] = root;
  parts.forEach((part, position) => {
    const last = position === parts.length - 1;
    const next = parts[position + 1];
    if (part.key !== undefined) {
      if (Array.isArray(current)) throw new Error('Android model bridge returned a malformed object path.');
      if (last) {
        current[part.key] = value;
        return;
      }
      let child = current[part.key];
      if (child === undefined) {
        child = next?.index !== undefined ? [] : {};
        current[part.key] = child;
      }
      if (typeof child !== 'object' || child === null) {
        throw new Error('Android model bridge returned a conflicting field path.');
      }
      current = child as JsonRecord | JsonValue[];
      return;
    }
    if (part.index === undefined || !Array.isArray(current)) {
      throw new Error('Android model bridge returned a malformed array path.');
    }
    if (last) {
      current[part.index] = value;
      return;
    }
    let child = current[part.index];
    if (child === undefined) {
      child = next?.index !== undefined ? [] : {};
      current[part.index] = child;
    }
    if (typeof child !== 'object' || child === null) {
      throw new Error('Android model bridge returned a conflicting array path.');
    }
    current = child as JsonRecord | JsonValue[];
  });
}

function inflate(records: readonly string[]): JsonValue {
  if (records.length % 3 !== 0) throw new Error('Android model bridge returned malformed records.');
  const root: JsonRecord = {};
  for (let index = 0; index < records.length; index += 3) {
    const path = records[index]!;
    const kind = records[index + 1]!;
    const raw = records[index + 2]!;
    let value: JsonValue;
    if (kind === 'array') {
      const length = Number(raw);
      if (!Number.isInteger(length) || length < 0 || length > 262_144) {
        throw new Error('Android model bridge returned an invalid array length.');
      }
      value = new Array<JsonValue>(length);
    } else if (kind === 'string') value = raw;
    else if (kind === 'number') {
      const number = Number(raw);
      if (!Number.isFinite(number)) throw new Error('Android model bridge returned a non-finite numeric record.');
      value = number;
    } else if (kind === 'boolean' && (raw === 'true' || raw === 'false')) value = raw === 'true';
    else throw new Error('Android model bridge returned an unsupported record kind.');
    assign(root, path, value);
  }
  return root;
}

function requestRecords(input: JsonObject): string {
  const records = flatten(input);
  const text = JSON.stringify(records);
  if (new TextEncoder().encode(text).length > MAX_REQUEST_BYTES) {
    throw new Error('Android model workbench request exceeds its transport limit.');
  }
  return text;
}

function parseReply(reply: NativeRecordsReply): ModelWorkbenchReply {
  if (!Array.isArray(reply.records) || reply.records.some(item => typeof item !== 'string')) {
    throw new Error('Android model workbench returned an invalid native envelope.');
  }
  const bytes = new TextEncoder().encode(JSON.stringify(reply.records)).length;
  if (bytes > MAX_REPLY_BYTES || reply.records.length < 1) {
    throw new Error('Android model workbench reply exceeds its transport limit.');
  }
  if (reply.records[0] === 'ok-null' && reply.records.length === 1) {
    return { version: MODEL_WORKBENCH_CONVENTION, ok: true, value: null };
  }
  if (reply.records[0] === 'ok') {
    return {
      version: MODEL_WORKBENCH_CONVENTION,
      ok: true,
      value: inflate(reply.records.slice(1)),
    };
  }
  if (reply.records[0] === 'error' && (reply.records.length === 3 || reply.records.length === 5)) {
    const code = Number(reply.records[1]);
    const reason = reply.records[2]!;
    if (!Number.isInteger(code) || code < 1 || code > 16 || !reason) {
      throw new Error('Android model workbench returned an invalid error envelope.');
    }
    const field = reply.records[4];
    const validationCode = reply.records[3];
    return {
      version: MODEL_WORKBENCH_CONVENTION,
      ok: false,
      error: {
        code: code as Code,
        reason,
        ...(reply.records.length === 5 && validationCode
          ? { validation: { version: 'MPMC/model/validation-detail/v1', code: validationCode, ...(field ? { field } : {}) } }
          : {}),
      },
    };
  }
  throw new Error('Android model workbench returned an unsupported native envelope.');
}

async function transport(operation: Promise<NativeRecordsReply>): Promise<ModelWorkbenchReply> {
  try {
    return parseReply(await operation);
  } catch {
    return {
      version: MODEL_WORKBENCH_CONVENTION,
      ok: false,
      error: { code: Code.Unavailable, reason: 'workbench.failed' },
    };
  }
}

export function createAndroidModelWorkbenchBridge(): ModelWorkbenchBridge {
  return Object.freeze({
    convention: MODEL_WORKBENCH_CONVENTION,
    apply(requestId: string, input: JsonObject) {
      return transport(plugin.modelApply({ requestId, recordsJson: requestRecords(input) }));
    },
    solve(requestId: string, input: JsonObject) {
      return transport(plugin.modelSolve({ requestId, recordsJson: requestRecords(input) }));
    },
    release(requestId: string) {
      return transport(plugin.modelRelease({ requestId }));
    },
    cancel(requestId: string) {
      void plugin.modelCancel({ requestId }).catch(() => {});
    },
  });
}
