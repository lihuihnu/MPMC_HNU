import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { runInNewContext } from 'node:vm';

import { describe, expect, it, vi } from 'vitest';
import { MODEL_DESKTOP_CONVENTION, MODEL_DESKTOP_CHANNEL, MODEL_DESKTOP_CANCEL_CHANNEL, type ModelDesktopBridge } from '../src/api/modelDesktopContract';
import type { PtDesktopBridge } from '../src/api/desktopBridgeContract';

function source(relativePath: string): string {
  return readFileSync(resolve(process.cwd(), relativePath), 'utf8');
}

describe('PT desktop architecture boundary', () => {
  it('keeps the Electron host, IPC, and native gateway model-neutral', () => {
    const files = [
      'desktop/main.ts',
      'desktop/window.ts',
      'desktop/hostSession.ts',
      'desktop/desktopIpc.ts',
      'desktop/desktopIpcPolicy.ts',
      'desktop/ptGateway.ts',
      'desktop/preload.cjs',
    ];
    for (const file of files) {
      const text = source(file).toLowerCase();
      expect(text, file).not.toMatch(/pr76|sw92|cpa/u);
      expect(text, file).not.toMatch(/modules\/(flash|thermodynamics)/u);
      expect(text, file).not.toMatch(/solve_(pr76|sw92|cpa)/u);
    }

    const gateway = source('desktop/ptGateway.ts');
    expect(gateway.match(/client\.solvePtFlash\(/gu)).toHaveLength(1);
    expect(gateway).toContain('toWireSolveRequest(request)');
    expect(gateway).toContain('mapSolveResponse(wire, expected, request)');
  });

  it('locks the renderer behind a narrow versioned preload bridge', () => {
    const window = source('desktop/window.ts');
    expect(window).toContain('contextIsolation: true');
    expect(window).toContain('nodeIntegration: false');
    expect(window).toContain('sandbox: true');
    expect(window).toContain('webSecurity: true');
    expect(window).toContain("setWindowOpenHandler(() => ({ action: 'deny' }))");
    expect(window).toContain("on('will-navigate'");
    expect(window).toContain("on('will-attach-webview'");

    const preload = source('desktop/preload.cjs');
    expect(preload).toContain("convention: 'MPMC/PT/desktop-bridge/v1'");
    expect([...preload.matchAll(/ipcRenderer\.(?:invoke|send)\(\s*'([^']+)'/gu)].map(match => match[1]))
      .toEqual(['mpmc:pt:discover:v1', 'mpmc:pt:solve:v1', 'mpmc:pt:cancel:v1',
        MODEL_DESKTOP_CHANNEL, MODEL_DESKTOP_CANCEL_CHANNEL]);
    const exposed: Record<string, unknown> = {};
    const invoke = vi.fn(async () => null); const send = vi.fn();
    runInNewContext(preload, { require(name: string) {
      expect(name).toBe('electron');
      return { contextBridge: { exposeInMainWorld(key: string, value: unknown) { exposed[key] = value; } },
        ipcRenderer: { invoke, send } };
    } });
    expect(Object.keys(exposed)).toEqual(['mpmcPtDesktop', 'mpmcModelDesktop']);
    const legacy = exposed.mpmcPtDesktop as PtDesktopBridge;
    expect(Object.keys(legacy)).toEqual(['convention', 'discoverPtCapabilities', 'solvePtFlash', 'cancel']);
    expect(Object.isFrozen(legacy)).toBe(true);
    const model = exposed.mpmcModelDesktop as ModelDesktopBridge;
    expect(Object.isFrozen(model)).toBe(true);
    expect(Object.keys(model)).toEqual(['convention', 'connect', 'reconnect', 'create', 'describe', 'solve', 'release', 'cancel']);
    expect(model.convention).toBe(MODEL_DESKTOP_CONVENTION);
    void model.connect('open'); void model.reconnect('again');
    void model.create('create', { presetId: 'explicit-test' });
    void model.describe('describe', 'local'); void model.solve('solve', 'local', { feed: [1] });
    void model.release('release', 'local'); model.cancel('cancel');
    expect(invoke.mock.calls).toEqual([
      [MODEL_DESKTOP_CHANNEL, { version: MODEL_DESKTOP_CONVENTION, operation: 'connect', requestId: 'open' }],
      [MODEL_DESKTOP_CHANNEL, { version: MODEL_DESKTOP_CONVENTION, operation: 'reconnect', requestId: 'again' }],
      [MODEL_DESKTOP_CHANNEL, { version: MODEL_DESKTOP_CONVENTION, operation: 'create', requestId: 'create', input: { presetId: 'explicit-test' } }],
      [MODEL_DESKTOP_CHANNEL, { version: MODEL_DESKTOP_CONVENTION, operation: 'describe', requestId: 'describe', model: 'local' }],
      [MODEL_DESKTOP_CHANNEL, { version: MODEL_DESKTOP_CONVENTION, operation: 'solve', requestId: 'solve', model: 'local', input: { feed: [1] } }],
      [MODEL_DESKTOP_CHANNEL, { version: MODEL_DESKTOP_CONVENTION, operation: 'release', requestId: 'release', model: 'local' }],
    ]);
    expect(send.mock.calls).toEqual([[MODEL_DESKTOP_CANCEL_CHANNEL, 'cancel']]);
    expect(preload).not.toContain('sendSync');
  });

  it('starts only an authenticated ephemeral loopback child session', () => {
    const session = source('desktop/hostSession.ts');
    expect(session).toContain("randomBytes(32).toString('base64url')");
    expect(session).toContain("['--desktop-session-token-stdin']");
    expect(session).toContain("shell: false");
    expect(session).toContain("child.stdin.write(`${token}\\n`)");
    expect(session).not.toContain("bearerToken: process.env");
  });

  it('keeps installed smoke generic and unreachable from the renderer API', () => {
    const smoke = source('desktop/installSmoke.ts').toLowerCase();
    expect(smoke).not.toMatch(/pr76|sw92|cpa/u);
    expect(smoke).not.toMatch(/modules\/(flash|thermodynamics)/u);
    expect(smoke).toContain('globalthis.mpmcptdesktop.solveptflash');
    expect(smoke).toContain('validateptflashrequest');

    const preload = source('desktop/preload.cjs');
    expect(preload).not.toContain('install-smoke');
    expect(preload).not.toContain('INSTALL_SMOKE');
  });
});
