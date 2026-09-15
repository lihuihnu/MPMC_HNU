import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { runInNewContext } from 'node:vm';

import { describe, expect, it, vi } from 'vitest';
import {
  MODEL_DESKTOP_CONVENTION,
  MODEL_DESKTOP_CHANNEL,
  MODEL_DESKTOP_CANCEL_CHANNEL,
  MODEL_DESKTOP_V2_CONVENTION,
  MODEL_DESKTOP_V2_CHANNEL,
  MODEL_DESKTOP_V2_CANCEL_CHANNEL,
  type ModelDesktopBridge,
} from '../src/api/modelDesktopContract';
import {
  MODEL_WORKBENCH_CANCEL_CHANNEL,
  MODEL_WORKBENCH_CHANNEL,
  MODEL_WORKBENCH_CONVENTION,
  type ModelWorkbenchBridge,
} from '../src/api/modelWorkbenchContract';
import type { PtDesktopBridge } from '../src/api/desktopBridgeContract';

function source(relativePath: string): string {
  return readFileSync(resolve(process.cwd(), relativePath), 'utf8');
}

function evaluatePreload(debug = false) {
  const preload = source('desktop/preload.cjs');
  const exposed: Record<string, unknown> = {};
  const invoke = vi.fn(async () => null);
  const send = vi.fn();
  runInNewContext(preload, {
    ...(debug ? { process: { env: { MPMC_MODEL_DESKTOP_DEBUG_BRIDGE: '1' } } } : {}),
    require(name: string) {
      expect(name).toBe('electron');
      return {
        contextBridge: {
          exposeInMainWorld(key: string, value: unknown) { exposed[key] = value; },
        },
        ipcRenderer: { invoke, send },
      };
    },
  });
  return { preload, exposed, invoke, send };
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
      'desktop/modelWorkbenchIpc.ts',
      'desktop/modelWorkbenchSession.ts',
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

  it('exposes only PT plus the handle-free editable-model workbench by default', () => {
    const window = source('desktop/window.ts');
    expect(window).toContain('contextIsolation: true');
    expect(window).toContain('nodeIntegration: false');
    expect(window).toContain('sandbox: true');
    expect(window).toContain('webSecurity: true');
    expect(window).toContain("setWindowOpenHandler(() => ({ action: 'deny' }))");
    expect(window).toContain("on('will-navigate'");
    expect(window).toContain("on('will-attach-webview'");

    const { preload, exposed, invoke, send } = evaluatePreload(false);
    expect(Object.keys(exposed)).toEqual(['mpmcPtDesktop', 'mpmcModelWorkbench']);

    const pt = exposed.mpmcPtDesktop as PtDesktopBridge;
    expect(Object.keys(pt)).toEqual(['convention', 'discoverPtCapabilities', 'solvePtFlash', 'cancel']);
    expect(Object.isFrozen(pt)).toBe(true);

    const workbench = exposed.mpmcModelWorkbench as ModelWorkbenchBridge;
    expect(Object.isFrozen(workbench)).toBe(true);
    expect(Object.keys(workbench)).toEqual(['convention', 'apply', 'solve', 'release', 'cancel']);
    expect(workbench.convention).toBe(MODEL_WORKBENCH_CONVENTION);
    void workbench.apply('apply', { definition: {} });
    void workbench.solve('solve', { pressurePa: 1, temperatureK: 1, feed: [1] });
    void workbench.release('release');
    workbench.cancel('cancel');
    expect(invoke.mock.calls).toEqual([
      [MODEL_WORKBENCH_CHANNEL, {
        version: MODEL_WORKBENCH_CONVENTION,
        requestId: 'apply', operation: 'apply', input: { definition: {} },
      }],
      [MODEL_WORKBENCH_CHANNEL, {
        version: MODEL_WORKBENCH_CONVENTION,
        requestId: 'solve', operation: 'solve',
        input: { pressurePa: 1, temperatureK: 1, feed: [1] },
      }],
      [MODEL_WORKBENCH_CHANNEL, {
        version: MODEL_WORKBENCH_CONVENTION,
        requestId: 'release', operation: 'release',
      }],
    ]);
    expect(send.mock.calls).toEqual([[MODEL_WORKBENCH_CANCEL_CHANNEL, 'cancel']]);
    expect(preload).not.toContain('sendSync');
  });

  it('keeps low-level model handles behind the dedicated regression flag', () => {
    const regular = evaluatePreload(false).exposed;
    expect(regular).not.toHaveProperty('mpmcModelDesktop');
    expect(regular).not.toHaveProperty('mpmcModelDesktopV2');

    const { exposed, invoke, send } = evaluatePreload(true);
    expect(Object.keys(exposed)).toEqual([
      'mpmcPtDesktop', 'mpmcModelWorkbench', 'mpmcModelDesktop', 'mpmcModelDesktopV2',
    ]);
    const legacy = exposed.mpmcModelDesktop as ModelDesktopBridge;
    const v2 = exposed.mpmcModelDesktopV2 as ModelDesktopBridge;
    expect(Object.keys(legacy)).toEqual([
      'convention', 'connect', 'reconnect', 'create', 'describe', 'solve', 'release', 'cancel',
    ]);
    expect(Object.isFrozen(legacy)).toBe(true);
    expect(legacy.convention).toBe(MODEL_DESKTOP_CONVENTION);
    expect(v2.convention).toBe(MODEL_DESKTOP_V2_CONVENTION);

    void legacy.create('c1', { presetId: 'test' });
    void legacy.solve('s1', 'local-token', { feed: [1] });
    legacy.cancel('x1');
    void v2.create('c2', { presetId: 'test' });
    void v2.solve('s2', 'local-token', { feed: [1] });
    v2.cancel('x2');
    expect(invoke.mock.calls.slice(-4)).toEqual([
      [MODEL_DESKTOP_CHANNEL, {
        version: MODEL_DESKTOP_CONVENTION, requestId: 'c1', operation: 'create', input: { presetId: 'test' },
      }],
      [MODEL_DESKTOP_CHANNEL, {
        version: MODEL_DESKTOP_CONVENTION, requestId: 's1', operation: 'solve', model: 'local-token', input: { feed: [1] },
      }],
      [MODEL_DESKTOP_V2_CHANNEL, {
        version: MODEL_DESKTOP_V2_CONVENTION, requestId: 'c2', operation: 'create', input: { presetId: 'test' },
      }],
      [MODEL_DESKTOP_V2_CHANNEL, {
        version: MODEL_DESKTOP_V2_CONVENTION, requestId: 's2', operation: 'solve', model: 'local-token', input: { feed: [1] },
      }],
    ]);
    expect(send.mock.calls.slice(-2)).toEqual([
      [MODEL_DESKTOP_CANCEL_CHANNEL, 'x1'],
      [MODEL_DESKTOP_V2_CANCEL_CHANNEL, 'x2'],
    ]);
  });

  it('starts only an ephemeral loopback child session with a private transport token', () => {
    const session = source('desktop/hostSession.ts');
    expect(session).toContain("randomBytes(32).toString('base64url')");
    expect(session).toContain("['--desktop-session-token-stdin']");
    expect(session).toContain('shell: false');
    expect(session).toContain("child.stdin.write(`${token}\\n`)");
    expect(session).not.toContain('bearerToken: process.env');
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
