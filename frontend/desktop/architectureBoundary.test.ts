import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';

import { describe, expect, it } from 'vitest';

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
    expect(preload.match(/ipcRenderer\.(invoke|send)\(/gu)).toHaveLength(3);
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
});
