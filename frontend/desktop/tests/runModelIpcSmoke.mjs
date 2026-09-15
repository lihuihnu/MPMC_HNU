import { readFileSync, writeFileSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import electron from 'electron';

// Electron may exit cleanly before an async test finishes; exit status alone is insufficient.
// The product preload has no low-level model handle API. This harness temporarily substitutes
// the dedicated test preload for its child process, then restores the product file exactly.
const productPreload = 'desktop/preload.cjs';
const originalPreload = readFileSync(productPreload);
const debugPreload = readFileSync('desktop/modelDebugPreload.cjs');
let child;
try {
  writeFileSync(productPreload, debugPreload);
  child = spawnSync(electron, ['desktop-model-test-dist/model-ipc-smoke.mjs'], {
    encoding: 'utf8', timeout: 200_000, maxBuffer: 16 * 1024 * 1024,
  });
} finally {
  writeFileSync(productPreload, originalPreload);
}
process.stdout.write(child?.stdout ?? '');
process.stderr.write(child?.stderr ?? '');
const complete = ['MODEL_SOLVE_VALIDATION_PATH_OK', 'MODEL_NESTED_VALIDATION_PATH_OK', 'MODEL_VALIDATION_DETAIL_OK', 'MODEL_RENDERER_CLIENT_OK', 'MODEL_DESKTOP_IPC_OK', 'MODEL_DESKTOP_IPC_COMPLETE']
  .every(marker => (child?.stdout ?? '').split(/\r?\n/u).some(line => line === marker || line.startsWith(`${marker} `)));
if (!child || child.error || child.status !== 0 || !complete) {
  console.error('Model IPC smoke did not complete all assertions and host teardown.');
  process.exitCode = 1;
}
