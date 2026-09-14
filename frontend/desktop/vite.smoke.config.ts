import { builtinModules } from 'node:module';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { defineConfig } from 'vite';

const here = dirname(fileURLToPath(import.meta.url));
const external = [
  'electron',
  ...builtinModules,
  ...builtinModules.map((name) => `node:${name}`),
];

export default defineConfig({
  build: {
    target: 'node24',
    ssr: resolve(here, 'tests/realBackendSmoke.ts'),
    outDir: resolve(here, '../desktop-test-dist'),
    emptyOutDir: true,
    minify: false,
    sourcemap: true,
    rollupOptions: {
      external,
      output: {
        entryFileNames: 'real-backend-smoke.mjs',
      },
    },
  },
  ssr: {
    noExternal: true,
  },
});
