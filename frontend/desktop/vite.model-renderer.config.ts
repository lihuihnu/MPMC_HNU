import { resolve } from 'node:path';
import { defineConfig } from 'vite';

export default defineConfig({
  build: {
    target: 'es2022',
    outDir: 'desktop-model-test-dist',
    emptyOutDir: false,
    minify: false,
    lib: {
      entry: resolve('desktop/tests/modelRendererSmoke.ts'),
      name: 'MpmcModelRendererSmoke',
      formats: ['iife'],
      fileName: () => 'model-renderer-smoke.js',
    },
  },
});
