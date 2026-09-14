import path from 'node:path';
import { fileURLToPath } from 'node:url';

const shellRoot = path.dirname(fileURLToPath(import.meta.url));

export default {
  root: shellRoot,
  base: './',
  esbuild: {
    jsx: 'automatic',
    jsxImportSource: 'react',
  },
  resolve: {
    alias: [
      {
        find: /^react$/u,
        replacement: path.join(shellRoot, 'node_modules/react/index.js'),
      },
      {
        find: /^react\/jsx-runtime$/u,
        replacement: path.join(shellRoot, 'node_modules/react/jsx-runtime.js'),
      },
      {
        find: /^react-dom\/client$/u,
        replacement: path.join(shellRoot, 'node_modules/react-dom/client.js'),
      },
    ],
    dedupe: ['react', 'react-dom'],
  },
  build: {
    outDir: path.join(shellRoot, 'dist'),
    emptyOutDir: true,
  },
};
