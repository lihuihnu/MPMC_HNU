import path from 'node:path';
import { fileURLToPath } from 'node:url';

const shellRoot = path.dirname(fileURLToPath(import.meta.url));
const repositoryRoot = path.resolve(shellRoot, '../../..');
const frontendRoot = path.join(repositoryRoot, 'frontend');

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
        replacement: path.join(frontendRoot, 'node_modules/react/index.js'),
      },
      {
        find: /^react\/jsx-runtime$/u,
        replacement: path.join(frontendRoot, 'node_modules/react/jsx-runtime.js'),
      },
      {
        find: /^react-dom\/client$/u,
        replacement: path.join(frontendRoot, 'node_modules/react-dom/client.js'),
      },
    ],
  },
  build: {
    outDir: path.join(shellRoot, 'dist'),
    emptyOutDir: true,
  },
};
