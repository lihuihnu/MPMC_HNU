declare module '*.css';

interface ImportMetaEnv {
  readonly VITE_MPMC_ANDROID_PRODUCT_SHELL_SMOKE?: string;
}

interface ImportMeta {
  readonly env: ImportMetaEnv;
}
