import { packager } from '@electron/packager';
import {
  access,
  cp,
  mkdir,
  readFile,
  readdir,
  rm,
  stat,
  writeFile,
} from 'node:fs/promises';
import { constants } from 'node:fs';
import { basename, join, resolve } from 'node:path';

const ELECTRON_VERSION = '44.3.0';
const allowedPlatforms = new Set(['linux', 'win32', 'darwin']);
const allowedArchitectures = new Set(['x64', 'arm64']);
const expectedNativeProfiles = new Map([
  ['linux/x64', 'linux-x86_64-release'],
  ['win32/x64', 'windows-x86_64-release'],
  ['darwin/arm64', 'macos-armv8-release'],
]);

function argumentsFrom(argv) {
  const result = new Map();
  for (let index = 0; index < argv.length; ++index) {
    const option = argv[index];
    if (!option.startsWith('--') || index + 1 >= argv.length) {
      throw new Error(`Invalid desktop package argument: ${option}`);
    }
    if (result.has(option)) {
      throw new Error(`Duplicate desktop package argument: ${option}`);
    }
    result.set(option, argv[++index]);
  }
  return result;
}

function required(options, name) {
  const value = options.get(name);
  if (typeof value !== 'string' || value.length === 0) {
    throw new Error(`Missing desktop package argument: ${name}`);
  }
  return value;
}

async function findNamed(root, wanted) {
  const entries = await readdir(root, { withFileTypes: true });
  for (const entry of entries) {
    const path = join(root, entry.name);
    if (entry.name === wanted) {
      return path;
    }
    if (entry.isDirectory()) {
      const nested = await findNamed(path, wanted);
      if (nested !== null) {
        return nested;
      }
    }
  }
  return null;
}

async function main() {
  const options = argumentsFrom(process.argv.slice(2));
  const nativeStage = resolve(required(options, '--native-stage'));
  const output = resolve(required(options, '--out'));
  const platform = required(options, '--platform');
  const arch = required(options, '--arch');
  const sourceRevision = required(options, '--source-revision');
  if (!allowedPlatforms.has(platform) || !allowedArchitectures.has(arch)) {
    throw new Error('Unsupported desktop preview platform or architecture.');
  }
  if (!/^[A-Za-z0-9][A-Za-z0-9._+-]{6,127}$/u.test(sourceRevision)) {
    throw new Error('Desktop preview source revision is not a stable identifier.');
  }

  const executable = platform === 'win32'
    ? 'mpmc_pt_service_host.exe'
    : 'mpmc_pt_service_host';
  const nativeHost = join(nativeStage, 'bin', executable);
  await access(nativeHost, constants.R_OK);
  if (platform !== 'win32') {
    const mode = (await stat(nativeHost)).mode;
    if ((mode & 0o111) === 0) {
      throw new Error('The staged PT host is not executable.');
    }
  }
  const stagingManifestPath = join(
    nativeStage,
    'share',
    'mpmc-pt',
    'product-staging-manifest.json',
  );
  const stagingManifest = JSON.parse(await readFile(stagingManifestPath, 'utf8'));
  if (stagingManifest.convention !== 'MPMC/PT/product-staging/v1') {
    throw new Error('The native PT staging convention changed.');
  }
  if (stagingManifest.product?.build_revision !== sourceRevision) {
    throw new Error('Desktop and native staging source revisions do not match.');
  }
  const expectedProfile = expectedNativeProfiles.get(`${platform}/${arch}`);
  if (
    expectedProfile === undefined ||
    stagingManifest.platform?.dependency_profile !== expectedProfile
  ) {
    throw new Error('The native staging profile does not match the desktop target.');
  }
  if (
    stagingManifest.entry?.transport !== 'native-grpc-mtls' ||
    stagingManifest.entry?.requires_explicit_tls_paths !== true
  ) {
    throw new Error('The production native host security contract changed.');
  }

  const frontendRoot = resolve(import.meta.dirname, '..');
  const buildRoot = join(frontendRoot, 'build', 'desktop-package');
  const appRoot = join(buildRoot, 'app');
  const copiedNative = join(buildRoot, 'desktop-native');
  const previewManifestPath = join(buildRoot, 'desktop-preview-manifest.json');
  await rm(buildRoot, { recursive: true, force: true });
  await rm(output, { recursive: true, force: true });
  await mkdir(appRoot, { recursive: true });

  await cp(join(frontendRoot, 'dist'), join(appRoot, 'renderer'), { recursive: true });
  await cp(join(frontendRoot, 'desktop-dist', 'main.mjs'), join(appRoot, 'main.mjs'));
  await cp(join(frontendRoot, 'desktop', 'preload.cjs'), join(appRoot, 'preload.cjs'));
  await cp(nativeStage, copiedNative, { recursive: true, preserveTimestamps: true });
  await writeFile(
    join(appRoot, 'package.json'),
    `${JSON.stringify({
      name: 'mpmc-pt-desktop-preview',
      productName: 'MPMC PT Desktop Preview',
      version: '0.1.0',
      author: 'MPMC_HNU contributors',
      private: true,
      type: 'module',
      main: 'main.mjs',
    }, null, 2)}\n`,
    'utf8',
  );
  const previewManifest = {
    convention: 'MPMC/PT/desktop-preview/v1',
    source_revision: sourceRevision,
    native_staging_convention: stagingManifest.convention,
    native_dependencies: stagingManifest.dependencies,
    electron_version: ELECTRON_VERSION,
    platform,
    architecture: arch,
    release_eligible: false,
    signature: { status: 'unsigned-preview' },
    production_identity: null,
    desktop_session: {
      transport: 'native-grpc-loopback-bearer',
      listener: '127.0.0.1:0',
      token_delivery: 'child-stdin',
      shutdown: 'child-stdin-or-parent-exit',
    },
  };
  await writeFile(previewManifestPath, `${JSON.stringify(previewManifest, null, 2)}\n`, 'utf8');

  const packagePaths = await packager({
    dir: appRoot,
    out: output,
    name: 'MPMC-PT-Desktop-Preview',
    executableName: 'MPMC-PT-Desktop-Preview',
    platform,
    arch,
    electronVersion: ELECTRON_VERSION,
    appVersion: '0.1.0',
    buildVersion: '0.1.0',
    appBundleId: 'invalid.mpmc-hnu.pt-desktop-preview',
    asar: true,
    overwrite: true,
    prune: false,
    extraResource: [
      copiedNative,
      previewManifestPath,
      join(frontendRoot, 'desktop', 'desktop-preview-notice.txt'),
    ],
  });
  if (packagePaths.length !== 1) {
    throw new Error('Electron packager returned an unexpected output set.');
  }
  const packageRoot = packagePaths[0];
  const packagedHost = await findNamed(packageRoot, executable);
  const packagedManifest = await findNamed(packageRoot, basename(previewManifestPath));
  const packagedNotice = await findNamed(packageRoot, 'desktop-preview-notice.txt');
  const appArchive = await findNamed(packageRoot, 'app.asar');
  if (
    packagedHost === null ||
    packagedManifest === null ||
    packagedNotice === null ||
    appArchive === null
  ) {
    throw new Error('Packaged desktop preview omitted a required runtime resource.');
  }
  console.info(`DESKTOP_PACKAGE_OK path=${packageRoot}`);
}

await main();
