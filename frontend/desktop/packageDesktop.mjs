import { createHash } from 'node:crypto';
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

import { packager } from '@electron/packager';

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

function requireString(value, label, pattern = null) {
  if (typeof value !== 'string' || value.length === 0) {
    throw new Error(`Windows desktop release identity is missing ${label}.`);
  }
  if (pattern !== null && !pattern.test(value)) {
    throw new Error(`Windows desktop release identity has invalid ${label}.`);
  }
  return value;
}

async function readReleaseIdentity(path) {
  const resolved = resolve(path);
  const text = await readFile(resolved, 'utf8');
  const identity = JSON.parse(text);
  if (
    identity?.convention !== 'MPMC/PT/windows-desktop-release-identity/v1' ||
    identity?.product_id !== 'mpmc-pt-desktop' ||
    identity?.platform !== 'windows-x86_64' ||
    identity?.release?.channel !== 'rc'
  ) {
    throw new Error('Windows desktop release identity convention changed.');
  }
  const msiVersion = requireString(
    identity.release?.msi_product_version,
    'MSI product version',
    /^[0-9]+\.[0-9]+\.[0-9]+$/u,
  );
  const displayVersion = requireString(
    identity.release?.display_version,
    'display version',
    /^[0-9]+\.[0-9]+\.[0-9]+-rc\.[1-9][0-9]*$/u,
  );
  const releaseExecutable = requireString(
    identity.application?.release_executable,
    'release executable',
    /^[A-Za-z0-9][A-Za-z0-9._+-]*\.exe$/u,
  );
  const productName = requireString(identity.product_name, 'product name');
  const bundleId = requireString(
    identity.application?.electron_bundle_id,
    'Electron bundle id',
    /^[A-Za-z0-9][A-Za-z0-9.-]+$/u,
  );
  const productGuid = requireString(
    identity.release?.product_guid,
    'RC product GUID',
    /^[0-9A-F]{8}(?:-[0-9A-F]{4}){3}-[0-9A-F]{12}$/u,
  );
  const upgradeGuid = requireString(
    identity.upgrade_guid,
    'upgrade GUID',
    /^[0-9A-F]{8}(?:-[0-9A-F]{4}){3}-[0-9A-F]{12}$/u,
  );
  return {
    path: resolved,
    sha256: createHash('sha256').update(text).digest('hex'),
    identity,
    msiVersion,
    displayVersion,
    releaseExecutable,
    executableStem: releaseExecutable.replace(/\.exe$/u, ''),
    productName,
    bundleId,
    productGuid,
    upgradeGuid,
  };
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

  const releaseIdentityArgument = options.get('--release-identity');
  const releaseIdentity = releaseIdentityArgument === undefined
    ? null
    : await readReleaseIdentity(releaseIdentityArgument);
  if (releaseIdentity !== null && (platform !== 'win32' || arch !== 'x64')) {
    throw new Error('The fixed Windows release identity may only package win32/x64.');
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
  const desktopManifestPath = join(buildRoot, 'desktop-preview-manifest.json');
  await rm(buildRoot, { recursive: true, force: true });
  await rm(output, { recursive: true, force: true });
  await mkdir(appRoot, { recursive: true });

  await cp(join(frontendRoot, 'dist'), join(appRoot, 'renderer'), { recursive: true });
  await cp(join(frontendRoot, 'desktop-dist', 'main.mjs'), join(appRoot, 'main.mjs'));
  await cp(join(frontendRoot, 'desktop', 'preload.cjs'), join(appRoot, 'preload.cjs'));
  await cp(nativeStage, copiedNative, { recursive: true, preserveTimestamps: true });

  const appVersion = releaseIdentity?.msiVersion ?? '0.1.0';
  const productName = releaseIdentity?.productName ?? 'MPMC PT Desktop Preview';
  const packageName = releaseIdentity === null
    ? 'mpmc-pt-desktop-preview'
    : 'mpmc-pt-desktop';
  const electronExecutable = releaseIdentity?.executableStem ?? 'MPMC-PT-Desktop-Preview';
  const electronPackageName = releaseIdentity === null
    ? 'MPMC-PT-Desktop-Preview'
    : 'MPMC-PT-Desktop';
  const bundleId = releaseIdentity?.bundleId ?? 'invalid.mpmc-hnu.pt-desktop-preview';
  const noticePath = releaseIdentity === null
    ? join(frontendRoot, 'desktop', 'desktop-preview-notice.txt')
    : join(frontendRoot, 'desktop', 'desktop-release-candidate-notice.txt');

  await writeFile(
    join(appRoot, 'package.json'),
    `${JSON.stringify({
      name: packageName,
      productName,
      version: appVersion,
      author: 'MPMC_HNU contributors',
      private: true,
      type: 'module',
      main: 'main.mjs',
    }, null, 2)}\n`,
    'utf8',
  );

  const baseManifest = {
    source_revision: sourceRevision,
    native_staging_convention: stagingManifest.convention,
    native_dependencies: stagingManifest.dependencies,
    electron_version: ELECTRON_VERSION,
    platform,
    architecture: arch,
    release_eligible: false,
    production_identity: null,
    desktop_session: {
      transport: 'native-grpc-loopback-bearer',
      listener: '127.0.0.1:0',
      token_delivery: 'child-stdin',
      shutdown: 'child-stdin-or-parent-exit',
    },
  };
  const desktopManifest = releaseIdentity === null
    ? {
        convention: 'MPMC/PT/desktop-preview/v1',
        ...baseManifest,
        signature: { status: 'unsigned-preview' },
      }
    : {
        convention: 'MPMC/PT/desktop-release-candidate-payload/v1',
        ...baseManifest,
        release_candidate: true,
        release_identity: {
          convention: releaseIdentity.identity.convention,
          product_id: releaseIdentity.identity.product_id,
          display_version: releaseIdentity.displayVersion,
          msi_product_version: releaseIdentity.msiVersion,
          product_guid: releaseIdentity.productGuid,
          upgrade_guid: releaseIdentity.upgradeGuid,
          identity_sha256: releaseIdentity.sha256,
        },
        signature: {
          status: 'awaiting-authenticode',
          required: 'desktop-executable-native-host-and-msi-rfc3161',
        },
      };
  await writeFile(desktopManifestPath, `${JSON.stringify(desktopManifest, null, 2)}\n`, 'utf8');

  const packagePaths = await packager({
    dir: appRoot,
    out: output,
    name: electronPackageName,
    executableName: electronExecutable,
    platform,
    arch,
    electronVersion: ELECTRON_VERSION,
    appVersion,
    buildVersion: appVersion,
    appBundleId: bundleId,
    asar: true,
    overwrite: true,
    prune: false,
    extraResource: [
      copiedNative,
      desktopManifestPath,
      noticePath,
    ],
  });
  if (packagePaths.length !== 1) {
    throw new Error('Electron packager returned an unexpected output set.');
  }
  const packageRoot = packagePaths[0];
  const packagedDesktopExecutable = releaseIdentity?.releaseExecutable ?? 'MPMC-PT-Desktop-Preview.exe';
  const packagedHost = await findNamed(packageRoot, executable);
  const packagedDesktop = await findNamed(packageRoot, packagedDesktopExecutable);
  const packagedManifest = await findNamed(packageRoot, basename(desktopManifestPath));
  const packagedNotice = await findNamed(packageRoot, basename(noticePath));
  const appArchive = await findNamed(packageRoot, 'app.asar');
  if (
    packagedDesktop === null ||
    packagedHost === null ||
    packagedManifest === null ||
    packagedNotice === null ||
    appArchive === null
  ) {
    throw new Error('Packaged desktop payload omitted a required runtime resource.');
  }
  console.info(
    `${releaseIdentity === null ? 'DESKTOP_PACKAGE_OK' : 'DESKTOP_RC_PAYLOAD_OK'} path=${packageRoot}`,
  );
}

await main();
