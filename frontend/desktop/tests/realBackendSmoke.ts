import { resolve } from 'node:path';

import { app, session, type BrowserWindow } from 'electron';

import type { DesktopBridgeReply } from '../../src/api/desktopBridgeContract';
import type { PtFlashRequest } from '../../src/domain/flash';
import type {
  PtCapabilityDiscovery,
  PtFlashResponse,
} from '../../src/domain/flash';
import { registerPtDesktopIpc } from '../desktopIpc';
import { PtHostSession } from '../hostSession';
import { PtDesktopGateway } from '../ptGateway';
import { createPtDesktopWindow, restrictPtDesktopSession } from '../window';

const expectedInventories = new Map<string, readonly string[]>([
  [
    'pr76.methane-ethane-propane.literature-r1',
    ['methane', 'ethane', 'propane'],
  ],
  [
    'sw92.carbon-dioxide-water.freshwater.literature-r1',
    ['carbon-dioxide', 'water'],
  ],
  [
    'cpa.methanol-water-333.15k.cr1.literature-r1',
    ['METHANOL', 'WATER'],
  ],
]);

const requests: readonly PtFlashRequest[] = [
  {
    configuredBackendId: 'pr76.methane-ethane-propane.literature-r1',
    pressurePa: 1e6,
    temperatureK: 350,
    feed: [
      { componentId: 'methane', moleFraction: 0.8 },
      { componentId: 'ethane', moleFraction: 0.1 },
      { componentId: 'propane', moleFraction: 0.1 },
    ],
  },
  {
    configuredBackendId: 'sw92.carbon-dioxide-water.freshwater.literature-r1',
    pressurePa: 3e6,
    temperatureK: 340,
    feed: [
      { componentId: 'carbon-dioxide', moleFraction: 0.7 },
      { componentId: 'water', moleFraction: 0.3 },
    ],
  },
  {
    configuredBackendId: 'cpa.methanol-water-333.15k.cr1.literature-r1',
    pressurePa: 48_852,
    temperatureK: 333.15,
    feed: [
      { componentId: 'METHANOL', moleFraction: 0.5 },
      { componentId: 'WATER', moleFraction: 0.5 },
    ],
  },
];

function requireCondition(condition: boolean, message: string): asserts condition {
  if (!condition) {
    throw new Error(message);
  }
}

function requireBridgeValue<T>(reply: DesktopBridgeReply<T>, operation: string): T {
  if (!reply.ok) {
    throw new Error(
      `${operation} failed at the desktop transport boundary: ${reply.error.code} ${reply.error.message}`,
    );
  }
  return reply.value;
}

async function waitForReactDiscovery(window: BrowserWindow): Promise<void> {
  const deadline = Date.now() + 20_000;
  while (Date.now() < deadline) {
    const ready = await window.webContents.executeJavaScript(
      `(() => {
        const shell = document.querySelector('[data-desktop-product-shell="true"]');
        const pt = document.querySelector('[data-product-mode="pt"]');
        const state = document.querySelector('.backend-state');
        return shell !== null && pt?.getAttribute('aria-pressed') === 'true' &&
          document.querySelector('[data-expert-session]') === null &&
          state?.getAttribute('data-configured') === 'true' &&
          state.textContent?.includes('3 backends discovered') === true;
      })()`,
      true,
    );
    if (ready === true) {
      console.info('DESKTOP_PRODUCT_DEFAULT_PT_OK shell=true expert_session=false');
      return;
    }
    await new Promise((resolveDelay) => setTimeout(resolveDelay, 100));
  }
  throw new Error('The embedded React renderer did not publish the default PT desktop shell and three-backend discovery.');
}

async function run(): Promise<void> {
  const rawHost = process.env.MPMC_PT_DESKTOP_HOST_BINARY;
  requireCondition(rawHost !== undefined, 'MPMC_PT_DESKTOP_HOST_BINARY is required');
  const hostLog: string[] = [];
  const host = new PtHostSession(resolve(rawHost), (line) => hostLog.push(line));
  const gateway = new PtDesktopGateway(host);
  restrictPtDesktopSession(session.defaultSession);
  const window = createPtDesktopWindow({
    preloadPath: resolve(process.cwd(), 'desktop', 'preload.cjs'),
    showWhenReady: false,
  });
  const removeIpc = registerPtDesktopIpc(gateway, () => window.webContents);
  try {
    await window.loadFile(resolve(process.cwd(), 'dist', 'index.html'));
    await waitForReactDiscovery(window);
    const discoveryReply = (await window.webContents.executeJavaScript(
      `globalThis.mpmcPtDesktop.discoverPtCapabilities('desktop-smoke-discovery')`,
      true,
    )) as DesktopBridgeReply<PtCapabilityDiscovery>;
    const discovery = requireBridgeValue(discoveryReply, 'desktop discovery');
    requireCondition(discovery.backends.length === 3, 'desktop discovery did not return three backends');
    for (const descriptor of discovery.backends) {
      const expected = expectedInventories.get(descriptor.configuredBackendId);
      requireCondition(expected !== undefined, 'desktop discovery returned an unknown backend');
      requireCondition(
        descriptor.componentInventory.components.map((component) => component.componentId).join('\0') ===
          expected.join('\0'),
        `desktop inventory changed for ${descriptor.configuredBackendId}`,
      );
    }

    for (const [index, request] of requests.entries()) {
      const solveReply = (await window.webContents.executeJavaScript(
        `globalThis.mpmcPtDesktop.solvePtFlash(
          'desktop-smoke-solve-${index}',
          ${JSON.stringify(request)}
        )`,
        true,
      )) as DesktopBridgeReply<PtFlashResponse>;
      const response = requireBridgeValue(
        solveReply,
        `desktop solve ${request.configuredBackendId}`,
      );
      requireCondition(response.kind === 'result', `desktop solve returned a service error for ${request.configuredBackendId}`);
      requireCondition(
        response.result.provenance.backend.configuredBackendId === request.configuredBackendId,
        `desktop solve provenance changed for ${request.configuredBackendId}`,
      );
      console.info(
        `DESKTOP_REAL_SOLVE_OK backend=${request.configuredBackendId} outcome=${response.result.outcome}`,
      );
    }
  } finally {
    removeIpc();
    if (!window.isDestroyed()) {
      window.destroy();
    }
    await gateway.stop();
  }
  requireCondition(
    hostLog.some((line) => line.includes('"event":"pt_process_stopped"')),
    'desktop host did not report graceful stdin shutdown',
  );
  console.info('DESKTOP_VERTICAL_SLICE_OK');
}

void app.whenReady().then(async () => {
  try {
    await run();
    app.exit(0);
  } catch (cause) {
    console.error(cause instanceof Error ? cause.stack ?? cause.message : String(cause));
    app.exit(1);
  }
});
