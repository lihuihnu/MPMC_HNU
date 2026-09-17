import { readFile, writeFile } from 'node:fs/promises';
import { isAbsolute } from 'node:path';

import type { BrowserWindow } from 'electron';

import type { DesktopBridgeReply } from '../src/api/desktopBridgeContract';
import {
  type PtBackendDescriptor,
  type PtCapabilityDiscovery,
  type PtFlashRequest,
  type PtFlashResponse,
  validatePtFlashRequest,
} from '../src/domain/flash';

export const PT_DESKTOP_INSTALL_SMOKE_PLAN_CONVENTION =
  'MPMC/PT/windows-desktop-install-smoke-plan/v1' as const;
export const PT_DESKTOP_INSTALL_SMOKE_RESULT_CONVENTION =
  'MPMC/PT/windows-desktop-install-smoke-result/v1' as const;

export interface PtDesktopInstallSmokePlan {
  convention: typeof PT_DESKTOP_INSTALL_SMOKE_PLAN_CONVENTION;
  expectedBackendCount: number;
  requests: readonly PtFlashRequest[];
}

function requireCondition(condition: boolean, message: string): asserts condition {
  if (!condition) {
    throw new Error(message);
  }
}

function requireBridgeValue<T>(
  reply: DesktopBridgeReply<T>,
  operation: string,
): T {
  if (!reply.ok) {
    throw new Error(
      `${operation} failed at the desktop transport boundary: ${reply.error.code} ${reply.error.message}`,
    );
  }
  return reply.value;
}

function requestFrom(value: unknown, index: number): PtFlashRequest {
  requireCondition(
    typeof value === 'object' && value !== null,
    `Installed smoke request ${index} is not an object.`,
  );
  const request = value as Record<string, unknown>;
  requireCondition(
    typeof request.configuredBackendId === 'string' &&
      request.configuredBackendId.length > 0,
    `Installed smoke request ${index} has no backend id.`,
  );
  requireCondition(
    typeof request.pressurePa === 'number' &&
      Number.isFinite(request.pressurePa),
    `Installed smoke request ${index} has invalid pressure.`,
  );
  requireCondition(
    typeof request.temperatureK === 'number' &&
      Number.isFinite(request.temperatureK),
    `Installed smoke request ${index} has invalid temperature.`,
  );
  requireCondition(
    Array.isArray(request.feed),
    `Installed smoke request ${index} has no feed array.`,
  );
  const feed = request.feed.map((entry, componentIndex) => {
    requireCondition(
      typeof entry === 'object' && entry !== null,
      `Installed smoke request ${index} feed ${componentIndex} is invalid.`,
    );
    const component = entry as Record<string, unknown>;
    requireCondition(
      typeof component.componentId === 'string' &&
        component.componentId.length > 0,
      `Installed smoke request ${index} feed ${componentIndex} has no component id.`,
    );
    requireCondition(
      typeof component.moleFraction === 'number' &&
        Number.isFinite(component.moleFraction),
      `Installed smoke request ${index} feed ${componentIndex} has invalid fraction.`,
    );
    return {
      componentId: component.componentId,
      moleFraction: component.moleFraction,
    };
  });
  return {
    configuredBackendId: request.configuredBackendId,
    pressurePa: request.pressurePa,
    temperatureK: request.temperatureK,
    feed,
  };
}

export function parseInstalledSmokePlan(value: unknown): PtDesktopInstallSmokePlan {
  requireCondition(
    typeof value === 'object' && value !== null,
    'Installed smoke plan is not an object.',
  );
  const plan = value as Record<string, unknown>;
  requireCondition(
    plan.convention === PT_DESKTOP_INSTALL_SMOKE_PLAN_CONVENTION,
    'Installed smoke plan convention changed.',
  );
  requireCondition(
    typeof plan.expectedBackendCount === 'number' &&
      Number.isInteger(plan.expectedBackendCount) &&
      plan.expectedBackendCount > 0,
    'Installed smoke expected backend count is invalid.',
  );
  requireCondition(
    Array.isArray(plan.requests) && plan.requests.length > 0,
    'Installed smoke plan has no solve requests.',
  );
  const requests = plan.requests.map((request, index) =>
    requestFrom(request, index),
  );
  return {
    convention: PT_DESKTOP_INSTALL_SMOKE_PLAN_CONVENTION,
    expectedBackendCount: plan.expectedBackendCount,
    requests,
  };
}

export function installedSmokeRequested(argv: readonly string[]): boolean {
  return argv.includes('--mpmc-install-smoke');
}

async function waitForReactProduct(window: BrowserWindow): Promise<void> {
  const deadline = Date.now() + 20_000;
  while (Date.now() < deadline) {
    const ready = await window.webContents.executeJavaScript(
      `(() => {
        const shell = document.querySelector('[data-desktop-product-shell="true"]');
        const workspace = document.querySelector('[data-expert-workbench="ready"]');
        return shell !== null && workspace !== null &&
          shell.textContent?.includes('Classic PR') === true;
      })()`,
      true,
    );
    if (ready === true) {
      return;
    }
    await new Promise((resolveDelay) => setTimeout(resolveDelay, 100));
  }
  throw new Error('The installed React renderer did not publish the Classic PR product shell.');
}

function backendFor(
  discovery: PtCapabilityDiscovery,
  request: PtFlashRequest,
): PtBackendDescriptor {
  const backend = discovery.backends.find(
    (candidate) =>
      candidate.configuredBackendId === request.configuredBackendId,
  );
  requireCondition(
    backend !== undefined,
    `Installed smoke requested an undiscovered backend: ${request.configuredBackendId}.`,
  );
  return backend;
}

export async function runInstalledDesktopSmoke(
  window: BrowserWindow,
  planPath: string,
  resultPath: string,
): Promise<void> {
  requireCondition(isAbsolute(planPath), 'Installed smoke plan path must be absolute.');
  requireCondition(
    isAbsolute(resultPath),
    'Installed smoke result path must be absolute.',
  );

  const plan = parseInstalledSmokePlan(
    JSON.parse(await readFile(planPath, 'utf8')) as unknown,
  );
  await waitForReactProduct(window);

  const discoveryReply = (await window.webContents.executeJavaScript(
    `globalThis.mpmcPtDesktop.discoverPtCapabilities('installed-desktop-smoke-discovery')`,
    true,
  )) as DesktopBridgeReply<PtCapabilityDiscovery>;
  const discovery = requireBridgeValue(
    discoveryReply,
    'installed desktop discovery',
  );
  requireCondition(
    discovery.backends.length === plan.expectedBackendCount,
    `Installed desktop discovered ${discovery.backends.length} backends; expected ${plan.expectedBackendCount}.`,
  );

  const solves: Array<{ configuredBackendId: string; outcome: string }> = [];
  for (const [index, request] of plan.requests.entries()) {
    const backend = backendFor(discovery, request);
    const validationErrors = validatePtFlashRequest(request, backend);
    requireCondition(
      validationErrors.length === 0,
      `Installed smoke request ${request.configuredBackendId} is invalid: ${validationErrors.join(' ')}`,
    );
    const solveReply = (await window.webContents.executeJavaScript(
      `globalThis.mpmcPtDesktop.solvePtFlash(
        'installed-desktop-smoke-solve-${index}',
        ${JSON.stringify(request)}
      )`,
      true,
    )) as DesktopBridgeReply<PtFlashResponse>;
    const response = requireBridgeValue(
      solveReply,
      `installed desktop solve ${request.configuredBackendId}`,
    );
    requireCondition(
      response.kind === 'result',
      `Installed desktop solve returned a service error for ${request.configuredBackendId}.`,
    );
    requireCondition(
      response.result.provenance.backend.configuredBackendId ===
        request.configuredBackendId,
      `Installed desktop solve provenance changed for ${request.configuredBackendId}.`,
    );
    solves.push({
      configuredBackendId: request.configuredBackendId,
      outcome: response.result.outcome,
    });
  }

  const result = {
    convention: PT_DESKTOP_INSTALL_SMOKE_RESULT_CONVENTION,
    discoveredBackendCount: discovery.backends.length,
    solves,
  };
  await writeFile(resultPath, `${JSON.stringify(result, null, 2)}\n`, 'utf8');
}
