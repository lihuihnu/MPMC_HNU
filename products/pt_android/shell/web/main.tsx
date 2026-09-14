import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';

import { App } from '../../../../frontend/src/App';
import type {
  PtBackendDescriptor,
  PtFlashRequest,
} from '../../../../frontend/src/domain/flash';
import '../../../../frontend/src/styles.css';
import { createAndroidFlashClient } from './androidFlashClient';

const rootElement = document.getElementById('root');
if (rootElement === null) {
  throw new Error('MPMC_HNU Android product shell root element is missing.');
}

const client = createAndroidFlashClient();

createRoot(rootElement).render(
  <StrictMode>
    <App client={client} />
  </StrictMode>,
);

function requiredBackend(
  backends: readonly PtBackendDescriptor[],
  configuredBackendId: string,
): PtBackendDescriptor {
  const backend = backends.find(
    (candidate) => candidate.configuredBackendId === configuredBackendId,
  );
  if (backend === undefined) {
    throw new Error(`Android product shell smoke is missing ${configuredBackendId}.`);
  }
  return backend;
}

async function acceptedSolve(
  backend: PtBackendDescriptor,
  request: PtFlashRequest,
): Promise<number> {
  const response = await client.solvePtFlash(request, backend);
  if (response.kind !== 'result' || response.result.outcome !== 'accepted') {
    throw new Error(
      `Android product shell smoke did not accept ${backend.configuredBackendId}.`,
    );
  }
  return response.result.phases.length;
}

async function runProductShellSmoke(): Promise<void> {
  const discovery = await client.discoverPtCapabilities();
  if (discovery.backends.length !== 3) {
    throw new Error(
      `Android product shell smoke expected 3 backends, found ${discovery.backends.length}.`,
    );
  }

  const pr76 = requiredBackend(
    discovery.backends,
    'pr76.methane-ethane-propane.literature-r1',
  );
  const sw92 = requiredBackend(
    discovery.backends,
    'sw92.carbon-dioxide-water.freshwater.literature-r1',
  );
  const cpa = requiredBackend(
    discovery.backends,
    'cpa.methanol-water-333.15k.cr1.literature-r1',
  );

  const phaseCounts = await Promise.all([
    acceptedSolve(pr76, {
      configuredBackendId: pr76.configuredBackendId,
      pressurePa: 1.0e6,
      temperatureK: 350.0,
      feed: [
        { componentId: 'methane', moleFraction: 0.8 },
        { componentId: 'ethane', moleFraction: 0.1 },
        { componentId: 'propane', moleFraction: 0.1 },
      ],
    }),
    acceptedSolve(sw92, {
      configuredBackendId: sw92.configuredBackendId,
      pressurePa: 3.0e6,
      temperatureK: 340.0,
      feed: [
        { componentId: 'carbon-dioxide', moleFraction: 0.7 },
        { componentId: 'water', moleFraction: 0.3 },
      ],
    }),
    acceptedSolve(cpa, {
      configuredBackendId: cpa.configuredBackendId,
      pressurePa: 48_852.0,
      temperatureK: 333.15,
      feed: [
        { componentId: 'METHANOL', moleFraction: 0.5 },
        { componentId: 'WATER', moleFraction: 0.5 },
      ],
    }),
  ]);

  document.documentElement.dataset.mpmcAndroidProductShellSmoke = 'ok';
  console.info(
    `ANDROID_PRODUCT_SHELL_WEB_OK backends=3 pr76_phases=${phaseCounts[0]} sw92_phases=${phaseCounts[1]} cpa_phases=${phaseCounts[2]}`,
  );
}

if (import.meta.env.VITE_MPMC_ANDROID_PRODUCT_SHELL_SMOKE === '1') {
  void runProductShellSmoke().catch((cause: unknown) => {
    console.error(
      `ANDROID_PRODUCT_SHELL_WEB_FAIL ${cause instanceof Error ? cause.message : String(cause)}`,
    );
  });
}
