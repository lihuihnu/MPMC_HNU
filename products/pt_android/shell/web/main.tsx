import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';

import { DesktopProductShell } from '../../../../frontend/src/components/DesktopProductShell';
import {
  addPr76Component,
  buildPr76CreateInput,
  editPr76ComponentScalar,
  editPr76ComponentText,
  editPr76Kij,
  newPr76ExpertDraft,
  type Pr76ExpertDraft,
} from '../../../../frontend/src/api/pr76ExpertDraft';
import { modelWorkbenchOwner } from '../../../../frontend/src/api/modelWorkbenchOwner';
import { PtComputationOutcome } from '../../../../frontend/src/gen/mpmc/runtime/v1/pt_service_pb';
import type {
  PtBackendDescriptor,
  PtFlashRequest,
} from '../../../../frontend/src/domain/flash';
import '../../../../frontend/src/styles.css';
import '../../../../frontend/src/prProduct.css';
import { createAndroidFlashClient } from './androidFlashClient';
import { createAndroidModelWorkbenchBridge } from './androidModelWorkbenchBridge';

const rootElement = document.getElementById('root');
if (rootElement === null) {
  throw new Error('MPMC_HNU Android product shell root element is missing.');
}

const client = createAndroidFlashClient();
const owner = modelWorkbenchOwner(createAndroidModelWorkbenchBridge());
if (owner === null) {
  throw new Error('MPMC_HNU Android Classic PR ownership adapter is unavailable.');
}

createRoot(rootElement).render(
  <StrictMode>
    <DesktopProductShell flashClient={client} expertOwner={owner} />
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

function addSmokeComponent(
  draft: Pr76ExpertDraft,
  componentId: string,
  displayName: string,
  criticalTemperatureK: number,
  criticalPressurePa: number,
  acentricFactor: number,
): Pr76ExpertDraft {
  let next = addPr76Component(draft);
  const component = next.components[next.components.length - 1]!;
  next = editPr76ComponentText(next, component.key, 'componentId', componentId);
  next = editPr76ComponentText(next, component.key, 'displayName', displayName);
  next = editPr76ComponentScalar(next, component.key, 'criticalTemperatureK', String(criticalTemperatureK));
  next = editPr76ComponentScalar(next, component.key, 'criticalPressurePa', String(criticalPressurePa));
  next = editPr76ComponentScalar(next, component.key, 'acentricFactor', String(acentricFactor));
  return next;
}

function classicPrSmokeDraft(): Pr76ExpertDraft {
  let draft: Pr76ExpertDraft = {
    ...newPr76ExpertDraft(),
    displayName: 'Android Classic PR integration regression',
    datasetId: 'android-classic-pr-repository-fixture',
    revision: 'android-classic-pr-smoke-r1',
    userRecord: {
      reference: 'modules/pt_process/src/parameter_snapshot_supplier.cpp',
      revision: 'repository-curated fixture in current source revision',
      locator: 'pr76_parameter_source_v1',
      acquisition: 'copied for Android Classic PR integration regression',
      usageTerms: 'test-only repository fixture; no external usage claim',
      note: 'Pure parameters cite https://doi.org/10.1002/aic.16730; explicit zero kij values cite https://doi.org/10.1021/i160057a011.',
    },
  };
  draft = addSmokeComponent(draft, 'methane', 'Methane', 190.555, 4.595e6, 0.0);
  draft = addSmokeComponent(draft, 'ethane', 'Ethane', 305.4, 4.88e6, 0.099);
  draft = addSmokeComponent(draft, 'propane', 'Propane', 369.825, 4.248e6, 0.15308);
  for (const pair of draft.pairs) draft = editPr76Kij(draft, pair.key, '0');
  return draft;
}

async function runClassicPrSmoke(): Promise<number> {
  const model = await owner.create(buildPr76CreateInput(classicPrSmokeDraft()));
  try {
    const result = await model.solve({
      pressurePa: 1.0e6,
      temperatureK: 350.0,
      feed: [0.8, 0.1, 0.1],
    });
    if (result.outcome !== PtComputationOutcome.ACCEPTED) {
      throw new Error('Android Classic PR model smoke did not return an accepted result.');
    }
    const phases = result.candidatePhaseSet?.phases.length ?? 0;
    if (phases < 1) throw new Error('Android Classic PR model smoke returned no accepted phase.');
    return phases;
  } finally {
    await model.release();
  }
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
  const classicPhases = await runClassicPrSmoke();

  document.documentElement.dataset.mpmcAndroidClassicPrSmoke = 'ok';
  document.documentElement.dataset.mpmcAndroidProductShellSmoke = 'ok';
  console.info(
    `ANDROID_PRODUCT_SHELL_WEB_OK backends=3 pr76_phases=${phaseCounts[0]} sw92_phases=${phaseCounts[1]} cpa_phases=${phaseCounts[2]} classic_pr_phases=${classicPhases}`,
  );
}

if (import.meta.env.VITE_MPMC_ANDROID_PRODUCT_SHELL_SMOKE === '1') {
  void runProductShellSmoke().catch((cause: unknown) => {
    console.error(
      `ANDROID_PRODUCT_SHELL_WEB_FAIL ${cause instanceof Error ? cause.message : String(cause)}`,
    );
  });
}
