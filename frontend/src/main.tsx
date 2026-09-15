import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';

import { App } from './App';
import {
  createGrpcWebFlashClient,
  unconfiguredFlashClient,
} from './api/flashClient';
import { desktopFlashClient } from './api/desktopBridge';
import './styles.css';

const rootElement = document.getElementById('root');
if (rootElement === null) {
  throw new Error('MPMC_HNU frontend root element is missing.');
}

const configuredBaseUrl = import.meta.env.VITE_MPMC_GRPC_WEB_BASE_URL;
const localDesktopClient = desktopFlashClient();
const flashClient = localDesktopClient ?? (
  typeof configuredBaseUrl === 'string'
    ? createGrpcWebFlashClient(configuredBaseUrl)
    : unconfiguredFlashClient
);
const root = createRoot(rootElement);

function renderPtOnly(expertShellUnavailable = false) {
  root.render(
    <StrictMode>
      {expertShellUnavailable ? (
        <div className="app-shell" data-desktop-product-shell="failed" role="alert">
          <div className="connection-banner connection-banner-error">
            <strong>Local model workbench unavailable.</strong>
            <span>PT Flash remains available; no editable model was applied.</span>
          </div>
        </div>
      ) : null}
      <App client={flashClient} />
    </StrictMode>,
  );
}

async function renderProduct() {
  // The editable-model bridge is a narrow Electron preload capability. Ordinary
  // Web/Android builds do not receive native model ownership or process transport.
  const modelBridge = window.mpmcModelWorkbench;
  if (!modelBridge) {
    renderPtOnly();
    return;
  }

  try {
    const [{ modelWorkbenchOwner }, { DesktopProductShell }] = await Promise.all([
      import('./api/modelWorkbenchOwner'),
      import('./components/DesktopProductShell'),
    ]);
    const expertOwner = modelWorkbenchOwner(modelBridge);
    if (!expertOwner) throw new Error('model workbench unavailable');
    root.render(
      <StrictMode>
        <DesktopProductShell flashClient={flashClient} expertOwner={expertOwner} />
      </StrictMode>,
    );
  } catch {
    // Never expose dynamic-loader or preload exception text. The PT path remains
    // usable, while desktop product smoke can detect this explicit fallback.
    renderPtOnly(true);
  }
}

void renderProduct();
