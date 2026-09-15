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
            <strong>Expert desktop shell unavailable.</strong>
            <span>PT Flash remains available; no model session was opened.</span>
          </div>
        </div>
      ) : null}
      <App client={flashClient} />
    </StrictMode>,
  );
}

async function renderProduct() {
  // The model bridge is an Electron preload capability. Ordinary Web/Android do
  // not import or initialize the Expert model runtime until Web identity/session
  // integration is separately reviewed.
  const modelBridge = window.mpmcModelDesktopV2 ?? window.mpmcModelDesktop;
  if (!modelBridge) {
    renderPtOnly();
    return;
  }

  try {
    const [{ rendererModelClient }, { DesktopProductShell }] = await Promise.all([
      import('./api/rendererModelClient'),
      import('./components/DesktopProductShell'),
    ]);
    const modelClient = rendererModelClient(modelBridge);
    if (!modelClient) throw new Error('model bridge unavailable');
    root.render(
      <StrictMode>
        <DesktopProductShell flashClient={flashClient} modelClient={modelClient} />
      </StrictMode>,
    );
  } catch {
    // Never expose dynamic-loader or preload exception text. The PT path remains
    // usable, while desktop product smoke requires the shell marker and will fail
    // packaging if this fallback is reached.
    renderPtOnly(true);
  }
}

void renderProduct();
