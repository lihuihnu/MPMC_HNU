import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';

import { App } from './App';
import {
  createGrpcWebFlashClient,
  unconfiguredFlashClient,
} from './api/flashClient';
import { desktopFlashClient } from './api/desktopBridge';
import type { HostedWebModelOwnershipCapability } from './api/hostedWebModelOwner';
import './styles.css';
import './prProduct.css';

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

async function renderHostedProduct(
  capability: HostedWebModelOwnershipCapability,
): Promise<boolean> {
  const [{ hostedWebExpertModelOwner }, { DesktopProductShell }] = await Promise.all([
    import('./api/hostedWebModelOwner'),
    import('./components/DesktopProductShell'),
  ]);
  const expertOwner = hostedWebExpertModelOwner(capability);
  if (expertOwner === null) return false;

  const dispose = () => {
    void expertOwner.dispose().catch(() => {});
  };
  window.addEventListener('pagehide', dispose, { once: true });
  try {
    root.render(
      <StrictMode>
        <DesktopProductShell
          flashClient={flashClient}
          expertOwner={expertOwner}
          environmentLabel="Hosted calculation · authenticated model session"
          workspaceAccessLabel="Hosted calculation · authenticated session"
        />
      </StrictMode>,
    );
    return true;
  } catch {
    window.removeEventListener('pagehide', dispose);
    await expertOwner.dispose().catch(() => {});
    return false;
  }
}

async function renderProduct() {
  // The editable-model bridge remains the authoritative local Electron path.
  const modelBridge = window.mpmcModelWorkbench;
  if (modelBridge) {
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
      return;
    } catch {
      // Never expose dynamic-loader or preload exception text. The PT path remains
      // usable, while desktop product smoke can detect this explicit fallback.
      renderPtOnly(true);
      return;
    }
  }

  // Hosted Classic PR is opt-in and fail-closed. The trusted application shell
  // supplies only a versioned same-origin endpoint + identity-provider function;
  // bearer values remain inside webModelSession.ts and are acquired on first apply.
  const hostedCapability = window.mpmcHostedWebModelOwnership;
  if (hostedCapability !== undefined) {
    try {
      if (await renderHostedProduct(hostedCapability)) return;
    } catch {
      // A malformed/unavailable hosted ownership integration must not turn into
      // anonymous model access or leak identity/transport details into the UI.
    }
  }

  renderPtOnly();
}

void renderProduct();
