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

createRoot(rootElement).render(
  <StrictMode>
    <App client={flashClient} />
  </StrictMode>,
);
