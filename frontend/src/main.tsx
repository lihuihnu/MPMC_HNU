import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';

import { App } from './App';
import { unconfiguredFlashClient } from './api/flashClient';
import './styles.css';

const rootElement = document.getElementById('root');
if (rootElement === null) {
  throw new Error('MPMC_HNU frontend root element is missing.');
}

createRoot(rootElement).render(
  <StrictMode>
    <App client={unconfiguredFlashClient} />
  </StrictMode>,
);
