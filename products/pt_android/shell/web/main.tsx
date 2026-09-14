import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';

import { App } from '../../../../frontend/src/App';
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
