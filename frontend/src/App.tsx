import { useState } from 'react';

import { BackendNotConfiguredError, type FlashClient } from './api/flashClient';
import { FlashForm } from './components/FlashForm';
import { ResultPanel } from './components/ResultPanel';
import type { ProfileCPtRequest, ProfileCPtResponse } from './domain/flash';

interface AppProps {
  client: FlashClient;
}

export function App({ client }: AppProps) {
  const [busy, setBusy] = useState(false);
  const [result, setResult] = useState<ProfileCPtResponse | null>(null);
  const [error, setError] = useState<string | null>(null);

  async function handleSubmit(request: ProfileCPtRequest) {
    setBusy(true);
    setError(null);
    try {
      const response = await client.solveProfileCPt(request);
      setResult(response);
    } catch (cause) {
      setResult(null);
      if (cause instanceof BackendNotConfiguredError) {
        setError(cause.message);
      } else if (cause instanceof Error) {
        setError(cause.message);
      } else {
        setError('The compute service returned an unknown error.');
      }
    } finally {
      setBusy(false);
    }
  }

  return (
    <div className="app-shell">
      <header className="app-header">
        <div>
          <p className="brand-mark">MPMC_HNU</p>
          <h1>Multiphase PT Flash</h1>
          <p>
            Profile-C topology orchestration for aqueous and nonaqueous phase
            candidates.
          </p>
        </div>
        <div className="backend-state" data-configured={client.configured}>
          <span className="backend-dot" />
          {client.configured ? 'Compute service connected' : 'Compute service not configured'}
        </div>
      </header>

      {!client.configured ? (
        <div className="connection-banner" role="status">
          <strong>Backend transport is intentionally not mocked.</strong>
          <span>
            This frontend build contains the validated Profile-C PT contract and UI
            only. Connect a real MPMC_HNU service implementation before calculations
            can run.
          </span>
        </div>
      ) : null}

      <main className="workspace">
        <section className="workspace-input">
          <FlashForm
            backendConfigured={client.configured}
            busy={busy}
            onSubmit={(request) => void handleSubmit(request)}
          />
        </section>
        <section className="workspace-result">
          <ResultPanel result={result} error={error} />
        </section>
      </main>

      <footer className="app-footer">
        <span>Internal SI contract · browser validation is advisory</span>
        <span>EOS and flash authority remain in the C++ backend</span>
      </footer>
    </div>
  );
}
