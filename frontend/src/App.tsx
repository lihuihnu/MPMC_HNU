import { useEffect, useRef, useState } from 'react';

import {
  BackendNotConfiguredError,
  PtTransportError,
  transportCodeLabel,
  type FlashClient,
} from './api/flashClient';
import { PtWireContractError } from './api/ptWire';
import { FlashForm } from './components/FlashForm';
import { ResultPanel, type PtClientFailure } from './components/ResultPanel';
import type {
  PtBackendDescriptor,
  PtCapabilityDiscovery,
  PtFlashRequest,
  PtFlashResponse,
} from './domain/flash';

interface AppProps {
  client: FlashClient;
}

type DiscoveryState =
  | { status: 'unconfigured' }
  | { status: 'loading' }
  | { status: 'ready'; discovery: PtCapabilityDiscovery }
  | { status: 'failed'; failure: PtClientFailure };

function clientFailure(cause: unknown, operation: string): PtClientFailure {
  if (cause instanceof BackendNotConfiguredError) {
    return { kind: 'configuration', title: 'Endpoint not configured', message: cause.message };
  }
  if (cause instanceof PtWireContractError) {
    return {
      kind: 'contract',
      title: 'Wire contract rejected',
      message: cause.message,
    };
  }
  if (cause instanceof PtTransportError) {
    return {
      kind: 'transport',
      title: `${operation} RPC ${transportCodeLabel(cause.code)}`,
      message: cause.message,
    };
  }
  return {
    kind: 'transport',
    title: `${operation} RPC failed`,
    message: cause instanceof Error ? cause.message : 'The client received an unknown failure.',
  };
}

export function App({ client }: AppProps) {
  const [discoveryState, setDiscoveryState] = useState<DiscoveryState>(
    client.configured ? { status: 'loading' } : { status: 'unconfigured' },
  );
  const [discoveryAttempt, setDiscoveryAttempt] = useState(0);
  const [busy, setBusy] = useState(false);
  const [response, setResponse] = useState<PtFlashResponse | null>(null);
  const [solveFailure, setSolveFailure] = useState<PtClientFailure | null>(null);
  const activeSolve = useRef<AbortController | null>(null);
  const mounted = useRef(true);

  useEffect(() => {
    if (!client.configured) {
      setDiscoveryState({ status: 'unconfigured' });
      return;
    }

    const controller = new AbortController();
    setDiscoveryState({ status: 'loading' });
    void client
      .discoverPtCapabilities({ signal: controller.signal })
      .then((discovery) => {
        if (!controller.signal.aborted) {
          setDiscoveryState({ status: 'ready', discovery });
        }
      })
      .catch((cause: unknown) => {
        if (!controller.signal.aborted) {
          setDiscoveryState({ status: 'failed', failure: clientFailure(cause, 'Discovery') });
        }
      });
    return () => controller.abort();
  }, [client, discoveryAttempt]);

  useEffect(() => {
    mounted.current = true;
    return () => {
      mounted.current = false;
      activeSolve.current?.abort();
    };
  }, []);

  async function handleSubmit(request: PtFlashRequest, backend: PtBackendDescriptor) {
    const controller = new AbortController();
    activeSolve.current?.abort();
    activeSolve.current = controller;
    setBusy(true);
    setResponse(null);
    setSolveFailure(null);
    try {
      const nextResponse = await client.solvePtFlash(request, backend, {
        signal: controller.signal,
      });
      if (!controller.signal.aborted) {
        setResponse(nextResponse);
      }
    } catch (cause) {
      if (mounted.current) {
        setResponse(null);
        setSolveFailure(clientFailure(cause, 'Solve'));
      }
    } finally {
      if (activeSolve.current === controller) {
        activeSolve.current = null;
        if (mounted.current) {
          setBusy(false);
        }
      }
    }
  }

  function cancelSolve() {
    activeSolve.current?.abort();
  }

  function resetResult() {
    setResponse(null);
    setSolveFailure(null);
  }

  const backends =
    discoveryState.status === 'ready' ? discoveryState.discovery.backends : [];
  const connectionLabel = (() => {
    switch (discoveryState.status) {
      case 'unconfigured':
        return 'gRPC-Web endpoint not configured';
      case 'loading':
        return 'Discovering service capabilities';
      case 'failed':
        return 'Capability discovery failed';
      case 'ready':
        return `${discoveryState.discovery.backends.length} backend${
          discoveryState.discovery.backends.length === 1 ? '' : 's'
        } discovered`;
    }
  })();

  return (
    <div className="app-shell">
      <header className="app-header">
        <div>
          <p className="brand-mark">MPMC_HNU</p>
          <h1>Model-neutral PT Flash</h1>
          <p>
            Runtime backend discovery, exact component inventories and variable phase
            publication through the versioned PT service boundary.
          </p>
        </div>
        <div
          className="backend-state"
          data-configured={discoveryState.status === 'ready'}
        >
          <span className="backend-dot" />
          {connectionLabel}
        </div>
      </header>

      {discoveryState.status === 'unconfigured' ? (
        <div className="connection-banner" role="status">
          <strong>No compute endpoint configured.</strong>
          <span>
            Set <code>VITE_MPMC_GRPC_WEB_BASE_URL</code> to a service implementing{' '}
            <code>mpmc.runtime.v1.PtFlashService</code>. No calculation data is mocked.
          </span>
        </div>
      ) : null}

      {discoveryState.status === 'failed' ? (
        <div className="connection-banner connection-banner-error" role="alert">
          <strong>{discoveryState.failure.title}</strong>
          <span>{discoveryState.failure.message}</span>
          <button
            className="secondary-button banner-button"
            type="button"
            onClick={() => setDiscoveryAttempt((value) => value + 1)}
          >
            Retry discovery
          </button>
        </div>
      ) : null}

      {discoveryState.status === 'ready' && backends.length === 0 ? (
        <div className="connection-banner" role="status">
          <strong>Service returned no configured PT backends.</strong>
          <span>Calculation stays disabled until the worker publishes an inventory.</span>
        </div>
      ) : null}

      <main className="workspace">
        <section className="workspace-input">
          <FlashForm
            backends={backends}
            discoveryLoading={discoveryState.status === 'loading'}
            busy={busy}
            onBackendChange={resetResult}
            onCancel={cancelSolve}
            onSubmit={(request, backend) => void handleSubmit(request, backend)}
          />
        </section>
        <section className="workspace-result">
          <ResultPanel response={response} clientFailure={solveFailure} />
        </section>
      </main>

      <footer className="app-footer">
        <span>Wire: mpmc.runtime.v1 · internal PT units are Pa and K</span>
        <span>No browser-side EOS, flash, normalization, retry or fallback</span>
      </footer>
    </div>
  );
}
