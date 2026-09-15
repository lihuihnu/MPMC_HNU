import { useMemo, useState } from 'react';
import { Code } from '@connectrpc/connect';

import { App } from '../App';
import type { FlashClient } from '../api/flashClient';
import { bindExpertModelOwner, type ExpertModelOwner } from '../api/expertModelOwner';
import { RendererModelError, type RendererModelClient } from '../api/rendererModelClient';
import type { ModelValidationDetail } from '../api/modelValidationDetail';
import { ExpertPr76Workspace } from './ExpertPr76Workspace';

export type DesktopProductMode = 'pt' | 'expert';
export interface DesktopExpertFailure {
  readonly reason: string;
  readonly code?: Code;
  readonly validation?: ModelValidationDetail;
}
export type DesktopExpertSessionState =
  | { readonly status: 'idle' }
  | { readonly status: 'loading' }
  | { readonly status: 'ready' }
  | { readonly status: 'failed'; readonly failure: DesktopExpertFailure };

export function desktopExpertFailure(cause: unknown): DesktopExpertFailure {
  if (cause instanceof RendererModelError) {
    return Object.freeze({
      reason: cause.reason,
      code: cause.code,
      ...(cause.validation === undefined ? {} : { validation: cause.validation }),
    });
  }
  return Object.freeze({ reason: 'desktop.expert_session_failed' });
}

export interface DesktopProductShellViewProps {
  flashClient: FlashClient;
  expertOwner: ExpertModelOwner;
  mode: DesktopProductMode;
  expertSession: DesktopExpertSessionState;
  onPtMode(): void;
  onExpertMode(): void;
  onReconnect(): void;
}

export function DesktopProductShellView({
  flashClient,
  expertOwner,
  mode,
  expertSession,
  onPtMode,
  onExpertMode,
  onReconnect,
}: DesktopProductShellViewProps) {
  return (
    <>
      <nav className="desktop-product-mode" data-desktop-product-shell="true" aria-label="Desktop product mode">
        <div>
          <strong>Desktop workspace</strong>
          <span>PT Flash stays the default; Expert models use the authenticated window session.</span>
        </div>
        <div className="submit-actions">
          <button
            type="button"
            className={mode === 'pt' ? 'primary-button' : 'secondary-button cancel-button'}
            data-product-mode="pt"
            aria-pressed={mode === 'pt'}
            onClick={onPtMode}
          >
            PT Flash
          </button>
          <button
            type="button"
            className={mode === 'expert' ? 'primary-button' : 'secondary-button cancel-button'}
            data-product-mode="expert"
            aria-pressed={mode === 'expert'}
            disabled={expertSession.status === 'loading'}
            onClick={onExpertMode}
          >
            PR76 Expert
          </button>
        </div>
      </nav>

      {mode === 'pt' ? <App client={flashClient} /> : null}

      {mode === 'expert' && expertSession.status === 'idle' ? (
        <section className="app-shell">
          <div className="result-panel result-empty" data-expert-session="idle" aria-live="polite">
            <p className="eyebrow">Expert session</p>
            <h2>Expert session not opened</h2>
            <p>Select PR76 Expert again to open the authenticated window-scoped model session.</p>
          </div>
        </section>
      ) : null}

      {mode === 'expert' && expertSession.status === 'loading' ? (
        <section className="app-shell">
          <div className="result-panel result-empty" data-expert-session="loading" aria-live="polite">
            <p className="eyebrow">Expert session</p>
            <h2>Opening authenticated model session</h2>
            <p>No model is created until the explicit Expert editor is submitted.</p>
          </div>
        </section>
      ) : null}

      {mode === 'expert' && expertSession.status === 'failed' ? (
        <section className="app-shell">
          <div className="result-panel result-empty client-failure" data-expert-session="failed" aria-live="polite">
            <p className="eyebrow">Expert session</p>
            <h2>Expert model session unavailable</h2>
            <p>
              Reason <code>{expertSession.failure.reason}</code>
              {expertSession.failure.code === undefined ? null : (
                <> · gRPC code <code>{expertSession.failure.code}</code></>
              )}
            </p>
            {expertSession.failure.validation ? (
              <div className="validation-panel">
                <strong>Structured validation detail</strong>
                <p>
                  <code>{expertSession.failure.validation.code}</code> · field{' '}
                  <code>{expertSession.failure.validation.field ?? 'not supplied'}</code>
                </p>
              </div>
            ) : null}
            <p className="failure-semantics">
              Reconnection is never automatic. The button below is an explicit request to replace the window model session; existing model references become stale.
            </p>
            <button type="button" className="secondary-button" onClick={onReconnect}>
              Reconnect Expert session
            </button>
          </div>
        </section>
      ) : null}

      {mode === 'expert' && expertSession.status === 'ready' ? (
        <div data-expert-session="ready">
          <ExpertPr76Workspace owner={expertOwner} />
        </div>
      ) : null}
    </>
  );
}

export interface DesktopProductShellProps {
  flashClient: FlashClient;
  modelClient: RendererModelClient;
}

/** Electron-only product shell. No model session is opened until Expert mode is explicitly selected. */
export function DesktopProductShell({ flashClient, modelClient }: DesktopProductShellProps) {
  const expertOwner = useMemo(() => bindExpertModelOwner(modelClient), [modelClient]);
  const [mode, setMode] = useState<DesktopProductMode>('pt');
  const [expertSession, setExpertSession] = useState<DesktopExpertSessionState>({ status: 'idle' });

  async function enterExpert() {
    setMode('expert');
    if (expertSession.status === 'ready' && !modelClient.requiresReconnect) return;
    if (modelClient.requiresReconnect) {
      setExpertSession({
        status: 'failed',
        failure: { reason: 'renderer.reconnect_required', code: Code.FailedPrecondition },
      });
      return;
    }
    setExpertSession({ status: 'loading' });
    try {
      await modelClient.connect();
      setExpertSession({ status: 'ready' });
    } catch (cause) {
      setExpertSession({ status: 'failed', failure: desktopExpertFailure(cause) });
    }
  }

  async function reconnectExpert() {
    setExpertSession({ status: 'loading' });
    try {
      await modelClient.reconnect();
      setExpertSession({ status: 'ready' });
    } catch (cause) {
      setExpertSession({ status: 'failed', failure: desktopExpertFailure(cause) });
    }
  }

  return (
    <DesktopProductShellView
      flashClient={flashClient}
      expertOwner={expertOwner}
      mode={mode}
      expertSession={expertSession}
      onPtMode={() => setMode('pt')}
      onExpertMode={() => void enterExpert()}
      onReconnect={() => void reconnectExpert()}
    />
  );
}
