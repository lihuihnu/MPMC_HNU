import { useEffect, useState } from 'react';
import { Code } from '@connectrpc/connect';

import type { ExpertModelInspection } from '../api/expertModelInspector';
import type { ExpertModelSource } from '../api/expertModelSource';
import { ModelClientError } from '../api/modelSessionClient';
import { RendererModelError } from '../api/rendererModelClient';
import {
  MODEL_VALIDATION_DETAIL_VERSION,
  type ModelValidationDetail,
} from '../api/modelValidationDetail';
import { ExpertModelInspector } from './ExpertModelInspector';

export interface ExpertModelSurfaceFailure {
  readonly reason: string;
  readonly code?: Code;
  readonly validation?: ModelValidationDetail;
}

export type ExpertModelSurfaceState =
  | { readonly status: 'loading' }
  | { readonly status: 'ready'; readonly inspection: ExpertModelInspection }
  | { readonly status: 'failed'; readonly failure: ExpertModelSurfaceFailure };

/** Keep user-facing failures on the existing typed error contracts; never show what()/message text. */
export function expertModelSurfaceFailure(cause: unknown): ExpertModelSurfaceFailure {
  if (cause instanceof ModelClientError || cause instanceof RendererModelError) {
    return Object.freeze({
      reason: cause.reason,
      code: cause.code,
      ...(cause.validation === undefined ? {} : { validation: cause.validation }),
    });
  }
  return Object.freeze({ reason: 'expert.describe_failed' });
}

interface ExpertModelSurfaceViewProps {
  state: ExpertModelSurfaceState;
  onRefresh(): void;
}

export function ExpertModelSurfaceView({ state, onRefresh }: ExpertModelSurfaceViewProps) {
  if (state.status === 'loading') {
    return (
      <section
        className="result-panel result-empty"
        data-expert-state="loading"
        aria-live="polite"
      >
        <p className="eyebrow">Expert model</p>
        <h2>Reading live model snapshot</h2>
        <p>
          The surface is calling the owning typed client&apos;s read-only describe path.
          No configuration is recreated from browser state.
        </p>
      </section>
    );
  }

  if (state.status === 'failed') {
    const { failure } = state;
    return (
      <section
        className="result-panel result-empty client-failure"
        data-expert-state="failed"
        aria-live="polite"
      >
        <p className="eyebrow">Expert model</p>
        <h2>Live model snapshot unavailable</h2>
        <p>
          Reason <code>{failure.reason}</code>
          {failure.code === undefined ? null : (
            <> · gRPC code <code>{failure.code}</code></>
          )}
        </p>
        {failure.validation ? (
          <div
            className="validation-panel"
            data-validation-version={MODEL_VALIDATION_DETAIL_VERSION}
          >
            <strong>Structured validation detail</strong>
            <ul>
              <li>
                Code <code>{failure.validation.code}</code>
              </li>
              <li>
                Field <code>{failure.validation.field ?? 'not supplied'}</code>
              </li>
            </ul>
          </div>
        ) : null}
        <p className="failure-semantics">
          No stale snapshot is shown and the surface does not reconnect, recreate or
          release the owning model reference automatically.
        </p>
        <button className="secondary-button" type="button" onClick={onRefresh}>
          Retry describe
        </button>
      </section>
    );
  }

  return (
    <section data-expert-state="ready" aria-live="polite">
      <div className="submit-row">
        <div className="submit-note">
          Live <code>describe()</code> snapshot · session/reference remain owned by the
          typed client
        </div>
        <button className="secondary-button cancel-button" type="button" onClick={onRefresh}>
          Refresh snapshot
        </button>
      </div>
      <ExpertModelInspector inspection={state.inspection} />
    </section>
  );
}

export interface ExpertModelSurfaceProps {
  source: ExpertModelSource;
}

/**
 * User-facing read-only shell around one already-owned model reference.
 * A source change or manual refresh aborts the prior describe request; no retry or
 * session mutation happens behind the caller's back.
 */
export function ExpertModelSurface({ source }: ExpertModelSurfaceProps) {
  const [attempt, setAttempt] = useState(0);
  const [state, setState] = useState<ExpertModelSurfaceState>({ status: 'loading' });

  useEffect(() => {
    const controller = new AbortController();
    setState({ status: 'loading' });
    void source.describe({ signal: controller.signal }).then(
      (inspection) => {
        if (!controller.signal.aborted) setState({ status: 'ready', inspection });
      },
      (cause: unknown) => {
        if (!controller.signal.aborted) {
          setState({ status: 'failed', failure: expertModelSurfaceFailure(cause) });
        }
      },
    );
    return () => controller.abort();
  }, [source, attempt]);

  return <ExpertModelSurfaceView state={state} onRefresh={() => setAttempt((value) => value + 1)} />;
}
