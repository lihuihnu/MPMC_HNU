import { useEffect, useRef, useState } from 'react';
import { Code } from '@connectrpc/connect';

import type { ExpertModelOwner, ExpertOwnedModel } from '../api/expertModelOwner';
import { ModelClientError } from '../api/modelSessionClient';
import { RendererModelError } from '../api/rendererModelClient';
import type { ModelValidationDetail } from '../api/modelValidationDetail';
import { ExpertModelSurface } from './ExpertModelSurface';
import { ExpertPr76Editor } from './ExpertPr76Editor';

export interface ExpertWorkspaceFailure {
  readonly reason: string;
  readonly code?: Code;
  readonly validation?: ModelValidationDetail;
}

export function expertWorkspaceFailure(cause: unknown): ExpertWorkspaceFailure {
  if (cause instanceof ModelClientError || cause instanceof RendererModelError) {
    return Object.freeze({
      reason: cause.reason,
      code: cause.code,
      ...(cause.validation === undefined ? {} : { validation: cause.validation }),
    });
  }
  return Object.freeze({ reason: 'expert.release_failed' });
}

/**
 * A successful create becomes current before the previous model is retired.
 * This preserves the old model when create itself fails, while still bounding
 * steady-state ownership to one model per workspace. No release is retried.
 */
export async function retirePreviousExpertModel(
  previous: ExpertOwnedModel | null,
  next: ExpertOwnedModel,
): Promise<ExpertWorkspaceFailure | null> {
  if (previous === null || previous === next || previous.released) return null;
  try {
    await previous.release();
    return null;
  } catch (cause) {
    return expertWorkspaceFailure(cause);
  }
}

export interface ExpertPr76WorkspaceViewProps {
  owner: ExpertModelOwner;
  current: ExpertOwnedModel | null;
  generation: number;
  retirementFailure: ExpertWorkspaceFailure | null;
  onCreated(model: ExpertOwnedModel): void;
}

export function ExpertPr76WorkspaceView({
  owner,
  current,
  generation,
  retirementFailure,
  onCreated,
}: ExpertPr76WorkspaceViewProps) {
  return (
    <div className="app-shell expert-pr76-workspace" data-expert-workspace="pr76">
      <header className="app-header">
        <div>
          <p className="brand-mark">MPMC_HNU · Expert</p>
          <h1>PR76 model workspace</h1>
          <p>
            Create immutable parameter/settings revisions on the authenticated model
            service, then inspect the live server-owned snapshot through describe().
          </p>
        </div>
        <div className="backend-state" data-configured="true">
          <span className="backend-dot" />
          Opaque session-owned model references
        </div>
      </header>

      {retirementFailure ? (
        <div className="connection-banner connection-banner-error" role="alert">
          <strong>Previous model cleanup was not confirmed.</strong>
          <span>
            Reason <code>{retirementFailure.reason}</code>
            {retirementFailure.code === undefined ? null : (
              <> · gRPC code <code>{retirementFailure.code}</code></>
            )}
            {retirementFailure.validation?.field === undefined ? null : (
              <> · field <code>{retirementFailure.validation.field}</code></>
            )}
          </span>
        </div>
      ) : null}

      <main className="workspace">
        <section className="workspace-input">
          <ExpertPr76Editor
            key={generation}
            owner={owner}
            {...(current === null ? {} : { seedSnapshot: current.snapshot })}
            onCreated={onCreated}
          />
        </section>
        <section className="workspace-result">
          {current === null ? (
            <section className="result-panel result-empty" aria-live="polite">
              <p className="eyebrow">Live model</p>
              <h2>No Expert model created yet</h2>
              <p>
                Complete the explicit PR76 draft on the left. A successful create will
                return a session-local opaque reference and this panel will read the
                authoritative snapshot through the existing typed describe path.
              </p>
            </section>
          ) : (
            <ExpertModelSurface source={current.source} />
          )}
        </section>
      </main>

      <footer className="app-footer">
        <span>Expert mutations create immutable models; existing handles are never edited in place</span>
        <span>No browser-side EOS, parameter fallback, implicit kij, retry or reconnect</span>
      </footer>
    </div>
  );
}

export interface ExpertPr76WorkspaceProps {
  owner: ExpertModelOwner;
}

/** Owns at most one current Expert model under normal successful cleanup. */
export function ExpertPr76Workspace({ owner }: ExpertPr76WorkspaceProps) {
  const [current, setCurrent] = useState<ExpertOwnedModel | null>(null);
  const [generation, setGeneration] = useState(0);
  const [retirementFailure, setRetirementFailure] = useState<ExpertWorkspaceFailure | null>(null);
  const currentRef = useRef<ExpertOwnedModel | null>(null);
  const alive = useRef(true);

  useEffect(() => {
    alive.current = true;
    return () => {
      alive.current = false;
      const owned = currentRef.current;
      currentRef.current = null;
      if (owned && !owned.released) {
        void owned.release().catch(() => {
          // Teardown cannot publish UI state. The typed client/session remains the
          // authority for ambiguous release and eventual session cleanup.
        });
      }
    };
  }, []);

  function handleCreated(next: ExpertOwnedModel) {
    const previous = currentRef.current;
    currentRef.current = next;
    setCurrent(next);
    setGeneration((value) => value + 1);
    setRetirementFailure(null);

    void retirePreviousExpertModel(previous, next).then((failure) => {
      if (failure && alive.current) setRetirementFailure(failure);
    });
  }

  return (
    <ExpertPr76WorkspaceView
      owner={owner}
      current={current}
      generation={generation}
      retirementFailure={retirementFailure}
      onCreated={handleCreated}
    />
  );
}
