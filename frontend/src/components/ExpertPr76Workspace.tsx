import { useEffect, useRef, useState } from 'react';
import { Code } from '@connectrpc/connect';

import type { ExpertModelOwner, ExpertOwnedModel } from '../api/expertModelOwner';
import { ModelClientError } from '../api/modelSessionClient';
import { RendererModelError } from '../api/rendererModelClient';
import type { ModelValidationDetail } from '../api/modelValidationDetail';
import { ExpertModelSurface } from './ExpertModelSurface';
import { ExpertPr76Editor } from './ExpertPr76Editor';
import { ExpertPtSolvePanel } from './ExpertPtSolvePanel';

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
 * Replace the applied snapshot only after create succeeds. A failed mutation
 * may invalidate the typed client's lease; keeping an old snapshot does not
 * promise that its reference remains usable. No release is retried.
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

/**
 * Defer irreversible teardown release until the current JS turn completes.
 * React StrictMode's development setup→cleanup→setup replay can invalidate the
 * cleanup token before release admission; a real unmount still releases only
 * the model that is current when this deferred callback executes.
 */
export function deferExpertWorkspaceRelease(
  readCurrent: () => ExpertOwnedModel | null,
  stillUnmounted: () => boolean,
  enqueue: (task: () => void) => void = queueMicrotask,
): void {
  enqueue(() => {
    if (!stillUnmounted()) return;
    const current = readCurrent();
    if (current === null || current.released) return;
    void current.release().catch(() => {
      // Teardown cannot publish UI state. The typed client/session remains the
      // authority for ambiguous release and eventual session cleanup.
    });
  });
}

export interface ExpertPr76WorkspaceViewProps {
  owner: ExpertModelOwner;
  current: ExpertOwnedModel | null;
  generation: number;
  retirementFailure: ExpertWorkspaceFailure | null;
  onCreated(model: ExpertOwnedModel): void;
  unapplied?: boolean;
  onDraftChanged?(): void;
}

export function ExpertPr76WorkspaceView({
  owner,
  current,
  generation,
  retirementFailure,
  onCreated,
  unapplied = false,
  onDraftChanged,
}: ExpertPr76WorkspaceViewProps) {
  return (
    <div className="app-shell expert-pr76-workspace" data-expert-workspace="pr76">
      <header className="app-header">
        <div>
          <p className="brand-mark">MPMC_HNU · Expert</p>
          <h1>PR76 model workspace</h1>
          <p>
            Edit component parameters and kij, apply the model, then enter P, T and z
            for C++ phase-stability analysis and flash calculations up to three phases.
          </p>
        </div>
        <div className="backend-state" data-configured="true">
          <span className="backend-dot" />
          Local computation · no account or login required
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
            {...(onDraftChanged === undefined ? {} : { onDraftChanged })}
            onCreated={onCreated}
          />
        </section>
        <section className="workspace-result">
          {current === null ? (
            <section className="result-panel result-empty" aria-live="polite">
              <p className="eyebrow">Live model</p>
              <h2>No Expert model created yet</h2>
              <p>
                Complete and apply the explicit PR76 draft on the left. The calculation
                panel will use that model, not a preconfigured example fluid.
              </p>
            </section>
          ) : (
            <>
              <ExpertPtSolvePanel key={generation} model={current} blocked={unapplied || retirementFailure !== null} />
              <details className="contract-details"><summary>Applied model parameters and provenance</summary>
                <ExpertModelSurface source={current.source} />
              </details>
            </>
          )}
        </section>
      </main>

      <footer className="app-footer">
        <span>Parameters are editable data; EOS and flash algorithms remain in C++</span>
        <span>No implicit kij, feed normalization, model fallback or account setup</span>
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
  const [unapplied, setUnapplied] = useState(false);
  const [retirementFailure, setRetirementFailure] = useState<ExpertWorkspaceFailure | null>(null);
  const currentRef = useRef<ExpertOwnedModel | null>(null);
  const alive = useRef(true);
  const mountToken = useRef<symbol | null>(null);

  useEffect(() => {
    const token = Symbol('expert-pr76-workspace-mount');
    mountToken.current = token;
    alive.current = true;
    return () => {
      alive.current = false;
      deferExpertWorkspaceRelease(
        () => currentRef.current,
        () => mountToken.current === token && !alive.current,
      );
    };
  }, []);

  function handleCreated(next: ExpertOwnedModel) {
    const previous = currentRef.current;
    currentRef.current = next;
    setCurrent(next);
    setGeneration((value) => value + 1);
    setUnapplied(false);
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
      unapplied={unapplied}
      onDraftChanged={() => setUnapplied(true)}
    />
  );
}
