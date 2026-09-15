import { clone } from '@bufbuild/protobuf';
import { Code } from '@connectrpc/connect';

import {
  ModelClientError,
  type ModelCallOptions,
  type ModelCreateInput,
  type ModelSolveInput,
} from './modelSessionClient';
import type { ExpertModelSource } from './expertModelSource';
import { bindExpertModelSource } from './expertModelSource';
import type { ModelSnapshotDescriber } from './expertModelInspector';
import {
  FullPtResultSchema,
  ModelSnapshotSchema,
  type FullPtResult,
  type ModelSnapshot,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';

export interface ExpertModelCreationClient<ModelReference>
  extends ModelSnapshotDescriber<ModelReference> {
  create(
    input: ModelCreateInput,
    options?: ModelCallOptions,
  ): Promise<{ model: ModelReference; snapshot: ModelSnapshot }>;
  // The shared client returns the full message; the renderer also supplies a
  // presentation outcome. Both must retain exactly the same native message.
  solve(
    model: ModelReference,
    input: ModelSolveInput,
    options?: ModelCallOptions,
  ): Promise<FullPtResult | { readonly result: FullPtResult }>;
  release(model: ModelReference, options?: ModelCallOptions): Promise<void>;
}

export interface ExpertOwnedModel {
  readonly snapshot: ModelSnapshot;
  readonly source: ExpertModelSource;
  readonly released: boolean;
  /** Solve this immutable model, not a repository-curated backend selected by ID. */
  solve(input: ModelSolveInput, options?: ModelCallOptions): Promise<FullPtResult>;
  /** Idempotent local ownership release. A failed remote release is not retried. */
  release(options?: ModelCallOptions): Promise<void>;
}

export interface ExpertModelOwner {
  create(input: ModelCreateInput, options?: ModelCallOptions): Promise<ExpertOwnedModel>;
}

/**
 * Adapt the existing typed clients without exposing their branded session-local
 * reference. No account, login, retry or reconnect policy belongs in this adapter.
 */
export function bindExpertModelOwner<ModelReference>(
  client: ExpertModelCreationClient<ModelReference>,
): ExpertModelOwner {
  return Object.freeze({
    async create(
      input: ModelCreateInput,
      options: ModelCallOptions = {},
    ): Promise<ExpertOwnedModel> {
      const made = await client.create(input, options);
      const snapshot = clone(ModelSnapshotSchema, made.snapshot);
      const source = bindExpertModelSource(client, made.model);
      let released = false;
      const owned: ExpertOwnedModel = {
        snapshot,
        source,
        get released() {
          return released;
        },
        async solve(solveInput: ModelSolveInput, solveOptions: ModelCallOptions = {}) {
          if (released) {
            throw new ModelClientError(Code.NotFound, 'model.stale_reference');
          }
          const solved = await client.solve(made.model, solveInput, solveOptions);
          const result = 'result' in solved ? solved.result : solved;
          // Preserve candidates, diagnostics, hints-derived results and every
          // native field. Acceptance is never inferred from phase count here.
          return clone(FullPtResultSchema, result);
        },
        async release(releaseOptions: ModelCallOptions = {}) {
          if (released) return;
          // Both typed clients invalidate the reference at release admission, so
          // retrying a failed release would be semantically wrong.
          released = true;
          await client.release(made.model, releaseOptions);
        },
      };
      return Object.freeze(owned);
    },
  });
}
