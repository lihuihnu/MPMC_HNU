import { clone } from '@bufbuild/protobuf';

import type {
  ModelCallOptions,
  ModelCreateInput,
} from './modelSessionClient';
import type { ExpertModelSource } from './expertModelSource';
import { bindExpertModelSource } from './expertModelSource';
import type { ModelSnapshotDescriber } from './expertModelInspector';
import {
  ModelSnapshotSchema,
  type ModelSnapshot,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';

export interface ExpertModelCreationClient<ModelReference>
  extends ModelSnapshotDescriber<ModelReference> {
  create(
    input: ModelCreateInput,
    options?: ModelCallOptions,
  ): Promise<{ model: ModelReference; snapshot: ModelSnapshot }>;
  release(model: ModelReference, options?: ModelCallOptions): Promise<void>;
}

export interface ExpertOwnedModel {
  readonly snapshot: ModelSnapshot;
  readonly source: ExpertModelSource;
  readonly released: boolean;
  /** Idempotent local ownership release. A failed remote release is not retried. */
  release(options?: ModelCallOptions): Promise<void>;
}

export interface ExpertModelOwner {
  create(input: ModelCreateInput, options?: ModelCallOptions): Promise<ExpertOwnedModel>;
}

/**
 * Adapt either typed Web or Electron model ownership without exposing its branded
 * session-local reference. Creation remains a mutation on the existing client;
 * this adapter performs no retry or reconnect.
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
