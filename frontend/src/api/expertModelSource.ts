import type {
  ExpertDescribeOptions,
  ExpertModelInspection,
  ModelSnapshotDescriber,
} from './expertModelInspector';
import { describeExpertModel } from './expertModelInspector';

/**
 * Read-only closure over one already-owned, session-local model reference.
 *
 * The opaque reference never enters React state/DOM and no wire handle is exposed.
 * The upstream typed client remains the lifetime/session authority; this source
 * deliberately cannot create, mutate, release, reconnect or recover a model.
 */
export interface ExpertModelSource {
  describe(options?: ExpertDescribeOptions): Promise<ExpertModelInspection>;
}

/**
 * Hand an existing typed model reference to the read-only Expert surface without
 * weakening its nominal ownership. Both ModelSessionClient and RendererModelClient
 * satisfy ModelSnapshotDescriber with their own branded reference types.
 */
export function bindExpertModelSource<ModelReference>(
  client: ModelSnapshotDescriber<ModelReference>,
  model: ModelReference,
): ExpertModelSource {
  return Object.freeze({
    describe(options: ExpertDescribeOptions = {}) {
      return describeExpertModel(client, model, options);
    },
  });
}
