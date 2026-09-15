import { clone } from '@bufbuild/protobuf';

import {
  ModelSnapshotSchema,
  type ModelSnapshot,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';

export interface ExpertDescribeOptions {
  signal?: AbortSignal;
  timeoutMs?: number;
}

/**
 * Small structural boundary shared by ModelSessionClient (Web) and
 * RendererModelClient (Electron). The inspector never needs a wire handle and
 * never creates or mutates a model.
 */
export interface ModelSnapshotDescriber<ModelReference> {
  describe(
    model: ModelReference,
    options?: ExpertDescribeOptions,
  ): Promise<ModelSnapshot>;
}

export interface ExpertModelInspection {
  /** Exact validated service snapshot retained for read-only presentation. */
  readonly snapshot: ModelSnapshot;
}

/**
 * Isolate the presentation layer from caller mutations without projecting away
 * any service fields. The complete ModelSnapshot remains the source of truth.
 */
export function inspectionFromSnapshot(snapshot: ModelSnapshot): ExpertModelInspection {
  return Object.freeze({ snapshot: clone(ModelSnapshotSchema, snapshot) });
}

/**
 * Read-only describe path used by both browser and Electron clients. No retry,
 * normalization, defaulting or backend interpretation is added here.
 */
export async function describeExpertModel<ModelReference>(
  client: ModelSnapshotDescriber<ModelReference>,
  model: ModelReference,
  options: ExpertDescribeOptions = {},
): Promise<ExpertModelInspection> {
  return inspectionFromSnapshot(await client.describe(model, options));
}
