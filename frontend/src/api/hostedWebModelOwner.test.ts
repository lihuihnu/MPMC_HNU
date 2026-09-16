import { create } from '@bufbuild/protobuf';
import { describe, expect, it, vi } from 'vitest';

import {
  ModelSnapshotSchema,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';
import {
  MODEL_SESSION_CONTRACT,
  MODEL_WIRE_CONTRACT,
  ModelClientError,
  ModelSessionClient,
  type ModelCreateInput,
  type ModelSessionConnector,
} from './modelSessionClient';
import {
  HOSTED_WEB_MODEL_OWNERSHIP_CONVENTION,
  hostedWebExpertModelOwner,
  type HostedWebModelOwnershipCapability,
} from './hostedWebModelOwner';

function capability(): HostedWebModelOwnershipCapability {
  return {
    convention: HOSTED_WEB_MODEL_OWNERSHIP_CONVENTION,
    baseUrl: '/model-api',
    identity: {
      getAccessToken: vi.fn(async () => ({ accessToken: 'opaque.web.test.token' })),
    },
  };
}

function liveSession(order: string[]) {
  const createModel = vi.fn(async () => {
    order.push('create');
    return {
      wireContract: MODEL_WIRE_CONTRACT,
      modelHandle: 'mh1_web_test',
      snapshot: create(ModelSnapshotSchema, {}),
    };
  });
  const releaseModel = vi.fn(async () => {
    order.push('release');
    return { wireContract: MODEL_WIRE_CONTRACT };
  });
  const connector: ModelSessionConnector = async () => ({
    models: { createModel, releaseModel } as never,
    sessions: {
      openModelSession: vi.fn(async function* (
        _request: unknown,
        options: { signal?: AbortSignal },
      ) {
        order.push('open');
        yield { wireContract: MODEL_SESSION_CONTRACT, sessionId: 'ms1_web_test' };
        await new Promise<void>((resolve) => {
          options.signal?.addEventListener('abort', () => resolve(), { once: true });
        });
      }),
    } as never,
    close: vi.fn(),
  });
  return { client: new ModelSessionClient(connector), createModel, releaseModel };
}

describe('hosted Web Expert model owner', () => {
  it('fails closed unless the trusted host supplies the exact runtime capability', async () => {
    expect(hostedWebExpertModelOwner(undefined)).toBeNull();
    expect(hostedWebExpertModelOwner({
      convention: 'MPMC/model/hosted-web-ownership/v0',
      baseUrl: '/model-api',
      identity: { getAccessToken: async () => ({ accessToken: 'opaque' }) },
    })).toBeNull();
    expect(hostedWebExpertModelOwner({
      convention: HOSTED_WEB_MODEL_OWNERSHIP_CONVENTION,
      baseUrl: '/model-api',
      identity: { getAccessToken: async () => ({ accessToken: 'opaque' }) },
      accessToken: 'must-not-be-accepted-as-runtime-capability-shape',
    })).toBeNull();

    const owner = hostedWebExpertModelOwner(capability());
    expect(owner).not.toBeNull();
    await owner?.dispose();
  });

  it('opens the authenticated session before the first typed create and reuses the shared owner', async () => {
    const order: string[] = [];
    const { client, createModel, releaseModel } = liveSession(order);
    const owner = hostedWebExpertModelOwner(capability(), () => client);
    expect(owner).not.toBeNull();

    const model = await owner!.create({} as ModelCreateInput);
    expect(order.slice(0, 2)).toEqual(['open', 'create']);
    expect(createModel).toHaveBeenCalledTimes(1);
    await model.release();
    expect(releaseModel).toHaveBeenCalledTimes(1);
    await owner!.dispose();
  });

  it('cancels a still-opening session instead of dispatching Create after the caller aborts', async () => {
    const createModel = vi.fn();
    const connector: ModelSessionConnector = async () => ({
      models: { createModel } as never,
      sessions: {
        openModelSession: vi.fn(async function* (
          _request: unknown,
          options: { signal?: AbortSignal },
        ) {
          await new Promise<void>((resolve) => {
            options.signal?.addEventListener('abort', () => resolve(), { once: true });
          });
        }),
      } as never,
      close: vi.fn(),
    });
    const client = new ModelSessionClient(connector);
    const owner = hostedWebExpertModelOwner(capability(), () => client)!;
    const controller = new AbortController();
    const pending = owner.create({} as ModelCreateInput, { signal: controller.signal });
    controller.abort();

    await expect(pending).rejects.toMatchObject<ModelClientError>({
      reason: 'web.owner_cancelled',
    });
    expect(createModel).not.toHaveBeenCalled();
    await owner.dispose();
  });
});
