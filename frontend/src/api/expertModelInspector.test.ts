import { equals, fromJson, toJson } from '@bufbuild/protobuf';
import { describe, expect, expectTypeOf, it, vi } from 'vitest';

import type { ModelReference, ModelSessionClient } from './modelSessionClient';
import type { RendererModelClient, RendererModelReference } from './rendererModelClient';
import { readModelSnapshot } from './rendererModelWire';
import {
  describeExpertModel,
  inspectionFromSnapshot,
  type ModelSnapshotDescriber,
} from './expertModelInspector';
import { ModelSnapshotSchema } from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { expertSnapshotJson } from '../test/modelInspectorFixtures';

describe('read-only Expert snapshot boundary', () => {
  it('accepts both existing typed client describe shapes without a new transport abstraction', () => {
    expectTypeOf<ModelSessionClient>().toMatchTypeOf<ModelSnapshotDescriber<ModelReference>>();
    expectTypeOf<RendererModelClient>().toMatchTypeOf<ModelSnapshotDescriber<RendererModelReference>>();
  });

  it('retains the complete Web service snapshot and isolates caller mutation', async () => {
    const snapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
    const describe = vi.fn(async () => snapshot);
    const inspection = await describeExpertModel({ describe }, Symbol('web-model'));

    expect(describe).toHaveBeenCalledTimes(1);
    expect(equals(ModelSnapshotSchema, inspection.snapshot, snapshot)).toBe(true);
    expect(toJson(ModelSnapshotSchema, inspection.snapshot)).toEqual(toJson(ModelSnapshotSchema, snapshot));

    snapshot.definition!.displayName = 'mutated by caller';
    expect(inspection.snapshot.definition?.displayName).toBe('Read-only PR76 fixture');
  });

  it('retains the complete Electron renderer snapshot after its contract decoder', async () => {
    const serviceSnapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
    const rendererSnapshot = readModelSnapshot(toJson(ModelSnapshotSchema, serviceSnapshot));
    const inspection = await describeExpertModel(
      { describe: vi.fn(async () => rendererSnapshot) },
      Symbol('electron-model'),
    );

    expect(equals(ModelSnapshotSchema, inspection.snapshot, serviceSnapshot)).toBe(true);
    expect(toJson(ModelSnapshotSchema, inspection.snapshot)).toEqual(
      toJson(ModelSnapshotSchema, serviceSnapshot),
    );
  });

  it('preserves one-sided applicability presence instead of inventing an opposite endpoint', () => {
    const snapshot = fromJson(ModelSnapshotSchema, expertSnapshotJson());
    const inspection = inspectionFromSnapshot(snapshot);
    expect(inspection.snapshot.definition?.applicability?.temperatureLowerK).toBe(250);
    expect(inspection.snapshot.definition?.applicability?.temperatureUpperK).toBeUndefined();
    expect(inspection.snapshot.definition?.applicability?.pressureLowerPa).toBeUndefined();
    expect(inspection.snapshot.definition?.applicability?.pressureUpperPa).toBe(30_000_000);
  });
});
