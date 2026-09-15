import { fromJson } from '@bufbuild/protobuf';
import { renderToStaticMarkup } from 'react-dom/server';
import { describe, expect, it } from 'vitest';

import { MODEL_VALIDATION_DETAIL_VERSION } from '../api/modelValidationDetail';
import { inspectionFromSnapshot } from '../api/expertModelInspector';
import {
  FullPtResultSchema,
  ModelSnapshotSchema,
} from '../gen/mpmc/model_configuration/v1/model_service_pb';
import { expertResultJson, expertSnapshotJson } from '../test/modelInspectorFixtures';
import { ExpertModelInspector } from './ExpertModelInspector';

describe('Expert model inspector presentation', () => {
  it('renders ordered components, PR76 kij provenance, settings identity and exact applicability presence', () => {
    const inspection = inspectionFromSnapshot(fromJson(ModelSnapshotSchema, expertSnapshotJson()));
    const html = renderToStaticMarkup(<ExpertModelInspector inspection={inspection} />);

    expect(html).toContain('Expert model inspector');
    expect(html).toContain('Read-only PR76 fixture');
    expect(html).not.toContain('Søreide–Whitson');
    expect(html).toContain('0: methane');
    expect(html).toContain('1: ethane');
    expect(html.indexOf('0: methane')).toBeLessThan(html.indexOf('1: ethane'));
    expect(html).toContain('kij methane');
    expect(html).toContain('0.0123');
    expect(html).toContain('Doe 2026');
    expect(html).toContain('mpmc-balanced-default/v1');
    expect(html).toContain('250 K · open');
    expect(html).toContain('30000000 Pa · closed');
    expect(html).toContain('Temperature upper');
    expect(html).toContain('Pressure lower');
    expect(html.match(/unknown/gu)?.length).toBe(2);
    expect(html).not.toContain('Infinity');
    expect(html).not.toContain('unlimited');
  });

  it('renders structured validation and complete-result diagnostics without parsing error text', () => {
    const inspection = inspectionFromSnapshot(fromJson(ModelSnapshotSchema, expertSnapshotJson()));
    const result = fromJson(FullPtResultSchema, expertResultJson());
    const validation = {
      version: MODEL_VALIDATION_DETAIL_VERSION,
      code: 'configuration.invalid_range',
      field: 'definition.applicability.temperature_lower_k',
    } as const;
    const html = renderToStaticMarkup(
      <ExpertModelInspector inspection={inspection} validation={validation} result={result} />,
    );

    expect(html).toContain('Structured service state');
    expect(html).toContain('configuration.invalid_range');
    expect(html).toContain('definition.applicability.temperature_lower_k');
    expect(html).toContain('indeterminate');
    expect(html).toContain('Candidate phases:');
    expect(html).toContain('Opaque backend diagnostic preserved verbatim.');
    expect(html).toContain('Transition evidence:');
  });
});
