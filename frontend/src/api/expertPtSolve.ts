import { Code } from '@connectrpc/connect';

import type { ExpertOwnedModel } from './expertModelOwner';
import { ModelClientError, type ModelCallOptions, type ModelSolveInput } from './modelSessionClient';
import { MODEL_VALIDATION_DETAIL_VERSION } from './modelValidationDetail';
import type { FullPtResult, ModelSnapshot } from '../gen/mpmc/model_configuration/v1/model_service_pb';

export interface ExpertPtInput {
  readonly pressurePa: string;
  readonly temperatureK: string;
  readonly feed: ReadonlyArray<{ readonly componentId: string; readonly fraction: string }>;
}

function invalid(field: string): never {
  throw new ModelClientError(Code.InvalidArgument, 'expert.invalid_pt_input', {
    version: MODEL_VALIDATION_DETAIL_VERSION,
    code: 'request.rejected',
    field,
  });
}

function componentIds(snapshot: ModelSnapshot): string[] {
  const ids = snapshot.definition?.components.map(component => component.componentId) ?? [];
  if (!ids.length || ids.some(id => !id) || new Set(ids).size !== ids.length) {
    throw new ModelClientError(Code.DataLoss, 'expert.invalid_model_components');
  }
  return ids as string[];
}

/** A new applied revision starts with missing feed, never guessed fractions. */
export function emptyExpertPtInput(snapshot: ModelSnapshot): ExpertPtInput {
  return {
    pressurePa: '',
    temperatureK: '',
    feed: componentIds(snapshot).map(componentId => ({ componentId, fraction: '' })),
  };
}

function number(text: string, field: string): number {
  const value = text.trim();
  // Decimal/scientific notation only: Number('') and Number('0x10') are not input.
  if (!/^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/iu.test(value)) invalid(field);
  const parsed = Number(value);
  if (!Number.isFinite(parsed)) invalid(field);
  return parsed;
}

/** No EOS, normalization tolerance or implicit unit conversion in presentation. */
export function buildExpertPtRequest(snapshot: ModelSnapshot, input: ExpertPtInput): ModelSolveInput {
  const ids = componentIds(snapshot);
  const pressurePa = number(input.pressurePa, 'pressure_pa');
  const temperatureK = number(input.temperatureK, 'temperature_k');
  if (pressurePa <= 0) invalid('pressure_pa');
  if (temperatureK <= 0) invalid('temperature_k');
  if (input.feed.length !== ids.length) invalid('feed');
  const feed = input.feed.map((entry, index) => {
    if (entry.componentId !== ids[index]) invalid(`feed[${index}]`);
    const fraction = number(entry.fraction, `feed[${index}]`);
    if (fraction < 0 || fraction > 1) invalid(`feed[${index}]`);
    return fraction;
  });
  // The native model decides normalization tolerance, applicability and support.
  // Send the entered values unchanged rather than silently repairing them.
  return { pressurePa, temperatureK, feed };
}

/**
 * One admitted calculation, no queue or retry. Input/model edits invalidate only
 * publication: cancelling a typed transport can invalidate its whole model lease.
 * The existing backend deadline bounds work; the native lease owns admitted work.
 */
export class ExpertPtSolveGate {
  #revision = 0;
  #busy = false;
  get busy(): boolean { return this.#busy; }

  invalidate(): void { ++this.#revision; }

  async solve(
    model: Pick<ExpertOwnedModel, 'solve'>,
    input: ModelSolveInput,
    options: ModelCallOptions = {},
  ): Promise<FullPtResult | null> {
    if (this.#busy) throw new ModelClientError(Code.ResourceExhausted, 'expert.solve_busy');
    const revision = this.#revision;
    this.#busy = true;
    try {
      const result = await model.solve(input, options);
      return revision === this.#revision ? result : null;
    } catch (cause) {
      if (revision !== this.#revision) return null;
      throw cause;
    } finally {
      this.#busy = false;
    }
  }
}
