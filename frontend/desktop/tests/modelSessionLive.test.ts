import { readFileSync } from 'node:fs';
import { fromJson, type JsonValue } from '@bufbuild/protobuf';
import { Code, ConnectError } from '@connectrpc/connect';
import { expect, it } from 'vitest';
import { CreateModelRequestSchema, ModelServiceErrorSchema, SolveModelRequestSchema } from '../../src/gen/mpmc/model_configuration/v1/model_service_pb';
import { MODEL_SESSION_HEADER, MODEL_WIRE_CONTRACT, ModelSessionClient } from '../../src/api/modelSessionClient';
import { PtHostSession } from '../hostSession';
import { createDesktopModelConnection } from '../modelConnection';
import { PtDesktopGateway } from '../ptGateway';

const binary = process.env.MPMC_MODEL_CLIENT_HOST;
const fixture = process.env.MPMC_MODEL_CLIENT_FIXTURE;

it.skipIf(!binary || !fixture)('shared client reconnect/release against the real native host', async () => {
  // C++ exports the existing attributed nitrogen/ethane fixture and PT state.
  const data = JSON.parse(readFileSync(fixture!, 'utf8')) as { create: JsonValue; solve: JsonValue };
  const definition = fromJson(CreateModelRequestSchema, data.create);
  const state = fromJson(SolveModelRequestSchema, data.solve);
  const logs: string[] = [];
  const host = new PtHostSession(binary!, line => logs.push(line), { enableModelSessions: true });
  const gateway = new PtDesktopGateway(host);
  const connection = await host.start();
  const observer = createDesktopModelConnection(connection);
  const observed: { id: string; handle: string }[] = [];
  const closed: boolean[] = [];
  const client = new ModelSessionClient(async signal => {
    const hostConnection = await host.start(); signal.throwIfAborted();
    const transport = createDesktopModelConnection(hostConnection);
    const entry = { id: '', handle: '' }; observed.push(entry);
    const index = closed.length; closed.push(false);
    return {
      ...transport,
      sessions: {
        async *openModelSession(request, options) {
          for await (const opened of transport.sessions.openModelSession(request, options)) {
            entry.id = opened.sessionId!; yield opened;
          }
        },
      },
      models: {
        ...transport.models,
        async createModel(request, options) {
          const made = await transport.models.createModel(request, options);
          entry.handle = made.modelHandle!; return made;
        },
      },
      close() { closed[index] = true; transport.close(); },
    };
  });
  function options(id: string) {
    const headers = new Headers(observer.headers); headers.set(MODEL_SESSION_HEADER, id);
    return { headers, timeoutMs: 5_000 };
  }
  async function gone(entry: { id: string; handle: string }) {
    const end = Date.now() + 8_000;
    while (true) {
      try {
        await observer.models.describeModel({ wireContract: MODEL_WIRE_CONTRACT, modelHandle: entry.handle }, options(entry.id));
      } catch (cause) {
        const error = ConnectError.from(cause);
        expect(error.code).toBe(Code.NotFound);
        expect(error.findDetails(ModelServiceErrorSchema)).toMatchObject([
          { wireContract: MODEL_WIRE_CONTRACT, code: 'session.not_found' },
        ]);
        return;
      }
      expect(Date.now(), 'server retained a disconnected model/session').toBeLessThan(end);
      await new Promise(resolve => setTimeout(resolve, 10));
    }
  }
  try {
    await client.connect();
    const initial = await client.create(definition);
    const solved = await client.solve(initial.model, state);
    expect(solved.pressurePa).toBe(state.pressurePa);
    expect(solved.temperatureK).toBe(state.temperatureK);
    expect(solved.feed).toEqual(state.feed);
    expect(solved.candidatePhaseSet?.phases).toHaveLength(2);
    expect(await client.describe(initial.model)).toEqual(initial.snapshot);
    // Four models fill the actual per-session host cap. Release restores a slot.
    const others = [];
    for (let i = 0; i < 3; ++i) others.push(await client.create(definition));
    const releasedHandle = observed.at(-1)!.handle;
    await client.release(others.at(-1)!.model);
    const rejected = observer.models.describeModel({ wireContract: MODEL_WIRE_CONTRACT, modelHandle: releasedHandle },
      options(observed.at(-1)!.id));
    await expect(rejected).rejects.toMatchObject({ code: Code.NotFound });
    await rejected.catch(cause => {
      expect(ConnectError.from(cause).findDetails(ModelServiceErrorSchema)).toMatchObject([
        { wireContract: MODEL_WIRE_CONTRACT, code: 'registry.model_not_found' },
      ]);
    });
    await client.create(definition);
    let old = initial.model;
    // More reconnects than the host's 16-session cap expose unreclaimed slots.
    for (let i = 0; i < 20; ++i) {
      const previous = { ...observed.at(-1)! };
      await client.reconnect(); await gone(previous);
      await expect(client.solve(old, state)).rejects.toMatchObject({ reason: 'model.stale_reference' });
      await expect(observer.models.describeModel({ wireContract: MODEL_WIRE_CONTRACT, modelHandle: previous.handle },
        options(observed.at(-1)!.id))).rejects.toMatchObject({ code: Code.NotFound });
      const made = await client.create(definition);
      expect(await client.solve(made.model, state)).toEqual(solved);
      old = made.model;
    }
    const last = { ...observed.at(-1)! };
    await client.dispose(); await client.dispose(); await gone(last);
    expect(closed.every(Boolean)).toBe(true);
    await expect(client.reconnect()).rejects.toMatchObject({ reason: 'client.disposed' });

    // Exercise the actual desktop gateway composition, still without renderer IPC.
    await gateway.models.connect(); const made = await gateway.models.create(definition);
    expect(await gateway.models.solve(made.model, state)).toEqual(solved);
    await gateway.models.reconnect();
    await expect(gateway.models.describe(made.model)).rejects.toMatchObject({ code: Code.NotFound });
    await gateway.models.create(definition); // Owned until gateway shutdown.
    await gateway.stop();
    expect(gateway.models.connected).toBe(false);
    expect(logs.some(line => line.includes('"event":"pt_process_stopped"'))).toBe(true);
    for (const entry of observed) {
      expect(logs.join('\n')).not.toContain(entry.id);
      expect(logs.join('\n')).not.toContain(entry.handle);
    }
    expect(logs.join('\n')).not.toContain(connection.bearerToken);
    console.info('MODEL_CLIENT_LIVE_OK create/solve/release, 20 reconnects, stale handles, disposal and gateway shutdown');
  } finally {
    observer.close(); await client.dispose(); await gateway.stop();
  }
}, 90_000);
