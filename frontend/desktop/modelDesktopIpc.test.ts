import { EventEmitter } from 'node:events';
import { Code } from '@connectrpc/connect';
import type { WebContents } from 'electron';
import { expect, it, vi } from 'vitest';
import { MODEL_DESKTOP_CHANNEL, MODEL_DESKTOP_CANCEL_CHANNEL, MODEL_DESKTOP_CONVENTION, MODEL_DESKTOP_V2_CHANNEL, MODEL_DESKTOP_V2_CANCEL_CHANNEL, MODEL_DESKTOP_V2_CONVENTION, type ModelDesktopReply } from '../src/api/modelDesktopContract';
import { ModelSessionClient } from '../src/api/modelSessionClient';

const main = vi.hoisted(() => ({
  handlers: new Map<string, (event: unknown, request: unknown) => Promise<ModelDesktopReply>>(),
  listeners: new Map<string, (event: unknown, request: unknown) => void>(),
  handle: vi.fn(), on: vi.fn(), removeHandler: vi.fn(), removeListener: vi.fn(),
}));
vi.mock('electron', () => ({ ipcMain: main }));
import { registerModelDesktopIpc } from './modelDesktopIpc';

const url = 'file:///packaged/renderer/index.html';
class Contents extends EventEmitter {
  mainFrame = { url };
  destroyed = false;
  isDestroyed() { return this.destroyed; }
  getURL() { return this.mainFrame.url; }
  get native() { return this as unknown as WebContents; }
  event(frame: unknown = this.mainFrame) { return { sender: this, senderFrame: frame }; }
}
function setup() {
  main.handlers.clear(); main.listeners.clear();
  main.handle.mockImplementation((channel, handler) => main.handlers.set(channel, handler));
  main.on.mockImplementation((channel, handler) => main.listeners.set(channel, handler));
  main.removeHandler.mockImplementation(channel => main.handlers.delete(channel));
  main.removeListener.mockImplementation(channel => main.listeners.delete(channel));
  const clients: ModelSessionClient[] = [];
  const ipc = registerModelDesktopIpc(() => {
    const client = new ModelSessionClient(async () => { throw new Error('No native connection expected.'); });
    vi.spyOn(client, 'disconnect'); vi.spyOn(client, 'dispose');
    clients.push(client); return client;
  }, url);
  const invoke = (contents: Contents, frame: unknown = contents.mainFrame) => main.handlers.get(MODEL_DESKTOP_CHANNEL)!(
    contents.event(frame), { version: 'unknown', operation: 'connect', requestId: 'probe' },
  );
  return { ipc, clients, invoke };
}

it('authorizes only attached current main frames at the exact packaged URL', async () => {
  const { ipc, clients, invoke } = setup(); const a = new Contents(); const rogue = new Contents();
  ipc.attach(a.native);
  try {
    expect(await invoke(a)).toMatchObject({ error: { code: Code.InvalidArgument, reason: 'ipc.unsupported_version' } });
    for (const reply of [await invoke(rogue), await invoke(a, null), await invoke(a, { url }), await invoke(a, { url, parent: a.mainFrame })]) {
      expect(reply).toMatchObject({ error: { code: Code.PermissionDenied } });
    }
    for (const other of ['https://untrusted.invalid/', `${url}?injected=1`, 'file:///packaged/renderer/other.html']) {
      a.mainFrame.url = other;
      expect(await invoke(a)).toMatchObject({ error: { code: Code.PermissionDenied } });
    }
    a.mainFrame.url = `${url}#section`;
    expect(await invoke(a)).toMatchObject({ error: { code: Code.InvalidArgument } });
    expect(clients[0]!.connected).toBe(false);
  } finally { await ipc.dispose(); }
});

it('revokes on main-document navigation, renderer loss and destruction; keeps another window independent', async () => {
  const { ipc, clients, invoke } = setup(); const a = new Contents(); const b = new Contents();
  ipc.attach(a.native); ipc.attach(b.native);
  try {
    a.emit('did-start-navigation', { isMainFrame: false, isSameDocument: false });
    a.emit('did-start-navigation', { isMainFrame: true, isSameDocument: true });
    expect(clients[0]!.disconnect).not.toHaveBeenCalled();
    a.emit('did-start-navigation', { isMainFrame: true, isSameDocument: false });
    expect(clients[0]!.disconnect).toHaveBeenCalledTimes(1);
    expect(await invoke(a)).toMatchObject({ error: { code: Code.PermissionDenied } });
    expect(await invoke(b)).toMatchObject({ error: { code: Code.InvalidArgument } });
    a.emit('did-frame-navigate', {}, url, 200, 'OK', true);
    expect(await invoke(a)).toMatchObject({ error: { code: Code.InvalidArgument } });
    a.emit('render-process-gone', {}, { reason: 'crashed' });
    expect(await invoke(a)).toMatchObject({ error: { code: Code.PermissionDenied } });
    a.destroyed = true; a.emit('destroyed');
    expect(clients[0]!.disconnect).toHaveBeenCalledTimes(3);
    expect(clients[1]!.disconnect).not.toHaveBeenCalled();
    expect(a.listenerCount('did-start-navigation')).toBe(0);
    expect(a.listenerCount('render-process-gone')).toBe(0);
  } finally { await ipc.dispose(); await ipc.dispose(); }
  expect(main.handlers.has(MODEL_DESKTOP_CHANNEL)).toBe(false);
  expect(main.handlers.has(MODEL_DESKTOP_V2_CHANNEL)).toBe(false);
  expect(main.listeners.has(MODEL_DESKTOP_V2_CANCEL_CHANNEL)).toBe(false);
  expect(main.listeners.has(MODEL_DESKTOP_CANCEL_CHANNEL)).toBe(false);
  expect(clients.every(client => !client.connected)).toBe(true);
});

it('bounds attached and draining windows and refuses attachment after shutdown', async () => {
  const { ipc, clients } = setup(); const windows = Array.from({ length: 8 }, () => new Contents());
  for (const window of windows) ipc.attach(window.native);
  expect(() => ipc.attach(new Contents().native)).toThrow('limit');
  let release!: () => void;
  vi.mocked(clients[0]!.disconnect).mockImplementationOnce(() => new Promise<void>(yes => { release = yes; }));
  windows[0]!.destroyed = true; windows[0]!.emit('destroyed');
  expect(() => ipc.attach(new Contents().native)).toThrow('limit');
  release();
  await vi.waitFor(() => expect(clients[0]!.dispose).toHaveBeenCalled());
  await vi.waitFor(() => ipc.attach(new Contents().native));
  await ipc.dispose();
  expect(() => ipc.attach(new Contents().native)).toThrow('rejected');
});

it('rejects malformed calls from a trusted sender without exposing a generic Electron API', async () => {
  const { ipc } = setup(); const a = new Contents(); ipc.attach(a.native);
  try {
    const reply = await main.handlers.get(MODEL_DESKTOP_CHANNEL)!(a.event(), {
      version: MODEL_DESKTOP_CONVENTION, operation: 'arbitrary-channel', requestId: 'probe',
    });
    expect(reply).toMatchObject({ error: { code: Code.InvalidArgument, reason: 'ipc.unknown_operation' } });
    expect([...main.handlers.keys()]).toEqual([MODEL_DESKTOP_CHANNEL, MODEL_DESKTOP_V2_CHANNEL]);
  } finally { await ipc.dispose(); }
});

it('applies the same sender and lifecycle authorization to both protocol channels', async () => {
  const { ipc, clients } = setup(); const a = new Contents(); const rogue = new Contents(); ipc.attach(a.native);
  try {
    for (const [channel, version] of [[MODEL_DESKTOP_CHANNEL, MODEL_DESKTOP_CONVENTION], [MODEL_DESKTOP_V2_CHANNEL, MODEL_DESKTOP_V2_CONVENTION]]) {
      const invoke = main.handlers.get(channel!)!;
      const request = { version, operation: 'connect', requestId: 'x' };
      for (const event of [rogue.event(), a.event(null), a.event({ url })]) {
        expect(await invoke(event, request)).toEqual({ version, ok: false, error: { code: Code.PermissionDenied, reason: 'ipc.sender_rejected' } });
      }
    }
    expect(clients).toHaveLength(1);
    a.emit('render-process-gone');
    expect(await main.handlers.get(MODEL_DESKTOP_V2_CHANNEL)!(a.event(), { version: MODEL_DESKTOP_V2_CONVENTION }))
      .toMatchObject({ error: { code: Code.PermissionDenied } });
    main.listeners.get(MODEL_DESKTOP_V2_CANCEL_CHANNEL)!(rogue.event(), 'x');
    expect(clients[0]!.disconnect).toHaveBeenCalledTimes(1);
  } finally { await ipc.dispose(); }
});
