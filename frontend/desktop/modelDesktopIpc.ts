import { Code } from '@connectrpc/connect';
import { ipcMain, type IpcMainEvent, type IpcMainInvokeEvent, type WebContents } from 'electron';
import { MODEL_DESKTOP_CONVENTION, MODEL_DESKTOP_V2_CONVENTION, MODEL_DESKTOP_CANCEL_CHANNEL, MODEL_DESKTOP_CHANNEL, MODEL_DESKTOP_V2_CHANNEL, MODEL_DESKTOP_V2_CANCEL_CHANNEL } from '../src/api/modelDesktopContract';
import type { ModelSessionClient } from '../src/api/modelSessionClient';
import { ModelDesktopSession, modelDesktopFailure } from './modelDesktopSession';

/** Explicit main-process attachment; renderer IDs/URLs never grant ownership. */
export function registerModelDesktopIpc(createClient: () => ModelSessionClient, entryUrl: string) {
  interface Binding { owner: ModelDesktopSession; ready: boolean; closing: boolean; detach(): void }
  const windows = new Map<WebContents, Binding>();
  const drains = new Set<Promise<void>>();
  let disposed = false;
  function trusted(url: string): boolean {
    try { const parsed = new URL(url); parsed.hash = ''; return parsed.href === entryUrl; }
    catch { return false; }
  }
  function authorized(event: IpcMainEvent | IpcMainInvokeEvent): Binding | undefined {
    const binding = windows.get(event.sender);
    if (disposed || !binding?.ready || binding.closing || event.sender.isDestroyed()) return undefined;
    const frame = event.senderFrame;
    if (!frame || frame !== event.sender.mainFrame || !trusted(frame.url)) return undefined;
    return binding;
  }
  function track(work: Promise<void>): void {
    if (drains.has(work)) return;
    drains.add(work);
    void work.finally(() => drains.delete(work)).catch(() => {});
  }
  for (const [channel, version] of [[MODEL_DESKTOP_CHANNEL, MODEL_DESKTOP_CONVENTION],
    [MODEL_DESKTOP_V2_CHANNEL, MODEL_DESKTOP_V2_CONVENTION]] as const) {
    ipcMain.handle(channel, (event, request: unknown) => {
      const binding = authorized(event);
      return binding ? binding.owner.invoke(request, version)
        : modelDesktopFailure(Code.PermissionDenied, 'ipc.sender_rejected', version);
    });
  }
  const cancel = (event: IpcMainEvent, id: unknown) => authorized(event)?.owner.cancel(id);
  ipcMain.on(MODEL_DESKTOP_CANCEL_CHANNEL, cancel);
  ipcMain.on(MODEL_DESKTOP_V2_CANCEL_CHANNEL, cancel);
  return {
    attach(contents: WebContents): void {
      if (disposed || contents.isDestroyed() || windows.has(contents)) throw new Error('Model IPC attachment rejected.');
      // Closing bindings retain their reservation until their client drains.
      if (windows.size >= 8) throw new Error('Model IPC window limit reached.');
      const owner = new ModelDesktopSession(createClient());
      const binding: Binding = { owner, ready: trusted(contents.getURL()), closing: false, detach };
      function navigation(details: { isMainFrame: boolean; isSameDocument: boolean }) {
        if (details.isMainFrame && !details.isSameDocument) {
          binding.ready = false; track(owner.reset());
        }
      }
      function committed(_event: Electron.Event, url: string, _status: number, _text: string, isMainFrame: boolean) {
        if (isMainFrame && !binding.closing) binding.ready = trusted(url);
      }
      function gone() { binding.ready = false; track(owner.reset()); }
      function detach() {
        if (binding.closing) return;
        binding.closing = true; binding.ready = false;
        contents.removeListener('did-start-navigation', navigation);
        contents.removeListener('did-frame-navigate', committed);
        contents.removeListener('render-process-gone', gone);
        contents.removeListener('destroyed', detach);
        track(owner.dispose().finally(() => windows.delete(contents)));
      }
      contents.on('did-start-navigation', navigation);
      contents.on('did-frame-navigate', committed);
      contents.on('render-process-gone', gone);
      contents.once('destroyed', detach);
      windows.set(contents, binding);
    },
    async dispose(): Promise<void> {
      if (!disposed) {
        disposed = true;
        ipcMain.removeHandler(MODEL_DESKTOP_CHANNEL);
        ipcMain.removeHandler(MODEL_DESKTOP_V2_CHANNEL);
        ipcMain.removeListener(MODEL_DESKTOP_CANCEL_CHANNEL, cancel);
        ipcMain.removeListener(MODEL_DESKTOP_V2_CANCEL_CHANNEL, cancel);
        for (const binding of windows.values()) binding.detach();
      }
      await Promise.all([...drains]);
    },
  };
}
