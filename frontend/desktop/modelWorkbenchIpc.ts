import { Code } from '@connectrpc/connect';
import { ipcMain, type IpcMainEvent, type IpcMainInvokeEvent, type WebContents } from 'electron';

import {
  MODEL_WORKBENCH_CANCEL_CHANNEL,
  MODEL_WORKBENCH_CHANNEL,
  MODEL_WORKBENCH_CONVENTION,
} from '../src/api/modelWorkbenchContract';
import type { ModelSessionClient } from '../src/api/modelSessionClient';
import { ModelWorkbenchSession } from './modelWorkbenchSession';

function rejected() {
  return {
    version: MODEL_WORKBENCH_CONVENTION,
    ok: false as const,
    error: { code: Code.PermissionDenied, reason: 'workbench.sender_rejected' },
  };
}

/** Main-process attachment for the renderer-safe editable-model workbench. */
export function registerModelWorkbenchIpc(createClient: () => ModelSessionClient, entryUrl: string) {
  interface Binding {
    owner: ModelWorkbenchSession;
    ready: boolean;
    closing: boolean;
    detach(): void;
  }
  const windows = new Map<WebContents, Binding>();
  const drains = new Set<Promise<void>>();
  let disposed = false;

  function trusted(url: string): boolean {
    try {
      const parsed = new URL(url);
      parsed.hash = '';
      return parsed.href === entryUrl;
    } catch {
      return false;
    }
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

  ipcMain.handle(MODEL_WORKBENCH_CHANNEL, (event, request: unknown) => {
    const binding = authorized(event);
    return binding ? binding.owner.invoke(request) : rejected();
  });
  const cancel = (event: IpcMainEvent, id: unknown) => authorized(event)?.owner.cancel(id);
  ipcMain.on(MODEL_WORKBENCH_CANCEL_CHANNEL, cancel);

  return {
    attach(contents: WebContents): void {
      if (disposed || contents.isDestroyed() || windows.has(contents)) {
        throw new Error('Model workbench IPC attachment rejected.');
      }
      if (windows.size >= 8) throw new Error('Model workbench window limit reached.');
      const owner = new ModelWorkbenchSession(createClient());
      const binding: Binding = { owner, ready: trusted(contents.getURL()), closing: false, detach };

      function navigation(details: { isMainFrame: boolean; isSameDocument: boolean }) {
        if (details.isMainFrame && !details.isSameDocument) {
          binding.ready = false;
          track(owner.reset());
        }
      }
      function committed(
        _event: Electron.Event,
        url: string,
        _status: number,
        _text: string,
        isMainFrame: boolean,
      ) {
        if (isMainFrame && !binding.closing) binding.ready = trusted(url);
      }
      function gone() {
        binding.ready = false;
        track(owner.reset());
      }
      function detach() {
        if (binding.closing) return;
        binding.closing = true;
        binding.ready = false;
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
        ipcMain.removeHandler(MODEL_WORKBENCH_CHANNEL);
        ipcMain.removeListener(MODEL_WORKBENCH_CANCEL_CHANNEL, cancel);
        for (const binding of windows.values()) binding.detach();
      }
      await Promise.all([...drains]);
    },
  };
}
