import { join, resolve } from 'node:path';

import { app, session, type BrowserWindow } from 'electron';

import { registerPtDesktopIpc } from './desktopIpc';
import { PtHostSession } from './hostSession';
import { PtDesktopGateway } from './ptGateway';
import { createPtDesktopWindow, restrictPtDesktopSession } from './window';

function nativeHostPath(): string {
  const developmentOverride = process.env.MPMC_PT_DESKTOP_HOST_BINARY;
  if (!app.isPackaged && developmentOverride !== undefined) {
    return resolve(developmentOverride);
  }
  const executable = process.platform === 'win32'
    ? 'mpmc_pt_service_host.exe'
    : 'mpmc_pt_service_host';
  return join(process.resourcesPath, 'desktop-native', 'bin', executable);
}

const gotSingleInstanceLock = app.requestSingleInstanceLock();
if (!gotSingleInstanceLock) {
  app.quit();
} else {
  const host = new PtHostSession(nativeHostPath(), (line) => {
    console.info(`[pt-host] ${line}`);
  });
  const gateway = new PtDesktopGateway(host);
  let mainWindow: BrowserWindow | null = null;
  let removeIpc: (() => void) | null = null;
  let quitAfterStop = false;

  function createWindow(): void {
    const preload = join(app.getAppPath(), 'preload.cjs');
    mainWindow = createPtDesktopWindow({
      preloadPath: preload,
      showWhenReady: true,
    });
    removeIpc = registerPtDesktopIpc(gateway, () => mainWindow?.webContents ?? null);
    mainWindow.on('closed', () => {
      mainWindow = null;
    });
    void mainWindow.loadFile(join(app.getAppPath(), 'renderer', 'index.html'));
  }

  app.on('second-instance', () => {
    if (mainWindow !== null) {
      if (mainWindow.isMinimized()) {
        mainWindow.restore();
      }
      mainWindow.focus();
    }
  });

  app.on('before-quit', (event) => {
    if (quitAfterStop) {
      return;
    }
    event.preventDefault();
    quitAfterStop = true;
    removeIpc?.();
    void gateway.stop().finally(() => app.quit());
  });

  app.on('window-all-closed', () => app.quit());

  void app.whenReady().then(() => {
    restrictPtDesktopSession(session.defaultSession);
    createWindow();
  });
}
