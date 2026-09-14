import { join, resolve } from 'node:path';

import { app, session, type BrowserWindow } from 'electron';

import { registerPtDesktopIpc } from './desktopIpc';
import { PtHostSession } from './hostSession';
import {
  installedSmokeRequested,
  runInstalledDesktopSmoke,
} from './installSmoke';
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
  const installSmoke = installedSmokeRequested(process.argv);
  let mainWindow: BrowserWindow | null = null;
  let removeIpc: (() => void) | null = null;
  let quitAfterStop = false;

  async function stopAndExit(exitCode: number): Promise<void> {
    quitAfterStop = true;
    removeIpc?.();
    removeIpc = null;
    if (mainWindow !== null && !mainWindow.isDestroyed()) {
      mainWindow.destroy();
    }
    mainWindow = null;
    try {
      await gateway.stop();
    } finally {
      app.exit(exitCode);
    }
  }

  async function createWindow(): Promise<void> {
    const preload = join(app.getAppPath(), 'preload.cjs');
    mainWindow = createPtDesktopWindow({
      preloadPath: preload,
      showWhenReady: !installSmoke,
    });
    removeIpc = registerPtDesktopIpc(gateway, () => mainWindow?.webContents ?? null);
    mainWindow.on('closed', () => {
      mainWindow = null;
    });
    await mainWindow.loadFile(join(app.getAppPath(), 'renderer', 'index.html'));

    if (!installSmoke) {
      return;
    }

    const planPath = process.env.MPMC_PT_DESKTOP_INSTALL_SMOKE_PLAN;
    const resultPath = process.env.MPMC_PT_DESKTOP_INSTALL_SMOKE_RESULT;
    if (planPath === undefined || resultPath === undefined) {
      throw new Error(
        'Installed desktop smoke requires absolute plan and result paths.',
      );
    }
    await runInstalledDesktopSmoke(mainWindow, planPath, resultPath);
    await stopAndExit(0);
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

  void app.whenReady().then(async () => {
    restrictPtDesktopSession(session.defaultSession);
    try {
      await createWindow();
    } catch (cause) {
      console.error(cause instanceof Error ? cause.stack ?? cause.message : String(cause));
      await stopAndExit(1);
    }
  });
}
