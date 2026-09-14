import {
  BrowserWindow,
  type BrowserWindowConstructorOptions,
  type Session,
} from 'electron';

export interface PtDesktopWindowOptions {
  readonly preloadPath: string;
  readonly showWhenReady: boolean;
}

export function restrictPtDesktopSession(target: Session): void {
  target.setPermissionCheckHandler(() => false);
  target.setPermissionRequestHandler((_webContents, _permission, callback) => {
    callback(false);
  });
  target.webRequest.onBeforeRequest((details, callback) => {
    const allowed = details.url.startsWith('file:') || details.url.startsWith('devtools:');
    callback({ cancel: !allowed });
  });
}

export function createPtDesktopWindow(options: PtDesktopWindowOptions): BrowserWindow {
  const browserOptions: BrowserWindowConstructorOptions = {
    width: 1280,
    height: 820,
    minWidth: 900,
    minHeight: 640,
    show: false,
    webPreferences: {
      preload: options.preloadPath,
      contextIsolation: true,
      nodeIntegration: false,
      sandbox: true,
      webSecurity: true,
      allowRunningInsecureContent: false,
    },
  };
  const window = new BrowserWindow(browserOptions);
  window.webContents.setWindowOpenHandler(() => ({ action: 'deny' }));
  window.webContents.on('will-navigate', (event) => event.preventDefault());
  window.webContents.on('will-attach-webview', (event) => event.preventDefault());
  if (options.showWhenReady) {
    window.once('ready-to-show', () => window.show());
  }
  return window;
}
