'use strict';

const { contextBridge, ipcRenderer } = require('electron');

const bridge = Object.freeze({
  convention: 'MPMC/PT/desktop-bridge/v1',
  discoverPtCapabilities(requestId) {
    return ipcRenderer.invoke('mpmc:pt:discover:v1', requestId);
  },
  solvePtFlash(requestId, request) {
    return ipcRenderer.invoke('mpmc:pt:solve:v1', requestId, request);
  },
  cancel(requestId) {
    ipcRenderer.send('mpmc:pt:cancel:v1', requestId);
  },
});

contextBridge.exposeInMainWorld('mpmcPtDesktop', bridge);

// Production renderer capability for editable models. Main owns the session and
// native model reference; no handle, credential, connect/reconnect or transport
// primitive crosses contextBridge.
const workbenchVersion = 'MPMC/model/workbench-bridge/v1';
const workbenchInvoke = (requestId, operation, fields = {}) => ipcRenderer.invoke(
  'mpmc:model-workbench:invoke:v1',
  { version: workbenchVersion, requestId, operation, ...fields },
);
contextBridge.exposeInMainWorld('mpmcModelWorkbench', Object.freeze({
  convention: workbenchVersion,
  apply: (requestId, input) => workbenchInvoke(requestId, 'apply', { input }),
  solve: (requestId, input) => workbenchInvoke(requestId, 'solve', { input }),
  release: (requestId) => workbenchInvoke(requestId, 'release'),
  cancel: (requestId) => ipcRenderer.send('mpmc:model-workbench:cancel:v1', requestId),
}));

// Legacy low-level model bridges exist only for the dedicated native IPC regression
// harness. Packaged/development product windows do not receive them by default.
const exposeModelDebugBridge = typeof process !== 'undefined'
  && process.env?.MPMC_MODEL_DESKTOP_DEBUG_BRIDGE === '1';
if (exposeModelDebugBridge) {
  const modelVersion = 'MPMC/model/desktop-bridge/v1';
  const modelInvoke = (requestId, operation, fields = {}) => ipcRenderer.invoke(
    'mpmc:model:invoke:v1', { version: modelVersion, requestId, operation, ...fields },
  );
  contextBridge.exposeInMainWorld('mpmcModelDesktop', Object.freeze({
    convention: modelVersion,
    connect: (requestId) => modelInvoke(requestId, 'connect'),
    reconnect: (requestId) => modelInvoke(requestId, 'reconnect'),
    create: (requestId, input) => modelInvoke(requestId, 'create', { input }),
    describe: (requestId, model) => modelInvoke(requestId, 'describe', { model }),
    solve: (requestId, model, input) => modelInvoke(requestId, 'solve', { model, input }),
    release: (requestId, model) => modelInvoke(requestId, 'release', { model }),
    cancel: (requestId) => ipcRenderer.send('mpmc:model:cancel:v1', requestId),
  }));

  // Additive v2 carries sanitized validation details for the low-level regression.
  const modelVersionV2 = 'MPMC/model/desktop-bridge/v2';
  const modelInvokeV2 = (requestId, operation, fields = {}) => ipcRenderer.invoke(
    'mpmc:model:invoke:v2', { version: modelVersionV2, requestId, operation, ...fields },
  );
  contextBridge.exposeInMainWorld('mpmcModelDesktopV2', Object.freeze({
    convention: modelVersionV2,
    connect: (requestId) => modelInvokeV2(requestId, 'connect'),
    reconnect: (requestId) => modelInvokeV2(requestId, 'reconnect'),
    create: (requestId, input) => modelInvokeV2(requestId, 'create', { input }),
    describe: (requestId, model) => modelInvokeV2(requestId, 'describe', { model }),
    solve: (requestId, model, input) => modelInvokeV2(requestId, 'solve', { model, input }),
    release: (requestId, model) => modelInvokeV2(requestId, 'release', { model }),
    cancel: (requestId) => ipcRenderer.send('mpmc:model:cancel:v2', requestId),
  }));
}
