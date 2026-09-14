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


// Fixed channels/methods only. Never expose ipcRenderer or host credentials.
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

// Additive v2 carries sanitized validation details; v1 remains byte-shape compatible.
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
