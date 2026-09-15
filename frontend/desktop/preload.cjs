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

// Product renderer capability for editable models. Main owns the session and
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
