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
