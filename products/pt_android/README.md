# PT Android Product Shell v1

This directory is an isolated Android product path for the existing C++ PT
runtime. It does **not** change the Windows/Linux/macOS products, Electron
process, gRPC transport, EOS/flash equations, parameter values, solver
thresholds, or publication rules.

## Native boundary

`libmpmc_pt_android_core.so` is built with the Android NDK for `arm64-v8a` and
`x86_64`. It links the model-neutral `mpmc::runtime` chain and directly compiles
the existing repository-curated PT snapshot/backend assembly translation units:

- `modules/pt_process/src/configured_backends.cpp`;
- `modules/pt_process/src/parameter_snapshot_supplier.cpp`.

Android does **not** link `runtime_grpc` or the process host. A product-private
`composition_root.hpp` shim mirrors only `OwnedConfiguredPtBackend` and its
aliasing-owner helper. This keeps the scientific parameter source singular while
leaving the desktop/process composition root untouched.

The native library now exposes three JNI surfaces:

1. the existing repository-curated emulator smoke entry point;
2. model-neutral capability discovery as versioned JSON; and
3. model-neutral PT solve as versioned JSON.

The product bridge convention is `MPMC/PT/android-product-bridge/v1`. The JNI
bridge owns no PR76, SW92, or CPA special cases: it constructs the existing
`PtService`, transports the service contracts, and preserves service-side
validation and scientific outcomes.

## Product Shell v1

`products/pt_android/shell/` wraps the proven JNI runtime with an isolated
Capacitor shell while reusing the existing React/Vite PT UI directly:

```text
existing frontend/src/App + styles
        ↓
AndroidFlashClient (FlashClient)
        ↓
Capacitor app-local MpmcPt plugin
        ↓
JNI discoverJson / solveJson
        ↓
libmpmc_pt_android_core.so
        ↓
PtService
        ↓
PR76 / SW92 / CPA
```

The shell pins Capacitor `8.5.2` and React `19.2.3`. The generated Android Studio
project is deliberately **not** committed: CI generates it from the pinned
Capacitor input, sets the native compatibility floor to API 26, overlays the
small Java plugin, embeds the validated x86_64 native core, and builds a debug
engineering APK. This keeps generated Gradle/template churn out of the main
repository.

The app-local Java plugin is transport only. It runs native work on one dedicated
executor thread and contains no EOS equations, parameter data, phase-selection
logic, retry policy, or publication thresholds. The React client implements the
same `FlashClient` interface already consumed by `App`.

## Gates

Three independent Android gates protect the path:

- `PT Android NDK portability`: compile/audit the native runtime for
  `arm64-v8a` and `x86_64`, including the product JNI exports;
- `PT Android JNI emulator`: retain the minimal no-UI JNI/PtService
  discovery/solve regression;
- `PT Android Product Shell`: build the shared React UI, generate the Capacitor
  project, build the debug APK, boot an API 35 x86_64 emulator, render the real
  React UI, and exercise JavaScript -> Capacitor -> JNI -> `PtService` discovery
  plus one accepted PR76, SW92, and CPA solve.

Android's native accessibility hierarchy exposes a WebView as one opaque node,
so it cannot by itself prove which React DOM content rendered. For the debug
engineering APK only, `MainActivity` enables WebView DevTools when Android marks
the app debuggable. The Product Shell gate uses an ADB-forwarded local DevTools
socket and Chrome DevTools Protocol to inspect the live DOM without OCR. It
requires a populated React root, the CI smoke marker, and the rendered shared UI
text `MPMC_HNU` plus `Model-neutral PT Flash`, then removes the app at the end of
the run. This debugging surface is not enabled by this code for non-debuggable
builds.

## Deliberate exclusions

Product Shell v1 is an **engineering/debug APK**, not a public Android release.
It intentionally has no Play/App signing identity, release keystore, AAB/Play
track, updater, persistent project storage, telemetry, account system, or cloud
transport. JavaScript cancellation/deadline prevents stale UI publication, but
v1 does not yet hard-cancel an already-running C++ solve. Those are later,
separate product/release increments and must not be used to change the scientific
solver boundary.
