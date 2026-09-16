# PT Android Product Shell v1

This directory is the isolated Android product path for the existing C++ PT
runtime and the shared Classic PR React product workflow. It does **not** change
the Windows/Linux/macOS products, Electron process, gRPC transport, EOS/flash
equations, parameter values, solver thresholds, or publication rules.

## Native boundary

`libmpmc_pt_android_core.so` is built with the Android NDK for `arm64-v8a` and
`x86_64`. It retains the model-neutral `mpmc::runtime` path and directly compiles
the existing repository-curated PT snapshot/backend assembly translation units:

- `modules/pt_process/src/configured_backends.cpp`;
- `modules/pt_process/src/parameter_snapshot_supplier.cpp`.

For editable Classic PR, the same library also links the existing
`mpmc::pr76_solver_configuration` interface and uses
`model_configuration::Pr76ExecutableModel`. Android therefore reuses the same
parameter preparation, solver-setting validation, stability/split/max3 flash and
result contracts as the other local model-workbench path. No PR EOS or flash
formula is reimplemented in Java, TypeScript, or the Android bridge.

Android does **not** link `runtime_grpc`, Protobuf C++ transport, or the process
host. A product-private `composition_root.hpp` shim remains limited to the
configured-backend ownership primitive needed by the model-neutral path.

The native library exposes two transport families:

1. repository-curated smoke plus model-neutral discovery/solve for the configured
   PR76/SW92/CPA compatibility path; and
2. typed Classic PR ownership operations `apply`, `solve`, `release`, and
   `cancel`, backed by one active immutable `Pr76ExecutableModel`.

The configured product bridge remains `MPMC/PT/android-product-bridge/v1`.
Classic PR uses the same renderer `MPMC/model/workbench-bridge/v1` semantics as
the shared `ModelWorkbenchOwner`: Java/Capacitor transports bounded records only;
native configuration/solve errors are mapped to the existing structured model
validation categories; raw native handles are never exposed to React.

## Product Shell v1

`products/pt_android/shell/` wraps the native library with the pinned Capacitor
shell and now renders the same Classic PR product shell used by Electron:

```text
shared DesktopProductShell / ExpertPr76Workspace
        ↓                         ↓
Android ModelWorkbenchBridge      AndroidFlashClient
        ↓                         ↓
Capacitor app-local MpmcPt plugin
        ↓
JNI apply/solve/release/cancel + discover/solve
        ↓
libmpmc_pt_android_core.so
        ↓
Pr76ExecutableModel + PtService
```

The default Android UI is therefore **Classic PR**: the same component editor,
explicit `kij`, solver-settings disclosure, P/T/z form, accepted phase-count /
phase-fraction / phase-composition charts, exact table and non-accepted candidate
isolation are reused directly from `frontend/`. `Configured PT` remains the same
secondary compatibility tab and still exercises repository-curated PR76, SW92
and CPA backends.

The shell pins Capacitor `8.5.2` and React `19.2.3`. The generated Android Studio
project is deliberately **not** committed: CI generates it from the pinned
Capacitor input, keeps API 26 as the native compatibility floor, overlays the
small Java plugin, embeds validated `arm64-v8a` and `x86_64` native cores, and
builds one universal debug engineering APK. The package gate fails closed unless
both native payloads are present at the standard Android paths.

The app-local Java plugin is transport only. Apply/solve/release run on one
native executor thread. `cancel` marks native ownership immediately rather than
waiting behind that executor, so an abandoned apply cannot later publish a model.
The plugin contains no EOS equations, parameter constants, phase-selection logic,
retry policy, or publication thresholds.

## Gates

Three independent Android gates protect the path:

- `PT Android NDK portability`: compile/audit the native runtime for
  `arm64-v8a` and `x86_64`, including both configured-PT and Classic PR JNI
  exports;
- `PT Android JNI emulator`: retain the minimal no-UI JNI/PtService
  discovery/solve regression;
- `PT Android Product Shell`: build both native ABIs, require both native cores
  inside one universal APK, typecheck/build the shared React shell, boot an API
  35 x86_64 emulator, retain the configured PR76/SW92/CPA solves, and additionally
  execute a real Classic PR `apply -> solve -> release` through
  React -> Capacitor -> JNI -> `Pr76ExecutableModel`.

The Classic PR smoke copies the repository-curated methane/ethane/propane PR76
fixture already traced in `parameter_snapshot_supplier.cpp`; it does not invent a
new physical reference dataset. The model is solved at the existing accepted
Android regression state `1 MPa / 350 K / z=(0.8,0.1,0.1)`.

Android's native accessibility hierarchy exposes a WebView as one opaque node.
For the debuggable engineering APK only, `MainActivity` enables WebView DevTools.
The Product Shell gate uses an ADB-forwarded local Chrome DevTools Protocol
endpoint to verify the live DOM without OCR. It requires the populated React
root, both smoke markers, and the shared `Flash workspace`, `Classic PR`, and
`Define PR fluid` text before uninstalling the app.

## Deliberate exclusions

Product Shell v1 is an **engineering/debug APK**, not a public Android release.
It intentionally has no Play/App signing identity, release keystore, AAB/Play
track, updater, persistent project storage, telemetry, account system, or cloud
transport. Native cancellation protects model ownership/publication; it does not
interrupt the mathematical kernel mid-call. Hard interruption of C++ numerical
work remains a separate concern and must not be used to alter scientific solver
semantics.
