# PT Android native runtime v2

This directory remains an isolated Android product path for the existing C++ PT
core. It does **not** change the Windows/Linux/macOS products, Electron process,
gRPC transport, EOS/flash equations, parameter values, solver thresholds, or
publication rules.

## Native boundary

`libmpmc_pt_android_core.so` is built with the Android NDK for `arm64-v8a` and
`x86_64`. It links the model-neutral `mpmc::runtime` chain and directly compiles
the existing repository-curated PT snapshot/backend assembly translation units:

- `modules/pt_process/src/configured_backends.cpp`;
- `modules/pt_process/src/parameter_snapshot_supplier.cpp`.

Android does **not** link `runtime_grpc` or the process host. A product-private
`composition_root.hpp` shim mirrors only `OwnedConfiguredPtBackend` and its
aliasing-owner helper, which are the two canonical translation units' ownership
requirements. This keeps the scientific parameter source singular while leaving
the desktop/process composition root untouched.

The original C portability probes remain exported. v2 additionally exports the
JNI entry point used by the emulator smoke application.

## Emulator acceptance path

The separate `PT Android JNI emulator` gate builds an x86_64 native library and
a deliberately minimal, offline APK without introducing Gradle/AndroidX/UI
framework dependencies. The APK contains two Java classes only:

1. `NativeBridge`, which loads `libmpmc_pt_android_core.so` and calls JNI;
2. `SmokeActivity`, which runs the native smoke off the Android UI thread and
   reports the result through logcat.

Inside the Android app process, JNI performs the real runtime sequence:

1. load the canonical repository-curated PR76, SW92 Profile-C, and CPA snapshot
   bundle;
2. construct the existing `PtService` over the three owned configured backends;
3. call `discover_capabilities()` and require all three configured backend IDs;
4. solve the established PR76 methane/ethane/propane state;
5. solve the established SW92 CO2/water fresh-water state;
6. solve the established CPA methanol/water 333.15 K state; and
7. require a structurally valid accepted service result from every solve.

The hosted-runner gate then boots an Android x86_64 emulator, installs the APK,
launches the activity, waits for:

`ANDROID_JNI_PT_SERVICE_OK backends=3`

and finally uninstalls the smoke application.

## Deliberate exclusions

This is still not an Android end-user product. There is no Capacitor/React shell,
no production Android package identity, no release signing, no persistent data,
no network transport, and no public JNI solve API. The next product increment
should wrap this proven in-process JNI path with the Android application shell
without changing the native scientific boundary.
