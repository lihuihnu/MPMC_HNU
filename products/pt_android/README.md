# PT Android NDK portability v1

This directory is an isolated Android portability probe for the existing C++ PT
core. It does **not** change the Windows/Linux/macOS product, desktop process,
gRPC service, scientific models, parameters, solver thresholds, or publication
rules.

## Scope of v1

The v1 target builds `libmpmc_pt_android_core.so` with the Android NDK and links
only the existing model-neutral `mpmc::runtime` dependency chain. Its translation
unit includes the established PR76, SW92 Profile-C, and CPA PT backend types so
that Android/Clang must parse and instantiate their public class relationships.

The exported ABI is deliberately limited to a version and convention probe. It
is **not** a solve API and carries no component data, parameter snapshots, JNI,
gRPC, process-management, React, Electron, or Android UI behavior.

CI cross-compiles the same source for:

- `arm64-v8a` — primary modern Android-device ABI;
- `x86_64` — emulator/CI-oriented ABI.

The gate verifies the expected ELF machine architecture and the two exported C
symbols. A green gate means the current C++ PT core is NDK-compilable for those
ABIs. It does not yet mean an APK can execute a flash calculation.

## Next isolated increment

After this gate is stable, add a separate JNI bridge that consumes this native
core and exercises a repository-curated `PtService` discovery/solve path without
modifying the existing desktop or gRPC boundaries.
