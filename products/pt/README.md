# PT cross-platform product staging

This directory defines the relocatable native PT product staging gate. It is
not an end-user installer and contains no desktop shell; a downstream Electron
job consumes its exact staged tree without rebuilding the native product.

Hosted product CI consumes the checked-in Conan lock and one fixed Release
profile. All packages must already be published: both dependency seeding and
product staging pass `--build=never`, and staging additionally passes
`--no-remote`. The retained vcpkg manifest/triplets remain available as a
source-build path for local development, but they are not used by the
30-minute hosted staging gate.

The build produces a staging tree containing:

```text
bin/mpmc_pt_service_host[.exe]
lib/<non-system runtime libraries, when required>
share/mpmc-pt/product-staging-manifest.json
share/mpmc-pt/dependency-manifest.json
share/mpmc-pt/third-party/*/copyright
```

The host contains the exact repository-curated parameter bundle
`MPMC/PT/repository-curated-literature-snapshots/v1@r1`. It still exposes only
native gRPC with mandatory mutual TLS and registers all configured backends in
one `PtService`. Product staging does not add EOS, flash, retry, fallback or
result-interpretation logic.

## Restore-only CI contract

`PT cross-platform product staging` has three separate matrix jobs on the
official Ubuntu 24.04, Windows Server 2022 and macOS 15 runners:

1. `seed` restores the exact cache key. On a miss it downloads locked published
   binaries with source builds disabled, verifies the graph offline, and saves
   the cache only after that verification succeeds.
2. `stage` requires an exact cache hit, resolves the graph offline with source
   builds disabled, validates that every dependency node came from the local
   cache, then configures, builds, installs and runs the product tests. This job
   never saves or repairs a dependency cache.
3. `desktop` downloads that job's permission-preserving native staging archive,
   builds the shared React/Electron layer, loads the real renderer, discovers
   all three backends, requests one real solve from each, and assembles an
   unsigned engineering preview. It performs no Conan resolution or C++ build.

Seed and stage have 30-minute timeouts; desktop has a 20-minute timeout. Cache identity includes `conanfile.txt`,
`conan.lock` and all fixed profiles, so source-only changes reuse the verified
dependency seed while dependency changes create a new immutable cache entry.
The generated dependency manifest records every recipe/package revision, the
Conan client and profile, and the lockfile SHA-256; all seven dependency license
notices are copied into the staged tree.

## Native restore-only build

Install the pinned Conan client and choose the checked-in profile matching the
target (`linux-x86_64-release`, `windows-x86_64-release`, or
`macos-armv8-release`). For example, from a shell on Linux:

```text
python -m pip install conan==2.32.0
export CONAN_HOME="$PWD/build/conan-home"
conan install products/pt \
  --lockfile=products/pt/conan.lock \
  --profile:all=products/pt/conan/profiles/linux-x86_64-release \
  --build=never \
  --output-folder=build/conan-generators \
  --format=json > build/conan-graph.json
python products/pt/scripts/prepare_conan_metadata.py \
  --graph=build/conan-graph.json \
  --lockfile=products/pt/conan.lock \
  --output=build/conan-metadata \
  --profile-id=linux-x86_64-release \
  --conan-version=2.32.0
cmake -C build/conan-metadata/dependency-metadata.cmake \
  -S products/pt -B build/pt-product \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE="$PWD/build/conan-generators/conan_toolchain.cmake" \
  -DMPMC_PT_PRODUCT_BUILD_REVISION=<git-commit> \
  -DBUILD_TESTING=ON
cmake --build build/pt-product --config Release --parallel 2
ctest --test-dir build/pt-product -C Release --output-on-failure
```

Adding `--no-remote` to the Conan command reproduces the staging job after the
exact dependency cache has already been seeded.

## Retained local source-build path

Use the vcpkg checkout whose commit matches `builtin-baseline` in
`vcpkg.json`. This path is intentionally outside the hosted 30-minute gate:

```text
cmake -S products/pt -B build/pt-product \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake \
  -DVCPKG_OVERLAY_TRIPLETS=$PWD/products/pt/triplets \
  -DVCPKG_HOST_TRIPLET=x64-linux-release \
  -DVCPKG_TARGET_TRIPLET=x64-linux-release \
  -DMPMC_PT_PRODUCT_BUILD_REVISION=<git-commit> \
  -DBUILD_TESTING=ON
cmake --build build/pt-product --config Release --parallel 2
ctest --test-dir build/pt-product -C Release --output-on-failure
```

The staging smoke test runs `cmake --install` into
`build/pt-product/stage`, generates job-local test authorities, launches the
installed executable from that tree, discovers the exact PR76/SW92/CPA
inventory, checks that a service error remains outside the scientific result
arm, and requests graceful shutdown. Synthetic adapter tests separately retain
the accepted/indeterminate/service-error and deadline/cancel mappings.

## Boundary and release status

- The staged host must always receive explicit absolute TLS paths. Its
  `/run/secrets` defaults remain a server deployment convention, not a desktop
  filesystem assumption.
- Envoy, public DNS, CORS and production collectors belong to hosted Web
  deployment and are not staged here.
- Test CA private keys are generated only in a temporary test directory and
  are deleted after the smoke test.
- GitHub artifacts from this gate are unsigned engineering evidence. They are
  not public releases; project licensing, platform signing, notarization,
  desktop per-install identity and update policy remain later gates.
