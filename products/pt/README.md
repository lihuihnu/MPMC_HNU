# PT cross-platform product staging

This directory defines the relocatable native PT product staging gate. It is
not an end-user installer and does not add a desktop shell.

The build consumes the checked-in `vcpkg.json` baseline and the matching
release-only overlay triplet for its target platform. It produces a staging
tree containing:

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

## Native build

Use the vcpkg checkout whose commit matches `builtin-baseline` in
`vcpkg.json`:

```text
cmake -S products/pt -B build/pt-product \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake \
  -DVCPKG_OVERLAY_TRIPLETS=$PWD/products/pt/triplets \
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
