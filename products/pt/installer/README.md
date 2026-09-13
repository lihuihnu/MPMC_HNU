# PT signed-installer staging gate

This directory is the model-neutral packaging boundary after PT product
staging. It consumes only the installed staging tree and produces native
installer **candidates**:

| Target | Candidate | Release trust required later |
| --- | --- | --- |
| Linux x86_64 | Debian `.deb` | signed APT archive `Release` metadata and an approved public-key distribution path |
| Windows x86_64 | WiX `.msi` | Authenticode on the executable/runtime payload and MSI, with SHA-256 RFC 3161 timestamps |
| macOS arm64 | productbuild `.pkg` | Developer ID Application payload signing, Developer ID Installer package signing, notarization and stapling |

The packager copies the staging tree as opaque payload. It does not construct
`PtService`, register a backend, parse a flash result or call an EOS/flash
solver. The installed host continues to require explicit mTLS paths and remains
loopback-only; this gate does not invent a desktop identity or weaken transport
security.

## CI contract

`PT signed installer gate` runs independently on official Ubuntu 24.04,
Windows Server 2022 and macOS 15 runners. Each job:

1. requires the exact dependency cache previously seeded by product staging;
2. resolves the Conan graph offline with `--no-remote --build=never`;
3. rebuilds the small product payload from the current revision and reruns the
   installed mTLS/contract tests;
4. packages that staging tree with CPack DEB, WiX or productbuild;
5. inspects the native package payload and checks that the platform's release
   signing tools are present; and
6. writes `MPMC/PT/installer-candidate/v1` provenance with package and input
   SHA-256 values.

The Windows inventory is fail-closed. The current reviewed payload allows the
host and Microsoft Visual C++ redistributable runtime only; a Windows system
DLL or any new DLL fails the gate until its redistribution and provenance are
reviewed. Credential-like files are rejected on every platform.

Every job is capped at 30 minutes and never saves or repairs the dependency
cache. Uploaded evidence is named `unsigned-installer-candidate`, retained for
seven days, contains a prominent non-distribution notice, and records
`release_eligible=false`. A green candidate job is therefore not evidence of a
production-trusted signature.

## Production signing is intentionally blocked

No private key, certificate or signing password belongs in this repository.
Production signing must be a separate protected-environment workflow and must
fail unless all of these are approved and configured:

- the project license and release version/channel;
- a Windows code-signing identity usable by the hosted runner plus an approved
  RFC 3161 timestamp service, and final approval of the compiler-runtime
  redistribution notices;
- Apple Developer ID Application and Installer identities plus App Store
  Connect notarization credentials and team identifier; and
- an APT repository origin/suite/component, archive signing identity, public
  key publication path and rotation policy.

Linux `.deb` files are not treated as individually platform-trusted packages:
APT trust is established by signed repository metadata. Detached test
signatures would not substitute for that release chain.

## Current product coverage

This first gate covers Linux x86_64 DEB, Windows x86_64 MSI and macOS arm64 PKG.
It does not yet claim universal one-click support. Windows arm64, macOS x86_64
or universal packages, non-Debian Linux formats, desktop launcher/service
lifecycle, per-install mTLS identity enrollment and automatic updates remain
separate product gates.

## Authoritative packaging references

- [CPack WiX generator](https://cmake.org/cmake/help/latest/cpack_gen/wix.html)
- [CPack productbuild generator](https://cmake.org/cmake/help/latest/cpack_gen/productbuild.html)
- [CPack DEB generator](https://cmake.org/cmake/help/latest/cpack_gen/deb.html)
- [Microsoft SignTool](https://learn.microsoft.com/en-us/windows-hardware/drivers/devtest/signtool)
- [Apple notarization](https://developer.apple.com/documentation/security/notarizing-macos-software-before-distribution)
- [Debian package-signing model](https://www.debian.org/doc/manuals/securing-debian-manual/deb-pack-sign.en.html)
