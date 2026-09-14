# Windows Electron desktop installer gates

This directory owns the Windows x64 packaging boundary after the audited
Electron desktop payload and embedded native PT staging tree. It does not add
EOS, flash, fallback, retry, or result-interpretation logic.

## Engineering installability gate

`PT Windows desktop installer` remains the ordinary pull-request gate. It:

1. rebuilds the exact Windows x64 native staging payload from the locked Conan
   graph with `--no-remote --build=never`;
2. builds the shared React/Electron renderer and the engineering desktop
   preview;
3. wraps the complete directory in a per-machine WiX MSI under
   `Program Files\MPMC PT Desktop`;
4. performs a real silent `msiexec /i` installation;
5. launches the **installed** Electron executable in the opt-in hidden smoke
   mode;
6. feeds a versioned external smoke plan through the existing preload bridge
   and requires one scientific result from each configured PR76, SW92 and CPA
   backend;
7. verifies the Start Menu shortcut resolves to the installed executable; and
8. performs `msiexec /x` and requires both the executable and shortcut to be
   removed.

The installed-smoke plan is supplied only by CI. The production Electron bundle
contains generic plan parsing and bridge orchestration but no model-specific
backend IDs, compositions, EOS branching, or scientific acceptance thresholds.

The engineering MSI remains explicitly unsigned and not release-eligible.

## Fixed Windows release identity

`release-identity.json` is the single authoritative identity for the first
Windows desktop release candidate:

- logical product: `mpmc-pt-desktop` / `MPMC PT Desktop`;
- RC display version: `0.1.0-rc.1`;
- Windows Installer ProductVersion: `0.1.0`;
- per-machine install directory: `Program Files\MPMC PT Desktop`;
- permanent UpgradeCode: `40656E02-E6EE-59D5-A5A6-73F3C3C452D6`;
- RC1 ProductCode: `68080E81-97BF-5EDE-A618-FE9DE579A101`.

The UpgradeCode is inherited from the already-validated installable MSI and
must remain stable across future Windows upgrades. The ProductCode is frozen
for this exact RC identity; a later RC/version must receive a reviewed new
ProductCode instead of silently mutating RC1.

The normal engineering package path remains unchanged. Only a build that
explicitly passes the release identity to `frontend/desktop/packageDesktop.mjs`
produces `MPMC-PT-Desktop.exe` and the release-candidate payload manifest.

## Protected Windows release-candidate gate

`PT Windows desktop release candidate` is deliberately **manual only**
(`workflow_dispatch`) and is bound to the `windows-desktop-release` GitHub
environment. It refuses to run from anything except `main`.

The environment must provide all four protected settings:

- `WINDOWS_AUTHENTICODE_PFX_BASE64`;
- `WINDOWS_AUTHENTICODE_PFX_PASSWORD`;
- `WINDOWS_AUTHENTICODE_TIMESTAMP_URL`;
- `WINDOWS_AUTHENTICODE_EXPECTED_SUBJECT`.

There is no unsigned fallback. If any setting is missing, malformed, expired,
has the wrong subject/EKU, cannot reach the RFC 3161 timestamp service, or
fails verification, the RC gate stops before publishing an artifact.

The signing interface imports the supplied PFX only into the ephemeral runner's
CurrentUser certificate store and selects the code-signing certificate by
thumbprint. The PFX password is never passed to SignTool. The gate requires
SHA-256 file digests plus RFC 3161/SHA-256 timestamps and signs exactly these
release-owned roles:

1. `MPMC-PT-Desktop.exe`;
2. `resources/desktop-native/bin/mpmc_pt_service_host.exe`;
3. `mpmc-pt-desktop-0.1.0-rc.1-windows-x86_64.msi`.

Each signature is checked with SignTool and `Get-AuthenticodeSignature`; the
signer subject must equal the protected expected subject and a timestamp
certificate must be present. Only after payload signing does the gate build the
fixed ProductCode/UpgradeCode MSI, then it signs the MSI itself.

`windows-desktop-release-candidate-provenance.json` is emitted only when the
three signed roles, their SHA-256 values, the fixed identity, native/dependency
provenance, and the source revision all agree. At that point it records
`release_candidate_eligible=true` but still keeps `release_eligible=false`.
Final public-release approval, automatic-update policy, support policy, and
distribution policy remain later gates.

The signed RC is then subjected to the same real install -> Start Menu ->
three-backend solve -> uninstall path, with an additional check that the
**installed** desktop executable and native host still have valid timestamped
Authenticode signatures.

Until real production signing material is configured in the protected
environment, this RC workflow is expected to fail closed at its signing-input
preflight. That failure is intentional and is not replaced with a test
certificate or fabricated release evidence.
