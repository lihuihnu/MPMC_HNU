# Windows Desktop RC Gate acceptance

The protected RC gate is accepted only when all of the following are true:

1. source ref is `main`;
2. the fixed `0.1.0-rc.1` identity matches `release-identity.json`;
3. restore-only Conan dependency provenance is unchanged;
4. frontend and desktop packaging tests pass;
5. `MPMC-PT-Desktop.exe` and `mpmc_pt_service_host.exe` are signed by the
   configured production Authenticode identity with SHA-256 and RFC 3161
   timestamps;
6. the fixed WiX ProductCode/UpgradeCode MSI is built and itself receives the
   same verified timestamped Authenticode signature;
7. signed RC provenance covers all three roles and exact SHA-256 values;
8. the signed MSI installs silently, the installed payload signatures remain
   valid, the installed Electron app completes PR76/SW92/CPA smoke solves, and
   MSI uninstall removes the application and shortcut; and
9. only then is the signed RC artifact uploaded.

Missing production signing material is an expected hard failure. There is no
unsigned RC fallback.
