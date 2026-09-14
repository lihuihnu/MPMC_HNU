# Windows desktop release identity change policy

The values in `release-identity.json` are release-control inputs, not build
conveniences.

- `upgrade_guid` is the permanent Windows upgrade family identity and must not
  change for future MPMC PT Desktop upgrades.
- `release.product_guid` identifies the exact RC/version package identity. A
  later RC or product version must receive a reviewed new ProductCode.
- `release.msi_product_version`, `release.display_version`, and
  `release.sequence` move together under an explicit release increment.
- `install_directory` and `application.release_executable` are user-visible
  installation identities and require migration review before change.
- Authenticode requirements may only become stricter; removing SHA-256,
  RFC 3161 timestamps, or one of the required signed roles is not compatible
  with the current RC contract.

The protected release-candidate workflow independently pins the current RC1
version and GUIDs so an accidental identity edit fails closed before any
signing material is used.
