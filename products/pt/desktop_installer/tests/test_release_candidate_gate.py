from __future__ import annotations

import hashlib
import importlib.util
import json
import pathlib
import sys
import tempfile
import unittest

INSTALLER_ROOT = pathlib.Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = INSTALLER_ROOT.parents[2]
IDENTITY_PATH = INSTALLER_ROOT / "release-identity.json"
RC_SCRIPT_PATH = INSTALLER_ROOT / "scripts" / "write_release_candidate_manifest.py"
SPEC = importlib.util.spec_from_file_location("desktop_rc_manifest", RC_SCRIPT_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load desktop RC manifest module")
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class WindowsDesktopReleaseCandidateGateTest(unittest.TestCase):
    revision = "1" * 40
    lock_id = "2" * 64
    signer_subject = "CN=MPMC_HNU Test Signing Identity"

    def _fixture(self, root: pathlib.Path) -> tuple[pathlib.Path, pathlib.Path, pathlib.Path]:
        identity = json.loads(IDENTITY_PATH.read_text(encoding="utf-8"))
        release = identity["release"]
        application = identity["application"]
        package = root / "desktop"
        resources = package / "resources"
        native_share = resources / "desktop-native" / "share" / "mpmc-pt"
        native_bin = resources / "desktop-native" / "bin"
        resources.mkdir(parents=True)
        native_share.mkdir(parents=True)
        native_bin.mkdir(parents=True)

        desktop_executable = package / application["release_executable"]
        native_host = native_bin / "mpmc_pt_service_host.exe"
        desktop_executable.write_bytes(b"signed-electron" * 256)
        native_host.write_bytes(b"signed-native-host" * 256)

        (resources / "desktop-preview-manifest.json").write_text(
            json.dumps(
                {
                    "convention": "MPMC/PT/desktop-release-candidate-payload/v1",
                    "source_revision": self.revision,
                    "platform": "win32",
                    "architecture": "x64",
                    "release_candidate": True,
                    "release_eligible": False,
                    "desktop_session": {
                        "transport": "native-grpc-loopback-bearer",
                    },
                    "release_identity": {
                        "display_version": release["display_version"],
                        "msi_product_version": release["msi_product_version"],
                        "product_guid": release["product_guid"],
                        "upgrade_guid": identity["upgrade_guid"],
                        "identity_sha256": sha256(IDENTITY_PATH),
                    },
                    "signature": {
                        "status": "authenticode-signed-timestamped",
                        "verified": True,
                        "signer_subject": self.signer_subject,
                        "file_digest": "SHA256",
                        "timestamp_protocol": "RFC3161",
                        "timestamp_digest": "SHA256",
                    },
                }
            ),
            encoding="utf-8",
        )
        (native_share / "product-staging-manifest.json").write_text(
            json.dumps(
                {
                    "convention": "MPMC/PT/product-staging/v1",
                    "product": {"build_revision": self.revision},
                    "platform": {"dependency_profile": "windows-x86_64-release"},
                    "dependencies": {"provider": "conan", "lock_id": self.lock_id},
                }
            ),
            encoding="utf-8",
        )
        (native_share / "dependency-manifest.json").write_text(
            json.dumps(
                {
                    "convention": "MPMC/PT/conan-dependency-manifest/v1",
                    "profile": "windows-x86_64-release",
                    "lock_sha256": self.lock_id,
                }
            ),
            encoding="utf-8",
        )

        artifact = root / f"mpmc-pt-desktop-{release['display_version']}-windows-x86_64.msi"
        artifact.write_bytes(b"signed-msi" * 256)
        evidence = root / "authenticode-evidence.json"
        evidence.write_text(
            json.dumps(
                {
                    "convention": "MPMC/PT/windows-authenticode-evidence/v1",
                    "signer_subject": self.signer_subject,
                    "signer_thumbprint": "a" * 40,
                    "file_digest": "SHA256",
                    "timestamp_protocol": "RFC3161",
                    "timestamp_digest": "SHA256",
                    "files": [
                        {
                            "role": "desktop-executable",
                            "status": "Valid",
                            "timestamped": True,
                            "signer_subject": self.signer_subject,
                            "sha256": sha256(desktop_executable),
                        },
                        {
                            "role": "native-host",
                            "status": "Valid",
                            "timestamped": True,
                            "signer_subject": self.signer_subject,
                            "sha256": sha256(native_host),
                        },
                        {
                            "role": "msi",
                            "status": "Valid",
                            "timestamped": True,
                            "signer_subject": self.signer_subject,
                            "sha256": sha256(artifact),
                        },
                    ],
                }
            ),
            encoding="utf-8",
        )
        return package, artifact, evidence

    def test_release_identity_is_frozen_and_upgrade_compatible(self):
        identity = json.loads(IDENTITY_PATH.read_text(encoding="utf-8"))
        self.assertEqual(identity["convention"], "MPMC/PT/windows-desktop-release-identity/v1")
        self.assertEqual(identity["product_id"], "mpmc-pt-desktop")
        self.assertEqual(identity["release"]["display_version"], "0.1.0-rc.1")
        self.assertEqual(identity["release"]["msi_product_version"], "0.1.0")
        self.assertEqual(
            identity["upgrade_guid"],
            "40656E02-E6EE-59D5-A5A6-73F3C3C452D6",
        )
        self.assertEqual(
            identity["release"]["product_guid"],
            "68080E81-97BF-5EDE-A618-FE9DE579A101",
        )
        self.assertEqual(identity["install_directory"], "MPMC PT Desktop")
        self.assertEqual(identity["authenticode"]["file_digest"], "SHA256")
        self.assertEqual(identity["authenticode"]["timestamp_protocol"], "RFC3161")

        cmake = (INSTALLER_ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
        self.assertIn('set(CPACK_WIX_UPGRADE_GUID "${mpmc_upgrade_guid}")', cmake)
        self.assertIn('set(CPACK_WIX_PRODUCT_GUID "${mpmc_rc_product_guid}")', cmake)
        self.assertIn("MPMC_PT_DESKTOP_RELEASE_CANDIDATE", cmake)
        self.assertIn("release-candidate-notice.txt", cmake)

    def test_release_payload_mode_does_not_replace_engineering_preview_mode(self):
        packager = (
            REPOSITORY_ROOT / "frontend" / "desktop" / "packageDesktop.mjs"
        ).read_text(encoding="utf-8")
        self.assertIn("--release-identity", packager)
        self.assertIn("MPMC/PT/desktop-preview/v1", packager)
        self.assertIn("MPMC/PT/desktop-release-candidate-payload/v1", packager)
        self.assertIn("MPMC-PT-Desktop-Preview", packager)
        self.assertIn("MPMC-PT-Desktop", packager)
        self.assertIn("desktop-release-candidate-notice.txt", packager)

    def test_authenticode_interface_is_fail_closed_and_timestamped(self):
        signing = (
            INSTALLER_ROOT / "scripts" / "sign_release_candidate.ps1"
        ).read_text(encoding="utf-8")
        for required in (
            "MPMC_WINDOWS_CODESIGN_PFX_PATH",
            "MPMC_WINDOWS_CODESIGN_PFX_PASSWORD",
            "MPMC_WINDOWS_CODESIGN_TIMESTAMP_URL",
            "MPMC_WINDOWS_CODESIGN_EXPECTED_SUBJECT",
            "Import-PfxCertificate",
            "1.3.6.1.5.5.7.3.3",
            "'/fd', 'SHA256'",
            "'/tr', $timestampUrl",
            "'/td', 'SHA256'",
            "verify /pa /all /v",
            "TimeStamperCertificate",
            "authenticode-signed-timestamped",
        ):
            self.assertIn(required, signing)
        self.assertNotIn("'/p'", signing)
        self.assertNotIn("-Exportable:$true", signing)

    def test_provenance_requires_all_three_valid_timestamped_roles(self):
        with tempfile.TemporaryDirectory() as raw_directory:
            root = pathlib.Path(raw_directory)
            package, artifact, evidence = self._fixture(root)
            manifest = MODULE.build_manifest(
                identity_path=IDENTITY_PATH,
                package_dir=package,
                artifact=artifact,
                signing_evidence_path=evidence,
                expected_revision=self.revision,
                workflow_run_id=123,
                workflow_run_attempt=1,
            )
            self.assertEqual(
                manifest["convention"],
                "MPMC/PT/windows-desktop-release-candidate/v1",
            )
            self.assertIs(manifest["release_candidate_eligible"], True)
            self.assertIs(manifest["release_eligible"], False)
            self.assertEqual(
                manifest["authenticode"]["covered_roles"],
                ["desktop-executable", "msi", "native-host"],
            )

            evidence_data = json.loads(evidence.read_text(encoding="utf-8"))
            evidence_data["files"] = evidence_data["files"][:-1]
            evidence.write_text(json.dumps(evidence_data), encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "full RC surface"):
                MODULE.build_manifest(
                    identity_path=IDENTITY_PATH,
                    package_dir=package,
                    artifact=artifact,
                    signing_evidence_path=evidence,
                    expected_revision=self.revision,
                    workflow_run_id=123,
                    workflow_run_attempt=1,
                )

    def test_release_workflow_is_manual_protected_and_fail_closed(self):
        workflow = (
            REPOSITORY_ROOT
            / ".github"
            / "workflows"
            / "pt_windows_desktop_release_candidate.yml"
        ).read_text(encoding="utf-8")
        self.assertIn("workflow_dispatch:", workflow)
        self.assertNotIn("pull_request:", workflow)
        self.assertIn("environment: windows-desktop-release", workflow)
        self.assertIn("refs/heads/main", workflow)
        for secret in (
            "WINDOWS_AUTHENTICODE_PFX_BASE64",
            "WINDOWS_AUTHENTICODE_PFX_PASSWORD",
            "WINDOWS_AUTHENTICODE_TIMESTAMP_URL",
            "WINDOWS_AUTHENTICODE_EXPECTED_SUBJECT",
        ):
            self.assertIn(f"secrets.{secret}", workflow)
        self.assertIn("Protected Windows RC signing inputs are incomplete", workflow)
        payload_sign = workflow.index("-Mode payload")
        package_build = workflow.index("Build fixed Windows RC MSI")
        msi_sign = workflow.index("-Mode msi")
        provenance = workflow.index("write_release_candidate_manifest.py")
        upload = workflow.index("Upload signed Windows desktop release candidate")
        self.assertLess(payload_sign, package_build)
        self.assertLess(package_build, msi_sign)
        self.assertLess(msi_sign, provenance)
        self.assertLess(provenance, upload)
        self.assertIn("Get-AuthenticodeSignature", workflow)
        self.assertIn("TimeStamperCertificate", workflow)
        self.assertIn("--mpmc-install-smoke", workflow)
        self.assertNotIn("if (Test-Path $installedExe -or Test-Path $installedHost)", workflow)
        self.assertIn("(Test-Path $installedExe) -or (Test-Path $installedHost)", workflow)


if __name__ == "__main__":
    unittest.main()
