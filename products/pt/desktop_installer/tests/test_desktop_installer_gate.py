from __future__ import annotations

import importlib.util
import json
import pathlib
import sys
import tempfile
import unittest

INSTALLER_ROOT = pathlib.Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = INSTALLER_ROOT.parents[2]
SCRIPT_PATH = INSTALLER_ROOT / "scripts" / "write_candidate_manifest.py"
SPEC = importlib.util.spec_from_file_location("desktop_installer_manifest", SCRIPT_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load desktop installer manifest module")
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class WindowsDesktopInstallerGateTest(unittest.TestCase):
    revision = "1" * 40
    lock_id = "2" * 64

    def _package(self, root: pathlib.Path) -> pathlib.Path:
        package = root / "desktop"
        resources = package / "resources"
        native_share = resources / "desktop-native" / "share" / "mpmc-pt"
        resources.mkdir(parents=True)
        native_share.mkdir(parents=True)
        (package / MODULE.EXPECTED_EXECUTABLE).write_bytes(b"electron-exe")
        (resources / "desktop-preview-manifest.json").write_text(
            json.dumps(
                {
                    "convention": "MPMC/PT/desktop-preview/v1",
                    "source_revision": self.revision,
                    "electron_version": "44.3.0",
                    "platform": "win32",
                    "architecture": "x64",
                    "release_eligible": False,
                    "desktop_session": {
                        "transport": "native-grpc-loopback-bearer",
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
                    "dependencies": {
                        "provider": "conan",
                        "lock_id": self.lock_id,
                    },
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
        return package

    def test_candidate_provenance_is_explicitly_unsigned(self):
        with tempfile.TemporaryDirectory() as raw_directory:
            root = pathlib.Path(raw_directory)
            package = self._package(root)
            artifact = root / "mpmc-pt-desktop-0.1.0-windows-x86_64.msi"
            artifact.write_bytes(b"msi-candidate" * 256)
            manifest = MODULE.build_manifest(
                package_dir=package,
                artifact=artifact,
                expected_revision=self.revision,
                workflow_run_id=123,
                workflow_run_attempt=1,
            )
        self.assertEqual(
            manifest["convention"],
            "MPMC/PT/windows-desktop-installer-candidate/v1",
        )
        self.assertIs(manifest["release_eligible"], False)
        self.assertEqual(manifest["signature"]["status"], "unsigned-candidate")
        self.assertIs(manifest["signature"]["verified"], False)

    def test_revision_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as raw_directory:
            root = pathlib.Path(raw_directory)
            package = self._package(root)
            artifact = root / "mpmc-pt-desktop-0.1.0-windows-x86_64.msi"
            artifact.write_bytes(b"msi-candidate" * 256)
            with self.assertRaisesRegex(RuntimeError, "desktop package identity"):
                MODULE.build_manifest(
                    package_dir=package,
                    artifact=artifact,
                    expected_revision="3" * 40,
                    workflow_run_id=123,
                    workflow_run_attempt=1,
                )

    def test_workflow_performs_real_install_smoke_and_uninstall(self):
        workflow = (
            REPOSITORY_ROOT / ".github/workflows/pt_windows_desktop_installer.yml"
        ).read_text(encoding="utf-8")
        self.assertNotIn("self-hosted", workflow)
        self.assertIn("windows-2022", workflow)
        self.assertIn("fail-on-cache-miss: true", workflow)
        self.assertIn("--no-remote", workflow)
        self.assertIn("--build=never", workflow)
        self.assertIn("npm run desktop:package", workflow)
        self.assertIn("'/i'", workflow)
        self.assertIn("'/x'", workflow)
        self.assertIn("--mpmc-install-smoke", workflow)
        self.assertIn("MPMC_PT_DESKTOP_INSTALL_SMOKE_PLAN", workflow)
        self.assertIn("MPMC_PT_DESKTOP_INSTALL_SMOKE_RESULT", workflow)
        self.assertIn("Start Menu", workflow)
        self.assertIn("unsigned-desktop-msi-candidate", workflow)

    def test_packager_remains_model_neutral(self):
        cmake = (INSTALLER_ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
        script = SCRIPT_PATH.read_text(encoding="utf-8")
        for forbidden in (
            "solve_pr76",
            "solve_sw92",
            "solve_cpa",
            "fugacity",
            "common_tangent",
            "material_balance",
        ):
            self.assertNotIn(forbidden, cmake.lower())
            self.assertNotIn(forbidden, script.lower())
        self.assertIn('set(CPACK_GENERATOR "WIX")', cmake)
        self.assertIn("CPACK_PACKAGE_EXECUTABLES", cmake)
        self.assertIn("CPACK_CREATE_DESKTOP_LINKS", cmake)


if __name__ == "__main__":
    unittest.main()
