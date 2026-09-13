from __future__ import annotations

import importlib.util
import json
import pathlib
import sys
import tempfile
import unittest


INSTALLER_ROOT = pathlib.Path(__file__).resolve().parents[1]
PRODUCT_ROOT = INSTALLER_ROOT.parent
REPOSITORY_ROOT = PRODUCT_ROOT.parents[1]
SCRIPT_PATH = INSTALLER_ROOT / "scripts" / "write_candidate_manifest.py"
SPEC = importlib.util.spec_from_file_location("write_candidate_manifest", SCRIPT_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load installer manifest module")
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class InstallerGateTest(unittest.TestCase):
    revision = "1" * 40
    lock_id = "2" * 64

    def _stage(self, root: pathlib.Path, target_id: str) -> pathlib.Path:
        target = MODULE.TARGETS[target_id]
        stage = root / "stage"
        bin_directory = stage / "bin"
        share_directory = stage / "share" / "mpmc-pt"
        bin_directory.mkdir(parents=True)
        share_directory.mkdir(parents=True)
        entry_name = pathlib.PurePosixPath(target.entry).name
        (bin_directory / entry_name).write_bytes(b"host")
        staging = {
            "convention": "MPMC/PT/product-staging/v1",
            "product": {
                "id": "mpmc-pt",
                "version": "0.1.0",
                "build_revision": self.revision,
            },
            "platform": {
                "system_name": target.system_name,
                "system_processor": target.processor,
                "dependency_profile": target.profile,
            },
            "dependencies": {
                "provider": "conan",
                "provider_version": "2.32.0",
                "lock_kind": "conan-lock-sha256",
                "lock_id": self.lock_id,
            },
            "entry": {
                "path": f"bin/{entry_name}",
                "transport": "native-grpc-mtls",
                "loopback_only": True,
                "requires_explicit_tls_paths": True,
            },
            "snapshot_bundle": {
                "id": "MPMC/PT/repository-curated-literature-snapshots/v1",
                "revision": "r1",
                "location": "builtin://mpmc/test",
            },
        }
        dependency = {
            "convention": "MPMC/PT/conan-dependency-manifest/v1",
            "provider": "conan",
            "client_version": "2.32.0",
            "profile": target.profile,
            "lock_sha256": self.lock_id,
            "packages": [
                {
                    "name": name,
                    "reference": f"{name}/1.0#recipe",
                    "package_id": "package-id",
                    "package_revision": "package-revision",
                    "license": "test-only",
                }
                for name in sorted(MODULE.EXPECTED_DEPENDENCIES)
            ],
        }
        (share_directory / "product-staging-manifest.json").write_text(
            json.dumps(staging), encoding="utf-8"
        )
        (share_directory / "dependency-manifest.json").write_text(
            json.dumps(dependency), encoding="utf-8"
        )
        return stage

    def test_candidate_manifest_is_explicitly_not_release_eligible(self):
        with tempfile.TemporaryDirectory() as raw_directory:
            root = pathlib.Path(raw_directory)
            stage = self._stage(root, "linux-x86_64")
            artifact = root / "mpmc-pt-0.1.0-linux-x86_64.deb"
            artifact.write_bytes(b"candidate" * 256)
            manifest = MODULE.build_manifest(
                stage=stage,
                artifact=artifact,
                target_id="linux-x86_64",
                expected_revision=self.revision,
                workflow_run_id=123,
                workflow_run_attempt=2,
            )
        self.assertEqual(manifest["convention"], "MPMC/PT/installer-candidate/v1")
        self.assertIs(manifest["release_eligible"], False)
        self.assertEqual(manifest["signature"]["status"], "unsigned-candidate")
        self.assertIs(manifest["signature"]["verified"], False)
        self.assertEqual(
            manifest["signature"]["required_for_release"],
            "apt-archive-release-openpgp",
        )
        self.assertEqual(manifest["payload"]["dependency_package_count"], 7)

    def test_revision_and_platform_mismatch_fail_closed(self):
        with tempfile.TemporaryDirectory() as raw_directory:
            stage = self._stage(pathlib.Path(raw_directory), "macos-arm64")
            with self.assertRaisesRegex(RuntimeError, "source revision"):
                MODULE.validate_stage(stage, "macos-arm64", "3" * 40)
            with self.assertRaisesRegex(RuntimeError, "installer target"):
                MODULE.validate_stage(stage, "linux-x86_64", self.revision)

    def test_windows_system_or_unreviewed_dll_is_rejected(self):
        with tempfile.TemporaryDirectory() as raw_directory:
            stage = self._stage(pathlib.Path(raw_directory), "windows-x86_64")
            (stage / "bin" / "kernel32.dll").write_bytes(b"system")
            with self.assertRaisesRegex(RuntimeError, "system binaries"):
                MODULE.validate_stage(stage, "windows-x86_64", self.revision)

    def test_workflow_is_separate_bounded_and_restore_only(self):
        workflow = (
            REPOSITORY_ROOT / ".github/workflows/pt_signed_installer_gate.yml"
        ).read_text(encoding="utf-8")
        self.assertNotIn("pull_request_target", workflow)
        self.assertIn("contents: read", workflow)
        self.assertNotIn("self-hosted", workflow)
        self.assertEqual(workflow.count("timeout-minutes: 30"), 1)
        self.assertIn("fail-on-cache-miss: true", workflow)
        self.assertIn("--no-remote", workflow)
        self.assertIn("--build=never", workflow)
        self.assertNotIn("actions/cache/save@", workflow)
        self.assertIn("unsigned-installer-candidate", workflow)
        self.assertIn("pkgutil --expand-full", workflow)
        for runner in ("ubuntu-24.04", "windows-2022", "macos-15"):
            self.assertIn(f"os: {runner}", workflow)

    def test_packager_has_no_model_or_transport_implementation(self):
        implementation = "\n".join(
            (
                (INSTALLER_ROOT / "CMakeLists.txt").read_text(encoding="utf-8"),
                SCRIPT_PATH.read_text(encoding="utf-8"),
                (
                    REPOSITORY_ROOT
                    / ".github/workflows/pt_signed_installer_gate.yml"
                ).read_text(encoding="utf-8"),
            )
        )
        for forbidden in (
            "solve_pr76",
            "solve_sw92",
            "solve_cpa",
            "PtService(",
            "fugacity",
            "material_balance",
            "common_tangent",
        ):
            self.assertNotIn(forbidden, implementation)
        for generator in ('"DEB"', '"WIX"', '"productbuild"'):
            self.assertIn(generator, implementation)
        self.assertEqual(
            implementation.count("CPACK_RESOURCE_FILE_LICENSE"), 1
        )


if __name__ == "__main__":
    unittest.main()
