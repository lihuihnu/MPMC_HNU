import json
import pathlib
import unittest


PRODUCT_ROOT = pathlib.Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PRODUCT_ROOT.parents[1]
LOCKED_VCPKG_BASELINE = "a1cae005c39be7b18ba319fced856b68d7276271"


class ProductBoundaryTest(unittest.TestCase):
    def test_dependency_manifest_and_workflow_share_one_baseline(self):
        manifest = json.loads(
            (PRODUCT_ROOT / "vcpkg.json").read_text(encoding="utf-8")
        )
        self.assertEqual(manifest["builtin-baseline"], LOCKED_VCPKG_BASELINE)
        dependencies = {
            item if isinstance(item, str) else item["name"]
            for item in manifest["dependencies"]
        }
        self.assertEqual(dependencies, {"grpc", "openssl", "protobuf"})
        feature_sets = {
            item["name"]: set(item.get("features", []))
            for item in manifest["dependencies"]
            if isinstance(item, dict)
        }
        self.assertEqual(feature_sets["grpc"], {"codegen"})
        self.assertEqual(feature_sets["openssl"], {"tools"})

        workflow = (
            REPOSITORY_ROOT / ".github/workflows/pt_product_staging.yml"
        ).read_text(encoding="utf-8")
        self.assertIn(f"ref: {LOCKED_VCPKG_BASELINE}", workflow)
        self.assertNotIn("self-hosted", workflow)
        self.assertIn("uses: actions/cache/restore@", workflow)
        self.assertIn("uses: actions/cache/save@", workflow)
        self.assertIn(
            "if: always() && steps.vcpkg-cache.outputs.cache-hit != 'true'",
            workflow,
        )
        self.assertIn("-DVCPKG_OVERLAY_TRIPLETS=", workflow)
        self.assertIn("-DVCPKG_HOST_TRIPLET=${{ matrix.triplet }}", workflow)
        self.assertIn("-DVCPKG_TARGET_TRIPLET=${{ matrix.triplet }}", workflow)
        for runner in ("ubuntu-24.04", "windows-2022", "macos-15"):
            self.assertIn(f"os: {runner}", workflow)
        for triplet in (
            "x64-linux-release",
            "x64-windows-release",
            "arm64-osx-release",
        ):
            self.assertIn(f"triplet: {triplet}", workflow)
            triplet_text = (PRODUCT_ROOT / "triplets" / f"{triplet}.cmake").read_text(
                encoding="utf-8"
            )
            self.assertIn("set(VCPKG_BUILD_TYPE release)", triplet_text)

    def test_staging_does_not_absorb_web_or_model_logic(self):
        cmake = (PRODUCT_ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
        smoke = (PRODUCT_ROOT / "tests/staged_smoke.py").read_text(
            encoding="utf-8"
        )
        for forbidden in (
            "deploy/pt-grpc-web",
            "frontend/",
            "solve_pr76",
            "solve_sw92",
            "solve_cpa",
        ):
            self.assertNotIn(forbidden, cmake)
            self.assertNotIn(forbidden, smoke)

        staging_manifest = (
            PRODUCT_ROOT / "product-staging-manifest.json.in"
        ).read_text(encoding="utf-8")
        self.assertIn('"transport": "native-grpc-mtls"', staging_manifest)
        self.assertIn('"loopback_only": true', staging_manifest)
        self.assertIn(
            '"requires_explicit_tls_paths": true', staging_manifest
        )
        self.assertIn(
            '"builtin_baseline": "@MPMC_PT_VCPKG_BASELINE@"',
            staging_manifest,
        )


if __name__ == "__main__":
    unittest.main()
