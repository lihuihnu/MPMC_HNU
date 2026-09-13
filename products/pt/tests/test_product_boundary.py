import pathlib
import unittest


PRODUCT_ROOT = pathlib.Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PRODUCT_ROOT.parents[1]
LOCKED_RECIPES = {
    "abseil/20260107.1": "883e95a7be2b767999f64669ac642e5d",
    "c-ares/1.34.8": "074b97a9e38159813ed5e2dcfeb290ca",
    "grpc/1.81.1": "b1374b9a466d214d491948c17c77acc3",
    "openssl/3.6.4": "bc460bf37686ccae86a15bfdce71be11",
    "protobuf/6.33.5": "ca5ff466767b31a1b496ec60247e105c",
    "re2/20251105": "8579cfd0bda4daf0683f9e3898f964b4",
    "zlib/1.3.2": "1cb806da49011867778ffb6ac7190fcb",
}


class ProductBoundaryTest(unittest.TestCase):
    def test_dependency_seed_and_staging_are_strictly_separated(self):
        conanfile = (PRODUCT_ROOT / "conanfile.txt").read_text(encoding="utf-8")
        lockfile = (PRODUCT_ROOT / "conan.lock").read_text(encoding="utf-8")
        for reference, revision in LOCKED_RECIPES.items():
            self.assertIn(reference, conanfile)
            self.assertIn(f'"{reference}#{revision}"', lockfile)
        self.assertIn("grpc/*:secure=True", conanfile)

        workflow = (
            REPOSITORY_ROOT / ".github/workflows/pt_product_staging.yml"
        ).read_text(encoding="utf-8")
        cmake = (PRODUCT_ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
        self.assertNotIn("self-hosted", workflow)
        self.assertNotIn("timeout-minutes: 90", workflow)
        self.assertGreaterEqual(workflow.count("timeout-minutes: 30"), 2)
        self.assertIn("uses: actions/cache/restore@", workflow)
        self.assertIn("uses: actions/cache/save@", workflow)
        self.assertIn("fail-on-cache-miss: true", workflow)
        self.assertIn("--build=never", workflow)
        self.assertIn("--no-remote", workflow)
        self.assertNotIn("if: always()", workflow)
        stage = workflow.split("  stage:", maxsplit=1)[1]
        self.assertNotIn("actions/cache/save@", stage)
        self.assertIn("--no-remote", stage)
        self.assertIn("RESTORE_ONLY_DEPENDENCIES_OK", workflow)
        self.assertIn("/external:W0", cmake)
        self.assertNotIn("/wd4996", cmake)
        grpc_headers = (
            REPOSITORY_ROOT
            / "modules/runtime_grpc/include/mpmc/runtime_grpc/grpc_headers.hpp"
        ).read_text(encoding="utf-8")
        self.assertIn("#pragma warning(push)", grpc_headers)
        self.assertIn("#pragma warning(disable : 4996)", grpc_headers)
        self.assertIn("#pragma warning(pop)", grpc_headers)
        for runner in ("ubuntu-24.04", "windows-2022", "macos-15"):
            self.assertIn(f"os: {runner}", workflow)
        for profile in (
            "linux-x86_64-release",
            "windows-x86_64-release",
            "macos-armv8-release",
        ):
            self.assertIn(f"profile: {profile}", workflow)
            profile_text = (
                PRODUCT_ROOT / "conan" / "profiles" / profile
            ).read_text(encoding="utf-8")
            self.assertIn("build_type=Release", profile_text)

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
            '"lock_id": "@MPMC_PT_PRODUCT_DEPENDENCY_LOCK_ID@"',
            staging_manifest,
        )


if __name__ == "__main__":
    unittest.main()
