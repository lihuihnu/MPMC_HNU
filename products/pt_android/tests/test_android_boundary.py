from __future__ import annotations

import pathlib
import unittest

ANDROID_ROOT = pathlib.Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ANDROID_ROOT.parents[1]
WORKFLOW = REPOSITORY_ROOT / ".github" / "workflows" / "pt_android_ndk_portability.yml"


class AndroidPortabilityBoundaryTest(unittest.TestCase):
    def test_android_product_stays_out_of_transport_and_ui_layers(self) -> None:
        texts = "\n".join(
            path.read_text(encoding="utf-8")
            for path in (
                ANDROID_ROOT / "CMakeLists.txt",
                ANDROID_ROOT / "src" / "portability.cpp",
                ANDROID_ROOT / "include" / "mpmc" / "pt_android" / "portability.h",
            )
        ).lower()
        for forbidden in (
            "runtime_grpc",
            "pt_process",
            "grpc",
            "electron",
            "capacitor",
            "jni.h",
            "android/log.h",
        ):
            self.assertNotIn(forbidden, texts)

    def test_portability_translation_unit_covers_current_pt_backend_types(self) -> None:
        source = (ANDROID_ROOT / "src" / "portability.cpp").read_text(
            encoding="utf-8"
        )
        self.assertIn("Pr76PtFlashBackend", source)
        self.assertIn("Sw92ProfileCPtFlashBackend", source)
        self.assertIn("CpaPtFlashBackend", source)
        self.assertIn("runtime::PtService", source)

    def test_workflow_pins_ndk_and_builds_two_android_abis(self) -> None:
        workflow = WORKFLOW.read_text(encoding="utf-8")
        self.assertNotIn("self-hosted", workflow)
        self.assertIn("ubuntu-24.04", workflow)
        self.assertIn("30.0.16248370", workflow)
        self.assertIn("arm64-v8a", workflow)
        self.assertIn("x86_64", workflow)
        self.assertIn("android.toolchain.cmake", workflow)
        self.assertIn("ANDROID_PLATFORM=android-26", workflow)
        self.assertIn("libmpmc_pt_android_core.so", workflow)

    def test_android_cmake_only_links_model_neutral_runtime_chain(self) -> None:
        cmake = (ANDROID_ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
        self.assertIn("mpmc::runtime", cmake)
        self.assertNotIn("mpmc::runtime_grpc", cmake)
        self.assertNotIn("mpmc::pt_process", cmake)
        self.assertNotIn("products/pt", cmake)


if __name__ == "__main__":
    unittest.main()
