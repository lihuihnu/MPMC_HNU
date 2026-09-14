from __future__ import annotations

import pathlib
import unittest

ANDROID_ROOT = pathlib.Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ANDROID_ROOT.parents[1]
PORTABILITY_WORKFLOW = (
    REPOSITORY_ROOT / ".github" / "workflows" / "pt_android_ndk_portability.yml"
)
EMULATOR_WORKFLOW = (
    REPOSITORY_ROOT / ".github" / "workflows" / "pt_android_jni_emulator.yml"
)


class AndroidPortabilityBoundaryTest(unittest.TestCase):
    def test_android_product_stays_out_of_existing_transport_and_ui_layers(self) -> None:
        texts = "\n".join(
            path.read_text(encoding="utf-8")
            for path in (
                ANDROID_ROOT / "CMakeLists.txt",
                ANDROID_ROOT / "src" / "portability.cpp",
                ANDROID_ROOT / "src" / "jni_smoke.cpp",
                ANDROID_ROOT
                / "private_include"
                / "mpmc"
                / "pt_process"
                / "composition_root.hpp",
            )
        ).lower()
        for forbidden_token in (
            "mpmc::runtime_grpc",
            "<mpmc/runtime_grpc/",
            "<grpc/",
            "<electron/",
            "<capacitor/",
            "#include <android/log.h>",
        ):
            self.assertNotIn(forbidden_token, texts)

    def test_canonical_repository_curated_sources_are_reused_without_process_host(self) -> None:
        cmake = (ANDROID_ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
        self.assertIn("modules/pt_process/src/configured_backends.cpp", cmake)
        self.assertIn("modules/pt_process/src/parameter_snapshot_supplier.cpp", cmake)
        self.assertIn('modules/runtime"', cmake)
        self.assertIn("mpmc::runtime", cmake)
        self.assertIn("private_include", cmake)
        self.assertNotIn("add_subdirectory(\"${mpmc_repository_root}/modules/pt_process", cmake)
        self.assertNotIn("mpmc::runtime_grpc", cmake)

        shim = (
            ANDROID_ROOT
            / "private_include"
            / "mpmc"
            / "pt_process"
            / "composition_root.hpp"
        ).read_text(encoding="utf-8")
        self.assertIn("OwnedConfiguredPtBackend", shim)
        self.assertIn("retain_configured_pt_backend", shim)
        self.assertNotIn("class PtCompositionRoot", shim)
        self.assertNotIn("<mpmc/runtime_grpc/", shim)

    def test_jni_smoke_uses_real_pt_service_discovery_and_three_backends(self) -> None:
        source = (ANDROID_ROOT / "src" / "jni_smoke.cpp").read_text(
            encoding="utf-8"
        )
        self.assertIn("#include <jni.h>", source)
        self.assertIn("load_repository_curated_pt_parameter_snapshots_v1", source)
        self.assertIn("rt::PtService service", source)
        self.assertIn("discover_capabilities", source)
        self.assertIn("pr76.methane-ethane-propane.literature-r1", source)
        self.assertIn("sw92.carbon-dioxide-water.freshwater.literature-r1", source)
        self.assertIn("cpa.methanol-water-333.15k.cr1.literature-r1", source)
        self.assertIn("ANDROID_JNI_PT_SERVICE_OK", source)
        self.assertIn("PtServiceOutcome::accepted", source)

    def test_emulator_app_is_minimal_offline_and_loads_native_library(self) -> None:
        manifest = (ANDROID_ROOT / "emulator" / "AndroidManifest.xml").read_text(
            encoding="utf-8"
        )
        bridge = (
            ANDROID_ROOT
            / "emulator"
            / "java"
            / "org"
            / "mpmc"
            / "ptandroid"
            / "NativeBridge.java"
        ).read_text(encoding="utf-8")
        activity = (
            ANDROID_ROOT
            / "emulator"
            / "java"
            / "org"
            / "mpmc"
            / "ptandroid"
            / "SmokeActivity.java"
        ).read_text(encoding="utf-8")
        self.assertNotIn("android.permission.INTERNET", manifest)
        self.assertIn('android:minSdkVersion="26"', manifest)
        self.assertIn('android:targetSdkVersion="35"', manifest)
        self.assertIn('System.loadLibrary("mpmc_pt_android_core")', bridge)
        self.assertIn("native String runSmoke", bridge)
        self.assertIn("new Thread", activity)
        self.assertIn("ANDROID_JNI_PT_SERVICE_FAIL", activity)

    def test_portability_workflow_builds_real_jni_runtime_on_two_abis(self) -> None:
        workflow = PORTABILITY_WORKFLOW.read_text(encoding="utf-8")
        self.assertNotIn("self-hosted", workflow)
        self.assertIn("ubuntu-24.04", workflow)
        self.assertIn("30.0.16248370", workflow)
        self.assertIn("arm64-v8a", workflow)
        self.assertIn("x86_64", workflow)
        self.assertIn("android.toolchain.cmake", workflow)
        self.assertIn("ANDROID_PLATFORM=android-26", workflow)
        self.assertIn("ANDROID_STL=c++_static", workflow)
        self.assertIn("Java_org_mpmc_ptandroid_NativeBridge_runSmoke", workflow)

    def test_emulator_workflow_builds_installs_launches_and_checks_result(self) -> None:
        workflow = EMULATOR_WORKFLOW.read_text(encoding="utf-8")
        build_script = (
            ANDROID_ROOT / "emulator" / "build_smoke_apk.sh"
        ).read_text(encoding="utf-8")
        run_script = (
            ANDROID_ROOT / "emulator" / "run_emulator_smoke.sh"
        ).read_text(encoding="utf-8")
        self.assertNotIn("self-hosted", workflow)
        self.assertIn("ubuntu-24.04", workflow)
        self.assertIn("30.0.16248370", workflow)
        self.assertIn("system-images;android-", workflow)
        self.assertIn("google_apis;x86_64", workflow)
        self.assertIn("ANDROID_STL=c++_static", workflow)
        for token in ("aapt2", "d8", "zipalign", "apksigner"):
            self.assertIn(token, build_script)
        self.assertIn('install -r "$apk"', run_script)
        self.assertIn("am start -W", run_script)
        self.assertIn("ANDROID_JNI_PT_SERVICE_OK backends=3", run_script)
        self.assertIn("ANDROID_EMULATOR_JNI_PT_SERVICE_OK", run_script)


if __name__ == "__main__":
    unittest.main()
