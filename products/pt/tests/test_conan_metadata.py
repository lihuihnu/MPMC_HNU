import json
import pathlib
import subprocess
import sys
import tempfile
import unittest


PRODUCT_ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = PRODUCT_ROOT / "scripts" / "prepare_conan_metadata.py"
PACKAGES = (
    "abseil",
    "c-ares",
    "grpc",
    "openssl",
    "protobuf",
    "re2",
    "zlib",
)


class ConanMetadataTest(unittest.TestCase):
    def _graph(self, root: pathlib.Path, binary: str = "Cache") -> pathlib.Path:
        nodes: dict[str, dict[str, object]] = {
            "0": {"ref": "conanfile", "context": "host"}
        }
        for index, name in enumerate(PACKAGES, start=1):
            package = root / name / "package"
            recipe = root / name / "recipe"
            (package / "licenses").mkdir(parents=True)
            recipe.mkdir()
            (package / "licenses" / "LICENSE.txt").write_text(
                f"license for {name}\n", encoding="utf-8"
            )
            if name == "openssl":
                executable = package / "bin" / "openssl"
                executable.parent.mkdir()
                executable.write_text("test executable\n", encoding="utf-8")
            nodes[str(index)] = {
                "ref": f"{name}/1.0#recipe-revision",
                "name": name,
                "context": "host",
                "binary": binary,
                "package_folder": str(package),
                "recipe_folder": str(recipe),
                "package_id": f"package-{index}",
                "prev": f"revision-{index}",
                "license": "test-license",
            }
        graph = root / "graph.json"
        graph.write_text(json.dumps({"graph": {"nodes": nodes}}), encoding="utf-8")
        return graph

    def test_emits_relocatable_provenance_and_all_notices(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            output = root / "metadata"
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--graph",
                    str(self._graph(root)),
                    "--lockfile",
                    str(PRODUCT_ROOT / "conan.lock"),
                    "--output",
                    str(output),
                    "--profile-id",
                    "test-profile",
                    "--conan-version",
                    "2.32.0",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertIn("RESTORE_ONLY_DEPENDENCIES_OK", result.stdout)
            manifest = json.loads(
                (output / "dependency-manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                {package["name"] for package in manifest["packages"]}, set(PACKAGES)
            )
            for name in PACKAGES:
                self.assertTrue(
                    (output / "third-party" / name / "copyright").is_file()
                )
            generated = (output / "dependency-metadata.cmake").read_text(
                encoding="utf-8"
            )
            self.assertIn('DEPENDENCY_PROVIDER "conan"', generated)
            self.assertIn("MPMC_PT_PRODUCT_OPENSSL_EXECUTABLE", generated)

    def test_rejects_a_download_or_source_build_in_staging_graph(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--graph",
                    str(self._graph(root, binary="Download")),
                    "--lockfile",
                    str(PRODUCT_ROOT / "conan.lock"),
                    "--output",
                    str(root / "metadata"),
                    "--profile-id",
                    "test-profile",
                    "--conan-version",
                    "2.32.0",
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("non-cache binary", result.stderr)


if __name__ == "__main__":
    unittest.main()
