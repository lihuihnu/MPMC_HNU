#!/usr/bin/env python3
"""Validate a restore-only Conan graph and stage dependency provenance."""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import re
import shutil


EXPECTED_PACKAGES = {
    "abseil",
    "c-ares",
    "grpc",
    "openssl",
    "protobuf",
    "re2",
    "zlib",
}
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._+-]*$")
LICENSE_PREFIXES = ("license", "licence", "copying", "copyright", "notice")


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", required=True, type=pathlib.Path)
    parser.add_argument("--lockfile", required=True, type=pathlib.Path)
    parser.add_argument("--output", required=True, type=pathlib.Path)
    parser.add_argument("--profile-id", required=True)
    parser.add_argument("--conan-version", required=True)
    return parser.parse_args()


def _license_files(node: dict[str, object]) -> list[pathlib.Path]:
    package_folder = pathlib.Path(str(node["package_folder"]))
    recipe_folder = pathlib.Path(str(node["recipe_folder"]))
    candidates: list[pathlib.Path] = []
    license_directory = package_folder / "licenses"
    if license_directory.is_dir():
        candidates.extend(path for path in license_directory.rglob("*") if path.is_file())
    for folder in (package_folder, recipe_folder):
        if not folder.is_dir():
            continue
        candidates.extend(
            path
            for path in folder.iterdir()
            if path.is_file()
            and path.name.casefold().startswith(LICENSE_PREFIXES)
        )
    return sorted(set(candidates), key=lambda path: str(path).casefold())


def _cmake_bracket(value: pathlib.Path) -> str:
    serialized = str(value.resolve())
    if "]==]" in serialized:
        raise RuntimeError("dependency path cannot be represented safely in CMake")
    return f"[==[{serialized}]==]"


def main() -> int:
    args = _arguments()
    if not SAFE_ID.fullmatch(args.profile_id):
        raise RuntimeError("profile id must be a stable non-secret identifier")
    if not re.fullmatch(r"[0-9]+(?:[.][0-9]+){2}", args.conan_version):
        raise RuntimeError("Conan client version must be an exact semantic version")

    graph = json.loads(args.graph.read_text(encoding="utf-8"))
    nodes = graph.get("graph", {}).get("nodes", {})
    if not isinstance(nodes, dict):
        raise RuntimeError("Conan graph does not contain a node map")

    packages: dict[str, dict[str, object]] = {}
    for node in nodes.values():
        if (
            not isinstance(node, dict)
            or node.get("context") != "host"
            or not node.get("package_folder")
            or node.get("ref") == "conanfile"
        ):
            continue
        if node.get("binary") != "Cache":
            raise RuntimeError(
                f"restore-only graph contains non-cache binary {node.get('ref')}: "
                f"{node.get('binary')}"
            )
        name = str(node.get("name", ""))
        if name not in EXPECTED_PACKAGES:
            raise RuntimeError(f"unexpected host dependency in locked graph: {name}")
        if name in packages:
            raise RuntimeError(f"duplicate host dependency in locked graph: {name}")
        packages[name] = node

    if set(packages) != EXPECTED_PACKAGES:
        missing = sorted(EXPECTED_PACKAGES - set(packages))
        raise RuntimeError(f"locked host dependency set is incomplete: {missing}")

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    notices_root = output / "third-party"
    notices_root.mkdir()

    manifest_packages: list[dict[str, str]] = []
    openssl_executable: pathlib.Path | None = None
    for name in sorted(packages):
        node = packages[name]
        package_folder = pathlib.Path(str(node["package_folder"]))
        notices = _license_files(node)
        if not notices:
            raise RuntimeError(f"locked dependency {name} has no license notice")
        destination = notices_root / name
        destination.mkdir()
        combined = bytearray()
        for index, source in enumerate(notices, start=1):
            target_name = f"{index:02d}-{source.name}"
            shutil.copyfile(source, destination / target_name)
            combined.extend(f"===== {target_name} =====\n".encode("utf-8"))
            combined.extend(source.read_bytes())
            if not combined.endswith(b"\n"):
                combined.extend(b"\n")
        (destination / "copyright").write_bytes(combined)

        manifest_packages.append(
            {
                "name": name,
                "reference": str(node["ref"]),
                "package_id": str(node["package_id"]),
                "package_revision": str(node["prev"]),
                "license": str(node.get("license") or "unspecified"),
            }
        )
        if name == "openssl":
            executable_name = "openssl.exe" if package_folder.drive else "openssl"
            candidate = package_folder / "bin" / executable_name
            if candidate.is_file():
                openssl_executable = candidate

    if openssl_executable is None:
        discovered = shutil.which("openssl")
        if discovered is None:
            raise RuntimeError("no OpenSSL executable is available for the mTLS smoke")
        openssl_executable = pathlib.Path(discovered)

    lock_digest = hashlib.sha256(args.lockfile.read_bytes()).hexdigest()
    dependency_manifest = {
        "convention": "MPMC/PT/conan-dependency-manifest/v1",
        "provider": "conan",
        "client_version": args.conan_version,
        "profile": args.profile_id,
        "lock_sha256": lock_digest,
        "packages": manifest_packages,
    }
    (output / "dependency-manifest.json").write_text(
        json.dumps(dependency_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (output / "dependency-metadata.cmake").write_text(
        "\n".join(
            (
                "# Generated from the restore-only Conan graph.",
                'set(MPMC_PT_PRODUCT_DEPENDENCY_PROVIDER "conan" CACHE STRING "")',
                "set(MPMC_PT_PRODUCT_DEPENDENCY_METADATA_DIR "
                f"{_cmake_bracket(output)} CACHE PATH \"\")",
                "set(MPMC_PT_PRODUCT_OPENSSL_EXECUTABLE "
                f"{_cmake_bracket(openssl_executable)} CACHE FILEPATH \"\")",
                f'set(MPMC_PT_PRODUCT_CONAN_PROFILE_ID "{args.profile_id}" CACHE STRING "")',
                f'set(MPMC_PT_PRODUCT_CONAN_LOCK_SHA256 "{lock_digest}" CACHE STRING "")',
                f'set(MPMC_PT_PRODUCT_CONAN_VERSION "{args.conan_version}" CACHE STRING "")',
                "",
            )
        ),
        encoding="utf-8",
    )
    print(
        "RESTORE_ONLY_DEPENDENCIES_OK "
        f"profile={args.profile_id} packages={len(manifest_packages)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
