#!/usr/bin/env python3
"""Write fail-closed provenance for one Windows desktop MSI candidate."""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import re

FULL_SHA = re.compile(r"^[0-9a-f]{40}$")
EXPECTED_EXECUTABLE = "MPMC-PT-Desktop-Preview.exe"


def _read_object(path: pathlib.Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise RuntimeError(f"cannot read JSON object {path}: {error}") from error
    if not isinstance(value, dict):
        raise RuntimeError(f"JSON document is not an object: {path}")
    return value


def _sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def build_manifest(
    *,
    package_dir: pathlib.Path,
    artifact: pathlib.Path,
    expected_revision: str,
    workflow_run_id: int,
    workflow_run_attempt: int,
) -> dict[str, object]:
    if not FULL_SHA.fullmatch(expected_revision):
        raise RuntimeError("expected source revision must be one full lowercase SHA")
    if workflow_run_id < 1 or workflow_run_attempt < 1:
        raise RuntimeError("workflow run identity must be positive")

    package_dir = package_dir.resolve()
    if not package_dir.is_dir():
        raise RuntimeError("desktop package directory is missing")
    executable = package_dir / EXPECTED_EXECUTABLE
    desktop_manifest_path = package_dir / "resources" / "desktop-preview-manifest.json"
    native_manifest_path = (
        package_dir
        / "resources"
        / "desktop-native"
        / "share"
        / "mpmc-pt"
        / "product-staging-manifest.json"
    )
    dependency_manifest_path = (
        package_dir
        / "resources"
        / "desktop-native"
        / "share"
        / "mpmc-pt"
        / "dependency-manifest.json"
    )
    for required in (
        executable,
        desktop_manifest_path,
        native_manifest_path,
        dependency_manifest_path,
    ):
        if not required.is_file():
            raise RuntimeError(f"desktop package is incomplete: {required}")

    desktop = _read_object(desktop_manifest_path)
    native = _read_object(native_manifest_path)
    dependencies = _read_object(dependency_manifest_path)
    session = desktop.get("desktop_session")
    native_product = native.get("product")
    native_platform = native.get("platform")
    native_dependencies = native.get("dependencies")
    if not all(
        isinstance(value, dict)
        for value in (session, native_product, native_platform, native_dependencies)
    ):
        raise RuntimeError("desktop or native staging manifest is incomplete")
    assert isinstance(session, dict)
    assert isinstance(native_product, dict)
    assert isinstance(native_platform, dict)
    assert isinstance(native_dependencies, dict)

    if (
        desktop.get("convention") != "MPMC/PT/desktop-preview/v1"
        or desktop.get("source_revision") != expected_revision
        or desktop.get("platform") != "win32"
        or desktop.get("architecture") != "x64"
        or desktop.get("release_eligible") is not False
        or session.get("transport") != "native-grpc-loopback-bearer"
    ):
        raise RuntimeError("desktop package identity or session contract changed")
    if (
        native.get("convention") != "MPMC/PT/product-staging/v1"
        or native_product.get("build_revision") != expected_revision
        or native_platform.get("dependency_profile") != "windows-x86_64-release"
        or native_dependencies.get("provider") != "conan"
    ):
        raise RuntimeError("embedded native staging identity changed")
    if (
        dependencies.get("convention")
        != "MPMC/PT/conan-dependency-manifest/v1"
        or dependencies.get("profile") != "windows-x86_64-release"
        or dependencies.get("lock_sha256") != native_dependencies.get("lock_id")
    ):
        raise RuntimeError("embedded dependency provenance changed")

    artifact = artifact.resolve()
    expected_name = "mpmc-pt-desktop-0.1.0-windows-x86_64.msi"
    if artifact.name != expected_name or not artifact.is_file():
        raise RuntimeError(f"expected desktop MSI is missing: {expected_name}")
    if artifact.stat().st_size < 1024:
        raise RuntimeError("desktop MSI is unexpectedly small")

    return {
        "convention": "MPMC/PT/windows-desktop-installer-candidate/v1",
        "release_eligible": False,
        "source": {
            "revision": expected_revision,
            "workflow": "PT Windows desktop installer",
            "run_id": workflow_run_id,
            "run_attempt": workflow_run_attempt,
        },
        "desktop": {
            "executable": EXPECTED_EXECUTABLE,
            "executable_sha256": _sha256(executable),
            "manifest_sha256": _sha256(desktop_manifest_path),
            "electron_version": desktop.get("electron_version"),
            "session_transport": session.get("transport"),
        },
        "native": {
            "manifest_sha256": _sha256(native_manifest_path),
            "dependency_manifest_sha256": _sha256(dependency_manifest_path),
            "dependency_lock_id": native_dependencies.get("lock_id"),
        },
        "package": {
            "format": "msi",
            "file": artifact.name,
            "size_bytes": artifact.stat().st_size,
            "sha256": _sha256(artifact),
        },
        "signature": {
            "status": "unsigned-candidate",
            "verified": False,
            "required_for_release": "authenticode-payload-and-msi-rfc3161",
        },
    }


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-dir", required=True, type=pathlib.Path)
    parser.add_argument("--artifact", required=True, type=pathlib.Path)
    parser.add_argument("--expected-revision", required=True)
    parser.add_argument("--workflow-run-id", required=True, type=int)
    parser.add_argument("--workflow-run-attempt", required=True, type=int)
    parser.add_argument("--output", required=True, type=pathlib.Path)
    return parser.parse_args()


def main() -> int:
    args = _arguments()
    manifest = build_manifest(
        package_dir=args.package_dir,
        artifact=args.artifact,
        expected_revision=args.expected_revision,
        workflow_run_id=args.workflow_run_id,
        workflow_run_attempt=args.workflow_run_attempt,
    )
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        "WINDOWS_DESKTOP_INSTALLER_PROVENANCE_OK "
        f"package={manifest['package']['file']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
