#!/usr/bin/env python3
"""Validate one staged PT payload and write installer-candidate provenance."""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import pathlib
import re


@dataclasses.dataclass(frozen=True)
class Target:
    profile: str
    system_name: str
    processor: str
    entry: str
    package_format: str
    extension: str
    release_signature: str


TARGETS = {
    "linux-x86_64": Target(
        profile="linux-x86_64-release",
        system_name="Linux",
        processor="x86_64",
        entry="bin/mpmc_pt_service_host",
        package_format="deb",
        extension=".deb",
        release_signature="apt-archive-release-openpgp",
    ),
    "windows-x86_64": Target(
        profile="windows-x86_64-release",
        system_name="Windows",
        processor="AMD64",
        entry="bin/mpmc_pt_service_host.exe",
        package_format="msi",
        extension=".msi",
        release_signature="authenticode-payload-and-msi-rfc3161",
    ),
    "macos-arm64": Target(
        profile="macos-armv8-release",
        system_name="Darwin",
        processor="arm64",
        entry="bin/mpmc_pt_service_host",
        package_format="pkg",
        extension=".pkg",
        release_signature=(
            "developer-id-payload-and-installer-notarized-stapled"
        ),
    ),
}

SAFE_VERSION = re.compile(r"^[0-9]+[.][0-9]+[.][0-9]+(?:[-+][A-Za-z0-9.-]+)?$")
FULL_SHA = re.compile(r"^[0-9a-f]{40}$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
ALLOWED_TOP_LEVEL = {"bin", "lib", "share"}
EXPECTED_DEPENDENCIES = {
    "abseil",
    "c-ares",
    "grpc",
    "openssl",
    "protobuf",
    "re2",
    "zlib",
}
FORBIDDEN_SECRET_SUFFIXES = {".key", ".p12", ".pfx", ".jks", ".keystore"}
WINDOWS_ALLOWED_BINARIES = {
    "mpmc_pt_service_host.exe",
    "msvcp140.dll",
    "vcruntime140.dll",
    "vcruntime140_1.dll",
}


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


def _child(root: pathlib.Path, relative: str) -> pathlib.Path:
    relative_path = pathlib.PurePosixPath(relative)
    if relative_path.is_absolute() or ".." in relative_path.parts:
        raise RuntimeError(f"staging manifest contains unsafe path: {relative}")
    candidate = root.joinpath(*relative_path.parts)
    if not candidate.is_file():
        raise RuntimeError(f"staging manifest entry is missing: {relative}")
    return candidate


def validate_stage(
    stage: pathlib.Path, target_id: str, expected_revision: str
) -> tuple[dict[str, object], dict[str, object], pathlib.Path, pathlib.Path]:
    target = TARGETS[target_id]
    if not FULL_SHA.fullmatch(expected_revision):
        raise RuntimeError("expected source revision must be one full lowercase SHA")
    stage = stage.resolve()
    if not stage.is_dir():
        raise RuntimeError("staging input is not a directory")

    entries = list(stage.iterdir())
    unexpected_roots = sorted(path.name for path in entries if path.name not in ALLOWED_TOP_LEVEL)
    if unexpected_roots:
        raise RuntimeError(f"unexpected staging roots: {unexpected_roots}")
    for path in stage.rglob("*"):
        if path.is_symlink() and not path.resolve().is_relative_to(stage):
            raise RuntimeError(f"staging symlink escapes its root: {path}")
        if path.is_file() and path.suffix.casefold() in FORBIDDEN_SECRET_SUFFIXES:
            raise RuntimeError(f"staging tree contains a credential-like file: {path}")

    staging_path = stage / "share" / "mpmc-pt" / "product-staging-manifest.json"
    dependency_path = stage / "share" / "mpmc-pt" / "dependency-manifest.json"
    staging = _read_object(staging_path)
    dependencies = _read_object(dependency_path)
    product = staging.get("product")
    platform = staging.get("platform")
    dependency_ref = staging.get("dependencies")
    entry = staging.get("entry")
    snapshot = staging.get("snapshot_bundle")
    if not all(
        isinstance(value, dict)
        for value in (product, platform, dependency_ref, entry, snapshot)
    ):
        raise RuntimeError("staging manifest structure is incomplete")
    assert isinstance(product, dict)
    assert isinstance(platform, dict)
    assert isinstance(dependency_ref, dict)
    assert isinstance(entry, dict)
    assert isinstance(snapshot, dict)

    version = product.get("version")
    if (
        staging.get("convention") != "MPMC/PT/product-staging/v1"
        or product.get("id") != "mpmc-pt"
        or product.get("build_revision") != expected_revision
        or not isinstance(version, str)
        or not SAFE_VERSION.fullmatch(version)
    ):
        raise RuntimeError("staged product identity or source revision changed")
    if (
        platform.get("dependency_profile") != target.profile
        or platform.get("system_name") != target.system_name
        or platform.get("system_processor") != target.processor
    ):
        raise RuntimeError("staged platform does not match the installer target")
    if (
        dependency_ref.get("provider") != "conan"
        or dependency_ref.get("lock_kind") != "conan-lock-sha256"
        or not isinstance(dependency_ref.get("lock_id"), str)
        or not SHA256.fullmatch(str(dependency_ref.get("lock_id")))
        or dependencies.get("convention")
        != "MPMC/PT/conan-dependency-manifest/v1"
        or dependencies.get("provider") != "conan"
        or dependencies.get("profile") != target.profile
        or dependencies.get("lock_sha256") != dependency_ref.get("lock_id")
        or dependencies.get("client_version")
        != dependency_ref.get("provider_version")
    ):
        raise RuntimeError("staged dependency provenance is inconsistent")
    packages = dependencies.get("packages")
    if not isinstance(packages, list):
        raise RuntimeError("staged dependency package inventory is missing")
    indexed_packages = {
        package.get("name"): package
        for package in packages
        if isinstance(package, dict) and isinstance(package.get("name"), str)
    }
    if set(indexed_packages) != EXPECTED_DEPENDENCIES:
        raise RuntimeError("staged dependency package inventory changed")
    for name, package in indexed_packages.items():
        if not all(
            isinstance(package.get(field), str) and package.get(field)
            for field in ("reference", "package_id", "package_revision", "license")
        ):
            raise RuntimeError(f"staged dependency provenance is incomplete: {name}")
    if (
        entry.get("transport") != "native-grpc-mtls"
        or entry.get("loopback_only") is not True
        or entry.get("requires_explicit_tls_paths") is not True
        or entry.get("path") != target.entry
    ):
        raise RuntimeError("staged process-boundary contract changed")
    _child(stage, str(entry["path"]))
    if (
        snapshot.get("id")
        != "MPMC/PT/repository-curated-literature-snapshots/v1"
        or snapshot.get("revision") != "r1"
    ):
        raise RuntimeError("staged parameter snapshot identity changed")

    if target_id == "windows-x86_64":
        binary_names = {
            path.name.casefold()
            for path in (stage / "bin").iterdir()
            if path.is_file()
        }
        unexpected = sorted(binary_names - WINDOWS_ALLOWED_BINARIES)
        if unexpected:
            raise RuntimeError(
                "Windows staging contains unreviewed or system binaries: "
                f"{unexpected}"
            )

    return staging, dependencies, staging_path, dependency_path


def build_manifest(
    *,
    stage: pathlib.Path,
    artifact: pathlib.Path,
    target_id: str,
    expected_revision: str,
    workflow_run_id: int,
    workflow_run_attempt: int,
) -> dict[str, object]:
    target = TARGETS[target_id]
    if workflow_run_id < 1 or workflow_run_attempt < 1:
        raise RuntimeError("workflow run identity must be positive")
    staging, dependencies, staging_path, dependency_path = validate_stage(
        stage, target_id, expected_revision
    )
    product = staging["product"]
    dependency_ref = staging["dependencies"]
    entry = staging["entry"]
    snapshot = staging["snapshot_bundle"]
    assert isinstance(product, dict)
    assert isinstance(dependency_ref, dict)
    assert isinstance(entry, dict)
    assert isinstance(snapshot, dict)

    version = str(product["version"])
    expected_name = f"mpmc-pt-{version}-{target_id}{target.extension}"
    artifact = artifact.resolve()
    if artifact.name != expected_name or not artifact.is_file():
        raise RuntimeError(f"expected installer artifact is missing: {expected_name}")
    if artifact.stat().st_size < 1024:
        raise RuntimeError("installer artifact is unexpectedly small")

    return {
        "convention": "MPMC/PT/installer-candidate/v1",
        "release_eligible": False,
        "product": {
            "id": "mpmc-pt",
            "version": version,
            "build_revision": expected_revision,
        },
        "source": {
            "workflow": "PT signed installer gate",
            "run_id": workflow_run_id,
            "run_attempt": workflow_run_attempt,
        },
        "platform": {
            "target": target_id,
            "dependency_profile": target.profile,
        },
        "payload": {
            "entry": entry["path"],
            "transport": entry["transport"],
            "staging_manifest_sha256": _sha256(staging_path),
            "dependency_manifest_sha256": _sha256(dependency_path),
            "dependency_lock_id": dependency_ref["lock_id"],
            "dependency_package_count": len(dependencies.get("packages", [])),
            "snapshot_id": snapshot["id"],
            "snapshot_revision": snapshot["revision"],
        },
        "package": {
            "format": target.package_format,
            "file": artifact.name,
            "size_bytes": artifact.stat().st_size,
            "sha256": _sha256(artifact),
        },
        "signature": {
            "status": "unsigned-candidate",
            "verified": False,
            "required_for_release": target.release_signature,
        },
    }


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, type=pathlib.Path)
    parser.add_argument("--artifact", required=True, type=pathlib.Path)
    parser.add_argument("--target", required=True, choices=sorted(TARGETS))
    parser.add_argument("--expected-revision", required=True)
    parser.add_argument("--workflow-run-id", required=True, type=int)
    parser.add_argument("--workflow-run-attempt", required=True, type=int)
    parser.add_argument("--output", required=True, type=pathlib.Path)
    return parser.parse_args()


def main() -> int:
    args = _arguments()
    manifest = build_manifest(
        stage=args.stage,
        artifact=args.artifact,
        target_id=args.target,
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
        "INSTALLER_CANDIDATE_PROVENANCE_OK "
        f"target={args.target} package={manifest['package']['file']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
