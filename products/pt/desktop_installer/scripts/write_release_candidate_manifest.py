#!/usr/bin/env python3
"""Validate and write provenance for one signed Windows desktop release candidate."""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import re

FULL_SHA = re.compile(r"^[0-9a-f]{40}$")
GUID = re.compile(r"^[0-9A-F]{8}(?:-[0-9A-F]{4}){3}-[0-9A-F]{12}$")
DISPLAY_VERSION = re.compile(r"^[0-9]+\.[0-9]+\.[0-9]+-rc\.[1-9][0-9]*$")
MSI_VERSION = re.compile(r"^[0-9]+\.[0-9]+\.[0-9]+$")
EXPECTED_ROLES = {"desktop-executable", "native-host", "msi"}


def _read_object(path: pathlib.Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8-sig"))
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


def _identity(path: pathlib.Path) -> dict[str, object]:
    identity = _read_object(path)
    release = identity.get("release")
    application = identity.get("application")
    authenticode = identity.get("authenticode")
    if not all(isinstance(value, dict) for value in (release, application, authenticode)):
        raise RuntimeError("release identity structure is incomplete")
    assert isinstance(release, dict)
    assert isinstance(application, dict)
    assert isinstance(authenticode, dict)
    if (
        identity.get("convention") != "MPMC/PT/windows-desktop-release-identity/v1"
        or identity.get("product_id") != "mpmc-pt-desktop"
        or identity.get("platform") != "windows-x86_64"
        or release.get("channel") != "rc"
        or release.get("sequence") != 1
        or not isinstance(release.get("display_version"), str)
        or not DISPLAY_VERSION.fullmatch(str(release.get("display_version")))
        or not isinstance(release.get("msi_product_version"), str)
        or not MSI_VERSION.fullmatch(str(release.get("msi_product_version")))
        or not isinstance(release.get("product_guid"), str)
        or not GUID.fullmatch(str(release.get("product_guid")))
        or not isinstance(identity.get("upgrade_guid"), str)
        or not GUID.fullmatch(str(identity.get("upgrade_guid")))
        or authenticode.get("file_digest") != "SHA256"
        or authenticode.get("timestamp_protocol") != "RFC3161"
        or authenticode.get("timestamp_digest") != "SHA256"
        or set(authenticode.get("required_roles", [])) != EXPECTED_ROLES
    ):
        raise RuntimeError("fixed Windows desktop release identity is invalid")
    return identity


def build_manifest(
    *,
    identity_path: pathlib.Path,
    package_dir: pathlib.Path,
    artifact: pathlib.Path,
    signing_evidence_path: pathlib.Path,
    expected_revision: str,
    workflow_run_id: int,
    workflow_run_attempt: int,
) -> dict[str, object]:
    if not FULL_SHA.fullmatch(expected_revision):
        raise RuntimeError("expected source revision must be one full lowercase SHA")
    if workflow_run_id < 1 or workflow_run_attempt < 1:
        raise RuntimeError("workflow run identity must be positive")

    identity_path = identity_path.resolve()
    identity = _identity(identity_path)
    release = identity["release"]
    application = identity["application"]
    assert isinstance(release, dict)
    assert isinstance(application, dict)

    package_dir = package_dir.resolve()
    if not package_dir.is_dir():
        raise RuntimeError("desktop RC package directory is missing")
    desktop_executable = package_dir / str(application["release_executable"])
    native_host = (
        package_dir
        / "resources"
        / "desktop-native"
        / "bin"
        / "mpmc_pt_service_host.exe"
    )
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
        desktop_executable,
        native_host,
        desktop_manifest_path,
        native_manifest_path,
        dependency_manifest_path,
    ):
        if not required.is_file():
            raise RuntimeError(f"desktop RC payload is incomplete: {required}")

    desktop = _read_object(desktop_manifest_path)
    native = _read_object(native_manifest_path)
    dependencies = _read_object(dependency_manifest_path)
    session = desktop.get("desktop_session")
    release_identity = desktop.get("release_identity")
    desktop_signature = desktop.get("signature")
    native_product = native.get("product")
    native_platform = native.get("platform")
    native_dependencies = native.get("dependencies")
    if not all(
        isinstance(value, dict)
        for value in (
            session,
            release_identity,
            desktop_signature,
            native_product,
            native_platform,
            native_dependencies,
        )
    ):
        raise RuntimeError("desktop or native RC manifest is incomplete")
    assert isinstance(session, dict)
    assert isinstance(release_identity, dict)
    assert isinstance(desktop_signature, dict)
    assert isinstance(native_product, dict)
    assert isinstance(native_platform, dict)
    assert isinstance(native_dependencies, dict)

    if (
        desktop.get("convention")
        != "MPMC/PT/desktop-release-candidate-payload/v1"
        or desktop.get("source_revision") != expected_revision
        or desktop.get("platform") != "win32"
        or desktop.get("architecture") != "x64"
        or desktop.get("release_candidate") is not True
        or desktop.get("release_eligible") is not False
        or session.get("transport") != "native-grpc-loopback-bearer"
        or release_identity.get("display_version") != release.get("display_version")
        or release_identity.get("msi_product_version")
        != release.get("msi_product_version")
        or release_identity.get("product_guid") != release.get("product_guid")
        or release_identity.get("upgrade_guid") != identity.get("upgrade_guid")
        or release_identity.get("identity_sha256") != _sha256(identity_path)
        or desktop_signature.get("status") != "authenticode-signed-timestamped"
        or desktop_signature.get("verified") is not True
        or desktop_signature.get("file_digest") != "SHA256"
        or desktop_signature.get("timestamp_protocol") != "RFC3161"
        or desktop_signature.get("timestamp_digest") != "SHA256"
    ):
        raise RuntimeError("desktop RC payload identity or signature state changed")
    if (
        native.get("convention") != "MPMC/PT/product-staging/v1"
        or native_product.get("build_revision") != expected_revision
        or native_platform.get("dependency_profile") != "windows-x86_64-release"
        or native_dependencies.get("provider") != "conan"
        or dependencies.get("convention")
        != "MPMC/PT/conan-dependency-manifest/v1"
        or dependencies.get("profile") != "windows-x86_64-release"
        or dependencies.get("lock_sha256") != native_dependencies.get("lock_id")
    ):
        raise RuntimeError("embedded native/dependency provenance changed")

    artifact = artifact.resolve()
    expected_name = (
        f"mpmc-pt-desktop-{release['display_version']}-windows-x86_64.msi"
    )
    if artifact.name != expected_name or not artifact.is_file():
        raise RuntimeError(f"expected signed desktop RC MSI is missing: {expected_name}")
    if artifact.stat().st_size < 1024:
        raise RuntimeError("signed desktop RC MSI is unexpectedly small")

    evidence = _read_object(signing_evidence_path.resolve())
    if (
        evidence.get("convention") != "MPMC/PT/windows-authenticode-evidence/v1"
        or evidence.get("file_digest") != "SHA256"
        or evidence.get("timestamp_protocol") != "RFC3161"
        or evidence.get("timestamp_digest") != "SHA256"
        or not isinstance(evidence.get("signer_subject"), str)
        or not str(evidence.get("signer_subject")).strip()
        or desktop_signature.get("signer_subject") != evidence.get("signer_subject")
        or not isinstance(evidence.get("files"), list)
    ):
        raise RuntimeError("Authenticode evidence is incomplete or inconsistent")
    entries = {
        entry.get("role"): entry
        for entry in evidence["files"]
        if isinstance(entry, dict) and isinstance(entry.get("role"), str)
    }
    if set(entries) != EXPECTED_ROLES:
        raise RuntimeError("Authenticode evidence does not cover the full RC surface")
    expected_paths = {
        "desktop-executable": desktop_executable,
        "native-host": native_host,
        "msi": artifact,
    }
    for role, path in expected_paths.items():
        entry = entries[role]
        if (
            entry.get("status") != "Valid"
            or entry.get("timestamped") is not True
            or entry.get("signer_subject") != evidence.get("signer_subject")
            or entry.get("sha256") != _sha256(path)
        ):
            raise RuntimeError(f"Authenticode evidence is invalid for {role}")

    return {
        "convention": "MPMC/PT/windows-desktop-release-candidate/v1",
        "release_candidate_eligible": True,
        "release_eligible": False,
        "identity": {
            "product_id": identity["product_id"],
            "product_name": identity["product_name"],
            "display_version": release["display_version"],
            "msi_product_version": release["msi_product_version"],
            "release_channel": release["channel"],
            "release_sequence": release["sequence"],
            "product_guid": release["product_guid"],
            "upgrade_guid": identity["upgrade_guid"],
            "install_directory": identity["install_directory"],
            "identity_sha256": _sha256(identity_path),
        },
        "source": {
            "revision": expected_revision,
            "workflow": "PT Windows desktop release candidate",
            "run_id": workflow_run_id,
            "run_attempt": workflow_run_attempt,
        },
        "desktop": {
            "executable": application["release_executable"],
            "executable_sha256": _sha256(desktop_executable),
            "manifest_sha256": _sha256(desktop_manifest_path),
            "session_transport": session["transport"],
        },
        "native": {
            "host_sha256": _sha256(native_host),
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
        "authenticode": {
            "status": "signed-and-rfc3161-timestamped",
            "verified": True,
            "signer_subject": evidence["signer_subject"],
            "signer_thumbprint": evidence.get("signer_thumbprint"),
            "file_digest": "SHA256",
            "timestamp_protocol": "RFC3161",
            "timestamp_digest": "SHA256",
            "covered_roles": sorted(EXPECTED_ROLES),
        },
        "remaining_release_blocks": [
            "final-public-release-approval",
            "automatic-update-policy",
            "support-and-distribution-policy",
        ],
    }


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--identity", required=True, type=pathlib.Path)
    parser.add_argument("--package-dir", required=True, type=pathlib.Path)
    parser.add_argument("--artifact", required=True, type=pathlib.Path)
    parser.add_argument("--signing-evidence", required=True, type=pathlib.Path)
    parser.add_argument("--expected-revision", required=True)
    parser.add_argument("--workflow-run-id", required=True, type=int)
    parser.add_argument("--workflow-run-attempt", required=True, type=int)
    parser.add_argument("--output", required=True, type=pathlib.Path)
    return parser.parse_args()


def main() -> int:
    args = _arguments()
    manifest = build_manifest(
        identity_path=args.identity,
        package_dir=args.package_dir,
        artifact=args.artifact,
        signing_evidence_path=args.signing_evidence,
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
        "WINDOWS_DESKTOP_RELEASE_CANDIDATE_PROVENANCE_OK "
        f"package={manifest['package']['file']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
