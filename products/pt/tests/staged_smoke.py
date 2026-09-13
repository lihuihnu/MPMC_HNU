#!/usr/bin/env python3
"""Install and exercise the relocatable PT product staging tree."""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import shutil
import signal
import subprocess
import tempfile
import time


def _run(command: list[str], *, timeout: int = 120) -> None:
    subprocess.run(command, check=True, timeout=timeout)


def _issue_ca(
    openssl: str, directory: pathlib.Path, name: str
) -> tuple[pathlib.Path, pathlib.Path]:
    certificate = directory / f"{name}.crt"
    key = directory / f"{name}.key"
    configuration = directory / f"{name}.cnf"
    configuration.write_text(
        "[req]\n"
        "distinguished_name=subject\n"
        "x509_extensions=ca_extensions\n"
        "prompt=no\n"
        "[subject]\n"
        f"CN={name}\n"
        "[ca_extensions]\n"
        "basicConstraints=critical,CA:TRUE\n"
        "keyUsage=critical,keyCertSign,cRLSign\n"
        "subjectKeyIdentifier=hash\n",
        encoding="ascii",
    )
    _run(
        [
            openssl,
            "req",
            "-x509",
            "-newkey",
            "rsa:2048",
            "-nodes",
            "-sha256",
            "-days",
            "1",
            "-config",
            str(configuration),
            "-keyout",
            str(key),
            "-out",
            str(certificate),
        ]
    )
    return certificate, key


def _issue_identity(
    openssl: str,
    directory: pathlib.Path,
    name: str,
    common_name: str,
    authority: tuple[pathlib.Path, pathlib.Path],
    extended_key_usage: str,
    subject_alternative_name: str,
) -> tuple[pathlib.Path, pathlib.Path]:
    certificate = directory / f"{name}.crt"
    key = directory / f"{name}.key"
    request = directory / f"{name}.csr"
    request_configuration = directory / f"{name}.cnf"
    extensions = directory / f"{name}.ext"
    request_configuration.write_text(
        "[req]\n"
        "distinguished_name=subject\n"
        "prompt=no\n"
        "[subject]\n"
        f"CN={common_name}\n",
        encoding="ascii",
    )
    extensions.write_text(
        "basicConstraints=critical,CA:FALSE\n"
        "keyUsage=critical,digitalSignature,keyEncipherment\n"
        f"extendedKeyUsage={extended_key_usage}\n"
        f"subjectAltName={subject_alternative_name}\n",
        encoding="ascii",
    )
    _run(
        [
            openssl,
            "req",
            "-newkey",
            "rsa:2048",
            "-nodes",
            "-sha256",
            "-config",
            str(request_configuration),
            "-keyout",
            str(key),
            "-out",
            str(request),
        ]
    )
    _run(
        [
            openssl,
            "x509",
            "-req",
            "-in",
            str(request),
            "-CA",
            str(authority[0]),
            "-CAkey",
            str(authority[1]),
            "-CAcreateserial",
            "-days",
            "1",
            "-sha256",
            "-extfile",
            str(extensions),
            "-out",
            str(certificate),
        ]
    )
    return certificate, key


def _wait_for_ready(
    process: subprocess.Popen[str], output_path: pathlib.Path, timeout: float
) -> int:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if output_path.exists():
            for line in output_path.read_text(encoding="utf-8").splitlines():
                try:
                    event = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if event.get("event") == "pt_process_ready":
                    port = event.get("selected_port")
                    if isinstance(port, int) and 1 <= port <= 65535:
                        return port
                    raise RuntimeError("staged host published an invalid port")
        if process.poll() is not None:
            raise RuntimeError("staged host exited before readiness")
        time.sleep(0.05)
    raise RuntimeError("timed out waiting for staged host readiness")


def _request_stop(process: subprocess.Popen[str]) -> None:
    if os.name == "nt":
        process.send_signal(signal.CTRL_BREAK_EVENT)
    else:
        process.send_signal(signal.SIGTERM)


def _logs(output: pathlib.Path, error: pathlib.Path) -> str:
    stdout = output.read_text(encoding="utf-8") if output.exists() else ""
    stderr = error.read_text(encoding="utf-8") if error.exists() else ""
    return f"\n--- host stdout ---\n{stdout}\n--- host stderr ---\n{stderr}"


def _validate_stage(stage: pathlib.Path, host_relative_path: pathlib.Path) -> pathlib.Path:
    manifest_path = stage / "share" / "mpmc-pt" / "product-staging-manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("convention") != "MPMC/PT/product-staging/v1":
        raise RuntimeError("staged product manifest convention changed")
    if manifest.get("entry", {}).get("path") != host_relative_path.as_posix():
        raise RuntimeError("staged product entry path does not match the manifest")
    snapshot = manifest.get("snapshot_bundle", {})
    if (
        snapshot.get("id")
        != "MPMC/PT/repository-curated-literature-snapshots/v1"
        or snapshot.get("revision") != "r1"
    ):
        raise RuntimeError("staged snapshot identity changed")
    dependencies = manifest.get("dependencies", {})
    dependency_manifest = json.loads(
        (stage / "share" / "mpmc-pt" / "dependency-manifest.json").read_text(
            encoding="utf-8"
        )
    )
    if (
        dependencies.get("provider") != "conan"
        or dependency_manifest.get("provider") != "conan"
        or dependency_manifest.get("convention")
        != "MPMC/PT/conan-dependency-manifest/v1"
        or dependencies.get("provider_version")
        != dependency_manifest.get("client_version")
        or dependencies.get("lock_kind") != "conan-lock-sha256"
        or dependencies.get("lock_id")
        != dependency_manifest.get("lock_sha256")
        or manifest.get("platform", {}).get("dependency_profile")
        != dependency_manifest.get("profile")
    ):
        raise RuntimeError("staged dependency provenance changed")
    packages = {
        package.get("name"): package
        for package in dependency_manifest.get("packages", [])
        if isinstance(package, dict)
    }
    if set(packages) != {
        "abseil",
        "c-ares",
        "grpc",
        "openssl",
        "protobuf",
        "re2",
        "zlib",
    }:
        raise RuntimeError("staged dependency package set changed")
    for name, package in packages.items():
        if not all(
            package.get(field)
            for field in ("reference", "package_id", "package_revision", "license")
        ):
            raise RuntimeError(f"staged {name} provenance is incomplete")
    serialized = manifest_path.read_text(encoding="utf-8")
    if str(stage.parent) in serialized:
        raise RuntimeError("staged product manifest contains a build-tree path")

    for package in packages:
        copyright_path = (
            stage / "share" / "mpmc-pt" / "third-party" / package / "copyright"
        )
        if not copyright_path.is_file() or copyright_path.stat().st_size == 0:
            raise RuntimeError(f"missing staged {package} copyright file")

    host = stage / host_relative_path
    if not host.is_file():
        raise RuntimeError("staged PT host executable is missing")
    return host


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cmake", required=True)
    parser.add_argument("--build-directory", required=True, type=pathlib.Path)
    parser.add_argument("--stage-directory", required=True, type=pathlib.Path)
    parser.add_argument("--configuration", required=True)
    parser.add_argument("--host-relative-path", required=True, type=pathlib.Path)
    parser.add_argument("--probe", required=True, type=pathlib.Path)
    parser.add_argument("--openssl", required=True)
    return parser


def main() -> int:
    args = _parser().parse_args()
    build = args.build_directory.resolve()
    stage = args.stage_directory.resolve()
    if stage.parent != build or stage.name != "stage":
        raise RuntimeError("staging test may only replace BUILD_DIRECTORY/stage")
    if stage.exists():
        shutil.rmtree(stage)
    _run(
        [
            args.cmake,
            "--install",
            str(build),
            "--config",
            args.configuration,
            "--prefix",
            str(stage),
        ]
    )
    host = _validate_stage(stage, args.host_relative_path)

    with tempfile.TemporaryDirectory(prefix="mpmc-pt-staged-") as raw_temp:
        temporary = pathlib.Path(raw_temp).resolve()
        backend_ca = _issue_ca(args.openssl, temporary, "backend-test-ca")
        edge_ca = _issue_ca(args.openssl, temporary, "edge-test-ca")
        backend = _issue_identity(
            args.openssl,
            temporary,
            "backend",
            "pt-backend.product.test",
            backend_ca,
            "serverAuth",
            "DNS:pt-backend.product.test",
        )
        edge = _issue_identity(
            args.openssl,
            temporary,
            "edge",
            "pt-product-probe",
            edge_ca,
            "clientAuth",
            "URI:spiffe://mpmc.test/product-probe",
        )
        output_path = temporary / "host.stdout.jsonl"
        error_path = temporary / "host.stderr.jsonl"
        creation_flags = (
            subprocess.CREATE_NEW_PROCESS_GROUP if os.name == "nt" else 0
        )
        with output_path.open("w", encoding="utf-8") as output, error_path.open(
            "w", encoding="utf-8"
        ) as error:
            process = subprocess.Popen(
                [
                    str(host),
                    "--listen-address",
                    "127.0.0.1:0",
                    "--tls-certificate-chain",
                    str(backend[0]),
                    "--tls-private-key",
                    str(backend[1]),
                    "--trusted-edge-client-ca",
                    str(edge_ca[0]),
                    "--shutdown-grace-seconds",
                    "5",
                ],
                cwd=stage,
                stdout=output,
                stderr=error,
                text=True,
                creationflags=creation_flags,
            )
            try:
                port = _wait_for_ready(process, output_path, 15.0)
                probe = subprocess.run(
                    [
                        str(args.probe.resolve()),
                        f"127.0.0.1:{port}",
                        "pt-backend.product.test",
                        str(backend_ca[0]),
                        str(edge[0]),
                        str(edge[1]),
                    ],
                    capture_output=True,
                    text=True,
                    timeout=30,
                )
                if (
                    probe.returncode != 0
                    or "STAGED_PRODUCT_OK" not in probe.stdout
                ):
                    raise RuntimeError(
                        "staged product probe failed\n"
                        f"probe stdout: {probe.stdout}\n"
                        f"probe stderr: {probe.stderr}"
                    )
                _request_stop(process)
                if process.wait(timeout=15) != 0:
                    raise RuntimeError("staged host returned a nonzero status")
            except Exception as error_value:
                if process.poll() is None:
                    process.kill()
                    process.wait(timeout=10)
                raise RuntimeError(
                    str(error_value) + _logs(output_path, error_path)
                ) from error_value

        events = [
            json.loads(line)
            for line in output_path.read_text(encoding="utf-8").splitlines()
            if line.startswith("{")
        ]
        if not any(event.get("event") == "pt_process_stopped" for event in events):
            raise RuntimeError("staged host did not publish graceful shutdown")

    print("STAGED_PRODUCT_SMOKE_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
