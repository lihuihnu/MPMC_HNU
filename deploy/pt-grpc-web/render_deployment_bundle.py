#!/usr/bin/env python3
"""Render a versioned, non-secret PT production deployment bundle."""

from __future__ import annotations

import argparse
import json
import pathlib
import re
from dataclasses import dataclass
from typing import Any

from render_production_config import (
    ProductionEdgeConfig,
    _validate_dns_name,
    _validate_origin,
    _validate_public_dns_name,
    render_production_config,
)


DEPLOYMENT_CONVENTION = "MPMC/PT/production-deployment/v1"
SNAPSHOT_BUNDLE_ID = "MPMC/PT/repository-curated-literature-snapshots/v1"
SNAPSHOT_BUNDLE_REVISION = "r1"
SECRET_MOUNT_ROOT = "/run/secrets/mpmc-pt"
ROTATION_MODE = "overlapping-trust-bundles-and-rolling-restart"
PROMETHEUS_SCRAPE_URL = "http://127.0.0.1:9901/stats/prometheus"
JSON_LOG_TARGET = "stdout-json-lines"
_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:/@+-]{0,255}$")


def _object(value: Any, field: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{field} must be an object")
    return value


def _exact_keys(
    value: dict[str, Any], expected: set[str], field: str
) -> dict[str, Any]:
    missing = expected - value.keys()
    extra = value.keys() - expected
    if missing or extra:
        raise ValueError(
            f"{field} fields mismatch; missing={sorted(missing)}, "
            f"unexpected={sorted(extra)}"
        )
    return value


def _identifier(value: Any, field: str) -> str:
    if not isinstance(value, str) or not _IDENTIFIER.fullmatch(value):
        raise ValueError(f"{field} must be an explicit stable identifier")
    if value.startswith("REQUIRED_"):
        raise ValueError(f"{field} still contains a required-value marker")
    return value


def _text(value: Any, field: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{field} must be nonempty text")
    return value


def _positive_days(value: Any, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{field} must be an integer day count")
    if not 1 <= value <= 3660:
        raise ValueError(f"{field} must be in [1,3660]")
    return value


def _loopback_listener(value: Any) -> tuple[str, int]:
    if not isinstance(value, str) or value.count(":") != 1:
        raise ValueError(
            "network.backend_listen_address must be 127.0.0.1:PORT"
        )
    host, port_text = value.split(":", 1)
    if host != "127.0.0.1":
        raise ValueError(
            "network.backend_listen_address must remain loopback-only"
        )
    try:
        port = int(port_text)
    except ValueError as error:
        raise ValueError(
            "network.backend_listen_address has an invalid port"
        ) from error
    if not 1 <= port <= 65535:
        raise ValueError(
            "network.backend_listen_address port must be in [1,65535]"
        )
    return host, port


@dataclass(frozen=True)
class DeploymentBundle:
    deployment_id: str
    api_dns: str
    frontend_origins: tuple[str, ...]
    backend_listen_address: str
    backend_server_dns: str
    ca_ids: tuple[str, str, str, str]
    edge_client_identity: str
    secret_generation: str
    certificate_rotation_owner: str
    leaf_max_validity_days: int
    leaf_renew_before_days: int
    ca_overlap_days: int
    prometheus_collector_id: str
    json_log_collector_id: str

    @classmethod
    def from_mapping(cls, document: dict[str, Any]) -> "DeploymentBundle":
        top = _exact_keys(
            _object(document, "deployment"),
            {
                "convention",
                "deployment_id",
                "snapshot_bundle",
                "network",
                "identity",
                "observability",
            },
            "deployment",
        )
        if top["convention"] != DEPLOYMENT_CONVENTION:
            raise ValueError("unsupported deployment convention")

        snapshot = _exact_keys(
            _object(top["snapshot_bundle"], "snapshot_bundle"),
            {"id", "revision", "location"},
            "snapshot_bundle",
        )
        if (
            snapshot["id"] != SNAPSHOT_BUNDLE_ID
            or snapshot["revision"] != SNAPSHOT_BUNDLE_REVISION
            or snapshot["location"]
            != "builtin://mpmc/pt/repository-curated-literature-snapshots/v1"
        ):
            raise ValueError(
                "deployment snapshot bundle does not match the compiled host"
            )

        network = _exact_keys(
            _object(top["network"], "network"),
            {
                "api_dns",
                "frontend_origins",
                "backend_listen_address",
                "backend_server_dns",
            },
            "network",
        )
        api_dns = _validate_public_dns_name(
            _text(network["api_dns"], "network.api_dns"),
            "network.api_dns",
        )
        if not isinstance(network["frontend_origins"], list):
            raise ValueError("network.frontend_origins must be an array")
        origins = tuple(
            dict.fromkeys(
                _validate_origin(
                    _text(origin, "network.frontend_origins[]")
                )
                for origin in network["frontend_origins"]
            )
        )
        if not origins:
            raise ValueError(
                "network.frontend_origins requires at least one exact origin"
            )
        _, backend_port = _loopback_listener(network["backend_listen_address"])
        backend_dns = _validate_dns_name(
            _text(
                network["backend_server_dns"],
                "network.backend_server_dns",
            ),
            "network.backend_server_dns",
        )
        if backend_dns.startswith("required-"):
            raise ValueError(
                "network.backend_server_dns still contains a required-value marker"
            )

        identity = _exact_keys(
            _object(top["identity"], "identity"),
            {
                "secret_mount_root",
                "public_server_ca_id",
                "browser_client_ca_id",
                "backend_server_ca_id",
                "edge_client_ca_id",
                "edge_client_identity",
                "secret_generation",
                "rotation_mode",
                "certificate_rotation_owner",
                "leaf_max_validity_days",
                "leaf_renew_before_days",
                "ca_overlap_days",
            },
            "identity",
        )
        if identity["secret_mount_root"] != SECRET_MOUNT_ROOT:
            raise ValueError(
                f"identity.secret_mount_root must be {SECRET_MOUNT_ROOT}"
            )
        if identity["rotation_mode"] != ROTATION_MODE:
            raise ValueError("unsupported identity.rotation_mode")
        ca_ids = tuple(
            _identifier(identity[name], f"identity.{name}")
            for name in (
                "public_server_ca_id",
                "browser_client_ca_id",
                "backend_server_ca_id",
                "edge_client_ca_id",
            )
        )
        if len(set(ca_ids)) != len(ca_ids):
            raise ValueError("the four CA roles require distinct identifiers")
        max_validity = _positive_days(
            identity["leaf_max_validity_days"],
            "identity.leaf_max_validity_days",
        )
        renew_before = _positive_days(
            identity["leaf_renew_before_days"],
            "identity.leaf_renew_before_days",
        )
        overlap = _positive_days(
            identity["ca_overlap_days"], "identity.ca_overlap_days"
        )
        if renew_before >= max_validity:
            raise ValueError(
                "leaf renewal must begin before the maximum validity interval"
            )

        observability = _exact_keys(
            _object(top["observability"], "observability"),
            {
                "prometheus_scrape_url",
                "process_log_target",
                "edge_log_target",
                "prometheus_collector_id",
                "json_log_collector_id",
            },
            "observability",
        )
        if observability["prometheus_scrape_url"] != PROMETHEUS_SCRAPE_URL:
            raise ValueError(
                "observability.prometheus_scrape_url must remain loopback-only"
            )
        if (
            observability["process_log_target"] != JSON_LOG_TARGET
            or observability["edge_log_target"] != JSON_LOG_TARGET
        ):
            raise ValueError("production logs must use stdout JSON lines")

        return cls(
            deployment_id=_identifier(top["deployment_id"], "deployment_id"),
            api_dns=api_dns,
            frontend_origins=origins,
            backend_listen_address=f"127.0.0.1:{backend_port}",
            backend_server_dns=backend_dns,
            ca_ids=ca_ids,
            edge_client_identity=_identifier(
                identity["edge_client_identity"],
                "identity.edge_client_identity",
            ),
            secret_generation=_identifier(
                identity["secret_generation"], "identity.secret_generation"
            ),
            certificate_rotation_owner=_identifier(
                identity["certificate_rotation_owner"],
                "identity.certificate_rotation_owner",
            ),
            leaf_max_validity_days=max_validity,
            leaf_renew_before_days=renew_before,
            ca_overlap_days=overlap,
            prometheus_collector_id=_identifier(
                observability["prometheus_collector_id"],
                "observability.prometheus_collector_id",
            ),
            json_log_collector_id=_identifier(
                observability["json_log_collector_id"],
                "observability.json_log_collector_id",
            ),
        )

    def edge_config(self) -> ProductionEdgeConfig:
        backend_port = int(self.backend_listen_address.rsplit(":", 1)[1])
        return ProductionEdgeConfig(
            public_host=self.api_dns,
            public_origins=self.frontend_origins,
            downstream_certificate_chain=(
                f"{SECRET_MOUNT_ROOT}/edge/public.crt"
            ),
            downstream_private_key=f"{SECRET_MOUNT_ROOT}/edge/public.key",
            downstream_client_ca=(
                f"{SECRET_MOUNT_ROOT}/edge/browser-client-ca.pem"
            ),
            upstream_host="127.0.0.1",
            upstream_port=backend_port,
            upstream_server_name=self.backend_server_dns,
            upstream_client_certificate_chain=(
                f"{SECRET_MOUNT_ROOT}/edge/upstream-client.crt"
            ),
            upstream_client_private_key=(
                f"{SECRET_MOUNT_ROOT}/edge/upstream-client.key"
            ),
            upstream_server_ca=(
                f"{SECRET_MOUNT_ROOT}/edge/backend-server-ca.pem"
            ),
        )

    def host_launch(self) -> dict[str, Any]:
        return {
            "executable": "mpmc_pt_service_host",
            "argv": [
                "--listen-address",
                self.backend_listen_address,
                "--tls-certificate-chain",
                f"{SECRET_MOUNT_ROOT}/backend/tls.crt",
                "--tls-private-key",
                f"{SECRET_MOUNT_ROOT}/backend/tls.key",
                "--trusted-edge-client-ca",
                f"{SECRET_MOUNT_ROOT}/backend/edge-client-ca.pem",
            ],
            "snapshot_bundle": {
                "id": SNAPSHOT_BUNDLE_ID,
                "revision": SNAPSHOT_BUNDLE_REVISION,
            },
        }

    def metadata(self) -> dict[str, Any]:
        return {
            "convention": DEPLOYMENT_CONVENTION,
            "deployment_id": self.deployment_id,
            "api_dns": self.api_dns,
            "frontend_origins": list(self.frontend_origins),
            "identity": {
                "public_server_ca_id": self.ca_ids[0],
                "browser_client_ca_id": self.ca_ids[1],
                "backend_server_ca_id": self.ca_ids[2],
                "edge_client_ca_id": self.ca_ids[3],
                "edge_client_identity": self.edge_client_identity,
                "secret_generation": self.secret_generation,
                "rotation_mode": ROTATION_MODE,
                "certificate_rotation_owner": self.certificate_rotation_owner,
                "leaf_max_validity_days": self.leaf_max_validity_days,
                "leaf_renew_before_days": self.leaf_renew_before_days,
                "ca_overlap_days": self.ca_overlap_days,
            },
            "observability": {
                "prometheus_scrape_url": PROMETHEUS_SCRAPE_URL,
                "prometheus_collector_id": self.prometheus_collector_id,
                "process_log_target": JSON_LOG_TARGET,
                "edge_log_target": JSON_LOG_TARGET,
                "json_log_collector_id": self.json_log_collector_id,
            },
        }


def render_bundle(document: dict[str, Any], output_directory: pathlib.Path) -> None:
    deployment = DeploymentBundle.from_mapping(document)
    output_directory.mkdir(parents=True, exist_ok=False)
    (output_directory / "envoy.yaml").write_text(
        render_production_config(deployment.edge_config()), encoding="utf-8"
    )
    (output_directory / "pt-host-launch.json").write_text(
        json.dumps(deployment.host_launch(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (output_directory / "deployment-metadata.json").write_text(
        json.dumps(deployment.metadata(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=pathlib.Path)
    parser.add_argument("--output-directory", required=True, type=pathlib.Path)
    return parser


def main() -> int:
    args = _parser().parse_args()
    try:
        document = json.loads(args.manifest.read_text(encoding="utf-8"))
        render_bundle(document, args.output_directory)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        _parser().error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
