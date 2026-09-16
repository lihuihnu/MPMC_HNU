#!/usr/bin/env python3
"""Render the fail-closed production Envoy edge for PT and hosted model services."""

from __future__ import annotations

import argparse
import ipaddress
import json
import os
import re
from dataclasses import dataclass
from urllib.parse import urlsplit


_HOST_LABEL = re.compile(r"^[a-z0-9](?:[a-z0-9-]{0,61}[a-z0-9])?$")


def _yaml(value: str) -> str:
    return json.dumps(value, ensure_ascii=True)


def _validate_port(value: int, field: str) -> int:
    if not 1 <= value <= 65535:
        raise ValueError(f"{field} must be in [1, 65535]")
    return value


def _validate_dns_name(value: str, field: str) -> str:
    candidate = value.rstrip(".").lower()
    if not candidate or len(candidate) > 253 or any(
        not _HOST_LABEL.fullmatch(label) for label in candidate.split(".")
    ):
        raise ValueError(f"{field} must be an ASCII DNS name")
    return candidate


def _validate_public_dns_name(value: str, field: str) -> str:
    candidate = _validate_dns_name(value, field)
    try:
        ipaddress.ip_address(candidate)
    except ValueError:
        pass
    else:
        raise ValueError(f"{field} must be a DNS name, not an IP address")
    if "." not in candidate or candidate == "localhost":
        raise ValueError(f"{field} must be an explicit deployed DNS name")
    return candidate


def _validate_socket_host(value: str, field: str) -> str:
    candidate = value.strip().lower()
    try:
        return str(ipaddress.ip_address(candidate))
    except ValueError:
        return _validate_dns_name(candidate, field)


def _validate_origin(value: str) -> str:
    parsed = urlsplit(value)
    if (
        parsed.scheme != "https"
        or not parsed.hostname
        or parsed.username is not None
        or parsed.password is not None
        or parsed.path
        or parsed.query
        or parsed.fragment
    ):
        raise ValueError(
            "public origins must be exact HTTPS origins without path, query, or fragment"
        )
    host = _validate_public_dns_name(parsed.hostname, "public origin host")
    try:
        port = parsed.port
    except ValueError as error:
        raise ValueError("public origin has an invalid port") from error
    canonical = f"https://{host}"
    if port is not None:
        canonical += f":{_validate_port(port, 'public origin port')}"
    return canonical


def _validate_absolute_path(value: str, field: str) -> str:
    if (
        not os.path.isabs(value)
        or "\x00" in value
        or any(ord(character) < 0x20 for character in value)
    ):
        raise ValueError(f"{field} must be an absolute non-control path")
    return value


@dataclass(frozen=True)
class ProductionEdgeConfig:
    public_host: str
    public_origins: tuple[str, ...]
    downstream_certificate_chain: str
    downstream_private_key: str
    downstream_client_ca: str
    upstream_host: str
    upstream_port: int
    upstream_server_name: str
    upstream_client_certificate_chain: str
    upstream_client_private_key: str
    upstream_server_ca: str
    public_port: int = 8443
    admin_port: int = 9901
    model_authz_loopback_port: int | None = None

    def validated(self) -> "ProductionEdgeConfig":
        public_host = _validate_public_dns_name(self.public_host, "public host")
        origins = tuple(dict.fromkeys(_validate_origin(x) for x in self.public_origins))
        if not origins:
            raise ValueError("at least one deployed public origin is required")
        upstream_host = _validate_socket_host(self.upstream_host, "upstream host")
        model_authz_port = (
            None
            if self.model_authz_loopback_port is None
            else _validate_port(
                self.model_authz_loopback_port, "model authz loopback port"
            )
        )
        public_port = _validate_port(self.public_port, "public port")
        admin_port = _validate_port(self.admin_port, "admin port")
        if model_authz_port is not None and model_authz_port in {
            public_port,
            admin_port,
        }:
            raise ValueError(
                "model authz loopback port must differ from public/admin ports"
            )
        return ProductionEdgeConfig(
            public_host=public_host,
            public_origins=origins,
            downstream_certificate_chain=_validate_absolute_path(
                self.downstream_certificate_chain, "downstream certificate chain"
            ),
            downstream_private_key=_validate_absolute_path(
                self.downstream_private_key, "downstream private key"
            ),
            downstream_client_ca=_validate_absolute_path(
                self.downstream_client_ca, "downstream client CA"
            ),
            upstream_host=upstream_host,
            upstream_port=_validate_port(self.upstream_port, "upstream port"),
            upstream_server_name=_validate_dns_name(
                self.upstream_server_name, "upstream server name"
            ),
            upstream_client_certificate_chain=_validate_absolute_path(
                self.upstream_client_certificate_chain,
                "upstream client certificate chain",
            ),
            upstream_client_private_key=_validate_absolute_path(
                self.upstream_client_private_key, "upstream client private key"
            ),
            upstream_server_ca=_validate_absolute_path(
                self.upstream_server_ca, "upstream server CA"
            ),
            public_port=public_port,
            admin_port=admin_port,
            model_authz_loopback_port=model_authz_port,
        )


def _model_routes(enabled: bool) -> str:
    if not enabled:
        return ""
    return """                        - match:
                            prefix: /model-api/mpmc.model_configuration.v1.ModelSessionService/
                          route:
                            cluster: pt_service_native_grpc
                            prefix_rewrite: /mpmc.model_configuration.v1.ModelSessionService/
                            timeout: 1805s
                            max_stream_duration:
                              grpc_timeout_header_max: 1800s
                        - match:
                            prefix: /model-api/mpmc.model_configuration.v1.ModelConfigurationService/
                          route:
                            cluster: pt_service_native_grpc
                            prefix_rewrite: /mpmc.model_configuration.v1.ModelConfigurationService/
                            timeout: 125s
                            max_stream_duration:
                              grpc_timeout_header_max: 120s
"""


def _ext_authz_filter(config: ProductionEdgeConfig) -> str:
    if config.model_authz_loopback_port is None:
        return ""
    return f"""                  - name: envoy.filters.http.ext_authz
                    typed_config:
                      \"@type\": type.googleapis.com/envoy.extensions.filters.http.ext_authz.v3.ExtAuthz
                      failure_mode_allow: false
                      status_on_error:
                        code: ServiceUnavailable
                      http_service:
                        server_uri:
                          uri: http://127.0.0.1:{config.model_authz_loopback_port}
                          cluster: web_identity_authz
                          timeout: 1s
                        authorization_request:
                          allowed_headers:
                            patterns:
                              - exact: authorization
                              - exact: x-mpmc-model-session
                        authorization_response:
                          allowed_upstream_headers:
                            patterns:
                              - exact: x-mpmc-authenticated-principal
                          allowed_client_headers:
                            patterns:
                              - exact: www-authenticate
"""


def _authz_cluster(config: ProductionEdgeConfig) -> str:
    if config.model_authz_loopback_port is None:
        return ""
    return f"""
    - name: web_identity_authz
      connect_timeout: 0.25s
      type: STATIC
      lb_policy: ROUND_ROBIN
      load_assignment:
        cluster_name: web_identity_authz
        endpoints:
          - lb_endpoints:
              - endpoint:
                  address:
                    socket_address:
                      address: 127.0.0.1
                      port_value: {config.model_authz_loopback_port}
"""


def render_production_config(config: ProductionEdgeConfig) -> str:
    config = config.validated()
    model_enabled = config.model_authz_loopback_port is not None
    origins = "\n".join(
        f"                          - exact: {_yaml(origin)}"
        for origin in config.public_origins
    )
    allow_headers = (
        "content-type,x-grpc-web,grpc-timeout,x-user-agent,"
        "grpc-encoding,grpc-accept-encoding"
    )
    if model_enabled:
        allow_headers += ",authorization,x-mpmc-model-session"
    legacy_filter_config = (
        """                          typed_per_filter_config:
                            envoy.filters.http.ext_authz:
                              \"@type\": type.googleapis.com/envoy.extensions.filters.http.ext_authz.v3.ExtAuthzPerRoute
                              disabled: true
"""
        if model_enabled
        else ""
    )
    return f"""admin:
  access_log_path: /dev/null
  address:
    socket_address:
      address: 127.0.0.1
      port_value: {config.admin_port}

static_resources:
  listeners:
    - name: pt_grpc_web_production
      address:
        socket_address:
          address: 0.0.0.0
          port_value: {config.public_port}
      filter_chains:
        - transport_socket:
            name: envoy.transport_sockets.tls
            typed_config:
              \"@type\": type.googleapis.com/envoy.extensions.transport_sockets.tls.v3.DownstreamTlsContext
              common_tls_context:
                tls_params:
                  tls_minimum_protocol_version: TLSv1_2
                alpn_protocols: [\"h2\", \"http/1.1\"]
                tls_certificates:
                  - certificate_chain:
                      filename: {_yaml(config.downstream_certificate_chain)}
                    private_key:
                      filename: {_yaml(config.downstream_private_key)}
                validation_context:
                  trusted_ca:
                    filename: {_yaml(config.downstream_client_ca)}
              require_client_certificate: true
          filters:
            - name: envoy.filters.network.http_connection_manager
              typed_config:
                \"@type\": type.googleapis.com/envoy.extensions.filters.network.http_connection_manager.v3.HttpConnectionManager
                codec_type: AUTO
                stat_prefix: pt_grpc_web
                max_request_headers_kb: 16
                access_log:
                  - name: envoy.access_loggers.stdout
                    typed_config:
                      \"@type\": type.googleapis.com/envoy.extensions.access_loggers.stream.v3.StdoutAccessLog
                      log_format:
                        json_format:
                          event: pt_edge_request
                          start_time: \"%START_TIME%\"
                          request_id: \"%REQ(X-REQUEST-ID)%\"
                          method: \"%REQ(:METHOD)%\"
                          path: \"%REQ(X-ENVOY-ORIGINAL-PATH?:PATH)%\"
                          response_code: \"%RESPONSE_CODE%\"
                          grpc_status: \"%GRPC_STATUS%\"
                          duration_ms: \"%DURATION%\"
                          request_bytes: \"%BYTES_RECEIVED%\"
                          response_bytes: \"%BYTES_SENT%\"
                          response_flags: \"%RESPONSE_FLAGS%\"
                route_config:
                  name: pt_service_routes
                  virtual_hosts:
                    - name: pt_service
                      domains: [{_yaml(config.public_host)}]
                      cors:
                        allow_origin_string_match:
{origins}
                        allow_methods: \"POST, OPTIONS\"
                        allow_headers: \"{allow_headers}\"
                        expose_headers: \"grpc-status,grpc-message,grpc-status-details-bin\"
                        max_age: \"600\"
                        allow_credentials: true
                      routes:
{_model_routes(model_enabled)}                        - match:
                            prefix: /mpmc.runtime.v1.PtFlashService/
{legacy_filter_config}                          route:
                            cluster: pt_service_native_grpc
                            timeout: 125s
                            max_stream_duration:
                              grpc_timeout_header_max: 120s
                http_filters:
                  - name: envoy.filters.http.cors
                    typed_config:
                      \"@type\": type.googleapis.com/envoy.extensions.filters.http.cors.v3.Cors
{_ext_authz_filter(config)}                  - name: envoy.filters.http.buffer
                    typed_config:
                      \"@type\": type.googleapis.com/envoy.extensions.filters.http.buffer.v3.Buffer
                      max_request_bytes: 65541
                  - name: envoy.filters.http.grpc_web
                    typed_config:
                      \"@type\": type.googleapis.com/envoy.extensions.filters.http.grpc_web.v3.GrpcWeb
                  - name: envoy.filters.http.router
                    typed_config:
                      \"@type\": type.googleapis.com/envoy.extensions.filters.http.router.v3.Router

  clusters:
    - name: pt_service_native_grpc
      connect_timeout: 0.25s
      type: STATIC
      lb_policy: ROUND_ROBIN
      health_checks:
        - timeout: 1s
          interval: 5s
          unhealthy_threshold: 2
          healthy_threshold: 1
          grpc_health_check:
            service_name: {pt_service_name()}
      typed_extension_protocol_options:
        envoy.extensions.upstreams.http.v3.HttpProtocolOptions:
          \"@type\": type.googleapis.com/envoy.extensions.upstreams.http.v3.HttpProtocolOptions
          explicit_http_config:
            http2_protocol_options: {{}}
      load_assignment:
        cluster_name: pt_service_native_grpc
        endpoints:
          - lb_endpoints:
              - endpoint:
                  address:
                    socket_address:
                      address: {_yaml(config.upstream_host)}
                      port_value: {config.upstream_port}
      transport_socket:
        name: envoy.transport_sockets.tls
        typed_config:
          \"@type\": type.googleapis.com/envoy.extensions.transport_sockets.tls.v3.UpstreamTlsContext
          sni: {_yaml(config.upstream_server_name)}
          common_tls_context:
            tls_params:
              tls_minimum_protocol_version: TLSv1_2
            alpn_protocols: [\"h2\"]
            tls_certificates:
              - certificate_chain:
                  filename: {_yaml(config.upstream_client_certificate_chain)}
                private_key:
                  filename: {_yaml(config.upstream_client_private_key)}
            validation_context:
              trusted_ca:
                filename: {_yaml(config.upstream_server_ca)}
              match_typed_subject_alt_names:
                - san_type: DNS
                  matcher:
                    exact: {_yaml(config.upstream_server_name)}
{_authz_cluster(config)}"""


def pt_service_name() -> str:
    return "mpmc.runtime.v1.PtFlashService"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--public-host", required=True)
    parser.add_argument("--public-origin", action="append", required=True)
    parser.add_argument("--downstream-certificate-chain", required=True)
    parser.add_argument("--downstream-private-key", required=True)
    parser.add_argument("--downstream-client-ca", required=True)
    parser.add_argument("--upstream-host", default="127.0.0.1")
    parser.add_argument("--upstream-port", type=int, default=50051)
    parser.add_argument("--upstream-server-name", required=True)
    parser.add_argument("--upstream-client-certificate-chain", required=True)
    parser.add_argument("--upstream-client-private-key", required=True)
    parser.add_argument("--upstream-server-ca", required=True)
    parser.add_argument("--public-port", type=int, default=8443)
    parser.add_argument("--admin-port", type=int, default=9901)
    parser.add_argument("--model-authz-loopback-port", type=int)
    return parser


def main() -> int:
    args = _parser().parse_args()
    try:
        rendered = render_production_config(
            ProductionEdgeConfig(
                public_host=args.public_host,
                public_origins=tuple(args.public_origin),
                downstream_certificate_chain=args.downstream_certificate_chain,
                downstream_private_key=args.downstream_private_key,
                downstream_client_ca=args.downstream_client_ca,
                upstream_host=args.upstream_host,
                upstream_port=args.upstream_port,
                upstream_server_name=args.upstream_server_name,
                upstream_client_certificate_chain=args.upstream_client_certificate_chain,
                upstream_client_private_key=args.upstream_client_private_key,
                upstream_server_ca=args.upstream_server_ca,
                public_port=args.public_port,
                admin_port=args.admin_port,
                model_authz_loopback_port=args.model_authz_loopback_port,
            )
        )
    except ValueError as error:
        _parser().error(str(error))
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
