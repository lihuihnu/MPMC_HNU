import importlib.util
import pathlib
import sys
import unittest


MODULE_PATH = pathlib.Path(__file__).parents[1] / "render_production_config.py"
SPEC = importlib.util.spec_from_file_location("pt_edge_renderer", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
RENDERER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RENDERER
SPEC.loader.exec_module(RENDERER)


def valid_config(**changes):
    values = {
        "public_host": "api.example.test",
        "public_origins": ("https://ui.example.test",),
        "downstream_certificate_chain": "/run/secrets/public.crt",
        "downstream_private_key": "/run/secrets/public.key",
        "downstream_client_ca": "/run/secrets/browser-ca.crt",
        "upstream_host": "127.0.0.1",
        "upstream_port": 50051,
        "upstream_server_name": "pt-backend.internal",
        "upstream_client_certificate_chain": "/run/secrets/edge.crt",
        "upstream_client_private_key": "/run/secrets/edge.key",
        "upstream_server_ca": "/run/secrets/backend-ca.crt",
    }
    values.update(changes)
    return RENDERER.ProductionEdgeConfig(**values)


class ProductionConfigTest(unittest.TestCase):
    def test_fail_closed_tls_cors_health_and_observability(self):
        rendered = RENDERER.render_production_config(
            valid_config(
                public_origins=(
                    "https://ui.example.test",
                    "https://ops.example.test:4443",
                )
            )
        )
        self.assertIn("require_client_certificate: true", rendered)
        self.assertEqual(rendered.count("envoy.transport_sockets.tls"), 2)
        self.assertIn("exact: \"https://ui.example.test\"", rendered)
        self.assertIn("exact: \"https://ops.example.test:4443\"", rendered)
        self.assertNotIn("exact: \"*\"", rendered)
        self.assertNotIn("localhost:5173", rendered)
        self.assertIn("allow_credentials: true", rendered)
        self.assertIn("service_name: mpmc.runtime.v1.PtFlashService", rendered)
        self.assertIn("match_typed_subject_alt_names:", rendered)
        self.assertIn("address: 127.0.0.1\n      port_value: 9901", rendered)
        self.assertIn("json_format:", rendered)
        self.assertIn("grpc_status: \"%GRPC_STATUS%\"", rendered)
        self.assertNotIn("authorization", rendered.lower())
        self.assertNotIn("/model-api/", rendered)

    def test_model_routes_require_external_authorization_and_preserve_pt_route(self):
        rendered = RENDERER.render_production_config(
            valid_config(model_authz_loopback_port=10003)
        )
        self.assertIn(
            "prefix: /model-api/mpmc.model_configuration.v1.ModelSessionService/",
            rendered,
        )
        self.assertIn(
            "prefix_rewrite: /mpmc.model_configuration.v1.ModelSessionService/",
            rendered,
        )
        self.assertIn(
            "prefix: /model-api/mpmc.model_configuration.v1.ModelConfigurationService/",
            rendered,
        )
        self.assertIn(
            "prefix_rewrite: /mpmc.model_configuration.v1.ModelConfigurationService/",
            rendered,
        )
        self.assertIn("envoy.filters.http.ext_authz", rendered)
        self.assertIn("failure_mode_allow: false", rendered)
        self.assertIn("uri: http://127.0.0.1:10003", rendered)
        self.assertIn("- exact: authorization", rendered)
        self.assertIn("- exact: x-mpmc-model-session", rendered)
        self.assertIn("- exact: x-mpmc-authenticated-principal", rendered)
        self.assertIn(
            "allow_headers: \"content-type,x-grpc-web,grpc-timeout,x-user-agent,"
            "grpc-encoding,grpc-accept-encoding,authorization,x-mpmc-model-session\"",
            rendered,
        )
        self.assertNotIn(
            "grpc-accept-encoding,x-mpmc-authenticated-principal", rendered
        )
        self.assertIn("prefix: /mpmc.runtime.v1.PtFlashService/", rendered)
        self.assertIn("disabled: true", rendered)
        self.assertIn("cluster_name: web_identity_authz", rendered)

    def test_requires_https_origin_without_url_suffix(self):
        for origin in (
            "http://ui.example.test",
            "https://ui.example.test/",
            "https://ui.example.test/path",
            "https://user@ui.example.test",
            "https://localhost",
            "https://127.0.0.1",
        ):
            with self.subTest(origin=origin), self.assertRaises(ValueError):
                RENDERER.render_production_config(
                    valid_config(public_origins=(origin,))
                )

    def test_requires_explicit_origin_and_absolute_secret_paths(self):
        with self.assertRaises(ValueError):
            RENDERER.render_production_config(valid_config(public_origins=()))
        with self.assertRaises(ValueError):
            RENDERER.render_production_config(
                valid_config(downstream_private_key="relative.key")
            )

    def test_model_authz_port_must_be_valid_and_not_collide_with_edge_ports(self):
        for value in (0, 65536):
            with self.subTest(value=value), self.assertRaises(ValueError):
                RENDERER.render_production_config(
                    valid_config(model_authz_loopback_port=value)
                )
        for value in (8443, 9901):
            with self.subTest(value=value), self.assertRaises(ValueError):
                RENDERER.render_production_config(
                    valid_config(model_authz_loopback_port=value)
                )


if __name__ == "__main__":
    unittest.main()
