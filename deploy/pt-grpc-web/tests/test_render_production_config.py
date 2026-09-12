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


if __name__ == "__main__":
    unittest.main()
