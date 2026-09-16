import copy
import importlib.util
import json
import pathlib
import sys
import tempfile
import unittest


DEPLOY_DIRECTORY = pathlib.Path(__file__).parents[1]
sys.path.insert(0, str(DEPLOY_DIRECTORY))
MODULE_PATH = DEPLOY_DIRECTORY / "render_deployment_bundle.py"
SPEC = importlib.util.spec_from_file_location(
    "pt_deployment_renderer", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
RENDERER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RENDERER
SPEC.loader.exec_module(RENDERER)


def valid_manifest():
    return {
        "convention": "MPMC/PT/production-deployment/v1",
        "deployment_id": "production-us-east-1",
        "snapshot_bundle": {
            "id": "MPMC/PT/repository-curated-literature-snapshots/v1",
            "revision": "r1",
            "location": (
                "builtin://mpmc/pt/"
                "repository-curated-literature-snapshots/v1"
            ),
        },
        "network": {
            "api_dns": "api.example.test",
            "frontend_origins": ["https://ui.example.test"],
            "backend_listen_address": "127.0.0.1:50051",
            "backend_server_dns": "pt-backend.internal",
            "model_authz_loopback_port": 10003,
        },
        "identity": {
            "secret_mount_root": "/run/secrets/mpmc-pt",
            "public_server_ca_id": "pki/public-server-ca/v7",
            "browser_client_ca_id": "pki/browser-client-ca/v4",
            "backend_server_ca_id": "pki/backend-server-ca/v3",
            "edge_client_ca_id": "pki/edge-client-ca/v5",
            "edge_client_identity": "envoy-edge/production-us-east-1",
            "secret_generation": "2026-09-12.1",
            "rotation_mode": (
                "overlapping-trust-bundles-and-rolling-restart"
            ),
            "certificate_rotation_owner": "platform-security",
            "leaf_max_validity_days": 90,
            "leaf_renew_before_days": 30,
            "ca_overlap_days": 14,
        },
        "observability": {
            "prometheus_scrape_url": (
                "http://127.0.0.1:9901/stats/prometheus"
            ),
            "process_log_target": "stdout-json-lines",
            "edge_log_target": "stdout-json-lines",
            "prometheus_collector_id": "prometheus/production-primary",
            "json_log_collector_id": "logs/production-primary",
        },
    }


class DeploymentBundleTest(unittest.TestCase):
    def test_renders_fail_closed_hosted_model_edge_and_observability_contract(self):
        document = valid_manifest()
        with tempfile.TemporaryDirectory() as temporary:
            output = pathlib.Path(temporary) / "rendered"
            RENDERER.render_bundle(document, output)
            envoy = (output / "envoy.yaml").read_text(encoding="utf-8")
            launch = json.loads(
                (output / "pt-host-launch.json").read_text(encoding="utf-8")
            )
            metadata = json.loads(
                (output / "deployment-metadata.json").read_text(
                    encoding="utf-8"
                )
            )

        self.assertIn('domains: ["api.example.test"]', envoy)
        self.assertIn('exact: "https://ui.example.test"', envoy)
        self.assertIn("require_client_certificate: true", envoy)
        self.assertIn("envoy.filters.http.ext_authz", envoy)
        self.assertIn("uri: http://127.0.0.1:10003", envoy)
        self.assertIn(
            "prefix: /model-api/mpmc.model_configuration.v1.ModelSessionService/",
            envoy,
        )
        self.assertIn(
            "prefix: /model-api/mpmc.model_configuration.v1.ModelConfigurationService/",
            envoy,
        )
        self.assertEqual(
            launch["snapshot_bundle"],
            {
                "id": (
                    "MPMC/PT/repository-curated-literature-snapshots/v1"
                ),
                "revision": "r1",
            },
        )
        self.assertEqual(launch["executable"], "mpmc_pt_service_host")
        self.assertIn(
            "/run/secrets/mpmc-pt/backend/edge-client-ca.pem",
            launch["argv"],
        )
        self.assertIn("--enable-model-sessions", launch["argv"])
        self.assertEqual(
            metadata["hosted_model"],
            {"enabled": True, "authz_loopback_port": 10003},
        )
        self.assertEqual(
            metadata["observability"]["prometheus_scrape_url"],
            "http://127.0.0.1:9901/stats/prometheus",
        )
        serialized = json.dumps({"launch": launch, "metadata": metadata})
        self.assertNotIn("PRIVATE KEY", serialized)
        self.assertNotIn("certificate_payload", serialized)

    def test_legacy_manifest_without_model_authz_remains_pt_only(self):
        document = valid_manifest()
        del document["network"]["model_authz_loopback_port"]
        deployment = RENDERER.DeploymentBundle.from_mapping(document)
        envoy = RENDERER.render_production_config(deployment.edge_config())
        launch = deployment.host_launch()
        metadata = deployment.metadata()
        self.assertNotIn("/model-api/", envoy)
        self.assertNotIn("envoy.filters.http.ext_authz", envoy)
        self.assertNotIn("--enable-model-sessions", launch["argv"])
        self.assertEqual(
            metadata["hosted_model"],
            {"enabled": False, "authz_loopback_port": None},
        )

    def test_rejects_snapshot_identity_and_trust_domain_aliasing(self):
        wrong_snapshot = valid_manifest()
        wrong_snapshot["snapshot_bundle"]["revision"] = "latest"
        with self.assertRaises(ValueError):
            RENDERER.DeploymentBundle.from_mapping(wrong_snapshot)

        aliased_ca = valid_manifest()
        aliased_ca["identity"]["edge_client_ca_id"] = (
            aliased_ca["identity"]["backend_server_ca_id"]
        )
        with self.assertRaises(ValueError):
            RENDERER.DeploymentBundle.from_mapping(aliased_ca)

    def test_rejects_nonloopback_backend_unknown_fields_and_authz_collision(self):
        public_backend = valid_manifest()
        public_backend["network"]["backend_listen_address"] = "0.0.0.0:50051"
        with self.assertRaises(ValueError):
            RENDERER.DeploymentBundle.from_mapping(public_backend)

        typo = valid_manifest()
        typo["observability"]["prometheus_target"] = "ignored"
        with self.assertRaises(ValueError):
            RENDERER.DeploymentBundle.from_mapping(typo)

        wrong_domain_type = valid_manifest()
        wrong_domain_type["network"]["api_dns"] = 42
        with self.assertRaises(ValueError):
            RENDERER.DeploymentBundle.from_mapping(wrong_domain_type)

        collision = valid_manifest()
        collision["network"]["model_authz_loopback_port"] = 50051
        with self.assertRaises(ValueError):
            RENDERER.DeploymentBundle.from_mapping(collision)

        invalid_authz = valid_manifest()
        invalid_authz["network"]["model_authz_loopback_port"] = 0
        with self.assertRaises(ValueError):
            RENDERER.DeploymentBundle.from_mapping(invalid_authz)

    def test_checked_in_template_is_deliberately_not_deployable(self):
        template = json.loads(
            (DEPLOY_DIRECTORY / "production-deployment.template.json").read_text(
                encoding="utf-8"
            )
        )
        with self.assertRaises(ValueError):
            RENDERER.DeploymentBundle.from_mapping(template)

    def test_refuses_to_overwrite_an_existing_output_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = pathlib.Path(temporary) / "rendered"
            output.mkdir()
            with self.assertRaises(FileExistsError):
                RENDERER.render_bundle(copy.deepcopy(valid_manifest()), output)


if __name__ == "__main__":
    unittest.main()
