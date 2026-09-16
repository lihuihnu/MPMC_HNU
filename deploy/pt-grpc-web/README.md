# PT gRPC-Web edge configuration

[`envoy.yaml`](envoy.yaml) is the audited local-development and CI edge for the
native C++ `mpmc.runtime.v1.PtFlashService` adapter:

```text
browser :8080 -> Envoy gRPC-Web/CORS -> native gRPC :50051
```

Envoy owns browser transport concerns; the C++ process implements ordinary gRPC
and does not duplicate HTTP, CORS or gRPC-Web framing.

## Development edge policy

The checked-in `envoy.yaml` remains deliberately narrow and PT-only:

- only `/mpmc.runtime.v1.PtFlashService/` is routed;
- binary gRPC-Web request bodies are buffered and capped at 65,541 bytes: the
  64-KiB Protobuf limit plus the five-byte gRPC data-frame header;
- request headers are capped at 16 KiB;
- `grpc-timeout` is capped at 120 seconds and the route has a 125-second outer
  timeout;
- CORS accepts `POST` and preflight `OPTIONS` only;
- the allowlist contains only Vite development origins
  `http://localhost:5173` and `http://127.0.0.1:5173`;
- wildcard origins and credentials are disabled; and
- JavaScript may read `grpc-status`, `grpc-message`, and
  `grpc-status-details-bin`.

The adapter independently configures gRPC's pre-deserialization receive limit,
post-decode size checks, response limit, memory quota, solve concurrency gate,
and deadline/cancellation publication checks. Defense at one layer does not
replace the other.

## Production edge

`envoy.yaml` is development-only. Production configuration is rendered by
[`render_production_config.py`](render_production_config.py), which refuses an
empty CORS allowlist, non-HTTPS or non-origin URLs, relative secret paths, invalid
ports/domains and conflicting control-plane ports.

All production variants configure:

- downstream TLS with required, CA-verified client certificates;
- exact deployed frontend origins and one exact API virtual host;
- upstream HTTP/2 mutual TLS to the native process, including DNS SAN matching;
- active gRPC health checks for `mpmc.runtime.v1.PtFlashService`;
- JSON access logs without bodies, bearer values or certificate payloads; and
- an Envoy admin/Prometheus endpoint bound only to `127.0.0.1`.

### PT-only default

Without `--model-authz-loopback-port`, the renderer preserves the historical
production contract: only `PtFlashService` is exposed. No `/model-api/` route,
external-authorization filter or model-session host enablement is implied.

```bash
python3 deploy/pt-grpc-web/render_production_config.py \
  --public-host api.example.com \
  --public-origin https://app.example.com \
  --downstream-certificate-chain /run/secrets/public.crt \
  --downstream-private-key /run/secrets/public.key \
  --downstream-client-ca /run/secrets/browser-client-ca.crt \
  --upstream-server-name pt-backend.internal \
  --upstream-client-certificate-chain /run/secrets/envoy-client.crt \
  --upstream-client-private-key /run/secrets/envoy-client.key \
  --upstream-server-ca /run/secrets/backend-server-ca.crt \
  > /etc/envoy/envoy.yaml
```

### Hosted-model opt-in

Hosted Classic PR model ownership is enabled only when the deployment explicitly
supplies `--model-authz-loopback-port PORT`. That switch adds exactly the model
routes used by the browser client:

```text
/model-api/mpmc.model_configuration.v1.ModelSessionService/
/model-api/mpmc.model_configuration.v1.ModelConfigurationService/
```

and installs a fail-closed Envoy `ext_authz` filter backed by a loopback-only
HTTP authorization service. The model authorization request forwards only the
browser bearer/session inputs required by that policy. A successful authorization
response may provide one internal `x-mpmc-authenticated-principal` header for the
native upstream. The legacy `PtFlashService` route explicitly disables this model
authorization filter and keeps its previous behavior.

The browser is not trusted to choose the native principal. The native production
listener still requires the dedicated edge-client mTLS certificate. Only after
that peer is verified may the host consume a single bounded, character-validated
`x-mpmc-authenticated-principal` value and bind the model session to
`web:<principal>`. Duplicate or invalid values fail closed. If no internal
principal is present, the existing direct-mTLS certificate identity contract is
retained for native/host regressions and non-Web callers.

A model-enabled native host must be launched with `--enable-model-sessions`.
The external authorization service remains a deployment responsibility: the
repository does not invent a JWT issuer, JWKS URL, account database, roles or a
third-party identity provider. Certificate enrollment, revocation and PKI
rotation are also deployment responsibilities. Do not replace exact CORS origins
with `*`, and do not expose the native listener directly to browsers.

## Hosted Web product regression

The `Hosted Web model product` workflow verifies the real product path on a
GitHub-hosted Ubuntu runner:

```text
Chrome
  -> same-origin test hosting shell
  -> production-style Envoy gRPC-Web edge
  -> external authz contract
  -> edge-client mTLS
  -> native ModelSessionService / ModelConfigurationService
  -> shared Classic PR workspace
```

The test authorization sidecar recognizes only two synthetic opaque repository
tokens and returns stable test principals. It is **not** a production identity
provider. The browser regression deliberately injects a spoofed internal-principal
header, verifies another principal cannot use the first principal's session, and
then exercises shared Classic PR apply, native solve and release. The success
marker is:

```text
HOSTED_WEB_PRODUCT_SMOKE_OK shared_classic_pr_ui=true cross_principal_denied=true release=true
```

The same workflow validates the model-enabled Envoy configuration with the pinned
Envoy image before running Chrome. Existing PT-only renderer tests and the normal
PT gRPC process-adapter workflow remain separate regression gates.

## Versioned deployment bundle

[`production-deployment.template.json`](production-deployment.template.json) is
the v1 fail-closed input contract. Its `REQUIRED_*` values deliberately make the
checked-in template non-deployable. Copy it outside the source tree, replace them
with actual deployment values, and render into a new directory:

```bash
python3 deploy/pt-grpc-web/render_deployment_bundle.py \
  /secure/config/mpmc-pt-production.json \
  --output-directory /secure/rendered/mpmc-pt-v1
```

The renderer refuses unknown fields, a snapshot other than the compiled
`MPMC/PT/repository-curated-literature-snapshots/v1@r1`, a non-loopback native
listener, aliased CA roles, a changed secret root, placeholder identities,
non-HTTPS origins, incomplete observability ownership, invalid authz ports and an
authz port colliding with the native listener.

`network.model_authz_loopback_port` is optional for compatibility:

- omitted: the bundle remains PT-only and the host argv does not contain
  `--enable-model-sessions`;
- present: the bundle renders the two model routes plus external authorization,
  adds `--enable-model-sessions` to `pt-host-launch.json`, and records
  `hosted_model.enabled=true` plus the non-secret authz loopback port in
  `deployment-metadata.json`.

The checked-in template includes the field so the hosted-model production shape is
visible, but it remains non-deployable until all `REQUIRED_*` values are replaced
with authoritative deployment inputs.

The bundle emits:

- `envoy.yaml`, with the mTLS/CORS/health/resource policy and optional hosted-model
  authorization routes;
- `pt-host-launch.json`, an argv array for `mpmc_pt_service_host`; and
- `deployment-metadata.json`, containing only non-secret identity, hosted-model,
  rotation and collection routing.

The fixed mount contract separates native-backend and edge material:

| Role | Read-only path |
| --- | --- |
| backend server leaf/key | `/run/secrets/mpmc-pt/backend/tls.crt`, `tls.key` |
| backend trust for edge clients | `/run/secrets/mpmc-pt/backend/edge-client-ca.pem` |
| public API leaf/key | `/run/secrets/mpmc-pt/edge/public.crt`, `public.key` |
| edge trust for browser clients | `/run/secrets/mpmc-pt/edge/browser-client-ca.pem` |
| edge upstream leaf/key | `/run/secrets/mpmc-pt/edge/upstream-client.crt`, `upstream-client.key` |
| edge trust for backend server | `/run/secrets/mpmc-pt/edge/backend-server-ca.pem` |

The four CA identifiers must be distinct. The edge-client CA is dedicated to the
Envoy-to-backend role and must not issue unrelated client identities. Rotation is
`overlapping-trust-bundles-and-rolling-restart`: publish old+new CA bundles,
rotate leaves and roll both processes, then remove the old CA in a later secret
generation. The manifest records the responsible owner, maximum leaf lifetime,
renewal lead time and CA overlap; issuance and automation remain external PKI
responsibilities.

Envoy metrics remain available only at
`http://127.0.0.1:9901/stats/prometheus`. Steady-state process and edge logs go to
stdout as JSON lines; process startup failures use stderr with the same JSON
discipline. The deployment manifest must name the local or sidecar agents that
scrape and forward those streams; the service does not embed a vendor client or
remote credentials.
