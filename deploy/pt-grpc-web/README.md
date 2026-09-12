# PT gRPC-Web edge configuration

[`envoy.yaml`](envoy.yaml) is the audited local-development and CI edge for the
native C++ `mpmc.runtime.v1.PtFlashService` adapter:

```text
browser :8080 -> Envoy gRPC-Web/CORS -> native gRPC :50051
```

Envoy is used because gRPC-Web has a distinct browser transport and the upstream
project identifies Envoy as its official server-side proxy. The C++ process
therefore implements ordinary gRPC only; it does not duplicate HTTP, CORS, or
gRPC-Web framing code.

## Enforced edge policy

- only `/mpmc.runtime.v1.PtFlashService/` is routed;
- binary gRPC-Web request bodies are buffered and capped at 65,541 bytes: the
  64-KiB Protobuf limit plus the five-byte gRPC data-frame header;
- request headers are capped at 16 KiB;
- `grpc-timeout` is capped at 120 seconds and the route has a 125-second outer
  timeout;
- CORS accepts `POST` and preflight `OPTIONS` only;
- the checked-in allowlist contains only Vite development origins
  `http://localhost:5173` and `http://127.0.0.1:5173`;
- wildcard origins and credentials are disabled;
- JavaScript may read `grpc-status`, `grpc-message`, and
  `grpc-status-details-bin`.

The adapter independently configures gRPC's pre-deserialization receive limit,
post-decode size checks, response limit, memory quota, solve concurrency
gate, and deadline/cancellation publication checks. Defense at one layer does
not replace the other.

## Production deployment

`envoy.yaml` remains development-only. Production configuration is rendered by
[`render_production_config.py`](render_production_config.py), which refuses an
empty CORS allowlist, non-HTTPS or non-origin URLs, relative secret paths, and
invalid ports/domains. It configures:

- downstream TLS with required, CA-verified browser/client certificates;
- exact deployed frontend origins and one exact API virtual host (credentialed
  CORS is enabled only for those origins because client certificates are used);
- upstream HTTP/2 mutual TLS to the native process, including DNS SAN matching;
- active gRPC health checks for `mpmc.runtime.v1.PtFlashService`;
- JSON access logs without bodies or authorization/certificate fields; and
- an Envoy admin/Prometheus endpoint bound only to `127.0.0.1`.

Example (paths are inside the Envoy runtime):

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

The native process should bind loopback/private-only and trust the CA that issued
the Envoy client certificate. Public mTLS provides concrete client
authentication. Certificate enrollment, revocation, and authorization policy are
deployment responsibilities; no JWT issuer, JWKS URL, identity claims, or role
mapping is invented here. Add such a policy only when those authoritative inputs
exist. Do not replace the exact CORS entries with `*`.

CI validates this configuration and then exercises the generated TypeScript
gRPC-Web client against a synthetic C++ golden server. The golden backend is
explicitly a wire-structure fixture without EOS or physical-validation meaning.

## Versioned deployment bundle

[`production-deployment.template.json`](production-deployment.template.json)
is the v1 fail-closed input contract. Its `REQUIRED_*` values deliberately make
the checked-in template non-deployable. Copy it outside the source tree, replace
them with actual deployment values, and render into a new directory:

```bash
python3 deploy/pt-grpc-web/render_deployment_bundle.py \
  /secure/config/mpmc-pt-production.json \
  --output-directory /secure/rendered/mpmc-pt-v1
```

The renderer refuses unknown fields, a snapshot other than the compiled
`MPMC/PT/repository-curated-literature-snapshots/v1@r1`, a non-loopback native
listener, aliased CA roles, a changed secret root, placeholder identities,
non-HTTPS origins, and incomplete observability ownership. It emits:

- `envoy.yaml`, with the existing mTLS/CORS/health/resource policies;
- `pt-host-launch.json`, an argv array for `mpmc_pt_service_host`; and
- `deployment-metadata.json`, containing only non-secret identity, rotation
  and collection routing.

The fixed mount contract separates native-backend and edge material:

| Role | Read-only path |
| --- | --- |
| backend server leaf/key | `/run/secrets/mpmc-pt/backend/tls.crt`, `tls.key` |
| backend trust for edge clients | `/run/secrets/mpmc-pt/backend/edge-client-ca.pem` |
| public API leaf/key | `/run/secrets/mpmc-pt/edge/public.crt`, `public.key` |
| edge trust for browser clients | `/run/secrets/mpmc-pt/edge/browser-client-ca.pem` |
| edge upstream leaf/key | `/run/secrets/mpmc-pt/edge/upstream-client.crt`, `upstream-client.key` |
| edge trust for backend server | `/run/secrets/mpmc-pt/edge/backend-server-ca.pem` |

The four CA identifiers must be distinct. The edge-client CA is dedicated to
the Envoy-to-backend role; it must not issue unrelated client identities.
Rotation is explicitly
`overlapping-trust-bundles-and-rolling-restart`: publish old+new CA bundles,
rotate leaves and roll both processes, then remove the old CA in a later secret
generation. The manifest records the responsible owner, maximum leaf lifetime,
renewal lead time and CA overlap; issuance and automation remain external PKI
responsibilities.

Envoy metrics remain available only at
`http://127.0.0.1:9901/stats/prometheus`. Steady-state process and edge logs
go to stdout as JSON lines; process startup failures use stderr with the same
JSON discipline. The deployment manifest must name the local or sidecar agents
that scrape and forward those streams; the service does not embed a vendor
client or remote credentials.
