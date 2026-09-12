# PT process composition root and secure host

`mpmc::pt_process` is the executable-application boundary around the frozen PT
service contract. It owns configured backend object graphs, registers their
`PtFlashBackend` views in one `PtService`, exposes that service through one gRPC
adapter, and hosts the adapter with mandatory mutual TLS.

## Assembly boundary

The PR76, SW92 Profile-C, and CPA factories accept parameter sets that were
already created and provenance-validated by the caller. Each factory only invokes
the established model/evaluator/backend constructors and retains their complete
lifetime. It does not parse or invent parameters, choose a model from a request,
call a solver, alter tolerances, or reinterpret a result. This repository has no
audited production parameter database, so no test fixture or implicit default
dataset is wired into the host.

```cpp
std::vector<mpmc::pt_process::OwnedConfiguredPtBackend> backends;
backends.push_back(mpmc::pt_process::make_pr76_configured_backend(
    "pr76.production-r1", pr76_parameters));
backends.push_back(mpmc::pt_process::make_sw92_configured_backend(
    "sw92.production-r1", sw92_parameters, sw92_options));
backends.push_back(mpmc::pt_process::make_cpa_configured_backend(
    "cpa.production-r1", cpa_parameters));

auto observer = std::make_shared<mpmc::pt_process::PtJsonLineObserver>(std::clog);
mpmc::pt_process::PtCompositionRoot root(std::move(backends), {}, {}, observer);

mpmc::pt_process::PtProcessHostOptions host_options;
host_options.listen_address = "127.0.0.1:50051";
host_options.tls = mpmc::pt_process::load_pt_process_tls_identity({
    "/run/secrets/pt-backend.crt",
    "/run/secrets/pt-backend.key",
    "/run/secrets/edge-client-ca.crt"});
mpmc::pt_process::PtProcessHost host(root, std::move(host_options));
host.start();
host.wait();
```

The default adapter admits one solve at a time because the currently configured
PR76 and CPA evaluators own reusable sequential scratch. Raising concurrency is
permitted only after every reachable backend object graph has been audited for
concurrent calls.

## Security and operations

The native host has no insecure credential mode. It requires a server certificate,
private key, and trusted client CA, then asks gRPC to require and verify every
client certificate. In the recommended deployment, it listens on loopback and
trusts only the Envoy edge client CA. Certificate issuance and identity
authorization remain deployment policy; this module does not fabricate a JWT
issuer, JWKS, user directory, or role model.

The default gRPC health service publishes both the overall server and
`mpmc.runtime.v1.PtFlashService`. Graceful shutdown marks both non-serving before
the bounded shutdown grace period. The JSON-lines observer records method,
transport completion, byte counts, duration, and the service outcome. In
particular, scientific `indeterminate` and `service_error` are different values;
no pressure, temperature, composition, backend parameters, diagnostic text, or
PEM content is logged.

Use the production Envoy renderer in
[`deploy/pt-grpc-web`](../../deploy/pt-grpc-web/README.md) for public TLS,
client-certificate authentication, exact deployed-origin CORS, active native
gRPC health checks, JSON access logs, and loopback-only Prometheus/admin access.
