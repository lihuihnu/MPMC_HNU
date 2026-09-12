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
post-decode size checks, response limit, memory/thread quota, solve concurrency
gate, and deadline/cancellation publication checks. Defense at one layer does
not replace the other.

## Deployment boundary

This file binds cleartext local ports and is not a production TLS/auth policy.
Before deployment, replace the development CORS entries with explicit deployed
frontend origins, terminate TLS, add authentication/authorization, bind the
native gRPC port to a private interface, and preserve the size/deadline filters.
Do not replace the allowlist with `*` when credentials or authorization headers
are enabled.

CI validates this configuration and then exercises the generated TypeScript
gRPC-Web client against a synthetic C++ golden server. The golden backend is
explicitly a wire-structure fixture without EOS or physical-validation meaning.
