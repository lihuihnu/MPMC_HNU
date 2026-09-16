# Authenticated Web model session connector

This module connects an **authoritative browser identity provider** to the existing
`ModelSessionClient`. It does not make Electron credentials usable from Web code,
and it does not itself authenticate or authorize a browser principal at the
network edge.

## Identity contract

The caller supplies a `WebIdentityProvider` that returns one opaque access token.
The token is treated strictly as bearer material:

- the browser client does not decode JWT claims or make authorization decisions;
- the token is never stored by this module in `localStorage`, session IDs, model
  objects, DOM state or user-facing errors;
- provider exception text is sanitized;
- missing, malformed or already-expired credentials fail closed;
- identity acquisition is cancellation-bounded even if the provider ignores its
  `AbortSignal`;
- a session epoch snapshots one token. The same Authorization header is used for
  `OpenModelSession` and its model RPCs;
- token refresh is **not** automatic. Explicit `ModelSessionClient.reconnect()` is
  required, which invalidates old model references before acquiring a new token.

This prevents a session that owns model references from silently changing its
principal midway through its lifetime.

## Endpoint contract

The v1 browser connector accepts only a root-relative same-origin gRPC-Web base
path such as `/model-api`. Absolute and protocol-relative URLs, query strings,
fragments and malformed paths reject locally. This is a deliberate token-leakage
boundary; cross-origin bearer forwarding requires a separately reviewed edge
policy.

The connector uses the existing pinned `@connectrpc/connect-web` transport and
canonical generated `ModelConfigurationService` / `ModelSessionService` types.
There is no new dependency and no duplicate wire schema.

## Hosted product ownership adapter

[`hostedWebModelOwner.ts`](hostedWebModelOwner.ts) is the product-facing adapter.
The trusted hosting shell may install one runtime capability on `window`:

```ts
window.mpmcHostedWebModelOwnership = {
  convention: 'MPMC/model/hosted-web-ownership/v1',
  baseUrl: '/model-api',
  identity: authoritativeIdentityProvider,
};
```

The capability contains a provider **function**, never an access-token value. Its
shape is exact and versioned. Missing, malformed or unknown-version capabilities
fail closed: the normal hosted configured-PT application remains active and no
model session is opened.

When the capability is valid, the frontend renders the existing shared Classic PR
product shell and `ExpertPr76Workspace`. The first model `create` lazily opens the
authenticated `ModelSessionClient`, then delegates model ownership to the existing
`bindExpertModelOwner`. No handle/session DTO, EOS logic, parameter conversion or
result interpretation is reimplemented. Page teardown disposes the session; an
aborted opening closes it rather than dispatching a late Create. There is no
automatic Create retry, reconnect or principal change.

The runtime capability must be supplied by trusted host integration code. Do not
put a bearer token in `VITE_*`, generated HTML, DOM attributes, browser storage or
a static runtime-config object.

## Production trusted-edge contract

The product adapter is usable in production only when the deployment provides an
authoritative trusted edge. The repository production renderer remains **PT-only
by default**. Model-session/configuration exposure is an explicit opt-in through
`model_authz_loopback_port` (CLI: `--model-authz-loopback-port`); without that
value no `/model-api/` routes or external-authorization filter are emitted.

When model routing is enabled, the rendered edge exposes only the versioned
`ModelSessionService` and `ModelConfigurationService` prefixes required by this
client. Those routes pass through a fail-closed external authorization service on
a loopback-only port. The edge sends that service only the bearer/session inputs
needed for authorization and accepts one internal
`x-mpmc-authenticated-principal` value back for the native upstream. The legacy
`PtFlashService` route explicitly bypasses this model authorization filter and
keeps its existing contract.

The browser is **not** authoritative for that principal header. Production native
model sessions are accepted only over the existing verified edge-client mTLS
channel. Once that peer is authenticated, the native host accepts at most one
bounded, character-validated `x-mpmc-authenticated-principal` value and binds the
model session to `web:<principal>`. Duplicate or invalid values fail closed. A
browser-supplied spoofed value therefore cannot establish session ownership; the
trusted edge authorization result is the value that reaches the native model
service.

A production deployment that enables this path must therefore provide the
external identity policy that:

1. validates the authoritative Web bearer at the trusted edge;
2. maps it to a stable end-user principal;
3. returns that principal to Envoy in the internal header contract;
4. keeps the native listener private and mutually authenticated to the dedicated
   edge-client CA;
5. starts `mpmc_pt_service_host` with model sessions enabled; and
6. publishes `mpmcHostedWebModelOwnership` only from trusted hosting integration
   after those conditions are satisfied.

No desktop launch token, model session ID or model handle may substitute for this
mapping. The repository intentionally does **not** invent a JWT issuer, JWKS URL,
account database, role model or third-party identity provider.

The versioned production bundle carries the same fail-closed behavior: omitting
`network.model_authz_loopback_port` preserves the historical PT-only edge/host;
providing it renders the model routes, adds `--enable-model-sessions` to the host
argv and records the non-secret hosted-model state in deployment metadata.

## Verification boundary

Focused connector tests cover same-origin endpoint enforcement, one-token-per-session
snapshot behavior, explicit token rotation on the next connector invocation,
cancellation, malformed/expired credentials and provider-error redaction. Hosted
owner tests additionally require exact runtime-capability shape, session opening
before the first Create, no token-valued shortcut field, and cancellation of a
still-opening session before Create dispatch. The existing `ModelSessionClient`
suite remains the authority for session metadata, reconnect, stale-reference and
mutation-failure behavior.

The `Hosted Web model product` workflow additionally builds the production native
host, validates the model-enabled Envoy configuration and drives a real Chrome /
ChromeDriver product regression through a test hosting shell, trusted edge and
native host. Its authorization sidecar recognizes only repository test tokens and
exists solely to exercise the contract: the regression deliberately sends a
spoofed internal-principal header, verifies a second principal cannot use the
first principal's session, then completes shared Classic PR apply/solve/release.
That fixture is **not** a production identity provider and carries no account or
physical-data meaning.
