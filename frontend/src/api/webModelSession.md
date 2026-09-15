# Authenticated Web model session connector

This module is the first Web-identity increment for the versioned model service. It
connects an **authoritative browser identity provider** to the existing
`ModelSessionClient`; it does not expose model routes publicly and it does not make
Electron credentials usable from Web code.

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

## What this does not provide

This increment does **not** make the current native model service a Web-user
service by itself. The existing production Envoy route intentionally exposes only
the legacy PT service, and the native host currently authenticates its private
bearer or mTLS edge identity. An Envoy certificate is an edge identity, not an end
user.

Before the browser connector can be enabled in the product, a later edge increment
must:

1. validate the authoritative Web access token at the trusted edge;
2. map it to a stable end-user principal without accepting browser-supplied subject
   headers;
3. propagate that principal to the native host through an authenticated trusted
   channel;
4. expose only the versioned model-session/configuration routes required by this
   client;
5. preserve session-owner isolation, CORS/CSRF policy, quotas, deadlines and error
   redaction already frozen by the model-service contract.

No desktop launch token, model session ID or model handle may substitute for that
identity mapping.

## Verification boundary

Focused tests cover same-origin endpoint enforcement, one-token-per-session
snapshot behavior, explicit token rotation on the next connector invocation,
cancellation, malformed/expired credentials and provider-error redaction. The
existing `ModelSessionClient` suite remains the authority for session metadata,
reconnect, stale-reference and mutation-failure behavior.
