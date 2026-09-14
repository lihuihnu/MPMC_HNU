# Versioned model validation errors

The native model service already emits `ModelServiceError` with its existing
`mpmc.model_configuration.v1/model-service/v1` contract. The shared session client
now preserves eligible details in `ModelClientError.validation`; the typed
renderer exposes the same immutable DTO in `RendererModelError.validation`:

```ts
{
  version: 'MPMC/model/validation-detail/v1',
  code: 'configuration.missing_parameter',
  field: 'parameters.pure[nitrogen].critical_temperature_k'
}
```

`code` is a native domain code, separate from the RPC status/category and fixed
client reason. `field` is optional. It preserves native path spelling, including
component key, pair and index selectors; it is advisory, not a JSON Pointer.
The existing service sometimes reports only `request`, `scalar.value`, `version`
or an options block. No finer location is inferred, and no field is fabricated
when none is provided. Parameter keys remain keys even when they contain dots
or happen to look numeric. Native validation remains authoritative.

| Boundary/version | Behavior |
| --- | --- |
| Native service v1 | Existing rich error bytes and numerical validation unchanged. |
| Desktop bridge v1 | Existing `mpmcModelDesktop`, invoke/cancel channels and exact `{code, reason}` error shape unchanged. |
| Desktop bridge v2 | `mpmcModelDesktopV2`, `mpmc:model:invoke:v2` and `mpmc:model:cancel:v2`; errors may add `validation`. |
| Validation detail v1 | Independent version, allowlisted domain/status pair and optional safe native field path. |

Both desktop versions share the same owner, reference map, pending-call cap and
native session per window. Version switching does not grant extra capacity or
ownership. They enforce identical attached-main-frame authorization, navigation,
crash/close cleanup and disposal. A request must match the version of its channel.

The renderer factory prefers the v2 preload when present and otherwise supports
an old v1 preload. A selected bridge's reply must match its version. Failed calls
are never retried on v1 and no implicit model recreation is introduced. Old
renderers using v1 receive no additional error keys. New renderers with old
preloads receive the existing status/category without field details.

Only known `configuration.*` validation codes, `wire.missing_field`,
`wire.unknown_field`, `wire.feed_limit` and `request.rejected` are eligible, with
exact expected RPC statuses. Unknown versions/codes, status mismatches, duplicate,
malformed or oversized details produce status-only errors. Auth, session,
registry, transport and internal failures do not gain validation details.
The client accepts at most eight incoming binary details, at most 4096 bytes
each and 8192 bytes total, with exactly one matching ModelServiceError before
decoding. No raw metadata or native exception object is retained.

Paths are at most 256 ASCII characters with recognized schema/native names and
restricted component selectors. Unrecognized paths, control characters,
credential/handle/session-like text and current private header/session/handle
values are omitted while the known domain code remains. Some legal native
component IDs contain punctuation outside the safe selector subset; their field
is deliberately omitted, not rewritten. Main and renderer recheck the DTO and
make an immutable copy. Future incompatible optional details are ignored without
changing the authoritative RPC status. Neither exceptions, stack traces,
request bodies nor private routing values are included.

Tests cover bounded binary decoding, exact paths, privacy filtering, immutable
copies, protocol compatibility and shared version ownership. Official real
sandboxed renderer regressions retain v1 traffic and execute v2 native parameter,
settings, preset and coarse solve validation failures, explicit recovery, full
result equality and independent window reclamation. `MODEL_VALIDATION_DETAIL_OK`
is required by the smoke launcher in addition to all existing completion markers.
UI field rendering and configuration forms remain separate work.
