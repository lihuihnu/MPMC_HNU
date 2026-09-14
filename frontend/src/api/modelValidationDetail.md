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
Native wire decoding now reports enclosing snake_case request paths, for example
`definition.pr76.pure[1].critical_temperature_k.provenance.kind` or
`definition.components[1].molar_mass_kg_per_mol.value`. These indices are zero-based
positions in the submitted records; unvalidated component IDs are never inserted.
Later semantic validation retains its existing keyed paths and may still report
only `request`, `version` or an options block. No finer location is inferred, and
no field is fabricated when none is provided. Older native services can still
return coarse aliases such as `scalar.value`. Parameter keys remain keys even when they contain dots
or happen to look numeric. Native validation remains authoritative.

For solve validation, `request.rejected` can identify `pressure_pa`,
`temperature_k`, `feed` (dimension or aggregate normalization), or `feed[1]`
(an invalid item in submitted component order). Missing P/T uses the same fields
with `wire.missing_field`. Only canonical nonnegative decimal feed indices are
accepted: keys, pairs, signs, fractions and leading zeroes are not feed indices.
Native declared-interval violations also identify P/T. This metadata is added
only after native rejection; internal failures with valid input remain coarse.
A rejected solve does not invalidate ownership or require reconnect. Old v1
clients keep their status-only shape, and older services can still return `PT`
or `request`. Protocol/detail versions are unchanged.

| Boundary/version | Behavior |
| --- | --- |
| Native service v1 | Same rich-error envelope/status/domain codes; decode fields now retain enclosing records. Numerical validation unchanged. |
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
result equality and independent window reclamation. `MODEL_VALIDATION_DETAIL_OK`,
`MODEL_NESTED_VALIDATION_PATH_OK` and `MODEL_SOLVE_VALIDATION_PATH_OK` are required by the smoke launcher alongside
all existing completion markers. Nested missing scalar values and provenance
kinds are checked through the real native host, with explicit reconnect/recreate
and recovered full-result equality. Solve coverage includes missing, zero,
negative and nonfinite P/T, feed shape/item/sum failures, Protobuf JSON nonfinite
roundtrips and full-result recovery after every failure on the same reference.
UI field rendering and configuration forms remain separate work.
