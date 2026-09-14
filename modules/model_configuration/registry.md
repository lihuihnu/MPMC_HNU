# Bounded PR76 model registry

`Pr76ModelRegistry` owns immutable executable PR76 models for one native session.
It adds no global cache, TTL, transport endpoint or background worker. Link the
separate optional compiled target (the previous preparation/factory targets keep
their existing dependencies):

```cmake
add_subdirectory(modules/model_configuration/pr76_registry)
target_link_libraries(my_consumer PRIVATE mpmc::pr76_model_registry)
```

```cpp
namespace mc = mpmc::model_configuration;
mc::Pr76ModelRegistryLimits limits;
limits.max_models = 16;
mc::Pr76ModelRegistry registry(limits);
auto handle = registry.create(definition,
    mc::resolve_pt_solver_preset(mc::mpmc_balanced_default_v1));
auto snapshot = registry.describe(handle);
auto result = registry.solve(handle, {pressure_pa, temperature_k, mole_feed});
registry.release(handle);
```

The caller supplies a complete definition and explicit settings. The registry
uses the existing factory/validators; construction never substitutes data,
changes tolerances or creates a mutable configuration alias. `describe` returns
owning copies of definition, resolved settings, host limits and backend capability.
Changing those copies cannot change the model. Model changes require creation of
a new handle; existing configured PR/SW/CPA backends and PtService v1 are untouched.

## Capacity and ownership

The default `max_models` is 64 and must be positive. It bounds **all model slots**:
creation reservations, published models and released models still retained by
admitted solves or in-progress description copies. A slot is reserved before
entropy generation, validation, parameter matrices or model construction. Every
failure rolls the reservation back; publication happens only after full creation.
Close racing with creation prevents publication and reclaims the reservation.

The inherited parameter and solver limits remain host-owned. Component/pair,
matrix, identifier, provenance/text and solver ceilings are enforced by the
existing adapters, including the factory's cross-layer component check. The
registry capacity is a model-count bound, not a measured byte quota. Returned
results/descriptions belong to callers. A future service must separately bound
serialized bytes **before parsing**, concurrent requests and caller-owned copies.

`status()` reports `registered_models`, `resident_models` (including reservations)
and `closed`. After release, the first may decrease before the second. Capacity
is reclaimed only after the retained model is destroyed, not when its handle is
removed. Registry state and capacity accounting have separate owners to avoid a
map/entry reference cycle. No released-handle tombstone list accumulates.

## Admission and release

| Operation | Lifetime / concurrency contract |
| --- | --- |
| `create` | Reserve under a short mutex; construct outside it; publish under it. |
| `describe` | Pin the entry under the mutex, then copy immutable snapshots outside it. |
| `acquire_solve` | Validate handle and acquire an exclusive, move-only solve lease under the mutex. |
| `solve` | Convenience acquire-and-execute; uses the same lease path. No registry lock during numerical work. |
| `release` | Remove the handle under the mutex. Later lookup/admission fails immediately. Existing leases remain valid. |
| `close` | Idempotently invalidate all handles, reject future operations and pending publication. Existing leases remain valid. |

The acquisition/removal mutex defines which operation wins a release/admission
race. An acquired lease is an **already admitted operation**, even before the
numerical kernel starts. It can finish after release, close or registry destruction.
It cannot admit a second operation. Call `std::move(lease).solve(request)` once;
empty, moved-from or consumed leases throw `invalid_lease`. Moving a lease transfers
ownership; destruction or replacement of an unused lease cancels its admission.
Return and exception both release admission and ownership after numerical work ends.

Only one lease may occupy a model at a time; another acquisition throws the
existing `Pr76ModelBusyError`. Immutable inspection can proceed concurrently.
Different models own independent workspaces. Neither registry mutex nor a global
solver lock serializes different models. Every successful solve returns the
existing complete native `PtFlashBackendResult` unchanged; native exceptions and
indeterminate outcomes propagate unchanged.

Keep the registry C++ object alive while its public member functions are running;
use `close()` for concurrent shutdown. Detached leases can outlive the registry.
Registry destruction retires its registered entries without waiting for leases.
There is no cancellation of a numerical operation and no API session authorization
layer in this increment.

## Handles and errors

Handles are bearer capabilities: opaque to callers, neither model data nor memory
addresses. Each creation uses 256 fresh random bits plus a checked monotonic
sequence. The sequence forbids reuse within a registry even if entropy repeats;
exhaustion rejects creation instead of wrapping. New entropy prevents prediction
of another handle from an observed one. Cross-registry/session separation relies
on the independent secure random draws. Handle internals are not a wire contract.

The built-in Linux adapter requests 32 bytes from
[`getrandom`](https://man7.org/linux/man-pages/man2/getrandom.2.html) with
`GRND_NONBLOCK`; incomplete/unavailable entropy fails creation. Windows uses
[`BCryptGenRandom`](https://learn.microsoft.com/en-us/windows/win32/api/bcrypt/nf-bcrypt-bcryptgenrandom)
with the system-preferred RNG. Neither adapter falls back to a counter-only token,
`std::random_device`, a clock or a user-space seeded PRNG. Other platforms must
supply a host-owned `ModelHandleEntropySource`; the default fails closed there.
Native platform entropy support is verified on Linux/Windows only; this does not
claim Android/macOS product integration. A supplied callback must provide fresh
cryptographic entropy, be thread-safe and propagate failure. Deterministic test
sources are explicitly unsuitable for production and cannot come from a draft.

Malformed handles are rejected at fixed length/character bounds before lookup;
valid-shaped unknown, foreign or released handles yield `model_not_found`, as
does duplicate release. A closed registry rejects operations with `registry_closed`.
Quota, empty-lease and entropy failures have distinct `ModelRegistryErrorCode`s.
Existing model configuration and scientific exceptions retain their original types.
Errors do not echo handles; hosts must not log/display complete bearer tokens.

## Audit and verification

Pre-write baseline: PR #82 `4a921f453d71b9bc99b0e83744074b49ea870196`.
Reviewed root AGENTS, uploaded Gate registry/resource/concurrency requirements,
factory/evaluator ownership, immutable snapshots and PtService's borrowed static
backend registry. Extending PtService would change v1 identity/lifetime semantics;
counting only map entries would permit release to bypass the resident-model cap.
A separate compiled owner with reservation RAII and one-shot leases is the minimal
extension. OS entropy linkage stays outside the header-only factory module.
No EOS, flash formula, tolerance, parameter fixture or product call site changes.

The independent registry project contains 17 cases plus public-header ownership
assertions. It reuses the existing complete-result comparator and attributed or
explicitly synthetic native fixtures. Coverage includes creation/description/solve
parity; malformed/foreign/stale handles and repeated release; capacity and host
quotas; configuration/entropy failure rollback; concurrent reservation and close;
release or registry destruction while an admitted solve is paused; close with a
lease; single-use/move/cancel/busy behavior; solve exception cleanup; parallel A/B;
release/acquire races; simultaneous creators competing for capacity.

Lifecycle tests synchronize at the admission boundary and then execute real PR76
**after** invalidation, covering the ownership gap without sleeps or claims of
pausing inside a scientific kernel. Repeated deterministic test entropy exercises
non-reuse separately from production OS entropy tests. These are ownership and
adapter tests, not physical validation, RNG certification or a TSan result.
The focused GitHub-hosted GCC ASan/UBSan, Clang and MSVC jobs instrument/build the
new library and tests, retaining all earlier model-configuration regressions.

The [additive configuration service](../model_configuration_grpc/README.md) now
provides bounded protobuf/gRPC mapping and create/describe/solve/release tests.
Native host session routing is opt-in. Web identity routing, shared clients/UI, public hints, one-sided applicability
and product vertical slices remain open Gate items. Existing `configured_backend_id`
wire fields retain their old meanings.
