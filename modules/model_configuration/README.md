# Public thermodynamic model configuration — executable PR76 models

The first five increments of the **unified model configuration + PR76 Expert
parameter interface Gate** prepare immutable parameter and numerical-settings
snapshots and combine them in an [executable PR76 model](executable_model.md).
The [versioned solver contract](solver_settings.md) resolves an explicit preset
and maps complete public settings to the existing PR76 solver options.
A [bounded native registry](registry.md) adds opaque model handles and release.
An [optional versioned gRPC service](../model_configuration_grpc/README.md) now
provides native RPC creation/description/solve/release. Web/Electron/Android
client/UI integration remains later work; native host sessions are now opt-in.
Parameter preparation and executable creation never silently select settings.

## Current boundary

```text
ThermodynamicModelDefinition (owning, unvalidated draft)
  -> version/family/shape/text checks
  -> PR76 field mapping
  -> existing OrderedComponents / PrParameterSet::create
  -> Pr76ModelParameters (owning, immutable parameter snapshot)
```

`model_definition.hpp` includes only the C++20 standard library. The separate
`pr76_parameters.hpp` adapter depends on the existing thermodynamics data
contract, with no EOS evaluation, flash, runtime, Protobuf, gRPC, or UI dependency.
The optional header-only target is `mpmc::model_configuration`:

```cmake
add_subdirectory(modules/model_configuration)
target_link_libraries(my_consumer PRIVATE mpmc::model_configuration)
```

Create with `Pr76ModelParameters::create(definition, host_limits)`. Validation
failures throw `ModelConfigurationError` with a code and field context. Allocation
failures and programming errors propagate. The returned object keeps both the
original public definition and the validated PR parameter set; no mutable getters
or assignment are provided. Editing or destroying the draft cannot change an
existing snapshot. Concurrent const reads are supported; callers must not mutate
the input draft during creation.

## Parameter contract

| Field | Semantics |
| --- | --- |
| `version` | Required exact `thermodynamic-model/parameter-definition/v1`; versions are not inferred. This versions the parameter layer, not the eventual full creation RPC. |
| `family` | `peng_robinson_1976` accepts the PR parameter variant. `soreide_whitson_1992` and `cubic_plus_association` are reserved and reject custom configuration. |
| `components` | All selected components, in canonical model order; stable, unique, visible ASCII IDs, display names, pure/pseudo kind, definition provenance, optional kg/mol molar mass. |
| `parameters.pure` | ID-keyed mandatory Tc in K, Pc in Pa, and dimensionless omega, each with provenance and original-unit/conversion history. No implicit missing values. |
| `parameters.binary` | Explicit constant, symmetric dimensionless `kij` for every unordered pair. Missing, duplicate/reversed, self and unknown pairs fail; diagonal zero is structural. |
| `dataset_id`, `revision`, `provenance` | Explicit dataset and model-definition identity; provenance is retained rather than upgraded or invented. |
| `applicability` | Optional independent T/P endpoints, closed by default with explicit open-endpoint flags, plus a required declaration source. Missing endpoints mean **unknown**, never unlimited validity; see [endpoint semantics](applicability.md). |

The existing thermodynamics contract remains the authority for physical domains,
component identity, provenance completeness, parameter completeness and range
validity. In particular, Tc/Pc must be finite positive values; omega and `kij`
have no invented universal range. Public values are already canonical SI; the
adapter records original-unit history without performing a second conversion.

Each provenance record preserves kind, reference, revision, locator, note,
acquisition and usage terms. `user_supplied` is supported and still needs an
explicit user record identity and revision; it does not require a fabricated
literature citation. Ordinary creation rejects `synthetic_test`; a host-owned
test policy can explicitly admit artificial software fixtures.

**Applicability contract:** the public definition supports lower-only, upper-only,
complete and absent bounds independently for temperature and pressure. A missing
endpoint remains unknown rather than becoming an implicit infinity. Declared
endpoints are closed unless their corresponding `*_exclusive` flag is true;
exclusive markers without endpoints and empty/inverted intervals are rejected.
The existing thermodynamics `Applicability` still stores complete closed intervals,
so the adapter maps only complete numeric envelopes and never fabricates a missing
endpoint. `Pr76ExecutableModel` supplies the public-only enforcement gap for
one-sided violations and exact open-endpoint equality, while strict violations of
complete intervals continue through the unchanged native PR76 property path.
No PR76 formula, root, stability, split, three-phase tolerance or acceptance rule
is changed. See [the full endpoint contract](applicability.md).

## Resource policy

Host limits default to 256 components/pure records, 32,640 pair records, 65,536
dense matrix entries, 128 bytes per identifier, 512 per display name, 4,096 per
provenance/history string and 1 MiB total text. All string occurrences count,
including repeated sources. Bounds run before mapping/copying; matrix products
are guarded by division. These are logical shape/text limits, not an exact heap
or serialized-byte quota. Future transports must reject oversized serialized
input **before parsing**. The optional registry bounds published, creating and
released-but-retained model slots together; see [registry policy](registry.md).
Limits and synthetic-data policy are host arguments, never fields of the draft.

## Audit and staged integration

Audit baseline: main `3a1c3c1c09e2176341fcf55a84525503de50c5d9` (PR #81).
The existing parameter, PR pure/mixing/root/phase, stability/split/max3/publication,
runtime, gRPC, process host, React clients/components, Electron and Android paths
were inspected before this increment. Key findings and decisions:

- Existing `PrParameterSet::create` already implements the required missing-data,
  symmetric-pair, provenance and ordering rules. This adapter reuses it and the
  existing provenance/text helpers instead of duplicating physical validation.
- A new optional module is the smallest shared entry for native and transport
  consumers. Extending only the process factory would couple Android to its
  gRPC composition root (Android currently uses a private ownership shim).
- Process factories own a stable model/evaluator/backend object graph;
  `PtService` borrows registered backend pointers. The dynamic factory must own
  that graph for the entire lifetime of every admitted solve.
- PR evaluators own mutable workspaces; the gRPC adapter defaults to one admitted
  solve. Immutable parameter snapshots do not make a shared evaluator concurrent.
  Preserve serialized calls until a separate workspace design is validated.
- `mpmc.runtime.v1` freezes existing field meanings and currently exposes only
  discovery and solve. Add a separately versioned configuration service later;
  do not reinterpret `configured_backend_id` as a model handle.
- Electron lives in `frontend/desktop`, not a root `desktop` directory. Its
  gateway and Android's Capacitor/JNI bridge both need later transport increments.
  No existing configured PR76/SW92/CPA call site imports this new module yet.

The solver-settings increment adds `PtSolverSettings`, `mpmc-balanced-default/v1`,
independent initial/final stability controls, host ceilings and direct-C++ option
and solve parity. Its optional `mpmc::pr76_solver_configuration` target imports
flash; the original parameter target still does not.

The executable increment adds `make_pr76_executable_model`: one stable-address
owner retains both snapshots and an evaluator/backend, rejects overlapping solves
on that model, and returns the complete native result by value. Independent
models have independent workspaces. It uses the same optional solver target.

The optional `mpmc::pr76_model_registry` target adds capacity reservations, opaque
handles, release/close and one-shot solve leases. Releasing a handle prevents new
access while preserving admitted work and its capacity charge until destruction.

The separately versioned configuration service adds bounded gRPC mapping and
real native roundtrip tests. The opt-in native host now connects authenticated
session routing; shared clients precede Expert UI, Electron/Android integration and final
regression. The full Gate additionally requires end-to-end dynamic-vs-direct flash parity,
session lifecycle and all three existing configured-backend product regressions.

SW/CPA custom editing, PR morphology classification, wide-range PR validation,
workspace performance refactoring and persistent fluid-project storage remain
outside this increment. No global stability or arbitrary P/T validity is claimed.

## Verification

The parameter test project has the original 11 cases plus an independent
applicability-endpoint contract case covering one-sided presence, open/closed
boundaries, singleton/empty intervals, exact field errors and conservative native
mapping. Values are explicitly synthetic software fixtures, not physical validation.

The focused GitHub-hosted workflow runs GCC Debug + ASan/UBSan, Clang Release and
MSVC Release; the GCC job also runs the unchanged thermodynamics contract suite.
The registry project adds 17 capacity, handle and concurrent-lifetime cases.
The separate executable test project has the original 17 cases plus an independent
applicability-endpoint solve regression covering one-sided and exact open-boundary
rejection for both cold and hinted paths. The solver test project adds 14 cases,
including exhaustive field presence, distinct custom option mapping, preset identity,
quotas, immutable settings and actual one-/three-phase and exhausted-budget PR76
parity. The versioned gRPC project separately exercises applicability endpoint
presence/snapshot roundtrip and structured solve error locations. Its commands include:

```bash
cmake -S tests/model_configuration/parameters -B build/model-parameters \
  -DCMAKE_BUILD_TYPE=Debug -DMPMC_MODEL_ENABLE_SANITIZERS=ON
cmake --build build/model-parameters --parallel 2
ctest --test-dir build/model-parameters -R '^model[.]parameters[.]' --verbose --no-tests=error
```

No existing scientific algorithm, threshold, fixture, or `mpmc.runtime.v1`
`PtFlashService` field is changed. The versioned model-configuration service gains
only additive applicability endpoint flags with closed defaults. UI code is not
changed by this increment; product workflows may still run when their dependency
filters include the shared configuration/proto paths.
