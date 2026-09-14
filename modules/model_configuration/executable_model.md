# Executable PR76 model factory

`pr76_executable_model.hpp` joins the validated parameter and settings snapshots
with the existing PR76 phase/evaluator/max3 backend. Link the optional
`mpmc::pr76_solver_configuration` target; the parameter-only target still has no
flash dependency. No new runtime, transport or third-party dependency is added.

```cpp
#include <mpmc/model_configuration/pr76_executable_model.hpp>

namespace mc = mpmc::model_configuration;
// definition is a complete, owning ThermodynamicModelDefinition draft.
auto settings = mc::resolve_pt_solver_preset(mc::mpmc_balanced_default_v1);
auto model = mc::make_pr76_executable_model(definition, settings);
auto result = model->solve({pressure_pa, temperature_k, mole_feed});
```

Creation requires both inputs explicitly. Settings version/preset/domain and host
ceilings are validated before preparing parameters. Component count must fit both
parameter and solver ceilings; reject a mismatch before allocating the parameter
matrix. The existing parameter validator remains the source of physical-domain,
provenance, units, component-order and pair-completeness rules. Construction
failures publish no model. Synthetic fixture policy remains a host-only opt-in.

## Ownership and concurrency

`make_pr76_executable_model` returns `unique_ptr<Pr76ExecutableModel>`. The object
also implements the existing `flash::PtFlashBackend` coarse dispatch interface.
Its constructor is available for native stack ownership, but the object cannot
be copied or moved: its backend borrows its evaluator by reference. Moving the
owning pointer preserves that address. Reverse destruction releases the backend
before the evaluator. The evaluator owns its own phase-model copy and workspace.

`parameter_snapshot()`, `solver_configuration()` and `parameter_limits()` expose
const references. Both public definitions and prepared native options are owned;
editing or destroying either input draft after creation cannot alter execution.
Do not mutate drafts during construction. Read-only inspection of snapshots and
capability is supported while solving. Const references remain valid only while
the owning model lives.

Each model admits one active `solve`; overlapping/reentrant calls throw
`Pr76ModelBusyError` immediately. An atomic RAII guard releases admission on both
normal return and exceptions. Rejection is a scheduling error, not an EOS result.
There is no waiting queue, cancellation or asynchronous lifetime management.
Callers must keep the model alive until all calls, including rejected calls, end.
Independent models can solve concurrently because they share no mutable workspace.

## Result contract

`solve` passes P [Pa], T [K] and mole feed in snapshot component order to the
existing backend once. It does not rebuild parameters, change settings, inject
hints, repair input or cache a previous result. Native acceptance/rejection and
numerical outcomes remain authoritative. After a native input exception, the
model can attach the `Pr76SolveRequestError` location interface while preserving
the original `std::invalid_argument`, `std::domain_error` or `std::length_error`
base and `what()` text. Catch those standard categories as before, or catch the
location interface first to read `field()`. Located errors have a derived dynamic
type; code must not require exact `typeid` equality with the standard base.

Diagnosis executes only on a rejected call and reuses native PT/composition
validation, including its compensated sum and normalization tolerance. It emits
`pressure_pa`, `temperature_k`, `feed` (length/aggregate normalization), or
`feed[index]` (a nonfinite or out-of-range mole fraction). Indices refer to the
submitted snapshot component order. Declared dataset intervals use the prepared
native snapshot and inclusive endpoints. Absent ranges remain unknown.
For multiple invalid fields the diagnostic priority is feed dimension, PT,
individual feed values, aggregate normalization, then declared P/T bounds.
This is an advisory invalid-request location, not a claim about an internal
throw site or an aggregate error report. No exception-message parsing is used.
A valid request that hits an internal search/resource/property failure keeps
its original exception without a fabricated request field. Result envelopes,
per-model admission release and model isolation remain unchanged.

The returned **entire existing `PtFlashBackendResult`** owns capability/model and
dataset identity, solution state/status/feed, candidate phases (also when not
accepted), mole fractions, compositions, Z, ln(phi), branch/smoothness flags,
diagnostic, stability limitations, phase-transition evidence, provider convention,
phase metadata and morphology flag. It survives the model's destruction.
`accepted_phase_set()` remains the existing publication guard. This does not add
new residual/iteration arrays or serialize the settings into the result: the full
resolved settings and their version/preset identity are queried on the model.
No global stability or PR phase morphology claim is added.

Public hints are still unsupported. The existing max3 backend may fail to close a
three-phase topology without starts; the factory must preserve that outcome.
Seeded accepted-three-phase option parity remains covered in the previous settings
suite, where identical native hints are explicitly supplied to both paths.

## Audit and verification scope

Baseline: PR #82 at `174e004e366c4b31f2dc476fb43feeaafabd365f`.
Read the root `AGENTS.md`, parameter/settings adapters, phase/evaluator ownership,
backend result/publication/transition contracts, existing process ownership graph
and focused workflow before editing. Leaving separate preparation APIs alone
would not provide safe execution lifetime; extending the process host would couple
native callers to transport composition. A small optional owner reuses all existing
scientific algorithms and avoids a premature registry or workspace refactor.

The independent `tests/model_configuration/executable_model` project has 14 cases:
full-envelope direct-C++ parity for single/binary/ternary cold search; changed-kij
A/B isolation; independent settings/root-budget exhaustion; draft/owner/result
lifetime; construction failure and host limits; input exceptions and declared
bounds; precise request locations, normalization boundary and coarse internal-failure fallback; deterministic admission/unwind; parallel A/B solves. The public header
is compiled separately with ownership/const-interface static assertions.

Expected results come from separately constructed native PR76 models/options,
not the public mapper. Same-platform exact field equality is required. Data reuse
is limited to the already attributed Hua nitrogen/ethane numerical fixture and
the explicitly synthetic max3 structural fixture; this is adapter/ownership
regression, not independent physical validation or a new parameter dataset.

The focused GitHub-hosted workflow adds these cases to GCC ASan/UBSan, Clang and
MSVC, retaining parameter/settings and the existing GCC thermodynamic-contract
checks. No local compilation is used as official verification. Existing algorithms,
curated product backends and wire/UI call sites do not change, so unrelated
product suites are not selected. Per-model admission uses standard acquire/release
atomics; thread tests check actual independent execution, without claiming a TSan
run or a general thread-local-workspace redesign.

The separate [bounded registry](registry.md) now owns these models, issues opaque
handles and retains admitted solves through release/close. It accounts for model
slots until destruction. Standalone factory users still own lifetime themselves.
The versioned model service and desktop typed client now transport safe error
locations; UI configuration forms remain deferred.
