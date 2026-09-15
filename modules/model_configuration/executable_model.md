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

## Result and request-error contract

The coarse `solve(request)` passes P [Pa], T [K] and mole feed in snapshot
component order to the existing backend once. It does not rebuild parameters,
change settings, inject hints, repair input or cache a previous result. Native
acceptance/rejection and numerical outcomes remain authoritative.

After a native input exception, the model can attach the `Pr76SolveRequestError`
location interface while preserving the original `std::invalid_argument`,
`std::domain_error` or `std::length_error` base and `what()` text. Catch those
standard categories as before, or catch the location interface first to read
`field()`. Located errors have a derived dynamic type; code must not require exact
`typeid` equality with the standard base.

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

## Public per-solve initialization and continuation hints

`pt_solve_hints.hpp` adds the transport-neutral, stdlib-only
`pt-solve-hints/v1` DTO. Hints belong to **one solve call**; they are not stored
in the immutable model/settings snapshot and never constitute phase-count
or stability evidence.

```cpp
auto hints = mc::make_pt_solve_hints_v1();
hints.initial_stability_starts = initial_compositions;
hints.final_two_phase_stability_starts = final_review_compositions;

mc::PtThreePhaseContinuationHint three;
three.compositions = {phase0, phase1, phase2};
three.phase_fraction_seed = {beta1, beta2};
hints.three_phase_continuation_starts.push_back(three);

auto result = model->solve({pressure_pa, temperature_k, mole_feed}, hints);
```

`solve(request, hints)` copies the validated model's prepared native options into
a temporary backend-options snapshot and maps only the supplied starts into the
existing `initial_starts`, `final_starts` and `three_phase_starts` fields. The
model's persistent coarse backend and settings snapshot remain unchanged, so a
later `solve(request)` is still a cold/native-default solve.

Version and host storage ceilings are checked before copying untrusted hint
storage. State-dependent scientific validation remains in the existing solver:
all stability starts must match the current feed dimension, be normalized under
the native tolerance and have exactly the current feed's active support. Final
two-phase review still reserves its two freshly solved phase starts. Three-phase
continuation seeds additionally require three valid compositions and a feasible
phase-fraction simplex seed, while native count/entry/attempt budgets remain
independent. A hint can initialize a search; it cannot force an accepted phase
set or bypass final stability/equilibrium gates.

When the native solver rejects a hint, the executable model can attach a precise
state-related location such as
`hints.initial_stability_starts[0][2]`,
`hints.final_two_phase_stability_starts[0]`,
`hints.three_phase_continuation_starts[0].compositions[1][2]`, or
`hints.three_phase_continuation_starts[0].phase_fraction_seed`, while preserving
the native standard exception category and message. Unsupported DTO versions and
pre-copy host ceilings use `ModelConfigurationError` with their public field.

`make_pr76_continuation_hints(previous_result)` converts only an **accepted** PR76
publication into numerical starts. It copies accepted phase compositions into
initial/final stability starts and, for an accepted three-phase result, also
copies the three compositions plus independent phase fractions into one
three-phase continuation seed. An unresolved/non-accepted result deliberately
produces an empty v1 hint set so stale state is not carried across a failed point.
The next solve is always fresh and must re-establish the phase set.

This increment is direct-C++ only. The model-neutral `PtFlashBackend` virtual
interface, bounded registry/service wire request and desktop/UI request types do
not yet serialize these optional hints. Their current coarse solve behavior is
unchanged; UI configuration remains deferred.

## Audit and verification scope

Baseline: PR #82 through `6db1ae25e4922f38fbac707bd6b2046a6748e309`.
Read the root `AGENTS.md`, parameter/settings adapters, phase/evaluator ownership,
backend result/publication/transition contracts, existing process ownership graph
and focused workflow before editing. Leaving separate preparation APIs alone
would not provide safe execution lifetime; extending the process host in the same
increment would couple native semantics to transport composition. A small optional
owner plus per-call DTO reuses all existing scientific algorithms without a
premature UI or workspace refactor.

The independent `tests/model_configuration/executable_model` project has 17 cases:
full-envelope direct-C++ parity for single/binary/ternary cold search; changed-kij
A/B isolation; independent settings/root-budget exhaustion; draft/owner/result
lifetime; construction failure and host limits; input exceptions and declared
bounds; precise request locations, normalization boundary and coarse internal-
failure fallback; public seeded-hint direct-C++ parity; accepted-result
continuation parity plus unresolved-state clearing; nested hint/state/quota
validation; deterministic admission/unwind; and parallel A/B solves. The public
hint header is compiled first in a separate translation unit to enforce stdlib-
only/self-contained use, with ownership/const-interface static assertions.

Expected results come from separately constructed native PR76 models/options,
not the public mapper. Same-platform exact field equality is required. Data reuse
is limited to the already attributed Hua nitrogen/ethane numerical fixture and
the explicitly synthetic max3 structural fixture; this is adapter/ownership
regression, not independent physical validation or a new parameter dataset.

The focused GitHub-hosted workflow runs these cases on GCC Debug + ASan/UBSan,
Clang Release and MSVC Release, retaining parameter/settings, registry and GCC
thermodynamic-contract checks. Commit `6db1ae25` passed all three compiler jobs;
the executable suite was 17/17 and the registry suite remained 17/17 under the
GCC sanitizer job. No local compilation is used as official verification.
Existing EOS/flash algorithms and UI call sites were not changed.

The separate [bounded registry](registry.md) owns coarse-dispatch models, issues
opaque handles and retains admitted solves through release/close. It accounts for
model slots until destruction. Standalone factory users still own lifetime
themselves. The versioned model service and desktop typed client continue their
existing safe request/error transport; public hint transport and UI configuration
remain deferred.
