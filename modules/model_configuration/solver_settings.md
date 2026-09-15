# Versioned public PT solver settings

`PtSolverSettings` is a transport-neutral C++20 numerical configuration with 57
explicitly present number/boolean fields. The DTO/preset header
[`pt_solver_settings.hpp`](include/mpmc/model_configuration/pt_solver_settings.hpp)
uses only the standard library and public configuration types. Its version is
`pt-solver-settings/v1`. The PR76-specific adapter is separate:
[`pr76_solver_settings.hpp`](include/mpmc/model_configuration/pr76_solver_settings.hpp).

## Preset and custom snapshots

`resolve_pt_solver_preset("mpmc-balanced-default/v1")` returns **MPMC balanced
default v1**, including every actual resolved number. It uses frozen literal
values audited against PR #82 at `9d856b8acd3162d0b8a38ef2cd9b4f9d47e8e056`, not
the default constructors of current backend option types. Unknown IDs fail;
there is no `latest`, implicit preset or field-by-field fallback.

A preset-labeled snapshot must exactly match its declared preset. To edit it,
declare a complete custom configuration and clear the preset ID:

```cpp
#include <mpmc/model_configuration/pr76_solver_settings.hpp>

namespace config = mpmc::model_configuration;
auto settings = config::resolve_pt_solver_preset(config::mpmc_balanced_default_v1);
settings.kind = config::PtSolverSettingsKind::custom;
settings.preset_id.clear();
settings.final_two_phase_stability.max_property_evaluations = 50000;
const auto prepared = config::Pr76SolverConfiguration::create(settings);
// prepared.settings(): owning, immutable resolved values and identity.
// prepared.safety_limits(): the separately applied host policy.
// prepared.root_options()/backend_options(): native factory handoff.
```

Link the optional `mpmc::pr76_solver_configuration` target by adding
`modules/model_configuration/pr76_solver`. Existing parameter-only consumers
can continue to link `mpmc::model_configuration` without a flash dependency.

All omitted fields fail, including missing booleans and counters. Numeric zero
and boolean false are never interpreted as absence. Counts use explicit signed
64-bit input, checked for nonnegativity, native representability and host
ceilings before narrowing. The snapshot cannot be assigned or edited in place;
modifications require a new preparation. There is no mutable global preset.

## Frozen values and native mapping

All tolerances and step/progress coefficients below are dimensionless. Limits
count the operations specified by the existing algorithms, not wall-clock time.
Mapping names in this document are for C++ integration, not Standard/Expert UI
labels.

`eos_root.max_iterations = 2048` maps to the existing PR root iteration limit.

Each of `initial_stability`, `final_two_phase_stability` and
`final_three_phase_stability` has its **own independent** complete group:

| Public field | Default v1 | Native stability option |
| --- | ---: | --- |
| `tpd_tolerance` | 1e-10 | `tpd_tolerance` |
| `stationarity_tolerance` | 1e-8 | `stationarity_tolerance` |
| `line_search_armijo_coefficient` | 1e-4 | `armijo` |
| `max_log_composition_step` | 4 | `log_step_limit` |
| `max_iterations` | 512 | `max_iterations` |
| `max_backtracks` | 32 | `max_backtracks` |
| `max_property_evaluations` | 100000 | `max_evaluations` |
| `automatic_multistart` | true | `automatic_starts` |
| `max_starts` | 1024 | `max_starts` |

The groups map respectively to initial split stability, final split stability,
and final three-phase stability. TPD tolerance is absolute D/(RT); stationarity
uses the max unweighted chemical-potential residual/(RT). The existing final
common-reference allowances remain inside the solver: resolved settings record
the configured base tolerances and do not remove those scientific semantics.

| `two_phase` field | Default v1 | Native option |
| --- | ---: | --- |
| `fugacity_equilibrium_tolerance` | 1e-11 | `iteration.fugacity_tolerance` |
| `absolute_mass_balance_tolerance` | 1e-12 | `iteration.mass_absolute_tolerance` |
| `relative_mass_balance_tolerance` | 1e-10 | `iteration.mass_relative_tolerance` |
| `minimum_phase_fraction` | 1e-10 | `iteration.minimum_phase_fraction` |
| `minimum_log_composition_separation` | 1e-7 | `iteration.log_k_separation` |
| `minimum_relative_z_separation` | 1e-8 | `iteration.relative_z_separation` |
| `max_log_equilibrium_ratio_step` | 2 | `iteration.max_log_step` |
| `residual_progress_coefficient` | 1e-4 | `iteration.residual_decrease` |
| `gibbs_progress_coefficient` | 1e-4 | `iteration.gibbs_armijo` |
| `max_iterations` | 512 | `iteration.max_iterations` |
| `max_backtracks` | 24 | `iteration.max_backtracks` |
| `max_property_evaluations` | 20000 | `iteration.max_evaluations` |
| `rachford_rice_max_iterations` | 192 | `iteration.rr.max_iterations` |
| `max_split_attempts` | 16 | `max_split_attempts` |

Composition separation is max absolute log composition ratio, not a raw mole
fraction difference. Relative mass balance applies per positive feed component,
without an invented floor. Progress coefficients do not replace acceptance
tolerances. The two-phase property budget is accumulated across split attempts
by the existing implementation.

| `three_phase` field | Default v1 | Native option |
| --- | ---: | --- |
| `chemical_potential_tolerance` | 1e-11 | `three_phase.chemical_potential_tolerance` |
| `absolute_mass_balance_tolerance` | 1e-12 | `three_phase.mass_absolute_tolerance` |
| `relative_mass_balance_tolerance` | 1e-10 | `three_phase.mass_relative_tolerance` |
| `generalized_rr_balance_tolerance` | 2e-13 | `three_phase.balance_tolerance` |
| `minimum_phase_fraction` | 1e-10 | `three_phase.minimum_phase_fraction` |
| `minimum_log_composition_separation` | 1e-7 | `three_phase.log_composition_separation` |
| `max_log_step` | 2 | `three_phase.max_log_step` |
| `residual_progress_coefficient` | 1e-4 | `three_phase.residual_decrease` |
| `max_iterations` | 512 | `three_phase.max_iterations` |
| `max_line_search_backtracks` | 32 | `three_phase.max_backtracks` |
| `max_balance_iterations` | 192 | `three_phase.max_balance_iterations` |
| `max_balance_backtracks` | 48 | `three_phase.max_balance_backtracks` |
| `max_property_evaluations` | 30000 | `three_phase.max_evaluations` |
| `max_three_phase_attempts` | 16 | `max_three_phase_attempts` |
| `new_phase_seed_fraction` | 0.1 | `new_phase_seed_fraction` |

Three-phase property evaluations are bounded per attempt. Attempt limits remain
independent for caller-supplied hints and automatic negative-TPD witnesses; the
settings mapper does not combine these budgets or reinterpret hints as phase
evidence.

## Validation and host safety ceilings

The adapter checks presence, finite values, integer conversions and host policy.
It delegates scientific option-domain checks to the existing stability, split,
three-phase and max3 validators. Root iterations use the existing positive-only
evaluator rule. No EOS evaluation is performed to validate settings. Domain
errors identify the public group; missing/nonfinite/count errors identify the
individual field. Unknown versions/presets and edited-but-preset-labeled input
fail explicitly. Allocation failures and programming errors propagate.

Existing distinctions remain: TPD tolerance and iteration counts may be zero;
root iterations/backtracks must be positive; two-phase evaluation and attempt
budgets can be zero; stability and three-phase evaluation budgets must be
positive. Exhaustion is delegated to the solver and can produce indeterminate,
never automatic success or repaired defaults.

`PtSolverSafetyLimits` is supplied by the **host**, separately from the public
settings. Its defaults are resource policy, not scientifically fitted values:

| Host ceiling | Default | Application |
| --- | ---: | --- |
| Root iterations | 8192 | Root setting upper bound |
| Iterations | 4096 | Stability, split, RR and three-phase/balance iteration upper bounds |
| Backtracks | 256 | Stability, split and three-phase/balance upper bounds |
| Property evaluations | 1000000 | Upper bound for each configured evaluation budget |
| Stability starts | 1024 | Upper bound for each group's `max_starts` |
| Split attempts | 64 | Two-phase attempt upper bound |
| Three-phase attempts | 64 | Each initialization-source class |
| Components | 256 | Injected into stability, RR and three-phase native options |
| Stability start entries | 262144 | Injected separately into all three stability groups |
| Three-phase hint count | 16 | Injected into max3 options |
| Three-phase hint entries | 12288 | Injected into max3 options |

Zero host ceilings are invalid. A stricter host may reject the standard preset;
it never silently clips it. The snapshot exposes the applied host limits so the
effective configuration is queryable. Nested ceilings are not an aggregate
flash-time budget or cancellation facility. The executable model factory checks
the model component count against this policy, and service admission must retain
the existing serialized evaluator-workspace contract.

## Per-solve hints and state-dependent validation

`Pr76SolverConfiguration` deliberately keeps the prepared backend's
initial/final/three-phase start lists empty. Starts are **not settings**: they
belong to a particular solve state and must not make the immutable settings
snapshot depend on a previous P/T/z point.

The executable PR76 model now exposes the separate stdlib-only
`pt-solve-hints/v1` DTO and `solve(request, hints)` overload documented in
[the executable model contract](executable_model.md). That per-call adapter copies
the already prepared native options, injects the supplied starts into the existing
native fields, and leaves the stored settings/coarse backend unchanged.

Version and host storage ceilings are checked at the public boundary before hint
storage is copied. Feed support, exact component dimension, composition
normalization, generated-start capacity, final-two-phase reservation, three-phase
simplex feasibility, and native attempt budgets remain state dependent and are
validated by the existing solver against the **current request**. Preparing
settings therefore never certifies that every possible P/T/z request or hint set
is solvable. Disabling automatic initial multistart still requires a valid
per-solve initial start under the existing empty-search rule.

This division is intentional: public settings describe persistent numerical
policy; public solve hints describe ephemeral initialization. Hints remain search
initializers only and do not become thermodynamic evidence or bypass equilibrium,
material-balance or final stability gates.

## Scope and verification

The independent test project `tests/model_configuration/solver_settings` covers
all 57 fields for missing-value rejection; a separately constructed direct C++
configuration with distinct values tests every numeric mapping and each host
injection. Frozen literal and current-native-default comparisons detect default
drift. Other tests cover preset identity, domain/quotas, zero/false preservation,
independent stability groups, snapshot isolation and public-header independence.

Actual PR76 solves compare outcome, phase count, phase fractions, compositions,
Z, ln(phi), diagnostic and provenance for single-/three-phase results and root,
initial stability, final two-phase and final three-phase budget exhaustion.
They reuse `tests/flash/pr76_three_phase/synthetic_fixture.hpp`; seeded comparisons
apply identical native hints to both paths. Exact same-platform equality is
required. This validates settings mapping, not independent EOS accuracy.

The separate executable-model suite adds direct-C++ public-hint parity,
accepted-result continuation parity, unresolved-state clearing, and nested
state-dependent hint validation against independently assembled native options.
At commit `6db1ae25`, the official model-configuration workflow passed GCC Debug
+ ASan/UBSan, Clang Release and MSVC Release; GCC reported 17/17 executable-model
cases and the existing registry remained 17/17.

The focused official-hosted workflow keeps the solver-settings 14-case suite on
all three compilers alongside parameter/executable/registry coverage. Commands:

```bash
cmake -S tests/model_configuration/solver_settings -B build/solver-settings \
  -DCMAKE_BUILD_TYPE=Debug -DMPMC_SOLVER_SETTINGS_ENABLE_SANITIZERS=ON
cmake --build build/solver-settings --parallel 2
ctest --test-dir build/solver-settings -R '^model[.]solver_settings[.]' --verbose --no-tests=error
```

No existing EOS/flash formula, tolerance default, v1 wire meaning, configured
PR/SW/CPA backend or product call site is changed. An additive settings module
was chosen over exposing native Options types or changing existing defaults.
Executable model ownership and per-solve hints are provided by
[the PR76 factory](executable_model.md). The optional [bounded registry](registry.md)
adds handle/release lifetime. Hint serialization through registry/service/desktop
transport, API/UI controls, one-sided applicability and full product-facing
configuration work remain later Gate items; UI stays deferred.
