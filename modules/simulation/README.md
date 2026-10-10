# Simulation module

`mpmc::simulation` owns simulation-lifecycle semantics above the existing
mesh, thermodynamics/flash, flow, well and solver layers. It decides **when an
already materialized physical system is advanced**; it does not reimplement
residuals, Jacobians, phase stability, well physics, timestep adaptation or
PETSc nonlinear/linear solves.

The first production slice is deliberately narrow: a PETSc-independent hard
timeline plus a multi-timestep PETSc runner that repeatedly consumes the
existing
`flow_discretization_petsc::advance_one_physical_timestep_3d()` contract.
This closes the gap between one accepted physical timestep and a simulation
that reaches exact report/end-time boundaries.

## Ownership and dependency direction

The intended dependency direction is:

```text
simulation_petsc
    |
    +--> simulation
    +--> flow_discretization_petsc
            |
            +--> flow / discretization / mesh / thermodynamics

flow, well, mesh and thermodynamics never depend on simulation.
```

The core `simulation` layer remains free of PETSc/MPI types. The PETSc bridge
may consume the accepted reservoir system, accepted physical-time clock and
phase-transition bindings already owned by `flow_discretization_petsc`.

Simulation does **not** own:

- component/energy conservation equations or constitutive models;
- SNES/KSP/PC configuration or adaptive cutback/growth policy;
- post-SNES phase scanning or same-dt phase rebuilds;
- well connection physics or rate/BHP arbitration;
- mesh topology, geometry, transmissibility or partitioning.

Those responsibilities remain in their existing modules.

## Time model

All production simulation time is SI seconds.

A simulation timeline is an ordered set of finite, strictly increasing **hard
time boundaries**. The first version uses these boundaries only for report/end
times; later control/output/checkpoint payloads may bind to the same accepted
time semantics without changing the core clock contract.

For an accepted time `t_n` and next hard boundary `t_h > t_n`, the simulation
runner supplies

```text
initial_timestep_cap_seconds <= t_h - t_n
```

to the existing physical-timestep driver. Any adaptive retry/cutback remains
inside that driver. A trial solve, rejected attempt or same-dt phase-transition
restart cannot advance the simulation cursor.

A hard boundary is consumed only after the authoritative
`AcceptedPhysicalTimeClock3D` has committed an accepted step at that
boundary. The simulation layer never mutates the accepted reservoir state or
physical clock independently.

If a caller already supplies a stricter physical-timestep cap, the runner must
preserve it by using the minimum of the caller cap and the remaining distance
to the next hard boundary.

## Implemented T1 core API

The PETSc-independent core now exposes:

- `SimulationTimeline`: an immutable value object holding finite,
  non-negative, numerically distinct and strictly increasing hard-boundary
  times;
- `SimulationCursor`: accepted lifecycle state with exactly one
  next-unconsumed-boundary index;
- `SimulationHardBoundary`: stable boundary index/time identity.

A cursor initialized from an already accepted time treats boundaries at or
before that accepted time as historical and emits no synthetic reach event.
During forward execution, `accept_time()` advances at most one boundary and
rejects an accepted time that would skip an unconsumed boundary.
`cap_initial_timestep_seconds()` returns the stricter of the caller cap and
the remaining distance to the next hard boundary.

Time-point identity uses the same `64 * epsilon * max(1, |t_a|, |t_b|)`
scale already used by the physical-time layer. This tolerance is only an
identity guard for hard times; it is not a nonlinear, conservation or adaptive
timestep tolerance.

The core is an INTERFACE C++20 target `mpmc::simulation` with standard-library
dependencies only. Its standalone regression project is
`tests/simulation/core`.

```sh
cmake -S tests/simulation/core -B build/simulation-core -DCMAKE_BUILD_TYPE=Release
cmake --build build/simulation-core --parallel 2
ctest --test-dir build/simulation-core --output-on-failure
```

The current regression owns seven behavior cases plus independent
self-contained probes for both public headers. Central CI ownership is not yet
claimed in this Draft slice; it remains an explicit T6 completion item for the
full PR.

## Implemented T2 report/status API

The PETSc-independent lifecycle report is intentionally a **bounded projection**
rather than a solver-history container.

- `AcceptedSimulationStepSummary` carries accepted step index, accepted time
  interval, accepted/next timestep, retry count, attempt count and phase-restart
  count. It does not contain residual vectors, nonlinear norms or KSP details.
- `ReachedSimulationBoundary` binds one hard-boundary identity to the accepted
  physical step that actually reached it.
- `SimulationRunStatus` distinguishes `running`, `completed`, terminal
  physical-timestep rejection, lower-layer execution error and lifecycle
  contract violation without importing PETSc enums.
- `SimulationReport` stores fixed-size lifecycle diagnostics, the **last**
  accepted-step summary, the **last** reached boundary and an optional signed
  lower-layer error code. It never accumulates an unbounded vector of steps,
  attempts or solver diagnostics.

A report starts from `SimulationReportAnchor` containing authoritative accepted
time, accepted-step count and already-consumed boundary count. New accepted
steps must be sequential in both step index and accepted time. Reached
boundaries must be sequential from the anchor and match the endpoint of the
last accepted step. Every accepted step must correspond to exactly one new
physical-timestep call.

The lower-layer error code is valid only for `execution_error`; other terminal
states do not smuggle backend-specific diagnostics into the core contract.
After any terminal status, further lifecycle records are rejected.

T2 extends the standalone core suite to 18 CTests: seven timeline/cursor
behaviors, seven report/status behaviors and four public-header self-contained
probes.

## Implemented T3 PETSc runner

The first PETSc bridge is `mpmc::simulation_petsc`, an INTERFACE C++20 target
that depends on `mpmc::simulation` and the existing
`mpmc::flow_discretization_petsc` bridge.

`advance_simulation_timeline_3d()` accepts an already materialized
`PhaseTransitionRebuiltNaturalVariableSystem3D`, authoritative
`AcceptedPhysicalTimeClock3D`, `SimulationCursor`, transition bindings and
the existing physical-timestep options. For each internal step it:

1. reads the next unconsumed hard boundary;
2. composes the boundary distance with any existing caller
   `initial_timestep_cap_seconds` using the stricter cap;
3. lowers only the per-call adaptive minimum to the strictest of the caller
   floor, the current hard-boundary/caller cap and the authoritative
   `clock.next_timestep_seconds()`; this admits a sub-minimum next-dt inherited
   from a prior exact boundary clip without allowing a retry below that inherited
   floor, and does not mutate the caller's persistent adaptive policy;
4. calls `advance_one_physical_timestep_3d()` exactly once;
5. validates the accepted record, adaptive report and authoritative clock
   against each other;
6. projects the lower report into the bounded
   `AcceptedSimulationStepSummary` and discards the detailed lower report;
7. advances the `SimulationCursor` only when the accepted clock actually
   reaches the next hard boundary;
8. repeats until the terminal boundary is consumed.

A normal lower-level terminal timestep rejection maps to
`physical_timestep_rejected` only if accepted time, accepted step count and
next suggested timestep are unchanged. A nonzero PETSc/MPI return maps to
`execution_error` only when the accepted clock remains unchanged. Any
contradiction between the lower report, accepted clock and hard-boundary
contract terminates as `contract_violation`.

The existing Flow `physical_time_loop.hpp` remains unchanged. It is a
lower-layer single-target helper with an unbounded vector of detailed physical
step reports; Simulation does not reuse that report container because the
lifecycle contract is intentionally bounded and supports multiple hard
boundaries.

The current integration consumer extends the already-owned real PR76 transient
fixture: from accepted `t=2.5`, hard boundaries `3.0/3.5` and caller cap
`0.4` require four accepted internal physical steps, two reached boundaries,
terminal clock `3.5`, and accepted-history rebase on the final state. The
source is wired into the existing PETSc test target; runtime execution still
requires the pinned PETSc/MPI environment.

## T5 transactional regression — retry/cutback slice

The first T5 slice uses the existing real PR76/PETSc production fixture rather
than a mock solver. A test-only scanner wrapper makes the first converged trial
return `indeterminate`, forcing the existing physical-timestep driver to
discard that trial and execute its normal adaptive cutback.

The wrapper directly observes the authoritative physical clock and
`SimulationCursor` at three scanner points:

1. the first trial before the forced rejection;
2. the cutback retry before its accepted commit;
3. the next physical step after the cutback step was accepted but before the
   hard boundary is reached.

The first two observations must retain the entry accepted time, entry accepted
step count and unconsumed hard-boundary cursor. The third must show exactly one
accepted physical step of progress while the boundary is still unconsumed.

For the frozen regression, accepted `t=3.5`, hard boundary `3.6`, initial cap
`0.1`, minimum dt `0.025` and cutback factor `0.5` produce:

- first trial `dt=0.1`: forced recoverable rejection, no accepted commit;
- retry `dt=0.05`: accepted at `t=3.55`, cursor still unconsumed;
- next physical step `dt=0.05`: accepted at `t=3.6`, cursor advances once.

The bounded simulation report must contain two physical calls, two accepted
steps, one accepted retry and one reached boundary. Accepted component/energy
history must match the terminal committed state.

## Accepted-state transaction

The first implementation has three accepted authorities:

1. the accepted natural-variable system;
2. `AcceptedPhysicalTimeClock3D`;
3. the simulation timeline cursor.

The first two continue to be committed exclusively by the existing
`PhysicalTimestepDriver`. The timeline cursor may advance only after the
returned report contains a valid accepted timestep record consistent with the
accepted clock.

On terminal rejection or PETSc/MPI error:

- accepted reservoir history remains unchanged according to the existing flow
  contract;
- accepted physical time remains unchanged;
- the timeline cursor remains unchanged;
- no hard boundary is reported as reached.

Phase-transition restarts inside one physical timestep are invisible to the
timeline except through the final accepted record.

## Initial v1 result

The v1 runner is expected to advance an already materialized reservoir system
from its current accepted time through one or more hard boundaries until a
requested terminal boundary is reached or the existing physical-timestep
driver returns a terminal failure.

The report should be bounded and audit-friendly: accepted physical-step
records, reached boundary indices/times and the terminal status. It must not
copy full field state or duplicate lower-level SNES/KSP diagnostics already
owned by the physical-timestep reports.

## Explicit non-goals for v1

This first slice does not add:

- YAML/JSON/ECLIPSE-style input parsing;
- case materialization or `p/T/z` initialization;
- boundary-condition physics;
- well schedule changes or multi-well scheduling;
- summary/VTU field output;
- checkpoint/restart serialization;
- a generic event bus or plugin system;
- new EOS/flash models, flow equations or solver settings.

Those capabilities require separate contracts after the timeline/runner
transaction is proven.

## Verification boundary

The PETSc-independent timeline tests must cover invalid/non-finite time input,
strict ordering, exact boundary selection, cap composition and cursor
progression without requiring PETSc.

The PETSc integration regression must prove at minimum:

- an adaptive timestep cannot cross the next hard boundary;
- multiple accepted internal timesteps can land exactly on a boundary;
- retry/cutback does not advance the simulation cursor;
- terminal failure leaves accepted time and cursor unchanged;
- an accepted same-dt phase-transition restart still advances the cursor only
  once, after the final phase-stable commit;
- serial and the existing affected MPI ownership path preserve the same
  accepted boundary sequence where that path is applicable.

Tests reuse existing physical-timestep scientific fixtures and solver
contracts; this module does not create a new physical oracle.
