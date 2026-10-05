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
