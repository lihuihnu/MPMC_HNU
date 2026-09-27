# Mixed-cardinality PETSc regression

This directory contains private source fragments of
`../mixed_cardinality_physical_snes_assembly_test.cpp`. Start with that file to
read the assembled-system assertions and the scenario call order.

| File | Responsibility |
| --- | --- |
| `fixture.inc` | Collective assertions, target states, cell evaluators, mesh and phase identities |
| `thermodynamic_adapter.inc` | Manufactured PR data, absent-phase extension and rebuild inputs |
| `phase_transition.inc` | Controller fixture, scanner, rebuild and lifecycle checks |
| `well_fixture.inc` | Shared well timestep setup, conserved totals and rate probes |
| `fixed_bhp.inc` | Frozen/rebound single- and multi-connection BHP cases |
| `rate_control.inc` | Rate/BHP arbitration, hysteresis and transactional phase transitions |
| `fixed_bhp_transition.inc` | Explicit BHP source rebind across phase appearance/disappearance |

The parent includes these files in dependency order inside its anonymous namespace.
They are not standalone headers or additional translation units. Helpers and state
remain private and are defined once, shared by the existing cases; other test suites
must not include these fragments. This preserves compiler visibility, symbol linkage,
fixture lifetimes, collective call order and optimization opportunities. It does not
claim reduced compilation time or faster solver execution.

Owner: `mpmc_flow_discretization_petsc_tests`.
CTest: `flow_discretization.petsc.distributed_owner_targeted_conservation` (two ranks).
The existing `mixed_cardinality_physical_snes_assembly_test()` entry remains called
from the executable main. All assertions, tolerances and manufactured provenance
are retained. Changes anywhere under this directory route to the owning PETSc Gate.

Validation of the initial split expanded the includes and compared the result
byte-for-byte with the previous monolithic source. Runtime validation uses the same
Gate; no additional test executable or workflow is introduced.
