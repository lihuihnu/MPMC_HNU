# Flow: natural variables, conservation and time integration

`mpmc::flow` supplies the local contracts for fully implicit, non-isothermal,
multicomponent flow with one, two or three active fluid phases. Spatial assembly
and PETSc integration live in separate subdirectories. Phase slots are numerical
coordinates; a physical phase identity must come from explicit model evidence.

This page is the development entry point. The complete equations, assumptions,
references and historical implementation records are preserved in
[the contract reference](contracts.md). Its numbered sections describe the
scope of their original increments; they are not new implementation requests.

## Follow a calculation through the code

| Responsibility | Start here | Boundary |
| --- | --- | --- |
| Natural-variable layout and state | [natural_variable_cell_state.hpp](include/mpmc/flow/natural_variable_cell_state.hpp), [1P](include/mpmc/flow/single_phase_natural_variable.hpp), [2P](include/mpmc/flow/two_phase_natural_variable.hpp) | Component order, active phases, pivots, strict positive support |
| EOS and property adapters | [PR76 closure](include/mpmc/flow/pr76_selected_phase_property_closure.hpp), [SW92 closure](include/mpmc/flow/sw92_selected_phase_property_closure.hpp) | Frozen selected phases, sourced properties, value/derivative consistency |
| Local storage | [component accumulation](include/mpmc/flow/component_accumulation.hpp), [backward Euler](include/mpmc/flow/component_accumulation_time.hpp), [energy accumulation](include/mpmc/flow/energy_accumulation.hpp) | Frozen accepted history; local residual and Jacobian |
| Mobility and face transport | [phase transport](include/mpmc/flow/phase_transport.hpp), [phase potential](include/mpmc/flow/phase_potential_upwind.hpp), [component flux](discretization/include/mpmc/flow_discretization/tpfa_component_molar_flux.hpp) | TPFA admissibility, frozen upwind choice, owner-to-neighbour orientation |
| Cell conservation | [component residual](discretization/include/mpmc/flow_discretization/local_component_conservation_residual.hpp), [energy residual](discretization/include/mpmc/flow_discretization/local_energy_conservation_residual.hpp) | Bulk-volume normalization and matching state identities |
| Distributed nonlinear assembly | [mixed-cardinality assembly](discretization/petsc/include/mpmc/flow_discretization_petsc/mixed_cardinality_physical_snes_assembly.hpp) | Frozen 1P/2P/3P chart during each SNES solve; owner-only contributions |
| Physical timestep | [timestep driver](discretization/petsc/include/mpmc/flow_discretization_petsc/physical_timestep_driver.hpp), [outer rebuild](discretization/petsc/include/mpmc/flow_discretization_petsc/phase_transition_outer_rebuild.hpp) | Candidate solve, phase scan, same-dt restart, final accepted-history commit |
| Well physics and controls | [Well module](../well/README.md) | Well index, completion aggregation and rate/BHP arbitration are owned by Well |

Mesh owns topology, geometry and fields. Discretization owns TPFA admissibility
and transmissibility. Flow consumes these results; its local core does not
introduce PETSc/MPI types. AD and thermodynamics keep their existing independent
module boundaries. See [cross-module ownership](contracts.md#14-cross-module-ownership).

## State checks and numerical policies

The internal [validation helpers](include/mpmc/flow/detail/validation.hpp) contain
shared comparison algorithms. Call sites retain their numerical policy and
context-specific failure messages:

| Comparison | Existing policy retained |
| --- | --- |
| Local transport and storage consistency | `4096 * epsilon` with the caller's scale |
| Component flux, scatter and local component closure | `8192 * epsilon` with the caller's scale |
| Energy/distributed closure where previously selected | `16384 * epsilon` with the caller's scale |
| Well component-rate comparison | Separate implementation with its `1e-30` scale floor |

These constants are existing consistency allowances, not a universal solver
convergence tolerance. The shared unit-floor comparison rejects non-finite input
and invalid extra scales; it does not clip, renormalize or repair a state.

- A layout identity includes component count, active phase count, unknown count
  and dependent-component pivots.
- Component names are compared as ordered, borrowed sequences. Comparing a
  snapshot with a cell state requires no temporary vector of strings.
- Active pressure, temperature, saturation and composition use the caller's
  comparison policy. Inactive phase sidecars remain **exactly zero/empty**.
- Equality does not prove physical validity: positive support, material balance,
  energy conservation, admissibility and stability retain their own checks.
- Local wrappers remain available to existing callers; the helper is an internal
  implementation detail, not a new application API.

## Scientific and capability boundaries

The natural-variable sizes are `Nc + 1`, `2*Nc + 1` and `3*Nc + 1` for
one, two and three phases. A solve never changes its active chart inside a
residual/Jacobian callback. Phase changes are handled outside SNES, rebuilding
from stable cell/phase identities and the same accepted history.

The default nonlinear path remains Newton with GMRES and ASM. The
variable-cardinality ASM subblock LU explicitly uses MUMPS without factor shift;
this refactor does not change solver settings, conservation tolerances or
strict-positive acceptance.

The SW92 transactional restart currently rejects authoritative faces and local
ghost overlap because no validated family/root-aware absent-phase transport
extension is available. The CO2/H2O production-property regression and Sample-6
serve different purposes: Sample-6 uses real SW92 thermodynamics with
**manufactured transport/caloric properties** to test restart and conservation
software. It is not an independently validated Sample-6 transport model.
See [SW92 restart details](contracts.md#52-sw92-transactional-same-dt-phase-transition-restart).

## Build and verify the layer you change

Run these commands from the repository root:

```sh
cmake -S tests/flow/core -B build/flow-core -DCMAKE_BUILD_TYPE=Release
cmake --build build/flow-core --parallel 2
ctest --test-dir build/flow-core --output-on-failure

cmake -S tests/flow_discretization/core -B build/flow-discretization -DCMAKE_BUILD_TYPE=Release
cmake --build build/flow-discretization --parallel 2
ctest --test-dir build/flow-discretization --output-on-failure
```

The second project includes Well discretization tests. PETSc/MPI assembly and
restart tests belong to `tests/flow_discretization/petsc`; use the pinned toolchain
and commands in [the existing PETSc workflow](../../.github/workflows/flow_discretization_petsc.yml).
The root CMake project configures AD only, so a successful root build is not a
Flow verification.

The shared comparison regression is `flow.core.shared_validation`, owned by
`mpmc_flow_core_tests`. It covers threshold boundaries, non-finite values, ordered
identities, pivots and inactive sidecars. Existing constitutive, AD/Jacobian,
conservation and restart regressions remain their original owners.

A change to the shared validation header selects Flow core, Flow discretization
and Flow PETSc through the [central impact router](../../.github/ci/impact_rules.py).
No additional workflow is needed. Follow [AGENTS.md](../../AGENTS.md) for the
project's audit, provenance, testing and submission requirements.
