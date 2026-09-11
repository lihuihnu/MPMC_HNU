# Physics thermodynamic-closure consumption contract

## Scope

The `physics` layer consumes already accepted thermodynamic/flash states and exposes owned thermodynamic snapshots for future conservation-law and discretization code. It does **not** create mesh, flux, saturation, Darcy flow, time stepping, global residual assembly, a nonlinear solver, viscosity, enthalpy or an energy equation.

Two adapters are currently implemented:

- `mpmc/physics/pr76_thermodynamic_closure.hpp`: accepted PR76 internal two-phase VLE primal plus its validated local linearization;
- `mpmc/physics/sw92_thermodynamic_closure.hpp`: accepted SW92 Profile-C authoritative 1/2/3-phase PT primal plus the validated fixed-phase-set local implicit linearization when available.

The generic header `mpmc/physics/thermodynamic_closure.hpp` contains two payload shapes:

- `ThermodynamicClosureSnapshot`: the original fixed liquid/vapor PR76 VLE contract;
- `PtPhaseSetThermodynamicClosureSnapshot`: the variable-cardinality phase-set reduced-feed v2 contract used by SW92 and intended for later phase-set models.

The CMake target remains `mpmc::physics`. Dependency direction is `physics -> flash_sensitivity -> ad/thermodynamics`; no lower layer depends on physics.

## Atomic snapshot and two-axis availability

Primal availability and linearization availability are independent decisions.

A valid primal may be used by a future conservation residual or line-search merit evaluation. A Newton/Jacobian base point may consume a derivative only when the **same owned closure snapshot** publishes a supported linearization.

No adapter may replace an unavailable derivative with zeros, a previous-cell or previous-iterate derivative, clipping, or a hidden finite-difference fallback.

### PR76

For the existing PR76 two-phase adapter, an accepted two-phase flash can remain a valid primal when sensitivity is unavailable because of a phase-boundary guard, local ill-conditioning, unsupported feed support, derivative property failure or derivative arithmetic failure. In those cases `residual_available()==true`, `can_seed_newton()==false`, and `linearization` is empty.

If the sensitivity path reports `solution_not_accepted`, or the flash itself is unstable/indeterminate, the adapter does not publish a usable primal. Stable PR76 single-phase closure remains outside that v1 adapter.

### SW92 Profile-C

The SW92 adapter now follows the same two-axis policy for authoritative 1/2/3-phase states.

When the fixed-phase-set flash sensitivity succeeds:

```text
primal_status        = valid
linearization_status = available
linearization_reason = none
residual_available() = true
can_seed_newton()    = true
```

When the authoritative primal is valid but the local derivative is unavailable because of a phase-boundary guard, local ill-conditioning, unsupported positive-support chart, derivative property failure or derivative arithmetic failure:

```text
primal_status        = valid
linearization_status = unavailable
residual_available() = true
can_seed_newton()    = false
```

If the sensitivity path reports `solution_not_accepted`, the derivative path has failed to reproduce the supposedly authoritative owned base state. That is treated as an atomic consistency failure: both primal and linearization are withdrawn and the snapshot is `indeterminate`.

## Generic variable-cardinality phase-set v2

`PtPhaseSetThermodynamicClosureSnapshot` owns:

- pressure `p` and temperature `T`;
- ordered overall feed and ordered component IDs;
- thermodynamic model / dataset / revision identity;
- a vector of `ThermodynamicPhaseState` objects;
- an optional `PtPhaseSetThermodynamicLinearization`.

Each phase primal state contains:

```text
mole phase fraction beta_alpha
ordered mole composition x_alpha
compressibility factor Z_alpha
molar density c_alpha = p/(Z_alpha R T)
```

The phase fraction is a **mole phase fraction, not pore-volume saturation**.

The generic linearization columns are:

```text
q = (p_Pa, T_K, z_0, ..., z_{N-2})
z_{N-1} = 1 - sum(z_0, ..., z_{N-2})
```

and it stores, in the exact primal representation order:

```text
d beta_alpha / dq
d x_alpha,i / dq
d Z_alpha / dq
d c_alpha / dq
```

plus the local equilibrium-Jacobian reciprocal-condition estimate, implicit-solve backward error and accepted-base residual norm.

`can_seed_newton()` checks not only status and payload presence but also component/phase/input dimensions, all Jacobian shapes, finite entries and finite local-solve diagnostics. The generic vector deliberately has no `liquid`, `vapor`, `aqueous`, AQ or NA label. Model-specific role/family identity remains in the model adapter sidecar so generic physics code cannot silently turn a numerical phase slot into physical morphology.

## SW92 authoritative Profile-C adapter

Public entry:

```cpp
#include <mpmc/physics/sw92_thermodynamic_closure.hpp>

const auto closure = mpmc::physics::build_sw92_profile_c_thermodynamic_closure(
    authoritative_phase_set, sw92_model);
```

The input must be an authoritative `Sw92ProfileCPtPhaseSetResult` from the Profile-C publication contract. The adapter does not run another flash, TPD search, material-balance solve, topology search, composition normalization or morphology classifier.

### Provenance and phase identity

The adapter preserves and checks:

- p/T and ordered feed;
- dataset ID / revision / ordered component IDs;
- prescribed NaCl molality;
- SW92 model and phase convention;
- Profile-C equilibrium/orchestration/boundary/publication conventions;
- the SW92 fixed-phase-set sensitivity convention;
- per-phase physical role and AQ/NA family;
- per-phase selected cubic-root branch as numerical provenance.

The SW sidecar retains only the roles already authorized upstream:

- `aqueous` with AQ family;
- `nonaqueous_unclassified` with NA family.

H0/H1 remain `nonaqueous_unclassified`; branch index or Z ordering is not promoted to liquid/vapor or LV/LL morphology. `morphology_resolved=false` and `global_stability_proven=false` remain explicit.

### Primal reproduction and density

Before publishing the physics primal, the adapter evaluates the existing SW92 phase-property kernel at each **exact accepted** composition, family and root branch and requires the reproduced `ln(phi)` to match the activity retained by the authoritative flash result within a roundoff guard.

For accepted 2/3-phase states whose authoritative publication already owns Z, the reproduced Z must agree and the physics state preserves the exact accepted Z.

The authoritative dry no-W single-H publication intentionally stores reference activity but no selected Z, because the publication adapter remains a pure projection. The physics closure needs a molar volume/density, so for that one case it obtains Z from the already identified family/root branch during the same property-reproduction call. It does **not** rerun phase selection or flash logic.

Molar density introduces no new empirical model:

```text
Z = p v_m/(R T)
c = 1/v_m = p/(Z R T)   [mol/m^3]
```

No mass density, saturation, mobility or transport property is implied.

### Atomic SW92 linearization handoff

After the primal is reproduced, the adapter invokes:

```cpp
differentiate_sw92_profile_c_phase_set(authoritative_phase_set, sw92_model, ...)
```

with the same selected-root options used for primal reproduction. No physics-layer differentiation is reimplemented.

A successful sensitivity is accepted only if it still matches the owned primal snapshot in:

- p/T and ordered feed;
- dataset/revision/component order;
- model/phase/equilibrium/orchestration/boundary/publication conventions;
- phase count and SW physical-role / thermodynamic-family sidecar;
- reduced-feed input count and every Jacobian shape;
- finite Jacobian entries, positive finite reciprocal condition estimate, finite backward error and finite base residual.

The validated flash arrays are then copied unchanged into `PtPhaseSetThermodynamicLinearization`. This is an ownership boundary, not a second derivative algorithm.

The differentiated domain remains exactly the flash sensitivity domain: fixed accepted phase count, fixed AQ/NA family assignment, fixed selected cubic-root branch and fixed representation slots. It does not differentiate TPD, SSI, generalized RR, line search, topology routing, phase appearance/disappearance, family/root switching or H morphology classification.

### Failure policy

A non-accepted authoritative status, malformed phase set, inconsistent publication convention, invalid role/family pair, phase-fraction/composition shape error, or unreproducible activity/Z yields an `indeterminate` closure with no residual primal.

Passing a different ordered SW92 model snapshot (dataset/revision/component order mismatch) is an input-contract error and throws `std::invalid_argument`, consistent with the PR76 adapter.

Derivative-only statuses `phase_boundary`, `ill_conditioned_equilibrium`, `unsupported_feed_support`, `property_failure` and `arithmetic_failure` preserve the valid primal but publish no linearization. `solution_not_accepted` invalidates the entire snapshot because the sensitivity and owned primal disagree on the base state.

## PR76 coordinates and local linearization

The existing PR76 derivative columns remain:

```text
q = (p_Pa, T_K, z_0, ..., z_{N-2})
z_{N-1} = 1 - sum(z_0, ..., z_{N-2})
```

The PR76 primal is anchored to the exact accepted flash `beta/x/y/Z` snapshot. Before publication, the adapter re-evaluates the selected branches with the accepted full phase compositions and requires reproduced Z agreement.

For an available linearization, the adapter obtains local selected-branch Z partials from the existing AD/simple-root IFT phase kernel and combines them with converged flash composition derivatives. Molar-density derivatives use the exact `c=p/(ZRT)` identity.

No flash iteration, RR bisection, line search, TPD search or root iteration is differentiated by the closure adapter itself.

## Validation

The SW92 physics-closure regression covers:

1. traceable dry CO2/H2O authoritative single-H state, including reconstruction of the source-missing Z and successful local linearization;
2. traceable wet CO2/H2O authoritative W+H state with exact accepted primal and exact handoff of the validated flash Jacobian;
3. physical Mortezazadeh–Rasaei Sample-6 authoritative W+H0+H1 primal plus successful fixed-phase-set local linearization;
4. `c=p/(ZRT)` for every published phase;
5. AQ/NA and `aqueous` / `nonaqueous_unclassified` metadata preservation;
6. direct equality between physics linearization arrays and the independently called flash-sensitivity payload;
7. derivative phase-boundary guard: primal remains available while `can_seed_newton()==false`;
8. sensitivity/base inconsistency: `solution_not_accepted` invalidates the entire atomic snapshot;
9. publication-convention and phase-activity tamper rejection;
10. ordered model/source mismatch rejection;
11. runtime component permutation for both primal and p/T derivative columns;
12. public-header self containment and reduced-feed v2 convention.

The underlying flash-sensitivity regression independently validates the fixed-phase-set implicit Jacobian with forward AD, fresh re-solves, the physical Sample-6 three-phase state, permutation/H-slot invariance and derivative boundary guards; see the [SW92 Profile-C sensitivity contract](../flash/sw92_profile_c_sensitivity.md).

A dedicated GitHub-hosted GCC Debug+ASan/UBSan / Clang Release / MSVC Release workflow regenerates the existing independent Sample-6 Decimal(80) primal oracle before the focused C++ suite. Changes to the shared generic closure header also select the existing PR76 closure workflow.

## Current non-capabilities

The physics layer still does not provide:

- a validated SW92 H0/H1 liquid/vapor or LV/LL morphology resolver;
- derivatives across phase appearance/disappearance, family/root switching or critical/near-multiple-root states;
- NaCl inventory conservation (Profile-C molality remains prescribed model input);
- pore-volume saturation, mass density, viscosity, mobility, enthalpy/internal energy;
- mass/energy fluxes, capillary pressure or relative permeability;
- conservation residual assembly, mesh/discretization, time integration or a nonlinear solver;
- CPA closure.

Those are separate model/numerical increments and must not be inferred from a valid local thermodynamic linearization.
