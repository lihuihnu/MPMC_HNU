# Physics thermodynamic-closure consumption contract

## Scope

The `physics` layer consumes already accepted thermodynamic/flash states and exposes owned thermodynamic snapshots for future conservation-law and discretization code. It does **not** create mesh, flux, saturation, Darcy flow, time stepping, global residual assembly, a nonlinear solver, viscosity, enthalpy or an energy equation.

Two adapters are currently implemented:

- `mpmc/physics/pr76_thermodynamic_closure.hpp`: accepted PR76 internal two-phase VLE primal plus its validated local linearization;
- `mpmc/physics/sw92_thermodynamic_closure.hpp`: accepted SW92 Profile-C authoritative 1/2/3-phase PT **primal only**.

The generic header `mpmc/physics/thermodynamic_closure.hpp` therefore contains two separate payload shapes:

- `ThermodynamicClosureSnapshot`: the original fixed liquid/vapor PR76 VLE contract;
- `PtPhaseSetThermodynamicClosureSnapshot`: a variable-cardinality phase-set primal contract that does not invent model-specific phase labels.

The CMake target remains `mpmc::physics`. Dependency direction is `physics -> flash/thermodynamics` (and, for the existing PR76 derivative adapter, `flash_sensitivity -> ad`); no lower layer depends on physics.

## Atomic snapshot and two-axis availability

Primal availability and linearization availability are independent decisions.

A valid primal may be used by a future conservation residual or line-search merit evaluation. A Newton/Jacobian base point may consume a derivative only when the same closure snapshot explicitly publishes a supported linearization.

No adapter may replace an unavailable derivative with zeros, a previous-cell or previous-iterate derivative, clipping, or a hidden finite-difference fallback.

### PR76

For the existing PR76 two-phase adapter, an accepted two-phase flash can remain a valid primal when sensitivity is unavailable because of a phase-boundary guard, local ill-conditioning, unsupported feed support, derivative property failure or derivative arithmetic failure. In those cases `residual_available()==true`, `can_seed_newton()==false`, and `linearization` is empty.

If the sensitivity path reports `solution_not_accepted`, or the flash itself is unstable/indeterminate, the adapter does not publish a usable primal. Stable PR76 single-phase closure remains outside that v1 adapter.

### SW92 Profile-C

The SW92 physics adapter intentionally still publishes a **usable 1/2/3-phase primal only**:

```text
primal_status        = valid
linearization_status = unavailable
linearization_reason = not_implemented
residual_available() = true
can_seed_newton()    = false
```

A separately validated flash-layer Jacobian is now available through [`differentiate_sw92_profile_c_phase_set(...)`](../flash/sw92_profile_c_sensitivity.md) for a fixed authoritative phase set, fixed AQ/NA family assignment and fixed selected root branches. This physics adapter does **not** consume that Jacobian yet. Wiring it into an atomic physics closure snapshot is the next downstream adapter gate; until then `can_seed_newton()` must remain false. This separation prevents a newly available flash derivative from silently changing the established physics consumption contract.

## Generic variable-cardinality phase-set primal

`PtPhaseSetThermodynamicClosureSnapshot` owns:

- pressure `p` and temperature `T`;
- ordered overall feed and ordered component IDs;
- thermodynamic model / dataset / revision identity;
- a vector of `ThermodynamicPhaseState` objects.

Each phase state contains:

```text
mole phase fraction beta_alpha
ordered mole composition x_alpha
compressibility factor Z_alpha
molar density c_alpha = p/(Z_alpha R T)
```

The phase fraction is a **mole phase fraction, not pore-volume saturation**.

The generic vector deliberately has no `liquid`, `vapor`, `aqueous`, AQ or NA label. Model-specific role/family identity is carried by the model adapter sidecar so generic physics code cannot silently turn a numerical phase slot into physical morphology.

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
- per-phase physical role and AQ/NA family;
- per-phase selected cubic-root branch as numerical provenance.

The SW sidecar retains only the roles already authorized upstream:

- `aqueous` with AQ family;
- `nonaqueous_unclassified` with NA family.

H0/H1 remain `nonaqueous_unclassified`; branch index or Z ordering is not promoted to liquid/vapor or LV/LL morphology. `morphology_resolved=false` and `global_stability_proven=false` remain explicit.

### Property reproduction and density

Before publishing the physics primal, the adapter evaluates the existing SW92 phase-property kernel at each **exact accepted** composition, family and root branch and requires the reproduced `ln(phi)` to match the activity retained by the authoritative flash result within a roundoff guard.

For accepted 2/3-phase states whose authoritative publication already owns Z, the reproduced Z must agree and the physics state preserves the exact accepted Z.

The authoritative dry no-W single-H publication intentionally stores reference activity but no selected Z, because the publication adapter was required to remain a pure projection. The physics closure has a different responsibility: it needs a molar volume/density. For that one case it obtains Z from the already identified family/root branch during the same property-reproduction call. It still does **not** rerun phase selection or flash logic.

Molar density introduces no new empirical model:

```text
Z = p v_m/(R T)
c = 1/v_m = p/(Z R T)   [mol/m^3]
```

No mass density, saturation, mobility or transport property is implied.

### Failure policy

A non-accepted authoritative status, malformed phase set, inconsistent publication convention, invalid role/family pair, phase-fraction/composition shape error, or unreproducible activity/Z yields an `indeterminate` closure with no residual primal.

Passing a different ordered SW92 model snapshot (dataset/revision/component order mismatch) is an input-contract error and throws `std::invalid_argument`, consistent with the existing PR76 closure adapter's model/snapshot mismatch behavior.

## PR76 coordinates and local linearization

The existing PR76 derivative columns remain:

```text
q = (p_Pa, T_K, z_0, ..., z_{N-2})
z_{N-1} = 1 - sum(z_0, ..., z_{N-2})
```

The PR76 primal is anchored to the exact accepted flash `beta/x/y/Z` snapshot. Before publication, the adapter re-evaluates the selected branches with the accepted full phase compositions and requires reproduced Z agreement.

For an available linearization, the adapter obtains local selected-branch Z partials from the existing AD/simple-root IFT phase kernel and combines them with converged flash composition derivatives:

```text
dZ_alpha/dq = partial Z_alpha/partial(p,T,w_alpha) * d(p,T,w_alpha)/dq
```

Molar-density derivatives use:

```text
dc_alpha/c_alpha = dp/p - dZ_alpha/Z_alpha - dT/T
```

No flash iteration, RR bisection, line search, TPD search or root iteration is differentiated by the closure adapter itself.

## Validation

The SW92 primal-consumption regression covers:

1. traceable dry CO2/H2O authoritative single-H state, including reconstruction of the source-missing Z;
2. traceable wet CO2/H2O authoritative W+H state with exact accepted phase-fraction/composition/Z preservation;
3. physical Mortezazadeh–Rasaei Sample-6 authoritative W+H0+H1 state;
4. `c=p/(ZRT)` for every published phase;
5. AQ/NA and `aqueous` / `nonaqueous_unclassified` metadata preservation;
6. publication-convention and phase-activity tamper rejection;
7. ordered model/source mismatch rejection;
8. runtime component permutation;
9. public-header self containment;
10. valid-primal / unavailable-physics-linearization policy.

The separate SW92 flash-sensitivity regression validates the fixed-phase-set implicit Jacobian with forward AD, fresh re-solves, the physical Sample-6 three-phase state, permutation/H-slot invariance and derivative boundary guards; see the [SW92 Profile-C sensitivity contract](../flash/sw92_profile_c_sensitivity.md).

A dedicated GitHub-hosted GCC Debug+ASan/UBSan / Clang Release / MSVC Release workflow regenerates the existing independent Sample-6 Decimal(80) oracle before the focused C++ suite. Changes to the shared generic closure header also rerun the existing PR76 closure workflow and its independent reference regeneration.

## Current non-capabilities

The physics layer still does not provide:

- consumption/publication of the SW92 flash Jacobian in `Sw92ProfileCThermodynamicClosureSnapshot`;
- a validated SW92 H0/H1 liquid/vapor or LV/LL morphology resolver;
- NaCl inventory conservation (Profile-C molality remains prescribed model input);
- pore-volume saturation, mass density, viscosity, mobility, enthalpy/internal energy;
- mass/energy fluxes, capillary pressure or relative permeability;
- conservation residual assembly, mesh/discretization, time integration or nonlinear solvers;
- CPA closure.

Those are separate model/numerical increments and must not be inferred from a valid SW92 thermodynamic primal or a standalone flash sensitivity result.
