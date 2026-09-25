# Three-phase multicomponent porous-media flow

## 1. Scope and goal

`mpmc::flow` is the orchestration layer for the first production porous-media flow model in MPMC_HNU. Its target is a **fully implicit, non-isothermal, up-to-three-fluid-phase, multicomponent compositional model** for subsurface porous-media flow.

The model combines:

- component-wise molar conservation;
- local thermodynamic phase equilibrium written in a **natural-variable formulation**;
- one local-thermal-equilibrium energy conservation equation;
- multiphase Darcy velocity with gravity;
- optional capillary pressure;
- configurable relative-permeability constitutive laws;
- phase appearance/disappearance with a variable active phase set;
- a first well model based on Peaceman well-block pressure / well index.

The fluid phases are represented by generic phase slots. The core flow equations do **not** assume that the three phases are necessarily named oil/gas/water. A configured thermodynamic backend may provide physically meaningful phase-family or role metadata, but numerical phase ordering, compressibility-factor ordering, or a root index must never be promoted to a physical phase identity.

All configured components are allowed to partition among all active fluid phases **when the selected thermodynamic model and parameter set support that mutual solubility**. No component is silently made immiscible and no fixed “water phase” is appended after a hydrocarbon flash.

This PR does not add chemical reactions, mineral precipitation/dissolution, solid deformation/geomechanics, poroelasticity, fractures, adsorption, kinetic interphase mass transfer, non-Darcy flow, or compositional diffusion/dispersion. Those require separate model contracts and evidence.

## 2. Pre-implementation audit

Baseline: `main@256b24edffb5b66433f5dcc6984d61072ac11328` after merge of PR #117.

Existing reusable capabilities:

- `thermodynamics / flash`: PR76, SW92 and CPA PT phase-equilibrium paths, including existing one/two/three-phase phase-set machinery and stability logic within their documented validity;
- `physics`: model-neutral variable-cardinality PT phase-set closure snapshots, phase molar fractions/compositions, `Z`, molar density and local reduced-feed Jacobians for already supported backends; component inventory is per total fluid volume, not pore volume;
- `mesh`: topology/geometry/fields/DoF/partition plus scalar/full-symmetric permeability and conductivity property contracts;
- `discretization`: TPFA admissibility and transmissibility coefficients;
- `mesh_petsc / discretization_petsc`: DMPlex/Section/SF migration, owned connection schedules, sparsity and MPIAIJ symbolic preallocation;
- `ad`: reusable forward automatic differentiation.

Missing capabilities that must be added before a scientifically complete flow residual can be claimed:

- pore-volume accumulation and saturation;
- phase mass density on a common explicit basis;
- phase viscosity;
- phase enthalpy and internal energy;
- relative permeability and capillary-pressure constitutive contracts;
- Darcy phase/component/energy flux;
- accumulation + flux + well source residuals;
- natural-variable equilibrium rows inside the global nonlinear system;
- time integration and global nonlinear solve;
- well controls and Peaceman source terms.

No existing PT flash result may be relabeled as a completed flow solution. No missing viscosity, enthalpy, capillary, relative-permeability or well property may be invented, defaulted to a convenient physical value, or inferred from a phase name.

## 3. Physical assumptions

The first production baseline freezes the following assumptions.

1. **Continuum porous medium.** A representative elementary volume exists and Darcy-scale conservation laws apply.
2. **Rigid geometry.** The rock skeleton does not deform and no mechanical equilibrium equation is solved. A stationary rock thermal-energy term is still required by the non-isothermal energy balance; this is thermal storage, not geomechanical coupling.
3. **Local thermal equilibrium.** All fluid phases and the stationary rock share one cell temperature `T`.
4. **Local phase equilibrium.** Active fluid phases satisfy equality of component fugacity/chemical potential at the local phase pressures and common temperature.
5. **Up to three fluid phases.** The valid active phase set can contain one, two or three phases. Three phases are not forced when the equilibrium state is single- or two-phase.
6. **Mutual solubility is model-driven.** Every configured component may occur in every active phase if permitted by the chosen EOS/phase model and parameter set.
7. **Darcy flow.** Inertial/Forchheimer effects are excluded.
8. **No reactions.** Component source terms arise from wells/boundaries only in this PR.
9. **No molecular diffusion or mechanical dispersion** in the first baseline. Component transport is advective. Adding diffusion/dispersion later requires a separate flux and entropy/energy audit.
10. **No capillary hysteresis in the first baseline.** Capillary pressure may be enabled, but hysteresis/scanning curves require a later stateful constitutive contract.

## 4. Conserved quantities and governing equations

Let:

- `i = 0..Nc-1` denote components;
- `alpha in A` denote the current active fluid phases, with `1 <= |A| <= 3`;
- `phi` be porosity;
- `S_alpha` be pore-volume saturation;
- `c_alpha [mol/m3]` be phase molar density;
- `x_alpha,i` be phase mole fraction;
- `rho_alpha [kg/m3]` be phase mass density;
- `mu_alpha [Pa s]` be phase viscosity;
- `u_alpha [J/kg]` and `h_alpha [J/kg]` be phase specific internal energy and enthalpy;
- `K [m2]` be the absolute permeability tensor;
- `kr_alpha` be relative permeability;
- `v_alpha [m/s]` be Darcy volumetric flux;
- `lambda_eff [W/(m K)]` be effective porous-medium thermal conductivity.

### 4.1 Component molar conservation

For every component `i`,

```text
d/dt [ phi * sum_alpha( S_alpha * c_alpha * x_alpha,i ) ]
+ div [ sum_alpha( c_alpha * x_alpha,i * v_alpha ) ]
= q_i
```

where `q_i` is the external component molar source per bulk volume.

The accumulation basis is therefore **moles per bulk porous-medium volume**. It must not be confused with the existing physics `ComponentInventory`, which is expressed per total fluid volume.

### 4.2 Multiphase Darcy velocity

For each active phase,

```text
v_alpha =
  - K * (kr_alpha / mu_alpha)
      * ( grad(p_alpha) - rho_alpha * g )
```

with one consistent sign convention for the gravity vector `g`.

The model must distinguish:

- absolute permeability `K`;
- dimensionless relative permeability `kr_alpha`;
- dynamic viscosity `mu_alpha`;
- phase mass density `rho_alpha`.

Molar density may not be substituted for mass density in the gravity term.

### 4.3 Saturation and capillary closure

```text
sum_alpha S_alpha = 1
S_alpha >= 0
```

One active phase is selected as the explicit pressure reference. When capillary pressure is enabled,

```text
p_alpha = p_ref + pc_alpha(S, constitutive_state)
```

with `pc_ref = 0`.

When capillary pressure is disabled, all active phase pressures equal `p_ref`.

Capillary models must declare their reference phase, phase-pair meaning, pressure sign convention, input saturation definition, endpoint/residual saturations, units and validity domain. A water/oil/gas-specific correlation must not be applied to generic unclassified phase slots without an explicit physical-role mapping.

### 4.4 Local thermodynamic equilibrium

For each component `i` and every non-reference active phase `alpha`,

```text
f_i,alpha(p_alpha, T, x_alpha)
=
f_i,ref(p_ref, T, x_ref)
```

with

```text
sum_i x_alpha,i = 1
x_alpha,i >= 0
```

for every active phase.

The equality is evaluated through the selected thermodynamic model at the **actual phase pressure**. If capillary pressure is enabled, equilibrium must not silently evaluate every phase at the same pressure.

Zero/trace-component handling, logarithmic fugacity residuals, active positive support and phase disappearance require explicit numerical contracts. A zero composition must not be repaired by an undocumented epsilon.

### 4.5 Energy conservation

The first non-isothermal baseline uses one local-thermal-equilibrium energy equation. Neglecting kinetic energy and mechanical deformation, a consistent internal-energy / enthalpy form is required:

```text
d/dt [
    phi * sum_alpha( S_alpha * rho_alpha * u_alpha )
  + (1 - phi) * rho_r * u_r(T)
]
+ div [
    sum_alpha( rho_alpha * h_alpha * v_alpha )
  - lambda_eff * grad(T)
]
= q_E
```

where `rho_r * u_r(T)` is stationary-rock thermal energy per solid volume and `q_E` is the external heat/energy source per bulk volume.

The exact gravitational-potential-energy convention must be frozen before the first executable energy residual so that gravitational work is neither omitted inconsistently nor double counted. The same convention must be used by well/source energy terms.

No constant heat capacity, rock density, thermal conductivity, enthalpy or internal energy is supplied implicitly. Synthetic constants are allowed only in explicitly synthetic verification fixtures.

## 5. Natural-variable nonlinear formulation

For a fixed active phase set of size `P` and `Nc` components, the non-isothermal natural-variable state uses independent unknowns:

```text
p_ref
T
S_0 ... S_(P-2)
x_0,0 ... x_0,Nc-2
...
x_(P-1),0 ... x_(P-1),Nc-2
```

with the final saturation and final component fraction reconstructed from the exact sum constraints.

Unknown count:

```text
2 + (P - 1) + P * (Nc - 1) = P * Nc + 1
```

Equation count:

```text
Nc component balances
+ 1 energy balance
+ (P - 1) * Nc fugacity-equality constraints
= P * Nc + 1
```

Thus the local equilibrium constraints and transport conservation equations form one square fully implicit nonlinear system.

This is the core architectural requirement of this PR:

- phase saturations and phase compositions are global nonlinear unknowns;
- equilibrium equations are explicit nonlinear residual rows;
- a standalone PT flash may provide initialization, phase-stability evidence and phase-set transition candidates;
- a standalone PT flash must **not** replace the equilibrium rows while the code is described as a natural-variable formulation.

For an isothermal verification problem, remove `T` and the energy row consistently; do not retain an unused temperature unknown.

## 6. Phase appearance, disappearance and identity

The solver must support `1 <-> 2 <-> 3` phase-set changes without forcing three positive saturations everywhere.

The first implementation may validate an interior fixed-three-phase positive-support state before phase switching is enabled. However, production acceptance of this PR ultimately requires:

- explicit phase-set transition rules;
- stability re-check using the selected thermodynamic backend;
- conservative state transfer across phase appearance/disappearance;
- no stale phase composition or saturation carried into an inactive phase;
- no assignment of oil/gas/water identity based only on slot index, `Z`, density or EOS root order.

Existing model-specific family/role metadata may be preserved as sidecar provenance.

## 7. Relative permeability contract

Relative permeability is a configurable constitutive model owned above the mesh geometry layer.

Every model must declare:

- phase-role requirements, if any;
- saturation coordinates and normalization;
- residual/end-point saturations;
- parameter units and bounds;
- whether it is two-phase, three-phase, pairwise-combined or tabulated;
- differentiability and endpoint behavior;
- hysteresis support status;
- provenance/reference.

Allowed implementation directions include reference-backed Corey/Brooks-Corey-type two-phase laws, Stone I/II-type three-phase constructions and explicit tabulated curves.

Stone-type models are **not generic phase-slot models**: their oil/water/gas and wettability assumptions must be satisfied by an explicit configured role map before they can be selected.

A tabulated model must define interpolation/extrapolation and derivative behavior; silent clipping is forbidden.

## 8. Capillary-pressure contract

Capillary pressure is optional and disabled by an explicit `none` model, not by absent/guessed data.

Reference-backed options may include Brooks-Corey- or van-Genuchten-type laws where their porous-medium assumptions are appropriate. Petroleum reservoir table inputs may also be supported later.

Each law must define:

- wetting/non-wetting or explicit phase-pair roles;
- `pc` sign convention;
- pressure units;
- effective saturation definition;
- residual saturations and endpoints;
- derivative behavior near endpoints;
- whether hysteresis is unsupported or explicitly modeled.

## 9. Transport and thermal property boundary

The current thermodynamic closure does not provide all properties needed by the flow equations.

Before the first production Darcy/energy residual, a model-neutral phase-property contract must provide, for the exact natural-variable phase state:

- molar density;
- mass density;
- viscosity;
- enthalpy;
- internal energy;
- the derivatives actually required by the Newton linearization;
- model/dataset/revision/component-order provenance;
- explicit failure status when a requested property is unsupported.

A flow layer must not infer viscosity from density, infer enthalpy from `Z`, or substitute a constant property unless the selected configured property model explicitly defines that constant.

Existing `mesh::conductivity` is only a typed material-property storage contract. The flow energy model must additionally declare the thermal meaning and require a compatible conductivity unit such as `W/(m*K)`; it must not reinterpret electrical conductivity as thermal conductivity.

## 10. Spatial and temporal discretization

The first production discretization is cell-centered finite volume with fully implicit backward Euler time integration.

The existing TPFA path may be used only when its geometry/K admissibility contract allows materialization. A face requiring non-orthogonal treatment must not be silently forced through the direct TPFA coefficient.

Initial acceptance therefore allows:

- supported grids/faces for which the existing TPFA admissibility gate succeeds;
- explicit unsupported/failure status outside that domain.

Adding MPFA or another non-orthogonal discretization is a separate later increment.

The residual/Jacobian layer consumes authoritative owned connection rows from the existing discretization/PETSc ownership chain. It must not rebuild mesh ownership or invent a second stable-cell numbering system.

## 11. Wells: initial Peaceman scope

The first well model is based on Peaceman well-block pressure / well-index theory.

Initial scope:

- cell-connected wells;
- BHP control and total/component rate controls as separately audited control modes;
- phase/component source terms derived from the same phase mobility, density and composition conventions as the reservoir residual;
- energy carried by injected/produced fluid when the non-isothermal model is active;
- skin and well radius only when explicitly configured;
- hydrostatic well-depth correction with a documented sign convention.

The original Peaceman result and the nonsquare/anisotropic extension have restricted grid and continuum assumptions. A Peaceman formula must not be generalized to arbitrary unstructured cells by analogy. Unsupported geometry must fail or wait for a separately validated well-index model.

This PR does not implement multi-segment wells, wellbore storage, tubing hydraulics, fractures, inflow-control devices or facility networks.

## 12. Nonlinear solve and Jacobian requirements

The target solver is fully implicit Newton-based.

Requirements:

- residual and Jacobian use the same physical state and phase set;
- analytic/AD derivatives are preferred; hidden production finite differences are forbidden;
- local thermodynamic Jacobians must be evaluated on the selected phase/root/family branch;
- phase-boundary and ill-conditioned local derivatives remain explicit failure/transition states;
- damping/line search must not repair negative saturation/composition by unreported clipping;
- component conservation and energy conservation are checked independently after convergence;
- a small nonlinear residual alone is insufficient if phase stability or physical-domain checks fail.

PETSc may be used in the solver adapter, but PETSc types must not enter thermodynamics, flash, mesh core or model-neutral flow-domain contracts.

## 13. Verification plan and acceptance gates

The PR is complete only when the implemented scope has independent, incremental evidence for all of the following:

1. **0D natural-variable algebra:** exact equation/unknown count, fugacity constraints, saturation/composition closure and Jacobian checks on a fixed three-phase positive-support state.
2. **Single-phase reduction:** the same conservation/energy code reduces correctly when only one phase is active.
3. **Two/three-phase local transitions:** accepted phase appearance/disappearance preserves component inventory and never fabricates phase identity.
4. **Darcy sign/unit tests:** pressure, gravity, permeability, viscosity, density and transmissibility conventions are independently checked.
5. **Relative-permeability/capillary tests:** reference curves or independently calculated values, endpoint/domain failures and derivatives.
6. **Energy tests:** static conduction, advective enthalpy transport and coupled accumulation against analytical/manufactured or independent reference cases.
7. **Finite-volume conservation:** closed-domain component and energy conservation to the declared nonlinear/time-discretization tolerance.
8. **Peaceman well test:** independently reproduce a supported analytical/reference well-block relation before coupled well controls are accepted.
9. **Three-phase multicomponent end-to-end case:** use a thermodynamic/transport dataset with traceable parameters and an independent reference program or literature result. Synthetic data may validate software structure but may not be presented as physical validation.
10. **Parallel consistency:** serial and MPI/PETSc runs preserve stable ownership, global conservation and converged state within declared numerical tolerances.

No tolerance may be relaxed after seeing a failure merely to make a gate pass.

## 14. Cross-module ownership

Expected dependency direction:

```text
thermodynamics / flash
          |
          v
       physics
          |
          v
        flow  <----- constitutive transport/thermal models
          |
          +------> mesh
          |
          +------> discretization
          |
          v
     solver adapter
          |
          v
        PETSc
```

Exact CMake dependency edges will be introduced only when the first executable slice requires them.

Rules:

- thermodynamics and flash do not depend on flow;
- mesh does not own relative permeability, capillary pressure, mobility, Darcy flux or wells;
- discretization owns spatial coefficient/discretization policy, not thermodynamic phase equilibrium;
- flow owns conservation equations, phase pressure/saturation coupling, constitutive orchestration and well/source coupling;
- solver/PETSc code owns global nonlinear/linear algebra mechanics, not physical formulas.

Any required change to an existing module is audited before writing. Existing code may be refactored only when needed for a verified contract or when there is strong measured evidence of a material performance/maintenance benefit; behavior and scientific invariants must be preserved or the approved change stated explicitly.

## 15. Explicit non-goals for this PR

- chemical reactions or reactive transport;
- mineral precipitation/dissolution;
- geomechanics, poroelasticity or solid deformation;
- adsorption/desorption;
- kinetic nonequilibrium interphase transfer;
- molecular diffusion/dispersion;
- capillary hysteresis;
- non-Darcy / Forchheimer flow;
- fractures/DFM;
- multi-segment well/facility network;
- AMR/dynamic repartition;
- GPU-specific kernels;
- frontend/UI work.

## 16. Scientific and software references

Primary literature:

1. K. H. Coats, “An Equation of State Compositional Model,” *Society of Petroleum Engineers Journal*, 20(5), 363–376, 1980. DOI: https://doi.org/10.2118/8284-PA
2. D. W. Peaceman, “Interpretation of Well-Block Pressures in Numerical Reservoir Simulation,” *Society of Petroleum Engineers Journal*, 18(3), 183–194, 1978. DOI: https://doi.org/10.2118/6893-PA
3. D. W. Peaceman, “Interpretation of Well-Block Pressures in Numerical Reservoir Simulation With Nonsquare Grid Blocks and Anisotropic Permeability,” *Society of Petroleum Engineers Journal*, 23(3), 531–543, 1983. DOI: https://doi.org/10.2118/10528-PA
4. H. L. Stone, “Probability Model for Estimating Three-Phase Relative Permeability,” *Journal of Petroleum Technology*, 22(2), 214–218, 1970. DOI: https://doi.org/10.2118/2116-PA
5. H. L. Stone, “Estimation of Three-Phase Relative Permeability and Residual Oil Data,” *Journal of Canadian Petroleum Technology*, 12(4), 53–61, 1973. DOI: https://doi.org/10.2118/73-04-06
6. R. H. Brooks and A. T. Corey, *Hydraulic Properties of Porous Media*, Hydrology Paper No. 3, Colorado State University, 1964.
7. M. Th. van Genuchten, “A Closed-form Equation for Predicting the Hydraulic Conductivity of Unsaturated Soils,” *Soil Science Society of America Journal*, 44(5), 892–898, 1980. DOI: https://doi.org/10.2136/sssaj1980.03615995004400050002x
8. Z. Chen, G. Huan, and Y. Ma, *Computational Methods for Multiphase Flows in Porous Media*, SIAM, 2006. DOI: https://doi.org/10.1137/1.9780898718942

Reference programs/documentation:

- MRST compositional module, including `NaturalVariablesCompositionalModel`: https://www.sintef.no/projectweb/mrst/modules/compositional/
- MRST source repository: https://github.com/SINTEF-AppliedCompSci/MRST
- DuMuX non-isothermal porous-medium model and energy balance: https://dumux.org/docs/doxygen/master/group___n_i_model.html
- DuMuX generalized multiphase/multicomponent model documentation: https://dumux.org/docs/doxygen/master/topics.html
- OPM Flow reference manual: https://opm-project.org/?page_id=955

Reference software is used for formula/behavior comparison and independent test design. Source code is not copied into MPMC_HNU without a separate license audit.

## 17. First implementation slice

The first code slice after this scope document is intentionally smaller than the full model:

> establish a **model-neutral natural-variable cell-state and phase-property prerequisite contract** for an interior fixed-three-phase, positive-support state.

It must freeze:

- ordered component identity and generic three-phase slot identity;
- `p_ref`, `T`, two independent saturations and reconstructed third saturation;
- `Nc-1` independent compositions per phase and exact normalization;
- phase-pressure construction with optional `pc=none` baseline;
- required property payload: molar density, mass density, viscosity, enthalpy and internal energy;
- equation/unknown indexing for `Nc` component rows + one energy row + `2*Nc` fugacity rows;
- finite/domain validation and explicit unsupported-property failures.

The first slice must **not** create a mesh residual, Darcy face flux, time stepper, Newton solver, PETSc matrix values, phase switching or a well model.

Only after this local contract is independently tested should the next cross-module consumption step be selected.


## 18. First implementation slice: local natural-variable contract

The first executable slice is implemented as the header-only `mpmc::flow` target with public header:

```text
mpmc/flow/natural_variable_cell_state.hpp
```

It remains independent of thermodynamics, flash, physics, mesh, discretization, MPI and PETSc.

The fixed interior contract freezes:

- generic `phase0/phase1/phase2` numerical slots with no inferred morphology;
- `phase0` as the local fugacity-reference slot for deterministic equation indexing;
- ordered unique component IDs with `Nc >= 2`;
- strictly positive finite `p_ref` and `T`;
- two independent strictly positive saturations and dependent `S2 = 1-S0-S1 > 0`;
- `Nc-1` strictly positive independent mole fractions per phase and a dependent final fraction `x_last = 1-sum(x_independent) > 0`;
- `pc=none` only: all three local phase pressures equal `p_ref`;
- mandatory per-phase molar density, mass density and dynamic viscosity, all finite and strictly positive;
- mandatory finite specific enthalpy and specific internal energy; these may be negative because their zero depends on the thermodynamic reference convention;
- exact deterministic unknown and equation indexing for `3*Nc+1` local unknowns/equations.

Missing required phase properties are represented explicitly by absent input optionals and are rejected as unsupported prerequisites. The validated state contains no optional property values and performs no density/viscosity/energy inference.

Positive support is a strict domain contract for this slice. Zero saturations, zero component fractions, hidden epsilon insertion, clipping and silent renormalization are all rejected rather than interpreted as phase appearance/disappearance.

The dedicated `tests/flow/core` regression owns:

- unknown/equation count and index mapping;
- dependent saturation/composition reconstruction;
- ordered component identity and `pc=none` phase pressures;
- complete validated property payload;
- duplicate IDs, non-finite/non-positive p/T, boundary/invalid simplex inputs, malformed composition shape, missing property and invalid property rejection;
- public-header self containment.

`.github/workflows/flow_core.yml` is the required GitHub-hosted GCC Debug+ASan/UBSan, Clang Release and MSVC Release gate for this contract.

That first slice intentionally did **not** evaluate fugacity or create equilibrium residual values. Those capabilities begin in the next local thermodynamic-equilibrium slice; pore-volume accumulation, mesh/face fluxes, global Jacobian/PETSc matrix assembly, phase switching, time advancement and wells remain outside both slices.


## 19. Flow -> thermodynamics audit and fixed-three-phase fugacity residual

### 19.1 Cross-module audit

The EoS layer now exposes a common selected-phase fugacity façade in `mpmc/thermodynamics/selected_phase_fugacity.hpp`. It is deliberately a fixed-branch property interface, not a phase/root selector.

- **PR76:** selected algebraic root + exact `p,T,x` -> scalar-generic `ln_phi`; derivative capability is `scalar_generic_first_order`.
- **SW92:** selected family/molality/root + exact `p,T,x` -> scalar-generic `ln_phi`; derivative capability is `scalar_generic_first_order`.
- **CPA:** selected PT density-root index + exact `p,T,x` -> scalar-preserving `ln_phi`; derivative capability is now `scalar_generic_first_order`. The CPA path differentiates association site fractions and the selected density root by explicit IFTs and uses second-directional residual-Helmholtz derivatives; near-multiple/ill-conditioned branches remain explicit failures.

The façade never chooses a different root, family or phase on behalf of flow. CPA near-multiple/tangent topology is rejected, and an invalid selected root index is an error. Flow still consumes the model-neutral evaluator boundary below; a concrete flow-to-EoS adapter can wrap this thermodynamics façade without moving branch selection into the residual.

The evaluator contract is:

```text
(slot, actual_phase_pressure, T, full_positive_composition)
    -> ln(phi_i), i=0..Nc-1
```

with the **same scalar type** on inputs and outputs. The adapter may capture model-specific root/family provenance, but the flow residual performs no root search, stability search, family selection, fallback or finite-difference differentiation.

### 19.2 Residual definition

For `phase0` as the local fugacity-reference slot and each `alpha in {phase1, phase2}`, the positive-support equilibrium row for component `i` is

```text
R_(alpha,i) =
    log(x_0,i / x_alpha,i)
  + ln(phi_0,i) - ln(phi_alpha,i)
  + log(p_0 / p_alpha)
```

which is

```text
R_(alpha,i) = log(f_0,i / f_alpha,i)
```

when `f_i = x_i phi_i p`.

The ratio form is deliberate: every logarithm argument is dimensionless. It is equivalent to fugacity equality while avoiding a logarithm of a dimensional pressure. The reference-minus-other sign matches the repository's existing generic PT three-phase `mu0-mu1` / `mu0-mu2` chemical-potential residual convention.

Local row order is exactly:

```text
phase1 vs phase0: components 0..Nc-1
phase2 vs phase0: components 0..Nc-1
```

so the compact residual has `2*Nc` entries and maps directly to the final `2*Nc` equation rows in `NaturalVariableLayout3P`.

### 19.3 Actual phase pressure

`FugacityEquilibriumStateView3P<Number>` contains three explicit `phase_pressures_pa` values. The residual calls the thermodynamic evaluator separately with each phase's own supplied pressure.

The current `NaturalVariableCellState3P` helper still produces equal pressures because the first state slice implements `pc=none`. That helper is only a convenience adapter. The residual itself does **not** assume equal phase pressure; a future capillary-pressure layer must resolve `p_alpha` before constructing the fugacity state view.

### 19.4 Differentiability contract

`evaluate_fugacity_equilibrium_residual_3p()` is templated on the scalar `Number`. It preserves that scalar through

- actual phase pressures;
- temperature;
- full phase compositions;
- thermodynamic `ln(phi)`;
- composition-ratio logarithms;
- pressure-ratio logarithms;
- the final `2*Nc` residual values.

The dedicated regression instantiates the full residual with the existing `mpmc::ad::Dual<double,2>` type and verifies analytic pressure/composition derivative propagation, including the dependent-last composition direction. Production finite differences are not used.

The evaluator must return exactly `Nc` finite `ln(phi)` values. Model-specific thermodynamic exceptions are propagated unchanged so flow cannot silently replace a failed selected branch with another root, zeros, stale values or a fallback model.

### 19.5 Positive-support and failure policy

This slice remains interior fixed-three-phase only:

- all actual phase pressures are finite and strictly positive;
- `T` is finite and strictly positive;
- every phase has the same `Nc >= 2`;
- every `x_alpha,i` is finite and strictly positive;
- each phase composition is normalized within a machine-roundoff structural tolerance;
- no epsilon insertion, clipping or renormalization is performed.

The contract intentionally does not handle phase disappearance, zero-support logarithms or root/family switching.

### 19.6 Validation ownership

`tests/flow/core` now additionally verifies:

- exact `2*Nc` residual shape and row ordering;
- independent reconstruction of the logarithmic fugacity-ratio formula;
- thermodynamic evaluator calls receiving three distinct supplied phase pressures;
- AD propagation through actual pressure and composition coordinates;
- wrong evaluator result size, non-finite `ln(phi)`, invalid pressure/T/composition and normalization rejection;
- unchanged propagation of a thermodynamic selected-root/property failure;
- public-header self containment.

The flow core CI listens to the AD public headers because the differentiability compatibility test is an owned downstream consumer. It still has no mesh, discretization, MPI or PETSc dependency.

This slice still does **not** evaluate a concrete PR76/SW92/CPA adapter in production, assemble component/energy conservation, create Darcy flux, construct a global Newton system, insert PETSc values, switch phase sets or create wells.

## 20. Concrete thermodynamics fugacity adapters

The first concrete flow-to-EoS bridge is now owned by
`mpmc/flow/thermodynamics_fugacity_adapters.hpp`.

It exposes three fixed-branch evaluators:

- `Pr76SelectedPhaseFugacityEvaluator3P`;
- `Sw92SelectedPhaseFugacityEvaluator3P`;
- `CpaSelectedPhaseFugacityEvaluator3P`.

Each evaluator stores only a model reference and the three caller-supplied selected-phase records. A phase-slot evaluation forwards the exact supplied `p,T,x` to the thermodynamics-owned selected-phase fugacity façade. The adapters do not perform root searches, family selection, Gibbs ranking, density/Z sorting or phase-identity inference.

The base `mpmc::flow` target remains thermodynamics-independent. Concrete
EoS consumption is isolated in the `mpmc::flow_thermodynamics` interface target,
which depends one-way on `mpmc::flow + mpmc::thermodynamics`; thermodynamics and
flash remain independent of flow.

### Accepted-state integration evidence

The adapter boundary is validated inside the existing owners of the accepted three-phase fixtures rather than by copying fixture data into flow tests:

- **PR76:** Li–Firoozabadi 2012 sour-gas literature benchmark after accepted three-phase max3 closure and final stability review;
- **SW92:** authoritative Profile-C Sample-6 three-phase publication, preserving AQ/NA family metadata and the selected root branch carried by each published phase. The regression now keeps the authoritative component order; the natural-variable pivot selects a well-conditioned dependent component per phase instead of reconstructing the ~7.45e-12 aqueous trace component;
- **CPA:** the repository's accepted max3 symmetric structural fixture. It remains explicitly synthetic and is used only for equation/software validation, not physical validation.

For each state, the integration regression transfers `activity.branch` directly into the adapter selection. SW92 additionally transfers the published thermodynamic family and configured NaCl molality. No root is reconstructed from density, Z or slot ordering.

### Residual and Jacobian regression

The shared test helper builds the full non-isothermal natural-variable coordinate layout:

```text
p_ref, T, S0, S1, 3*(Nc-1) phase-composition coordinates
```

and evaluates the complete `2*Nc` fugacity-equilibrium residual with runtime forward AD. Therefore the checked Jacobian has shape

```text
(2*Nc) x (3*Nc + 1).
```

The regressions require:

- the accepted-state residual to remain within its declared equilibrium tolerance;
- every Jacobian entry to remain finite;
- the two saturation columns to be exactly zero while `pc=none`;
- pressure, temperature, one phase-1 composition column and one phase-2 composition column to agree with fixed-branch fresh central perturbations.

Finite differences are test-only cross-checks. Production residual/Jacobian evaluation remains AD/analytic through the selected EoS branch.

The Sample-6 audit exposed the numerical limitation of a fixed-last dependent
composition. That limitation is removed by the pivot contract below; no component
permutation, clipping or tolerance relaxation is required.

This slice still does not construct component or energy conservation rows, Darcy face fluxes, time-discretization terms, a global Newton system, PETSc matrix values, phase switching or wells.

## 21. Natural-variable dependent-component pivot contract

A local natural-variable chart now owns an explicit
`NaturalVariableCompositionPivot3P`. The pivot is selected from the **current full
positive-support phase compositions** before a local linearization begins.

For each phase independently:

1. choose the component with the largest mole fraction as the dependent component;
2. if two or more components are exactly tied, choose the lowest canonical component
   index;
3. freeze that choice for the lifetime of the local chart/Jacobian evaluation;
4. repivot only by explicitly constructing a new chart from a new accepted/current
   composition state.

The selection rule makes the reconstructed value

```text
x_dep = 1 - sum(x_independent)
```

the largest available composition coordinate, which avoids forcing a trace component
through subtractive reconstruction. The pivot selector requires normalized, finite,
strictly positive compositions and never clips or renormalizes them.

### Stable identity and ordering

Pivoting changes **coordinates**, not component identity.

- external component IDs and their canonical order are unchanged;
- each phase's independent composition block is ordered by canonical component index
  with only that phase's dependent component omitted;
- `independent_composition_unknown_index(phase, component)` maps a physical component
  identity to its frozen Jacobian column, returning no column only for the dependent
  component;
- `composition_unknown_identity(column)` provides the inverse phase/component
  identity for every composition column;
- fugacity-equilibrium rows remain ordered by non-reference phase and canonical
  component index, independent of the pivot;
- `fugacity_equilibrium_row_identity(row)` exposes that invariant row identity.

Therefore repivoting may change Jacobian **column coordinates**, but never silently
reorders components or residual equations. A solver that accepts a new pivot must also
accept the new explicit chart metadata; a pivot is never changed inside one residual or
Jacobian evaluation.

### Backward-compatible fixed-last input

The existing count-only `NaturalVariableLayout3P(Nc)` constructor and cell-state input
without an explicit pivot retain the original fixed-last chart for compatibility.
Production code that starts from full phase compositions should select and pass an
explicit pivot. `NaturalVariableCellStateInput3P::composition_pivot` preserves that
chart when reconstructing the validated full composition.

### Integration evidence

The PR76, SW92 and CPA accepted-three-phase adapter/Jacobian regressions now all select
their pivot directly from the accepted full phase compositions.

Most importantly, SW92 Profile-C Sample-6 is again tested in its authoritative normal
component order. Its aqueous final component is approximately `7.45e-12`, but the
pivot selects the dominant aqueous component instead, so the trace fraction remains an
independent physical coordinate rather than the subtraction remainder.

The regressions continue to require the same `2*Nc` residual closure, finite AD
Jacobian, zero saturation columns under `pc=none`, and fixed-branch central-perturbation
agreement. This pivot slice still does not add conservation, Darcy flux, time stepping,
global Newton/PETSc assembly, phase switching or wells.

## 22. Pore-volume component accumulation contract

The first local conservation-law quantity is now defined without introducing a
face flux or time-discrete residual.

For canonical component `i`:

```text
N_i = phi * sum_alpha(S_alpha * c_alpha * x_alpha,i)
```

where:

- `phi` is the dimensionless pore-volume fraction of bulk porous-medium volume;
- `S_alpha` is phase saturation, i.e. phase volume / pore volume;
- `c_alpha` is phase molar density in `mol / phase-fluid m^3`;
- `x_alpha,i` is the phase mole fraction;
- `N_i` is component accumulation in **`mol / bulk-m^3`**.

This is intentionally distinct from `physics::PtComponentInventory`, which reports
component moles per **total fluid volume** and contains no porosity or saturation.

### Current and previous snapshots

`build_pore_volume_component_accumulation()` consumes one validated
`NaturalVariableCellState3P` plus explicit porosity and returns an owned
`PoreVolumeComponentAccumulationSnapshot3P` in canonical component order.

`make_pore_volume_component_accumulation_pair()` owns the current and previous
snapshots needed by a later time-discrete balance. The pair does **not** contain a time
step, divide by `dt`, or create a residual.

Because the current PR assumes rigid porous-medium geometry, current and previous
porosity must agree to roundoff. A future pressure-dependent/poroelastic porosity model
must be introduced as an explicit constitutive extension rather than hidden in this
pair.

Current and previous component IDs/order must match exactly.

### Conservation identity

Every snapshot independently enforces

```text
sum_i N_i
  = phi * sum_alpha(S_alpha * c_alpha)
  [mol / bulk-m^3]
```

using the exact normalized phase compositions already owned by the cell state.

No mole-phase fraction is used here. `S_alpha` is the pore-volume saturation natural
variable; this is therefore not the existing total-fluid-volume inventory formula.

### Current-state Jacobian

The current accumulation Jacobian is defined on the exact frozen
`NaturalVariableLayout3P` chart:

```text
q = p_ref, T, S0, S1, pivoted phase-composition coordinates
```

with shape

```text
Nc x (3*Nc + 1).
```

For fixed porosity:

```text
dN_i/dq =
    phi * sum_alpha [
        (dS_alpha/dq) * c_alpha * x_alpha,i
      + S_alpha * (dc_alpha/dq) * x_alpha,i
      + S_alpha * c_alpha * (dx_alpha,i/dq)
    ].
```

The saturation and composition derivatives are exact chart identities:

- `S2 = 1 - S0 - S1`;
- each phase's frozen dependent composition has derivative `-1` with respect to
  every independent composition coordinate in that phase;
- canonical independent composition identity and Jacobian column mapping come directly
  from `NaturalVariableLayout3P`.

Molar density is **not** treated as constant. The caller must supply a
`PhaseMolarDensityNaturalVariableLinearization3P` containing the exact current
phase-density primal and `dc_alpha/dq` for every natural-variable column. The
accumulation layer rejects a different pivot chart, mismatched density primal, malformed
gradient shape or non-finite derivative. It never inserts zero density derivatives as a
fallback.

The Jacobian separately checks the differentiated conservation identity

```text
sum_i dN_i/dq
  = d/dq [phi * sum_alpha(S_alpha c_alpha)].
```

### Validation ownership

The dedicated `flow.accumulation.*` regression covers:

- direct `mol / bulk-m^3` evaluation from porosity, saturation, molar density and
  canonical phase composition;
- current/previous ordered component identity and rigid-medium porosity consistency;
- primal component-to-total conservation;
- a pivoted three-phase current chart;
- analytic Jacobian assembly including `dS`, `dc` and dependent-component `dx`;
- every Jacobian column against fresh perturbation of an independent synthetic density
  law;
- differentiated component-to-total conservation;
- invalid porosity, chart mismatch and non-finite density derivative rejection;
- public-header self containment.

The synthetic density law is only a structural derivative oracle. It is not physical
validation and supplies no scientific property data.

This slice still does **not** add `dt`, accumulation residual assembly, face fluxes,
Darcy velocity, mobility, gravity, capillary pressure, wells, global Newton or PETSc
value insertion.

## 23. Backward-Euler local accumulation residual contract

The first time-discrete conservation contribution is now defined strictly from the
already validated current/previous pore-volume accumulation snapshots:

```text
R_i^acc = (N_i^(n+1) - N_i^n) / dt
```

with

```text
R_i^acc [mol / (bulk-m^3 s)].
```

This remains a **local cell accumulation/time operator only**. It does not contain a
face flux, source, well, Darcy velocity, mobility, gravity or global residual assembly.

### Inputs and frozen history

`build_backward_euler_component_accumulation_residual()` consumes:

- a validated current/previous `PoreVolumeComponentAccumulationPair3P`;
- the current `PoreVolumeComponentAccumulationLinearization3P`;
- an explicit finite `dt > 0` in seconds.

The previous snapshot and `dt` are frozen inputs. They are not natural-variable
unknowns and are not differentiated.

The builder requires:

- current/previous canonical component IDs/order to match;
- current/previous porosity to match under the rigid-medium contract;
- current linearization component identity to match the current snapshot;
- current linearization porosity to match the current snapshot;
- current Jacobian shape to match its frozen `NaturalVariableLayout3P`;
- finite Jacobian entries and differentiated component closure.

### Current-state Jacobian

For the frozen current chart:

```text
d R_i^acc / d q^(n+1)
    = (1/dt) * d N_i^(n+1) / d q^(n+1).
```

No derivative of `N_i^n` appears. Changing the previous snapshot may change the
residual value but must leave the current Jacobian unchanged.

The published result retains:

- current `NaturalVariableLayout3P` and its composition pivot;
- canonical component row identity;
- `dt [s]`;
- porosity;
- per-component residual `mol / (bulk-m^3 s)`;
- total accumulation residual;
- the complete current-state Jacobian;
- the corresponding total-residual gradient.

### Conservation checks

The primal time-discrete component residuals must satisfy

```text
sum_i R_i^acc
  = [N_total^(n+1) - N_total^n] / dt.
```

For every current natural-variable column:

```text
sum_i dR_i^acc/dq
  = (1/dt) * dN_total^(n+1)/dq.
```

The builder verifies both identities and rejects non-finite scaled values. Because
backward differencing can subtract two large nearly equal inventories, the floating-point
closure check is scaled by the original current/previous inventory magnitude rather than
the already-cancelled residual magnitude. This preserves the upstream machine-roundoff
integrity contract without relaxing any physical convergence tolerance.

### Validation ownership

The independent `flow.accumulation_time.*` regression covers:

- exact backward-Euler component residual values;
- canonical component-row identity;
- exact `1/dt` scaling of every current Jacobian entry;
- total component/total accumulation closure;
- changing previous accumulation changes the residual but not the Jacobian;
- invalid/NaN/non-positive `dt`;
- current-linearization component/porosity mismatch;
- malformed Jacobian shape;
- non-finite derivative and differentiated-closure rejection;
- public-header self containment.

This slice still does **not** add face flux, Darcy velocity, source/well terms,
mobility, gravity, capillary pressure, global Newton/PETSc assembly or any nonlinear
solve.

## 24. Model-neutral three-phase saturation constitutive contract

The first Darcy prerequisites are now frozen without constructing a face flux.

Public entry:

```cpp
#include <mpmc/flow/saturation_constitutive.hpp>
```

The local positive-support saturation chart is exactly the fixed-three-phase
natural-variable chart:

```text
q_S = (S0, S1)
S2  = 1 - S0 - S1
```

with all three saturations finite and strictly positive. No clipping, endpoint
regularization or hidden effective-saturation transform is performed by the generic
wrapper.

### Relative permeability

A relative-permeability evaluator consumes the full three-phase saturation state and
returns

```text
kr_alpha  [dimensionless], alpha = phase0, phase1, phase2.
```

The model-neutral layer requires only finite nonnegative values. A stricter declared
bound such as `kr <= 1`, residual saturations, endpoint values, interpolation policy
or a particular Corey/Stone normalization belongs to the configured model and its
reference/provenance contract.

No default Corey, Brooks-Corey, Stone or tabulated model is installed by this slice.
In particular, Stone-type oil/water/gas assumptions must later be bound through an
explicit physical-role map; numerical phase slots are never promoted to those roles by
the generic wrapper.

### Optional capillary pressure

Capillary pressure is represented by offsets relative to the fixed natural-variable
reference phase `phase0`:

```text
pc_0 = 0
pc_1 = p_phase1 - p_phase0
pc_2 = p_phase2 - p_phase0
```

and resolved phase pressures are

```text
p_phase0 = p_ref
p_phase1 = p_ref + pc_1
p_phase2 = p_ref + pc_2.
```

The sign convention is therefore explicit and model-independent. Capillary offsets may
be positive or negative, but every resolved phase pressure must remain finite and
strictly positive.

Capillary pressure is disabled only by the explicit built-in
`NoCapillaryPressure3P`, which returns identically zero offsets and therefore zero
saturation derivatives. Missing/absent capillary data is not interpreted as `none`.

The first contract remains stateless: capillary hysteresis/scanning curves are not
supported here.

### Differentiable saturation interface

Both relative-permeability and capillary evaluators are scalar-generic callable
contracts. The scalar type supplied for `S0/S1` is preserved through

- reconstructed `S2`;
- `kr_alpha`;
- capillary offsets;
- resolved phase pressures.

Consequently an analytic scalar-generic constitutive model can be evaluated directly
with the existing forward-AD scalar and expose exact derivatives with respect to

```text
(S0, S1).
```

The generic layer contributes the exact chart derivative

```text
dS2/dS0 = -1
dS2/dS1 = -1.
```

It does not finite-difference a constitutive law, does not differentiate clipping and
does not switch models/roles inside one evaluation.

### Phase-role/provenance boundary

A concrete law may capture residual saturations, endpoints, saturation normalization,
wettability/phase-role mapping, tables, parameter provenance and reference revision in
its configured evaluator. If a law requires water/oil/gas or wetting/non-wetting
identity, construction of that evaluator must reject missing/incompatible role mapping
before it is passed to the model-neutral flow wrapper.

The generic contract itself only consumes numerical `phase0/phase1/phase2` slots.

### Validation ownership

The dedicated `flow.constitutive.*` regression covers:

- positive-support `S0/S1 -> S2` reconstruction;
- three dimensionless relative-permeability outputs;
- explicit capillary-offset sign convention and actual phase pressure resolution;
- forward-AD derivatives of `kr`, `pc` and resolved pressure with respect to both
  `S0` and `S1`;
- the exact zero value/zero derivative behavior of `NoCapillaryPressure3P`;
- zero/invalid saturation, negative/non-finite relative permeability, non-finite
  capillary pressure and non-positive resolved pressure rejection;
- public-header self containment;
- a reference-backed Kenyon-Behie SPE3 three-phase saturation point using the
  published/table-2 SWOF/SGOF data with explicit `phase0=oil`,
  `phase1=water`, `phase2=gas` role mapping. The regression checks source-table
  node values, interior piecewise-linear `kr` and nonzero water-oil capillary
  pressure, SI/sign conversion, exact segment slopes through forward AD, and
  explicit rejection outside the registered source brackets.

The synthetic polynomial `kr` and linear capillary law remain structural derivative
fixtures only. Physical constitutive evidence is provided separately by the SPE3
regression, which cites Kenyon & Behie (SPE 12278 / DOI 10.2118/12278-PA) and a
fixed OPM/opm-data revision of `SPE3CASE2.DATA`. The regression does not turn that
table into a new production constitutive model and does not extrapolate beyond its
registered brackets.

This slice still does **not** divide by viscosity, construct phase mobility, evaluate
gravity, create a Darcy face flux, apply transmissibility, assemble a spatial residual,
add source/well terms or insert PETSc values.

## 25. Phase transport-property linearization and local mobility contract

The local Darcy prerequisites now include exact-state phase mass density, viscosity and
mobility without constructing any face flux.

Public entry:

```cpp
#include <mpmc/flow/phase_transport.hpp>
```

### Exact-state transport-property linearization

`PhaseTransportPropertyNaturalVariableLinearization3P` binds one exact
`NaturalVariableCellState3P` to:

```text
rho_alpha [kg / m^3]
mu_alpha  [Pa s]
d rho_alpha / dq
d mu_alpha  / dq
```

for every current natural-variable column

```text
q = p_ref, T, S0, S1, pivoted phase-composition coordinates.
```

The primal values are not re-invented by this layer. They must exactly reproduce the
validated `mass_density_kg_per_m3` and `dynamic_viscosity_pa_s` already stored in the
current phase-property payload.

The transport linearization owns an exact-state identity containing the frozen layout,
ordered component IDs, p/T, three saturations and all three canonical phase
compositions. A derivative payload from another state/chart is rejected.

Mass density and viscosity each also carry explicit nonempty
`model / dataset_id / revision` provenance. A missing property/derivative is not
replaced by a constant or zero fallback.

### Saturation constitutive promotion to the full chart

The relative-permeability/capillary layer remains natively defined on

```text
(S0, S1).
```

`ThreePhaseSaturationCoordinateDerivatives3P` supplies the already analytic/AD-backed
two-column derivatives of `kr` and capillary offsets. The flow transport layer promotes
them to the full natural-variable chart:

```text
d kr_alpha / dq
d p_alpha  / dq
```

with

```text
d p_alpha / d p_ref = 1
```

for every phase, saturation derivatives inherited from `pc_alpha(S0,S1)`, and zero
T/composition columns under this v1 saturation-only constitutive contract.

The promoted payload verifies

```text
p_alpha = p_ref + pc_alpha
```

for every phase and preserves `pc_0 = 0` exactly.

### Local phase mobility

`LocalPhaseMobilityLinearization3P` combines the two validated local payloads:

```text
lambda_alpha = kr_alpha / mu_alpha
```

with

```text
lambda_alpha [1 / (Pa s)].
```

Its analytic derivative is

```text
d lambda_alpha / dq =
    [ mu_alpha * dkr_alpha/dq
      - kr_alpha * dmu_alpha/dq ]
    / mu_alpha^2.
```

The published local package contains, for each phase:

- actual phase pressure `p_alpha [Pa]`;
- mass density `rho_alpha [kg/m^3]`;
- dynamic viscosity `mu_alpha [Pa s]`;
- relative permeability `kr_alpha`;
- mobility `lambda_alpha [1/(Pa s)]`;
- the full natural-variable gradient of every one of those quantities;
- density/viscosity provenance;
- the exact current state identity/chart.

Mass density is carried now because the later gravity term requires it, but gravity is
not evaluated in this slice.

### Validation ownership

The dedicated `flow.transport.*` regression uses an explicitly synthetic
state-dependent density/viscosity law and the existing synthetic saturation constitutive
law. It checks:

- exact-state p/T/saturation/composition/pivot identity;
- required density/viscosity provenance;
- exact primal `rho` and `mu` agreement with the current cell-state payload;
- full-chart `d rho/dq` and `d mu/dq`;
- promotion of `dkr/d(S0,S1)` and `dpc/d(S0,S1)` into the complete current chart;
- `dp_alpha/dp_ref = 1`;
- exact no-capillary zero saturation derivative of phase pressure;
- `lambda=kr/mu`;
- every column of `dp`, `d rho`, `d mu`, `dkr` and `d lambda` against fresh
  central perturbations of the same synthetic model;
- state/chart/primal/provenance/non-finite derivative mismatch rejection;
- public-header self containment.

The synthetic transport law is only a structural derivative oracle. It supplies no
physical viscosity/density data and is not a reservoir-property validation.

This slice still does **not** multiply by absolute permeability or transmissibility,
evaluate `rho*g`, construct a pressure potential difference, upwind mobility, create
a Darcy face flux, assemble a spatial residual, add source/well terms or insert PETSc
values.

## 26. Two-cell / one-face phase-potential and upwind contract

The first internal-face flow-direction contract is now defined without multiplying by
static transmissibility and without creating a Darcy flux.

Public entry:

```cpp
#include <mpmc/flow/phase_potential_upwind.hpp>
```

### Geometry and direction

The potential direction is frozen as **owner -> neighbour**:

```text
d_on = x_neighbour - x_owner [m].
```

The model-neutral flow core consumes only that displacement. It does not own mesh
indices and does not recompute cell centers. A later mesh/flow adapter must use the
existing mesh owner-to-neighbour cell-center displacement from
`CellFaceGeometricOperator3D`.

`g [m/s^2]` is the physical gravitational-acceleration vector. Thus in a coordinate
system whose z-axis is positive upward, ordinary gravity has a negative z component.

The phase potential difference is

```text
DeltaPhi_alpha =
    (p_alpha,n - p_alpha,o)
  - rho_alpha,* * g . d_on
```

in Pa.

This convention is consistent with the already frozen Darcy form

```text
v_alpha = -K lambda_alpha (grad(p_alpha) - rho_alpha g).
```

A hydrostatic state therefore satisfies `DeltaPhi_alpha = 0`.

### Face phase-density policy

The v1 policy is explicit and unique:

```text
rho_alpha,* = 0.5 * (rho_alpha,o + rho_alpha,n).
```

It is represented by
`FacePhaseDensityPolicy3P::arithmetic_mean_owner_neighbour`; unsupported enum values
are rejected rather than silently mapped to another policy.

Its split two-cell Jacobian is retained:

```text
d rho_* / d q_o = 0.5 d rho_o / d q_o
d rho_* / d q_n = 0.5 d rho_n / d q_n.
```

This arithmetic face-density gravity baseline follows the standard form demonstrated
in the SINTEF/MRST reservoir-simulation material, where the discrete gravity term uses
`avg(rho)`; MRST multiphase flow also uses single-point upstream mobility weighting.

Reference:
- K.-A. Lie / SINTEF MRST book material, discrete Darcy flux with arithmetic
  `avg(rho)`: https://www.sintef.no/contentassets/8af8db2e42614f7fb94fb0c68f5bc256/mrst-book-2015.pdf
- SINTEF MRST multiphase implementation description: standard upstream mobility
  weighting.

### Potential Jacobian

Geometry and gravity are frozen face data. With arithmetic face density,

```text
d DeltaPhi_alpha / d q_o =
    - d p_alpha,o / d q_o
    - 0.5 * (d rho_alpha,o / d q_o) * (g . d_on)

d DeltaPhi_alpha / d q_n =
    + d p_alpha,n / d q_n
    - 0.5 * (d rho_alpha,n / d q_n) * (g . d_on).
```

Owner and neighbour gradients are stored separately because the two cells may have
different dependent-component pivots. This contract does not invent a face-global
natural-variable column numbering.

### Upstream selection

The future owner->neighbour phase flux is already frozen as

```text
F_alpha = -T_f * lambda_alpha,up * DeltaPhi_alpha
```

for positive materialized `T_f`, but this slice does **not** evaluate that expression.

Therefore the single-point upstream rule is:

```text
DeltaPhi_alpha < 0  -> owner upstream
DeltaPhi_alpha > 0  -> neighbour upstream
DeltaPhi_alpha = 0  -> deterministic owner tie branch
```

The exact-zero state is a non-differentiable switching surface. It is exposed as
`owner_exact_zero_tie`; selecting owner there is a deterministic branch convention,
not a claim that the upwind selector has a classical derivative at zero.

For a frozen nonzero branch, the upwind mobility Jacobian is exactly the selected
cell's local mobility Jacobian and the opposite cell block is zero. The exact-zero tie
publishes the same owner-branch Jacobian under the explicit tie status.

### Identity boundary

Owner and neighbour must use the same canonical component identity/order. Their local
composition pivots may differ. Numerical phase slots are assumed to have already been
aligned by the flow/phase-set orchestration; this face contract does not infer physical
phase identity from density, Z or root ordering.

### Validation ownership

The dedicated `flow.potential.*` regression covers:

- arithmetic face density and its split owner/neighbour Jacobians;
- fresh owner/neighbour natural-variable perturbations against the published face-density,
  phase-potential and branch-frozen upwind-mobility Jacobians for nonzero-drive phases;
- the physical-gravity-vector sign convention;
- owner->neighbour pressure and gravity potential direction;
- one phase with neighbour upstream, one with owner upstream and one exact-zero tie;
- branch-frozen upwind mobility values/Jacobians;
- owner/neighbour orientation reversal: `DeltaPhi` changes sign and the same physical
  upstream state remains selected for nonzero-drive phases;
- zero-gravity reduction to `p_n-p_o`;
- mismatched component identity, zero displacement, non-finite gravity, malformed local
  Jacobian and unsupported density-policy rejection;
- public-header self containment.

This slice still does **not** consume `T_f`, multiply transmissibility, construct
phase/component Darcy flux, choose a spatial upwind stencil beyond this one face,
assemble a residual, add source/well terms or insert PETSc values.
## 27. Materialized-TPFA internal-face phase Darcy flux contract

The first actual spatial phase-flow quantity is isolated in a cross-module adapter:

```text
modules/flow/discretization
target: mpmc::flow_discretization
```

This keeps `mpmc::flow` independent of mesh/discretization ownership. The adapter
depends on both `mpmc::flow` and `mpmc::discretization`.

Public entry:

```cpp
#include <mpmc/flow_discretization/tpfa_phase_darcy_flux.hpp>
```

### Admissibility-gate requirement

The builder consumes an existing
`TpfaInternalFaceTransmissibilityEntry3D`. It accepts the face only when all of the
following already hold:

```text
entry.disposition == materialized
entry.admissibility.disposition ==
    direct_normal_projection_k_orthogonal_candidate
entry.static_transmissibility exists
static disposition == positive_harmonic_combination
T_f > 0 and finite
face area > 0 and finite
```

A blocked/non-orthogonal/degenerate face cannot manufacture a phase flux by supplying a
number manually. Zero static transmissibility is also rejected in this first
materialized-positive contract.

The adapter does not recompute geometry, permeability, half transmissibility,
K-orthogonality or the harmonic combination.

### Phase volumetric Darcy flux

For each aligned phase slot:

```text
F_alpha =
    -T_f * lambda_alpha,up * DeltaPhi_alpha
```

with units:

```text
T_f          [m^3]
lambda       [1/(Pa s)]
DeltaPhi     [Pa]
F_alpha      [m^3/s]
```

The sign remains the owner-to-neighbour convention:

```text
F_alpha > 0  -> owner -> neighbour
F_alpha < 0  -> neighbour -> owner.
```

The input `DeltaPhi` and upwind state are not recomputed here. They come from the
already-validated `TwoCellPhasePotentialUpwindLinearization3P` contract.

### Frozen static transmissibility

Under the current rigid mesh / static permeability contract, `T_f` is frozen for this
local nonlinear derivative. Therefore:

```text
d T_f = 0
```

and

```text
dF_alpha =
  -T_f [
      lambda_alpha,up * dDeltaPhi_alpha
    + DeltaPhi_alpha * dlambda_alpha,up
  ].
```

Owner and neighbour Jacobian blocks remain separate because their natural-variable
charts/pivots may differ.

A future poromechanics/permeability-dependent transmissibility model would require a new
explicit derivative contract; it must not be smuggled into this static-TPFA path.

### Exact-zero upwind tie

If the upstream contract reports `owner_exact_zero_tie`, then

```text
DeltaPhi_alpha = 0
F_alpha = 0
```

and the published derivative is the already-declared frozen owner branch:

```text
dF_alpha =
    -T_f * lambda_alpha,owner * dDeltaPhi_alpha.
```

The `DeltaPhi * dlambda` term vanishes at the exact zero primal. This remains a branch
derivative on a non-smooth switching surface, not a claim that the global upwind map is
classically differentiable there.

### Validation ownership

The dedicated `flow_discretization.tpfa_phase_flux.*` suite checks:

- only materialized/direct/positive-harmonic/positive-`T_f` entries are accepted;
- phase volumetric flux value and owner->neighbour sign convention;
- owner and neighbour Jacobian blocks against the analytic product rule;
- fresh perturbation of `DeltaPhi` and branch-frozen `lambda_up` against the
  published flux Jacobian;
- exact-zero tie flux and owner-branch derivative;
- blocked, mismatched-admissibility, missing-`T_f`, zero-`T_f` and malformed
  phase-potential payload rejection;
- face/state identity propagation;
- public-header self containment.

This slice still does **not** multiply by `c_alpha x_alpha,i`, construct component
molar face flux, assemble owner/neighbour conservation residual rows, add source/well
terms, or insert PETSc matrix/vector values.

## 28. Upwind phase molar-content and component molar face-flux contract

The component-transport layer now consumes the already-materialized phase volumetric
flux and applies the **same per-phase upstream cell** to phase molar density and phase
composition:

```text
m_alpha,i,up = c_alpha,up * x_alpha,i,up   [mol/m^3]

n_dot_i^f =
    sum_alpha(
        m_alpha,i,up * F_alpha
    )                                      [mol/s]
```

No second, independent component upwind decision is made.

### Upwind coupling

For each phase:

- `owner_negative_phase_potential` -> owner `c_alpha, x_alpha,i`;
- `neighbour_positive_phase_potential` -> neighbour `c_alpha, x_alpha,i`;
- `owner_exact_zero_tie` -> the same deterministic owner branch already frozen by the
  phase-potential/phase-flux contract.

Thus mobility, molar density and composition all come from one consistent phase-upstream
state.

### Molar-content Jacobian

The existing `PhaseMolarDensityNaturalVariableLinearization3P` supplies
`c_alpha` and `dc_alpha/dq` on each cell's frozen natural-variable chart. Composition
derivatives come directly from that cell's dependent-component pivot:

```text
d(c_alpha x_alpha,i)
  = x_alpha,i dc_alpha
  + c_alpha dx_alpha,i.
```

For an independent phase-composition coordinate, its own component has `dx=+1` and
the frozen dependent component has `dx=-1`; all other local coordinates have the chart
derivative implied by `NaturalVariableLayout3P`.

The molar-content Jacobian exists only in the selected upstream cell block. The opposite
cell block is exactly zero.

### Component molar face-flux Jacobian

For either owner or neighbour current-state block:

```text
d n_dot_i^f =
  sum_alpha [
      (c_alpha,up x_alpha,i,up) dF_alpha
    + F_alpha d(c_alpha,up x_alpha,i,up)
  ].
```

Owner/neighbour Jacobian blocks remain separate because the cells may use different
dependent-component pivots.

At an exact-zero phase-potential tie, `F_alpha=0`, so that phase contributes

```text
d n_dot_i^f =
    (c_alpha,o x_alpha,i,o) dF_alpha
```

under the already-declared frozen owner branch.

### Component/total molar closure

Every upwind phase payload verifies

```text
sum_i c_alpha,up x_alpha,i,up
  = c_alpha,up
```

and, column by column,

```text
sum_i d(c_alpha,up x_alpha,i,up)
  = dc_alpha,up.
```

The final face result also verifies

```text
sum_i n_dot_i^f
  = sum_alpha(c_alpha,up F_alpha)
```

and the corresponding owner/neighbour differentiated closure.

The resulting canonical component rows therefore retain the same ordered component
identity as the owner/neighbour cell states.

### Validation ownership

The dedicated `flow_discretization.component_flux.*` suite covers:

- owner-, neighbour- and exact-zero-owner phase upstream selections;
- upwind `c_alpha x_alpha,i` values and selected-side-only Jacobians;
- dependent-component pivot derivatives;
- component molar flux values and component-to-total molar closure;
- owner and neighbour full-chart Jacobians against fresh reconstructed cell/density/
  phase-flux perturbations;
- exact-zero tie reusing owner molar density/composition;
- mismatched exact-state identity, molar-density primal/shape, inconsistent phase-flux
  primal and invalid upwind selection rejection;
- public-header self containment.

The density law and face-driving law in these tests are explicitly synthetic structural
fixtures; they are not physical property validation.

This slice still does **not** scatter flux into owner/neighbour conservation residual
rows, add accumulation and spatial terms together, introduce source/well terms, or
insert PETSc matrix/vector values.

## 29. Two-cell conservative component face-rate scatter contract

The component molar face flux is now mapped conservatively to the two adjacent cell
**face-rate contributions**, without combining it with accumulation.

For a component flux whose positive direction is owner -> neighbour:

```text
R_i,o^face = + n_dot_i^f
R_i,n^face = - n_dot_i^f
```

Both quantities remain in:

```text
mol/s
```

They are deliberately named face-rate contributions rather than a complete cell
residual. The existing backward-Euler accumulation contribution remains in
`mol/(bulk-m^3 s)`; the two are not added in this slice.

### Four explicit Jacobian blocks

The scatter publishes all four local row/column blocks:

```text
d R_owner^face     / d q_owner
d R_owner^face     / d q_neighbour
d R_neighbour^face / d q_owner
d R_neighbour^face / d q_neighbour
```

and freezes strict antisymmetry:

```text
d R_owner^face / d q_owner
  = + d n_dot^f / d q_owner

d R_neighbour^face / d q_owner
  = - d n_dot^f / d q_owner

d R_owner^face / d q_neighbour
  = + d n_dot^f / d q_neighbour

d R_neighbour^face / d q_neighbour
  = - d n_dot^f / d q_neighbour.
```

The same exact sign relation is retained for the total molar face-rate diagnostic.

Owner and neighbour keep their independent frozen natural-variable charts and
dependent-component pivots. No global matrix numbering is introduced here.

### Exact conservation

For every canonical component:

```text
R_i,o^face + R_i,n^face = 0
```

and for every owner/neighbour natural-variable column:

```text
d R_i,o^face + d R_i,n^face = 0.
```

The scatter uses exact sign negation of the already-validated component face flux and
its Jacobians. It does not recompute phase flux, upwind state, molar density,
composition or transmissibility.

A zero primal face rate does **not** imply a zero Jacobian. The contract preserves
nonzero derivatives even when the current face-rate value is exactly zero.

### Validation ownership

The dedicated `flow_discretization.scatter.*` regression covers:

- positive owner->neighbour and negative component face-rate signs;
- exact component-wise owner/neighbour conservation;
- all four owner/neighbour Jacobian blocks;
- exact Jacobian antisymmetry for both owner and neighbour columns;
- fresh perturbation of the source component face flux against both scattered row
  derivatives;
- zero primal face rate with nonzero Jacobian preservation;
- source component/total molar-flux closure and differentiated-closure validation;
- malformed identity/shape, non-finite derivative and inconsistent total rejection;
- public-header self containment.

This slice still does **not** introduce cell bulk volume, divide spatial face rates by
cell volume, combine face-rate and accumulation units, assemble multi-face cell
residuals, add source/well terms, assign global matrix rows/columns, or insert PETSc
values.

## 30. Rigid cell bulk-volume normalized spatial contribution contract

The conservative component face-rate scatter can now be normalized by explicit owner
and neighbour **cell bulk volumes**:

```text
V_b,o [m^3]
V_b,n [m^3]
```

The authoritative volume remains mesh geometry. Existing 3D mesh contracts already
publish cell volumes in cubic metres, for example
`CornerPointGeometry3D::cell_volume_m3()` and the linear/cartesian
`cell_volumes_m3` snapshots. The flow-discretization contract consumes explicit
positive finite values and does not recompute geometry.

### Rigid-grid derivative contract

The current v1 geometry is rigid:

```text
d V_b,o = 0
d V_b,n = 0
```

This is exposed in
`NormalizedComponentFaceContributionLinearization3D::bulk_volume_derivative_is_zero`.

A future poromechanics/deforming-grid model must introduce an explicit volume
linearization rather than silently reusing this contract.

### Normalized spatial contribution

For each canonical component:

```text
R_i,o^(face,V) = R_i,o^face / V_b,o
R_i,n^(face,V) = R_i,n^face / V_b,n
```

with units:

```text
mol / (bulk-m^3 s)
```

which now match the existing backward-Euler accumulation residual units.

Because `dV_b=0`, the four Jacobian blocks are simply:

```text
d R_i,o^(face,V) / dq = (1/V_b,o) d R_i,o^face / dq
d R_i,n^(face,V) / dq = (1/V_b,n) d R_i,n^face / dq.
```

Owner and neighbour retain their own frozen natural-variable charts/pivots.

### Weighted conservation after normalization

If cell volumes differ, the normalized spatial values are generally **not** direct
negatives. The correct invariant is volume weighted:

```text
V_b,o R_i,o^(face,V)
  + V_b,n R_i,n^(face,V)
  = 0.
```

For every owner/neighbour natural-variable column:

```text
V_b,o dR_i,o^(face,V)
  + V_b,n dR_i,n^(face,V)
  = 0.
```

The same invariant is checked for the total molar diagnostic.

The input conservative scatter is also revalidated for canonical identity, exact
owner/neighbour sign conservation, component-to-total closure and differentiated
closure before normalization.

### Validation ownership

The dedicated `flow_discretization.normalized_spatial.*` suite covers:

- unequal positive owner/neighbour bulk volumes;
- `mol/s -> mol/(bulk-m^3 s)` scaling;
- the fact that unequal volumes destroy direct normalized antisymmetry while preserving
  volume-weighted conservation;
- all four owner/neighbour Jacobian blocks under frozen `dV_b=0`;
- volume-weighted Jacobian conservation;
- fresh perturbation of the source face-rate contribution against normalized
  Jacobians;
- zero primal face rate with nonzero normalized Jacobian preservation;
- zero/non-finite volume and malformed source-scatter rejection;
- public-header self containment and the explicit rigid-volume derivative flag.

This slice still does **not** sum multiple faces into a cell, combine normalized spatial
contributions with the backward-Euler accumulation residual, add source/well terms,
assign global rows/columns, or insert PETSc matrix/vector values.

## 31. Local multi-face component conservation residual contract

The first complete **local component-conservation residual** now combines the existing
backward-Euler accumulation contribution with all already-normalized internal-face
spatial contributions incident on one cell:

```text
R_i,c =
    R_i,c^acc
  + sum_{f in internal(c)} R_i,c^(face,V)
```

Every term is in:

```text
mol / (bulk-m^3 s)
```

so no unit conversion is performed inside this builder.

Public entry:

```cpp
#include <mpmc/flow_discretization/local_component_conservation_residual.hpp>
```

### Explicit local-cell binding

The caller supplies:

- the frozen local `NaturalVariableStateIdentity3P`;
- the explicit local cell bulk volume `V_b,c [m^3]`;
- one validated `BackwardEulerComponentAccumulationResidual3P`;
- zero or more borrowed normalized internal-face contributions, each with an explicit
  statement that the local cell is the face owner or neighbour.

The side is never inferred by comparing numerical state values. Two different cells may
legitimately have identical `p/T/S/x`, so equal state payloads are not promoted to cell
identity.

Every selected incident face must match the local component order, frozen natural-variable
chart/pivot and local bulk volume. Duplicate face identities are rejected.

### Jacobian structure

The local diagonal block is assembled as:

```text
dR_c/dq_c =
    dR_c^acc/dq_c
  + sum_f dR_c^(face,V)/dq_c.
```

For every incident internal face, the opposite-cell derivative is retained as a separate
off-diagonal block:

```text
dR_c/dq_neighbour(f).
```

These blocks remain keyed by the local face identity. This slice deliberately does not
invent global cell-column numbering and does not use numerical state equality to guess a
neighbour cell. A later topology/global-assembly layer may map each face to the
authoritative neighbour and coalesce duplicate cell-pair blocks if required.

Component-to-total closure is revalidated for the accumulation input, every normalized
face payload, the final local residual, the final diagonal Jacobian and every neighbour
block. The upstream rigid-grid contract remains unchanged:

```text
dV_b = 0.
```

### Validation ownership

The dedicated `flow_discretization.local_conservation.*` regression covers:

- two incident internal faces with the local cell appearing once as owner and once as
  neighbour;
- accumulation + multi-face residual summation;
- local diagonal and per-face neighbour Jacobian blocks;
- fresh perturbation of the accumulation and face source payloads for every local and
  neighbour natural-variable column;
- zero primal residual with nonzero structural Jacobian;
- duplicate face, null binding, wrong side/state, inconsistent bulk volume, malformed
  accumulation/face closure and unsupported side rejection;
- public-header self containment.

This slice still does **not** add boundary fluxes, source/well terms, energy residuals,
fugacity-row global assembly, global row/column numbering, PETSc Mat/Vec insertion,
Newton/time-step orchestration or phase switching.

## 32. Serial closed-owned multi-cell component conservation snapshot

The first multi-cell conservation layer now consumes the authoritative connection
schedule semantics without importing PETSc into the flow/discretization public API.

Public entry:

```cpp
#include <mpmc/flow_discretization/owned_multi_cell_component_conservation.hpp>
```

The schedule adapter reads only:

```text
schedule.assembly_rows()
```

and copies each authoritative face/local+stable cell identity into a PETSc-free view.
`schedule.ghost_rows()` are not consumed.

### Strict serial v1 boundary

This v1 builder intentionally requires a one-rank `PartitionSnapshot`.

That restriction is scientific/assembly correctness, not a portability limitation:
under MPI the rank that owns an authoritative face may differ from the owner rank of one
or both endpoint cells. A face-owner rank can therefore compute both conservative
spatial sides, but it cannot claim that a remote owned cell row is complete until those
residual/Jacobian contributions have been routed to the endpoint cell owner.

The serial baseline avoids pretending that this owner-targeted exchange already exists.
Distributed residual/Jacobian exchange is a later explicit contract.

### Complete closed-patch rows

For every serial cell the snapshot publishes:

- local and stable cell identity;
- the complete existing local component residual;
- its accumulation + multi-face diagonal natural-variable Jacobian;
- off-diagonal blocks coalesced by **stable neighbour cell identity**;
- every stable face contributing to each coalesced cell-pair block.

Every authoritative internal face is matched exactly once to one normalized face
contribution and is scattered to both geometric endpoints.

The final spatial invariant is checked after restoring bulk volume:

```text
sum_c V_b,c * R_i,c^spatial = 0
```

for every component and for the total molar diagnostic.

The differentiated invariant is also checked for every source cell and every column of
that cell's frozen natural-variable chart:

```text
sum_target V_b,target
    * dR_target^spatial / dq_source = 0.
```

Thus unequal cell volumes, different dependent-component pivots and two-sided face
Jacobians are all retained without inventing global scalar numbering.

### Schedule compatibility gate

The core header is PETSc-free. The existing PETSc integration test separately compiles
the adapter against the real
`discretization_petsc::ParallelOwnedConnectionSchedule3D`, proving that its
`assembly_rows()` payload satisfies the production schedule interface.

No PETSc matrix/vector value insertion occurs in this slice.

### Validation ownership

The dedicated `flow_discretization.multi_cell_conservation.*` tests cover a
three-cell/two-authoritative-face closed patch with unequal cell bulk volumes and
different composition pivots:

- stable/local row identity;
- owner/neighbour two-sided scatter;
- middle-cell multi-face residual accumulation;
- diagonal and stable cell-pair off-diagonal Jacobian mapping;
- volume-weighted component/total spatial conservation;
- per-source-cell spatial Jacobian conservation;
- duplicate/missing authoritative face, stable-ID, volume and non-serial rejection;
- public-header self containment.

This slice still does **not** implement distributed owner-targeted residual/Jacobian
exchange, boundary/source/well terms, energy/fugacity rows, global scalar numbering,
PETSc Mat/Vec insertion, Newton iteration or phase switching.

## 33. Distributed owner-targeted component conservation exchange

The first distributed finite-volume component-row exchange is isolated in a PETSc/MPI
adapter module:

```text
modules/flow/discretization/petsc
target: mpmc::flow_discretization_petsc
```

The core `mpmc::flow_discretization` public API remains PETSc/MPI-free.

Public entry:

```cpp
#include <mpmc/flow_discretization_petsc/distributed_component_conservation.hpp>
```

### Authoritative-face evaluation and endpoint ownership

Only:

```text
ParallelOwnedConnectionSchedule3D::assembly_rows()
```

generate physical spatial contributions. A ghost schedule row never evaluates or
duplicates flux/residual physics.

For every authoritative internal face, the face-owning rank creates exactly two
endpoint-row payloads:

```text
authoritative face
    -> owner-cell row payload
    -> neighbour-cell row payload
```

Each payload is routed to the `PartitionSnapshot::owner_rank(cell)` of its target
endpoint. The face-owner rank is not assumed to own either endpoint cell.

No process-local index crosses MPI. The exchange transports stable target/column cell IDs
and stable face ID, then the receiver resolves its own local/ghost indices.

### Payload and validation

The owner-targeted payload carries:

- stable target cell ID;
- stable opposite/column cell ID;
- stable face ID;
- exact frozen-state/chart fingerprints for target and column cells;
- target bulk volume;
- normalized component and total spatial values;
- target-row diagonal face Jacobian;
- opposite-cell off-diagonal Jacobian.

State fingerprints include canonical component order, dependent-component pivots,
`p/T/S/x` and therefore reject stale or mismatched ghost charts before assembly.

The endpoint owner requires the target, column and face stable IDs to exist in its local
overlap and requires the local schedule copy to preserve the same canonical
owner/neighbour orientation.

### Complete owned rows

The endpoint owner starts from its authoritative backward-Euler accumulation row and adds
all received face-side payloads.

The result retains:

- complete owned-cell component residual;
- complete owned-cell diagonal natural-variable Jacobian;
- per-face neighbour blocks;
- off-diagonal blocks coalesced by stable neighbour cell identity;
- contributing stable face IDs for every cell-pair block.

No global scalar DoF number is assigned in this layer.

### Global conservation gate

After owner-targeted exchange, every rank contributes only its final owned rows to the
collective closed-patch audit:

```text
sum_owned_cells V_b,c * R_i,c^spatial = 0
```

for every component and for the total molar diagnostic.

The dedicated two-rank regression additionally reconstructs the final distributed owned
rows and, for each source cell and every column of its frozen natural-variable chart,
checks:

```text
sum_owned_target_cells
    V_b,target * dR_target^spatial / dq_source = 0.
```

The fixture deliberately reverses local cell ordering between the two ranks so passing
the test requires stable-ID routing rather than accidental LocalIndex agreement.

### Communication boundary

The implementation uses owner-targeted `MPI_Alltoallv` for fixed-width metadata and
numeric endpoint payloads. It does not create a PETSc matrix/vector and does not call
`MatSetValues`, `VecSetValues` or any nonlinear solver operation.

This slice still does **not** add boundary/source/well terms, energy/fugacity rows,
global scalar numbering, PETSc value insertion, Newton iteration or phase switching.

## 34. Global component-row / natural-variable scalar mapping bridge

The complete distributed owned component rows can now be converted into deterministic
scalar assembly entries without creating or modifying PETSc vectors or matrices.

Public entry:

```cpp
#include <mpmc/flow_discretization_petsc/global_component_assembly_mapping.hpp>
```

### Two numbering spaces remain explicit

The repository already has two valid but different global numbering concepts:

1. `DofNumberingSnapshot`: mesh-global scalar DoF identity ordered by global entity
   ordinals;
2. `PetscMpiAijSymbolicPreallocation3D`: PETSc rank-contiguous owned **cell** rows.

They are not assumed to be numerically equal.

For the fixed-three-phase natural-variable formulation, every cell has one fixed-width
block

```text
q = 3*Nc + 1
```

and the PETSc scalar system is expanded from the existing PETSc cell-row bridge:

```text
petsc_scalar(cell, slot)
    = petsc_cell_global_row(cell) * q + slot.
```

Because every rank owns an integer number of complete cell blocks with the same `q`,
multiplying the existing PETSc cell ownership range by `q` preserves a contiguous
rank-local scalar row range. A transient `PetscLayout` verifies this property.

### Row and column semantics

The selected `DofLayout` variable must:

- live on `EntityKind::cell`;
- have exactly `q` components;
- occupy the complete cell scalar block.

For a component conservation row:

```text
row slot =
    NaturalVariableLayout3P::component_conservation_equation_index(i)
```

and for a Jacobian column:

```text
column slot =
    frozen local/neighbor natural-variable unknown index.
```

Thus different dependent-component pivots remain local chart semantics; they do not
change the fixed scalar block width or silently reorder canonical component rows.

The current bridge emits only the first `Nc` equation rows. The later energy and
fugacity rows will occupy the remaining slots in the same square per-cell equation block.

### Assembly-ready output

The bridge publishes:

```text
residual:
    (PETSc global row, value)

Jacobian:
    (PETSc global row, PETSc global column, value)
```

with stable row/column cell identities, component-row identity, natural-variable column
index and diagonal/off-diagonal cell-block classification.

Numerically zero Jacobian values are deliberately retained. Structural sparsity comes
from the existing cell-pair pattern, not from testing whether a current coefficient is
zero.

Every output entry also keeps the independent `DofNumberingSnapshot` mesh-global DoF
index as provenance. The test deliberately reverses mesh-global cell ordinals relative
to PETSc row ownership and requires the two numbering values to differ, preventing an
accidental conflation.

### Structural and PetscSection gates

The bridge consumes the existing
`OwnedCellStructuralColumnPatternSnapshot3D` and requires the component Jacobian's
coalesced cell-pair blocks to match it exactly:

- self and same-rank neighbours -> diagonal cell columns;
- ghost/remote-owned neighbours -> off-diagonal cell columns;
- no missing or extra structural neighbour is accepted.

It also calls the existing `mesh_petsc::create_section_mapping()` to verify that the
selected local natural-variable scalar ordering and `DofNumberingSnapshot` mapping
agree with the PETSc Section representation.

No `Mat`, `Vec`, `MatSetValues`, `VecSetValues`, Newton iteration, energy/fugacity
rows, boundary/source/well terms or phase switching are introduced in this slice.

## 35. Fugacity-equilibrium local linearization and global assembly mapping

The fixed-three-phase fugacity-equilibrium rows now have a dedicated frozen
linearization carrier and a global scalar mapping bridge.

### Local carrier

Public entry:

```cpp
#include <mpmc/flow/fugacity_equilibrium_linearization.hpp>
```

`FugacityEquilibriumResidualLinearization3P` stores:

- the frozen `NaturalVariableLayout3P` chart;
- canonical component identity/order;
- the existing `2*Nc` fugacity residual values;
- an already-computed `(2*Nc) x q` Jacobian.

It does **not** evaluate an EOS, select roots/families/phases, or finite-difference the
thermodynamics. The derivative payload is supplied by the existing scalar-generic
fugacity residual path and selected-phase thermodynamic adapters.

The carrier validates shape, finiteness and local/global fugacity-row ordering before
the PETSc layer can consume it.

### Equation rows

For reference phase `phase0`, the two equilibrium blocks are:

```text
phase1 vs phase0 : Nc rows
phase2 vs phase0 : Nc rows
```

and map through the existing natural-variable equation contract:

```text
component rows : 0 ... Nc-1
energy row     : Nc
fugacity rows  : Nc+1 ... 3*Nc
```

For `Nc=3`, one cell block therefore uses:

```text
0,1,2 -> component conservation
3     -> reserved energy row
4..9  -> fugacity equilibrium
```

The energy slot remains untouched by this slice.

### Global fugacity mapping

Public entry:

```cpp
#include <mpmc/flow_discretization_petsc/fugacity_equilibrium_global_assembly_mapping.hpp>
```

The mapping reuses exactly the same:

- `DofLayout`;
- `DofNumberingSnapshot`;
- `PetscMpiAijSymbolicPreallocation3D`;
- `OwnedCellStructuralColumnPatternSnapshot3D`;
- natural-variable PETSc scalar block width `q=3*Nc+1`;

as component conservation.

Every fugacity row is local thermodynamics, so it emits only the cell-diagonal
natural-variable Jacobian block:

```text
row    -> owned cell fugacity equation slot
column -> same cell natural-variable slot 0..q-1
```

No spatial neighbour Jacobian is invented.

Different dependent-component pivots remain frozen per-cell chart semantics. The
binding must match the exact chart already carried by the distributed owned component
row for that cell.

### Numbering provenance

As with component rows, every fugacity residual/Jacobian entry stores both:

- PETSc rank-contiguous global scalar row/column;
- independent mesh-global `GlobalDofIndex` provenance.

The two numbering spaces are explicitly cross-checked but never assumed numerically
equal.

### Validation

Flow core adds dedicated carrier regressions for:

- `2*Nc` row ordering;
- pivot-aware equation indices;
- full `2*Nc x q` Jacobian indexing;
- malformed shape/component identity;
- non-finite residual/Jacobian rejection;
- public-header self containment.

The 2-rank PETSc regression additionally checks:

- fugacity rows occupy exactly slots `Nc+1 ... 3*Nc`;
- the energy row remains unoccupied;
- all `2*Nc*q` dense local Jacobian scalars map to the same cell block;
- current pc=none exact-zero saturation derivatives are retained rather than dropped;
- reversed mesh-global/PETSc numbering remains distinct;
- pivot/chart mismatch is collectively rejected;
- missing owned-cell fugacity linearization is collectively rejected.

No `MatSetValues`, `VecSetValues`, energy residual, boundary/source/well term,
Newton iteration or phase switching is introduced here.

## 36. Local non-isothermal energy conservation contract

The fixed-three-phase natural-variable model now has its first complete **local energy
equation** through accumulation plus internal-face advection/conduction.

### Bulk energy storage

Public entry:

```cpp
#include <mpmc/flow/energy_accumulation.hpp>
```

The bulk-volume internal-energy storage is:

```text
E_bulk =
    phi * sum_alpha(S_alpha * rho_alpha * u_alpha)
  + (1 - phi) * e_r
```

with units:

```text
J / bulk-m^3
```

where:

- `rho_alpha [kg/m^3]` is phase mass density;
- `u_alpha [J/kg]` is phase specific internal energy;
- `e_r [J/rock-m^3]` is explicit stationary-rock volumetric internal energy.

The rock model is intentionally caller-supplied. This layer does **not** invent a rock
density, heat capacity, reference temperature or reference internal energy.

The current rigid/static porous-medium contract freezes:

```text
d(phi) = 0
d(V_b) = 0
```

and therefore excludes geomechanical pore-volume work.

Backward Euler gives:

```text
R_E^acc =
    (E_bulk^(n+1) - E_bulk^n) / dt
```

in:

```text
W / bulk-m^3.
```

The previous storage and time step are frozen; the Jacobian contains only derivatives of
the current natural-variable state.

### Caloric derivative carrier

`PhaseCaloricPropertyNaturalVariableLinearization3P` carries exact current-state:

- specific enthalpy `h_alpha [J/kg]`;
- specific internal energy `u_alpha [J/kg]`;
- their full natural-variable gradients.

Primal values are copied from the validated `NaturalVariableCellState3P`; the carrier
therefore cannot silently mix caloric values from another thermodynamic state.

### Internal-face energy rate

Public entry:

```cpp
#include <mpmc/flow_discretization/energy_face_flux.hpp>
```

The advective energy rate reuses the already-frozen phase Darcy upwind branch:

```text
Q_adv^f =
    sum_alpha(
        rho_alpha,up
      * h_alpha,up
      * F_alpha)
```

with `F_alpha [m^3/s]`, so `Q_adv^f [W]`.

The same upstream side is used for density and enthalpy. The Jacobian differentiates both
the upstream `rho*h` payload and the existing phase volumetric flux.

Conductive heat is:

```text
Q_cond^f = G_f * (T_owner - T_neighbour)
```

with:

```text
G_f [W/K]
dG_f = 0
```

and total internal-face rate:

```text
Q_E^f = Q_adv^f + Q_cond^f.
```

Positive `Q_E^f` means owner -> neighbour.

`G_f` is supplied explicitly. Although mesh core already has scalar/full symmetric
conductivity field contracts, this slice does **not** guess a thermal interpretation or
derive a face conductance from those tensors. Conductivity-to-face discretization remains
a separate future contract.

### Bulk-volume normalization and conservative scatter

For explicit positive cell bulk volumes:

```text
R_E,o^(face,V) = +Q_E^f / V_b,o
R_E,n^(face,V) = -Q_E^f / V_b,n
```

in:

```text
W / bulk-m^3.
```

With rigid `dV_b=0`, the Jacobian uses the same fixed-volume scaling.

Unequal cell volumes do not preserve direct normalized antisymmetry. The exact physical
invariant is volume weighted:

```text
V_b,o * R_E,o^(face,V)
+ V_b,n * R_E,n^(face,V)
= 0
```

and the corresponding invariant is checked for every owner and neighbour
natural-variable Jacobian column.

### Local energy residual

Public entry:

```cpp
#include <mpmc/flow_discretization/local_energy_conservation_residual.hpp>
```

For one cell:

```text
R_E,c =
    R_E,c^acc
  + sum_f R_E,c^(face,V)
```

with a local diagonal Jacobian and one explicit neighbour Jacobian block per incident
internal face.

The current equation includes:

- fluid internal-energy storage;
- stationary-rock thermal storage;
- advective phase enthalpy transport;
- conductive heat transport.

The current v1 deliberately excludes kinetic-energy storage and gravitational potential
energy storage. Gravity can still affect the energy equation indirectly through the
already-defined phase Darcy flux used by the advective enthalpy term.

### Validation

Flow core energy tests cover:

- fluid + stationary-rock bulk storage;
- backward-Euler units and current-only Jacobian;
- fresh perturbation of every natural-variable column against the assembled storage
  Jacobian;
- zero accumulation residual with nonzero current-state Jacobian;
- snapshot closure, porosity, gradient shape and finite-value rejection;
- public-header self containment.

Flow-discretization energy tests cover:

- mixed owner/neighbour frozen upstream branches;
- `rho*h*F` advective rate;
- fixed-`G_f` conductive heat;
- fresh source perturbation of owner and neighbour face Jacobians;
- unequal-volume normalization;
- local energy residual assembly;
- two-cell closed-patch primal energy conservation;
- volume-weighted cross-cell Jacobian conservation;
- invalid conductance, volume and duplicate-face rejection;
- public-header self containment.

This slice still does **not** implement energy global/PETSc row mapping, conductivity
tensor -> thermal face conductance discretization, boundary/source/well terms,
`MatSetValues`, `VecSetValues`, Newton iteration or phase switching.

## 37. Distributed owner-targeted energy conservation exchange

The local non-isothermal energy equation now has the same owner-targeted MPI exchange
semantics as component conservation.

Public entry:

```cpp
#include <mpmc/flow_discretization_petsc/distributed_energy_conservation.hpp>
```

### Authoritative evaluation

Only:

```text
ParallelOwnedConnectionSchedule3D::assembly_rows()
```

produce physical energy-face contributions.

For every authoritative internal face, the face-owning rank starts from the already
validated normalized energy face payload and emits exactly two endpoint-row payloads:

```text
owner endpoint:
    + Q_E^f / V_b,o
    + owner-row Jacobian blocks

neighbour endpoint:
    - Q_E^f / V_b,n
    + neighbour-row Jacobian blocks
```

The target rank is always:

```text
PartitionSnapshot::owner_rank(endpoint cell)
```

and is not assumed to be the face-owner rank.

Ghost schedule rows are used only to validate the received stable face/cell orientation.
They never re-evaluate or duplicate the energy physics.

### Stable-ID communication

No process-local `LocalIndex` crosses MPI.

Each endpoint payload carries:

- stable target cell ID;
- stable opposite/column cell ID;
- stable face ID;
- target and column frozen-state/chart fingerprints;
- target cell bulk volume;
- normalized spatial energy residual contribution;
- target-cell diagonal face Jacobian;
- opposite-cell off-diagonal face Jacobian.

The receiver resolves stable IDs into its own local/ghost indices and revalidates:

- endpoint ownership;
- schedule owner/neighbour orientation;
- frozen component order/pivot/state fingerprint;
- local cell bulk volume;
- finite residual/Jacobian values.

### Complete owned energy rows

Every owned endpoint starts from its local
`BackwardEulerEnergyAccumulationResidual3P` and receives all incident authoritative
internal-face contributions.

The result publishes:

- complete owned energy residual [W/bulk-m^3];
- complete owned-cell diagonal natural-variable Jacobian;
- per-face neighbour Jacobian blocks;
- off-diagonal blocks coalesced by stable neighbour cell;
- contributing stable face IDs for each cell-pair block.

This is still a residual/Jacobian ownership contract. It does **not** assign the energy
equation to a PETSc scalar row yet.

### Distributed conservation gate

After exchange, only final owned rows contribute to the global spatial audit:

```text
sum_owned_cells
    V_b,c * R_E,c^spatial
= 0
```

in watts.

The dedicated two-rank regression additionally checks, for each source cell and every
column of its frozen natural-variable chart:

```text
sum_owned_target_cells
    V_b,target
  * d R_E,target^spatial / d q_source
= 0.
```

The two ranks deliberately use opposite local cell ordering, so passing the test requires
stable-ID routing rather than accidental LocalIndex agreement.

### Current scope

This slice reuses the already-established local energy physics:

- fluid + stationary-rock accumulation;
- upwind advective phase enthalpy;
- fixed face thermal conductance;
- rigid `dV_b=0`.

It does not derive `G_f` from conductivity tensors and does not add boundary/source/well
energy terms.

No PETSc scalar energy-row mapping, `MatSetValues`, `VecSetValues`, matrix/vector
creation, Newton iteration or phase switching is introduced here.

## 38. Energy global scalar assembly mapping

The complete distributed owned energy rows now map into the same fixed-width
natural-variable scalar block as component conservation and fugacity equilibrium.

Public entry:

```cpp
#include <mpmc/flow_discretization_petsc/energy_global_assembly_mapping.hpp>
```

### Energy equation slot

The mapping uses the existing authoritative natural-variable equation contract:

```text
component rows : 0 ... Nc-1
energy row     : Nc
fugacity rows  : Nc+1 ... 3*Nc
```

so the energy equation is always:

```text
NaturalVariableLayout3P::energy_equation_index() == Nc.
```

For `Nc=3`, energy therefore occupies scalar slot `3` inside each `q=10`
cell block.

### Same two numbering spaces

As with component and fugacity mappings, two different numbering spaces remain explicit:

1. PETSc rank-contiguous scalar row/column indices expanded from
   `PetscMpiAijSymbolicPreallocation3D`;
2. mesh-global `GlobalDofIndex` provenance from `DofNumberingSnapshot`.

The bridge uses:

```text
petsc_scalar(cell, slot)
    = petsc_cell_global_row(cell) * q + slot
```

and retains the independent mesh-global row/column identity for every emitted scalar.
The two values are cross-checked but never assumed numerically equal.

### Energy residual and Jacobian entries

For every owned energy row the bridge emits:

```text
Residual:
    (PETSc global energy row, value)

Jacobian:
    (PETSc global energy row,
     PETSc global natural-variable column,
     value)
```

The diagonal block contains all `q` columns of the owned cell. Every coalesced
neighbour energy block contributes another dense `q` columns using that neighbour's
frozen natural-variable chart.

Zero-valued Jacobian scalars are retained. Structural sparsity is defined by the
existing exact cell-pair structural column pattern, not by current coefficient values.

### Structural gates

The mapping consumes the same:

- `DofLayout`;
- `DofNumberingSnapshot`;
- `PetscMpiAijSymbolicPreallocation3D`;
- `OwnedCellStructuralColumnPatternSnapshot3D`;
- `mesh_petsc::create_section_mapping()`;

used by the existing component/fugacity assembly bridges.

The energy row's self and neighbour cell blocks must match the exact structural
cell-column pattern. Missing or extra neighbours are rejected collectively.

### Complete equation-block closure

The two-rank regression combines all three residual mappings for each owned cell and
requires the sorted PETSc-local equation slots to be exactly:

```text
0, 1, ..., 3*Nc
```

once each.

For `Nc=3`:

```text
component : 0,1,2
energy    : 3
fugacity  : 4,5,6,7,8,9
```

so the complete `3*Nc+1` natural-variable equation block now has no global row gap
and no row collision.

The regression continues to reverse local cell order and mesh-global entity ordinals
relative to PETSc ownership, verifies owned/ghost Jacobian columns, retains an
exact-zero energy Jacobian scalar, and collectively rejects wrong DoF block width or an
incomplete structural pattern.

No PETSc `Mat` or `Vec` is created or modified and no `MatSetValues` /
`VecSetValues` call is introduced in this slice.

## 39. Complete natural-variable assembly snapshot

The three equation families now combine into one validated assembly-ready snapshot before
any PETSc value insertion occurs.

Public entry:

```cpp
#include <mpmc/flow_discretization_petsc/complete_natural_variable_assembly_snapshot.hpp>
```

The snapshot consumes the already-validated global mappings for:

- component conservation;
- energy conservation;
- fugacity equilibrium.

It does not re-evaluate any physics, EOS, flux, accumulation or derivative.

### Native residual semantics are retained

The unified residual entry stores one equation-kind discriminator because the three row
families do not share one physical unit:

```text
component conservation : mol / (bulk-m^3 s)
energy conservation    : W / bulk-m^3
fugacity equilibrium   : dimensionless log residual
```

The snapshot unifies numbering and structure, not physical dimensions.

### Complete owned residual rows

For every rank, the combined residual entries must cover exactly:

```text
[petsc_scalar_row_start, petsc_scalar_row_end)
```

with one entry per owned scalar row.

Each row also stores:

- PETSc global scalar row;
- independent mesh-global `GlobalDofIndex`;
- stable row-cell identity;
- equation slot;
- equation family;
- native residual value.

The equation slot must satisfy:

```text
row % q == equation_slot
```

and must agree with the equation family:

```text
0 ... Nc-1    -> component
Nc            -> energy
Nc+1 ... 3Nc  -> fugacity
```

There are no row gaps or collisions.

### Complete Jacobian triplets

All component, energy and fugacity Jacobian entries are copied into one globally sorted
`(row,column,value)` set.

The snapshot rejects:

- duplicate `(row,column)` pairs;
- non-owned rows;
- columns outside the global scalar range;
- row/equation provenance mismatch;
- natural-variable column-slot mismatch;
- stable column cells that do not resolve through the local owned/ghost overlap.

Exact zero values are preserved.

### Symbolic compatibility is equation-aware

The existing cell-level symbolic pattern remains authoritative.

For component and energy rows, the exact scalar column set must equal all `q` columns
of:

```text
self cell
+ every structural neighbour cell.
```

For fugacity rows, the exact scalar column set is only the self-cell `q` columns,
because local thermodynamic equilibrium has no spatial neighbour derivative.

Thus the symbolic MPIAIJ cell pattern remains a safe superset for the complete system
without inventing fake fugacity neighbour coupling.

### Two-rank validation

The existing two-rank fixture continues to reverse both local cell ordering and
mesh-global entity ordinals relative to PETSc ownership.

The complete-snapshot regression checks:

- exactly `q=3*Nc+1` owned residual entries per owned cell;
- component / energy / fugacity row-family counts;
- contiguous row ownership with no gap or collision;
- globally unique `(row,column)` triplets;
- component and energy rows contain self + neighbour scalar blocks;
- fugacity rows contain only self scalar blocks;
- exact-zero Jacobian entries survive the merge;
- metadata mismatch is collectively rejected;
- a syntactically valid energy row deliberately collided with component row 0 is
  collectively rejected;
- an incomplete symbolic neighbour pattern is collectively rejected.

No PETSc `Mat` or `Vec` is created or modified and no `MatSetValues` /
`VecSetValues` call is introduced here.

## 40. PETSc complete-system materialization and value insertion

The validated complete natural-variable assembly snapshot can now be materialized into an
assembled PETSc residual vector and Jacobian matrix without introducing any solver logic.

Public entry:

```cpp
#include <mpmc/flow_discretization_petsc/complete_natural_variable_petsc_materialization.hpp>
```

### Scalar ownership

The materializer reuses the complete snapshot's scalar ownership:

```text
[petsc_scalar_row_start, petsc_scalar_row_end)
global scalar count
q = 3*Nc + 1
```

and verifies that those scalar ranges are exactly the existing cell-level PETSc ownership
ranges expanded by `q`.

### Residual vector

The materializer creates one MPI PETSc `Vec` with:

```text
local size  = owned cell count * q
global size = global cell count * q
```

and inserts only the complete snapshot's residual entries with `INSERT_VALUES`, followed
by `VecAssemblyBegin/End`.

The vector intentionally contains the three equation families in their native units. No
equation scaling, norm definition or solver convergence metric is introduced here.

### Scalar MPIAIJ preallocation

The existing cell-level symbolic counts remain the source of matrix capacity.

For every scalar row of one owned cell:

```text
d_nnz_scalar = q * d_nnz_cell
o_nnz_scalar = q * o_nnz_cell
```

These arrays are passed to `MatMPIAIJSetPreallocation()`.

The matrix is configured with:

```text
MAT_IGNORE_ZERO_ENTRIES          = false
MAT_NEW_NONZERO_ALLOCATION_ERR   = true
```

so snapshot locations cannot silently exceed the preallocated capacity.

### Equation-aware location freeze

Capacity preallocation is deliberately wider than the physical Jacobian for some rows,
especially fugacity rows.

The materializer therefore separates matrix structure from matrix values:

1. write explicit zero at every actual `(row,column)` location in the complete snapshot;
2. assemble the matrix;
3. enable `MAT_NEW_NONZERO_LOCATION_ERR = true`;
4. overwrite only those already-materialized locations with the snapshot's actual
   Jacobian values;
5. final-assemble the matrix.

Thus component and energy neighbour blocks are materialized, while a fugacity neighbour
location is **not** created merely because the cell-level MPIAIJ capacity could hold it.

Exact-zero physical Jacobian entries still obtain a real matrix location because the
structure-seeding pass uses `MAT_IGNORE_ZERO_ENTRIES=false`.

### Ownership and failure semantics

On success, the caller owns the returned `Vec` and `Mat` and must destroy them with
`VecDestroy()` and `MatDestroy()`.

On failure, both output handles remain null.

The materializer verifies:

- PETSc vector and matrix local/global sizes;
- row and column ownership ranges;
- `MATMPIAIJ` matrix type;
- snapshot/cell-bridge rank and scalar-range consistency.

It does not create a KSP or SNES object.

### Two-rank validation

The dedicated regression reads back:

- every owned residual value with `VecGetValues`;
- every explicit Jacobian triplet with `MatGetValues`.

It additionally requires:

- `MatInfo.mallocs == 0`, proving no dynamic nonzero growth was needed;
- an exact-zero snapshot triplet can be overwritten and restored after location freeze,
  proving that the zero location was materialized;
- inserting a remote-neighbour column into a fugacity row is rejected by
  `MAT_NEW_NONZERO_LOCATION_ERR`, even though the wider cell-level capacity exists;
- inconsistent scalar cell-bridge metadata is collectively rejected before output
  objects are returned.

No KSP/SNES/Newton solve, equation scaling, line search, phase switching,
boundary/source/well term or new physical model is added in this slice.

## 41. Frozen natural-variable Newton linear-system contract

The assembled natural-variable PETSc system now supports one algebraic Newton correction
without updating the physical state.

Public entry:

```cpp
#include <mpmc/flow_discretization_petsc/natural_variable_newton_linear_system.hpp>
```

The contract solves exactly one frozen linearization:

```text
J(q^k) * delta_q = -R(q^k)
```

and then stops.

### Input boundary

The solver consumes only:

- the already assembled PETSc Jacobian `Mat`;
- the already assembled residual `Vec`;
- the frozen complete-assembly numbering/provenance snapshot.

The API does not accept an EOS, flash backend, phase-property evaluator, flux operator,
accumulation evaluator or state-update callback. Therefore this layer cannot re-evaluate
thermodynamics or flow physics while solving the linear system.

### RHS ownership

The input residual is treated as read-only.

Internally:

```text
rhs = copy(R)
rhs = -rhs
```

using `VecDuplicate`, `VecCopy` and `VecScale(-1)`.

The original residual `Vec` is not modified.

### Small-system algebraic baseline

This gate intentionally uses a package-independent PETSc baseline:

```text
KSP = GMRES
PC  = NONE
restart = global scalar row count
rtol = 1e-12
atol = 1e-14
```

and rejects systems larger than 200 scalar rows.

For the current 2-cell / `Nc=3` regression, the global system has 20 scalar rows, so
GMRES retains the complete Krylov space. This is an algebraic-correctness baseline, not a
production reservoir solver or preconditioner policy.

Distributed direct LU is deliberately not frozen here because it would require selecting
an external package such as MUMPS or SuperLU_DIST, which is a separate solver-backend
decision.

### Correction provenance

The returned caller-owned `delta_q Vec` uses exactly the same PETSc scalar ownership as
the residual/Jacobian.

The solve report records every locally owned correction with:

- PETSc global scalar index;
- independent mesh-global `GlobalDofIndex`;
- stable cell identity;
- natural-variable slot;
- correction value.

Thus composition-slot interpretation remains tied to the already frozen per-cell
`NaturalVariableLayout3P` / pivot rather than being inferred from PETSc numbering.

### Linear residual audit

After `KSPSolve`, the contract independently forms:

```text
r_linear = J * delta_q + R
```

and computes:

```text
||r_linear||_2
```

The returned report includes the original residual norm, RHS norm, KSP convergence reason,
iteration count and the independently recomputed linear residual norm.

Because component, energy and fugacity equations retain different native physical units,
this combined Euclidean norm is an **unscaled algebraic verification norm only**. It is
not a physical Newton convergence criterion.

### Nonsingular manufactured regression

The earlier assembly tests intentionally use simple affine derivative fixtures whose
20-by-20 Jacobian is rank-deficient; those fixtures remain unchanged because they test
structure, ownership and numbering rather than solvability.

The linear-solve regression therefore keeps exactly the same:

- 2-rank ownership;
- `q=3*Nc+1=10` block width;
- PETSc scalar numbering;
- mesh-global provenance;
- actual Jacobian row/column locations;

but replaces only test numerical values with a strictly diagonally dominant manufactured
matrix.

A known correction `delta_q*` is chosen and the test residual is constructed as:

```text
R = -J * delta_q*
```

The regression first verifies the manufactured identity itself, then solves the assembled
PETSc system and checks:

- recovered `delta_q` against `delta_q*`;
- global/local scalar numbering;
- stable-cell and mesh-global DoF provenance;
- per-cell natural-variable slot semantics;
- unchanged input residual after the solve;
- positive PETSc convergence reason;
- independently recomputed `||J*delta_q + R||_2`;
- collective rejection of a residual Vec with incompatible size/ownership.

No state update, damping, line search, phase switching, SNES/Newton iteration, equation
scaling, well term or boundary/source term is introduced in this slice.

## 42. PETSc-owned nonlinear solve: Newton line search + GMRES + ASM

The solver adapter now delegates nonlinear and linear orchestration directly to PETSc.

Public entry:

```cpp
#include <mpmc/flow_discretization_petsc/natural_variable_snes_solver.hpp>
```

The default algorithm chain is fixed in this contract as:

```text
SNESNEWTONLS
  -> SNESLINESEARCHBT
  -> KSPGMRES
  -> PCASM
       overlap = 1
       type    = restricted ASM
```

This layer contains no custom Newton iteration loop and no custom linear solver.

### Evaluation boundary

PETSc state `x` reaches flow physics only through two explicit callbacks:

```text
function(x)  -> F(x)
jacobian(x)  -> J(x)
```

The SNES adapter itself does not evaluate EOS, flash, phase properties, accumulation,
Darcy flux or energy flux.

The function callback receives a PETSc residual `Vec`; the Jacobian callback receives an
already allocated MPIAIJ `Mat`. The adapter zeros those objects, lets the evaluator
insert the current values, and performs PETSc assembly.

The Jacobian matrix is duplicated from the existing audited natural-variable matrix
structure and keeps:

```text
MAT_NEW_NONZERO_LOCATION_ERR   = true
MAT_NEW_NONZERO_ALLOCATION_ERR = true
```

so nonlinear iterations cannot silently widen the frozen fixed-three-phase sparsity.

### PETSc 3.19 domain semantics

PETSc 3.19.6 distinguishes two cases that must not be conflated.

`SNESSetFunctionDomainError()` and `SNESSetJacobianDomainError()` are **hard**
domain-error signals: they cause the nonlinear solve to diverge with a negative SNES
convergence reason. They are therefore reserved for a state at which the authoritative
model cannot be evaluated.

Recoverable trial-step protection uses PETSc's line-search precheck instead:

```cpp
SNESLineSearchSetPreCheck(...)
```

The optional project callback receives the current state `X` and PETSc search direction
`Y`. PETSc's backtracking line search forms trial states as:

```text
X_trial = X - lambda * Y
```

The precheck may modify `Y` and mark it changed before BT evaluates any trial state.
Thus positive-support/simplex protection can remain a model-domain responsibility while
Newton and line-search orchestration remain PETSc-owned.

### Solver configuration boundary

This first audited solver contract deliberately does **not** call
`SNESSetFromOptions()`.

Unrestricted PETSc options could switch to matrix-free or finite-difference Jacobian
paths, which would violate the repository rule forbidding a production finite-difference
Jacobian fallback.

Future configurability must therefore enter through an explicit validated solver-settings
contract that whitelists supported SNES/KSP/PC choices.

### Result/provenance

The returned solution is a caller-owned PETSc `Vec`.

The solve report records, for every locally owned natural-variable scalar:

- PETSc global scalar index;
- independent mesh-global `GlobalDofIndex`;
- stable cell identity;
- natural-variable slot;
- final value.

It also records the actual PETSc solver types, nonlinear iteration count, evaluator call
counts, hard-domain-error counts, line-search precheck counts and the final native
function L2 norm.

As before, the combined component/energy/fugacity norm contains heterogeneous native
units and is only an algebraic SNES diagnostic until equation scaling is introduced.

### Two-rank nonlinear regression

The current regression retains:

- two MPI ranks;
- `Nc=3`;
- `q=3*Nc+1=10` natural-variable scalars per cell;
- existing PETSc scalar ownership;
- stable-cell and reversed mesh-global provenance;
- the same natural-variable Jacobian row/column locations.

The nonlinear values are test-only manufactured values:

```text
F_i(x)  = log(x_i / x_i*)
J_ii(x) = 1 / x_i
```

with:

```text
x0 = 10 * x*
```

The unconstrained full Newton direction would make the first trial negative. A
`SNESLineSearchSetPreCheck()` regression computes one global safe scaling for the PETSc
search direction, so the invalid trial is never passed to the function evaluator.

The regression requires:

- actual solver type `newtonls`;
- actual line-search type `bt`;
- actual outer KSP type `gmres`;
- actual PC type `asm`;
- default ASM overlap = 1;
- positive SNES convergence reason;
- at least one PETSc line-search precheck direction change;
- zero hard function/Jacobian domain errors;
- convergence to the known `x*`;
- final function norm within the audited algebraic tolerance;
- unchanged caller-owned initial-state vector;
- exact PETSc / mesh-global / stable-cell natural-variable provenance.

This is an orchestration and numbering regression, not a physical three-phase reference
case.

No phase switching, wells, boundary/source terms, equation scaling, custom ASM
subdomains or production KSP tuning are added in this slice.

## 43. Fixed-three-phase production assembly callbacks for PETSc SNES

The PETSc SNES adapter can now consume the existing fixed-three-phase flow assembly instead
of a manufactured nonlinear function.

Public entry:

```cpp
#include <mpmc/flow_discretization_petsc/fixed_three_phase_snes_assembly.hpp>
```

The assembly context owns no EOS formula and no alternative conservation equation. For
each SNES state it performs:

```text
PETSc owned x
  -> q-scaled PetscSF owned/ghost state exchange
  -> current-cell property/fugacity evaluator on each frozen local chart
  -> existing component + energy backward-Euler accumulation
  -> existing mobility / phase-potential / upwind path
  -> existing admissibility-gated TPFA phase flux
  -> existing component molar flux + conservative scatter + bulk-volume normalization
  -> existing advective enthalpy + conductive energy + bulk-volume normalization
  -> existing distributed owner-targeted component/energy conservation
  -> existing component/energy/fugacity global scalar mappings
  -> CompleteNaturalVariableAssemblySnapshot3D
```

Both the SNES function callback and Jacobian callback rebuild this chain from the current
PETSc state. No residual/Jacobian value is cached across nonlinear states and no
finite-difference path is added.

### MPI state exchange

The global SNES vector stores only owned scalar blocks. The context constructs one
q-scaled `PetscSF` from the existing cell ownership/global-row bridge and broadcasts
owned values to the local owned+depth-1-ghost cell blocks before any face evaluation.

This keeps authoritative owner/neighbour face evaluation on the existing overlap and does
not introduce a state `MPI_Allgather`.

### Frozen-chart positive-support precheck

The context supplies the SNES line-search precheck for the real natural-variable chart.
It keeps the per-cell dependent-component pivot frozen and limits PETSc's search
direction when a full trial would violate:

- positive reference pressure;
- positive temperature;
- positive `S0`, `S1`, or reconstructed `S2=1-S0-S1`;
- positive independent composition coordinates;
- positive reconstructed dependent composition in every phase.

The precheck only scales the PETSc search direction by one globally consistent safe
factor. It does not clip, renormalize, repivot, switch phase sets, or implement a second
line-search loop.

### Controlled two-rank regression

The regression uses the production accumulation, TPFA, component/energy face,
distributed-conservation and global-mapping builders. The only controlled part is the
cell closure providing simple analytic property/fugacity values and derivatives.

The controlled case sets relative permeability and thermal face conductance to zero.
Spatial component/energy flux is therefore physically zero, but the full authoritative
face/TPFA/distributed path is still executed. Previous component and energy storage are
generated at the known target with the same production accumulation builders.

The six local fugacity rows constrain the six composition coordinates, the energy row
constrains temperature, and the three component rows constrain pressure plus the two
independent saturations. The resulting per-cell 10-by-10 Jacobian is full rank; no
regularization or artificial diagonal shift is applied.

The regression verifies:

- PETSc owned-to-ghost state propagation through `PetscSF`;
- every function/Jacobian callback re-evaluates the current state;
- convergence through the existing `SNESNEWTONLS + BT + GMRES + ASM(1)` adapter;
- final state equals the known fixed-three-phase target;
- independent final production assembly has small residual;
- final production Jacobian entries are finite;
- frozen PETSc / mesh-global / stable-cell / natural-variable provenance remains intact.

This controlled closure is a software assembly/orchestration regression only. It is not
an independent physical three-phase validation case and does not replace later
PR76/SW92/CPA end-to-end validation.

No phase switching, wells, boundary/source terms, single-phase reduction, equation
scaling or custom ASM subdomains are introduced in this slice.



## 44. Variable-cardinality phase-set transition contract

The first explicit `1 <-> 2 <-> 3` transition contract is now separated from the
fixed-cardinality residual/Jacobian implementations.

A phase-count change is **not** authorized by clipping a small phase fraction or by
deleting/adding a numerical slot. The transition sequence is:

```text
current accepted natural-variable state
  -> stability / final phase-set / disappearance evidence
  -> target topology candidate
  -> fresh target-topology thermodynamic resolve/review
  -> model-neutral target phase set
  -> beta-to-saturation volume projection
  -> new target natural-variable chart
  -> fresh target natural-variable nonlinear solve
  -> final target review
```

The standalone PT flash remains a phase-set/stability candidate generator. It never
replaces the component, energy, or fugacity rows of the target natural-variable system.

### 44.1 Appearance and disappearance evidence

`PhaseSetTransitionCandidate` records the source/target active phase counts and the
reason that the current topology is being challenged. Appearance may be triggered by
initial stability, final phase-set instability, or provider topology evidence.
Disappearance requires explicit disappearance/boundary evidence.

A candidate with `target_resolve_required` cannot be projected into flow unknowns.
Projection is allowed only after the target topology is freshly resolved and published
as `target_resolved`. This preserves the existing flash-layer rule that a phase is
never added or removed merely from a source-state fraction endpoint.

The bridge from `PtFlashBackendResult` consumes only the generic accepted phase set
and `PtPhaseTransitionReport`. It does not infer oil/gas/water, liquid/vapor,
root identity, or phase family from slot order, density, or compressibility factor.

### 44.2 Mole phase fraction is not saturation

A flash candidate publishes mole phase fractions `beta_alpha`; the flow unknown is
pore-volume saturation `S_alpha`. They are not interchangeable.

For a freshly resolved target phase set with phase molar density
`c_alpha [mol / phase-m^3]`, the transition projection uses

```text
v_alpha ~ beta_alpha / c_alpha

S_alpha =
    (beta_alpha / c_alpha)
    / sum_gamma(beta_gamma / c_gamma)
```

and therefore

```text
c_mix =
    1 / sum_gamma(beta_gamma / c_gamma).
```

No density is guessed by flow. The caller must supply the target selected-phase molar
densities from the configured thermodynamic/property path.

Each target phase composition must remain strictly positive and normalized. The new
dependent composition component is selected deterministically as the largest target
mole fraction for that phase, with the lowest canonical component index winning an
exact tie. No epsilon floor, clipping, component permutation, or silent normalization
is introduced.

### 44.3 Material-balance guard

Before a target phase set is projected, its overall composition

```text
z_i(target) = sum_alpha beta_alpha * x_alpha,i
```

is checked against the normalized source pore-volume component inventory. The declared
transition material-balance tolerance is explicit. A phase-set candidate that changes
the cell's overall component composition is rejected rather than repaired.

The projection is an **initial guess for the target topology**, not an accepted flow
solution. The target `P*Nc+1` natural-variable residual must still be solved.

### 44.4 Conservative backward-Euler history migration

Previous-time component and energy accumulation are physical conserved history, not
phase-chart metadata. A phase-count transition therefore carries the previous

- component accumulation `N_i [mol / bulk-m^3]`; and
- total internal-energy accumulation `U [J / bulk-m^3]`

forward **without recomputing them from the new phase set**.

The transition history migration validates canonical component identity and frozen
porosity, then copies the previous snapshots unchanged. The previous energy snapshot
is allowed to retain its source phase-count/state identity while the current target
state uses a different active phase count. The backward-Euler energy/component
builders already consume the stored scalar history and do not require previous and
current topology to match.

The 1P and 2P PETSc assembly contexts no longer impose an extra previous-energy
phase-count equality check; the 3P context already had topology-neutral history
semantics.

### 44.5 Per-cell cardinality planning

`NaturalVariableActiveSetLayout` records per-cell

```text
q_c = P_c * Nc + 1
```

and deterministic prefix offsets for mixed `P_c in {1,2,3}` cells. This makes the
required heterogeneous scalar cardinality explicit. For example, with `Nc=3` and
cell phase counts `[1,2,3]`, local block widths are `[4,7,10]` and the packed
scalar count is `21`.

This slice does **not yet** replace the existing fixed-width mesh
`DofLayout/DofNumberingSnapshot` or PETSc MPIAIJ materialization with a
variable-per-cell section. Therefore fixed 1P, 2P, and 3P production systems remain
individually operational, while mixed-cardinality global PETSc numbering/assembly is
the remaining solver-structure step before automatic local switching can run in one
distributed solve.


## 45. Variable-cardinality cell DoF and PETSc scalar numbering bridge

The distributed active-set foundation now has a dedicated ragged natural-variable
numbering bridge in `flow_discretization_petsc`. It does not change the generic
fixed-width `mesh::DofLayout`; instead it consumes the existing cell
`PartitionSnapshot` plus the active phase count `P_c` for every local owned/ghost
cell.

For each cell,

```text
q_c = P_c * Nc + 1,  P_c in {1,2,3}.
```

Local packed DoFs are entity-major with variable width. PETSc scalar ownership is
rank-contiguous; within each rank, owned cells are ordered by stable
`GlobalEntityId`. Each owner publishes its stable cell ID, active phase count,
scalar width and global scalar start. Ghost copies must resolve to exactly the same
published width/start, otherwise construction fails collectively.

The bridge exposes:

- per-cell local packed scalar offset and width;
- per-cell PETSc global scalar start;
- local packed scalar -> PETSc global scalar mapping;
- owned PETSc scalar range and global scalar count;
- a ragged local `PetscSection` with one DoF block per cell;
- a distributed PETSc `Vec` whose ownership range matches the ragged numbering.

The 2-rank regression uses `Nc=3` with three stable cells simultaneously carrying
`P=[1,2,3]`. Their widths are `[4,7,10]`, the global scalar count is `21`,
and owner/ghost copies agree on every scalar index. A deliberately inconsistent ghost
phase count is rejected collectively.

This bridge establishes heterogeneous distributed numbering only. Complete
component/energy/fugacity residual assembly and MPIAIJ symbolic preallocation still
need to be generalized from uniform `q` to per-row-cell/per-column-cell widths before
mixed-cardinality SNES can assemble one coupled Jacobian.

## 46. Production post-SNES PT phase-transition scanner

The mixed-cardinality outer controller now has a production scanner bridge from one
converged frozen-cardinality flow state to the existing model-neutral
`PtFlashBackend`.

For every locally owned cell the scanner:

1. re-evaluates the same production cell closure used by SNES at the converged
   natural-variable state;
2. derives the overall feed from the current pore-volume component inventory,
   `z_i = N_i / sum_j N_j`, rather than substituting a phase composition;
3. calls the configured PT backend only after SNES convergence;
4. treats an accepted equal-cardinality result as a completed scan with no transition;
5. converts an accepted different-cardinality result through the existing
   `make_phase_set_transition_candidate_from_flash()` evidence guard; and
6. leaves target phase molar density to an explicit resolver callback because the
   generic PT publication does not promise a density or even a published `Z` for
   every provider/topology.

A missing target-density resolution, a non-accepted/indeterminate PT result, or a
three-phase state whose capillary offsets imply different active phase pressures is a
**scan-indeterminate** result. It is not silently interpreted as a stable phase set.
The post-SNES controller therefore accepts the timestep only when every MPI rank
reports a complete scan and the global accepted transition batch is empty.

The scanner never changes phase cardinality inside an SNES callback, never changes
EOS/flash tolerances, and never infers oil/gas/water identity from phase slots,
density, `Z`, or provider ordering. The physical SW92 Sample-6 PT backend remains
the real 2->3 transition regression for this bridge; controlled flow fixtures test
only orchestration and failure semantics.


## 47. PR76 selected-phase property closure

The first production property bridge from a fixed PR76 selected branch into the
natural-variable flow chart is explicit about model ownership.

PR76 supplies:

- selected-root `ln(phi_i)`;
- selected-root molar density `c [mol/m^3]`.

The flow-facing closure additionally requires every ordered component to carry an
explicit molar mass. It then derives

```text
rho_mass = c * sum_i x_i M_i
```

without a fitted conversion or phase-name heuristic.

PR76 does **not** define a viscosity correlation or an absolute caloric reference.
Therefore the closure requires an explicit scalar-generic transport/caloric provider
for dynamic viscosity and total specific enthalpy. The provider must carry provenance
and must preserve AD dependence. Missing molar mass, viscosity/enthalpy capability,
invalid selected-root topology or non-finite property values are hard failures; no
default constant property, clipping, finite-difference derivative or silent model
substitution is permitted.

Specific internal energy is derived from the supplied total specific enthalpy and the
selected PR76 mass density by the thermodynamic identity

```text
u = h - p / rho_mass .
```

The v1 bridge is explicitly `pc=none`: all active phases are evaluated at the
reference pressure. The production bridge now enforces that capability boundary:
every caller-owned three-phase saturation carrier must publish exactly zero
capillary offsets, `p_alpha == p_ref` for every active phase, and a phase-pressure
Jacobian containing only `dp_alpha/dp_ref = 1`. A nonzero capillary/phase-pressure
value, or a zero primal value with any nonzero capillary/phase-pressure derivative,
is a hard unsupported-capability error rather than a nonlinear domain failure.
Capillary-pressure-resolved phase pressures remain a later extension.

One topology-neutral chart evaluator now covers `P=1/2/3` with

```text
q = P*Nc + 1
```

and publishes, for each active selected phase:

- molar density;
- mass density;
- dynamic viscosity;
- specific enthalpy;
- specific internal energy;
- ordered `ln(phi_i)`;
- full first-order derivatives with respect to the frozen natural-variable chart.

For `P>1`, the same AD pass also publishes the local fugacity-equilibrium residual
and Jacobian using phase0 as the reference. Saturation columns remain exactly zero in
this `pc=none` property layer.

The `P=3` convenience bridge converts this publication directly into the existing
`NaturalVariableCellState3P`, molar-density, transport, caloric and fugacity
linearization carriers consumed by the production assembly. Relative permeability,
capillary pressure, stationary-rock thermal storage and mesh/PETSc assembly remain
separate caller-owned models.

Current validation uses a synthetic differentiable viscosity/enthalpy provider only
as a software oracle. It verifies the production contract, PR76 kernel reuse and AD
Jacobian plumbing; it is not a physical transport/caloric validation dataset. A real
case must supply independently traceable viscosity and enthalpy models before the
non-isothermal flow result can be called physically validated.


### 47.1 Repository-curated methane/ethane/propane provider

For the exact repository-curated PR76 ternary snapshot

```text
methane / ethane / propane
dataset = DeitersBell-aic16730-PengRobinson1976-ternary
```

the flow module provides a narrow production transport/caloric provider. It is
not a general hydrocarbon property database and rejects another dataset identity,
component count or component order.

The provider uses three independently attributed property layers:

1. **Molar mass and ideal-gas heat capacity:** NIST Chemistry WebBook SRD 69
   methane, ethane and propane species records. The common implemented temperature
   interval is `298.15 <= T <= 1500 K`; values outside that interval are rejected,
   not extrapolated.
2. **Dense-fluid viscosity:** Lohrenz-Bray-Clark with Stiel-Thodos dilute
   pure-component viscosity, Herning-Zipperer dilute-mixture blending and Kay
   pseudo-critical mixing. The LBC polynomial uses the audited coefficient
   `a5=0.0093324`. NIST critical molar-volume/density records supply the
   methane/ethane/propane critical-volume inputs. No C7+ critical-volume
   correlation and no tuning coefficient is enabled.
3. **Caloric departure:** the standard Peng-Robinson residual-enthalpy expression,
   using the same selected PR76 root as fugacity/density.

The ideal-gas enthalpy reference is deliberately fixed to

```text
h_i^ig(298.15 K) = 0
```

for each component. This is a reference choice for nonreactive flow, not a claim
that chemical formation enthalpies are zero. Any future boundary, well or coupled
energy model must use the same reference before energy terms are combined.

The PR residual enthalpy needs `da_mix/dT`. The implementation obtains it with one
nested forward-AD temperature direction through the explicit PR pure/mixing kernel.
The outer natural-variable AD dependence is therefore preserved. The selected
compressibility root itself remains on the existing first-order PR76 implicit-root
IFT path; nested AD is intentionally **not** applied to the root solver because that
solver does not currently publish a second-order implicit-root contract.

The convenience entry

```cpp
make_pr76_methane_ethane_propane_property_closure(model, selections)
```

binds this provider and its provenance to the existing P=1/2/3 selected-phase
property closure. Unsupported datasets or missing/mismatched molar masses fail
explicitly.


### 47.2 PR76 production cell evaluator bridge and first real transient step

The selected-phase property publication now has symmetric conversion helpers for
`P=1/2/3`. They convert the PR76 property chart into the existing single-, two-
and fixed-three-phase flow carriers without copying EOS, NIST or LBC equations
into the discretization layer.

The PETSc adapter `pr76_production_cell_evaluator.hpp` consumes those carriers
and adds only caller-owned constitutive pieces:

- stationary-rock thermal storage;
- two-phase relative permeability;
- three-phase relative permeability/capillary-pressure linearization.

The PR76 property backend remains responsible for density, viscosity, enthalpy,
internal energy and fugacity. The PETSc adapter does not infer phase roles, choose
a transport correlation, repivot compositions, or mutate phase cardinality inside
SNES callbacks.

The first end-to-end physical regression deliberately freezes a single-phase
topology on two distributed cells. It uses the repository-curated
methane/ethane/propane PR76 + NIST + LBC property provider, a positive internal
TPFA transmissibility, nonzero gravity, positive thermal conductance and a
positive backward-Euler timestep. The previous-state inventories are evaluated
with the same production property closure; the nonlinear state is not constructed
from a manufactured zero-residual target.

Acceptance requires:

- nonzero Darcy volumetric flux;
- nonzero component molar face flux;
- nonzero advective enthalpy flux;
- nonzero conductive heat flux;
- exact volume-weighted internal-face component/energy conservation to numerical
  roundoff;
- fresh central-perturbation agreement of the complete distributed residual
  Jacobian in pressure, temperature and one independent-composition direction;
- PETSc `SNESNEWTONLS + BT / GMRES / ASM(1)` convergence;
- independent final residual reassembly;
- closed-domain total component and energy inventory conservation.

Phase appearance/disappearance and outer phase-set rebuild are intentionally not
part of this regression. They remain independently validated by the existing
post-SNES transition controller tests.


### 47.3 Adaptive timestep retry / cutback / growth contract

The PETSc flow orchestration layer now owns an explicit adaptive timestep state
machine outside all SNES callbacks.

One physical timestep is attempted with a frozen accepted history and a candidate
`dt`. A rejected attempt never advances the accepted history. The attempt
callback may retain a pending converged state/system internally, but the separate
commit callback is invoked exactly once and only after a stable accepted attempt.

Recoverable rejection outcomes are:

- nonlinear solve divergence;
- nonlinear function/Jacobian domain error;
- indeterminate post-SNES phase-set scan;
- phase-transition restart-budget exhaustion;
- detected phase-set cycle.

All of those outcomes request a cutback while retry budget and the minimum
timestep permit it. PETSc/API/configuration errors are not converted into a
smaller timestep: they propagate as hard errors.

The timestep policy is explicit and versioned by
`AdaptiveTimestepControllerOptions3D`:

```text
minimum_timestep_seconds
maximum_timestep_seconds
cutback_factor              0 < f < 1
growth_factor               g > 1
maximum_retries
growth_nonlinear_iteration_limit
growth_line_search_direction_change_limit
growth_transition_restart_limit
```

On rejection:

```text
dt_retry = max(dt_min, cutback_factor * dt_current)
```

until either a stable attempt is accepted, the retry budget is exhausted, or
`dt_min` has already been reached.

On acceptance, growth is deliberately conservative. The next timestep grows only
when the accepted attempt is below all configured effort thresholds and has no
domain-error evidence. Otherwise the accepted `dt` is held. Growth is capped at
`dt_max`.

The contract provides adapters from:

- fixed-cardinality SNES success reports;
- fixed-cardinality SNES failure diagnostics, including domain-error evidence;
- post-SNES phase-transition controller outcomes.

The real CH4/C2H6/C3H8 PR76 transient regression feeds its actual converged SNES
report through this same adaptive decision path. Scripted two-rank policy
regressions independently verify cutback, retry-budget exhaustion, minimum-dt
termination, conservative growth/hold behavior, hard-error propagation, and the
rule that rejected attempts never commit history.

This contract does not yet advance physical simulation time over multiple
accepted steps, persist restart/checkpoint data, or choose model-specific
phase-transition tolerances. It establishes the deterministic timestep-control
boundary that those later features must consume.

### 47.4 Production single-phase adaptive attempt bridge

The fixed-cardinality single-phase production path now has a concrete adapter from
AdaptiveTimestepControllerBindings3D::attempt to the existing physical
backward-Euler/PETSc solve.

For every candidate timestep it freshly constructs the
SinglePhaseSnesAssemblyContext3D from the caller-owned accepted history, using the
candidate dt in the accumulation terms. It then freshly builds:

- the complete physical residual/Jacobian snapshot;
- frozen initial row equilibration;
- MPIAIJ structural materialization; and
- the existing SNESNEWTONLS + BT / GMRES + ASM(1) nonlinear solve.

SNES non-convergence and nonlinear function/Jacobian domain diagnostics are mapped
back into the adaptive-timestep outcome. PETSc/API/configuration failures remain hard
errors.

A converged SNES solve is not automatically accepted. A mandatory post-SNES review
callback owns the phase-set/stability decision. Only stable_phase_set retains a
pending solution/report. Scan-indeterminate, transition-budget and cycle outcomes
discard the trial state and request the normal adaptive cutback path. The adapter
never mutates the caller-owned accepted state or previous component/energy
accumulation snapshots.

The adaptive controller's separate commit callback must explicitly take the pending
stable result. A new attempt first destroys any still-uncommitted pending solution,
so a rejected candidate cannot leak into the next retry.

The real two-rank CH4/C2H6/C3H8 PR76 regression now runs the production adapter
twice in one physical timestep. The first real dt=0.1 s solve is deliberately
converted to a recoverable rejection by a test-only wrapper after the real solve;
its pending state is discarded. The controller cuts back to dt=0.05 s, rebuilds
the backward-Euler system, and runs a second real PR76/PETSc solve. Acceptance
requires:

- the state and previous component/energy history observed before both attempts to
  be exactly the original accepted history;
- exactly one commit, after the second stable solve;
- the committed history to equal the independently re-evaluated accepted-state
  component and energy accumulations;
- independent final residual reassembly using the original previous-time history
  and the accepted 0.05 s timestep;
- the existing closed-domain component/energy conservation checks.

The injected first rejection exists only in the regression harness; production code
does not manufacture a nonlinear or physical failure. This bridge remains
fixed-cardinality. Actual 1 <-> 2 <-> 3 topology rebuild continues to belong to
the existing post-SNES transition controller.


### 47.5 Real PR76 post-SNES equal-cardinality PT review

The real two-rank methane/ethane/propane transient no longer uses a callback that
unconditionally labels a converged frozen single-phase solve as stable.

The post-SNES review now constructs a model-neutral
`PostSnesPtFlashSourceCellSnapshot3D` from each rank's **converged production
cell state**:

```text
converged PETSc state
-> production PR76 selected-phase cell evaluator
-> p / T / overall composition
-> Pr76PtFlashBackend
-> scan_post_snes_pt_flash_source_cell_3d
```

The PT backend is built from the same frozen PR76 parameter snapshot as the flow
property closure. The scanner therefore performs the existing production PR76
initial-stability / max-three-phase / final phase-set publication path; the flow
test does not reproduce a second stability criterion.

For this fixed-cardinality slice, acceptance is deliberately narrow:

- every MPI rank must report scanner status `complete`;
- the accepted PT result must have the same one-phase cardinality, so the scanner
  publishes no transition proposal;
- an indeterminate scan or any unresolved different-cardinality path is **not**
  converted to stable and instead maps to the adaptive recoverable
  `phase_set_scan_indeterminate` outcome.

No target-density resolver is needed on the verified equal-cardinality path,
because the production scanner returns before transition projection when the
accepted phase count matches the source phase count.

The adaptive regression still keeps its deterministic test-only first rejection
after this real PT review, solely to exercise `0.1 s -> 0.05 s` cutback. Both
real SNES attempts must first pass the production PR76 PT scanner. The regression
checks two completed equal-cardinality scans and verifies that the second scan's
pressure, temperature and overall composition exactly match the state that is
subsequently committed.

This slice does not yet feed an accepted different-cardinality proposal into the
outer `1 <-> 2 <-> 3` rebuild controller. It closes only the
`SNES converged -> real PR76 PT scan -> equal-cardinality stable acceptance`
path.

### 47.6 Lossless adaptive phase-transition handoff

A resolved different-cardinality post-SNES proposal is no longer collapsed into
phase_set_scan_indeterminate or treated as a timestep cutback.

The adaptive state machine now distinguishes stable_phase_set,
phase_transition_proposed, and recoverable timestep rejection.
phase_transition_proposed maps to phase_transition_handoff_required.

That outcome does not commit accepted state/history, does not cut back dt,
does not consume retry budget, and leaves the converged candidate owned by the
attempt context for the outer topology controller.

The single-phase production attempt review now returns the exact local-owned
PostSnesPhaseTransitionProposal3D batch. A proposed transition retains the
attempt request, converged PETSc Vec, fixed-cardinality SNES report and local
resolved proposals. take_pending_transition() transfers that bundle exactly once.

Proposal validation is ownership-aware: local proposals must refer to locally
owned cells, preserve canonical component order, remain target_resolved, carry
evidence provenance, and define a valid positive-support target phase set.
MPI requires at least one proposal globally while allowing an empty local batch
on ranks that do not own a transitioning cell.

The production PR76 PT review now preserves scanner proposals directly. On the
current equal-cardinality CH4/C2H6/C3H8 case the proposal vectors remain empty.
A future resolved different-cardinality result stops at handoff instead of being
reclassified as indeterminate.

post_snes_phase_transition_handoff_scanner.hpp is the adapter to the existing
outer controller. Its first scan replays the preserved local proposal batch;
after the first topology rebuild, later generations delegate to the normal
production scanner. The outer controller keeps ownership validation, global
batch gathering, restart-budget and phase-set-cycle checks.

A controlled regression uses a real PR76/PETSc converged single-phase solve and
then publishes a synthetic test-provenance 1->2 proposal. It verifies no commit
or cutback, one-shot Vec/report/proposal transfer, preserved evidence/material
balance, and unchanged accepted component/energy history.

The controlled proposal is software orchestration evidence, not physical PR76
1->2 evidence. This slice does not yet construct the initial variable-cardinality
PhaseTransitionRebuiltNaturalVariableSystem3D from the fixed-cardinality handoff.

### 47.7 Fixed-1P handoff to initial outer-controller system

The fixed-cardinality handoff now has a production bridge to the source-topology
PhaseTransitionRebuiltNaturalVariableSystem3D used by the existing outer controller.

The bridge consumes the converged fixed-1P SNES report, accepted cell metadata and
unchanged previous-time histories. Only owned q values are copied into the new
ragged PETSc Vec. Ghost rebuild snapshots may omit q entirely; the rebuilt mixed
physical context obtains ghost state from the existing owner/PetscSF path. No
global solution Allgather is introduced at this control boundary.

The phase-transition rebuild cell contract was also corrected to separate two
different physical inventories:

- current converged component inventory: used only to verify candidate material
  balance and construct the target-topology initial guess;
- previous-time BE component/energy history: copied unchanged into the restarted
  system and never recomputed from the candidate.

This matches the earlier documented contract: PT scanning is based on current
pore-volume inventory while backward-Euler history remains the accepted t_n state.

The real CH4/C2H6/C3H8 regression now transfers a converged fixed-1P handoff into
an initial outer-controller system and checks global scalar count 8, both local
source cells remaining 1P, exact owned-q preservation, zero absent-phase provider
calls, and independently reassembled source residual <= 1e-6.

The preserved-proposal scanner is also exercised inside the existing controlled
mixed-cardinality outer-controller fixture. Its first generation replays the
preserved local 2P->3P proposal, the existing rebuild factory rebuilds the target
ragged topology, restarted SNES runs, and the controller reaches stable_phase_set
after exactly one transition restart.

The controlled 2P->3P proposal remains orchestration evidence only. This slice does
not claim a real PR76 1->2 physical transition; it closes the software boundary from
fixed handoff materialization through outer-controller proposal consumption and
topology restart.

### 47.8 PR76 resolved target-rebuild plan

The PR76 PT-to-flow transition contract now preserves the accepted phase's opaque
`activity.branch` alongside beta, composition, molar density and evidence. The
branch token is thermodynamic provider provenance only; it is explicitly not a
physical phase identity and is never interpreted as oil/gas/water or morphology.

New `pr76_single_phase_transition_target_rebuild.hpp` converts local-owned resolved
1P->2P/3P proposals into a distributed target-rebuild plan. For each owner proposal
it:

- reconstructs the current converged 1P q from the fixed SNES report;
- re-evaluates the production cell closure and builds the current component
  inventory used by the material-balance gate;
- keeps previous component/energy accumulation as unchanged t_n history;
- converts preserved activity.branch tokens into frozen PR76 selected-root
  selections using the audited root-options snapshot;
- requires an explicit phase-identity resolver; no slot/Z/density/root heuristic
  is allowed;
- validates source->target identity continuation; and
- builds the per-transition selected-phase branch registry.

Identity resolution may return unresolved without error. In that case the plan is
not published, so production orchestration can classify the transition as
indeterminate instead of inventing morphology.

Resolved owner topology metadata is serialized and MPI-allgathered by stable cell
ID. Every local ghost copy reconstructs the same target phase count, dependent
composition pivots, explicit phase identities, selected PR76 root bindings and
transition evidence. Owner target q/history remain owner-only; ghost target q and
BE histories are deliberately absent.

The two-rank controlled regression gives each owner one resolved 1P->2P proposal
with explicit test identities. Each rank then observes both local copies as target
2P, verifies its owned target q has width 7 while the ghost q remains empty, and
checks both synchronized transition records retain the same evidence and selected
root bindings. The identity resolver is test-only evidence and does not claim real
PR76 morphology.

This slice stops at the target-rebuild plan. The next step is to consume these
per-cell selected-root/identity records in a cell-dispatch PR76 production evaluator
and resolve the required frozen absent-phase coordinates before materializing the
actual mixed-cardinality restarted SNES system.

### 47.9 PR76 mixed-cardinality target materialization

The distributed PR76 target plan can now be materialized into the existing
`PhaseTransitionRebuiltNaturalVariableSystem3D` without assigning one global EOS
root selection per cardinality.

`CellScopedMixedCardinalityEvaluatorDispatcher3D` dispatches 1P/2P/3P production
cell evaluation by stable cell ID. Each target cell therefore consumes the selected
PR76 roots frozen in its own transition plan rather than reusing another cell's
root map.

Cross-cardinality absent-phase density uses a separate explicit branch resolver.
The active neighbour's selected root is supplied only as reference evidence; the
resolver must explicitly publish the branch to use on the absent host cell. No
neighbour-root, density, Z or phase-slot fallback is applied.

Materialization is deliberately two-stage:

1. build a bootstrap ragged system with the target cell dispatch;
2. use its existing PetscSF packed-state synchronization to evaluate every local
   target host state, including ghosts;
3. for every authoritative cross-cardinality face, freeze the absent phase at the
   actual target host state, using the adjacent resolved transition composition as
   the hypothetical composition reference and zero composition tangent;
4. attach an explicit cell-scoped selected-branch provenance key;
5. build `Pr76CellScopedAbsentPhasePotentialExtensionProvider`, then rebuild the
   final mixed system with the normal thermodynamic cross-cardinality TPFA bridge.

If the required active transition phase, host state, or absent-side selected branch
cannot be resolved, materialization returns no system rather than guessing. Multiple
faces requesting the same absent cell/identity must agree on the frozen composition
and selected branch.

The 2-rank CH4/C2H6/C3H8 regression now takes only cell10 through the controlled
1P->2P transition while cell20 remains 1P. The final ragged system therefore has
`7 + 4 = 11` global scalars and one real 2P/1P cross-cardinality face. It verifies
the cell20/phase-1 frozen coordinate entry, executes the cell-scoped PR76 absent
provider through the standard TPFA bridge, and assembles finite nonzero mixed
residual/Jacobian data. The synthetic transition is not required to converge and
remains orchestration evidence only.

### 47.10 Fixed real PR76 1P->2P discovery state

A fixed CH4/C2H6/C3H8 state is now promoted from an external discovery calculation
into the repository's own PR76 backend/scanner regression:

```text
T = 300 K
P = 2.0 MPa
z = [0.20, 0.15, 0.65]
component order = methane / ethane / propane
kij = 0
```

The regression does not store the external screening phase fraction as an oracle.
The repository backend must independently prove the state. Acceptance requires:

- a structurally valid `Pr76PtFlashBackendResult`;
- exactly two accepted phases;
- fresh closed `1P -> 2P` accepted-target transition evidence;
- nondegenerate phase fractions away from 0/1;
- two compositions with clear separation rather than a coalescence-boundary state;
- exact reconstruction of the feed from beta and phase compositions;
- positive target molar densities re-evaluated on each accepted phase's own
  `activity.branch`; and
- successful conversion by `scan_post_snes_pt_flash_source_cell_3d()` into a
  resolved `PostSnesPhaseTransitionProposal3D` carrying both provider branches.

This state is therefore the first fixed methane/ethane/propane operating point in
the flow PR that is required to produce real PR76 thermodynamic `1P -> 2P` evidence.
It is intentionally separate from the existing 450 K / approximately 8 MPa stable
single-phase transient fixture.

### 47.11 Real PR76-triggered fully implicit 1P->2P restart

The fixed 300 K / 2 MPa / z=[0.20,0.15,0.65] state now drives an end-to-end
two-rank restart regression. This path no longer injects a controlled transition
proposal.

The source is a closed two-cell MPI system with one structural connection whose physical flux is exactly zero. Cell10 uses
the real PR76 two-phase-trigger state; cell20 uses the already validated stable
450 K / 8 MPa one-phase state. Accepted t_n component/energy history is built from
exactly each initial one-phase state, so the fixed-cardinality source residual is
zero and source SNES must converge without moving q.

The normal adaptive post-SNES review then runs `Pr76PtFlashBackend` with selected-
branch molar-density resolution. Exactly one proposal is permitted globally:
cell10 1P->2P. Rank1 must publish no local proposal.

The retained real proposal flows through target-rebuild planning, per-cell PR76
selected-root closure construction, mixed 2P/1P ragged materialization and a
restarted PETSc SNES solve. Physical phase identity remains an explicit opaque
regression mapping because the current PR76 backend deliberately reports
`morphology_resolved=false`; EOS root indices are not renamed liquid/vapor.

The final restarted system must contain 11 global scalars (cell10 q=7, cell20 q=4),
converge with finite final SNES residual <= 1e-6, leave cell10 on strict-positive
two-phase saturation support, and pass a fresh production PT rescan with status
`complete` and no further transition proposal.

### 47.11.1 Production-range correction

The earlier 250 K flash point remains a valid PR76 thermodynamic two-phase test,
but it is below the production methane/ethane/propane caloric provider's sourced
NIST ideal-gas Cp lower bound of 298.15 K. It is therefore not used as a flow
restart source and no enthalpy extrapolation is introduced.

The production restart source is now the fixed 300 K / 2.0 MPa /
z=[0.20,0.15,0.65] state. This point lies inside the caloric provider temperature
range and retains an explicit metastable one-phase PR76 root for the frozen source
chart. The repository backend/scanner still decides whether it is accepted as the
real 1P->2P trigger; no external phase fraction is used as an oracle.

### 47.11.2 Preserve production ASM graph while freezing source flux

The source restart regression no longer removes the inter-cell connection. That
isolated algebraic graph was rejected by the audited GMRES+ASM solve contract before
SNES callbacks were entered (`PETSC_ERR_ARG_INCOMP`).

The source now retains the normal two-cell structural stencil and preallocation.
Both cells use 2.0 MPa reference pressure, the transition face uses zero gravity
and zero thermal conductance, so Darcy/component/energy face contributions are
exactly zero at the accepted source state while the off-diagonal Jacobian graph
remains present for ASM. The preflight independently requires the materialized
physical residual L2 norm to be <= 1e-12 before the source SNES is invoked.

For the restarted 2P/1P face, the absent-side PR76 branch is no longer a controlled
`root=0` fixture. A resolver evaluates `roots_full` at the absent host p/T and
hypothetical target composition and publishes a selection only when exactly one
positive-slope derivative-valid PR76 branch exists. Ambiguous root topology remains
unresolved instead of being guessed.

### 47.11.2 Frozen row scaling for restarted ragged SNES

The real PR76 1P->2P regression exposed a solver asymmetry: fixed-cardinality
production SNES already supported frozen left row equilibration, while the
variable-cardinality restart path solved the native component/energy/fugacity rows
without scaling.

`make_variable_cardinality_initial_row_equilibration_3d()` now independently
assembles the analytic ragged Jacobian at q0 and freezes
`D_i = 1 / max_j |J_ij(q0)|` for every owned scalar row. Rows without a finite
nonzero analytic coefficient are rejected.

`solve_variable_cardinality_natural_variable_snes_3d()` accepts the same optional
positive D contract as the fixed solver. The shared PETSc callbacks solve `D*R=0`
with `D*J`; physical residual/Jacobian evaluators remain unscaled. Newton,
backtracking, GMRES, ASM overlap, tolerances and analytic Jacobian policy are
unchanged.

`PhaseTransitionRebuiltNaturalVariableSystem3D::solve()` now creates and destroys
this frozen scaling automatically. The generic 1P/2P/3P ragged manufactured solver
regression also runs through the scaled path, so this is not a PR76 special case.

### 47.11.3 Ragged ASM stabilization and failure diagnostics

After frozen row scaling, the real restart still did not converge. Audit found the
variable-cardinality ASM sub-block policy had not yet mirrored the fixed solver's
scaled-real-SI safeguards. With D enabled it now uses the same `1e-10` weak-diagonal
reordering threshold and PETSc `MAT_SHIFT_NONZERO` stabilization on the private
ASM LU preconditioner matrix. The physical analytic Jacobian and nonlinear root are
unchanged.

The variable-cardinality solver also accepts optional
`NaturalVariableSnesFailureDiagnostics3D`. On a non-converged solve it records
SNES/KSP/top-level PC/ASM-sub-KSP/sub-PC reasons, iteration/evaluation/domain-error
counts and the PETSc function norm before returning `PETSC_ERR_NOT_CONVERGED`.
This is diagnostic evidence only and does not convert divergence into success.

### 47.12 Accepted physical-time clock and history advancement

The first real multi-timestep ownership boundary now sits above the existing
adaptive-timestep and variable-cardinality PETSc solvers.

`AcceptedPhysicalTimeClock3D` owns only accepted physical-time metadata:
accepted time, accepted-step count and the next proposed timestep. Time does not
advance on a rejected nonlinear attempt, cutback, or phase-transition handoff.
A phase-transition restart remains part of the same physical timestep.

After a stable solve has passed its production phase-set review,
`commit_accepted_physical_timestep_3d()` performs the accepted-step commit in
this order:

1. evaluate the converged mixed-cardinality cell state;
2. rebuild the owned backward-Euler component and energy history from that exact
   accepted state;
3. re-anchor every frozen absent-phase coordinate chart at the accepted host
   state while preserving the selected PR76 branch provenance;
4. copy the accepted PETSc state as the next timestep initial state;
5. install the adaptive controller's accepted next `dt`;
6. only then advance the physical clock.

The real two-rank CH4/C2H6/C3H8 regression now accepts two consecutive physical
steps. The first step consumes the production PR76 `1P -> 2P` proposal and
ragged `2P/1P` restart, so the adaptive policy holds `dt=1 s` because one
transition restart occurred. The committed history is then used by a second
fully implicit solve of the same physical system. Its fresh PT rescan remains
stable, the zero-effort accepted step requests growth to `dt=2 s`, and the
clock reaches `t=2 s` after exactly two accepted steps. A fresh residual
evaluation after the second commit checks that the rebased state/history is a
valid baseline for the next timestep.

This slice does not add a well, boundary/source control, checkpoint file,
long-duration schedule, or any new phase-transition physics.


### 47.13 Production one-physical-timestep driver

The PETSc flow layer now owns one explicit production orchestration boundary:

`advance_one_physical_timestep_3d()`.

The driver consumes an accepted `PhaseTransitionRebuiltNaturalVariableSystem3D`
and performs one physical backward-Euler step only. Its invariants are:

- each nonlinear retry starts from the same accepted state/history;
- trial `dt` may change for cutback without rebasing component/energy history;
- a phase-transition target is rebuilt on a disposable trial system, so failed
  ragged/topology restarts cannot replace the accepted topology;
- transition restarts remain inside the same physical timestep and are included
  in the adaptive effort/restart budget;
- only a stable post-SNES phase-set result reaches
  `commit_accepted_physical_timestep_3d()`;
- the accepted clock/history advance exactly once, after the final stable
  candidate has been selected.

`MixedCardinalityPhysicalSnesAssemblyContext3D::set_trial_timestep_seconds()`
and the corresponding rebuilt-system facade change only the trial backward-Euler
`dt`; they do not modify the stored initial state, accepted accumulation
history, physical phase identity or frozen absent-phase coordinates. Terminal
rejection restores the entry `dt`.

The real two-rank CH4/C2H6/C3H8 regression now uses this production driver for
the materialized PR76 2P/1P restart and the following physical step instead of
manually sequencing `solve -> PT scan -> adaptive decision -> history commit`.
The already-materialized 1P->2P restart is carried into the first driver call as
one restart belonging to that same physical timestep, so the first accepted
`dt=1 s` is held; the next stable step grows to `dt=2 s`. A forced
post-SNES-indeterminate attempt with zero retry budget additionally verifies that
terminal rejection leaves accepted time, step count, next proposed `dt`,
state and component/energy history unchanged.

This slice does not change PR76/flash physics or target materialization and does
not add wells, boundary/source controls, checkpointing or a long-duration
schedule.

### 47.14 Production physical-time interval loop

The PETSc flow layer now has a production interval-level orchestration entry:

`advance_physical_time_to_3d(..., target_time_seconds, ...)`.

It is deliberately a thin owner above `advance_one_physical_timestep_3d()`:
the loop never solves equations, scans phase stability or commits history by
itself. For every accepted step it delegates those responsibilities to the
one-timestep driver and stops immediately if that driver returns a terminal
timestep rejection.

The only interval-specific operation is end-time alignment. Before each step the
loop supplies a first-attempt `dt` cap equal to the remaining physical time.
The driver validates the accepted system against the clock's original proposal,
then uses `min(next_dt, remaining)` as the first trial. If a clipped trial is
rejected, the one-step rollback restores the entry proposal and leaves accepted
time/history unchanged. If the remaining terminal interval is below the ordinary
adaptive minimum, that remainder is allowed as the exact terminal first attempt
and becomes that step's retry floor; failure therefore exits instead of silently
stepping below the requested end-time boundary.

The real two-rank CH4/C2H6/C3H8 regression now keeps the production PR76
transition-restarted first step at `t=1 s`, then advances the same accepted
2P/1P system through this interval loop to `t=2.5 s`. The loop accepts one
ordinary `dt=1 s` step, observes adaptive growth, then truncates the proposed
`dt=2 s` to the remaining `0.5 s` without overshoot. The resulting accepted
step count is three and the next proposal is `1 s` from the clipped step's
growth policy. A second interval call forces an indeterminate post-SNES review
with zero retry budget and verifies immediate failure exit with no change to
accepted time, step count, state, history or the pre-call next-`dt` proposal.

This slice still does not introduce wells, boundary/source schedules,
checkpoint/restart files or long-duration case scheduling.

### 47.15 Model-neutral conservative cell source

The mixed-cardinality production assembly now accepts an optional borrowed
`MixedCardinalityPhysicalCellSourceEvaluatorBinding3D`. It is evaluated only
for cells owned by the current MPI rank; ghost cells never evaluate or insert
an external source.

The lower `flow_discretization` payload is
`CellSourceLinearization3D`. Positive component molar rate [mol/s] and positive
energy rate [W] mean injection into the control volume. With the existing
finite-volume residual convention, source contributions are
`R_i^src=-q_i/V_b` and `R_E^src=-Q_E/V_b`. Analytic source derivatives are
normalized by the same frozen rigid-grid bulk volume and inserted only into the
local diagonal natural-variable Jacobian block.

The source is external and therefore does not participate in internal-face MPI
exchange or closed-face zero-sum checks. The generic phase-transition rebuild
and PR76 target materializer preserve the borrowed source binding, while
compatibility overloads keep all existing no-source callers unchanged.

Tests cover exact source sign/unit/Jacobian normalization, owner-only 2-rank
PETSc residual/Jacobian insertion, and the real CH4/C2H6/C3H8 PR76 path: the
source stays disabled through the physical 1P->2P restart, is enabled on
accepted cell20, creates a nonzero production residual, and is then consumed by
the existing production multi-timestep loop.

This slice does not implement Peaceman well index, BHP/rate control, well
unknowns, boundary conditions, source schedules or facilities.


### 47.25 Li-Firoozabadi six-component sourced flow-property benchmark

The flow thermodynamics boundary now also exposes

`<mpmc/flow/pr76_li_firoozabadi_sour_gas_properties.hpp>`

for the repository's Li-Firoozabadi-2012 six-component acid-gas PR76 dataset.
It does not introduce a new transport or caloric law: it reuses the established
selected-PR76 + Stiel-Thodos/Herning-Zipperer/Lohrenz-Bray-Clark + sourced
ideal-gas-Cp + PR-departure structure, with NIST low-temperature data restricted
to the explicit 100-298.15 K interval.

The first consumer is a one-cell stationary non-isothermal serial PETSc
short-step at the independent 20-bar / 178.8-K three-phase equilibrium.  Its
fixed BHP equals phase pressure, so the independently expected well rate is
exactly zero and all six component inventories plus total internal energy are
invariants of the 1-s backward-Euler step.  NIST SRD 30 low-temperature quartz
data supply the stationary-rock thermal derivative; SPE1/Odeh supplies the
Cartesian cell/porosity/permeability/well-radius engineering inputs.

The regression's local `kr=S` callback is explicitly software-structural only.
Because the benchmark has no internal faces and zero well drawdown, that callback
does not enter any externally validated rate or inventory value and is not a
physical V-L1-L2 relative-permeability claim.  Full provenance and the independent
oracle contract are recorded in
`tests/flow_discretization/petsc/pr76_li_firoozabadi_sour_gas_short_step.md`.


## 48. SW92 selected-phase property -> flow natural-variable bridge

The flow thermodynamics boundary now exposes

`<mpmc/flow/sw92_selected_phase_property_closure.hpp>`

for a frozen SW92 selected family/root chart.  Each active phase owns an
explicit `Sw92SelectedPhase` carrying NaCl molality, AQ/NA family, algebraic
root index and root options.  The closure never changes those identities while
evaluating a natural-variable residual/Jacobian.

SW92 itself supplies selected-branch `ln(phi)` and molar density through the
existing scalar-generic thermodynamics façades.  Mass density is formed only
from explicit ordered component molar masses.  Dynamic viscosity and absolute
specific enthalpy remain an explicit caller-owned scalar-generic
`SelectedPhaseTransportCaloricValues` provider; this slice introduces no SW92
transport or caloric correlation.  Specific internal energy follows
`u=h-p/rho_mass`.

The existing topology-neutral P=1/2/3 property chart driver is now
closure-generic, so SW92 reuses the already audited natural-variable carriers
for molar density, transport, caloric properties and fugacity-equilibrium rows
without duplicating a second flow algebra implementation.  The public SW92
wrapper publishes value plus analytic forward-AD Jacobians with respect to
pressure, temperature and independent phase compositions.  Family, root and
molality are frozen discrete provenance, not differentiable unknowns.

This v1 bridge is explicitly `pc=none`: all active SW92 selected phases are
evaluated at `p_ref`.  A three-phase downstream saturation/pressure carrier
with any nonzero capillary/phase-pressure offset, or with any phase-pressure
Jacobian other than exact `dp_alpha/dp_ref=1`, is rejected as
`Sw92SelectedPhasePcNoneCapabilityError`.  No pressure offset is silently
ignored.

The flow-core regression covers 1P/2P/3P publication, direct equality with the
existing SW92 selected fugacity/density kernels, fresh-perturbation checks of
the forward-AD Jacobian, component permutation, frozen family/root/molality
identity, missing molar-mass and invalid-selection failures, root-resolution
failure propagation, and primal/Jacobian pc=none capability rejection.

This slice does not add an SW92 viscosity/enthalpy model, does not run an SW92
PETSc timestep, and does not change flash topology, wells, phase-transition
orchestration, solver tolerances or any EOS equation.


## 49. Source-backed zero-salinity SW92 CO2/H2O transport/caloric provider

The first non-synthetic SW92 transport/caloric provider is deliberately narrow:

`<mpmc/flow/sw92_co2_water_properties.hpp>`

supports only the repository-curated corrected-original **CO2/H2O binary at
NaCl molality = 0**.  Nonzero molality is rejected explicitly; the provider
does not reuse pure-water parameters as an unvalidated brine model.

Dynamic viscosity uses the Chung-Ajlan-Lee-Starling (1988) high-pressure
polar-mixture correlation (DOI `10.1021/ie00076a024`) with the standard Chung
mixture rules and Neufeld collision integral.  CO2 critical molar volume is the
NIST SRD 69 Li-Kiran record; water critical volume is obtained from the IAPWS
critical density `322 kg/m3` and the NIST molecular weight.  The water polar
inputs are `mu=1.8546 D` and Chung association factor `kappa=0.076`.
Chung transport binary `xi/zeta` parameters are unity, as in the standard
method.

Caloric enthalpy is the sum of:

1. CO2: NIST SRD 69 / Chase-1998 gas-phase Shomate sensible
   `H(T)-H(298.15 K)`;
2. H2O: IAPWS-95 R6-95(2018) ideal-gas Helmholtz contribution, evaluated as
   `h0(T)-h0(298.15 K)`; and
3. the standard Peng-Robinson departure expression evaluated using the **SW92
   family-specific** `a(T,x)`, `b(x)`, water alpha and BIPs for the frozen
   selected family/root.

The provider interval is now exactly `300 <= T <= 1200 K`. CO2 remains inside
its NIST 298--1200 K Shomate interval; H2O uses one continuous IAPWS-95
ideal-gas expression across the whole provider interval, so no 500 K caloric
splice is introduced. The reference is a nonreactive-flow sensible enthalpy,
not heat of formation.  Internal energy continues to be
derived by the selected-phase closure as `u=h-p/rho_mass`.

The provider validates the exact repository dataset/revision, component
identity, SW92 Table-3 `Tc/Pc/omega`, and NIST molar masses before use.
Component order may be permuted, but the physical source identity may not.

The flow-core regression separately freezes a source-formula oracle for the
NIST Shomate mixture enthalpy and Chung viscosity, checks component permutation
and range/salinity failures, then runs the real SW92 selected-root closure and
compares its pressure/temperature/composition forward-AD viscosity, enthalpy
and internal-energy derivatives against fresh central perturbations.

This slice does **not** add a brine transport model, thermal conductivity,
SW92 PETSc timestep, well coupling or phase-transition orchestration.


## 50. SW92 production cell evaluator and stationary PETSc short-step

The PETSc flow bridge now separates the model-neutral selected-phase production
cell evaluator core from EOS-specific wrappers:

- `selected_phase_production_cell_evaluator.hpp` owns the common 1P/2P/3P
  cell-evaluation pipeline, rock-storage callback contract, two-phase
  relative-permeability callback contract and three-phase saturation callback
  contract;
- the existing `pr76_production_cell_evaluator.hpp` remains source-compatible
  as a thin PR76 traits wrapper;
- `sw92_production_cell_evaluator.hpp` supplies the SW92 traits wrapper and
  maps SW92 `pc=none` capability failures to `PETSC_ERR_SUP`.

The SW92 wrapper also provides
`materialize_sw92_co2_water_profile_c_frozen_cell_3d(...)`.  It consumes an
already accepted authoritative Profile-C phase set, preserves its ordered
components, AQ/NA family and selected algebraic root, chooses a
well-conditioned dependent composition component per phase, converts mole
phase fractions to volume saturations from the selected SW92 molar densities,
and constructs the frozen 1P/2P/3P natural-variable chart.  It performs no
flash, stability search, phase transition or phase-identity inference.

The first production regression is intentionally stationary and minimal. The
existing authoritative zero-salinity CO2/H2O **W+H two-phase** Profile-C
regression state at 3 MPa / 340 K / z=[0.7,0.3] is materialized into the real
`Sw92Co2WaterPropertyProvider`, inserted through
`MixedCardinalityPhysicalSnesAssemblyContext3D`, and advanced as a 1 s
Backward-Euler system with no faces, source or well.  The initial accepted
state is therefore the exact nonlinear solution.

The regression checks component, energy and fugacity-equilibrium residuals at
the frozen state, compares the assembled analytic/AD Jacobian against fresh
pressure/temperature/saturation/composition central perturbations, solves
through the existing variable-cardinality PETSc
`SNESNEWTONLS -> GMRES -> restricted ASM` path, and verifies that accepted
component inventories and total internal energy are unchanged.

That production-cell slice remains `pc=none` and did not itself add SW92
post-SNES phase scanning, 1<->2<->3 restart, wells, brine transport, face fluxes
or long-time stepping; the following section adds only the scanner/materializer
boundary, not the restart.


## 51. SW92-aware post-SNES scan -> authoritative target materialization

The generic PT backend remains intentionally role-neutral: it does not publish
SW92 AQ/NA family labels, physical-role metadata or selected-root provenance.
Those fields are therefore **not** added to `PtFlashBackendResult`.

The PETSc bridge instead exposes

`<mpmc/flow_discretization_petsc/post_snes_sw92_profile_c_phase_transition_scanner.hpp>`

as a model-specific post-SNES boundary.  For each owned cell it runs the
boundary-aware Profile-C topology solve exactly once.  Two pure projections are
then taken from that same owned solve:

1. the existing role-neutral PT phase-set + transition evidence, used to build
   the model-neutral `PostSnesPhaseTransitionProposal3D`; and
2. the authoritative Profile-C publication, retaining ordered components,
   prescribed NaCl molality, AQ/NA family and selected algebraic root.

When the accepted cardinality differs from the frozen source cardinality, the
scanner evaluates each accepted SW92 selected-phase molar density on the
authoritative family/root, constructs the existing generic
`PhaseSetTransitionProjection` (dependent-component pivots, volume
saturations and target natural variables), and publishes a
`Sw92AuthoritativeTargetMaterialization3D` sidecar.  Proposal and sidecar are
one-to-one by stable cell ID and share the same transition evidence profile.
No liquid/vapor, oil/gas/water or cross-cell physical phase identity is inferred
from root, slot, density or compressibility factor.

The scanner context keeps only the current scan generation's local-owned
materializations so a later rebuild factory can consume them transactionally.
If any owned cell is indeterminate, the proposal list and sidecar cache are
both cleared.

The regression freezes two boundaries:

- the existing 3 MPa / 340 K zero-salinity CO2/H2O authoritative 2P state is a
  same-cardinality no-op and produces no proposal/sidecar;
- the sourced Sample-6 10 MPa / 350 K 2->3 topology path produces a resolved
  generic proposal plus an authoritative AQ/NA/NA sidecar whose family, root,
  density, component order and target natural-variable projection are checked.
  Generic material balance is checked against the source inventory, and the
  existing history-migration contract verifies component and total-energy
  histories are preserved exactly.

This slice does not rebuild/destroy an SNES system, restart a timestep, couple a
well, infer physical phase identities or add transport data.


## 52. SW92 transactional same-dt phase-transition restart

The authoritative SW92 target sidecar can now drive the existing model-neutral
post-SNES transition controller through a real transactional rebuild/re-solve
cycle.

`sw92_transactional_phase_transition_restart.hpp` adds a rebuild factory that
consumes the scanner generation's sidecars and the generic accepted transition
batch.  For each changed cell it:

- validates source/target cardinality and stable cell identity;
- obtains the target active-phase identity map from an explicit caller resolver
  (family/root/slot are **not** promoted to physical phase identity);
- reuses `make_accepted_phase_transition_rebuild_cell_3d` so target q,
  component material balance and frozen component/energy history use the same
  generic contract as existing phase-transition rebuilds;
- constructs a new SW92 selected-phase property closure directly from the
  authoritative target molality/family/root selections;
- installs the matching 1P/2P/3P production evaluator in a cell-scoped
  dispatcher; and
- calls the normal `rebuild_phase_transition_natural_variable_system_3d` with
  the **same** `time_step_seconds`.

The generic controller already destroys the converged source candidate after a
successful rebuild and immediately solves the rebuilt system in the same
physical timestep.  No accepted-history rebase occurs between generations.

The first production closure regression uses the sourced zero-salinity CO2/H2O
provider.  A deliberately frozen 1P NA cell at 3 MPa / 340 K /
z=[0.70,0.30] has zero backward-Euler residual on its own accepted baseline.
The authoritative scanner resolves the same inventory to W+H 2P; the source
candidate is discarded, the cell is rebuilt with authoritative AQ/NA
family/root selections, and the 2P system is solved again at the same 1 s dt.
The controller must report exactly one transition restart followed by a stable
generation, while final component inventory and total internal-energy history
still match the original pre-restart baseline.

### Current capability boundary

SW92 still has no validated family/root-aware absent-phase thermodynamic
extension provider.  Therefore this v1 transactional rebuild explicitly
requires **no authoritative faces and no local ghost overlap**.  It returns a
capability error rather than pretending cross-cardinality face transport is
available.  The rebuild factory accepts authoritative 1P/2P/3P targets.  In addition to
the sourced CO2/H2O 1P<->2P production-property regression, the PETSc test now
runs real Sample-6 SW92 thermodynamics through transactional 2P->3P and 3P->2P
same-dt rebuild/re-solve.  Sample-6 still has no source-complete eight-component
flow transport/caloric dataset in this repository, so those two algorithmic
regressions add only explicitly `synthetic_test` molar masses and a
manufactured topology-independent linear-u(T)/constant-viscosity provider.
Their rebuild contexts also explicitly opt into a conservation-storage anchor
that projects target q onto the frozen component/energy storage manifold before
SNES, and reuse the already-audited Sample-6 transition material-balance
tolerance of `2e-8` from the authoritative scanner regression.  Both are
explicit context settings: the production defaults remain no storage anchor and
the generic `1e-10` projection tolerance.  The sourced CO2/H2O
production-property regressions therefore pass the authoritative target
directly to SNES under the stricter default contract.  The manufactured values and opt-in projection never
enter production data or scientific claims; the phase equilibria,
AQ/NA families, selected roots, compositions and topology remain the real
SW92 Sample-6 solutions.

This slice does not commit physical time, does not rebase accepted history,
does not couple wells, and does not add brine transport.
