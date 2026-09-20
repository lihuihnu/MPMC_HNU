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
- public-header self containment.

The synthetic polynomial `kr` and linear capillary law used by the regression are
structural derivative fixtures only. They are not physical constitutive models and
provide no reservoir parameters or validation data.

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
modules/flow_discretization
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
