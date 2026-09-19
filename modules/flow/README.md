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
- **CPA:** selected PT density-root index + exact `p,T,x` -> `double` `ln_phi`; derivative capability is explicitly `value_only`. The current repository does not yet contain the association-state and density-root IFT derivative layer needed for a scientifically complete CPA `p,T,x` Jacobian.

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
