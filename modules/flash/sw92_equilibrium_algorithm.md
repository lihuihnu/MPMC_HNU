# SW92 equilibrium-algorithm contract

## 1. Purpose and scope

This document separates **thermodynamic model identity** from **equilibrium-algorithm identity** for the Søreide–Whitson path.

The thermodynamic profile already implemented in `modules/thermodynamics` is:

```text
SW92/corrected-original/PR76-base/NaCl-molality
```

That profile defines component parameters, the corrected-original SW92 water/brine alpha term, aqueous/non-aqueous binary-interaction rules, mixing, cubic roots, and `ln(phi)`. It does **not** by itself define how an unknown feed chooses an aqueous/non-aqueous phase family, how one or more phase families are combined in a PT equilibrium calculation, or how a final phase set is accepted.

This document defines two distinct equilibrium-algorithm profiles and their scientific status. It does not add a solver implementation, change any numerical tolerance, or claim a new phase-equilibrium validation result.

## 2. Evidence audit

### 2.1 Søreide and Whitson 1992: thermodynamic model basis

Primary source:

- I. Søreide and C. H. Whitson, *Peng-Robinson predictions for hydrocarbons, CO2, N2, and H2S with pure water and NaCl brine*, Fluid Phase Equilibria 77 (1992) 217-240, DOI `10.1016/0378-3812(92)85105-H`, including the authors' appended errata supplied to this project.

The paper states that two BIP sets are used for water/hydrocarbon mixtures and explicitly forms separate attraction parameters for the non-aqueous and aqueous phases, respectively (Eqs. 10a and 10b). The applications separately address non-aqueous-component solubility in aqueous phase and aqueous-component solubility in non-aqueous phase. The paper also states that fitted AQ interaction parameters were obtained by matching measured compositions with two-phase flash calculations.

The 1992 paper is therefore sufficient to establish **phase-family-specific thermodynamic parameterization**, but it does not specify a modern, complete, reproducible orchestration contract for initial phase-family selection, finite-search stability, arbitrary multicomponent phase-count determination, or a joint three-phase solver. MPMC must not invent those missing algorithm semantics and attribute them to the 1992 paper.

The corrected-original parameter profile remains the authoritative thermodynamic identity for the current implementation. Later correlations or tuning workflows are separate profiles.

### 2.2 Current Whitson Water bot: operational dual-model structure

Current Whitson documentation (`https://manual.whitson.com/methods/water-bot/`) explicitly says that **two EOS models are developed, one for the aqueous phase and one for the non-aqueous phase**. Its water-property procedure flashes a constructed water/reservoir-fluid mixture with the modified aqueous EOS model, then repeats the calculation with the non-aqueous EOS model for the complementary water-in-hydrocarbon quantities.

This is strong evidence for a practical **dual-model / dual-calculation** workflow. It is not evidence that the current Water bot uses the frozen MPMC `SW92/corrected-original` thermodynamic profile unchanged: the current manual also documents later CO2 parameterization and EOS tuning steps. Therefore MPMC may reuse the **algorithm structure** while preserving a separate thermodynamic-profile identifier.

A public 2026 Søreide–Whitson refresh repository (`mwburgoyne/SW_Framework_Refresh`) independently documents the same engineering pattern as two independent full-mixture flashes, one with AQ water-gas BIPs and one with NA water-gas BIPs. This repository is used only as corroborating implementation evidence; no source code or parameter table is copied into MPMC.

### 2.3 Xu, Haynes and Stadtherr 2005: generalized asymmetric-model route

Secondary/generalized theory:

- G. Xu, W. D. Haynes and M. A. Stadtherr, *Reliable Phase Stability Analysis for Asymmetric Models*, Fluid Phase Equilibria 235 (2005) 152-165, DOI `10.1016/j.fluid.2005.06.016`.

For distinct liquid- and vapor-phase thermodynamic models, Xu et al. define the physical Gibbs surface as the lower envelope

```text
g(x) = min[g^A(x), g^B(x)]
```

and require candidate-model tangent-plane distances to use the **same feed Gibbs value and feed tangent plane**, determined from the lower model at the feed. They introduce a binary model variable to handle the nondifferentiable model choice and solve the resulting problem with interval methods.

This provides a sound general framework for a future AQ/NA asymmetric-Gibbs interpretation, but MPMC's current `test_pt_stability` is a finite local/multistart search and does not reproduce the interval-Newton global guarantee of Xu et al. A future MPMC profile based on this theory must retain `global_stability_proven=false` unless a separate globally reliable method is implemented and validated.

### 2.4 Secondary caution on two-BIP cubic-EOS use

Later literature has noted that using two phase-specific BIP sets in a cubic EOS can improve mutual-solubility predictions while also raising thermodynamic-consistency questions if the phase-specific models are combined without a clear model identity and equilibrium construction. This is an additional reason to keep the thermodynamic profile, phase-family identity, and equilibrium algorithm explicit rather than treating `kij_AQ`/`kij_NA` as an implementation detail.

## 3. Architectural requirements

The following requirements are normative for MPMC:

1. **Thermodynamic profile and equilibrium algorithm are separate identities.** A result must not use the string `SW92` alone to imply formulas, parameter version, phase-selection policy, and flash algorithm.
2. **Phase family is not a cubic-root index.** `SwPhaseFamily::{aqueous,nonaqueous}` remains explicit; a root index is only a local algebraic branch diagnostic inside a selected family.
3. **Flash owns phase-set decisions.** `thermodynamics` supplies phase properties; it does not perform stability search or select the accepted phase set.
4. **No composition cutoff selects AQ/NA.** Rules such as `x_water > 0.8` are not part of the MPMC SW92 contract unless a separately sourced and validated algorithm profile explicitly requires one.
5. **No silent model upgrade.** Current Whitson/Yan/refreshed correlations, tuned water critical properties, or later BIPs may not replace the frozen corrected-original SW92 inputs without a new thermodynamic-profile identifier and new validation.
6. **Fixed NaCl molality remains an external state coordinate.** Neither algorithm profile defined here implies conserved salt inventory, explicit ions, precipitation, reactions, or finite-rate mass transfer.
7. **Finite search is not a global proof.** Existing generic stability statuses retain their current semantics.

These requirements follow the same separation seen in mature open-source thermodynamic software: ThermoPack distinguishes phase/root semantics, Clapeyron exposes explicit flash-method objects, DWSIM separates `PropertyPackage` from `FlashAlgorithm`, and OPM flash solvers consume a separate `FluidSystem`. NeqSim also separates Søreide–Whitson system/phase classes, but current code contains composition-based AQ/NA routing in parts of its mixing-rule implementation; MPMC deliberately does not copy that heuristic.

## 4. Algorithm profile A — primary compatibility path

### 4.1 Identity

```text
SW92-equilibrium/dual-model-directional/PT-v1
```

Status:

```text
specified; not yet implemented
```

Scientific interpretation:

- thermodynamic properties come from one explicit SW92 thermodynamic profile, initially `SW92/corrected-original/PR76-base/NaCl-molality`;
- the **dual-model calculation structure** follows current Whitson practice;
- this profile is an MPMC operational contract and must not be described as a verbatim algorithm printed in the 1992 paper.

### 4.2 Inputs

Both directional runs must use the same immutable calculation snapshot:

```text
pressure_pa
temperature_k
ordered overall composition z
thermodynamic_profile_id
dataset_id + revision
ordered component ids
fixed NaCl molality [mol/kg H2O]
root options
family-local stability/split options
resource limits
```

The two runs must not be built from different parameter datasets, component orderings, model revisions, molality values, or PT states.

### 4.3 Aqueous-direction run

The complete family-local equilibrium calculation uses the **aqueous thermodynamic parameterization consistently for all property evaluations in that run**:

```text
SwPhaseFamily::aqueous
    -> family-local stability
    -> family-local phase split if required
    -> family-local final stability review
```

The current `Sw92FamilyStabilityEvaluator` is the validated stability provider for the first step. A future `Sw92FamilyVleEvaluator` may provide the split properties, but must remain fixed to the same aqueous family for the complete run.

### 4.4 Non-aqueous-direction run

The second calculation independently uses:

```text
SwPhaseFamily::nonaqueous
    -> family-local stability
    -> family-local phase split if required
    -> family-local final stability review
```

It uses the same external state and parameter snapshot as the aqueous-direction run but the non-aqueous family BIPs.

### 4.5 Combined result semantics

The top-level result is a **pair of directional equilibrium calculations**, not one joint thermodynamic phase set.

It must preserve both complete inner results and at least:

```text
equilibrium_algorithm_id
thermodynamic_profile_id
dataset_id + revision
component ordering
p, T, fixed molality
AQ-run status/diagnostics
NA-run status/diagnostics
AQ-run evaluation counts
NA-run evaluation counts
```

A top-level convenience status may report whether both directional calculations were numerically resolved, but it must not use `accepted_phase_set` terminology for the combined pair.

The following claims are forbidden for this profile unless independently established by a later contract:

- one common phase fraction `beta` assembled from the two runs;
- one combined material balance assembled from phases taken from different runs;
- cross-family fugacity equality;
- a global Gibbs minimum across AQ and NA models;
- a joint one/two/three-phase equilibrium state;
- a conserved salt balance.

This restriction is essential: two independently equilibrated family-local calculations cannot be spliced into the existing generic `PtPhaseSetResult` and called an accepted multiphase state.

### 4.6 Directional target extraction

Current Whitson practice extracts aqueous quantities from the aqueous-model flash and non-aqueous quantities from the non-aqueous-model flash. MPMC does not yet define a general target-phase identification rule for arbitrary multicomponent `z`.

A future implementation increment must therefore separately audit how to identify the target phase inside each family-local result. It may not identify the target solely from cubic-root order. If composition-based phase identification is used, the criterion, tie behavior, and validity domain must be explicit and tested.

### 4.7 Intended use

This is the **first implementation target** because it has the strongest direct connection to current Whitson engineering practice while preserving the already validated corrected-original thermodynamic profile.

Its purpose is directional mutual-solubility / water-property computation and compatibility validation. It is not the final architecture for a fully coupled three-phase compositional-flow closure.

## 5. Algorithm profile B — generalized asymmetric-Gibbs path

### 5.1 Identity

```text
SW92-equilibrium/asymmetric-gibbs-envelope/PT-v1
```

Status:

```text
reserved/research; not implemented
```

Scientific interpretation:

This profile treats `aqueous` and `nonaqueous` as distinct thermodynamic phase models and applies the asymmetric-model stability construction of Xu et al. It is a generalized equilibrium interpretation, not an assertion that Søreide and Whitson 1992 used this exact stability algorithm.

### 5.2 Preconditions

Before implementation, MPMC must establish and regression-test that the two family models are comparable under a common chemical-potential reference convention for the chosen thermodynamic profile. For the current SW92 implementation this requires, at minimum:

- identical ordered component identities and pure-component parameter sources;
- identical PT state and fixed molality coordinate;
- a common fugacity definition `f_i = x_i p phi_i`;
- family differences confined to explicitly recorded family-dependent mixture rules/BIPs;
- no hidden phase-dependent pure-component standard-state shift;
- independently checked family Gibbs/chemical-potential consistency for representative states.

Passing these checks establishes a **model-internal** asymmetric equilibrium construction; it does not establish experimental validity throughout the whole composition simplex.

### 5.3 Feed lower-envelope reference

For the same `p,T,z,molality`, determine the mechanically admissible minimum-Gibbs candidate inside each family using the validated fixed-family root-selection semantics.

Away from an unresolved family tie, define the feed reference from the lower family Gibbs surface. With a common ideal-gas/reference convention, the family difference at identical feed composition reduces to a residual contribution such as

```text
Delta g/(RT) = sum_i z_i [ln(phi_i^AQ) - ln(phi_i^NA)]
```

with numerically careful summation and a **roundoff-only comparison guard**.

The ordinary TPD tolerance must not be reused as a family-selection tolerance.

If the AQ and NA feed surfaces are numerically tied within the comparison guard, the lower envelope is nonsmooth at the reference composition. The current generic single-reference TPD wrapper is insufficient to claim a resolved reference there. The top-level state must remain `indeterminate` until a dedicated asymmetric-envelope treatment equivalent in semantics to the Xu binary-model formulation is implemented.

### 5.4 Common tangent and candidate searches

For a resolved feed reference, construct one common reduced chemical-potential/tangent vector for all active components:

```text
d_i = ln(z_i) + ln(phi_i^reference)
```

with the existing finite placeholder convention for `z_i == 0`.

Then search the candidate models **against the same reference**:

```text
test_pt_stability_against(..., d, aqueous_evaluator, ...)
test_pt_stability_against(..., d, nonaqueous_evaluator, ...)
```

Each result retains its family identity and its local cubic-root branch diagnostics.

Aggregation for the current finite-search implementation is conservative:

- any robust negative TPD witness in either family -> `unstable`;
- both family searches `no_instability_found` -> `no_instability_found`, with `global_stability_proven=false`;
- no negative witness, but either family search unresolved -> `indeterminate`.

Independent per-family resource budgets must be reported; one shared evaluation counter must not create AQ-first/NA-first order bias.

### 5.5 Joint phase split

A future phase split for this profile must solve one joint material-balance/equilibrium problem with explicit family labels on phase instances. It may admit, for example:

```text
one aqueous + one non-aqueous
one aqueous + two non-aqueous
multiple non-aqueous phase instances at different compositions
```

as allowed by the future solver capability and thermodynamic phase rule.

For every transferable component present in coexisting phases, the accepted joint state must enforce the selected profile's common chemical-potential/fugacity equilibrium and component material balance. Two directional results from profile A may be useful as initial estimates or validation references, but may not be stitched together as the joint solution.

The existing generic `PtPhaseSetResult` is suitable as a publication boundary because phase count is independent of status and multiple phase instances need not have unique universal root labels. The numerical joint solver itself is not defined by this document.

### 5.6 Global-optimality boundary

Xu et al. obtain deterministic/global reliability using interval analysis. MPMC does not currently implement that solver. Reusing the Xu **thermodynamic formulation** does not transfer its global mathematical guarantee to MPMC's finite local/multistart search.

Therefore this profile must continue to publish:

```text
global_stability_proven = false
```

until a separately audited globally reliable search is implemented.

## 6. Algorithm/profile matrix

| Property | `dual-model-directional/PT-v1` | `asymmetric-gibbs-envelope/PT-v1` |
| --- | --- | --- |
| Scientific basis | Current Whitson dual-EOS operational structure + SW92 phase-specific BIPs | Xu 2005 asymmetric-model stability theory + SW92 family models |
| First implementation priority | **Yes** | No; research/reserved |
| Uses fixed family per property evaluation | Yes | Yes |
| Runs AQ and NA calculations | Independently | Under one common Gibbs-envelope interpretation |
| Cross-family feed Gibbs comparison | No | Yes |
| Common tangent across families | No | Yes |
| Combined result is a joint accepted phase set | **No** | Intended in future, after full solver validation |
| Suitable final basis for joint three-phase closure | No | Potentially, after validation |
| Existing finite-search global proof | No | No |
| Salt inventory conservation | No | No |

## 7. Required result identity

Any future public SW92 equilibrium result must retain enough metadata to distinguish model from algorithm. At minimum:

```text
thermodynamic_profile_id
equilibrium_algorithm_id
algorithm_revision
dataset_id
parameter_revision
ordered component ids
pressure_pa
temperature_k
fixed NaCl molality
root options
stability/split options and resource budgets
algorithm-specific status and diagnostics
```

A serialized result that contains only `model = SW92` is insufficient.

## 8. Forbidden shortcuts

The following are explicitly outside the contract:

- selecting AQ/NA from cubic-root order;
- selecting AQ/NA using an undocumented water-mole-fraction threshold;
- evaluating `min(AQ,NA)` at each composition but discarding which model produced the value;
- using a TPD convergence tolerance as a physical family-switch tolerance;
- silently replacing corrected-original SW92 correlations with current Whitson/Yan/refreshed parameters;
- calling two directional family-local flashes a joint VLE/VLLE state;
- composing a three-phase result from two independent two-phase calculations;
- treating fixed molality as a conserved salt inventory;
- claiming Xu-style global stability from the existing finite multistart search.

## 9. Mature-software architecture references

These references are **non-normative architecture precedents**. No code is copied from them.

- ThermoPack phase flags distinguish liquid, vapor, minimum-Gibbs root, and unidentified single phase rather than treating an algebraic root index as a universal phase identity: `https://thermotools.github.io/thermopack/vcurrent/phase_flags.html`.
- Clapeyron.jl exposes explicit flash-method objects (`RRTPFlash`, `MichelsenTPFlash`, `MultiPhaseTPFlash`, `DETPFlash`) and separates two-phase methods from multiphase phase-search algorithms: `https://clapeyronthermo.github.io/Clapeyron.jl/dev/properties/flash/`.
- DWSIM exposes a separate `PropertyPackage` and `FlashAlgorithm`, with stability, phase identification, and PT flash operations on the algorithm side: `https://dwsim.org/api_help/html/T_DWSIM_Thermodynamics_PropertyPackages_Auxiliary_FlashAlgorithms_FlashAlgorithm.htm`.
- OPM `NcpFlash<Scalar, FluidSystem>` separates the fluid-property system from the flash constraint solver and enforces component fugacity equilibrium and material balance at the solver layer: `https://opm-project.org/apidoc/latest/opm-material/html/class_opm_1_1_ncp_flash.php`.
- NeqSim has explicit Søreide–Whitson system/phase/mixing-rule types. Some current mixing-rule paths classify aqueous behavior by water composition; that is an engineering heuristic, not adopted here: `https://github.com/equinor/neqsim/blob/006727adc7cf0350f172a737024d21fdb581ceef/src/main/java/neqsim/thermo/mixingrule/EosMixingRuleHandler.java`.

## 10. Validation gates before code

### 10.1 First code increment: one-family VLE provider

Before the top-level dual-model wrapper, implement and validate a `Sw92FamilyVleEvaluator` that is fixed to one family and composes:

```text
Sw92FamilyStabilityEvaluator
+ explicit requested low-/high-density candidate root evaluation
+ existing generic PT split kernel
```

The family must remain fixed through initial stability, split iteration, and final family-local stability review. This increment must not add cross-family logic.

Tests should reuse traceable SW92 binary data and independent high-precision calculations. Existing PR76 split behavior and tolerances remain unchanged.

### 10.2 Second code increment: dual-model directional wrapper

Only after both family-local VLE paths are validated, add the `dual-model-directional/PT-v1` orchestration that runs AQ and NA calculations from one immutable input snapshot and returns both results separately.

The acceptance test is structural and numerical reproducibility of the two directional calculations; it is **not** joint material-balance or cross-family fugacity equality.

Where current Whitson output is used for comparison, its thermodynamic parameterization/tuning must be matched explicitly. Current Whitson results must not be used as direct golden values for the frozen corrected-original MPMC thermodynamic profile unless all relevant model inputs are demonstrated to be identical.

### 10.3 Asymmetric-Gibbs research gate

Before implementing `asymmetric-gibbs-envelope/PT-v1`, add an independent audit/reference increment that verifies:

- common AQ/NA chemical-potential reference-state convention;
- cross-family feed Gibbs differences with independent high-precision references;
- imposed-common-tangent TPD values for both families;
- model-envelope near-tie and nonsmooth semantics;
- family-local numerical-failure aggregation;
- component permutation and runtime component-count invariance;
- symmetric resource accounting;
- unchanged `global_stability_proven=false` semantics.

Only after those gates should a joint AQ/NA split be designed.

## 11. Current decision

The project decision at this stage is:

1. retain `SW92/corrected-original/PR76-base/NaCl-molality` as the current thermodynamic profile;
2. retain `Sw92FamilyStabilityEvaluator` as the validated fixed-family stability building block;
3. implement **`SW92-equilibrium/dual-model-directional/PT-v1` first**, because it has the strongest operational traceability to current Whitson practice;
4. keep **`SW92-equilibrium/asymmetric-gibbs-envelope/PT-v1` reserved** as the scientifically motivated route toward a true joint phase-set model;
5. do not call either route simply "the SW92 flash" without recording the explicit algorithm profile.

This decision may be revised if stronger primary-source evidence becomes available, but any revision requires a new audit record and must preserve existing validated profile identities.