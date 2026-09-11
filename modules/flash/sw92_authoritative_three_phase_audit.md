# SW92 Profile-C authoritative three-phase completion audit

## Purpose

This audit defines what must be true before MPMC_HNU may describe the current
SW92 Profile-C PT implementation as a complete three-phase flash rather than a
locally-closed candidate generator.

The frontend is frozen while these gates are completed.

## Current baseline

`solve_sw92_phase_assigned_pt(...)` already orchestrates:

- no-W one/two-NA candidates;
- `W(AQ)+H(NA)` C1;
- H-side NA phase-addition evidence C2a1;
- `W(AQ)+H0(NA)+H1(NA)` C2b.1;
- W-present final H-multiplicity review C2b.2.

The existing top-level result deliberately retains:

```text
global_stability_proven = false
accepted_phase_set_published = false
morphology_resolved = false
```

and the only full three-phase numerical golden is currently a synthetic
CH4/CO2/H2O structural fixture. That is sufficient for equation/routing
validation, but not for physical three-phase validation or authoritative
publication.

## Completion definition

For the current fixed-molality SW92 Profile-C scope, `complete three-phase
flash` means all of the following.

1. **One material inventory.** Every retained phase uses the same ordered EOS
   component inventory and reconstructs the normalized feed within the declared
   absolute/relative balance tolerances.
2. **Three positive phase fractions for an interior state.** No clipping or
   hidden phase deletion is permitted.
3. **Common component chemical potentials.** `W`, `H0`, and `H1` satisfy the
   same reduced chemical potentials/fugacities within the existing numerical
   gate.
4. **Family assignment is explicit.** `W -> AQ`; both H instances -> NA.
   Family identity is not a liquid/vapor classifier.
5. **Same-family root selection remains physical.** Every phase property call
   uses the mechanically admissible minimum-Gibbs root already established by
   the SW92 family evaluators.
6. **H-slot symmetry.** Exchanging H0/H1 cannot change the physical state;
   canonicalization is representation only.
7. **Role-constrained final review.** A retained W-present three-phase state
   must pass the finite NA phase-addition review without inventing unrestricted
   AQ-vs-NA lower-envelope competition.
8. **Rival topology routing is executable.** W disappearance and H
   disappearance route to re-solved neighboring topologies; an endpoint is not
   silently accepted merely because the three-phase Newton iteration reached a
   small fraction.
9. **Physical three-phase regression exists.** At least one literature-defined
   water/oil/gas system using Soreide-Whitson-style PR modifications must be
   reproduced independently. Published data may be used as a topology oracle;
   numerical golden values must match the exact MPMC_HNU model profile rather
   than being copied from a different EOS variant.
10. **Publication semantics are explicit.** A future accepted result may be an
    authoritative result under the declared finite-search/topology contract
    while still stating `global_stability_proven=false`. Finite numerical search
    is not a mathematical proof of a global minimum.
11. **No fake morphology.** Until an independently validated morphology
    classifier exists, hydrocarbon phases remain `nonaqueous_unclassified`.
    This does not invalidate their material/fugacity equilibrium solution.
12. **Fixed-molality scope is named.** NaCl molality is a prescribed model
    parameter in the current profile, not an independently conserved salt
    component. Salt-inventory coupling is a later thermodynamic model extension
    and must not be implied by this flash result.

## Physical regression anchor

Primary source:

- E. Mortezazadeh and M. R. Rasaei, *A robust procedure for three-phase
  equilibrium calculations of water-hydrocarbon systems using cubic equations
  of state*, Fluid Phase Equilibria 450 (2017) 160-174,
  DOI `10.1016/j.fluid.2017.07.007`.

The paper explicitly bases its water-hydrocarbon calculation on the
Soreide-Whitson PR modifications, supplies pure properties for synthetic samples
4-8, supplies sample feed compositions, states that non-water BIPs are assumed
zero, and publishes the gas-condensate Sample-6 phase envelope.

### Selected state

Use gas-condensate **Sample 6** at:

```text
P = 10 MPa
T = 350 K
NaCl molality = 0 mol/kg H2O
```

Figure 7 places this point inside the published oil+gas+water region and far
from the plotted phase boundaries. The rounded Table-5 feed is:

```text
H2O  0.5085
C1   0.4249
C2   0.0214
C3   0.0112
C4   0.0085
C5   0.0040
C6   0.0029
C7+  0.0185
```

The printed values sum to `0.9999`. Test construction must therefore perform an
**explicit, documented source-data normalization by the printed sum**; this is a
traceable conversion of rounded literature data, not production-input repair.

Table-4 `Pc`, `Tc`, and `omega` values are used as published. All non-water
BIPs are zero because the source explicitly makes that assumption. Water/non-
water NA constants and AQ correlations follow the SW92 correlation contract.

### Model-profile difference that must remain visible

Mortezazadeh-Rasaei Appendix A uses a later piecewise Peng-Robinson `kappa`
correlation for `omega > 0.49`. MPMC_HNU's declared
`SW92/corrected-original/PR76-base` profile intentionally uses the original PR76
quadratic `kappa` for every non-water component. Sample-6 `C7+` has
`omega=0.50879`, so exact numerical equality to the paper's plotted curves is
not an admissible test requirement.

Accordingly:

- the **published phase envelope is the independent physical topology oracle**;
- an independent stdlib-only Decimal(80) solver using the exact MPMC_HNU
  corrected-original profile generates the numerical phase-composition,
  fraction, Z and common-chemical-potential golden;
- the C++ top-level solver must reproduce that independent golden without
  changing formulas, source data or tolerances.

## Physical-anchor acceptance checks

The Sample-6 state is accepted only if the independent oracle and C++ solver
both obtain an interior `W+H0+H1` state with:

- all three phase fractions positive and greater than the existing disappearance
  threshold;
- W water-richer than both final H phases;
- H0/H1 compositionally distinct;
- feed reconstruction within the existing mass tolerances;
- common reduced chemical potentials within the existing C2b.1 tolerance;
- generalized-RR residual within the existing gate;
- deterministic H-slot canonicalization;
- C2b.2 W-present H-multiplicity review locally closed;
- runtime component permutation invariance;
- no claim that H0/H1 are liquid/vapor solely from their Z values.

## Implementation order

1. Add the independent physical Sample-6 Decimal(80) oracle and a focused C++
   top-level regression.
2. Fix only demonstrated equation/orchestration defects exposed by that physical
   case; never relax tolerances or change the source anchor to obtain a pass.
3. Add explicit re-solve routing for three-phase disappearance boundaries if any
   route is currently status-only.
4. Introduce a separate authoritative phase-set publication adapter only after
   the physical anchor and boundary routes pass. Keep
   `global_stability_proven=false` and hydrocarbon morphology unresolved.
5. Expand maximum-three-phase no-W topology only if a regression demonstrates
   that the current two-NA no-W cap can hide a physically admissible third NA
   phase in the supported model scope.

Frontend/service work remains frozen until these gates are complete.
