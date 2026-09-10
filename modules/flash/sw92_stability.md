# SW92 fixed-family PT stability adapter

## Scope

`Sw92FamilyStabilityEvaluator` connects the existing
`SW92/corrected-original/PR76-base/NaCl-molality` thermodynamics kernel to the
model-independent PT tangent-plane-distance search. One evaluator fixes:

- one validated `Sw92Phase<double>` snapshot;
- one NaCl molality in `mol NaCl / kg H2O`;
- exactly one `SwPhaseFamily::{aqueous,nonaqueous}`;
- one `Sw92RootOptions` value.

The family and molality remain fixed for every composition evaluated by that
search. This increment does **not** decide whether an unknown feed reference is
aqueous or nonaqueous, compare AQ and NA reference families, solve a phase split,
conserve salt inventory, or implement three-phase flash.

## Root selection inside one family

For each `(p,T,x)` requested by the generic stability driver, the adapter:

1. enumerates the selected family's Peng-Robinson cubic roots through
   `Sw92Phase::roots`;
2. retains only mechanically admissible roots with `H'(Z)>0`, equivalently the
   existing root diagnostic `slope_sign>0`;
3. evaluates `ln(phi)` for every retained root with the same family and fixed
   molality;
4. ranks those roots by `sum_i x_i ln(phi_i)`, for which ideal and component
   reference terms cancel at fixed `p,T,x,family`;
5. returns the lowest-Gibbs admissible root as `StabilityPhase`;
6. marks a same-family Gibbs near-tie `smooth=false` and rejects a selected root
   whose existing root-conditioning diagnostic is not derivative-valid.

`StabilityPhase::branch` remains the increasing-Z cubic root index **inside the
fixed family**. It never encodes aqueous/nonaqueous identity.

## Error and diagnostic mapping

Numerical property failures map onto the existing generic stability diagnostics:

| SW92 property condition | `StabilityPropertyIssue` |
| --- | --- |
| unresolved near-multiple roots | `root_topology` |
| cubic root iteration budget exhausted | `root_iteration_limit` |
| nonrepresentable root/property arithmetic | `root_range` |
| no mechanically admissible root | `no_admissible_branch` |
| selected root conditioning insufficient | `ill_conditioned_root` |
| nonfinite Gibbs comparison | `nonfinite_properties` |

Input/contract errors such as invalid family, invalid or out-of-declared-range
molality, wrong composition dimension, or invalid pressure/temperature continue
to throw; they are not relabeled as thermodynamic indeterminacy.

## Generic TPD bridge

`test_sw92_pt_family_stability` is a thin metadata-preserving wrapper around the
existing `test_pt_stability`. It records dataset/revision/component ordering,
model/phase conventions, fixed family, fixed molality, and root options. The TPD
search algorithm, tolerances, starts, finite-search semantics, and
`global_stability_proven=false` contract are unchanged.

The evaluator can also be passed directly to `test_pt_stability_against` when a
caller already owns a scientifically consistent common `log_activity` reference.
That capability does not itself define cross-family orchestration: deciding which
AQ/NA references must be tested and how their results form an initial phase set is
a separate audited increment.

## Model validity boundary

SW92 phase-family correlations are empirical. The generic TPD search may visit
compositions across the active simplex, including regions with limited or absent
experimental coverage. This adapter does not invent a water-fraction cutoff,
clip trial compositions, or infer family from composition/root order. Results are
therefore model-internal fixed-family numerical stability results, not an
experimental validity certificate or global stability proof.

## Verification

The focused stability integration includes:

- CO2/water nonaqueous three-root selection with the mechanically unstable middle
  root excluded and the minimum-Gibbs admissible root selected;
- CO2/water aqueous single-root routing at fixed brine molality;
- a fixed-family generic TPD negative-witness regression with an independent
  Decimal(80) implementation of SW92 correlations, PR mixing, the original-Z
  cubic, Gibbs root ranking, and TPD;
- phase-family/molality/root-budget/dimension failure semantics;
- component-order permutation invariance and public-header self-containment.

The reference uses the same traceable SW92 1992 corrected-original CO2/water data
already documented in `modules/thermodynamics/sw92.md`; it does not introduce a
new component database or experimental-validation claim.

## Equilibrium-algorithm boundary

A fixed-family stability result is not, by itself, a definition of how AQ and NA
parameterizations form one equilibrium calculation. The source audit now distinguishes
an evidence-aligned Whitson dual-model observable workflow from a generalized
Xu-style asymmetric-Gibbs equilibrium route. Their algorithm identities, publication
rules and implementation gates are defined in
[`sw92_equilibrium_algorithm.md`](sw92_equilibrium_algorithm.md).

In particular, this adapter must remain reusable by either route. No AQ/NA selection,
water-fraction heuristic, cross-family Gibbs policy or phase-split ownership is to be
moved into `Sw92FamilyStabilityEvaluator`.
