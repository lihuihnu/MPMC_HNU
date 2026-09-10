# SW92 Whitson dual-model observable orchestration

## Scope

This increment implements the documented algorithm profile
`SW92-equilibrium/whitson-dual-model-observables/v1` on top of the reviewed
`SW92-equilibrium/fixed-family-vle-primitive/v1` building block.

The profile is an **observable calculation composed of two independent equilibrium runs**:

1. one complete aqueous-family SW92 PT VLE run;
2. one complete nonaqueous-family SW92 PT VLE run.

It is deliberately **not** an AQ/NA joint thermodynamic phase split. The two runs never share a
Rachford-Rice material balance, phase fractions, a common tangent plane, or a joint Gibbs minimum.
The result type exposes both complete family-run results separately and provides no conversion to
an accepted `PtPhaseSetResult` containing one phase from each run.

Public entry point:

- `solve_sw92_whitson_dual_model_observables`

Public result/option types:

- `Sw92WhitsonDualModelResult`
- `Sw92DualModelOptions`
- `Sw92DualModelStarts`
- `Sw92DualModelObservables`

## Evidence and compatibility boundary

The thermodynamic model remains the repository's strict
`SW92/corrected-original/PR76-base/NaCl-molality` profile, sourced from Søreide and Whitson (1992)
plus the authors' errata already recorded by `modules/thermodynamics/sw92.md`.

The orchestration pattern is supported by current engineering evidence:

- the whitson+ Water Bot manual states that separate aqueous and non-aqueous EOS models are
  developed, performs an aqueous-model flash, and repeats the calculation with the non-aqueous
  EOS model: <https://manual.whitson.com/methods/water-bot/>;
- the public 2026 Søreide-Whitson framework refresh describes two independent full-mixture
  Rachford-Rice flashes, taking the aqueous result from the AQ-BIP pass and the non-aqueous result
  from the NA-BIP pass:
  <https://github.com/mwburgoyne/SW_Framework_Refresh/blob/main/shared/vle_engine/_lib_vle_engine.py>.

Those current sources include tuning, refreshed correlations and application-specific features
that are **not** imported into MPMC_HNU's corrected-original thermodynamic profile. This increment
uses the documented algorithmic pattern only; it copies no third-party implementation and changes
no SW92 formula or parameter.

## Same-snapshot invariant

The orchestration accepts one `Sw92Phase<double>` and constructs both family evaluators internally.
Therefore AQ and NA passes necessarily share:

- dataset and revision;
- ordered component snapshot;
- pure-component parameters;
- thermodynamic profile and phase convention;
- pressure, temperature, feed and fixed NaCl molality.

Only the explicit `SwPhaseFamily` changes. Each pass may have independent root, stability, split and
resource options, as required by the algorithm contract.

Each nested `Sw92FamilyPtSplitResult` keeps the Gate-1 identity
`SW92-equilibrium/fixed-family-vle-primitive/v1`. The outer result separately records
`SW92-equilibrium/whitson-dual-model-observables/v1`.

## Binary target-phase extraction

The existing equilibrium-algorithm contract authorizes automatic physical target-phase extraction
only for narrowly defined binary compatibility regressions. For an accepted binary two-phase run:

- the AQ observable is the **water-richer** phase from the AQ run;
- the NA observable is the **water-poorer** phase from the NA run.

No fixed water-fraction threshold is used. The two water mole fractions are ordered directly; if
their difference is unresolved at a machine-roundoff guard, the label is indeterminate. The stored
`source_role` only records whether the selected composition came from the generic lowest-Z
`liquid_candidate` or highest-Z `vapor_candidate`; it does not redefine those numerical roles as
physical AQ/NA identities.

For more than two components, both complete model runs are still performed and retained, but this
increment does **not** automatically extract physical target phases. A reusable multicomponent
phase-label rule remains a separate audit gate. In particular, no `x_water > constant` heuristic
is introduced.

## Cross-model equilibrium ratio

When binary target phases are available, the result reports componentwise

```text
cross_model_equilibrium_ratio_i = x_i^(NA target) / x_i^(AQ target)
```

The name is intentionally explicit. This ratio is **not** an ordinary `K` from one Rachford-Rice
problem and must not be used with either run's phase fraction to claim a joint material balance.
A zero/nonrepresentable denominator or overflow leaves the observable result indeterminate.

## Status and failure semantics

`Sw92DualModelStatus` describes whether the compatibility observables can be published. Regardless
of outer status, callers can inspect both nested family runs and their full initial stability,
split attempts, final stability, resource counters and diagnostics.

Input/parameter/domain errors retain existing exception semantics. Numerical inability of one
family run to produce an accepted two-phase target is not converted into success because the other
family run succeeds. Finite TPD searches in both runs retain `global_stability_proven=false`.

## Verification

The focused regression uses the same traceable Søreide-Whitson Table 3/Table 5 and corrected AQ
correlations as the existing family-VLE tests, but independently solves **both** family models at
the same state.

Two binary compatibility states are used:

1. CO2/H2O at 3 MPa, 340 K, fresh water, feed CO2 mole fraction 0.7;
2. CH4/H2O at 10 MPa, 350 K, NaCl molality 1 mol/kg H2O, feed methane mole fraction 0.5.

A stdlib-only `Decimal(80)` oracle independently recomputes both family VLE solutions and the
cross-model ratios. C++ regressions check both nested accepted runs, target extraction, ratio
values, component permutation, algorithm metadata, family-specific failure propagation and public
header self-containment.

These are model numerical anchors, not experimental validation and not evidence of a joint AQ/NA
Gibbs equilibrium.

## Remaining boundary

This profile is a compatibility/observable route. It does not satisfy the project's stronger
`mutual_solubility` goal requiring all coexisting phases to obey one component inventory and one
set of equilibrium constraints. If that joint result is required, the next scientific gate is the
separately identified `SW92-equilibrium/xu-asymmetric-gibbs/v1` reference-state/domain audit and
joint solver; Profile-A outputs may be used as comparison evidence or initial guesses only.
