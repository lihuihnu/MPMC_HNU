# SW92 one-family symmetric PT VLE baseline

## Scope

This increment implements Gate 1 from `sw92_equilibrium_algorithm.md`: a complete PT VLE run in
which **one and only one** SW92 phase family is used by every property evaluation in that run.

Public entry points:

- `Sw92FamilyVleEvaluator`
- `solve_sw92_pt_family_vle`
- `Sw92FamilyPtSplitResult`

One evaluator fixes one validated `Sw92Phase<double>` snapshot, one NaCl molality, one
`SwPhaseFamily::{aqueous,nonaqueous}`, and one `Sw92RootOptions` value. It is sequentially reusable
but owns scratch and is not concurrently callable.

This is an internally consistent **single-model VLE primitive**. It is not yet
`SW92-equilibrium/whitson-dual-model-observables/v1`, because that profile requires two separately
retained complete runs. It is also not `SW92-equilibrium/xu-asymmetric-gibbs/v1`, because no AQ/NA
lower envelope, cross-family common tangent, or joint asymmetric phase set is constructed.

## Provider contract

The evaluator supplies both generic flash provider roles without changing the generic algorithms.

### Stability call

```text
provider(p,T,w)
```

is delegated to `Sw92FamilyStabilityEvaluator`. It enumerates all mechanically admissible cubic
roots in the fixed family and returns the same-family minimum-Gibbs root. This preserves the
existing finite multistart TPD semantics and `global_stability_proven=false`.

### Split candidate call

```text
provider(p,T,w,PtPhaseRole)
```

uses the same model/family/molality but an explicitly requested density branch:

- `liquid_candidate`: lowest-Z mechanically admissible root;
- `vapor_candidate`: highest-Z mechanically admissible root.

These names are **requested numerical candidate roles**, not proof of aqueous/non-aqueous identity.
The fixed SW92 family remains a separate property of the evaluator/result. A single-root state can
therefore return the same algebraic branch for both requested roles; the generic split acceptance
logic still requires distinct converged phases and final stability.

Near-multiple, root-iteration, root-range, no-admissible-root and conditioning failures reuse the
existing `StabilityPropertyIssue` mapping. Input/parameter/domain errors remain contract errors and
are not relabeled as successful thermodynamic states.

## Full one-family calculation

`solve_sw92_pt_family_vle` is a metadata-preserving wrapper around `solve_pt_vle` using the same
`Sw92FamilyVleEvaluator` for both provider arguments:

```text
fixed-family feed TPD
        -> negative witness / stable single candidate
        -> material-balanced Rachford-Rice seed
        -> same-family two-phase fugacity iteration
        -> material-balance + phase-distinction checks
        -> same-family common-tangent final TPD review
```

No tolerance, Rachford-Rice rule, TPD rule, split iteration, or final acceptance rule is forked for
SW92.

An accepted two-phase result under this primitive is a valid phase set **for that one fixed family
model**. It must not be combined with a phase from an independently solved opposite-family run and
then presented as though both phases satisfied one common material balance.

## Fixed-molality semantics

NaCl molality remains the external state coordinate defined by
`SW92/corrected-original/PR76-base/NaCl-molality`. The same value is used for every feed, trial and
split composition in the run. No Na+/Cl- EOS components or conserved salt inventory are added.

## Verification

The focused regression is isolated from the established PR76 split tests and uses the existing
GitHub-hosted compiler matrix.

Independent `Decimal(80)` calculations use no production imports. They re-evaluate the corrected
SW92 alpha/BIP expressions, classical PR mixing, original-Z cubic roots, requested mechanical root
roles, and solve binary fugacity equality directly before obtaining the phase fraction from the
material balance.

Two traceable binary model anchors are used:

1. **CO2/H2O, nonaqueous family**, `p=3 MPa`, `T=340 K`, fresh water, feed CO2 mole fraction `0.7`.
   This exercises the Table 5 nonaqueous CO2/water BIP and a three-root intermediate state.
2. **CH4/H2O, aqueous family**, `p=10 MPa`, `T=350 K`, NaCl molality `1 mol/kg H2O`, feed methane
   mole fraction `0.5`. This exercises the corrected hydrocarbon aqueous correlation with nonzero
   molality.

The regression checks phase compositions, phase fraction, Z values, material balance, fugacity
residual, final same-family stability, component permutation, fixed family/molality metadata,
root-budget failure semantics, requested root roles, and public-header self-containment.

These anchors are independent **model numerical references**, not experimental validation and not
a global phase-stability proof.

## Next gate

Only after both fixed-family calculations are individually reliable should the project implement
`SW92-equilibrium/whitson-dual-model-observables/v1`: run AQ and NA calculations independently,
retain both full results, identify the requested application observables under an audited phase
labeling rule, and never fabricate a cross-run accepted phase set.
