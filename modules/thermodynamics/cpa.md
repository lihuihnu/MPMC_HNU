# CPA thermodynamics baseline

## Selected profile

This module begins CPA support with one explicit model profile:

`CPA/SRK-physical/simplified-rdf-1999/explicit-site-pairs/v1`

It combines:

- the SRK physical contribution used by the original CPA development;
- the CPA association contribution based on Wertheim association theory;
- the simplified radial-distribution expression proposed in the later CPA literature, `g = 1/(1 - 1.9 eta)`, with `eta = b*rho/4`;
- explicit site-pair association parameters. No hidden cross-association combining rule is selected in this profile.

Primary model anchor: G. M. Kontogeorgis et al., *Cubic equation of state for calculation of phase equilibria in association systems*, Fluid Phase Equilibria 118 (1996) 27-59, DOI `10.1016/0378-3812(95)02843-9`.

The simplified CPA radial-distribution form is the 1999 form commonly used in subsequent CPA implementations and is kept explicit in the profile identity so it is not confused with the original Carnahan-Starling contact expression.

## Parameter contract

`cpa_parameters.hpp` owns an ordered runtime component snapshot and requires:

- `Tc` [K];
- `a0` [Pa m^6 mol^-2];
- `b` [m^3 mol^-1];
- `c1` [dimensionless];
- declared association site classes and multiplicities;
- explicit physical-part `kij` for every selected distinct component pair;
- explicit allowed association site-pair records with association energy `epsilon` [J mol^-1] and association volume `beta` [dimensionless].

Every numerical datum carries repository provenance. Synthetic data require `DataPolicy::allow_synthetic_tests`.

Absence of an association site-pair record means that interaction is prohibited by the configured association scheme. It does **not** mean that an unknown cross-association parameter is silently set to a physical zero. If a future model profile adopts Elliott, CR-1, or another combining rule, that rule must receive its own profile/convention and validation.

## Association kernel

For a supplied `T`, molar density `rho`, composition `x`, and configured parameter set, `solve_cpa_association(...)` solves the site-fraction equations

```text
X_Ai = 1 / (1 + rho * sum_j x_j * sum_B n_Bj X_Bj Delta_AiBj)
```

with

```text
Delta_AiBj = g(rho) * [exp(epsilon_AiBj/(R T)) - 1] * b_ij * beta_AiBj
b_ij       = (b_i + b_j)/2
g(rho)     = 1/(1 - 1.9 eta)
eta        = b_mix rho / 4
```

The current numerical method is a bounded damped fixed-point iteration. It reports iteration-limit, radial-distribution singularity, and nonrepresentable-strength states explicitly. It does not clip site fractions or substitute an unconverged association state into a phase property.

The structural regression includes a one-component/one-site configuration for which the site-fraction equation reduces to a quadratic. Production output is checked against that independent analytic solution.

## Density-state phase kernel

`evaluate_cpa_phase_at_density(...)` is deliberately a lower-level density-state primitive. Given `T`, `rho`, and `x`, it computes

```text
a_i(T) = a0_i [1 + c1_i (1 - sqrt(T/Tc_i))]^2

a_mix = sum_i sum_j x_i x_j sqrt(a_i a_j) (1-kij)
b_mix = sum_i x_i b_i
```

and the SRK physical pressure plus CPA association pressure:

```text
P_physical = R T rho/(1-b rho) - a rho^2/(1+b rho)

P_assoc = -0.5 R T rho
          * (1 + rho d ln(g)/d rho)
          * sum_i x_i sum_A n_Ai (1-X_Ai)
```

For the simplified radial distribution,

```text
rho d ln(g)/d rho = 1.9 eta/(1-1.9 eta)
```

A component set with no association sites therefore reduces exactly to the configured SRK physical term; this limit is regression-tested.

## Current capability boundary

This gate establishes only:

1. runtime ordered CPA parameter/provenance contract;
2. association site-fraction kernel;
3. density-state pressure/property primitive;
4. component-permutation and analytic structural regression.

It does **not** yet provide:

- a `p,T,x -> density/root` solver;
- fugacity coefficients / chemical potentials;
- CPA TPD stability;
- CPA two-phase or three-phase flash;
- implicit sensitivity;
- water/hydrocarbon physical parameter database;
- an implicit cross-association combining rule;
- a claim of experimental validation.

The next CPA gate should derive and validate the full phase chemical-potential/fugacity kernel and density/root solve against an independent reference before connecting CPA to the generic PT stability and max-three-phase backend.
