# SW92/Xu max2 traceable CH4/CO2/H2O ternary regression

## Purpose

This validation increment adds the first **three-component SW92/Xu max2 regression whose complete
parameter snapshot is literature-traceable without synthetic binary-interaction data**.

It validates the already implemented algorithm

```text
SW92-equilibrium/xu-asymmetric-gibbs/max2-logK-SSI-RR/common-tangent/v1
```

and does not change any production thermodynamic or flash algorithm.

The regression is deliberately described as a **traceable hybrid-parameter model regression**, not as
a direct ternary experimental validation. All model inputs have literature provenance, but the
chosen `CH4/CO2/H2O` state is a numerical regression point assembled from separately sourced model
parameters.

## Parameter provenance

### SW92 terms

Soreide and Whitson (1992), *Peng-Robinson predictions for hydrocarbons, CO2, N2, and H2S with pure
water and NaCl brine*, Fluid Phase Equilibria 77, 217-240,
DOI `10.1016/0378-3812(92)85105-H`, supplies the already audited corrected-original SW92 portion:

- CH4, CO2 and H2O critical properties and acentric factors from Table 3;
- the corrected aqueous CH4/H2O correlation;
- the corrected aqueous CO2/H2O correlation;
- non-aqueous CH4/H2O `k_ij = 0.4850` from Table 5;
- non-aqueous CO2/H2O `k_ij = 0.1896` from Table 5;
- water alpha and the existing corrected-original PR76-base conventions.

The repository source audit and supplied-PDF digest remain recorded in
`modules/thermodynamics/sw92.md`.

### CH4/CO2 non-water pair

SW92 does not provide a universal non-water/non-water BIP database. The current
`Sw92ParameterSet` therefore requires every selected non-water pair to be supplied explicitly with
provenance; it never infers zero or copies a water-pair rule.

For this validation snapshot the missing CH4/CO2 pair is supplied from:

S.-E. K. Fateen, M. M. Khalil and A. O. Elnabawy (2013), *Semi-empirical correlation for binary
interaction parameters of the Peng-Robinson equation of state with the van der Waals mixing rules
for the prediction of high-pressure vapor-liquid equilibrium*, Journal of Advanced Research 4,
137-145, DOI `10.1016/j.jare.2012.03.004`, PMCID `PMC4195455`.

Table 1 row 34 reports for methane/carbon dioxide:

```text
constant Peng-Robinson + van-der-Waals-mixing-rule kij = 0.0919
12 data sets, 110 points
```

Figure 4 in that paper compares the constant-`k_ij` calculation with experimental CH4/CO2 VLE at
250 K and 270 K. The present regression uses the tabled constant as a sourced **gas/gas PR mixing
parameter**; it does not claim that Figure 4 validates the present 350 K ternary water-containing
state.

Because this external parameter is not phase-family-specific, the same sourced scalar `0.0919` is
supplied explicitly to both the `aqueous_kij` and `nonaqueous_kij` fields for the CH4/CO2 pair. This
is an explicit validation-snapshot data choice, not a new SW92 inference rule. The source metadata
continues to identify Fateen et al., not SW92.

## Regression state

Ordered components:

```text
CH4, CO2, H2O
```

State:

```text
p = 10 MPa
T = 350 K
NaCl molality = 0 mol/kg H2O
z = [0.35 CH4, 0.15 CO2, 0.50 H2O]
```

The temperature/feed are numerical model-regression coordinates. No direct ternary experimental
measurement is asserted for this exact state.

## Independent Decimal(80) result

`tests/flash/sw92_asymmetric_max2_physical_ternary/reference_decimal.py` is stdlib-only and imports no
production code. It independently rebuilds:

- PR76 pure-component terms;
- corrected SW92 water alpha and CH4/H2O, CO2/H2O family interactions;
- the sourced CH4/CO2 `k_ij=0.0919` in both family matrices;
- classical three-component PR mixing;
- cubic roots and mechanically admissible minimum-Gibbs root selection;
- a five-variable two-phase Newton solve (three common reduced chemical potentials plus two
  independent component balances);
- feed/pair reduced Gibbs and AQ-vs-NA family comparisons.

The independent AQ+AQ pair is:

```text
water-rich phase
  x_CH4 = 0.000926595565135305376744149839868598720...
  x_CO2 = 0.00410222546705065601534780964593278797...
  x_H2O = 0.994971178967814038607908040514198613...

other phase
  x_CH4 = 0.694825377331009477598127996194254247...
  x_CO2 = 0.294122280631855616926095936126350533...
  x_H2O = 0.0110523420371349054757760676793952200...

other-phase mole fraction
  beta = 0.503060984696532100395670072167533754...
```

Compressibility factors:

```text
Z_water-rich = 0.0758151675933028182659305846173684704...
Z_other      = 0.860286072014983298778271289543848281...
```

Reduced Gibbs evidence:

```text
G_feed,AQ / RT = -2.28631878230116056890700781726502666...
G_feed,NA / RT = -1.59857965464664726953666024786125066...
G_AQ - G_NA    = -0.687739127654513299370347569403776002...

G_pair / RT    = -3.09595922989048202539397969359502534...
G_pair-G_feed  = -0.809640447589321456486971876329998672...
```

Selected-pair common reduced tangent:

```text
d_CH4 = -0.462255240302127022411452093427274697...
d_CO2 = -1.50007205423714639252586253581859718...
d_H2O = -5.41831817529833121734218416104537923...
```

At the two converged phase compositions the independent family gaps are respectively

```text
g_AQ-g_NA = -0.0158605870907434694828533831379943984...
g_AQ-g_NA = -0.00607565267345640010750747291986237355...
```

so the Xu lower envelope selects AQ for both phases at this state. As in the earlier binary
validation, this is a property of the explicitly defined Xu-style generalized lower-envelope profile;
it must not be relabeled as proof that both phases are physically "aqueous" in the original SW92
phase-assigned engineering sense.

## C++ acceptance checks

The C++ regression requires:

- the hybrid parameter snapshot to build under ordinary (non-synthetic) `DataPolicy`;
- the CH4/CO2 non-water BIP to retain Fateen literature provenance in both family slots;
- missing CH4/CO2 pair data, or a missing family value, to be rejected instead of defaulted;
- full `solve_sw92_xu_asymmetric_max2` acceptance with two phases under the declared finite search;
- component material balance and common-chemical-potential residuals to remain within the existing
  production thresholds;
- the phase compositions, fractions, Z values, Gibbs evidence, common tangent and family gaps to
  match the independent Decimal(80) anchors;
- complete component permutation to preserve the same physical-parameter regression after remapping
  indices.

Every accepted result still records

```text
global_stability_proven = false
```

because the production AQ/NA stability searches are finite multistart searches.

## Scope boundary

This regression closes a previous validation-data gap: the repository now has one runtime
three-component max2 case with no synthetic thermodynamic parameter.

It does **not** establish:

- direct experimental validation of the complete CH4/CO2/H2O ternary equilibrium at 350 K;
- a validity range for the Fateen constant outside the evidence in its source;
- equivalence between the Xu lower-envelope family choice and original SW92 physical phase labels;
- global phase stability or proof that a third phase cannot exist;
- salt-inventory conservation;
- SW92 derivatives or physics closure;
- a three-phase solver.

Those claims require separate evidence and implementation gates.
