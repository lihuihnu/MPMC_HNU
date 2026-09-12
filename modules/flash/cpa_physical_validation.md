# CPA associating physical validation

## Scope

This validation gate tests the already-implemented CPA PT flash path against a traceable **associating** binary VLE system. It does not fit parameters, change production equations, relax flash/stability tolerances, or claim a three-phase physical oracle.

The validated system is methanol + water at `T = 333.15 K` using the repository profile

```text
CPA/SRK-physical/simplified-rdf-1.9eta/explicit-site-pairs/v1
```

with methanol represented by the 2B site scheme and water by the 4C site scheme.

## Literature chain

### CPA pure parameters and association schemes

G. M. Kontogeorgis et al., *Solvation Phenomena in Association Theories with Applications to Oil & Gas and Chemical Industries*, Oil & Gas Science and Technology 63 (2008) 305–319, DOI `10.2516/ogst:2008025`.

The test fixture transcribes the published CPA pure parameters for methanol and water and records the original table/units in every `SourceKind::literature` datum. Unit conversions are explicit:

- `a0`: `bar L^2 mol^-2 * 0.1 -> Pa m^6 mol^-2`;
- `b`: `L mol^-1 * 1e-3 -> m^3 mol^-1`;
- association energy: `bar L mol^-1 * 100 -> J mol^-1`.

The site schemes are explicit:

```text
methanol: H(1), e(1)        # 2B
water:    H(2), e(2)        # 4C
```

### Cross association

The CR-1 combining rule is **not** implemented as a hidden production rule. The fixture evaluates the published rule once and stores the resulting cross-site parameters as ordinary explicit `CpaAssociationPairInput` records:

```text
epsilon_cross = (epsilon_methanol + epsilon_water)/2
beta_cross    = sqrt(beta_methanol * beta_water)
```

Only complementary donor/acceptor cross pairs are configured. Absence of a site-pair record still means association is disabled for that pair, exactly as required by `CpaParameterSet`.

### Binary interaction parameter

G. K. Folas, *Modeling of Complex Mixtures Containing Hydrogen Bonding Molecules*, PhD thesis, Technical University of Denmark. Table 2.1 reports for water–methanol at `333.15 K`, CR-1:

```text
k12 = -0.055
Delta P = 0.6 %
Delta y * 100 = 0.8
```

The table identifies the experimental source as Kurihara et al. (1995). The repository does not refit `k12`.

### Experimental P-x-y data

K. Kurihara, T. Minoura, K. Takeda, K. Kojima, *Isothermal Vapor-Liquid Equilibria for Methanol + Ethanol + Water, Methanol + Water, and Ethanol + Water*, J. Chem. Eng. Data 40 (1995) 679–684, DOI `10.1021/je00019a033`.

Five interior 333.15 K methanol–water points are retained as the focused regression set:

| P / kPa | x(MeOH), liquid | y(MeOH), vapor |
| ---: | ---: | ---: |
| 39.223 | 0.1686 | 0.5714 |
| 48.852 | 0.3039 | 0.6943 |
| 56.652 | 0.4461 | 0.7742 |
| 63.998 | 0.6044 | 0.8383 |
| 72.832 | 0.7776 | 0.9141 |

The experimental compositions are observations, not exact model equilibrium points.

## Flash-validation semantics

For each literature point the test constructs a material-balanced feed containing 90 mol% of the reported liquid composition and 10 mol% of the reported vapor composition. The experimental phase compositions are supplied only as **initial split/stability starts**.

They are intentionally **not** supplied as extra final-common-tangent starts. `solve_pt_vle(...)` always appends its own converged liquid/vapor compositions to final stability. This separation is important: experimental measurement/model mismatch must be measured as accuracy error, not promoted into a topology-proof requirement.

Each accepted point must pass the unchanged production gates:

1. feed stability produces split evidence;
2. RR/log-K equations converge;
3. material balance satisfies the existing absolute/relative tolerances;
4. log-fugacity equality satisfies the existing `1e-11` maximum residual tolerance;
5. final common-tangent review from model-owned converged phases reports no sampled instability;
6. the accepted liquid root is re-evaluated and has a converged, nontrivial association state.

No TPD, stationarity, RR, fugacity, material-balance or phase-fraction acceptance tolerance is changed by this validation.

The physical-validation root search uses tighter numerical pressure-root tolerances than the standalone CPA PT default so that inner root error remains below the unchanged outer flash residual gates. This is test configuration only; no production default is changed.

## Frozen accuracy envelope

On the five retained points the hosted reference run gives approximately:

```text
mean |Delta x(MeOH)| = 0.00401
mean |Delta y(MeOH)| = 0.00666
max  |Delta x(MeOH)| = 0.00737
max  |Delta y(MeOH)| = 0.00937
```

The regression therefore freezes a stricter envelope than the initial audit threshold:

```text
mean |Delta x| <= 0.01
mean |Delta y| <= 0.01
max  |Delta x| <= 0.015
max  |Delta y| <= 0.015
```

The vapor mean threshold is consistent with the published CR-1 correlation scale (`Delta y * 100 = 0.8`) while still leaving room for platform-level floating-point differences.

Ordered-component permutation is independently required to preserve the physical solution.

## Explicit validation boundary

This gate establishes a traceable **associating two-phase CPA VLE** validation. It does not establish a physical three-phase CPA oracle.

Published CPA work exists for water/alcohol/alkane and other VLLE systems, but this repository will not create a three-phase physical regression until one source chain supplies, without guessing or incompatible model mixing:

- component identities and state (`p`, `T`, overall loading);
- all pure CPA parameters under the same CPA profile;
- association scheme and every required self/cross association rule or parameter;
- every binary physical interaction parameter used by the calculation;
- all three experimental phase compositions (and preferably phase fractions);
- enough methodological detail to map the literature equations onto the repository convention.

Until that gate is met, CPA max3 physical topology remains structurally validated with explicit synthetic fixtures, while the association equations and the complete two-phase path now have traceable physical validation.

`global_stability_proven` remains false: finite multistart TPD searches are not mathematical global proofs.
