# CPA associating physical validation

## Scope

This validation gate tests the implemented CPA PT flash path against a traceable **associating** binary VLE system. It does not fit parameters, relax flash/stability acceptance tolerances, or claim a three-phase physical oracle.

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

They are intentionally **not** supplied as extra final-common-tangent starts. The
final review instead receives the normalized feed as a model-input off-tangent
start, and `solve_pt_vle(...)` appends its own converged liquid/vapor
compositions. This separation is important: experimental measurement/model
mismatch must be measured as accuracy error, not promoted into a topology-proof
requirement. The feed trial is required to iterate to the unchanged stationarity
gate, so the review cannot silently degenerate into evaluating only its two
already-converged tangent points.

Each accepted point must pass the unchanged production gates:

1. feed stability produces split evidence;
2. RR/log-K equations converge;
3. material balance satisfies the existing absolute/relative tolerances;
4. log-fugacity equality satisfies the existing `1e-11` maximum residual tolerance;
5. final common-tangent review from the feed plus model-owned converged phases reports no sampled instability;
6. the accepted liquid root is re-evaluated and has a converged, nontrivial association state.

No TPD, stationarity, RR, fugacity, material-balance or phase-fraction acceptance tolerance is relaxed by this validation.

### 56.652 kPa final-TPD audit

The original validation also supplied both observed phase compositions as final
TPD starts. At 56.652 kPa the vapor-side trial reached `TPD=-5.06e-11`, which is not below
the existing `1e-10` negative-TPD gate, but the log-descent search terminated in
`line_search_failed` with stationarity about `3.66e-6`. The conservative
`indeterminate` publication was therefore correct: this was neither a robust
negative witness nor a completed stationary trial. It was not caused by the
global evaluation budget, root topology, or a property failure.

The regression now keeps observations out of topology search and uses the
normalized feed as the explicit off-tangent start. At 56.652 kPa all three final
trials (feed, converged liquid, converged vapor) must be stationary, the lowest
sampled TPD must remain above the unchanged effective negative gate, and the
existing fugacity and material-balance thresholds are asserted explicitly.

## CI tiers

- Clang Release and MSVC Release run the complete five-point physical/accuracy
  regression, ordered-component permutation, the focused final-TPD case,
  provenance, default precision, and public-header checks.
- GCC Debug with ASan/UBSan runs the representative 56.652 kPa production path,
  provenance/default precision, and public-header checks. It intentionally does
  not repeat the five-point accuracy sweep or the two full permutation solves;
  those are physical regressions rather than additional sanitizer safety paths.

## CPA flash-specific density-root default

The audit exposed a numerical-layer mismatch between the standalone CPA PT root default and the stricter outer flash equations. The thermodynamics-layer `CpaPtOptions{}` remains unchanged; `CpaVleEvaluator` now owns a **flash-specific** default:

```text
scan_intervals = 512                 # unchanged from standalone default
pressure_absolute_tolerance_pa = 3e-6
pressure_relative_tolerance    = 1e-12
```

The reason is not that a smaller root residual is always better. A focused physical probe showed:

- the old standalone root default (`1e-4 Pa`, `1e-10` relative) can close a representative physical VLE point, but leaves the resulting fugacity residual close to the `1e-11` outer gate and is unnecessarily expensive;
- an over-tightened `1e-7 Pa` flash candidate produced an indeterminate path at the low-pressure literature point under the same bounded finite root search;
- `3e-6 Pa / 1e-12` closes the full five-point associating literature set while preserving the original 512 scan intervals and unchanged outer flash/stability tolerances.

Therefore the validated value is frozen as a **numerical flash default**, not as a physical model parameter and not as a modification to standalone CPA phase-property semantics. Callers may still supply explicit `CpaPtOptions` when a different audited numerical regime is required.

## Frozen accuracy envelope

On the five retained points the hosted reference run gives approximately:

```text
mean |Delta x(MeOH)| = 0.00401
mean |Delta y(MeOH)| = 0.00666
max  |Delta x(MeOH)| = 0.00737
max  |Delta y(MeOH)| = 0.00937
```

The regression freezes:

```text
mean |Delta x| <= 0.01
mean |Delta y| <= 0.01
max  |Delta x| <= 0.015
max  |Delta y| <= 0.015
```

The vapor mean threshold is consistent with the published CR-1 correlation scale (`Delta y * 100 = 0.8`) while retaining a small cross-platform numerical margin.

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
