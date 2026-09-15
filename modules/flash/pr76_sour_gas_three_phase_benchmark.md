# PR76 sour-gas three-phase literature benchmark

## Scope

This benchmark replaces the previous PR76 max-three-phase validation gap in which the
full `2 -> 3` orchestration had only an explicitly synthetic structural three-phase
fixture. It adds a non-synthetic, literature-defined acid-gas mixture and a
repository-convention PR76 numerical reference.

It does **not** change any production EOS, root selection, TPD search, generalized
Rachford-Rice equation, convergence tolerance, phase-fraction gate, or max3 routing.
It is a validation increment only.

The benchmark is deliberately described as a **literature-defined engineering/model
benchmark**, not experimental validation and not a mathematical global-stability
certificate.

## Literature state

Primary benchmark source:

Z. Li and A. Firoozabadi,
“General Strategy for Stability Testing and Phase-Split Calculation in Two and Three
Phases,” *SPE Journal* 17(4), 1096–1107 (2012),
DOI `10.2118/129844-PA`.

The paper states that stability testing and phase-split calculations use the
Peng–Robinson EOS. Its Table 3 gives the acid-gas component properties and the
nonzero binary interaction coefficients. Table 9 selects the acid-gas condition

- `T = 178.8 K`,
- `P = 20 bar`,
- overall `n(CO2) = 0.5`,

and reports an actual three-phase split calculation for that point. Figure 1 places
the same acid-gas P–Z system in explicit `V-L1-L2` regions.

The recent independent algorithm study

J. Heringer et al.,
“New initialization procedures from phase stability testing in three-phase flash
calculations for CO2-hydrocarbon mixtures,”
*Fluid Phase Equilibria* 604 (2026) 114653,
DOI `10.1016/j.fluid.2025.114653`

reuses this sour-gas family. At `20 bar` and `178.8 K` it again reports a
three-phase interval and discusses the multiple TPD minima found on the feed
composition. That source is used only as topology/algorithm cross-evidence; the
numerical golden below is regenerated independently from the 2012 parameter table.

An earlier parameter source used by the same benchmark family is

K. B. Haugen, A. Firoozabadi and L. Sun,
“Efficient and robust three-phase split computations,”
*AIChE Journal* 57(9), 2555–2565 (2011),
DOI `10.1002/aic.12452`.

It explicitly describes Peng–Robinson with van-der-Waals mixing for its three-phase
examples and lists a near-identical rounded sour-gas parameter table. The present
reference nevertheless uses the more precise Table-3 numbers from Li–Firoozabadi
2012.

## Parameter snapshot

Ordered canonical components are

```text
CO2, N2, H2S, C1, C2, C3
```

with Table-3 values:

| component | Tc / K | Pc / bar | omega |
| --- | ---: | ---: | ---: |
| CO2 | 304.211 | 73.819 | 0.225 |
| N2 | 126.2 | 33.9 | 0.039 |
| H2S | 373.2 | 89.4 | 0.081 |
| C1 | 190.564 | 45.992 | 0.01141 |
| C2 | 305.322 | 48.718 | 0.10574 |
| C3 | 369.825 | 42.462 | 0.15813 |

Table 3 is explicitly headed as listing **nonzero** binary interaction
coefficients. Therefore every selected pair omitted from that list is represented
as an explicit zero in the repository parameter contract; missing data are never
silently treated as zero.

Nonzero symmetric constant `kij` values:

```text
CO2-N2  = -0.020
CO2-H2S =  0.120
N2-H2S  =  0.200
CO2-C1  =  0.125
N2-C1   =  0.031
H2S-C1  =  0.100
CO2-C2  =  0.135
N2-C2   =  0.042
H2S-C2  =  0.080
CO2-C3  =  0.150
N2-C3   =  0.091
H2S-C3  =  0.080
```

Hydrocarbon–hydrocarbon pairs not listed in the table are recorded explicitly as
zero.

The C++ fixture keeps source provenance on every pure and binary record. Critical
pressures are the only unit conversion: `bar * 100000 -> Pa`.

## Overall composition reconstruction

Table 9 publishes `n(CO2)=0.5` for a P–Z condition rather than printing all six
mole fractions again. Table 3 publishes the acid-gas starting vector.

For reproducibility, this benchmark defines the full P–Z composition explicitly:

1. set `z_CO2 = 0.5`;
2. keep the relative ratios of the five non-CO2 Table-3 entries;
3. normalize those five entries to the remaining `0.5`.

Thus, if `n_i^0` denotes the Table-3 initial entry,

```text
z_i = 0.5 * n_i^0 / sum_{j != CO2}(n_j^0),   i != CO2.
```

The resulting vector is

```text
CO2 = 0.5
N2  = 0.1195792770100074885969092518210906
H2S = 0.03346041255361154605487099189869971
C1  = 0.1167540336306079379127238069303560
C2  = 0.1797093062836135883994826060317244
C3  = 0.05049697052215943903601334331812921
```

This reconstruction is intentionally documented as a **derived P–Z composition
line**. The six numbers above are not claimed to be a row printed verbatim in the
paper.

## Independent PR76 numerical reference

`reference_sour_gas_decimal.py` uses only Python's standard-library `Decimal`.
It imports no production code and does not use C++ outputs.

The numerical model is fixed to the repository convention:

```text
PR76/printed-coefficients/R-SI-2019/exact-sqrt2/PT-v1
```

Specifically:

- `R = 8.31446261815324 J/(mol K)`;
- `a_c = 0.45724 (R Tc)^2/Pc`;
- `b = 0.07780 R Tc/Pc`;
- `kappa = 0.37464 + 1.54226 omega - 0.26992 omega^2`;
- classical symmetric `a_ij = sqrt(a_i a_j) (1-kij)`;
- exact `sqrt(2)` fugacity expression.

The reference solve uses 20 coupled unknowns:

- five log-ratio composition coordinates for each of three phases;
- two independent phase-fraction log-ratios;
- three positive compressibility factors.

It enforces

- two complete six-component chemical-potential equality sets;
- five independent overall material balances;
- one PR cubic equation per phase.

The Decimal Newton Jacobian is regenerated by symmetric high-precision numerical
differentiation inside the reference program. The initialization is intentionally
rounded to only about three decimal places. The emitted golden values therefore
come from the nonlinear solve, not from copying a production result.

Default CI recomputes the reference at Decimal(80). The sanitizer configuration
also recomputes Decimal(96); both must render the exact same fixed-digit header.

The independent result is:

```text
phase fractions =
  0.1200779322068530430482637509731291
  0.4839913617133007830070106657613211
  0.3959307060798461739447255832655498

Z =
  0.8681744728389315603807464806454163
  0.05774692992299103461273954832796937
  0.04735887844747790632380732337695566
```

The complete phase compositions are stored in `sour_gas_references.hpp`.

This three-phase numerical golden is an independent reconstruction of the
repository's exact PR76 profile at the literature-defined mixture state. It is not
a claim that the paper printed these phase compositions or used every floating
point convention of MPMC_HNU.

## C++ regression contract

`pr76_sour_gas_three_phase_test.cpp` verifies five separate surfaces.

### Exact fixed-three-phase candidate

The independent phase compositions are passed through the production
`Pr76ThreePhaseEvaluator` and `iterate_pt_three_phase`.

The test requires the unchanged production gates:

- max chemical-potential residual `<= 1e-11`;
- generalized-RR residual `<= 2e-13`;
- absolute material balance `<= 1e-12`;
- positive-feed relative material balance `<= 1e-10`;
- every accepted phase fraction `> 1e-10`.

Fresh PR76 phase properties at the reference compositions are also checked before
the nonlinear solve.

### `2 -> 3` max3 orchestration

The full `solve_pr76_pt_max3(...)` path must first retain evidence that the
two-phase set is incomplete/unstable, then run a fresh three-phase solve and final
common-tangent review.

The reference phase compositions are supplied as **continuation starts only**.
They cannot create phase-count evidence: the existing max3 contract does not read
them until the base two-phase final review has already rejected the pair.

The accepted three-phase attempt must have final stability status
`no_instability_found`.

### Interior `3 -> 3` reweighting

At fixed `P,T`, three coexisting equilibrium compositions define a three-phase
composition simplex. A second strictly interior overall feed is formed with phase
fractions

```text
0.20, 0.35, 0.45
```

using the same independently solved equilibrium phases. The full max3 path must
again close at three phases and reproduce the same equilibrium compositions.

This is a model-consistency/continuation-style regression derived from the
independent physical-mixture reference, not a second literature measurement.

### Runtime component permutation

The full parameter snapshot is rebuilt in reverse component order. Feed and starts
are permuted by component identity, then max3 must recover the same physical
three-phase result after mapping back to canonical IDs.

### Public backend publication

`Pr76PtFlashBackend` is configured with the same bounded continuation and
stability starts. The generic result must

- be structurally valid;
- publish exactly three accepted phases;
- retain `2 -> 3 accepted_target` evidence with a fresh target solve;
- keep `global_stability_proven=false`.

## Reference-comparison budget

The independent equilibrium is much more accurate than the production stopping
criteria, but composition sensitivity can amplify a small chemical-potential
residual. Before hosted C++ execution, this benchmark declares an absolute
`2e-6` comparison budget for phase composition, phase fraction and Z matching.

That comparison budget does **not** replace or relax any production residual gate.
All existing production equilibrium, mass-balance, TPD and phase-fraction
tolerances remain unchanged.

## Stability boundary

The production final common-tangent review is still a finite multi-start TPD
search. Acceptance means that the configured search found no additional
instability; it does not prove mathematical global stability.

The independent Decimal script establishes the equilibrium golden. It does not
pretend to certify every possible TPD minimum. Literature topology, the independent
equilibrium reconstruction, production final stability, and existing synthetic
robustness tests remain distinct evidence layers.

## What this increment does not claim

- no experimental phase-composition accuracy claim;
- no mathematical global-stability proof;
- no liquid/vapor/L1/L2 morphology contract;
- no change to automatic max3 initialization completeness;
- no new PR76/PR78 alpha branch;
- no EOS parameter fitting;
- no `2 -> 1` change;
- no CPA, SW92, physics or frontend change;
- no claim that the derived six-component P–Z vector was printed verbatim by the
  paper.

The next numerical increment should use this benchmark to audit a real
phase-boundary path around the literature-defined three-phase state, particularly
physical-mixture `3 -> 2` disappearance/fresh-neighbor behavior. Production changes
should be made only if that path exposes a concrete defect.
