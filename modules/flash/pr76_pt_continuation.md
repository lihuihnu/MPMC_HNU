# PR76 ordered PT continuation / phase-transition scan

## Scope

`mpmc/flash/pr76_pt_continuation.hpp` provides a deliberately small ordered-path layer over the established PR76 max-three-phase point solver. It does not introduce a new thermodynamic solve, predictor, pseudo-arclength method, Newton continuation system, morphology classifier, or global-stability proof.

The continuation convention is:

```text
PR76/PT/max3/ordered-path-continuation/v1
```

The input is one ordered component feed and an ordered list of `(pressure_pa, temperature_k)` states. Every point is solved independently through `solve_pr76_pt_max3(...)` at the current state. Previous-point information is only an initialization hint.

## Fresh-point rule

For every path point:

1. validate the current `p,T` state and feed;
2. build current initial/final stability starts from configured fallback starts plus, when available, the previous accepted phase compositions;
3. when the previous accepted point is three-phase, prepend its owned three-phase state as one `Pr76PtThreePhaseStart`;
4. run the full current-state `solve_pr76_pt_max3(...)`;
5. accept only the phase count owned by that current solve;
6. if the current point is unresolved, discard continuation state before the next point.

A previous phase count is never copied forward as a result. A previous three-phase state cannot create current 2→3 evidence: the underlying max3 contract still consumes three-phase starts only after the current point's two-phase final review has already established `phase_set_unstable`.

## Transition brackets

When two immediately adjacent scan points are both accepted and have different phase counts, the scan records a `Pr76PtTransitionBracket`.

The bracket means only:

```text
accepted left endpoint + accepted right endpoint + different phase count
```

It always publishes:

```text
exact_boundary_resolved = false
```

A phase-count difference of one (`1↔2` or `2↔3`) is marked as an adjacent topology step. A sampled `1↔3` jump is explicitly a path-resolution failure and requires refinement before interpreting the topology sequence.

The current synthetic structural regression intentionally exercises a path whose accepted endpoints are:

```text
(1.90 MPa, 325 K) -> 3 phases
(2.80 MPa, 331 K) -> 2 phases
(1.00 MPa, 347 K) -> 2 phases
(1.00 MPa, 348 K) -> 1 phase
```

and the exact reverse path. The resulting discrete sequences are therefore:

```text
3 -> 2 -> 2 -> 1
1 -> 2 -> 2 -> 3
```

The first `3↔2` segment spans a numerically unresolved belt observed during audit/refinement. The regression deliberately preserves that fact by treating the accepted endpoints as a bracket rather than pretending that the exact disappearance boundary was located. The fixture remains `SourceKind::synthetic_test`; these PT coordinates are algorithm/topology regression points, not experimental fluid data.

## Traceable literature critical-endpoint anchor

The continuation test suite also owns a separate literature data contract for the real ternary system:

```text
CO2 + 1-pentanol + n-tridecane
```

This contract is intentionally separate from the synthetic continuation path.

### Experimental anchor: Smits (1996)

J.C. Smits, *Phase Behaviour in Certain CO2 + n-Alkane + 1-Alkanol Systems: Experiments and Modelling*, Delft University of Technology, 1996, TU Delft repository UUID:

```text
ab0a04f3-3491-4314-a95c-a943aff2f955
```

Chapter 5, Table 5.1, system 4 reports for `CO2 + n-tridecane + 1-pentanol`:

```text
x_CO2                         = 0.9501
x*_C13 on a CO2-free basis    = 0.8603
(L=VL) UCEP pressure          = 9.22 MPa
(L=VL) UCEP temperature       = 317.43 K
```

For the project order `CO2, 1-pentanol, n-tridecane`, the isopleth loading is reconstructed explicitly as:

```text
z_CO2        = 0.9501
z_C13        = (1 - 0.9501) * 0.8603 = 0.04292897
z_1-pentanol = (1 - 0.9501) - z_C13 = 0.00697103
```

This is a real experimental **critical endpoint / K-point**. It is not a normal interior three-phase flash state: by definition, two coexisting phases become critical/coincident in the presence of another equilibrium phase. Consequently the repository marks this anchor as:

```text
critical_endpoint = true
suitable_as_distinct_three_phase_flash_oracle = false
```

and does not force it through the fixed-three-distinct-phase acceptance gate.

### PR76 pure data: Mushrif (2004), Appendix B Table B1

Samir H. Mushrif, *Determining Equation of State Binary Interaction Parameters Using K- and L-Points*, University of Saskatchewan, 2004, handle:

```text
hdl:10388/etd-10212004-233350
```

Appendix B, Table B1 (citing Yaws, 1999) gives:

| component | Tc (K) | Pc (bar) | omega |
| --- | ---: | ---: | ---: |
| CO2 | 304.1 | 73.8 | 0.2390 |
| 1-pentanol | 588.2 | 39.1 | 0.5784 |
| n-tridecane | 676.0 | 17.2 | 0.6203 |

The test contract stores pressure in SI Pa with the explicit conversion `bar * 1e5 -> Pa`.

### PR76 binary parameters: Mushrif (2004), Figure 3.4

Section 3.5.1, Figure 3.4 draws a Peng-Robinson Gibbs-energy surface for the same ternary system at `317.43 K` and `9.22 MPa` using:

```text
delta(CO2, 1-pentanol)        = 0.17
delta(CO2, n-tridecane)       = 0.15
delta(1-pentanol, n-tridecane)= 0.10
```

The thesis describes these nonzero values as a preliminary Gibbs-surface choice showing a possible phase split. They are therefore retained as an **illustrative PR76 model anchor**, not relabeled as a quantitatively fitted reproduction of the Smits experimental K-point.

The repository builds these values through the ordinary PR76 contract (`SourceKind::literature`, `DataPolicy::ordinary`) and verifies component order, units, conversions, every required binary pair, and provenance metadata.

## Why the K-point is not a three-phase flash oracle

The project's generic fixed-three-phase primitive requires three pairwise-distinct phase compositions before it can publish a candidate. That is appropriate for an interior three-phase flash state, but a K-point is a critical endpoint where two phases coalesce. Treating the Smits/Mushrif state as a three-distinct-phase numerical oracle would therefore contradict both the experimental meaning of the state and the solver's explicit separation gate.

A future quantitative literature flash regression should only be added when one source chain provides, for an interior three-phase state compatible with the current plain PR76 mixing contract:

- pressure and temperature;
- overall feed or enough information to reconstruct it;
- all required pure-component parameters;
- all constant symmetric PR76 `kij` values;
- independent phase fractions/compositions (and preferably densities or Z values);
- unambiguous units and source locators.

Missing values must not be inferred from unrelated correlations merely to manufacture a reference case.

## Non-capabilities

This layer does not provide:

- exact transition-boundary location;
- pseudo-arclength or derivative-based continuation;
- a guarantee that every finite start set resolves every PT point;
- liquid/vapor/L1/L2 morphology classification;
- mathematical global-stability proof;
- an experimental validation of the synthetic bidirectional phase sequence;
- a three-distinct-phase experimental flash oracle at the Smits K-point;
- frontend or physics coupling.
