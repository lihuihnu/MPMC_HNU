# SW92 Profile-C Gate C2a1 H-side NA phase-addition witness

## Scope

Algorithm identity:

```text
SW92-equilibrium/phase-assigned-aq-na-joint/h-side-na-phase-addition-witness/v1
```

This is the first production increment authorized by the literature-wide C2a topology audit. It
accepts an already admissible Gate-C1 `W(AQ)+H(NA)` candidate and asks only whether the retained
hydrocarbon-side phase `H` has finite-search evidence for adding a second **nonaqueous** phase.

It is intentionally not a final stability test, phase-number decision, hydrocarbon liquid/vapor
classifier, three-phase flash, or accepted phase-set publisher.

## Why the search is NA-only

The post-C1 C2 audit proved that a whole-simplex AQ search at the retained H composition evaluates

```text
D_AQ(y_H) = g_AQ(y_H) - g_NA(y_H),
```

so a negative value merely reintroduces Xu/Profile-B lower-envelope family reassignment. Profile C
instead fixes physical water phase `W -> AQ` and hydrocarbon phase `H -> NA`.

The C2a topology audit identified the next admissible topology question as

```text
W + H  ->  W + H0 + H1,
```

with both H-side phases evaluated by SW92 NA. A negative NA TPD can therefore be retained as an
**additional nonaqueous phase witness**. It is still only a seed for a future joint `W+H0+H1`
calculation.

## Source-candidate revalidation

The adapter does not rely only on the C1 status enum. Before running any stability search it checks
that:

- the C1 result and supplied `Sw92Phase<double>` have the same dataset, revision and ordered
  components;
- the result identifies the Profile-C C1 algorithm and corrected SW92 thermodynamic profile;
- the retained physical mapping is `W/aqueous -> AQ` and `H/nonaqueous -> NA`;
- both phase compositions and activities are finite, normalized, smooth and on exactly the feed
  support;
- phase fractions are positive, resolved and sum to one;
- material balance and cross-phase reduced chemical potentials still satisfy the C1 configured
  tolerances when recomputed from the retained phase state;
- the two phase compositions remain distinct;
- retained W is robustly more water-rich than retained H.

A failed C1 result returns `source_candidate_unavailable`. A nominal C1 candidate whose retained
state no longer passes those checks returns `source_candidate_inconsistent`. Supplying a different
model snapshot is a caller contract error.

## Reconstructed common tangent

For every active component the adapter reconstructs

```text
m_i^W = ln(x_i^W) + ln(phi_i^AQ)
m_i^H = ln(x_i^H) + ln(phi_i^NA)
d_i   = midpoint(m_i^W, m_i^H).
```

The C1 chemical-potential mismatch is retained as

```text
common_reference_allowance = 0.5 * max_i |m_i^W - m_i^H|.
```

The NA finite-search TPD threshold is

```text
effective_tpd_tolerance = base_tpd_tolerance + common_reference_allowance.
```

The allowance is numerical metadata for the imperfect C1 tangent; it is not a physical tolerance,
rigorous error bound, or global-stability allowance. The generic TPD roundoff guard remains active
separately.

The retained H phase is also evaluated algebraically against this reconstructed tangent using its
already retained NA activity. Its TPD must be zero within the common-reference allowance and the
generic arithmetic guard. This check costs no additional property evaluation.

## Diagnostic starts and finite search

The generic `test_pt_stability_against` driver is reused without modification. Its property provider
is one fixed

```text
Sw92FamilyStabilityEvaluator(..., SwPhaseFamily::nonaqueous, ...)
```

for the entire search.

In addition to caller starts and optional generic automatic starts, C2a1 always adds:

1. retained W composition evaluated under NA — diagnostic only;
2. retained H composition under NA — the trivial retained H start.

No AQ property provider exists in the adapter.

The normal generic guarantees remain unchanged:

- no composition floor;
- zero-feed components cannot be introduced;
- same-family mechanically admissible minimum-Gibbs root selection;
- root-envelope nonsmoothness/property failure/resource exhaustion remain explicit;
- finite multistart search never proves global stability.

## Witness topology gates

Every robustly negative NA trial is retained in `negative_witnesses`. A negative point is a usable
H-split seed only if both conditions hold.

### Composition distinction

On active support:

```text
max_i |ln(w_i) - ln(x_i^H)| > log_composition_separation.
```

This prevents the retained H phase, or a numerically equivalent representation of it, from being
promoted as an added phase.

### Water-role guard

Let `w_water` be the negative trial water fraction. The trial must satisfy

```text
x_water^W - w_water > 256*eps*(1 + |x_water^W| + |w_water|).
```

Thus a mathematical NA negative point that is W-like or at the W topology boundary remains
visible diagnostic evidence but cannot become a physically usable H-split seed.

This is only the relative two-role topology guard authorized by C1/C2a. It is not an absolute water
cutoff and is not presented as the empirical validity boundary of the SW92 parameterization.

## Status semantics

`Sw92PhaseAssignedHSideWitnessStatus` has five outcomes:

- `source_candidate_unavailable` — no admissible C1 `W+H` candidate exists;
- `source_candidate_inconsistent` — retained C1 evidence fails independent revalidation;
- `additional_nonaqueous_phase_witness_found` — at least one robust negative NA trial is distinct
  from H and water-poorer than W;
- `no_additional_nonaqueous_witness_found` — all required finite NA trials resolved without robust
  negative TPD;
- `indeterminate` — numerical search is unresolved, or negative mathematical NA trials exist but no
  negative point passes the physical H-split seed guards.

The last distinction is deliberate: a role-inadmissible negative trial cannot be converted to either
an accepted phase or a statement of H-side stability.

The result always states:

```text
global_stability_proven = false
accepted_phase_set_published = false
```

and exposes no `accepted_phase_set()` method.

## Validation

The focused regression suite covers:

1. traceable CO2/H2O C1 state: retained H is a zero-TPD NA point and finite NA search reports no
   additional witness;
2. traceable CH4/H2O, 1 molal NaCl C1 state: same check with nonzero molality;
3. the known negative `AQ-at-H` CO2 blocker is reproduced outside the adapter while C2a1 still
   reports no additional NA witness;
4. explicit `synthetic_test` CO2/water with only the NA water BIP changed to `k_ij=-0.1` produces a
   distinct, water-poorer robust negative NA witness near `x_CO2=0.01`, usable only as an H-split
   seed;
5. the same synthetic model produces a robust negative NA-at-W diagnostic that fails the
   water-role guard and leaves the outer result indeterminate;
6. finite evaluation exhaustion remains indeterminate;
7. NA root/property failure remains indeterminate with the generic family-search diagnostics;
8. failed/mismatched C1 source guards;
9. component permutation;
10. public-header self containment.

The standalone Decimal(80) script, independent of production C++, regenerates approximately

```text
x_CO2^W = 0.00593460227779478
y_CO2^H = 0.987169796047773
D_NA(W) = -0.002894806218964
D_NA(x_CO2=0.01) = -0.003774779298654
```

for that explicit synthetic fixture. The first negative point is intentionally role-inadmissible; the
second is compositionally distinct from retained H and water-poorer than retained W. These values
are structural test evidence only and are not SW92 physical parameter or phase-equilibrium claims.

The initially attempted synthetic `k_ij=0` fixture was independently rejected before acceptance:
its NA TPD did not become negative. No production tolerance or thermodynamic formula was changed
to force the intended test topology.

The existing C1 and post-C1 C2 blocker Decimal scripts are rerun unchanged in the C2a1 workflow.

## Explicitly not implemented

Gate C2a1 does not implement:

- AQ whole-simplex stability under Profile C;
- an `x_H2O=0.5` or any other absolute phase classifier;
- physical `H -> L/V` role assignment;
- a joint `W+H0+H1` / `W+L+V` solver;
- water appearance/disappearance orchestration from other topologies;
- autonomous phase-number selection;
- an authoritative Profile-C phase set;
- salt inventory conservation;
- SW92 flash derivatives or physics coupling.

A future negative C2a1 witness must be consumed only by a separately audited three-phase candidate
solver; it is not itself a new phase.
