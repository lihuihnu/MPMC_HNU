# SW92 Xu-style Gate 3B.3 maximum-two-phase final acceptance

## Scope

This increment closes the staged maximum-two-phase Gate 3B design under

```text
SW92-equilibrium/xu-asymmetric-gibbs/v1
```

with algorithm identity

```text
SW92-equilibrium/xu-asymmetric-gibbs/max2-logK-SSI-RR/common-tangent/v1
```

The chain is deliberately staged:

```text
Gate 3A
  lower-envelope feed reference + one common-tangent AQ/NA finite stability

Gate 3B.1
  one explicit family-pair joint material balance/common chemical potentials
  + per-phase lower-envelope family assignment

Gate 3B.2
  family-neutral witness plan + mandatory same-family alternatives
  + complete-plan fixed-pair attempts + candidate dedup/Gibbs selection

Gate 3B.3
  revalidated continuous pair evidence + pair-vs-feed Gibbs
  + selected-pair common tangent + final AQ/NA stability
  + authoritative family-aware publication guard
```

The thermodynamic profile remains

```text
SW92/corrected-original/PR76-base/NaCl-molality
```

and NaCl molality remains externally fixed. Only the ordered EOS-component inventory is conserved;
this is not a closed salt/electrolyte equilibrium.

## Public entry points

The full entry is

```cpp
solve_sw92_xu_asymmetric_max2(...)
```

and returns `Sw92AsymmetricMax2Result`.

The lower-level

```cpp
finalize_sw92_asymmetric_max2_selection(...)
```

accepts a Gate-3B.2 `Sw92AsymmetricPairSelectionResult` and performs final acceptance. Before using it,
the finalizer checks dataset/revision/component order and algorithm identities against the supplied
`Sw92Phase<double>` snapshot and retained Gate-3A evidence.

For a pair, it also rebuilds candidate classes and Gibbs choice from the retained attempts without new
thermodynamic calls. A truncated plan, a distinct candidate Gibbs tie, or loss of the selected
Gate-3B.1 admissible pair therefore cannot be hidden behind a stale cached top-level selection field.

For a single phase, the finalizer requires the retained Gate-3A result itself to be
`no_instability_found`, its lower family to be resolved, and the stored feed/reference candidate data
to agree.

## Family-aware candidate and publication contract

Gate 3B.3 introduces

```text
Sw92AsymmetricCandidatePhase
  family
  mole_phase_fraction
  composition
  StabilityPhase
  optional compressibility_factor

Sw92AsymmetricCandidatePhaseSet
  phases
```

The type deliberately says **candidate** because rejected/indeterminate final-review states retain the
same phase data for diagnostics. Only

```cpp
result.accepted_phase_set()
```

is an authoritative publication surface. It is non-null only on accepted one- or two-phase statuses
and additionally checks phase count, family identities, positive phase fractions, fraction sum,
composition normalization and exact active-feed support.

Phase-vector slot, cubic-root index and Z remain diagnostics; none is a universal liquid/vapor or
AQ/NA label.

## Top-level statuses

`Sw92AsymmetricMax2Status` has five outcomes:

- `single_phase_no_instability_found`: Gate 3A resolved the lower feed family and both finite family
  searches found no instability; the family-aware feed phase is publishable without forcing a split.
- `two_phase_no_instability_found`: a complete, unique Gate-3B.2 pair passed recomputed equation/balance
  evidence, pair-vs-feed Gibbs and both final AQ/NA searches.
- `phase_set_unstable`: a final family search found a robust negative TPD witness; the candidate is
  retained but not published.
- `pair_gibbs_above_feed`: the pair reduced Gibbs is resolved above the feasible lower-envelope feed
  state; the candidate is retained and final stability is not used to rescue it.
- `indeterminate`: required selection/equation/final-search evidence is unresolved or inconsistent.

Every path keeps

```text
global_stability_proven = false
```

because stability remains finite multistart search rather than interval/global proof.

## Single-phase acceptance

Gate 3B.2 retains a one-phase candidate only after Gate 3A returns `no_instability_found` with a
resolved lower family. That Gate-3A result already contains the required AQ and NA searches against
the lower-envelope feed tangent.

Gate 3B.3 therefore validates the retained Gate-3A/candidate consistency and publishes the one-phase
state directly. It does not force a split or repeat an equivalent final search. This remains
finite-search acceptance, not a proof of global stability.

## Recomputed pair evidence

A Gate-3B.2 selected pair is not accepted merely because its cached Gate-3B.1 fields say it converged.
Gate 3B.3 recomputes acceptance-critical continuous evidence from the retained phase fractions,
compositions and `ln(phi)` values:

```text
m_i^0 = ln(x_i^0) + ln(phi_i^F0)
m_i^1 = ln(x_i^1) + ln(phi_i^F1)

r_i = m_i^0 - m_i^1
z_i(recovered) = beta_0*x_i^0 + beta_1*x_i^1
```

The recomputed state must still satisfy the configured Gate-3B.1 chemical-potential, absolute/relative
material-balance, minimum-phase-fraction and composition-distinction gates. Same-family phase
activities must remain smooth.

This does not rerun the fixed-pair solver and does not change its tolerances. It prevents stale cached
`chemical_potential_norm`, mass residual, reduced Gibbs or common-tangent fields from becoming the
sole basis of authoritative publication.

The Gate-3B.1 per-phase lower-envelope family assignment remains part of the selected
`candidate_admissible()` evidence. The final AQ/NA searches then provide the independent phase-set
stability check, with both candidate compositions included as starts.

## Pair-vs-lower-feed Gibbs gate

The finalizer recomputes pair reduced Gibbs from the retained phase data:

```text
G_pair/(RT)
  = sum_alpha beta_alpha
      sum_i x_i^alpha [ln(x_i^alpha) + ln(phi_i^Falpha)].
```

It also recomputes the lower-feed reduced Gibbs from the retained Gate-3A lower-family feed reference:

```text
G_feed/(RT)
  = sum_i z_i [ln(z_i) + ln(phi_i^Fz(z))].
```

The omitted component reference terms are identical linear functions of the same total component
inventory and cancel under the already audited common gauge.

Both values receive arithmetic roundoff guards. If

```text
G_pair - G_feed > guard_pair + guard_feed
```

the result is `pair_gibbs_above_feed`; final stability is not allowed to promote that state.

The guards are floating-point arithmetic heuristics, not physical tolerances or interval error bounds.

## Recomputed selected-pair common tangent

From the same retained phase data, Gate 3B.3 recomputes

```text
d_i = midpoint(m_i^0, m_i^1)
```

and uses this one vector as the imposed reference for both final family searches. For an unmodified
Gate-3B.1 result this reproduces its stored midpoint tangent; the finalizer does not need to trust the
cached vector to obtain it.

The numerical reference allowance is

```text
common_reference_allowance
  = 0.5 * max_i |m_i^0 - m_i^1|.
```

AQ and NA base TPD tolerances are retained in result metadata. Their effective values are

```text
effective_tpd_tolerance_F
  = base_tpd_tolerance_F + common_reference_allowance.
```

Only the TPD threshold receives this allowance. Stationarity, root conditioning, iteration limits and
other stability semantics are unchanged. The generic contract permitting `tpd_tolerance == 0` is
preserved.

The allowance is a numerical reference-disagreement allowance, not an empirical model adjustment or a
rigorous enclosure.

## Mandatory final starts

Both selected phase compositions are appended as explicit starts to **both** required family searches:

```text
AQ: caller AQ starts + phase0 + phase1
NA: caller NA starts + phase0 + phase1
```

This directly checks the possibility that the opposite family lowers Gibbs at an already converged
phase composition and strengthens same-family phase-set review without a composition cutoff.

The existing asymmetric common-tangent helper preflights both final search contracts before its final
thermodynamic property calls. If configured quotas cannot accommodate the mandatory starts, the call
fails rather than dropping required evidence.

## Final two-family stability

Gate 3B.3 calls

```cpp
test_sw92_pt_asymmetric_stability_against(...)
```

with the recomputed selected-pair common tangent and effective final options. AQ and NA searches keep
separate options, resource accounting and family-tagged witnesses.

Mapping is strict:

```text
AQ unstable OR NA unstable
  -> phase_set_unstable

no robust negative witness, but AQ or NA indeterminate
  -> indeterminate

AQ no_instability_found AND NA no_instability_found
  -> two_phase_no_instability_found
```

No earlier equation/Gibbs result overrides a negative or unresolved final family search.

## Authoritative publication boundary

For the first time in the Xu-style path, Gate 3B.3 provides a family-aware `accepted_phase_set()`.
For two phases, publication requires the retained evidence chain to support all of the following:

- the configured family-neutral plan was complete rather than attempt-budget truncated;
- the selected fixed pair remains Gate-3B.1 admissible;
- candidate deduplication/Gibbs choice recomputed from attempts yields one unique class;
- recomputed common chemical potentials, material balance, phase fractions and composition distinction
  still satisfy Gate-3B.1 thresholds;
- recomputed pair Gibbs is not resolved above lower-feed Gibbs;
- one recomputed common tangent is supplied to both final family searches;
- both selected phase compositions are starts in both searches;
- both final searches return `no_instability_found`.

This is a **maximum-two-phase result under the declared finite search**. It is not a global phase-count
certificate; a third phase may still exist outside that finite evidence/capability.

## Independent reference

`tests/flash/sw92_asymmetric_max2/reference_decimal.py` imports no production code. It reuses the
previously independent Gate-3B.1 Decimal implementation of corrected SW92/PR equations, then
independently computes Gate-3B.3 quantities for the traceable binary state

```text
CO2/H2O
p = 3 MPa
T = 340 K
NaCl molality = 0
z_CO2 = 0.7
```

It independently anchors lower-feed AQ Gibbs, selected AQ+AQ pair Gibbs, pair-minus-feed Gibbs,
selected-pair common tangent and prescribed AQ/NA TPD values at both selected phase compositions and
the uniform composition.

A bounded 99-point-per-family coarse binary probe catches gross tangent/sign regressions only; it is
explicitly not used as a global stability claim. The Decimal path verifies the selected AQ phases lie
on their common tangent and the NA family is higher at those phase compositions for this reference.
These are model numerical anchors, not experimental validation.

## Verification and failure coverage

The dedicated hosted workflow uses:

- GCC Debug + ASan/UBSan;
- Clang Release;
- MSVC Release.

It runs Gate 3B.3 plus affected Gate 3B.2, Gate 3B.1 and Gate 3A regressions and regenerates independent
Decimal references for Gate 3A, 3B.1 and 3B.3.

Gate-3B.3 focused cases cover:

- accepted traceable two-phase reference and pair-vs-feed Gibbs anchor;
- recomputed common tangent and mandatory candidate starts in both final families;
- prescribed AQ/NA TPD anchors;
- accepted finite-search single-phase flow;
- Gate-3B.2 attempt-limit/indeterminate withholding;
- final-search indeterminate withholding;
- synthetic equation-preserving perturbation that forces a robust final negative witness and
  `phase_set_unstable`;
- synthetic equation-preserving Gibbs shift that forces `pair_gibbs_above_feed`;
- structural accepted-phase-set guard;
- component permutation;
- final-start quota preflight;
- zero-base-TPD-tolerance compatibility;
- public-header self containment.

Synthetic perturbations above test publication routing only and carry no physical/reference claim.

## Explicitly deferred

Gate 3B.3 still does not implement:

- three-phase flash or automatic phase addition after `phase_set_unstable`;
- interval/global stability certification;
- a universal water-rich/hydrocarbon-rich phase label;
- salt-inventory conservation/electrolyte components;
- SW92 flash derivatives or implicit sensitivity;
- SW92 physics closure;
- CPA.

The Whitson dual-model observable profile remains a separate engineering compatibility route. Its
independent AQ/NA runs must not be spliced into this joint material-balance result, and the Xu-style
max2 result must not be relabeled as Profile-A observables.
