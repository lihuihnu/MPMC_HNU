# SW92 Xu-style Gate 3B.3 maximum-two-phase final acceptance

## Scope

This increment closes the staged Gate 3B design for

```text
SW92-equilibrium/xu-asymmetric-gibbs/v1
```

with the maximum-two-phase algorithm identity

```text
SW92-equilibrium/xu-asymmetric-gibbs/max2-logK-SSI-RR/common-tangent/v1
```

It composes already reviewed capabilities rather than replacing them:

```text
Gate 3A
  lower-envelope feed reference + common-tangent AQ/NA finite stability

Gate 3B.1
  one explicit family-pair material-balance/common-chemical-potential solve
  + per-phase lower-envelope family assignment

Gate 3B.2
  family-neutral witness planning + same-family alternatives
  + complete-plan fixed-pair attempts + candidate deduplication/Gibbs selection

Gate 3B.3
  pair-vs-feed Gibbs gate + selected-pair common tangent
  + final AQ/NA stability + authoritative family-aware publication guard
```

The thermodynamic profile remains

```text
SW92/corrected-original/PR76-base/NaCl-molality
```

and NaCl molality remains an externally fixed coordinate. Gate 3B conserves only the ordered EOS
component inventory; it does not conserve or redistribute salt.

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

accepts an owned/copyable Gate-3B.2 `Sw92AsymmetricPairSelectionResult` and performs only final
acceptance. It validates that the selection and supplied `Sw92Phase<double>` refer to the same
ordered dataset/revision/component snapshot and algorithm profile before publishing anything.

The separate finalizer exists so final-search budgets can be tested/reviewed independently without
rerunning or reinterpreting Gate 3B.2 candidate selection. It is not permission to fabricate or
manually splice a selection result.

## Family-aware phase-set contract

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

The type deliberately says **candidate** because rejected and indeterminate final-review states retain
the same phase data for diagnostics. Publication is controlled only by

```cpp
result.accepted_phase_set()
```

which returns non-null only for the two accepted top-level statuses and performs a structural guard on
phase count, family identities, fractions, composition normalization and active feed support.

Vector position, cubic-root branch and Z remain diagnostics. They are not universal liquid/vapor or
AQ/NA selectors.

## Top-level statuses

`Sw92AsymmetricMax2Status` has five outcomes:

- `single_phase_no_instability_found`: Gate 3A resolved a lower feed family and both finite family
  searches found no instability. The family-aware feed phase is publishable without forcing a split.
- `two_phase_no_instability_found`: a unique Gate-3B.2 pair passed pair-vs-feed Gibbs and both final
  finite AQ/NA searches against the selected-pair common tangent.
- `phase_set_unstable`: the pair remained equation/Gibbs admissible, but a final family search found a
  robust negative TPD witness. The candidate is retained but not published.
- `pair_gibbs_above_feed`: the selected pair's reduced Gibbs is resolved above the feasible original
  lower-envelope feed state. The candidate is retained but no final search can promote it.
- `indeterminate`: initial selection or a required final family search did not support a reliable
  acceptance decision.

Every path keeps

```text
global_stability_proven = false
```

because all TPD searches remain finite multistart searches.

## Single-phase acceptance

Gate 3B.2 can retain a single-phase candidate only after Gate 3A returns
`no_instability_found` with a resolved lower feed family. That Gate-3A result already consists of two
required AQ/NA finite searches against the lower-envelope feed tangent.

Gate 3B.3 therefore publishes the one-phase candidate directly and does not perform a redundant second
final search. This is still qualified finite-search acceptance, not global stability certification.

## Pair-vs-lower-feed Gibbs gate

A selected Gate-3B.2 pair is not yet accepted. Gate 3B.3 first compares it with the feasible original
feed state.

For the selected lower feed family `Fz`:

```text
G_feed/(RT)
  = sum_i z_i [ln(z_i) + ln(phi_i^Fz(z))]
```

and Gate 3B.1 already stores

```text
G_pair/(RT)
  = sum_alpha beta_alpha
      sum_i x_i^alpha [ln(x_i^alpha) + ln(phi_i^Falpha(x^alpha))].
```

The component reference terms omitted from both expressions are identical linear functions of the
same total component inventory, so they cancel in this comparison under the audited Gate-3A common
gauge.

Both feed and pair carry arithmetic roundoff guards. If

```text
G_pair - G_feed > guard_pair + guard_feed
```

the pair is `pair_gibbs_above_feed` and cannot be accepted. A difference unresolved within the
combined arithmetic guard is not rejected by this test alone; final stability still decides the
remaining acceptance path.

The guards are floating-point arithmetic heuristics, not rigorous interval enclosures and not
physical tolerances.

## Selected-pair common tangent

Gate 3B.1 stores for every active component

```text
m_i^0 = ln(x_i^0) + ln(phi_i^F0)
m_i^1 = ln(x_i^1) + ln(phi_i^F1)

d_i = midpoint(m_i^0, m_i^1).
```

Gate 3B.3 uses that exact stored midpoint vector as the imposed reference for both final family
searches. It does not rebuild an AQ-only or NA-only tangent and does not pick either phase by slot or
family order.

Because the pair equations are converged only to a finite residual, define

```text
common_reference_allowance
  = 0.5 * max_i |m_i^0 - m_i^1|.
```

The user-supplied final AQ and NA base TPD tolerances are retained separately. The effective search
tolerances are

```text
effective_tpd_tolerance_F
  = base_tpd_tolerance_F + common_reference_allowance.
```

This allowance represents only the numerical disagreement of the two converged phase references. It
is not an empirical adjustment, model uncertainty, or rigorous error bound.

## Mandatory final starts

Both selected phase compositions are appended as explicit starts to **both** required final family
searches:

```text
AQ final starts: caller AQ starts + phase0 composition + phase1 composition
NA final starts: caller NA starts + phase0 composition + phase1 composition
```

This directly checks the discrete possibility that the opposite family is lower at an already
converged phase composition and also strengthens the same-family phase-set review.

Start quotas are preflighted by the existing asymmetric common-tangent helper before any final
thermodynamic property call. If the configured family search cannot accommodate the mandatory starts,
the call fails as an invalid resource contract instead of silently dropping required evidence.

## Final two-family stability

Gate 3B.3 calls the already reviewed

```cpp
test_sw92_pt_asymmetric_stability_against(...)
```

once with the selected-pair common tangent and the effective final options.

The helper performs independent AQ and NA finite searches with separate resource accounting and
family-tagged negative witnesses.

Final mapping is strict:

```text
AQ unstable OR NA unstable
  -> phase_set_unstable

no negative witness, but AQ or NA indeterminate
  -> indeterminate

AQ no_instability_found AND NA no_instability_found
  -> two_phase_no_instability_found
```

Equation convergence, lower-envelope family assignment, candidate Gibbs ranking or pair-vs-feed Gibbs
never override a negative or indeterminate required final search.

## Authoritative publication boundary

For the first time in the Xu-style path, Gate 3B.3 provides an authoritative family-aware
`accepted_phase_set()`.

A two-phase set is exposed only after the retained Gate-3B.2 evidence proves that:

- the complete configured family-neutral attempt plan was evaluated;
- the selected fixed pair is Gate-3B.1 `candidate_admissible()`;
- equivalent/slot-swapped candidates were deduplicated before Gibbs ranking;
- no distinct candidate Gibbs tie remained;
- pair Gibbs is not resolved above the lower-envelope feed Gibbs;
- both final family searches used the exact same selected-pair common tangent;
- both selected phase compositions were explicit starts in both searches;
- both final searches completed with `no_instability_found`.

The result remains a **maximum-two-phase candidate under the declared finite search**, not a global
phase-count or stability proof. A three-phase state can still exist outside the finite evidence or the
max2 capability.

## Independent reference

`tests/flash/sw92_asymmetric_max2/reference_decimal.py` imports no production code. It reuses the
previously independent Gate-3B.1 Decimal implementation of corrected SW92/PR equations and independently
computes the Gate-3B.3-only quantities for the traceable CO2/H2O state

```text
p = 3 MPa
T = 340 K
NaCl molality = 0
z_CO2 = 0.7
```

including:

- lower-envelope AQ feed reduced Gibbs;
- selected AQ+AQ pair reduced Gibbs and pair-minus-feed Gibbs;
- selected-pair common reduced tangent;
- prescribed AQ/NA TPD at both selected phase compositions;
- prescribed AQ/NA TPD at the uniform composition;
- a deterministic dense binary probe used only as a gross-regression check.

The Decimal path verifies that the two selected AQ phases lie on the common tangent and that changing
the same compositions to the NA family raises the prescribed TPD for this reference state. These are
model numerical anchors, not experimental validation or interval/global stability evidence.

## Verification and affected chain

The Gate-3B.3 workflow must run on GitHub-hosted runners for:

- GCC Debug + ASan/UBSan;
- Clang Release;
- MSVC Release.

It runs the new Gate-3B.3 cases plus affected Gate 3B.2, Gate 3B.1 and Gate 3A regressions, and
regenerates the independent Decimal references for Gate 3A, Gate 3B.1 and Gate 3B.3.

The dedicated Gate-3B.3 tests cover:

- accepted traceable two-phase reference;
- pair-vs-feed Gibbs anchor;
- exact selected-pair common tangent;
- mandatory candidate starts in both final family searches;
- prescribed AQ/NA TPD anchors;
- accepted finite-search single-phase flow;
- Gate-3B.2 indeterminate/attempt-limit propagation;
- final-search indeterminate publication rejection;
- synthetic pair-above-feed publication rejection;
- structural accepted-phase-set guard;
- component permutation;
- final-start quota preflight;
- public-header self containment.

## Explicitly deferred

Gate 3B.3 still does not implement:

- three-phase flash or automatic phase addition after `phase_set_unstable`;
- interval/global stability certification;
- a universal water-rich/hydrocarbon-rich phase label;
- salt-inventory conservation or electrolyte components;
- SW92 flash derivatives/implicit sensitivity;
- SW92 physics closure;
- CPA.

The Whitson dual-model observable profile remains a separate engineering compatibility route. The
Xu-style max2 result must not be substituted back into Profile A semantics, and Profile-A independent
AQ/NA runs must not be spliced into this joint material-balance result.
