# SW92 Xu-style max2 validation expansion

## Scope

This validation increment expands the evidence surface of

```text
SW92-equilibrium/xu-asymmetric-gibbs/max2-logK-SSI-RR/common-tangent/v1
```

without changing production algorithms, thermodynamic formulas, tolerances, or existing reference
anchors.

The requested validation order is covered as follows:

1. independent CH4/H2O max2 regression;
2. nonzero NaCl-molality/brine max2 regression;
3. one runtime three-component max2 snapshot;
4. phase-disappearance, family-envelope switch/tie, and root-envelope nonsmooth boundaries.

All successful stability decisions remain finite-search statements with
`global_stability_proven=false`.

## Primary SW92 source

The model constants and correlations remain those already audited against the project-supplied
Soreide-Whitson paper/errata:

- I. Soreide and C. H. Whitson, *Peng-Robinson predictions for hydrocarbons, CO2, N2, and H2S with
  pure water and NaCl brine*, Fluid Phase Equilibria 77 (1992) 217-240,
  DOI `10.1016/0378-3812(92)85105-H`;
- project PDF SHA-256 `cb5b1d5034d78d934e887449ce0c89692d43431835d1371da5606d95b2c6bf58`.

Relevant source scope is explicit:

- Table 3 supplies CH4/H2O critical properties and acentric factors;
- corrected Eq. (12) supplies the aqueous CH4/water interaction correlation, including salinity;
- Table 5 supplies the non-aqueous CH4/water value `k_ij^NA = 0.4850`;
- the paper states Eq. (12) was developed using hydrocarbon/NaCl-brine data from 0 to 5 molal;
- Figs. 13-14 report methane/NaCl-brine comparisons at 103 C and include 1 and 4 molal cases.

The regression points below are model-calculation states inside that documented temperature/salinity
context. Their overall feed compositions are numerical flash inputs, not experimental measurements.

## 1. Independent CH4/H2O fresh-water max2 regression

State:

```text
p = 10 MPa
T = 350 K
m_NaCl = 0 mol/kg H2O
z_CH4 = 0.5
```

`tests/flash/sw92_asymmetric_max2_validation/reference_decimal.py` imports no production code. It
reuses only the repository's already-independent stdlib Decimal SW92 equation implementation and
adds an independent minimum-Gibbs-root/fixed-family-pair/Xu-envelope calculation.

For the Xu lower-envelope model the independent Decimal(80) result selects an AQ+AQ pair:

```text
x_CH4(low)  = 0.0012488516673043369772690695862810869...
x_CH4(high) = 0.98976716870138545359588622097568157...
beta_high   = 0.50454416447146141705145186641410852...
G_feed/RT   = -1.7966801075584243682333600397004459...
G_pair/RT   = -2.7666079361625361393408003803537260...
G_pair-G_feed = -0.96992782860411177110744034065328018...
```

The final AQ/NA finite common-tangent searches must both complete with no robust negative TPD before
`accepted_phase_set()` publishes the pair.

## 2. Nonzero-molality CH4/brine regression

State:

```text
p = 20 MPa
T = 376.15 K (103 C)
m_NaCl = 4 mol/kg H2O
z_CH4 = 0.5
```

This deliberately exercises the salinity terms in both SW92 water alpha and the corrected aqueous
CH4/water Eq. (12) correlation. The fixed molality must propagate unchanged through Gate 3A,
3B.1/3B.2, final stability, and max2 result metadata.

Independent Decimal(80) anchors are:

```text
x_CH4(low)  = 0.00057927917850318821367216301214118614...
x_CH4(high) = 0.98563843743740236519943874461280656...
beta_high   = 0.50699566278255549394787817509225852...
G_feed/RT   = -1.7674270686613408746172242606525397...
G_pair/RT   = -2.6750517768359819185769968256602970...
G_pair-G_feed = -0.90762470817464104395977256500775723...
```

This remains the existing **fixed-molality conditional model**. It does not add salt as an EOS
component, conserve total salt inventory, or redistribute salt between phases.

## 3. Important profile-separation finding: Xu lower envelope is not original SW92 phase assignment

The expanded independent calculation intentionally also solves the fixed physical-family assignment

```text
phase0 = AQ
phase1 = NA
```

for both methane states.

Fresh-water reference:

```text
AQ+NA x_CH4 = 0.0012548893978171366881230600626601386...
AQ+NA y_CH4 = 0.99448254517749047477840934275427544...
```

4-molal reference:

```text
AQ+NA x_CH4 = 0.00058324135361841879751362867136907871...
AQ+NA y_CH4 = 0.99242009503855451670813326717549818...
```

Those equations converge, but the NA-assigned phase is above the AQ Gibbs surface at the same
composition under the Xu common gauge, so Gate 3B.1 correctly marks the fixed AQ+NA candidate as
`family_assignment_dominated`.

This is **not** a numerical failure. It records the algorithm-profile distinction already stated in
the Xu audit: the Xu profile is a new generalized lower-envelope model over the two SW92 family
surfaces; it is not the original 1992 phase-assigned engineering flash.

The validation therefore must not compare the Xu max2 AQ+AQ result against an AQ+NA engineering
calculation as though they were the same thermodynamic algorithm.

## 4. Family-envelope boundary audit

For both traceable CH4 states, the independent Decimal reference evaluates the family Gibbs gap at
99 interior gas fractions (`0.01` through `0.99`). Every sampled value has

```text
g_AQ(x) - g_NA(x) < 0.
```

The C++ validation independently repeats a smaller interior sample and verifies the pure-component
vertices meet (`Delta g -> 0`).

**Observed validation fact:** no interior AQ<->NA lower-envelope switch is observed in these sampled
CH4 binary scopes. The Xu profile therefore selects AQ+AQ for these accepted pairs.

This finite audit is not a theorem for every component/state. More importantly, `SwPhaseFamily`
remains a **thermodynamic model-family identity**, not a universal physical water-rich/non-aqueous
phase label.

A controlled synthetic BIP fixture is used only to verify the switch/tie machinery itself:

- one side resolves AQ lower;
- the other resolves NA lower;
- exact equality remains a guarded family-envelope tie and max2 publication is withheld.

The synthetic fixture is explicitly structural and carries no physical-validation claim.

## 5. Runtime three-component validation

The runtime chain is exercised with

```text
CH4 / N2 / H2O
p = 10 MPa
T = 350 K
m_NaCl = 1 mol/kg H2O
z = (0.49, 0.01, 0.50)
```

and a full component-order permutation.

The original SW92 source used here does not provide, in the current repository data contract, a
complete sourced CH4/N2 **non-water/non-water** pair suitable for claiming a physical ternary
reference. The test therefore uses a deliberately synthetic `k_CH4,N2 = 0` in both families under
`DataPolicy::allow_synthetic_tests`.

Consequently this case validates only:

- runtime component count 3;
- all-positive active support;
- family-aware max2 orchestration/final review at runtime dimension;
- material balance/common chemical potentials through the existing production gates;
- component permutation invariance and ordered metadata.

It is **not** a physical ternary SW92 validation. A physical ternary regression remains blocked on a
traceable non-water pair dataset compatible with the chosen SW92 parameter contract.

## 6. Phase-disappearance boundary

The fresh-water independent AQ+AQ coexistence anchor is reused without inventing new thermodynamic
numbers. A feed is constructed by material balance at a known small high-gas phase fraction, while
the test deliberately configures a larger numerical `minimum_phase_fraction`.

Expected behavior:

1. Gate 3A still identifies instability inside the coexistence region;
2. a fixed-pair attempt can satisfy chemical-potential and material-balance tolerances;
3. the small phase is then classified as `phase_disappearance`;
4. Gate 3B.2 has no admissible pair to publish;
5. max2 returns `indeterminate` and `accepted_phase_set()==nullptr`.

This checks numerical disappearance semantics. The configured threshold is not an experimental
physical phase-vanishing criterion.

## 7. Root-envelope nonsmooth boundary

The existing independently audited pure-CO2 PR saturation anchor at 280 K is propagated through the
full max2 entry. At the same-family minimum-Gibbs root tie, Gate 3A must remain
`family_reference_nonsmooth`; no pair plan or final search may start, and max2 publication must be
withheld.

This verifies that closing Gate 3B.3 did not bypass the earlier root-envelope nonsmooth contract.

## 8. Acceptance criteria

This validation expansion passes only if all of the following remain true:

- both independent CH4 binary Decimal(80) anchor sets match C++ max2 values;
- fresh and 4-molal cases reach finite-search two-phase acceptance without tolerance relaxation;
- fixed AQ+NA profile-separation solutions remain equation-converged but lower-envelope dominated;
- runtime ternary structure survives a full component permutation without a physical-data claim;
- disappearance, family tie/switch, and root-envelope nonsmooth cases never leak through
  `accepted_phase_set()`;
- the existing Gate 3B.3 regression remains green;
- no production thermodynamics/flash algorithm or existing tolerance/reference anchor is modified.
