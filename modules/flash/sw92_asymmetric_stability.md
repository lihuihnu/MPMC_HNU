# SW92 Xu-style asymmetric stability foundation

## Scope

This increment implements **Gate 3A** from `sw92_xu_asymmetric_audit.md` for the reserved parent
profile

```text
SW92-equilibrium/xu-asymmetric-gibbs/v1
```

The implemented numerical search identity is

```text
SW92-equilibrium/xu-asymmetric-gibbs/stability-foundation/v1
```

It stops at feed-family selection and finite common-reference stability search. It does **not**
solve a joint AQ/NA phase split, publish a joint phase set, conserve NaCl inventory, add SW92
derivatives, or implement three-phase equilibrium.

Public entry point:

```text
test_sw92_pt_asymmetric_stability
```

Public result:

```text
Sw92AsymmetricStabilityResult
```

## Audit refinement: pure-vertex gauge invariant

The implementation-time independent reference review refined one sentence from the pre-code audit.
At a pure-component simplex vertex, AQ and NA have the same pure-component EOS because all
cross-interaction terms vanish. Therefore they have the same root set and the same fugacity
coefficient for the **active pure component**.

The fugacity coefficient reported for a zero-composition *other* component is an infinite-dilution
quantity and can still depend on AQ/NA cross BIPs. The full `ln(phi)` vector is therefore not
required to be identical at a pure vertex.

This refinement does not weaken the common-gauge conclusion. Both families still use the same
ordered component identities, pure-component basis and PR chemical-potential convention; family
BIPs change physical residual chemical potentials, not the component reference-state zero. The
independent regression consequently checks active-component gauge equality, not equality of
irrelevant zero-component infinite-dilution values.

## Feed reference

For one validated `Sw92Phase<double>` snapshot and one fixed NaCl molality, the adapter constructs
one AQ and one NA `Sw92FamilyStabilityEvaluator`. Both evaluate the same normalized feed and each
returns its same-family mechanically admissible minimum-Gibbs cubic root.

Both required feed-family evaluations are attempted independently for numerical property failures.
Their family identity, selected root diagnostics, issues and evaluation counts remain separate.

If either feed family is unavailable or either same-family minimum envelope is nonsmooth, no
cross-family tangent is fabricated and the result is `indeterminate`.

For resolved smooth references, family ordering uses

```text
Delta g/(RT) = sum_i z_i [ln(phi_i^AQ(z)) - ln(phi_i^NA(z))].
```

The ideal mixing term is identical and cancels. The comparison uses compensated summation and an
arithmetic guard

```text
256 * epsilon * [1 + sum_i z_i (|ln(phi_i^AQ)| + |ln(phi_i^NA)|)].
```

This is a floating-point resolution guard, not a physical phase-selection tolerance. If the Gibbs
gap does not exceed the guard, the feed lower envelope is treated as nonsmooth and remains
`indeterminate`; execution order, root index and water fraction never break the tie.

If the lower family `F0` is resolved, the common reduced feed tangent is

```text
d_i = ln(z_i) + ln(phi_i^F0(z))
```

for every active feed component. Inactive feed components retain the finite placeholder required by
the existing generic imposed-reference API and cannot be introduced by trial compositions.

## Two required family searches

The AQ and NA family searches both call the existing generic `test_pt_stability_against` with the
**same** `d_i`. The SW92-specific adapter owns the proof that AQ and NA use a compatible reduced
chemical-potential gauge; the generic TPD algorithm itself is unchanged.

Each family has independent explicit root and stability options and independent extra starts. A
preflight validates both families' search resources and starts before any thermodynamic call, so an
invalid NA search contract cannot be hidden by earlier AQ work or vice versa.

Each family result retains:

- feed-reference evaluation count;
- feed minimum-Gibbs `StabilityPhase` or family-specific property issue;
- feed reduced Gibbs value;
- the full generic finite TPD `StabilityResult`, including its own options, trials, lowest sampled
  point and evaluation count;
- an overflow-safe helper for total feed-reference plus TPD property evaluations.

Every robust negative witness copied to the outer result carries explicit
`SwPhaseFamily::{aqueous,nonaqueous}` and the originating trial index. Family identity is never
encoded in `StabilityPhase::branch`; branch remains only a cubic-root diagnostic inside one family.

Combined status is exactly:

- either family has a robust negative TPD witness: `unstable`;
- otherwise, either required family search is indeterminate: `indeterminate`;
- otherwise both finite searches report no instability: `no_instability_found`.

A negative witness is sufficient instability evidence even if the opposite family search is
indeterminate. No finite absence of witnesses becomes a global proof:

```text
Sw92AsymmetricStabilityResult::global_stability_proven == false
```

## Fixed-molality semantics

NaCl molality remains one external conditional coordinate of
`SW92/corrected-original/PR76-base/NaCl-molality`. The same value is used by both feed references and
all AQ/NA trial properties.

Gate 3A does not add Na+/Cl- EOS components and does not claim:

- conservation of total salt inventory;
- salt redistribution among phases;
- automatic molality updates when phase amounts change;
- a closed H2O/NaCl/gas equilibrium problem.

## Independent validation

The focused regression uses existing sourced SW92 Table 3/5 data and corrected-original
correlations; it does not invent new experimental constants.

A stdlib-only `Decimal(80)` oracle independently rebuilds:

- corrected SW92 water alpha and AQ BIP correlations;
- sourced NA BIPs;
- classical PR mixing;
- the original-Z PR cubic and mechanically admissible roots;
- same-family minimum-Gibbs root ranking;
- feed-family Gibbs difference;
- common reduced feed tangent;
- prescribed AQ and NA TPD values against the same tangent.

The principal CO2/H2O anchor is `p=3 MPa`, `T=340 K`, fresh water and feed
`z_CO2=0.7`. The independent reference selects the AQ feed family and supplies a prescribed uniform
trial `w_CO2=0.5` with robust negative TPD in both family surfaces. The production finite search
therefore has explicit negative-witness evidence rather than relying only on a status enum.

Additional structural regressions cover:

- active-component pure-vertex AQ/NA gauge equality;
- exact and one-ULP-near cross-family Gibbs ties using a synthetic-test NA BIP chosen to coincide
  with the AQ correlation at one prescribed state;
- component permutation;
- AQ-only root-budget failure while the NA feed reference is still evaluated and retained;
- same-family root-envelope nonsmooth propagation using a high-precision PR pure-CO2 saturation
  pressure anchor;
- public-header self containment.

Synthetic data are used only to create controlled tie/routing states and carry explicit
`synthetic_test` provenance. They are not physical validation data.

## Explicit non-goals and next gate

This increment does not modify `pt_stability.hpp`, `pt_split.hpp`, Rachford-Rice, PR76, SW92
thermodynamic formulas, existing tolerances or existing reference anchors.

It also does not reuse `solve_pt_vle` as a joint asymmetric solver. Gate 3B requires a separate
phase-family-aware nonlinear/result contract because low/high-Z candidate roles are not AQ/NA
identity and the current generic `PtCandidatePhase` does not preserve thermodynamic family identity.

The next implementation gate, only after review of this foundation, is a maximum-two-phase joint
AQ/NA split with explicit phase families, one common EOS-component material balance, common reduced
chemical potentials, and final AQ **and** NA common-reference stability review. A converged pair
without that final asymmetric review must remain unaccepted.
