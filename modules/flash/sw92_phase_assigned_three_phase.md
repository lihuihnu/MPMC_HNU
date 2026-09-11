# SW92 Profile-C C2b.1 unordered three-phase candidate

## Scope

`sw92_phase_assigned_three_phase.hpp` implements the next numerical gate of
`SW92-equilibrium/phase-assigned-aq-na-joint/v1`:

```text
W(AQ) + H0(NA) + H1(NA)
```

The two `H` slots are **symmetric/unordered hydrocarbon phase instances**. They
are not liquid/vapor labels. C2b.1 is a candidate primitive only: it does not
perform autonomous phase-number selection, water appearance/disappearance
orchestration, final topology stability, hydrocarbon morphology resolution, or
publish an authoritative phase set.

Algorithm identity:

```text
SW92-equilibrium/phase-assigned-aq-na-joint/
unordered-w-h0-h1-three-phase-logK-SSI-generalized-RR/v1
```

`global_stability_proven=false`, `final_stability_checked=false`, and
`accepted_phase_set_published=false` are part of the result contract.

## Fixed thermodynamic mapping

For every property evaluation:

- `W -> SwPhaseFamily::aqueous`;
- `H0 -> SwPhaseFamily::nonaqueous`;
- `H1 -> SwPhaseFamily::nonaqueous`;
- all phases use the same ordered SW92 parameter snapshot and fixed NaCl
  molality supplied by the current Profile-C model;
- each family evaluation uses the existing same-family mechanically admissible
  minimum-Gibbs cubic-root selection.

C2b.1 does not alter SW92 formulas, BIPs, root tolerances, generic TPD, C1, or
C2a1.

## Generalized three-phase material balance

For active component `i`, use

```text
K0_i = x_i^H0 / x_i^W
K1_i = x_i^H1 / x_i^W
beta_W = 1 - beta_0 - beta_1
D_i = beta_W + beta_0 K0_i + beta_1 K1_i
x_i^W  = z_i / D_i
x_i^H0 = K0_i z_i / D_i
x_i^H1 = K1_i z_i / D_i
```

The two generalized Rachford-Rice equations are

```text
F0 = sum_i z_i (K0_i - 1) / D_i = 0
F1 = sum_i z_i (K1_i - 1) / D_i = 0
```

The implementation evaluates these without forming potentially overflowing
`K`: for each component it shifts the three exponentials by
`max(0, logK0_i, logK1_i)`. The 2x2 Newton Jacobian is evaluated from the raw
phase compositions,

```text
J00 = -sum_i (x_i^H0 - x_i^W)^2 / z_i
J01 = -sum_i (x_i^H0 - x_i^W)(x_i^H1 - x_i^W) / z_i
J11 = -sum_i (x_i^H1 - x_i^W)^2 / z_i
```

with a damped solve constrained to the **closed** phase-fraction simplex. The
closed boundary is required during phase addition: a valid C1 + negative-TPD
witness naturally starts with the new phase at an infinitesimal/zero-fraction
boundary. Final C2b.1 acceptance still requires all three fractions to be
strictly above `minimum_phase_fraction`.

After the generalized-RR solve, phase compositions are normalized only for
roundoff and the returned fractions/compositions are rechecked against the
original common component inventory. No clipping or composition floor is used.

## Chemical-potential iteration

Two independent log-ratio vectors are advanced:

```text
logK0_i = log(x_i^H0 / x_i^W)
logK1_i = log(x_i^H1 / x_i^W)

r0_i = [log x_i^W  + log phi_i^AQ]
     - [log x_i^H0 + log phi_i^NA]
r1_i = [log x_i^W  + log phi_i^AQ]
     - [log x_i^H1 + log phi_i^NA]
```

The outer SSI step updates both vectors together and uses a residual/Gibbs
line-search guard. A converged candidate must satisfy both chemical-potential
residual sets, generalized-RR residual, returned material balance, positive
phase fractions, pairwise composition distinction, and the Profile-C relative
water-richness topology (`W` water-richer than both `H` phases).

## C1/C2a1 source adapter

`solve_sw92_phase_assigned_c2b1_candidate(...)` accepts:

1. an admissible C1 `W(AQ)+H(NA)` result;
2. a matching C2a1 result from the same model snapshot;
3. one selected C2a1 robust negative NA trial.

Before solving, it independently revalidates the retained C1 equilibrium and
material balance, checks the C2a1 provenance/metadata, and reclassifies the
selected negative trial under the retained C2a1 composition-separation and
water-role rules. The retained H and additional witness become H0/H1 seeds.
The initial retained-H fraction is split only to place the new branch into the
numerical solve; this split is not a thermodynamic prior.

## Unordered H slots and canonical representation

The nonlinear equations are symmetric in H0/H1. After convergence only, the
result is canonicalized so equivalent slot swaps do not become different
candidate solutions. Relative selected `Z` is used only for canonical
representation, as permitted by the C2a2 audit; a component-ID-ordered
lexicographic composition fallback is used near a `Z` representation tie.
Neither rule is a physical `LV` versus `LL` classifier.

A later morphology gate must distinguish at least `LV`, `LL`, and unresolved.
C2b.1 deliberately publishes none of those labels.

## Phase-disappearance semantics

Equation convergence at a phase-fraction boundary is retained as evidence, not
silently converted to a different topology:

- W below the numerical threshold -> `aqueous_phase_disappearance`; later
  water-topology orchestration remains required;
- either H below the threshold -> `hydrocarbon_phase_disappearance`; C2b.1 does
  not accept a lower phase count.

Root/property failures, generalized-RR degeneracy/resource failures, outer
iteration/evaluation limits, and line-search failures remain explicit result
states.

## Independent structural reference

`tests/flash/sw92_phase_assigned_three_phase/reference_decimal.py` is a
stdlib-only Decimal(80) oracle and never imports production code. It reuses the
already-independent SW92 ternary equation transcription and changes only the
CH4/CO2 non-water pair to explicit `synthetic_test kij=0` in both AQ and NA.
The state is `p=3 MPa`, `T=260 K`, fresh water.

This synthetic non-water BIP is necessary because the currently audited SW92
source does not provide a general CH4/CO2 non-water pair. The fixture therefore
makes **no physical/experimental three-phase claim**. It independently solves:

- one W(AQ)+H0(NA)+H1(NA) common-chemical-potential state;
- strictly positive three-phase fractions for a common feed;
- a metastable C1 W+H tie-line for that feed;
- a robustly negative NA TPD toward the second H branch.

Focused C++ tests then cover the complete C1 -> C2a1 -> C2b.1 path, direct H
slot swap/canonicalization, component permutation, W/H disappearance,
root/property and resource failures, source guards, and public-header
self-containment. The dedicated hosted workflow also reruns the affected C1 and
C2a1 regressions and their independent Decimal references.
