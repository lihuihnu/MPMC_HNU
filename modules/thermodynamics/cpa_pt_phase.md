# CPA PT density-root and fugacity kernel

## Scope

This gate raises the CPA baseline from a supplied-density pressure primitive to a fixed-composition PT phase-property kernel:

`p, T, x -> detected density roots -> Z, ln(phi)`

The implementation remains below phase stability and flash. It does not select a final phase count and it does not claim a mathematical proof that a finite density scan found every possible CPA root.

The model profile remains:

`CPA/SRK-physical/simplified-rdf-1.9eta/explicit-site-pairs/v1`

and the PT numerical convention is:

`CPA/SRK-association/density-root-finite-scan/lnphi-helmholtz/v1`.

## Fugacity construction

The implementation uses the residual-Helmholtz definition

```text
ln(phi_i) = mu_i^r/(R T) - ln(Z)
```

instead of mixing a cubic-only compressibility factor with the total CPA compressibility factor.

For the SRK physical residual Helmholtz contribution,

```text
alpha_cub^r = -ln(1-b rho)
              - a/(b R T) ln(1+b rho)
```

and its fixed-`T,V,n_j` chemical-potential derivative is evaluated analytically with the classical one-fluid `a` and `b` rules.

For association, with converged site fractions `X_Ai`,

```text
alpha_assoc^r = sum_i x_i sum_A n_Ai
                [ln(X_Ai) - X_Ai/2 + 1/2]
```

and the configured `g=1/(1-1.9 eta)` gives the analytic chemical-potential contribution

```text
mu_i,assoc^r/(R T)
  = sum_A n_Ai ln(X_Ai)
    - rho b_i/(8 g) (dg/deta)
      sum_j x_j sum_A n_Aj (1-X_Aj)
```

with `(1/g) dg/deta = 1.9 g` for this profile.

The final total `Z` is always the requested PT pressure divided by `rho R T`; it is not the SRK physical-only `Z`.

Focused validation differentiates the extensive residual Helmholtz energy with respect to individual mole numbers at fixed `T,V,n_j`, re-solving the independent analytic one-site association relation after each perturbation. This test-only derivative is compared against production `ln(phi)` and is not used by production code.

## Density-root search

CPA pressure is not a cubic polynomial once association is active, so PR76/SW92 algebraic cubic-root code is deliberately not reused.

The root coordinate is

```text
u = b_mix rho,    0 < u < 1
```

which makes the SRK covolume boundary explicit. The finite search:

1. includes the analytic zero-density residual `P(0)-P_target = -P_target`;
2. samples a bounded Chebyshev-like grid `u_k = u_max sin^2(pi k/(2N))`;
3. brackets every detected sign-changing simple root;
4. refines each bracket by bisection to the declared absolute/relative pressure tolerance;
5. marks negative-to-positive crossings as mechanically stable simple roots and positive-to-negative crossings as mechanically unstable simple roots;
6. retains a scan point already within pressure tolerance with equal-sign neighbors as a tangent/near-multiple root, which makes the overall root set non-smooth rather than silently accepting it.

The scan has explicit interval, evaluation, root-count and bisection budgets. Exhaustion is a result status, not a successful root set.

This is a finite numerical search. A pair of very closely spaced roots can in principle lie between scan points, so `success` means all simple roots detected by the configured finite search were resolved; it is not a global root-topology proof. The later CPA stability adapter must preserve this limitation and must not publish mathematical global stability.

## Independent structural validation

All current CPA PT fixtures are explicitly synthetic.

### Non-associating SRK limit

A one-component non-associating state with three SRK roots is used to verify:

- exactly three roots are recovered by the CPA density search;
- stable/unstable/stable crossing order is preserved;
- each returned `Z` satisfies the independent SRK cubic polynomial;
- each `ln(phi)` agrees with the standard explicit SRK fugacity expression.

### Associating Helmholtz derivative

A binary fixture with one explicit self-associating site class uses an independent analytic quadratic for `X`. At a fixed density, the test independently constructs the total residual Helmholtz energy and forms centered fixed-`T,V,n_j` mole-number derivatives. These derivatives validate the complete production CPA `ln(phi)`, including the effect of density/composition on association.

### Runtime component permutation

The same associating state is rebuilt with reversed component order. Density roots and fugacity coefficients must map exactly back by component identity.

## Historical gate boundary

When this PT-property gate first landed, CPA stopped at the phase-property surface and
stability/flash were still future work. Those later layers now exist; current capability
is summarized by [the flash module](../flash/cpa_flash.md) and the repository root
README. The numerical root-search limitations above remain part of the current property
contract.


## PR #115 residual-Helmholtz source-of-truth completion record

This section is the minimal retained audit record for the production source-of-truth
switch completed by PR #115. It replaces the superseded chronological design/audit
document without discarding the scientific acceptance decisions needed to interpret the
current implementation.

For the frozen profile
`CPA/SRK-physical/simplified-rdf-1.9eta/explicit-site-pairs/v1`, production pressure
and residual chemical potentials now come from the canonical scalar-generic residual
Helmholtz terms in `cpa_residual_helmholtz.hpp`, differentiated by the
first-derivative adapter in `cpa_helmholtz_derivatives.hpp`. The association state is
solved once at the primal state and is held stationary during first differentiation.
The switch did not change the association solve, density-root algorithm or tolerances,
stability/split/max3 algorithms, parameters, applicability, failure semantics or fallback
behavior.

The final acceptance matrix was:

| Gate | Final status |
| --- | --- |
| A — scalar residual-Helmholtz value | **PASS** |
| B — pressure versus legacy analytic path | **PASS** |
| C — residual chemical potentials versus legacy analytic path | **PASS** |
| D — final `ln(phi)` plus pinned ThermoPack parity | **PASS** |
| E — named independent implementations | **WAIVED / NON-BLOCKING** |
| F — Helmholtz production source of truth | **PASS / COMPLETE** |

### Gate-E external numerical-defect witness

The named Clapeyron reference remains pinned to
`ClapeyronThermo/Clapeyron.jl@229b09452f36c2f812486150df0bb43b197bb4e5`.
Its raw ten-state regression is deliberately retained as a **known numerical-defect
witness**, not rewritten as a pass. The raw pinned implementation matched the frozen
formulation/parameters and met the scalar `F_res` and residual-chemical-potential
envelopes, but its exposed pressure / `ln(phi)` path missed the already-frozen
derivative-level envelopes.

The defect-attribution audit traced that miss to the pinned compressed
`X_exact2!` association path returning before configurable convergence settings were
consulted. Solving the same compressed association equations independently to high
stationarity, without changing the CPA formulation, parameters or acceptance thresholds,
reduced the complete ten-state discrepancies to:

- `max |Delta P| = 5.323330668403745e-7 Pa`;
- `max |Delta ln(phi)_Z-only| = 1.3571962033602398e-11`.

The project-owner decision therefore made Gate E **waived / non-blocking** for this
frozen profile; it did **not** claim that the raw pinned Clapeyron output passed. The raw
oracle and diagnostic scripts remain under
`tests/flash/cpa_clapeyron_oracle/`, including the pinned JSON witness and the
independent compressed-stationarity diagnostic.

### Gate-F production cross-check and downstream evidence

The pre-switch hand-written pressure and residual-chemical-potential formulas remain
test-only independent regression oracles. Across ten frozen liquid/vapor states, two
parameter snapshots and both component orders (40 comparisons), the final production
Helmholtz path versus those legacy analytic oracles reached:

- maximum absolute pressure delta: `2.1159648895263672e-6 Pa`;
- maximum absolute residual-chemical-potential delta: `7.5139894306630595e-13`;
- maximum absolute `ln(phi)` delta: `7.5139894306630595e-13`.

No numerical tolerance was widened for the switch. The reviewed PR #115 evidence also
passed the CPA physical-validation, density-root/fugacity, stability, two-phase split,
max-three-phase, model-neutral backend, ThermoPack-parity and standalone thermodynamics
contract suites on their declared hosted matrices.

The paired performance/structural audit (run `35303522137`) compared
`main@f9d65c9e5de03a6ac64fe96f6458faf311ae6fea` with the Gate-F production head
`5ad82d8520f5cf391cfc51d35baaf68fed04c094`. The five-state full-flash median moved
from `24.991 s` to `26.104 s` (ratio `1.0446`) and was classified
`no_clear_hosted_runner_regression_signal`; no speedup or hard hosted-wall-time gate is
claimed. High-level stability/split/root-search counts were unchanged. Lower-level
density/property evaluations, association solves and fixed-point sweeps moved by only
about `0.016-0.017%`, a separately reviewed root-tolerance-boundary sampling effect,
so Gate F is not described as bit-for-bit structural identity.

### Remaining scope limits

PR #115 authorized only the reviewed **first-derivative** Helmholtz source of truth at that time. The later selected-root derivative gate documented above adds the specific second-directional Helmholtz/association information required for local `d ln(phi)/d(p,T,x)` on a fixed smooth PT root. It still does not establish a general public Hessian API, caloric APIs, new CPA formulations or parameter sets, association continuation/caching, a performance optimization, or a physical three-phase VLLE oracle. Finite density-root and finite TPD
searches retain their existing non-global-proof semantics.

## Selected-root first derivatives for natural-variable flow

The selected PT phase now exposes first derivatives of `ln(phi_i)` with respect to arbitrary tangent directions in `(p,T,x)` through `cpa_pt_phase_ad.hpp`. This is a fixed-branch local derivative. It does not perform stability analysis, phase selection or root switching.

The implementation follows the stationary-`Q` treatment used for CPA association derivatives. For the site-fraction stationarity equations

```text
g(X,z) = 0
```

with `z in {T,V,n_i}`, site sensitivities are obtained from

```text
(dg/dX) dX/dz = - dg/dz .
```

The association Jacobian is factorized with partial pivoting and an explicit ill-conditioning guard. No derivative of the fixed-point iteration is taken.

The minimized residual Helmholtz second directional derivative uses

```text
d2F/dz1dz2
  = Q_z1z2 + Q_z1X X_z2 + X_z1^T Q_Xz2
    + X_z1^T Q_XX X_z2
```

at the stationary association state. The selected PT density root then satisfies

```text
P_EOS(T,V,n) - p = 0,
```

so its local volume/density derivative is obtained from the corresponding pressure IFT. Finally,

```text
ln(phi_i) = dF_res/dn_i - ln(Z)
```

is differentiated on the same selected root.

This route is consistent with the ThermoPack CPA memo's stationary-`Q` Hessian construction and Michelsen-style association derivative treatment:
`https://github.com/thermotools/thermopack/blob/main/docs/memo/CPA/cpa.tex`.

The derivative path rejects:

- non-tangent mole-fraction derivative directions;
- near-multiple/tangent PT roots;
- an ill-conditioned pressure root;
- singular/ill-conditioned association Jacobians;
- non-finite association, pressure or fugacity derivatives.

The production derivative contains no finite-difference fallback. Focused tests compare the analytic/IFT derivative against fresh independently re-solved PT perturbations of the associating binary structural fixture. Those perturbations are test-only numerical cross-checks, not production differentiation and not physical validation.
