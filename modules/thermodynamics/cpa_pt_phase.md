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

## Current boundary

After this gate CPA has the thermodynamic phase-property surface required by a future stability evaluator, but it still does **not** provide:

- a generic `StabilityPhase` adapter;
- TPD search acceptance;
- two-phase split;
- max-three-phase orchestration;
- phase-set publication/backend integration;
- implicit root/flash sensitivities;
- physical CPA parameter database or experimental regression.

The next gate should adapt only mechanically admissible resolved CPA roots into the existing generic PT stability contract, compare same-composition roots by Gibbs energy, and establish independent single-phase/unstable stability regressions before any CPA phase split is implemented.
