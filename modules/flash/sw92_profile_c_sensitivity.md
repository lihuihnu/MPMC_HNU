# SW92 Profile-C fixed-phase-set implicit sensitivity

## Scope

This increment differentiates an already authoritative SW92 Profile-C PT phase set under a **strictly local, fixed discrete state**:

- phase count is fixed;
- each phase keeps its accepted AQ/NA thermodynamic family;
- each phase keeps its accepted selected cubic-root branch;
- phase representation slots are fixed locally;
- prescribed NaCl molality is fixed.

It does **not** differentiate phase selection, TPD search, SSI, generalized Rachford-Rice iteration, line search, topology routing, phase appearance/disappearance, H-slot morphology classification, or a global-stability proof.

Public entry:

```cpp
#include <mpmc/flash/sw92_profile_c_sensitivity.hpp>

const auto derivative =
    mpmc::flash::differentiate_sw92_profile_c_phase_set(
        authoritative_phase_set, sw92_model);
```

The input must be a `Sw92ProfileCPtPhaseSetResult` already accepted by the authoritative Profile-C publication contract.

## External coordinates and outputs

The public independent coordinates are

```text
q = (p_Pa, T_K, z_0, ..., z_{N-2})
z_{N-1} = 1 - sum(z_0, ..., z_{N-2})
```

For every accepted phase `alpha`, the result publishes local derivatives of

```text
beta_alpha
x_{alpha,i}
Z_alpha
c_alpha = p/(Z_alpha R T)
```

with respect to every column of `q`.

The result also preserves the accepted dataset/revision/component ordering, fixed molality, model/equilibrium/orchestration/boundary/publication conventions, and SW92 Profile-C phase metadata. H0/H1 remain unordered `nonaqueous_unclassified` NA phase instances. A derivative never promotes them to liquid/vapor or LV/LL.

## Fixed-phase-set implicit system

For `P` accepted phases and `N` components, the local equilibrium manifold has `PN-1` internal unknowns:

- `N-1` independent composition coordinates for every phase;
- `P-1` independent phase-fraction coordinates.

The local residual contains exactly `PN-1` equations:

```text
(P-1) * N  common reduced chemical-potential equalities
N-1        component material balances
```

For phase `alpha > 0` and component `i`, the chemical-potential residual is

```text
F_mu(alpha,i)
  = log(x_0,i) + log(phi_0,i)
  - log(x_alpha,i) - log(phi_alpha,i)
```

and the material-balance rows are

```text
F_z(i) = sum_alpha(beta_alpha * x_alpha,i) - z_i,
         i = 0,...,N-2.
```

At the accepted base point, the implementation checks that these equations reproduce the authoritative state before any derivative is published.

## Internal log-ratio chart

The public derivative remains with respect to `q` above, but the IFT uses a numerically conditioned internal chart.

For each accepted phase, the largest accepted mole fraction is chosen as a local reference component `r_alpha`. Internal composition variables are

```text
eta_alpha,i = log(x_alpha,i / x_alpha,r_alpha), i != r_alpha.
```

A softmax-style reconstruction maps these variables back to a strictly positive, exactly normalized composition.

The largest accepted phase fraction is likewise chosen as a local reference phase `r_beta` and the internal phase-fraction variables are

```text
tau_alpha = log(beta_alpha / beta_r_beta), alpha != r_beta.
```

The phase fractions are reconstructed by the same normalized exponential map.

This chart is important for the physical Sample-6 state: the aqueous phase contains trace hydrocarbon components down to approximately `1e-11`. A raw mole-fraction chart would create both catastrophic subtraction when a trace component is dependent and `1/x` Jacobian scaling when a trace component is independent. The log-ratio chart avoids both pathologies without changing the physical equilibrium equations or public derivative coordinates.

## AD and selected-root derivative path

The existing corrected-original SW92 thermodynamic kernel is extended so the already implemented formulas can accept the same structural forward-mode AD number contract used by the PR76 thermodynamic path.

AD propagates through:

- SW92 water alpha;
- temperature-dependent AQ water-pair correlations;
- H2S non-aqueous Eq. (17);
- PR pure-component alpha for non-water components;
- AQ/NA mixing;
- `Z` and `ln(phi)` evaluation.

The cubic root iteration itself is **not** differentiated. The selected accepted simple root is differentiated by the local cubic implicit-function relation, reusing the existing PR root derivative contract. If the selected root is near multiple or has no reliable local derivative, sensitivity is unavailable.

No thermodynamic formula, BIP, root-selection rule, parameter, or primal flash tolerance is changed by making this arithmetic AD-capable. The historical floating-point `Sw92PhaseValues<double>` public spelling is retained.

## Implicit differentiation

Let the accepted local equations be

```text
F(u, q) = 0.
```

Forward AD evaluates both

```text
A = partial F / partial u
B = partial F / partial q
```

at the authoritative base point. The local state derivative is obtained from

```text
du/dq = -A^{-1} B.
```

The implementation reuses the repository's audited dense LU/rcond/backward-error machinery from the PR76 sensitivity infrastructure. It does not form an explicit matrix inverse.

The same AD pass also produces partial derivatives of `beta`, full phase compositions, `Z`, and molar density with respect to `(u,q)`. Their total derivatives are assembled by the chain rule using `du/dq`.

After assembly, differentiated identities are checked independently:

```text
sum_alpha d beta_alpha/dq = 0
sum_i d x_alpha,i/dq = 0
sum_alpha [d beta_alpha/dq * x_alpha,i
           + beta_alpha * d x_alpha,i/dq]
  = d z_i/dq.
```

A failed identity check is an arithmetic failure, not a silently clipped derivative.

## Availability and failure policy

Sensitivity is published only inside a smooth fixed discrete state. The adapter rejects or downgrades the derivative for conditions including:

- authoritative phase set absent or inconsistent;
- zero overall feed support in the current v1 reduced-feed chart;
- zero/invalid accepted phase composition support;
- a phase fraction inside the derivative boundary guard;
- phase compositions inside the local coalescence/separation guard;
- accepted activity/Z no longer reproducible from the supplied model snapshot;
- selected cubic root derivative unavailable;
- local equilibrium Jacobian singular or ill-conditioned;
- implicit linear solve failing its backward-error check;
- nonfinite total derivatives or differentiated identity failure.

No zero derivative, previous-state derivative, finite-difference production fallback, phase clipping, or hidden topology switch is used.

## Validation

The focused regression includes:

1. **AD phase-property cross-check.** Selected SW92 `Z` and `ln(phi)` derivatives with respect to `p`, `T`, and reduced composition are generated by AD and compared against an independent centered numerical perturbation used only in the test.
2. **Single-phase fixed-set identity.** The authoritative dry binary state verifies `beta=1`, the reduced-feed composition identity, and finite `Z`/density derivatives.
3. **Two-phase fresh re-solve cross-check.** A wet CO2/H2O state compares the implicit derivatives of phase fractions, phase compositions, `Z`, and molar densities against independently fresh authoritative solves at symmetric perturbations in `p`, `T`, and reduced feed.
4. **Physical Sample-6 three-phase state.** The 8-component Mortezazadeh-Rasaei state exercises the full `W(AQ)+H0(NA)+H1(NA)` implicit system with trace compositions, physical three-phase fractions, and the independently regenerated Decimal(80) primal oracle.
5. Runtime component permutation checks for the invariant `p/T` columns.
6. H0/H1 slot-swap symmetry.
7. Phase-boundary, tampered-property, and ordered-model mismatch guards.
8. Public-header self containment and historical `Sw92PhaseValues<double>` source compatibility.

The dedicated hosted matrix uses GCC Debug + ASan/UBSan, Clang Release, and MSVC Release and also reruns affected SW92 thermodynamics and authoritative Profile-C publication regressions.

## Remaining boundary

This increment creates the SW92 fixed-phase-set flash Jacobian contract. It still does not provide:

- derivatives across phase appearance/disappearance;
- derivatives across AQ/NA family changes or cubic-root switching;
- a critical/near-multiple-root derivative;
- a validated hydrocarbon LV/LL morphology derivative;
- NaCl inventory/salinity derivatives (molality remains prescribed);
- a mathematical global-stability derivative/proof;
- direct physics/Newton consumption of this Jacobian.

Wiring the successful sensitivity into the SW92 physics thermodynamic-closure snapshot is a separate downstream adapter gate so residual availability and Newton-seed availability remain independently auditable.
