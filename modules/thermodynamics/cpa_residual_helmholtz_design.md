# CPA residual Helmholtz kernel design

## Status and scope

This document resets the active CPA development line to the **thermodynamic model itself**.
The association-continuation / certificate work remains preserved as audit evidence but is
deferred and is not a prerequisite for completing CPA. No production continuation or
certificate optimization is enabled by this design.

This slice audits the current hand-written CPA pressure and fugacity-coefficient formulas
and defines the next production thermodynamic kernel around one canonical residual
Helmholtz potential.

Current production profile under audit:

```text
CPA/SRK-physical/simplified-rdf-1.9eta/explicit-site-pairs/v1
```

This design does **not** change the CPA formulation, parameter snapshot, association
scheme, root solver, stability solver, flash solver, numerical tolerances or failure
semantics.

## 1. Primary references

The model definition is taken from the CPA literature; external software is used as
independent implementation evidence rather than as a source of fitted golden values.

1. G. M. Kontogeorgis et al., *Ten Years with the CPA (Cubic-Plus-Association)
   Equation of State. Part 1. Pure Compounds and Self-Associating Systems*,
   Ind. Eng. Chem. Res. 45 (2006) 4855-4868, DOI `10.1021/ie051305v`.
2. G. M. Kontogeorgis et al., *Ten Years with the CPA ... Part 2. Cross-Associating
   and Multicomponent Systems*, Ind. Eng. Chem. Res. 45 (2006) 4869-4878,
   DOI `10.1021/ie051306n`.
3. ThermoPack pinned reference `thermotools/thermopack@d68c794c7342bfc6938eb424a1fbb88b7780b738`,
   especially `docs/memo/CPA/cpa.tex`. Its memo formulates CPA as
   `A^CPA = A^ideal + A^SRK + A^assoc` and derives first derivatives using the
   stationary Michelsen/CPA `Q` function.
4. Clapeyron.jl pinned source
   `ClapeyronThermo/Clapeyron.jl@229b09452f36c2f812486150df0bb43b197bb4e5`.
   Its CPA model is organized as `a_res(cubic) + a_assoc`; when association fractions
   are solved in primal space it explicitly evaluates the `Q` function to recover
   derivative information. Its general association solver also supports implicit AD.
5. NIST `teqp` CPA documentation as a secondary formula cross-check. It publishes the
   residual Helmholtz split and the SRK/PR cubic residual terms in a common framework.

No implementation code from these projects is copied into MPMC_HNU.

## 2. Current MPMC_HNU implementation under audit

The current production path is split across:

- `cpa_phase.hpp`
  - classic CPA alpha function `a_i(T)`;
  - classical `a_mix` and linear `b_mix`;
  - density-state association solve;
  - hand-written physical and association pressure.
- `cpa_association.hpp`
  - simplified radial distribution function;
  - explicit-site mass-action equations;
  - damped fixed-point solve.
- `cpa_pt_phase.hpp`
  - density roots;
  - hand-written cubic and association residual chemical potentials;
  - `ln(phi)` construction.

The audit question is not whether these routines reproduce the existing tests; they do.
The question is whether their equations are mutually derivable from one CPA thermodynamic
potential.

## 3. Canonical extensive variables

The new kernel must use **extensive mole numbers**, not normalized composition, as its
independent composition variables.

Define

```text
n       = sum_i n_i
x_i     = n_i / n
rho     = n / V
B(n)    = sum_i n_i b_i = n b_mix
A(T,n)  = sum_i sum_j n_i n_j a_ij(T) = n^2 a_mix
u       = B / V = b_mix rho
```

with

```text
a_ij(T) = sqrt(a_i(T) a_j(T)) (1-k_ij)
```

and

```text
a_i(T) = a0_i [1 + c1_i (1-sqrt(T/Tc_i))]^2.
```

Using `(V,n_i)` rather than `(rho,x_i)` is essential: chemical potentials are partial
derivatives at fixed `T,V,n_{j!=i}`. A Jacobian with respect to normalized mole fractions
cannot directly represent those derivatives.

For the existing density-state APIs, the extensive state can be chosen without loss of
generality as

```text
n_i = x_i mol
n   = 1 mol
V   = 1/rho m^3
```

because the residual Helmholtz energy is extensive.

## 4. Cubic residual Helmholtz contribution

For the current SRK physical contribution, define the dimensionless extensive residual
Helmholtz function

```text
F_cubic(T,V,n) = A_cubic^res/(R T)
```

as

```text
F_cubic = -n ln(1-B/V)
          - A(T,n)/(B R T) ln(1+B/V).
```

Dividing by total moles gives the familiar molar expression

```text
alpha_cubic^res = -ln(1-b rho)
                  - a_mix/(b_mix R T) ln(1+b_mix rho).
```

This agrees with the CPA literature and the independent `teqp` CPA formulation for SRK.

### 4.1 Pressure recovered from the potential

The total pressure is

```text
P = n R T / V - R T (partial F_res / partial V)_(T,n).
```

Applying this to `F_cubic` gives exactly

```text
P_physical = R T rho/(1-b_mix rho)
             - a_mix rho^2/(1+b_mix rho),
```

which is the expression currently implemented in `cpa_phase.hpp`.

### 4.2 Cubic chemical potential recovered from the potential

At fixed `T,V,n_{j!=i}`,

```text
mu_i^res/(R T) = partial F_res / partial n_i.
```

Differentiating `F_cubic` gives the current `cpa_pt_phase.hpp` expression

```text
mu_i,cubic^res/(R T)
  = (b_i/b_mix)(Z_physical-1)
    - ln(1-b_mix rho)
    - a_mix/(b_mix R T)
      [2 sum_j x_j a_ij/a_mix - b_i/b_mix]
      ln(1+b_mix rho),
```

where

```text
Z_physical = P_physical/(rho R T).
```

**Audit result:** the present hand-written cubic pressure and residual chemical-potential
formulas are mutually consistent with the same SRK residual Helmholtz potential.

## 5. Association residual Helmholtz contribution

For each association site class `A_i`, let `X_Ai` be the unbonded fraction. The CPA
association contribution is

```text
F_assoc = A_assoc^res/(R T)
        = sum_i n_i sum_Ai [ln X_Ai - X_Ai/2 + 1/2].
```

The mass-action equations are

```text
X_Ai = 1 /
       [1 + (1/V) sum_j n_j sum_Bj X_Bj Delta_AiBj].
```

For the current simplified/Kontogeorgis radial distribution function,

```text
g = 1 / (1 - 0.475 B/V)
  = 1 / (1 - 1.9 eta),

eta = B/(4V) = b_mix rho/4,
```

and

```text
Delta_AiBj = g [exp(epsilon_AiBj/(R T))-1] b_ij beta_AiBj.
```

This is algebraically the same formulation used by the current association kernel.
The current code's

```text
rho_dln_g_drho = (1.9 eta)/(1-1.9 eta)
```

is exactly `rho * d ln(g)/d rho` for this `g`.

## 6. The Q function is the correct first-derivative bridge

A direct AD evaluation of `F_assoc` while treating converged `X` as constants would be
wrong, because `F_assoc` is not stationary with respect to `X` in that representation.
Conversely, propagating Dual numbers through the iterative association solver is not
necessary for first derivatives.

Use the stationary CPA/Michelsen function

```text
Q(T,V,n,X)
 = sum_i sum_Ai n_i [ln X_Ai - X_Ai + 1]
   - (1/(2V))
     sum_i sum_j sum_Ai sum_Bj
       n_i n_j X_Ai X_Bj Delta_AiBj.
```

At the converged association state `X*`, the mass-action equations are equivalent to

```text
partial Q / partial X = 0,
```

and

```text
F_assoc(T,V,n) = Q(T,V,n,X*(T,V,n)).
```

Therefore the envelope/stationarity relation gives, for every first derivative
`z in {T,V,n_i}`,

```text
partial F_assoc / partial z
  = [partial Q / partial z]_(X=X*)
```

without requiring `partial X/partial z`.

This is exactly the structure documented by the pinned ThermoPack CPA memo. Clapeyron.jl
uses the same strategy when its association fractions are solved in primal space: its
source comments state that `X` does not carry derivative information and that Michelsen's
`Q` function is evaluated to recover the derivatives.

### Consequence for MPMC_HNU

The first residual-Helmholtz implementation should:

1. solve association once in ordinary `double` at the primal `(T,V,n)` state;
2. retain the converged `X*` as constants;
3. evaluate `F_cubic + Q_assoc` with AD-active `V` and `n_i`;
4. recompute `Delta(T,V,n)` inside `Q_assoc` with AD-active variables;
5. **not** differentiate through the fixed-point iteration.

This is both simpler and better aligned with the current independent AD module.

## 7. Audit of the current association pressure

The ThermoPack `Q` derivative gives

```text
partial F_assoc/partial V
 = (1/(2V)) sum n_i n_j X_Ai X_Bj
   [Delta/V - partial Delta/partial V].
```

For the current simplified `g`, this reduces using the mass-action identities to

```text
P_assoc = -0.5 R T rho
          [1 + rho d ln(g)/d rho]
          sum_i x_i sum_Ai (1-X_Ai).
```

That is exactly the current `cpa_phase.hpp` implementation:

```text
-0.5 * R*T * rho
* (1 + rho_dln_g_drho)
* association_sum.
```

**Audit result:** no model-form error was found in the current hand-written association
pressure for the present sCPA profile.

## 8. Audit of the current association chemical potential and ln(phi)

From the stationary `Q` function,

```text
partial F_assoc/partial n_k
 = sum_Ak ln X_Ak
   - (1/(2V)) sum n_i n_j X_Ai X_Bj
     partial Delta_AiBj/partial n_k.
```

For

```text
g = 1/(1-0.475 B/V),
```

and the present linear `B=sum n_i b_i`, this reduces to

```text
mu_k,assoc^res/(R T)
 = sum_Ak ln X_Ak
   - (1.9/8) rho b_k g
     sum_i x_i sum_Ai (1-X_Ai).
```

This is exactly the association term currently implemented by `cpa_fill_ln_phi`.

For the complete residual chemical potential,

```text
mu_i^res/(R T)
 = partial (F_cubic + F_assoc)/partial n_i.
```

At the solved total CPA pressure,

```text
Z = P/(rho R T)
```

and

```text
ln(phi_i) = mu_i^res/(R T) - ln Z.
```

The current code computes this as

```text
ln(phi_i) = mu_cubic + mu_association - ln(Z),
```

which is thermodynamically consistent with the same residual Helmholtz potential.

**Audit result:** no formula-level inconsistency was found in the current hand-written
`ln(phi)` for the present SRK+sCPA+explicit-site-pair model.

## 9. What is actually wrong with the current architecture

The current formulas are correct for the frozen profile, but the implementation has
multiple independent algebraic sources of truth:

- `pressure_physical_pa` is hand written;
- `pressure_association_pa` is hand written;
- `mu_cubic` is hand written separately;
- `mu_association` is hand written separately;
- the association model itself is implemented in a third location.

That is manageable for one formulation but fragile for future extensions:

- original CPA Carnahan-Starling RDF;
- simplified CPA RDF;
- alternate cubic physical terms such as PR-CPA;
- temperature-dependent `k_ij`;
- additional cross-association rules;
- solvation-specific association strengths;
- higher thermodynamic derivatives.

A new formulation can easily update pressure while forgetting the matching composition
derivative in `ln(phi)`.

The residual Helmholtz kernel should therefore become the canonical thermodynamic source;
the existing formulas should initially remain as independent regression oracles.

## 10. Proposed kernel boundary

The first production increment should introduce a small scalar-generic kernel, conceptually:

```cpp
template<class Scalar>
Scalar cpa_cubic_residual_helmholtz_reduced(
    Scalar temperature_k,
    Scalar volume_m3,
    std::span<const Scalar> mole_numbers,
    const CpaParameterSet& parameters);

template<class Scalar>
Scalar cpa_association_q_reduced(
    Scalar temperature_k,
    Scalar volume_m3,
    std::span<const Scalar> mole_numbers,
    const CpaParameterSet& parameters,
    const CpaAssociationResult& primal_association);

template<class Scalar>
Scalar cpa_residual_helmholtz_reduced(
    Scalar temperature_k,
    Scalar volume_m3,
    std::span<const Scalar> mole_numbers,
    const CpaParameterSet& parameters,
    const CpaAssociationResult& primal_association);
```

The name `reduced` here means

```text
F_res = A_res/(R T),
```

which is **extensive in moles**, not the per-mole Helmholtz energy.

Exact public naming is intentionally deferred to the implementation increment; the
important contract is the mathematical one.

### Scalar-generic, not AD-owned

The thermodynamics kernel should remain scalar-generic and should not make the parameter
contract depend on AD. Standard elementary operations should be written in a form that
can be instantiated by `double` and the existing `mpmc::ad::Dual` numeric type. The AD
driver belongs in a thin derivative adapter/test layer.

This avoids turning the thermodynamics module into an AD framework while still reusing
MPMC_HNU's existing block-forward runtime Jacobian machinery.

## 11. First-derivative adapter using the existing AD module

For a density state `(T,rho,x)`:

```text
n_i = x_i
V   = 1/rho
```

First solve the existing association equations at this primal state. Then capture the
converged `X*` and evaluate the scalar residual Helmholtz kernel with inputs

```text
[V, n_0, ..., n_(N-1)].
```

`value_and_jacobian_runtime<K>` already supports runtime component counts while keeping a
fixed `Dual<T,K>` propagation width.

For the returned scalar Jacobian:

```text
F_V   = dF_res/dV
F_ni  = dF_res/dn_i
```

compute

```text
P = n R T / V - R T F_V

mu_i^res/(R T) = F_ni

Z = P V/(n R T)

ln(phi_i) = F_ni - ln(Z).
```

This path has an important practical advantage: association is solved **once at the primal
state**, outside the repeated block-AD callbacks. The callback contains only deterministic
algebra with frozen `X*`, which satisfies the AD driver's repeated-evaluation contract.

## 12. Validation ladder before any replacement of production formulas

The current hand-written path must not be deleted immediately. It is useful as an
independent regression oracle while the new kernel is introduced.

### Gate A — scalar Helmholtz value

At fixed physical states verify independently:

- cubic `F_cubic` against its analytic formula;
- direct association
  `sum n_i sum_A (ln X - X/2 + 1/2)` against `Q(T,V,n,X*)`;
- total `F_res = F_cubic + F_assoc`.

### Gate B — AD pressure versus current analytic pressure

Across the existing physical methanol-water states and a wider density lattice:

```text
P_AD = nRT/V - RT dF_res/dV
```

must agree with the current `evaluate_cpa_phase_at_density(...).pressure_pa` at a
predeclared roundoff-scale tolerance.

Compare the decomposition as well:

- SRK physical pressure;
- association pressure;
- total pressure.

### Gate C — AD residual chemical potentials versus current analytic terms

Check every component separately:

```text
F_ni
```

against

```text
mu_cubic + mu_association
```

from the existing analytic implementation.

This is more diagnostic than comparing only final `ln(phi)`.

### Gate D — final ln(phi)

Using the same density root and target pressure,

```text
ln(phi_i)_AD = F_ni - ln Z
```

must match the current `cpa_fill_ln_phi` result and the already frozen ThermoPack
phase-kernel oracle.

### Gate E — independent external implementations

For matched formulation/parameters/states, compare the new kernel with:

- pinned ThermoPack;
- pinned Clapeyron.jl;
- literature formulas/fixtures.

External software parity must remain distinct from experimental VLE accuracy.

### Gate F — only then switch production source of truth

Only after A-E pass on GCC/Clang/MSVC should pressure and `ln(phi)` be redirected to the
new Helmholtz-derived path. The old analytic expressions can then remain as focused
cross-check tests or be removed in a later cleanup after a separate review.

## 13. Second derivatives are explicitly later

The stationarity trick removes `dX/dz` only for **first derivatives**. Second derivatives
require implicit association sensitivity, as the ThermoPack memo derives:

```text
d2F/dz2 dz1
 = Q_z1z2 - Q_z1X (Q_XX)^(-1) Q_Xz2.
```

Clapeyron also has an implicit-AD route for its association root solve.

MPMC_HNU should not add that complexity in the first residual-Helmholtz increment. Pressure
and fugacity require only first derivatives and can be made correct and unified first.

## 14. Explicit non-goals for the next increment

The first implementation must not:

- change the fixed-point association solver;
- enable association continuation/caching;
- alter `1e-12` association convergence tolerance;
- alter density-root scan/tolerances;
- change stability or flash algorithms;
- change pure/binary/cross-association parameters;
- introduce a new RDF or PR-CPA formulation;
- claim a performance improvement;
- fit external-software output.

## 15. Audit conclusion and recommended next increment

The current sCPA `pressure` and `ln(phi)` formulas are scientifically consistent with the
same CPA residual Helmholtz model. The issue to solve is **single-source thermodynamic
architecture**, not a discovered formula defect.

The next small production increment should therefore be:

> Implement a scalar-generic, extensive `F_res=A_res/(RT)` density-state kernel consisting
> of SRK residual Helmholtz plus the stationary association `Q` contribution, and add a
> test-only AD derivative adapter that cross-checks `P`, component residual chemical
> potentials and `ln(phi)` against the existing hand-written formulas. Do not redirect the
> production root/flash path yet.

That increment is small, scientifically auditable and directly reusable by later CPA
variants.