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

### 12.1 Gate-F readiness audit — PR #115 head `bf519bb40aa1d70fe09ab5b7a315c56909034373`

This audit maps the original A-E prerequisites to concrete current-head evidence. A gate
is marked `PASS` only when the requirement as written above is directly exercised; later
or stronger evidence does not waive an explicitly missing prerequisite.

| Gate | Status | Current-head evidence | Remaining blocker |
| --- | --- | --- | --- |
| A — scalar Helmholtz value | **PARTIAL / BLOCKED** | `cpa_residual_helmholtz_ad_test.cpp` checks stationary `Q` against the direct association Helmholtz expression and checks `F_res = F_cubic + Q`; `cpa_residual_helmholtz_limits_test.cpp` adds pure, full-binary pure endpoints, zero-association, dilute, permutation and extensivity limits. | There is no test-only scalar oracle that evaluates the SRK `F_cubic` analytic expression independently of `cpa_cubic_residual_helmholtz_reduced`. Derivative agreement does not substitute for the scalar-value gate. |
| B — pressure | **PASS** | The AD regression compares physical, association and total pressure against the current analytic density-state path over all ten frozen liquid/vapor phase states, both literature and ThermoPack-parity parameter snapshots, and normal/swapped component order. Limits additionally cover pure components, binary pure endpoints, zero association, dilute states and homogeneous scaling. | None for the frozen profile. |
| C — residual chemical potentials | **PASS** | Every component's cubic, association and total `mu^res/(RT)` is compared with the independent hand-written analytic path over the same ten-state/two-snapshot/two-order matrix; endpoint and extensivity limits add independent edge coverage. | None for the frozen profile. |
| D — `ln(phi)` | **PASS** | At the same density and target pressure, Helmholtz-derived `F_ni - ln Z` is checked against `cpa_fill_ln_phi`; the matched ThermoPack phase-kernel oracle additionally checks all ten phase states under the frozen external thresholds. | None for the frozen profile. |
| E — independent implementations | **PARTIAL / BLOCKED** | The pinned ThermoPack `d68c794...` phase-kernel oracle is frozen and reproduced; the Helmholtz AD test compares pressure decomposition, residual-chemical-potential decomposition and `ln(phi)` against it. The literature parameter/VLE fixture is source-traceable to Kontogeorgis/Folas/Kurihara and the design algebra is tied to the published CPA formulas. | The design explicitly requires a numeric comparison with pinned `ClapeyronThermo/Clapeyron.jl@229b09452f36c2f812486150df0bb43b197bb4e5`. No frozen Clapeyron oracle or current-head regression exists. `teqp` formula cross-checks do not satisfy this named implementation requirement. |

Supplemental first-derivative evidence is already stronger than the original A-D minimum:
`cpa_residual_helmholtz_temperature_test.cpp` independently re-solves association at
`T +/- 0.005 K` and confirms the stationary-Q temperature derivative from vapor-like,
moderate-density and liquid-like states. This supports the first-derivative architecture
but does not change the two blockers above.

**Readiness verdict:** Gate F is **NOT READY** at this head. The only open prerequisites
identified by this audit are (1) the missing independent scalar `F_cubic` value oracle in
Gate A and (2) the missing pinned Clapeyron numeric comparison in Gate E. Pressure,
chemical-potential and `ln(phi)` equivalence do not need additional formula work before
those two blockers are closed.

### 12.2 Frozen Gate-F production-switchover acceptance contract

The following contract must be satisfied **before** any production source-of-truth switch.
It applies only to
`CPA/SRK-physical/simplified-rdf-1.9eta/explicit-site-pairs/v1`; it must not be silently
reused as acceptance for PR-CPA, another RDF or another association rule.

#### A. Close the scalar-value blocker without changing production code

Add a test-only independent SRK scalar reference that evaluates

```text
F_cubic,ref = -n ln(1-B/V)
              - A(T,n)/(B R T) ln(1+B/V)
```

without calling `cpa_cubic_residual_helmholtz_reduced`. Exercise at minimum the ten frozen
ThermoPack liquid/vapor phase states in normal and swapped component order, plus the
existing full-binary pure endpoints. Freeze the scalar acceptance at the existing
roundoff convention

```text
roundoff(s) = 4096 * eps(double) * max(1, abs(s))

abs(F_cubic - F_cubic,ref)
  <= roundoff(F_cubic) + roundoff(F_cubic,ref).
```

The already-frozen association scalar gate remains unchanged:

```text
abs(Q_assoc - F_assoc,direct)
  <= 1e-10 + 1e-12 * max(abs(Q_assoc), abs(F_assoc,direct)).
```

No result obtained after implementation may be used to widen either contract.

#### B. Close the named Clapeyron leg of Gate E

Create a reproducible test oracle from the **unmodified** pinned revision
`ClapeyronThermo/Clapeyron.jl@229b09452f36c2f812486150df0bb43b197bb4e5`.
The generator must record the exact revision and explicitly align the same SRK physical
term, simplified CPA radial-distribution formulation, 2B methanol / 4C water association
topology, CR-1 cross association, pure parameters, `kij=-0.055`, units, gas constant and
component ordering. If the pinned implementation cannot represent the current MPMC_HNU
profile exactly, Gate E remains blocked; an approximately similar Clapeyron model must not
be accepted as a substitute.

Use the same ten frozen phase states `(T,V,n)` already carried by the ThermoPack
phase-kernel oracle. At `n=1 mol`, the predeclared cross-implementation envelope is:

```text
abs(F_res_MPMC - F_res_Clapeyron) <= 1e-10
abs(P_MPMC - P_Clapeyron)         <= 5e-6 Pa
abs(mu_i^res/(RT)_MPMC - mu_i^res/(RT)_Clapeyron) <= 1e-10
abs(ln(phi_i)_MPMC - ln(phi_i)_Clapeyron)          <= 1e-10
```

If Clapeyron does not expose one of these quantities directly, the oracle may derive it
from its own public residual-Helmholtz API, but it must not call MPMC_HNU code to create the
reference. Pressure comparison requires an explicitly reconciled gas constant; otherwise
that pressure row is not considered passed and the threshold is not relaxed.

The existing ThermoPack threshold contract remains independently frozen as
`MPMC_HNU/CPA/ThermoPack-parity-thresholds/v1` and is not replaced by the Clapeyron gate.

#### C. Production switch points and invariants

After A and E are fully green, one focused production change may redirect only these
existing observables:

1. `evaluate_cpa_phase_at_density(...)`:
   - retain the existing primal association solve and its failure semantics;
   - obtain `pressure_physical_pa`, `pressure_association_pa` and `pressure_pa` from
     derivatives of the canonical residual Helmholtz terms at the same `(T,V,n)` state.
2. `cpa_fill_ln_phi(...)`:
   - obtain every component's residual chemical potential from `F_ni`;
   - keep `Z = target_pressure/(rho R T)` and compute `ln(phi_i)=F_ni-ln(Z)`;
   - preserve `compressibility_factor`, `reduced_gibbs_offset` and all validation/error
     semantics.

Association must still be solved once at the primal state and `X*` held stationary during
first differentiation. The switch must not alter association convergence tolerance,
density-root scan or tolerances, root identity/status, parameter data, stability search,
flash algorithms, material-balance/fugacity/TPD acceptance, `global_stability_proven`
semantics, or introduce a silent fallback to the old hand-written formulas.

The old analytic pressure and chemical-potential expressions must remain available as
**test oracles** for the switchover regression. Removing them is a later cleanup task and
must not be bundled into the source-of-truth switch.

#### D. Required verification on the switchover commit

The switchover commit is acceptable only if all of the following hold without widening any
frozen threshold:

- Gate A through Gate E are all `PASS` on the exact switchover head;
- `CPA associating physical validation` passes on GCC, Clang and MSVC, including the
  residual-Helmholtz AD, limits and temperature regressions;
- the frozen ThermoPack phase and flash parity thresholds remain green;
- affected downstream CPA workflows pass: density-root/fugacity, stability, vapor-liquid
  split, max-three-phase orchestration and model-neutral PT-backend conformance;
- the existing pure/endpoints/zero-association/dilute/permutation/extensivity regressions
  remain green;
- no new fallback changes an association/root/property failure into success;
- for the fixed five-state performance workload, deterministic root-search calls,
  density-state phase evaluations, association-solve calls and association fixed-point
  sweep counts remain identical unless a separately reviewed algorithmic change explicitly
  changes them.

Performance remains a separate acceptance dimension. Re-run the existing same-runner
`main`/head CPA performance audit on the switchover head. The audit must use its current
interpretation rules: hosted timing is descriptive rather than a newly invented hard
wall-time threshold, but a reproducible regression signal must be investigated before the
switch is accepted. No speedup claim is implied by Gate F.

#### E. What Gate F does not authorize

Passing this contract authorizes only the pressure / residual-chemical-potential / `ln(phi)`
single-source-of-truth switch for the frozen SRK+sCPA profile. It does **not** authorize
second temperature derivatives, caloric-property APIs, association continuation/caching,
new formulations, new parameters, wider applicability claims, three-phase physical-oracle
claims or deletion of the analytic regression path.


### 12.3 Gate A / Gate E follow-up — current PR #115 evidence

This follow-up supersedes the readiness statuses recorded in section 12.1 while preserving
that earlier audit as historical evidence.

**Gate A is now PASS.** Commit
`cc95923f2a0fbcf0ffb6c452bc6c8ce4e64c32b4` added the independent analytic SRK scalar
reference required by the frozen contract. It exercises 40 frozen phase-state comparisons
and 8 full-binary pure-endpoint comparisons with normal/swapped component order. The
observed maximum scalar cubic difference is zero under the predeclared roundoff contract.

**The pinned Clapeyron formulation is exactly representable, but Gate E remains BLOCKED by
the frozen numerical envelope.** The external generator checks out the unmodified
`ClapeyronThermo/Clapeyron.jl@229b09452f36c2f812486150df0bb43b197bb4e5`
source and read-backs the instantiated model before producing any oracle values. The
successful audit confirms:

- RK/SRK cubic with `sCPAAlpha`;
- simplified/Kontogeorgis `g=1/(1-1.9 eta)`;
- `vdW1fRule` and no volume translation;
- methanol 2B (`H=1,e=1`) and water 4C (`H=2,e=2`);
- the frozen ThermoPack-parity `Tc/a0/b/c1` snapshot and `kij=-0.055`;
- `R=8.31446261815324 J mol^-1 K^-1`;
- source-complete CR-1 cross values
  `epsilon_cross=20623 J/mol` and
  `beta_cross=0.03337843615270194`.

The pinned Clapeyron revision has a relevant implementation detail: its native
`:cr1` recombination path calls the non-mutating `epsilon_assoc_mix` from
`assoc_mix!` and `recombine_assoc!` discards that returned mixed epsilon object. A
model constructed from pure association records alone therefore reads back zero
cross-association epsilon. The oracle does **not** patch the pinned source. Instead it
injects the already source-complete CR-1 cross records explicitly and uses
`combining=:nocombining`, exactly matching MPMC_HNU's explicit-site-pair parameter
contract. The read-back maximum parameter delta is
`6.245004513516506e-17`.

The final frozen ten-state oracle is
`tests/flash/cpa_clapeyron_oracle/clapeyron_229b094_meoh_h2o_phase_kernel_33315k.json`.
The generator uses Julia 1.10.12 and, for diagnostic purposes only, tightens its own
association numerics to `rtol=atol=1e-16`, `max_iters=4096`, and implicit AD. Even
with those oracle-only settings, the maximum reported Clapeyron association-equation
residual is `2.940869769929577e-12`.

The C++ Gate-E regression evaluates the MPMC_HNU canonical residual Helmholtz kernel at
the exact same ten `(T,V,n)` coordinates and applies the already-frozen thresholds
without widening them. The complete ten-state audit reports:

| Observable | Maximum absolute difference | Frozen threshold | Result |
| --- | ---: | ---: | --- |
| `F_res=A_res/(RT)` | `3.2996314014432926e-12` | `1e-10` | PASS |
| pressure | `4.9072827096097171e-05 Pa` | `5e-6 Pa` | **FAIL** |
| `mu_i^res/(RT)` | `2.6485480475457734e-11` | `1e-10` | PASS |
| `ln(phi_i)` | `6.7808336723373941e-10` | `1e-10` | **FAIL** |

This pattern is consistent with a pressure / `Z` derivative-level numerical floor in the
pinned Clapeyron oracle rather than a scalar-potential or residual-chemical-potential
formulation mismatch: the scalar `F_res` and `mu_i^res/(RT)` rows pass, while pressure
fails and the resulting `-ln Z` contribution causes `ln(phi)` to fail. This is a
diagnostic interpretation, not permission to relax the frozen envelope.

Therefore the current readiness matrix is:

| Gate | Status |
| --- | --- |
| A — scalar Helmholtz value | **PASS** |
| B — pressure vs MPMC analytic path | **PASS** |
| C — residual chemical potentials vs MPMC analytic path | **PASS** |
| D — final `ln(phi)` + ThermoPack parity | **PASS** |
| E — named independent implementations | **BLOCKED** |

**Gate F remains NOT READY.** No production pressure or `ln(phi)` source-of-truth switch
is authorized while the pinned Clapeyron leg exceeds the predeclared pressure and
`ln(phi)` thresholds. The failing Gate-E regression is intentionally retained as audit
evidence; the thresholds are not widened to make the gate green.

### 12.4 Gate-E adjudication — raw pinned witness versus corrected diagnostic

This adjudication reviews the frozen Gate-E acceptance wording against the complete
Clapeyron evidence chain now present on PR #115. It changes **no numerical threshold**,
does not replace the raw pinned oracle, and does not authorize a production
source-of-truth switch.

#### Frozen wording controls the evidence qualification

Section 12.2.B requires a reproducible oracle from the **unmodified** pinned revision
`ClapeyronThermo/Clapeyron.jl@229b09452f36c2f812486150df0bb43b197bb4e5`
and requires the ten-state cross-implementation comparison to satisfy the already-frozen
`F_res`, pressure, residual-chemical-potential and `ln(phi)` envelopes.

The contract contains one explicit derivation exception: **if Clapeyron does not expose a
quantity directly**, that quantity may be derived from Clapeyron's own public
residual-Helmholtz API. That exception does not apply to pressure or fugacity for this
pinned revision because those observables are exposed directly. The acceptance wording
therefore does not authorize replacing the pinned implementation's association solution
with a separately maintained nonlinear solve after the raw pinned result fails.

This distinction is also required by the original Gate-E purpose: the gate names
**independent external implementations**. A test-side independent Newton correction can
diagnose an external implementation defect, but it is no longer the numerical output of
the named pinned implementation itself.

#### Evidence now established

The raw pinned revision is exactly formulation-compatible with the frozen MPMC_HNU
profile, and its scalar / composition-derivative evidence is strong:

- exact formulation and parameter read-back is established;
- frozen ten-state `F_res` passes the `1e-10` envelope;
- frozen ten-state `mu_i^res/(RT)` passes the `1e-10` envelope;
- raw pressure and `ln(phi)` fail only the already-frozen derivative-level envelopes.

The subsequent diagnostic chain identifies the failure source without changing the model:

1. ordinary-double finite differences show a derivative floor but do not outperform the
   pinned internal derivative;
2. BigFloat scalar differentiation converges to the pinned internal derivative rather
   than to the frozen target;
3. configurable association `rtol/atol/max_iters` do not affect the problematic state;
4. all ten frozen states enter the same compressed 4-site -> 2-site exact-initializer
   association path;
5. the pinned `X_exact2!` path returns before configurable association convergence
   settings are consulted;
6. an independent BigFloat Newton solve of the same compressed equations to
   `max|R_X| <= 1e-30` removes the derivative discrepancy without changing formulation,
   parameters or thresholds;
7. over the complete ten-state set, the independently stationary diagnostic reaches
   `max |Delta P| = 5.323330668403745e-7 Pa` and
   `max |Delta ln(phi)_Z-only| = 1.3571962033602398e-11`, both inside the original
   frozen envelopes.

The raw pinned regression must therefore be retained as a **known numerical-defect
witness**. The independently stationary ten-state result is accepted as **supplemental
defect-attribution evidence** showing that the model/formulation itself is consistent
with the frozen targets.

#### Adjudication

**No — under the currently frozen Gate-E contract, the independently stationary
same-formulation ten-state diagnostic cannot substitute for the raw pinned Clapeyron
output and cannot by itself make Gate E PASS.**

Accepting it as the contract-satisfying named implementation would retroactively change
the identity of the reference being compared: from the unmodified pinned implementation
to a hybrid consisting of pinned Clapeyron thermodynamic algebra plus an MPMC-maintained
independent association solve. That is a change in evidence semantics even though the
numerical thresholds themselves are unchanged.

Accordingly:

- raw pinned Clapeyron remains a required regression and remains a documented
  numerical-defect witness;
- the independent-stationarity ten-state PASS remains a diagnostic proving that the raw
  failure is not evidence of an MPMC_HNU formulation mismatch;
- Gate E remains **BLOCKED** under the literal frozen acceptance contract;
- Gate F remains **NOT READY**;
- no threshold widening, oracle substitution or production switch is authorized by this
  adjudication.

Closing Gate E now requires a separately reviewed change to the external-reference
contract, not a reinterpretation of the existing frozen one. A scientifically clean
follow-up is to identify an **unmodified external revision/implementation** in which the
compressed exact-association numerical defect is fixed while the audited CPA formulation
and parameter mapping remain unchanged, then declare that reference change explicitly
before rerunning the same frozen ten-state numerical envelopes.

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

## 15. Current audit conclusion and recommended next increment

For the frozen SRK+sCPA profile, Gates A-D are now satisfied. The named Clapeyron leg has
also been investigated to completion at the current pinned revision: formulation and
parameter identity are established, the raw pinned scalar and residual-chemical-potential
rows pass, and the raw pressure / `ln(phi)` failure has been traced to the pinned
compressed `X_exact2!` numerical path. The independently stationary same-formulation
diagnostic passes the original ten-state pressure and Z-only `ln(phi)` envelopes without
changing parameters or thresholds.

That diagnostic is **not** contract-equivalent to the raw pinned implementation under the
wording frozen in section 12.2.B. Gate E therefore remains **BLOCKED** and Gate F remains
**NOT READY**. The raw pinned regression is retained as a known numerical-defect witness;
the corrected diagnostic is retained as defect-attribution evidence, not as a substituted
oracle.

The next small increment should therefore be:

> Audit later unmodified Clapeyron revisions for a fix to the compressed `X_exact2!`
> numerical path while holding the audited CPA formulation and parameter mapping fixed.
> If a source revision with the defect fixed is found, record the exact source change and
> propose a separately reviewed Gate-E reference-repin contract before rerunning the
> existing frozen ten-state thresholds. Do not modify MPMC_HNU production code or Gate-E
> numerical tolerances in that increment.
