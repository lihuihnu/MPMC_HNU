# SW92 Profile-C Gate C2a2 hydrocarbon H -> L/V physical-role resolution audit

## Status, scope, and controlling decision

This is the formal write-before-code Gate C2a2 scientific/architecture audit for

```text
SW92-equilibrium/phase-assigned-aq-na-joint/v1
```

against repository baseline

```text
main@6a46a82451f4c4f2744a22a263ba56c8733f4313
```

The baseline contains:

- Gate C1: a fixed physical `W(AQ)+H(NA)` joint two-phase candidate primitive;
- the post-C1 C2 blocker rejecting unrestricted whole-simplex AQ/NA TPD as Profile-C final
  stability;
- Gate C2a0: the discrete physical-topology contract;
- Gate C2a1: an H-side NA-only phase-addition witness adapter that may discover evidence for a
  second non-aqueous phase, but intentionally publishes no accepted phase set and does not classify
  hydrocarbon phases as liquid or vapor.

This increment is a **scientific/architecture audit only**. It changes no production C++, SW92
thermodynamics, parameter data, tolerance, cubic-root rule, stability implementation, reference
anchor, phase-set publication API, derivative, salt state, or three-phase solver.

### Controlling verdict

**GO for separating hydrocarbon physical morphology from SW92 family/root identity. NO-GO for an
authoritative single-H liquid/vapor classifier based only on Z, cubic-root index, or the
Mortezazadeh pseudo-critical-temperature rule. NO-GO for assuming that a second NA phase discovered
by C2a1 is necessarily vapor.**

The next globally useful numerical increment is not a single-H L/V labeler. It is an **unordered
three-phase candidate primitive**

```text
W(AQ) + H0(NA) + H1(NA)
```

seeded by C2a1 evidence, with both hydrocarbon phases intentionally left physically unclassified
inside the nonlinear equilibrium solve. Physical morphology may be resolved afterwards from
additional evidence. This ordering preserves the project's required `L+L` and `V+L+L` capability
and prevents an engineering phase-label heuristic from changing the thermodynamic equations.

Authoritative Profile-C phase-set publication and SW92 physics coupling remain blocked.

## 1. Evidence policy and source set

This audit distinguishes direct primary/full-text evidence, primary metadata/abstract evidence, and
project inference.

### Tier P1 — primary/full-text sources inspected

1. I. Søreide and C. H. Whitson, *Peng-Robinson predictions for hydrocarbons, CO2, N2, and H2S
   with pure water and NaCl brine*, Fluid Phase Equilibria 77 (1992) 217-240,
   DOI `10.1016/0378-3812(92)85105-H`, plus the project-recorded authors' errata.
2. M. L. Michelsen, *The isothermal flash problem. Part I. Stability*, Fluid Phase Equilibria 9
   (1982) 1-19, DOI `10.1016/0378-3812(82)85001-2`.
3. C. H. Whitson and M. L. Michelsen, *The negative flash*, Fluid Phase Equilibria 53 (1989)
   51-71, DOI `10.1016/0378-3812(89)80072-X`.
4. Z. Li and A. Firoozabadi, *General Strategy for Stability Testing and Phase-Split Calculation in
   Two and Three Phases*, SPE Journal 17 (2012) 1096-1107, DOI `10.2118/129844-PA`.
5. N. Sabet and H. R. Erfani Gahrooei, *A new robust stability algorithm for three phase flash
   calculations in presence of water*, Journal of Natural Gas Science and Engineering 35 (2016)
   382-391, DOI `10.1016/j.jngse.2016.08.068`.
6. E. Mortezazadeh and M. R. Rasaei, *A robust procedure for three-phase equilibrium calculations
   of water-hydrocarbon systems using cubic equations of state*, Fluid Phase Equilibria 450 (2017)
   160-174, DOI `10.1016/j.fluid.2017.07.007`.
7. M. Nazari, M. B. Asadi and S. Zendehboudi, *A new efficient algorithm to determine three-phase
   equilibrium conditions in the presence of aqueous phase: Phase stability and computational
   cost*, Fluid Phase Equilibria (2019), DOI `10.1016/j.fluid.2018.12.013`.

The wider C2a audit's SW92/Xu/Lapene/Ma/Wapperom evidence remains controlling for aqueous versus
non-aqueous model assignment and topology orchestration. C2a2 narrows the question to morphology
inside the non-aqueous family.

### Tier P2 — primary publication metadata/abstract available but full text not in the current evidence bundle

8. Y.-K. Li and L. X. Nghiem, *Phase equilibria of oil, gas and water/brine mixtures from a cubic
   equation of state and Henry's law*, Canadian Journal of Chemical Engineering 64 (1986)
   486-496, DOI `10.1002/cjce.5450640319`.
9. V. Venkatarathnam and L. R. Oellrich, *Identification of the phase of a fluid using partial
   derivatives of pressure, volume, and temperature without reference to saturation properties*,
   Fluid Phase Equilibria 301 (2011), DOI `10.1016/j.fluid.2010.12.001`.

Exact Li-Nghiem branch logic and an implementation of the Venkatarathnam-Oellrich phase
identification parameter are not authorized from abstract-level evidence alone.

### Supporting primary metadata / later description

- Nghiem and Li (1984), *Computation of multiphase equilibrium phenomena with an equation of
  state*, Fluid Phase Equilibria 17, DOI `10.1016/0378-3812(84)80013-8`, provides an important
  multiphase-EOS precedent including liquid-liquid-vapor behavior.
- Later PIP literature describes the Venkatarathnam-Oellrich criterion as a derivative-based phase
  identification method; that description is useful for candidate evaluation but does not substitute
  for auditing the 2011 primary paper before production adoption.

## 2. Four identities that C2a2 must keep separate

The C2a topology audit already separated physical role, SW92 family, cubic root, topology and model
applicability. C2a2 specializes that separation for the hydrocarbon side:

```text
phase instance
    H0 / H1 / ...

thermodynamic family
    SW92 NA

physical hydrocarbon morphology
    liquid / vapor / unresolved

cubic root / molar volume
    local property of that phase state
```

These are not synonyms.

In particular:

```text
SwPhaseFamily::nonaqueous != vapor
SwPhaseFamily::nonaqueous != liquid
low-Z root != universally liquid
high-Z root != universally vapor
```

The public/result architecture must be capable of storing multiple phase instances with the same
physical morphology, e.g. two distinct hydrocarbon liquids.

## 3. Original SW92 does not solve the L/V identity problem

Søreide-Whitson distinguishes aqueous and non-aqueous parameter sets. The non-aqueous set applies
to non-aqueous phase(s), but the 1992 formulation does not provide a general software contract that
classifies every non-aqueous equilibrium state as a unique physical liquid or vapor.

Therefore:

- `AQ/NA` assignment is a thermodynamic-family decision;
- `L/V` is a different physical morphology problem;
- no original-SW92 formula authorizes mapping NA directly to V or L;
- no original-SW92 rule authorizes using an algebraic cubic-root index as the final physical label.

## 4. Michelsen stability — strong basis for NA phase multiplicity, not a direct L/V label

Michelsen's TPD stability framework is naturally applicable when all candidate fluid phases share
one thermodynamic model. This is exactly the relevant structure inside the Profile-C hydrocarbon
side because both prospective `L` and `V` use the SW92 NA family.

Thus Michelsen-style same-model stability is strong scientific support for asking:

```text
Can H(NA) split into another distinct NA phase?
```

That question is already implemented as finite witness evidence in C2a1.

What Michelsen stability does **not** supply by itself is an immutable physical label `liquid` or
`vapor` for one isolated single-phase state. A TPD stationary point, a Gibbs minimum and a cubic
root are thermodynamic/numerical objects; the application-level morphology label remains a separate
contract, particularly near critical continuity.

## 5. Mortezazadeh 2017 — explicit PSD heuristic, useful but not a universal classifier

Mortezazadeh and Rasaei are the most directly relevant SW92-modified-PR precedent. Their algorithm
predefines the only possible phases as

```text
water, liquid hydrocarbon, gas.
```

After water-hydrocarbon stability/split calculations, the calculated hydrocarbon-phase composition
is sent to an oil-gas PST. If the hydrocarbon phase is found single, a phase-state-determination
(PSD) function estimates a mixture critical temperature using their Eq. (37):

```text
Tca = sum_i(Tci * Vci * zi) / sum_i(Vci * zi)
```

and assigns

```text
T >= Tca -> vapor
T <  Tca -> liquid.
```

This is explicit published engineering phase-identification logic, but it is **not** accepted as the
MPMC_HNU default for four reasons.

### 5.1 It belongs to a restricted topology model

The paper assumes at most one oil and one gas phase and ignores additional liquid-liquid
separation. That is narrower than MPMC_HNU's planned `L+L` and `V+L+L` capability.

### 5.2 It requires critical-volume data not present in the current SW92 contract

The current SW92 parameter contract contains `Tc`, `Pc`, acentric factor and family interaction data,
but does not carry the critical volume `Vc` required by Eq. (37). Fabricating `Vc` or silently deriving
it through an unaudited correlation would violate the project provenance rules.

### 5.3 It is a phase-identification heuristic, not a Gibbs/stability theorem

The criterion is not derived as a globally valid morphology discriminator for arbitrary
multicomponent mixtures, liquid-liquid states or near-critical states.

### 5.4 It would control physical labels before the project's broader topology is resolved

Using the PSD label as a solver constraint would prematurely force every H state into exactly one L
or V role.

**Decision:** Mortezazadeh-2017 PSD is rejected as the default authoritative resolver. A future named
compatibility/diagnostic implementation such as

```text
hydrocarbon-role/mortezazadeh2017-pseudocritical-psd/v1
```

is conditionally acceptable only after `Vc` provenance and independent validation are added.

## 6. Mortezazadeh root rules are numerical/engineering choices, not portable morphology identity

Mortezazadeh Appendix B states that, when three real roots exist, the smaller root is selected for
water-hydrocarbon PST and the conjugate root of the base phase is of interest for oil-gas PST.

These rules must not replace MPMC_HNU's current same-family mechanically admissible
minimum-Gibbs-root selector without a separate scientific/numerical audit.

More importantly for C2a2, even a reliable root-selection rule answers:

```text
Which local EOS branch represents this phase evaluation?
```

not:

```text
What physical morphology label must this phase instance carry in every topology?
```

Root branch remains diagnostic/property metadata.

## 7. Sabet 2016 — strongest evidence that phase-type diagnostics and final equilibrium are separate

Sabet and Gahrooei explicitly model physical `L`, `V` and `W` roles. Their algorithm adds an L-V
stability calculation **with water removed**, and the paper explicitly says this auxiliary
calculation is used only for distinguishing the number and type of phases, **not for calculating
final equilibrium properties**.

That is an important architecture precedent:

```text
phase-type evidence != authoritative equilibrium state
```

For MPMC_HNU a similar hydrocarbon-subsystem diagnostic could be useful, but it cannot be adopted
unchanged yet:

1. removing water changes the thermodynamic system;
2. the reduced composition must be renormalized explicitly;
3. non-water/non-water BIPs must remain traceable;
4. the current `Sw92ParameterSet` contract requires exactly one water component and is not a generic
   water-free NA-only model snapshot;
5. any reduced-system classification must remain metadata/initialization evidence and cannot replace
   the full-component W/H material balance.

**Decision:** a Sabet-style water-removed L/V diagnostic is conditionally admissible as a future
named auxiliary operation after a dedicated NA hydrocarbon-subsystem contract. It is rejected as a
substitute for full Profile-C equilibrium.

## 8. Sabet high-CO2 example directly refutes absolute-Z single-H classification

Sabet's high-CO2 example is especially diagnostic. At 15 MPa and 60 degrees C, Sabet's procedure
identifies a two-phase **L-W** state while the Li-Nghiem procedure identifies an approximately
compositionally equivalent **V-W** state. The reported hydrocarbon compressibility factors are
approximately

```text
Sabet:       Z_L = 0.4095
Li-Nghiem:  Z_V = 0.4094
```

while phase fractions and compositions are very similar.

Therefore an absolute rule such as

```text
Z < threshold -> liquid
Z >= threshold -> vapor
```

cannot resolve the physical dispute even in this published case. The authors also point out that an
incorrect L/V label can lead to large errors when downstream properties such as viscosity use
phase-specific correlations.

**Decision:** absolute Z is rejected as an authoritative single-H L/V classifier.

## 9. Nazari 2019 — topology labels pre-exist root selection

Nazari et al. explicitly formulate `L`, `V` and `W` as physical phase roles and enumerate candidate
topologies using negative flash plus TPD. Their EOS implementation chooses the best cubic roots by
Gibbs-energy minimization.

This ordering is decisive for MPMC_HNU:

```text
physical role/topology hypothesis
    -> root/property evaluation
```

not

```text
root index
    -> physical role identity.
```

Nazari also confirms the Sabet high-CO2 `L-W` versus Li-Nghiem `V-W` classification issue and warns
about phase-change/dew-point difficulties in EOS-based aqueous calculations. The work therefore
supports explicit unresolved/boundary semantics rather than forcing a label from a root ordering.

## 10. Whitson-Michelsen negative flash — continuation evidence is not morphology acceptance

Negative flash is useful for continuing algebraic phase-split solutions through phase appearance or
disappearance boundaries. A converged negative-flash branch may be a saddle-like continuation and
is not, by itself, the stable physical phase state.

For C2a2 this means phase fraction sign, continuation branch and L/V label must remain separate.
Negative flash may later help detect topology boundaries but cannot be the sole morphology
classifier.

## 11. Li-Firoozabadi staged stability — phase addition before final topology acceptance

Li-Firoozabadi's staged strategy reinforces that instability evidence creates new candidates and that
a split must subsequently be solved and validated. This matches C2a1's current semantics:

```text
negative NA witness -> seed another NA phase
```

not

```text
negative NA witness -> automatically call it vapor.
```

The role classifier must consume a converged phase state or phase set; it should not control the
thermodynamic equations from an unvalidated witness.

## 12. Nghiem-Li multiphase EOS evidence — a second NA phase need not be vapor

The broader multiphase-EOS literature includes liquid-liquid-vapor equilibrium under a single cubic
EOS. Nghiem and Li's general multiphase work is an explicit precedent for `L1-L2-V` behavior.

This matters directly because Profile C assigns **all hydrocarbon liquid and vapor phases to the
same NA thermodynamic family**. Consequently, after C2a1 discovers a second NA phase,

```text
H0(NA) + H1(NA)
```

may physically represent, depending on the system:

```text
L + V
L1 + L2
```

and near coalescence the morphology may remain unresolved.

Therefore the implication

```text
second NA phase == vapor
```

is scientifically invalid for the project's target state space.

## 13. Project requirement: L+L and V+L+L make a binary L/V-only resolver insufficient

The project planning contract explicitly requires future validation of

```text
V, L, V+L, L+L, V+L+L
```

and separately warns that a water-rich liquid versus hydrocarbon-rich liquid cannot be identified
from array ordering alone.

Thus C2a2 must support **morphology multiplicity**:

- there may be zero or one vapor instance;
- there may be more than one hydrocarbon liquid instance;
- phase-instance identity and physical-role enum must be separate;
- two phase instances may both carry `hydrocarbon_liquid` while retaining different compositions,
  fractions and identifiers.

Any API invariant of the form "there can be at most one phase with role liquid" is rejected.

## 14. What Z/molar-volume ordering can legitimately do

At common `p,T`, molar volume is proportional to compressibility factor:

```text
v = Z R T / p.
```

For two already converged and distinct NA phase instances, lower Z means lower molar volume and
higher molar density; higher Z means the opposite.

That ordering is useful for:

- deterministic representation/canonicalization;
- duplicate detection diagnostics;
- continuity tracking away from coalescence;
- labeling the denser/lighter member **after** an independent contract has established that the pair
  is specifically `L+V`.

It is not sufficient to distinguish `L+V` from `L1+L2`.

Near critical/coalescing states, the volume separation can become unresolved. In that region the
physical morphology result must be `indeterminate/merged`, not chosen by slot order.

**Decision:** relative Z/molar-volume ordering is GO as diagnostic/canonicalization only. It is
REJECTED as the sole LV-vs-LL classifier.

## 15. A single H phase may be intrinsically unclassifiable as liquid or vapor

For a single homogeneous hydrocarbon phase at a specified `p,T,z`, the terms liquid-like and
vapor-like can lose a unique thermodynamic meaning above or near the critical region. A software
interface that always returns either `liquid` or `vapor` therefore risks creating false precision.

Profile C should support an explicit state such as

```text
hydrocarbon_unclassified
```

or equivalent role-resolution status.

This is not a failure of the EOS property evaluation. It is a deliberate physical-role statement:
thermodynamic properties may be valid while a discrete morphology label is not authoritative.

Physics consumers must either support this state or refuse role-dependent constitutive relations.

## 16. Venkatarathnam-Oellrich PIP — promising candidate, not yet authorized

The 2011 Venkatarathnam-Oellrich paper is directly relevant because its stated purpose is to
identify the phase of a fluid from pressure-volume-temperature partial derivatives without relying
on saturation properties, including contexts involving liquid-liquid and vapor-liquid-liquid
behavior.

Later literature describes their phase-identification parameter as a derivative-based discriminator
between vapor and liquid/liquid-like states.

This is substantially more promising for MPMC_HNU than a fixed Z threshold or pseudo-critical
mixing rule because it is built as a thermodynamic phase-identification quantity rather than a
specific petroleum heuristic.

However it is **not authorized for production yet**:

- the project does not currently hold/inspect the complete 2011 primary paper;
- the precise definition, singular/critical behavior and mixture semantics must be verified from the
  primary source;
- implementation requires pressure derivatives with respect to volume and temperature that are not
  part of the current SW92 phase-role contract;
- independent mixture validation is needed for V, L, L+L and L+V states before using it as an
  authoritative resolver.

**Decision:** PIP-based morphology resolution is BLOCKED pending primary-source and derivative
validation.

## 17. Mortezazadeh PSD and PIP must not be conflated

The two approaches answer a superficially similar software question but have different scientific
status:

```text
Mortezazadeh PSD:
    compare T with an estimated mixture critical temperature

Venkatarathnam-Oellrich PIP:
    derivative-based local thermodynamic phase-identification quantity
```

They must have distinct algorithm identities, tests and provenance if both are ever implemented.
Neither may silently become a property of `SwPhaseFamily`.

## 18. Correct sequencing after C2a1

The earlier intuitive sequence was:

```text
C1 W+H
 -> classify H as L or V
 -> if unstable, add the opposite hydrocarbon role
```

This audit rejects that ordering as too restrictive for the project.

The preferred sequence is:

```text
C1 W(AQ)+H(NA)
        |
        v
C2a1 NA-only phase-addition witness
        |
        v
solve W(AQ)+H0(NA)+H1(NA) as an unordered/symmetric thermodynamic candidate
        |
        v
post-solve hydrocarbon morphology evidence
        |
        +-- L + V
        +-- L1 + L2
        +-- unresolved / near-coalescence
```

This prevents an L/V heuristic from controlling which equilibrium equations are solved.

## 19. Consequence for the next three-phase nonlinear primitive

A future `W+H0+H1` solver can be scientifically defined **before** an authoritative L/V classifier,
because both H phases use exactly the same NA thermodynamic family.

The nonlinear thermodynamic requirements are independent of their eventual morphology labels:

```text
one EOS-component inventory
phase fractions sum to one
z_i = beta_W*w_i + beta_0*h0_i + beta_1*h1_i
mu_i^W(AQ) = mu_i^H0(NA) = mu_i^H1(NA)
```

with:

- same p,T and ordered component snapshot;
- same externally fixed NaCl molality for every phase-property evaluation under the current profile;
- same-family minimum-Gibbs mechanically admissible roots;
- no composition floor;
- explicit phase-fraction disappearance;
- pairwise composition distinction;
- resource/failure semantics;
- `global_stability_proven=false` for finite searches.

The two H slots must be treated symmetrically. Slot 0/1 cannot mean liquid/vapor.

## 20. Canonicalization of an unordered H pair

A solver needs deterministic storage/deduplication even when H0/H1 are physically unordered.
Permitted representation rules include:

1. lexicographic composition ordering on a stable component identity; or
2. molar-volume/Z ordering with an explicit near-tie guard.

The selected rule is only **canonical representation**. It must not modify equations and must not
publish physical L/V labels.

Near representation ties, a permutation-equivalent phase set must be deduplicated rather than
reported as two candidate solutions.

## 21. Post-solve physical-role vocabulary

The long-term public role vocabulary should conceptually include:

```text
aqueous
hydrocarbon_liquid
hydrocarbon_vapor
hydrocarbon_unclassified
```

Role values are not unique keys. Multiple phase instances may legitimately have:

```text
phase A: hydrocarbon_liquid
phase B: hydrocarbon_liquid
```

Each phase therefore also requires a separate instance identity/slot/canonical index.

The thermodynamic mapping remains:

```text
aqueous                -> SW92 AQ
hydrocarbon_liquid     -> SW92 NA
hydrocarbon_vapor      -> SW92 NA
hydrocarbon_unclassified -> SW92 NA
```

## 22. Pair morphology needs more information than density ordering

For two converged H phases, a future classifier must distinguish at least:

```text
LV
LL
unresolved
```

A scientifically adequate classifier should use evidence beyond which phase is denser. Candidate
evidence may include:

- a validated derivative-based PIP or equivalent phase-identification quantity for each phase;
- topology continuation/history from a nearby trusted state;
- an explicitly audited auxiliary hydrocarbon-only L/V test;
- near-critical/coalescence diagnostics;
- possibly saturation/phase-envelope information when independently available.

No one of these is adopted by this audit as the default production rule.

## 23. Continuation can help but cannot create identity from history alone

For pressure/temperature paths, continuity from a previously resolved L or V phase can be useful
for tracking phase instances. However:

- continuation can cross critical endpoints;
- liquid-liquid branches can exchange density order;
- phase disappearance/reappearance can break identity;
- a stale historical label cannot override current thermodynamic evidence.

Therefore continuation is a supporting signal and requires an explicit unresolved state near
branch topology changes.

## 24. Single-H publication and physics remain blocked

C2a1 may report no additional NA witness for a C1 `W+H` candidate, but that does not tell whether
H is physically liquid or vapor. It also does not resolve rival W-only/H-only topologies.

Consequently:

- Profile-C `accepted_phase_set()` remains unavailable;
- role-dependent viscosity/density/capillary/relative-permeability models must not consume a forced
  L/V label;
- physics coupling remains blocked until a phase-role contract is authoritative for the states it
  consumes.

A generic hydrocarbon property consumer may use NA thermodynamic properties without an L/V label
only if its constitutive law is explicitly morphology-independent.

## 25. Status semantics required of a future role resolver

A future role-resolution result should distinguish at least:

```text
resolved_liquid
resolved_vapor
resolved_liquid_liquid_pair
resolved_liquid_vapor_pair
near_critical_or_coalescing
insufficient_evidence
out_of_applicability
numerical_indeterminate
```

The exact enum should be designed with the eventual API, but the scientific distinction is
mandatory.

A resolver must retain the evidence actually used so that downstream code can audit whether the
label came from PIP, a named PSD policy, continuation, or another method.

## 26. What may be implemented as named compatibility diagnostics

This audit conditionally permits two non-authoritative diagnostics after their prerequisites are
met.

### 26.1 Mortezazadeh-2017 pseudo-critical PSD

Allowed only as a named compatibility/initialization diagnostic after:

- traceable critical-volume data are added;
- Eq. (37) is independently reproduced;
- V/L and near-critical regression states are added;
- callers are told that the result is heuristic evidence, not authoritative morphology.

### 26.2 Sabet-2016 water-removed L/V diagnostic

Allowed only after a dedicated hydrocarbon-subsystem audit defines:

- component projection and normalization;
- parameter provenance in the projected model;
- how water interactions are removed rather than silently ignored;
- how the auxiliary result maps back to the full-component phase state;
- explicit statement that it does not determine final equilibrium properties.

Neither diagnostic is required before the unordered `W+H0+H1` thermodynamic primitive.

## 27. What is explicitly rejected

This audit rejects the following production semantics:

- `NA == vapor`;
- `NA == liquid`;
- low Z always means liquid and high Z always means vapor;
- selected cubic root index as a physical role enum;
- C2a1 additional-NA witness automatically becomes vapor;
- one and only one hydrocarbon liquid phase globally;
- Mortezazadeh Eq. (37) as an undocumented default;
- deleting water and using the resulting auxiliary flash as final full-system equilibrium;
- forcing a liquid/vapor label in a near-critical or phase-coalescence state;
- letting an inferred L/V role alter the `W+H0+H1` joint equations before that role is validated.

## 28. Revised Gate decisions

| Gate | Decision | Required meaning |
| --- | --- | --- |
| C1 `W(AQ)+H(NA)` joint equations | **PASS / implemented** | One imposed W+H candidate; H morphology unresolved. |
| C2a0 discrete topology contract | **PASS / established** | Physical topology controls orchestration. |
| C2a1 H-side NA phase-addition witness | **PASS / implemented** | Finite H-split evidence only; no new phase publication. |
| C2a2.0 role/family/root/topology separation | **GO — established by this audit** | H morphology is separate metadata/evidence. |
| Single-H L/V from root or absolute Z | **REJECTED** | Published evidence and near-critical semantics do not support it. |
| Mortezazadeh pseudo-critical PSD as default | **REJECTED** | Restricted-topology heuristic; missing `Vc` in current contract. |
| Mortezazadeh PSD as named diagnostic | **CONDITIONAL GO** | Requires sourced `Vc` + independent validation. |
| Sabet water-removed L/V test as final equilibrium | **REJECTED** | Auxiliary type evidence changes the system. |
| Sabet-style hydrocarbon-only diagnostic | **CONDITIONAL GO** | Requires dedicated projected-subsystem contract. |
| Relative Z/molar-volume ordering | **GO as diagnostic/canonicalization only** | Never sufficient for LV-vs-LL classification. |
| Venkatarathnam-Oellrich PIP classifier | **BLOCKED pending full primary paper + validation** | Promising derivative-based role evidence. |
| Authoritative arbitrary-mixture H->L/V resolver | **BLOCKED** | Insufficient validated morphology evidence. |
| `W+H0+H1` unordered joint candidate primitive | **CONDITIONAL GO** | Both H slots are symmetric NA phases; no L/V assumption needed. |
| `W+H0+H1 -> W+L+V / W+L1+L2` publication | **BLOCKED** | Requires post-solve pair morphology resolver. |
| Profile-C authoritative `accepted_phase_set()` | **BLOCKED** | Rival topology and morphology evidence incomplete. |
| SW92 physics coupling | **BLOCKED** | Physical phase roles not authoritative end-to-end. |

## 29. Required design for the next production increment: C2b.1

The next production increment should be

```text
C2b.1 — unordered W(AQ)+H0(NA)+H1(NA) joint three-phase candidate primitive
```

seeded by an admissible C2a1 additional-NA witness.

Minimum contract:

1. one `W` phase fixed to AQ family;
2. two **symmetric/unordered** H phases fixed to NA family;
3. no H slot is called liquid or vapor in the nonlinear solver API;
4. one common EOS-component inventory;
5. two independent phase fractions with explicit simplex/feasibility semantics;
6. common reduced chemical potentials across W, H0 and H1;
7. same fixed NaCl molality in all phase-property evaluations under the current profile;
8. same-family minimum-Gibbs root selection for every phase evaluation;
9. explicit material-balance, chemical-potential, phase-fraction and pair-distinction checks;
10. phase disappearance retained as boundary/unresolved evidence rather than automatic topology
    acceptance;
11. H0/H1 slot-swap equivalence and deterministic canonicalization;
12. no `accepted_phase_set()`;
13. no L/V physical role publication;
14. `global_stability_proven=false`.

The first implementation should be a candidate primitive only. It should not simultaneously solve
water appearance/disappearance orchestration, pair morphology or autonomous phase-number
selection.

## 30. Required validation before C2b.1 merge

At minimum:

1. independent high-precision three-phase equation reference for a state with one AQ and two NA
   phases, using only traceable parameters if a suitable physical source state is available;
2. otherwise an explicitly tagged synthetic structural fixture may validate the generalized
   material-balance solver, but must not be used for a physical-equilibrium claim;
3. all three phase fractions strictly positive in the interior anchor;
4. one common component inventory to numerical tolerance;
5. common reduced chemical potentials across all three phases;
6. pairwise composition distinction;
7. H0/H1 slot-swap invariance/deduplication;
8. W/H role topology guard retained without assigning H0/H1 as L/V;
9. one H-phase disappearance boundary;
10. W disappearance boundary retained as unresolved for later water-topology orchestration;
11. root/property failure and resource exhaustion remain explicit;
12. runtime component permutation;
13. public-header self containment;
14. C1 and C2a1 affected regressions rerun;
15. no existing SW92 thermodynamic formula/tolerance modified merely to obtain a three-phase anchor.

## 31. Evidence still required before an authoritative morphology resolver

Two primary-source gaps are important enough to close before productionizing a general L/V
resolver:

1. Venkatarathnam and Oellrich (2011), DOI `10.1016/j.fluid.2010.12.001` — needed to audit the
   exact derivative-based PIP definition, mixture behavior and singular/critical semantics.
2. Li and Nghiem (1986), DOI `10.1002/cjce.5450640319` — needed before reproducing any exact
   historical liquid/vapor/aqueous phase-detection branch attributed to that algorithm.

The absence of these papers does **not** block C2b.1, because C2b.1 intentionally avoids physical
L/V classification. It does block claiming that a future specific PIP or Li-Nghiem decision rule has
been faithfully implemented from primary evidence.

## 32. Final recommendation

Do **not** implement a single-H L/V classifier next.

Proceed first with C2b.1: solve an unordered, family-explicit `W(AQ)+H0(NA)+H1(NA)` three-phase
candidate seeded by C2a1. Keep both H phases physically unclassified throughout the equilibrium
solve. This closes the thermodynamic phase-addition path without prejudging whether the pair is
`L+V` or `L1+L2`.

In parallel, close the primary-source gap for Venkatarathnam-Oellrich 2011 and Li-Nghiem 1986. A
later post-solve morphology gate can then compare derivative-based PIP, a named hydrocarbon-only
L/V diagnostic, continuation evidence and near-critical/coalescence semantics on converged NA
phases. Only after that gate is independently validated should Profile C publish authoritative
`W+L`, `W+V`, `W+L+V`, or `W+L1+L2` physical phase labels or allow SW92 physics coupling.
