# SW92 Profile-C Gate C2a physical-phase-topology and stability-orchestration audit

## Status, scope, and controlling decision

This is the formal literature-driven Gate C2a audit for

```text
SW92-equilibrium/phase-assigned-aq-na-joint/v1
```

against repository baseline

```text
main@da6298742f8aeaeb68d8132501762d854562b6df
```

The baseline already contains:

- Gate C1, a fixed two-phase physical `aqueous + nonaqueous` joint primitive with one EOS-component
  material balance and cross-phase reduced-chemical-potential equality;
- the post-C1 C2 blocker proving that an unrestricted whole-simplex AQ/NA TPD review at the C1
  common tangent reintroduces Xu-style lower-envelope family competition.

This increment is a **scientific/architecture audit only**. It changes no production C++, EOS
parameter, tolerance, root-selection rule, reference anchor, phase-set publication API, salt state
variable, derivative, or three-phase solver.

### Controlling verdict

**GO for a discrete physical-phase-topology contract. NO-GO for treating a universal composition
partition `Omega_AQ/Omega_NA` as the default scientific definition of Profile C. Authoritative
Profile-C phase-set publication remains BLOCKED.**

The literature reviewed here converges on a stronger and more transferable architecture:

```text
physical phase hypothesis / topology
        -> thermodynamic model assigned to each physical role
        -> role-targeted stability/candidate generation
        -> joint phase-split equations for that topology
        -> topology-specific validation / phase addition or disappearance
        -> publication only after all required rival topologies are resolved
```

For Profile C the physical phase type, the SW92 thermodynamic family, and the cubic-root branch must
remain separate concepts.

A fixed absolute water-fraction cutoff is not part of original SW92. Ma et al. (2021) publish a
`x_H2O >= 0.5` engineering phase-identification policy, and it may be reproduced later under an
explicit named compatibility policy, but it is not adopted here as the Profile-C core definition.

## 1. Evidence policy and source set

The audit distinguishes direct primary-source statements from project inference.

### Tier P1: primary/full-text sources directly inspected

1. I. Søreide and C. H. Whitson, *Peng-Robinson predictions for hydrocarbons, CO2, N2, and H2S
   with pure water and NaCl brine*, Fluid Phase Equilibria 77 (1992) 217-240,
   DOI `10.1016/0378-3812(92)85105-H`, plus authors' errata already recorded by the project.
2. M. L. Michelsen, *The isothermal flash problem. Part I. Stability*, Fluid Phase Equilibria 9
   (1982) 1-19, DOI `10.1016/0378-3812(82)85001-2`.
3. C. H. Whitson and M. L. Michelsen, *The negative flash*, Fluid Phase Equilibria 53 (1989)
   51-71, DOI `10.1016/0378-3812(89)80072-X`.
4. G. Xu, W. D. Haynes and M. A. Stadtherr, *Reliable Phase Stability Analysis for Asymmetric
   Models*, Fluid Phase Equilibria 235 (2005) 152-165,
   DOI `10.1016/j.fluid.2005.06.016`.
5. A. Lapene, D. V. Nichita, G. Debenest and M. Quintard, *Three-phase free-water flash
   calculations using a new Modified Rachford-Rice equation*, Fluid Phase Equilibria 297 (2010)
   121-128, DOI `10.1016/j.fluid.2010.06.018`.
6. Z. Li and A. Firoozabadi, *General Strategy for Stability Testing and Phase-Split Calculation in
   Two and Three Phases*, SPE Journal 17 (2012) 1096-1107, DOI `10.2118/129844-PA`.
7. N. Sabet and H. R. Erfani Gahrooei, *A new robust stability algorithm for three phase flash
   calculations in presence of water*, Journal of Natural Gas Science and Engineering 35 (2016)
   382-391, DOI `10.1016/j.jngse.2016.08.068`.
8. E. Mortezazadeh and M. R. Rasaei, *A robust procedure for three-phase equilibrium calculations
   of water-hydrocarbon systems using cubic equations of state*, Fluid Phase Equilibria 450 (2017)
   160-174, DOI `10.1016/j.fluid.2017.07.007`.
9. M. Nazari, M. B. Asadi and S. Zendehboudi, *A new efficient algorithm to determine three-phase
   equilibrium conditions in the presence of aqueous phase: Phase stability and computational
   cost*, Fluid Phase Equilibria (2019), DOI `10.1016/j.fluid.2018.12.013`.
10. X. Ma, S. Wu, G. Huang and T. Fan, *Three-Phase Equilibrium Calculations of
    Water/Hydrocarbon/Nonhydrocarbon Systems Based on the Equation of State (EOS) in Thermal
    Processes*, ACS Omega 6 (2021) 34406-34415, DOI `10.1021/acsomega.1c04522`.
11. M. Wapperom, J. Heringer, D. V. Nichita and D. Voskov, *A hybrid-EoS approach for multiphase
    pressure-based equilibrium calculations of reservoir mixtures with brine*, Gas Science and
    Engineering 154 (2026) 205966, DOI `10.1016/j.jgsce.2026.205966`.

### Tier P2: primary publication metadata/abstract available, full article not in the current evidence bundle

12. Y.-K. Li and L. X. Nghiem, *Phase equilibria of oil, gas and water/brine mixtures from a cubic
    equation of state and Henry's law*, Canadian Journal of Chemical Engineering 64 (1986)
    486-496, DOI `10.1002/cjce.5450640319`.

Its publisher abstract directly establishes cubic-EOS oil/gas plus Henry-law aqueous modeling and
salinity treatment. Detailed stepwise topology claims in this audit are attributed instead to the
later full-text papers by Sabet/Gahrooei, Mortezazadeh/Rasaei, Nazari et al., and Ma et al., which
explicitly describe, adapt, compare, and in some cases correct the Li-Nghiem orchestration.

### Supporting symmetric-multiphase evidence

The audit also uses the phase-addition logic of Nghiem/Li-style general multiphase calculations and
the modern Li-Firoozabadi staged sequence as architectural precedent. These symmetric-EOS results
are not treated as authority for SW92 AQ/NA family assignment.

## 2. Five concepts that must not be conflated

Profile C needs five distinct axes:

```text
physical phase role
    aqueous / hydrocarbon-liquid / hydrocarbon-vapor

thermodynamic family
    SW92 AQ / SW92 NA

cubic root / local branch
    a property/root diagnostic inside the selected family

phase topology
    which physical roles coexist in the candidate state

model applicability
    p, T, salinity, component and possibly composition range supported by the empirical model/data
```

The source literature repeatedly separates at least the first, fourth, and fifth concepts. The
current MPMC_HNU architecture must keep all five explicit so that a numerical root or a lower Gibbs
model surface cannot silently become a physical phase identity.

## 3. What original Søreide-Whitson 1992 actually fixes

The 1992 paper is decisive about **phase-specific parameter assignment**:

- two BIP sets are introduced, one for the non-aqueous phase(s), one for the aqueous phase;
- Eqs. (10a)/(10b) construct separate `a^NA` and `a^AQ` using the corresponding compositions and
  BIPs;
- the symbol list defines `x_i` as normalized aqueous-phase mole fraction and `y_i` as normalized
  non-aqueous-phase mole fraction;
- the conclusions again assign `k_ij^AQ` and `k_ij^NA` to aqueous and non-aqueous phases,
  respectively;
- the authors warn about the thermodynamic inconsistency introduced by a phase-dependent EOS
  attraction parameter.

The same paper separately correlates hydrocarbon/gas solubility **in the aqueous phase** and water
solubility **in the non-aqueous phase**. This supports a physical water-rich versus non-aqueous
interpretation of the empirical parameter sets.

It does **not** provide:

- an autonomous `p,T,z -> {W,L,V}` topology algorithm;
- a universal `x_H2O` cutoff;
- an admissible-domain/KKT definition for AQ and NA trials;
- a rule that the lower of `g_AQ(x)` and `g_NA(x)` chooses physical phase identity;
- a modern finite-resource stability/publication contract.

Therefore original SW92 supports Profile-C phase assignment, but does not close Gate C2a by itself.

## 4. Michelsen 1982 — TPD is controlling mathematics for a symmetric fluid model, not a phase-role oracle

Michelsen's stability formulation is based on Gibbs tangent-plane distance and is explicitly
developed mainly for EOS calculations using **a single model for all fluid phases**. In that setting,
trial-phase composition changes do not change which thermodynamic model represents the phase.

This is the appropriate mathematical baseline for:

- same-model hydrocarbon liquid/vapor stability inside SW92 NA;
- generic stationary-point and negative-TPD semantics;
- using an instability witness to initialize an added phase.

It does not by itself authorize applying two phase-specific empirical models over the same entire
simplex and letting the lower one replace physical role. That is exactly the extra asymmetric-model
problem treated by Xu, and it is not the Profile-C definition.

## 5. Xu 2005 — valid Profile B, explicit counterexample to Profile-C semantics

Xu/Haynes/Stadtherr introduce a pseudo-TPD for asymmetric models with a binary model-choice
variable. The resulting generalized problem treats the model surfaces as competitors in one
stability problem. This provides a coherent basis for existing Profile B:

```text
SW92-equilibrium/xu-asymmetric-gibbs/v1
```

The post-C1 C2 audit already demonstrated why importing this whole-simplex competition into Profile
C is invalid: at an existing physical NA composition `y`, an unrestricted AQ TPD equals
`g_AQ(y)-g_NA(y)` and can be robustly negative. That event is model reassignment, not discovery of a
new physical phase topology.

**Decision:** Profile B remains unchanged and is not reused as Profile-C final stability.

## 6. Li & Nghiem 1986 — phase-specific models are established, but detailed orchestration remains a primary-source gap

The publisher abstract directly states:

- oil and gas are modeled by a cubic EOS;
- gas solubility in the aqueous phase is modeled by Henry's law;
- pressure/temperature dependent Henry constants and salt corrections are included.

This is an early clear example of `physical phase role -> thermodynamic model` assignment.

The detailed stepwise stability flow is not claimed here from inaccessible primary text. Later
full-text work consistently describes Li-Nghiem as a stepwise algorithm for identifying liquid,
vapor and aqueous phase sets and documents phase-type failures for some high-CO2 cases. Those later
comparisons are sufficient to motivate a better topology contract, but the Li-Nghiem original PDF
should still be added to the project evidence set before reproducing any of its exact decision rules.

## 7. Whitson-Michelsen negative flash 1989 — useful topology continuation, not acceptance by itself

Negative flash allows equal-fugacity/material-balance solutions to continue into regions where the
physical system is single phase. The authors show that the negative-flash solution can correspond to
a Gibbs saddle point.

For Profile C this establishes an important separation:

```text
algebraically converged phase-pair continuation
    != authoritative physical phase-set acceptance
```

Negative flash is therefore suitable for:

- phase appearance/disappearance continuation;
- topology boundary detection;
- initialization or diagnostic phase-state algorithms.

It cannot, by itself, replace stability or physical-role validation.

## 8. Lapene et al. 2010 — explicit physical topology, but under a free-water model

Lapene et al. formulate a three-phase hydrocarbon-water flash using the free-water assumption:
water-rich liquid is pure water and non-water components are absent from that phase. The work is a
strong precedent for explicit water/hydrocarbon topology and negative-flash handling.

It is **not** a direct model for Profile C because Profile C intentionally retains mutual solubility
of hydrocarbons and water. Its value is architectural: phase assignment and topology are explicit,
not recovered from global competition between aqueous and non-aqueous property surfaces.

## 9. Li & Firoozabadi 2012 — staged phase addition is the robust generic pattern

Li/Firoozabadi present a symmetric-model sequence:

```text
single-phase stability
    -> two-phase split when unstable
    -> two-phase stability
    -> three-phase split when the two-phase state is unstable
```

They emphasize multiple stability initializations, distinguish trivial/non-trivial stationary
solutions, and note that an apparent two-phase split may be an incorrect unstable solution even when
a corresponding three-phase solution does not exist; another candidate must then be tried and the
correct solution identified.

This supports three Profile-C design rules:

1. a stability witness is **candidate-generation evidence**, not automatically an accepted new
   topology;
2. a converged split must be validated as a whole phase set;
3. candidate multiplicity and retry/deduplication semantics are part of correctness.

The paper uses one PR-EOS model for the phases and therefore does not supply the missing
phase-specific AQ/NA role rule.

## 10. Sabet & Gahrooei 2016 — topology-first phase identification is explicit

Sabet/Gahrooei use Henry's law for the water-rich phase and an EOS for hydrocarbon liquid/vapor.
Their stability algorithm operates over explicit physical hypotheses. The paper states that a
water-rich phase can disappear at low overall water, high temperature, or low pressure.

Most importantly, the paper describes two-phase checks as flash calculations containing a stability
test, whose result identifies phase A, phase B, or A+B. Its flowchart enumerates `L-W`, `V-W`,
`L-V`, `L-V-W`, and single-phase outcomes.

The paper also introduces an `L-V` stability calculation with water removed **only for determining
phase number/type, not for final equilibrium properties**. That distinction is particularly useful
for MPMC_HNU: a phase-identification auxiliary calculation does not have to be the final
thermodynamic state.

Their high-CO2 example shows that two methods can produce similar phase compositions/fractions but
different `L-W` versus `V-W` labels; the authors stress that phase-type errors propagate to physical
properties such as viscosity.

**Decision for MPMC_HNU:** physical phase role must be authoritative metadata; `SwPhaseFamily::NA`
alone is not enough for physics coupling.

## 11. Mortezazadeh & Rasaei 2017 — strongest direct SW92 topology precedent

This paper directly uses the Søreide-Whitson modified PR model and explicitly describes different
water-hydrocarbon BIPs for aqueous and non-aqueous phases.

Its full algorithm treats all possible physical phases as:

```text
water
liquid hydrocarbon
gas
```

The sequence is topology-first:

1. perform water-hydrocarbon phase stability testing;
2. use a near-pure-water trial (`0.999` water, remaining `0.001` shared by other components) as an
   initialization strategy;
3. if a water/hydrocarbon split is obtained, use the calculated hydrocarbon composition for an
   oil-gas stability test;
4. if that hydrocarbon state is unstable, solve the full water-oil-gas three-phase problem;
5. otherwise retain a water+single-hydrocarbon topology.

The authors' `0.999` water number is explicitly an **initial trial composition**, not a universal
phase-domain boundary.

The paper also uses a temperature-versus-estimated-critical-temperature PSD heuristic to label a
single hydrocarbon phase as gas or oil and recommends small/conjugate cubic roots in particular PST
contexts. These are valuable literature comparison points but are **not adopted** as MPMC_HNU
scientific contracts without separate validation. The existing mechanically-admissible
same-family minimum-Gibbs root selection remains preferable as the numerical/root layer.

The paper reports broad phase-envelope tests and also acknowledges data limitations, including
assuming zero non-water BIPs when source information is unavailable in one case. Such examples may
be algorithmic regressions but must not be promoted to strict physical-validation anchors in this
project without traceable pair data.

## 12. Nazari et al. 2019 — explicit topology graph plus negative-flash/TPD cross-check

Nazari et al. formulate `L`, `V`, and `W` as distinct physical phases with common material balance
and fugacity equality. Their stability orchestration begins from an `L-V-W` flash/negative-flash
hypothesis and, if necessary, checks `V-L`, `V-W`, `L-W`, and single-phase alternatives.

They explicitly state that two-phase subsystems are checked with both negative flash and TPD; when
negative flash fails to converge or disagrees with TPD, the phase set predicted by TPD is used.
They also report that negative flash alone can mispredict near a dew point and that the full
flash/stability logic must resolve the phase set.

This reinforces:

- topology hypotheses are discrete objects;
- negative phase fraction is a topology diagnostic, not sufficient final proof;
- stability and flash evidence must be reconciled;
- aqueous property-model applicability matters near phase boundaries.

Their statement that Henry's law is valid for low dissolved-solute concentration is a
**property-model applicability condition**, not a universal aqueous/non-aqueous phase classifier.

## 13. Ma et al. 2021 — a published 50% policy, useful only as a named engineering option

Ma et al. explicitly identify the missing phase-identification problem when applying the two SW92
BIP sets. Their sequential algorithm uses water vapor pressure to guide whether an aqueous phase may
be present, starts the first stability analysis with a nearly pure-water trial, splits aqueous first
when found, then applies a second stability analysis to the hydrocarbon phase; a third phase is added
if that hydrocarbon state is unstable.

They also **stipulate** the rule:

```text
x_H2O >= 0.5 -> waterlike -> use aqueous BIPs
x_H2O <  0.5 -> use nonaqueous BIPs
```

This is real published evidence and must not be ignored. However:

- the rule is introduced by Ma et al., not by SW92 (1992);
- it is presented in the context of their thermal-process algorithm;
- it is an engineering phase-identification policy, not a general thermodynamic theorem;
- using it as the Profile-C default would silently change the model identity.

**Decision:** reserve it only as a future explicit policy, e.g.

```text
phase-domain-policy/ma2021-waterlike-50pct/v1
```

if independent validation demonstrates value. Do not hard-code it into `Sw92Phase` or the default
Profile-C topology contract.

## 14. Wapperom et al. 2026 — strongest modern evidence for explicit phase-type/model-validity separation

Wapperom et al. deliberately combine a cubic EOS for non-aqueous phases with a separate brine
thermodynamic model. The work describes the models as representing **different phase types** and
notes that the multiple-model construction is thermodynamically inconsistent most visibly near
water critical conditions.

The paper further distinguishes vapor-like and liquid-like cubic roots in water-rich compositions to
avoid interpreting two mathematical minima as two different physical phases when both describe
liquid water. This directly supports the MPMC_HNU rule that root identity is not physical phase
identity.

The aqueous model is intended for a water-rich regime and its validated range matters. This is
strong support for carrying an explicit **model applicability** contract alongside topology.
However, the work does not establish a universal composition threshold that can be transferred to
SW92 Profile C.

**Decision:** model-validity ranges can reject or mark a requested property evaluation as outside
validated applicability, but they must not silently become the universal topology classifier.

## 15. Literature synthesis matrix

| Source | Explicit physical topology | Phase-specific model assignment | Whole-simplex model competition | Phase addition / staged orchestration | Universal water cutoff |
| --- | --- | --- | --- | --- | --- |
| SW92 1992 | AQ vs non-AQ semantics | **Yes: AQ/NA BIPs** | No | Not specified as modern algorithm | **No** |
| Michelsen 1982 | Generic phases | Single model | No multi-model issue | Stability -> split | No |
| Li-Nghiem 1986 | W/L/V from later descriptions; abstract supports AQ vs HC models | **Yes: Henry vs EOS** | No evidence | Later papers describe stepwise route | No evidence |
| Whitson-Michelsen 1989 | V/L continuation | Single model | No | Negative-flash boundary utility | No |
| Xu 2005 | Model-defined V/L asymmetric problem | Different models | **Yes** | Global pseudo-TPD formulation | No |
| Lapene 2010 | W/O/V | Free-water vs HC | No | Explicit 3-phase flash/negative flash | Free-water assumption, not cutoff |
| Li-Firoozabadi 2012 | Generic 1/2/3 phases | Single model | No | **Yes** | No |
| Sabet 2016 | **W/L/V** | Henry W, EOS L/V | No | **Yes, topology graph** | No |
| Mortezazadeh 2017 | **W/L/V** | **SW92 AQ/NA** | No | **Yes, W-H then H L/V** | No; 0.999 is start |
| Nazari 2019 | **W/L/V** | EOS or Henry W variants | No | **Yes, negative-flash topology graph** | No |
| Ma 2021 | **W/HC then HC split** | **SW92-style AQ/NA** | No | **Yes** | **Yes, 0.5, author-stipulated** |
| Wapperom 2026 | **brine/non-AQ, V/L roots** | **Yes** | Avoids duplicate physical phase | Multiphase framework | No transferable universal cutoff |

The strongest common pattern among phase-specific water/hydrocarbon models is not a global
composition partition. It is **discrete phase-role assignment followed by topology-specific
stability/splitting**.

## 16. Revised Profile-C physical-role contract

The long-term public role concept should be at least:

```text
PhysicalPhaseRole::aqueous
PhysicalPhaseRole::hydrocarbon_liquid
PhysicalPhaseRole::hydrocarbon_vapor
```

with SW92 family mapping:

```text
aqueous             -> SwPhaseFamily::aqueous
hydrocarbon_liquid  -> SwPhaseFamily::nonaqueous
hydrocarbon_vapor   -> SwPhaseFamily::nonaqueous
```

This does **not** require changing `SwPhaseFamily` itself. The family remains a thermodynamic
parameterization identity.

The C1 enum value `Sw92PhysicalPhaseRole::nonaqueous` should be treated as a deliberate transitional
role for the unresolved single hydrocarbon phase `H`; it must not be consumed by physics as if it
already meant liquid or vapor.

### Root rule

Cubic root/branch remains separate. Profile C should retain the existing same-family selection of
the mechanically admissible minimum-Gibbs root and its nonsmooth/conditioning diagnostics.
Literature small-root, conjugate-root, or pseudo-critical-temperature heuristics are comparison
algorithms, not default identity.

## 17. Topology state space

For the targeted reservoir water-hydrocarbon problem, the explicit physical topology vocabulary is:

```text
W
L
V
W+L
W+V
L+V
W+L+V
```

where:

- `W` uses AQ;
- `L` and `V` use NA;
- each accepted multi-phase state must have one material balance and common component fugacities /
  reduced chemical potentials across all retained phases.

For staged development, use the coarser intermediate topology:

```text
W + H
```

where `H` means one unresolved hydrocarbon phase using NA. Gate C1 is exactly this `W+H`
coexistence primitive.

No baseline source reviewed here requires or validates two distinct aqueous phases. Therefore
`W+W`, `W+W+H`, etc. are **not** part of the initial Profile-C topology set. They would require a
separate physical-model audit.

## 18. Reinterpretation of Gate C1

C1 remains scientifically useful and unchanged:

```text
W(AQ) + H(NA)
```

It proves that, for an imposed `W+H` topology, one can solve:

- one component inventory;
- one phase-fraction balance;
- cross-role chemical-potential equality;
- same-family root selection;
- positive interior phase fractions;
- distinct compositions;
- relative water-richness `x_H2O^W > x_H2O^H`.

C1 does **not** determine:

- whether `H` is liquid or vapor;
- whether `H` is itself unstable to an additional NA phase;
- whether a `W` phase should appear from an H-only state;
- whether `W+H` is the authoritative topology among all alternatives.

This is the correct boundary for the existing `candidate()` API.

## 19. Why a composition-domain-first C2 is rejected as the default route

The post-C1 blocker originally suggested that a valid C2 might require explicit
`Omega_AQ/Omega_NA` and constrained KKT minimization. The broader literature now changes the
priority.

A composition domain may still be needed as **model applicability** or as a named engineering
policy, but the reviewed water/hydrocarbon algorithms do not derive phase identity primarily from a
universal partition of the composition simplex. Instead they target a physical phase, generate a
trial for that role, solve/validate a phase topology, and then test whether another physical phase
must be added.

Therefore this audit does not authorize an arbitrary default rule such as:

```text
Omega_AQ = {w : w_H2O >= boundary}
Omega_NA = {w : w_H2O < boundary}
```

nor the pair-relative variant as a complete C2 solution.

A future constrained optimizer may still be useful for a **specific named phase-domain policy**.
That would be a numerical feature subordinate to the topology contract, not the definition of
Profile C itself.

## 20. Correct interpretation of stability under Profile C

Under Profile C, stability is a question about **allowed physical topology changes**.

For a current `W+H` C1 state, the primary unresolved instability is not:

```text
Can AQ evaluated at the existing H composition beat NA?
```

That question is Profile B.

The primary topology question is:

```text
Can H split into two distinct non-aqueous phases while W remains?
```

Conceptually:

```text
W + H  ->  W + H0 + H1
```

and after physical hydrocarbon-role resolution:

```text
W + H0 + H1  ->  W + L + V
```

Both new hydrocarbon phases use the NA family. Thus a same-NA-family TPD search against the C1
common tangent is thermodynamically relevant as a **hydrocarbon-side phase-addition witness
search** and does not require comparing AQ and NA at the same existing H composition.

However, a negative TPD witness alone is not an accepted W+L+V state. It must seed a joint
three-phase candidate and that candidate must satisfy material balance, common chemical potentials,
phase distinction, physical L/V role validation and all final gates.

## 21. H-side NA stability: what can be implemented next, and what it may claim

### C2a1 conditional GO: witness adapter only

A small next increment may accept an existing C1 `W+H` candidate and run finite **NA-only**
stability searches against its common tangent.

Required semantics:

1. The retained W and H compositions are added as diagnostic starts; the H start is the trivial
   retained phase.
2. All thermodynamic trial evaluations use the NA family only.
3. A robust negative NA trial is retained as an `additional_nonaqueous_phase_witness`, not as an
   accepted phase.
4. To qualify as a physically usable H-split seed, the witness must be compositionally distinct from
   H and remain water-poorer than the retained W by a roundoff-scaled role guard.
5. A negative trial failing that role guard is diagnostic evidence about NA's mathematical extension,
   not permission to create another physical W-like phase; the outer status remains unresolved for
   authoritative publication.
6. Finite no-negative output may be called only `no_additional_nonaqueous_witness_found`; it is not
   a global proof.
7. `global_stability_proven=false` remains mandatory.

This adapter is useful because it separates a valid hydrocarbon-side instability mechanism from the
rejected whole-simplex AQ/NA family competition.

### Why C2a1 cannot publish `accepted_phase_set()`

An NA-only local search that finds no admissible negative witness does not close:

- boundary/role-domain minima missed by finite starts;
- W appearance/disappearance from other topologies;
- single H liquid/vapor role;
- rival W-only or H-only topologies;
- candidate multiplicity of W+H solutions;
- a full W+L+V solution when a witness is found.

Therefore C2a1 is a **candidate/rejection-evidence primitive**, not C2 final acceptance.

## 22. Water appearance/disappearance must be a separate topology operation

The reviewed literature repeatedly targets a water phase with a near-pure-water trial or an
explicit water-containing phase-pair calculation. Examples include Mortezazadeh's 0.999-water start,
Ma's nearly pure-water first stability analysis, and the water-specific branches of Sabet/Nazari.

For MPMC_HNU this motivates a later operation conceptually like:

```text
H-only state
   -> targeted W candidate generation
   -> solve W+H joint equations
   -> inspect phase fraction / topology feasibility
```

rather than an unrestricted AQ search whose local minimum is allowed to replace the existing H
model.

Water disappearance from an existing W+H branch should likewise be treated through phase-fraction
continuation/boundary semantics and negative-flash-like diagnostics, not through opposite-family
Gibbs reassignment.

No production implementation is authorized by this audit for the H-only starting case because its
physical liquid/vapor role and a complete single-phase publication contract remain unresolved.

## 23. Hydrocarbon liquid/vapor role requires its own contract

The literature offers several phase-identification heuristics:

- Mortezazadeh: temperature relative to an estimated mixture critical temperature;
- Sabet: an auxiliary L-V stability test with water removed, used only for phase type/number;
- cubic-root conventions in various papers;
- Wapperom: distinguish vapor-like and liquid-like EOS roots to avoid duplicate physical water
  minima.

None is adopted automatically.

MPMC_HNU should define a separate hydrocarbon phase-state layer with requirements:

- both L and V use the NA thermodynamic family;
- L/V identity is not `SwPhaseFamily`;
- L/V identity is not merely root index or smaller/larger Z;
- near-critical/indistinguishable states remain explicitly unresolved;
- if an auxiliary reduced/no-water calculation is used, its result is phase-identification evidence
  only and never substitutes for the final full-component material balance.

This audit reserves that work for a later gate before authoritative W+L/W+V/W+L+V publication.

## 24. Candidate topology graph for MPMC_HNU

The literature-consistent long-term orchestration graph is:

```text
                    unknown physical state
                           |
              +------------+------------+
              |                         |
        targeted W evidence       hydrocarbon-only evidence
              |                         |
          solve W+H                solve/classify H
              |                         |
       +------+-------+            +----+----+
       |              |            |         |
   W disappears    W+H branch      L         V
                      |
             NA-only H stability
                      |
             +--------+--------+
             |                 |
       no H-split witness   H-split witness
             |                 |
        W+H candidate       solve W+H0+H1
                               |
                           classify H0/H1
                               |
                             W+L+V
```

Every arrow denotes **candidate generation or topology transition**, not automatic acceptance.
Publication happens only after the required local/final evidence for the selected topology is
complete.

The graph deliberately avoids a top-level `AQ surface vs NA surface` comparison.

## 25. Publication contract implied by C2a

A future authoritative Profile-C result must expose at least:

```text
physical topology identity
physical phase role for every phase
thermodynamic family used for every phase
root/branch diagnostics for every phase
one joint material balance
common reduced chemical potentials/fugacities
phase fractions
phase-role/indistinguishability diagnostics
model-applicability diagnostics
stability/topology checks actually performed
unresolved rival topology diagnostics
resource accounting
global_stability_proven = false   // for finite searches
```

`accepted_phase_set()` must remain unavailable when a required rival topology is not checked or a
phase role remains indeterminate.

A C1 `W+H` candidate is therefore not authoritative merely because its equations converge.

## 26. Relationship to model applicability

Topology and applicability are complementary:

```text
physical role says what phase the model is asked to represent;
applicability says whether the model/data are scientifically supported there.
```

Examples:

- SW92 source p/T/salinity/component-data ranges;
- Wapperom's water-rich brine-model range;
- Henry-law dilute-solute limits in Sabet/Nazari-type models;
- Ma's 50% waterlike classifier as an explicit engineering policy.

An applicability failure should normally make the result `indeterminate/out_of_applicability`, not
silently switch thermodynamic family. Switching AQ <-> NA by lower Gibbs is Profile B.

## 27. Revised Gate decisions

| Gate | Decision | Required meaning |
| --- | --- | --- |
| C1 `W+H` fixed AQ/NA equations | **PASS / implemented** | One imposed W+H coexistence candidate; H physical L/V role unresolved. |
| C2 unrestricted AQ/NA TPD | **REJECTED** | Reintroduces Profile-B family competition. |
| C2a0 physical topology/type contract | **GO — established by this audit** | `W`, `L`, `V`, combinations; role/family/root/applicability separated. |
| C2a1 H-side NA phase-addition witness search | **CONDITIONAL GO** | NA-only finite witness search from C1 common tangent; no publication. |
| C2a2 targeted W appearance/disappearance orchestration | **BLOCKED pending dedicated design/validation** | Near-water trial/negative-flash-style candidate logic; no full-simplex AQ reassignment. |
| C2a3 hydrocarbon L/V physical-role contract | **BLOCKED pending dedicated audit** | Do not inherit pseudo-critical or root-index heuristic without validation. |
| C2 authoritative W+H `accepted_phase_set()` | **BLOCKED** | Rival topology coverage incomplete. |
| W+L+V three-phase joint solver/publication | **BLOCKED until C2a1/C2a3** | Needs two NA phases + one AQ phase, one inventory and role-resolved final review. |
| autonomous all-topology PT solver | **BLOCKED** | Requires W and H single-phase role/topology semantics. |
| SW92 physics coupling | **BLOCKED** | Physical phase roles must be authoritative before property/flow consumption. |

## 28. Required validation for C2a1 before production merge

The first post-audit implementation should remain deliberately small.

Minimum tests:

1. C1 CO2/H2O reference: reconstruct exact C1 common tangent and verify the NA retained H phase is a
   trivial zero-TPD start within numerical allowance.
2. C1 CH4/H2O/brine reference: same check under nonzero molality.
3. Confirm the previously fatal `AQ at retained NA composition` witness is **not evaluated as a C2a1
   rejection mechanism**, because C2a1 is NA-only.
4. At least one traceable/synthetic structural state with an additional distinct NA negative witness,
   retained explicitly as an H-split seed.
5. Role guard: a negative NA witness that is not water-poorer than retained W cannot become an
   admissible H-split seed; status remains unresolved for publication.
6. No-negative finite search reports only `no_additional_nonaqueous_witness_found`, with
   `global_stability_proven=false`.
7. Property/root failure and resource exhaustion remain family-tagged and cannot be converted to
   no-instability.
8. Component permutation.
9. Public-header self containment.
10. Existing C1 and C2-blocker Decimal/reference regressions rerun unchanged.

No new `accepted_phase_set()` is allowed in C2a1.

## 29. What is explicitly not authorized

This audit does not authorize:

- adding `x_H2O=0.5` as the default SW92 phase classifier;
- adding an arbitrary pair-dependent composition cutoff and claiming it is SW92;
- whole-simplex AQ+NA stability for Profile C;
- relabeling Profile-B AQ+AQ as two physical aqueous phases;
- treating `SwPhaseFamily::nonaqueous` as synonymous with vapor or liquid;
- selecting physical L/V solely by Z or root index;
- converting negative flash convergence into accepted equilibrium without stability/topology review;
- treating Mortezazadeh's zero non-water BIP assumption as traceable physical data;
- three-phase or physics publication before physical role identity is carried end-to-end.

## 30. Final recommendation

The next production increment should be **C2a1: a Profile-C hydrocarbon-side NA phase-addition
witness adapter** operating on an already converged C1 `W+H` candidate.

It should reuse the generic finite TPD machinery only with the NA family against the C1 common
tangent, retain explicit physical-role guards and diagnostic starts, and publish witnesses/statuses
only. It must not perform AQ/NA lower-envelope competition and must not expose an authoritative
phase set.

After C2a1, audit and implement the hydrocarbon `H -> L/V` role contract and the corresponding
`W+H0+H1 -> W+L+V` joint candidate path. Only after those topology transitions and water
appearance/disappearance semantics are validated should Profile C return to the question of an
authoritative `accepted_phase_set()`.

## Primary-source provenance gap to close

The only named controlling source in the current C2a literature set for which the project does not
yet hold/inspect a full primary text is Li & Nghiem (1986), DOI `10.1002/cjce.5450640319`.
Publisher metadata/abstract and several later full-text papers are sufficient for the conclusions above,
but any future implementation that reproduces a specific Li-Nghiem decision branch must first audit
the original article itself.
