# SW92 equilibrium algorithm-profile separation audit

## Status, scope, and controlling decision

This is the formal post-Gate-3B audit of equilibrium-algorithm identity for the thermodynamic
profile

```text
SW92/corrected-original/PR76-base/NaCl-molality
```

against repository baseline

```text
main@85869e9f7a3bb79e9d4966b5bf9504661d06c8b9
```

It is a scientific/architecture audit only. It does **not** change SW92 thermodynamics, flash
production code, numerical tolerances, reference anchors, salt state variables, derivatives, or
three-phase capability.

**Audit verdict:** the earlier two-profile taxonomy is scientifically incomplete for the project's
physical mutual-solubility goal. Three materially different algorithm semantics must be kept
separate:

1. `SW92-equilibrium/whitson-dual-model-observables/v1` — two independent complete family flashes
   used to publish compatibility observables; no joint phase set.
2. `SW92-equilibrium/xu-asymmetric-gibbs/v1` — a generalized joint lower-envelope model in which
   AQ/NA are competing thermodynamic surfaces and the lower surface may be selected at each
   composition.
3. **Reserved by this audit:**
   `SW92-equilibrium/phase-assigned-aq-na-joint/v1` — a constrained joint mutual-solubility model in
   which the physical aqueous phase uses the AQ parameterization and the physical non-aqueous phase
   uses the NA parameterization, while both phases obey one component inventory and one set of
   cross-phase equilibrium constraints.

The third identity is intentionally **not** named `original-SW92-flash`: the 1992 paper establishes
phase-assigned AQ/NA parameter semantics and mutual-solubility calculations, but does not uniquely
specify a modern autonomous PT phase-number/stability orchestration API.

Where the older `sw92_equilibrium_algorithm.md` says that only two profiles are reserved, this
later audit extends that taxonomy. Existing Profile-A and Profile-B identities and production
semantics remain unchanged.

## Sources reviewed

### Primary Søreide-Whitson source

I. Søreide and C. H. Whitson, *Peng-Robinson predictions for hydrocarbons, CO2, N2, and H2S with
pure water and NaCl brine*, Fluid Phase Equilibria 77 (1992) 217-240,
DOI `10.1016/0378-3812(92)85105-H`, with the project-recorded authors' errata PDF,
SHA-256 `cb5b1d5034d78d934e887449ce0c89692d43431835d1371da5606d95b2c6bf58`.

The source is decisive about the *parameter role*:

- the authors introduce two sets of BIPs, one for the non-aqueous phase(s) and one for the aqueous
  phase;
- Eqs. (10a) and (10b) define separate `a^NA` and `a^AQ` mixture parameters using the
  corresponding phase compositions and BIPs;
- the symbol list defines `x_i` as normalized aqueous-phase mole fraction and `y_i` as normalized
  non-aqueous-phase mole fraction;
- the conclusions again describe `k_ij^AQ` for the aqueous phase and `k_ij^NA` for the
  non-aqueous phase;
- the paper explicitly warns that this phase-dependent `a` construction may be thermodynamically
  inconsistent in the ordinary symmetric-EOS sense.

Therefore the original SW92 evidence does **not** support silently treating AQ and NA merely as two
interchangeable labels whose assignment is always selected by whichever whole-composition Gibbs
surface is lower. The parameterizations were developed and reported as physical phase-specific
empirical models.

At the same time, the paper does not define the complete modern software contract needed for an
autonomous `p,T,z -> phase set` solver: global stability search, finite-resource failure semantics,
root-envelope ties, phase disappearance, single-phase role identification, and multicomponent
phase-number orchestration require explicit project contracts.

### Current Whitson engineering evidence

The current whitson+ Water Bot procedure develops one aqueous EOS model and one non-aqueous EOS
model, performs an aqueous-model flash, then repeats the calculation with the non-aqueous model.
The current Søreide-Whitson framework refresh likewise describes a **dual-flash scheme** in which
AQ BIPs control dissolved-gas content and NA BIPs control gas-phase water content.

This supports the existing Profile-A compatibility route. It does not convert the outputs of the
two independent flashes into one material-balanced joint phase set.

### Xu, Haynes, and Stadtherr asymmetric model

G. Xu, W. D. Haynes and M. A. Stadtherr, *Reliable Phase Stability Analysis for Asymmetric Models*,
Fluid Phase Equilibria 235 (2005) 152-165, DOI `10.1016/j.fluid.2005.06.016`.

Xu et al. explicitly define, at a tested composition, the system Gibbs surface by whichever of the
model-specific surfaces is lower:

```text
g(x) = min[g^V(x), g^L(x)]
```

and form both model-specific tangent-plane-distance functions against the same reference tangent.
A binary variable represents model choice in the pseudo-TPD formulation. Their deterministic global
guarantee comes from interval analysis, not from the lower-envelope definition alone.

The current MPMC_HNU Profile B is therefore a coherent Xu-style generalized construction provided
its already-audited common reference gauge is retained. It is, however, a **different thermodynamic
construction from the phase-assigned empirical semantics of SW92 (1992)**.

## 1. Separation matrix

| Property | Profile A: dual-model observables | Profile B: Xu lower envelope | Profile C: phase-assigned joint |
| --- | --- | --- | --- |
| AQ and NA evaluations | Two complete independent runs | Competing surfaces in one generalized model | Fixed to physical phase role inside one joint problem |
| One component inventory | **No** across extracted AQ/NA observables | **Yes** | **Yes** |
| Common chemical potentials across retained phases | **No** across the two model runs | **Yes** | **Yes** |
| Model/family assignment | Fixed per complete pass | Chosen by lower Gibbs surface | Constrained by physical aqueous/non-aqueous role |
| AQ+AQ / NA+NA two-phase set | Not a joint-set concept | Allowed if lower-envelope solution selects it | Not an accepted two-role AQ+NA state |
| Opposite-family Gibbs at same phase composition | Diagnostic only | Acceptance-critical family-dominance test | Must **not** reassign the phase; only applicability/diagnostic evidence |
| Binary water-rich/water-poor extraction | Authorized for observables | Not family identity | Candidate physical-role ordering for two-phase primitive |
| Single-phase publication | Per independent family run only | Lower-envelope family may publish under Profile-B rules | **Blocked without a separate physical-role policy** |
| `global_stability_proven` | false for finite searches | false for finite searches | must remain false for finite searches |
| Original-SW92 phase semantics | Compatibility-oriented approximation | No; generalized Xu construction | Closest intended phase-role semantics, but not claimed as a fully specified 1992 flash algorithm |

The three profiles are not interchangeable modes of one numerical solver. The identity controls what
constitutes a valid thermodynamic state and what may be published.

## 2. Profile A audit — PASS, retain unchanged

Existing implementation:

```text
SW92-equilibrium/whitson-dual-model-observables/v1
```

`solve_sw92_whitson_dual_model_observables` correctly constructs AQ and NA evaluators from the same
ordered snapshot and runs two complete `solve_sw92_pt_family_vle` calculations. Both child results
are retained. The outer result exposes no joint `PtPhaseSetResult`, and the documented cross-model
ratio is explicitly not a Rachford-Rice `K` for a common material balance.

Binary target extraction uses relative water-rich/water-poor ordering without a hard composition
threshold. More-than-binary calls retain both complete runs but refuse automatic physical target
labeling.

**Decision:** no Profile-A production change is justified by this audit. Profile A remains the
engineering compatibility/observable route and may be used as comparison evidence or to generate
initial guesses for a later constrained joint primitive.

## 3. Profile B audit — PASS as the declared generalized model, retain identity

Existing parent identity:

```text
SW92-equilibrium/xu-asymmetric-gibbs/v1
```

Existing max2 convention:

```text
SW92-equilibrium/xu-asymmetric-gibbs/max2-logK-SSI-RR/common-tangent/v1
```

The current implementation is internally consistent with its declared Xu-style model:

- the AQ/NA common reduced-chemical-potential gauge was separately audited;
- feed reference uses the lower family Gibbs surface;
- both family stability searches use the same tangent/reference;
- fixed-pair candidate sets permit `AQ+AQ`, `AQ+NA`, and `NA+NA`;
- a converged phase must pass an acceptance-critical lower-envelope family-dominance check;
- final stability is two-family/common-tangent and finite;
- family identity is retained in the result and `global_stability_proven=false` remains explicit.

The latest traceable CH4/CO2/H2O regression is especially diagnostic: both accepted phases are
assigned AQ because `g_AQ-g_NA < 0` at both converged compositions. That is valid evidence for
Profile B, but it cannot be relabeled as evidence that both physical phases are aqueous under the
original SW92 phase-specific interpretation.

**Decision:** do not alter Profile B to force one AQ plus one NA phase. Doing so would no longer be
the implemented Xu lower-envelope model and would invalidate its existing algorithm identity and
regressions. The physical phase-assigned goal belongs in Profile C.

## 4. Profile C — reserve a separate physical phase-assigned joint identity

Reserved identity:

```text
SW92-equilibrium/phase-assigned-aq-na-joint/v1
```

This profile is intended to answer a different question:

> Given one physically aqueous phase and one physically non-aqueous phase at common `p,T`, can the
> SW92 AQ and NA parameterizations be used in one constrained material-balance/equilibrium problem
> to predict their mutual solubility without pretending that two independent flashes are one state
> and without replacing phase role by Xu lower-envelope family selection?

For a two-phase candidate let `x` denote the phase assigned physical role `aqueous`, `y` the phase
assigned physical role `nonaqueous`, and `beta` the non-aqueous mole fraction. The required
continuous equations are

```text
sum_i x_i = 1
sum_i y_i = 1
z_i = (1-beta) x_i + beta y_i
ln(x_i) + ln(phi_i^AQ(x)) = ln(y_i) + ln(phi_i^NA(y))
```

for every active transferable EOS component, at the same pressure, temperature, ordered component
snapshot, and externally fixed NaCl molality.

Within each assigned family, a phase property evaluation should continue to use the mechanically
admissible same-family minimum-Gibbs cubic root. Root index or Z remains a local numerical/property
diagnostic, not the AQ/NA selector.

### 4.1 What Profile C must *not* inherit from Profile B

Profile C must not use

```text
g_system(x) = min[g_AQ(x), g_NA(x)]
```

as a rule for reassigning an already identified physical phase. In particular, the existing
Profile-B `family_assignment_dominated` check is **not** an acceptance criterion for Profile C.
If the NA Gibbs surface is numerically lower at the composition of the physically aqueous phase,
that is evidence about the relationship between the two empirical models; it does not authorize
changing the phase's physical role to NA inside Profile C.

Likewise, Profile C may not publish `AQ+AQ` or `NA+NA` as the result of its first two-role primitive.
Those are valid Profile-B lower-envelope possibilities, not a solution of the constrained
AQ/non-AQ mutual-solubility question.

### 4.2 What Profile C may reuse

The following numerical ideas are reusable after role-neutral wrapping:

- Rachford-Rice as the private fixed-`K` material-balance engine;
- active-support/no-composition-floor contracts;
- log-composition-ratio SSI/Newton directions where independently verified;
- same-family minimum-Gibbs root evaluation;
- absolute/relative material-balance diagnostics;
- common reduced-chemical-potential residuals;
- phase-disappearance and resource-limit semantics;
- finite TPD search machinery against an externally supplied common tangent.

It must expose new Profile-C metadata and physical-role-aware results. Reusing implementation
utilities is not permission to inherit Profile-B model-selection semantics.

## 5. Physical phase-role contract — the main blocker

`SwPhaseFamily` currently identifies which thermodynamic family/model is evaluated. Under Profile B
that is exactly the discrete model variable. Profile C needs a **separate physical phase-role
concept** so that thermodynamic family identity is not silently overloaded.

Conceptually:

```text
physical_role = aqueous       -> evaluate SwPhaseFamily::aqueous
physical_role = nonaqueous    -> evaluate SwPhaseFamily::nonaqueous
```

The first implementation should keep this mapping inside the Profile-C result/solver layer rather
than changing the meaning of `SwPhaseFamily` globally.

### 5.1 Two-phase relative water ordering — conditional GO

For a distinct two-phase state with water active in the feed, the first constrained primitive may
use the relative ordering

```text
x_water(aqueous-role phase) > x_water(nonaqueous-role phase)
```

as its phase-role topology check, with a roundoff-scaled unresolved tie reported as
`phase_role_indeterminate`.

This uses **no absolute water-fraction threshold** and is invariant to component ordering and phase
slot ordering. It is appropriate as a two-phase mutual-solubility topology check, not as a global
statement of empirical model validity.

A result must still document the SW92 data applicability. Two numerically distinct phases that are
both outside the scientifically supported water-rich/water-poor regimes are not made experimentally
valid merely by relative ordering.

### 5.2 Single-phase role — BLOCKED

A single phase has no second phase against which relative water richness can be defined. Selecting
AQ or NA by lower Gibbs would import Profile-B semantics. Selecting by a fixed water threshold would
introduce an unaudited heuristic. Selecting by cubic-root order would confuse phase role with density.

Therefore Profile C must **not initially be an autonomous all-topology PT phase-set solver**.
Single-phase publication under Profile C remains blocked until a separate physical-role policy is
scientifically defined and validated. The first implementation should be a constrained two-phase
AQ+NA joint-split primitive whose success only means that such a two-role state was solved and
validated.

## 6. Stability semantics for the constrained two-phase primitive

Once a candidate AQ+NA pair satisfies one material balance and common reduced chemical potentials,
its common tangent `d_i` is well defined under the already audited compatible reference gauge.
A finite final review may therefore search:

1. AQ trial compositions evaluated with the AQ family against that common tangent;
2. NA trial compositions evaluated with the NA family against the same common tangent.

A robust negative TPD witness means the candidate is not stable within the constrained two-family
state space searched. If no negative witness exists but either required family search is
indeterminate, the candidate remains indeterminate. If both finite searches report no instability,
that is finite evidence only and must retain

```text
global_stability_proven = false
```

Crucially, final stability asks whether an additional allowed AQ or NA phase lowers the constrained
state Gibbs relative to the candidate common tangent. It does **not** compare AQ and NA surfaces at
an existing phase composition and reassign that phase by lower-envelope dominance.

## 7. Initialization and publication boundary

Profile-A outputs can provide candidate aqueous-oriented and non-aqueous-oriented compositions for
initialization, but the two independent Profile-A states are never themselves acceptance evidence
for Profile C.

The first Profile-C numerical increment must independently solve the joint equations and re-check:

- one component inventory;
- one phase-fraction sum;
- common reduced chemical potentials/fugacities;
- physical role ordering;
- same-family root smoothness/conditioning;
- phase disappearance;
- final constrained two-family stability;
- resource/failure status.

Only then may it expose an accepted **Profile-C two-phase candidate**. The result should retain both
`physical_role` and evaluated `SwPhaseFamily` so future audits can detect accidental remapping.

## 8. Validation gates before any production use

### C0 — contract/result metadata: GO

Reserve the identity and define a family/physical-role separation in documentation. No production
code is required in this audit.

### C1 — fixed AQ+NA joint two-phase numerical primitive: CONDITIONAL GO

Implement only after this audit is accepted. Required first evidence:

- direct independent Decimal(80) solution of the cross-family AQ+NA equations for traceable binary
  CO2/H2O and CH4/H2O cases;
- material balance and common reduced chemical potentials from the independent oracle;
- explicit water-rich/water-poor ordering and tie case;
- component/phase-slot permutation invariance;
- root-envelope, disappearance and resource-limit cases;
- no Profile-B family-dominance/reassignment logic.

Do not derive the reference by stitching the two Profile-A flashes.

### C2 — final constrained stability: CONDITIONAL GO after C1

Use a common candidate tangent and independent AQ/NA finite searches. Add negative-witness,
indeterminate-family, boundary, and final-publication regressions. Retain
`global_stability_proven=false`.

### C3 — autonomous phase-number/single-phase Profile-C solver: BLOCKED

Requires a separate phase-role/topology audit. No current evidence authorizes an automatic
single-phase AQ/NA classification rule without changing the algorithm semantics.

### C4 — physical three-component regression: BLOCKED until C1/C2

The existing traceable CH4/CO2/H2O parameter snapshot may be reused as **input data**, but a new
independent AQ+NA joint oracle is required. The current Profile-B AQ+AQ Decimal anchors are not
Profile-C references.

### C5 — three-phase and physics coupling: BLOCKED

A physically phase-assigned three-phase topology is not just a larger max2 problem. A likely
water/hydrocarbon topology can contain one aqueous phase and multiple non-aqueous phases, for
example `AQ+NA+NA`, requiring explicit physical-role multiplicity and phase-number stability.
It must not be inferred by extending the two-role max2 result mechanically.

SW92 physics coupling remains blocked because downstream conservation/closure must know which
equilibrium profile generated the phase set and what its physical roles mean.

## 9. Salinity boundary remains unchanged

All three profiles currently operate conditional on one externally fixed NaCl molality. NaCl is
not an EOS component. None of these profiles may claim closed salt inventory conservation, salt
redistribution between phases, or automatic molality change as water transfers.

A conserved-salt formulation is a separate state-model audit and must receive a different identity.

## 10. Required public metadata after Profile C exists

A downstream phase result must never be described only as `SW92`.

At minimum preserve:

```text
thermodynamics_profile = SW92/corrected-original/PR76-base/NaCl-molality
equilibrium_profile    = SW92-equilibrium/phase-assigned-aq-na-joint/v1
physical_role           = aqueous | nonaqueous
thermodynamic_family    = aqueous | nonaqueous
component_snapshot      = ordered IDs + dataset/revision
fixed_NaCl_molality
stability_search_scope  = finite/non-global
```

This is particularly important before physics coupling: density/root branch, thermodynamic family,
physical phase role, and equilibrium algorithm identity are separate concepts.

## 11. Formal audit conclusion

The repository should **not** choose between the existing Profile A and Profile B by renaming one of
them as the physical SW92 answer.

- Profile A is correct for documented dual-flash compatibility observables and must remain
  non-joint.
- Profile B is correct as the explicitly declared Xu lower-envelope generalized equilibrium and
  must continue to permit lower-envelope same-family solutions such as the validated AQ+AQ cases.
- The project's desired physical mutual-solubility joint calculation requires a third, explicit
  phase-assigned profile with one material balance and cross-family fugacity equality.

The first safe next numerical step is therefore **C1: a fixed AQ+NA joint two-phase primitive with
independent binary Decimal(80) references**, not three-phase, not SW92 physics coupling, and not a
modification of the existing Xu max2 acceptance logic.

## Verification policy for this audit

This is documentation-only scientific scope. Acceptance requires source/provenance review,
algorithm-identity consistency review, internal-link review, and final repository diff review.
C++/CTest/Decimal CI is required only when a later increment changes executable code, tests,
tolerances, or reference calculations.
