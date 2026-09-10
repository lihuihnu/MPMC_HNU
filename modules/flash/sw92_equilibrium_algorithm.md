# SW92 equilibrium-algorithm contract

## Status and scope

This document defines **equilibrium-algorithm identity and publication semantics** for the existing
`SW92/corrected-original/PR76-base/NaCl-molality` thermodynamics profile. It is a contract/audit
increment only: it does not implement a phase split, cross-family selector, three-phase solver,
salt-inventory balance, or SW92 derivatives.

The thermodynamics profile and the equilibrium algorithm are orthogonal identities. Future SW92
equilibrium entry points and results must record both. An unqualified label such as `SW92 flash`
is not sufficient because the literature supports more than one materially different way to use
the AQ/NA parameterizations.

The already implemented `Sw92FamilyStabilityEvaluator` remains a valid building block: one
instance fixes exactly one `SwPhaseFamily`, one NaCl molality, one ordered thermodynamic snapshot,
and one root-option set. This contract does not change its same-family root/Gibbs semantics.

## Evidence hierarchy

### 1. Primary model source: Søreide and Whitson (1992) plus authors' errata

Primary source: I. Søreide and C. H. Whitson, *Peng-Robinson predictions for hydrocarbons, CO2,
N2, and H2S with pure water and NaCl brine*, Fluid Phase Equilibria 77 (1992) 217-240,
DOI `10.1016/0378-3812(92)85105-H`, using the project-supplied PDF whose recorded SHA-256 is
`cb5b1d5034d78d934e887449ce0c89692d43431835d1371da5606d95b2c6bf58` and which contains the
authors' errata.

The paper establishes the following model facts relevant to equilibrium orchestration:

- it deliberately uses two BIP sets for water/brine binaries, one for the non-aqueous phase(s)
  and one for the aqueous phase;
- it explicitly warns that using different EOS `a` values in the two phases can be
  thermodynamically inconsistent, while arguing that the engineering approach should be useful
  away from critical conditions typical of hydrocarbon/water systems;
- Eqs. (10a) and (10b) define separate `a^NA` and `a^AQ` mixture parameters;
- the applications and conclusions treat `k_ij^AQ` and `k_ij^NA` as phase-specific empirical
  parameterizations for mutual-solubility prediction.

This source therefore supports **phase-specific thermodynamic parameterizations**. It does not by
itself define Xu-style global lower-envelope stability, nor does it uniquely prescribe a modern
software orchestration API.

### 2. Current Whitson engineering method

The current whitson+ Water Bot manual states that two EOS models are developed, one for the
aqueous phase and another for the non-aqueous phase, and that the calculation is repeated with the
non-aqueous EOS model after an aqueous-model flash:

`https://manual.whitson.com/methods/water-bot/`

That Water Bot workflow includes additional model tuning and application-specific property
calculations, so its numerical parameter state must not be silently equated with the repository's
strict `SW92/corrected-original` thermodynamics profile. Its value here is the **algorithmic
pattern**: AQ and NA are separate model passes, not a root-order heuristic.

A current public Whitson-associated refresh implementation documents an even more explicit
compatibility scheme:

`https://github.com/mwburgoyne/SW_Framework_Refresh/blob/main/shared/vle_engine/_lib_vle_engine.py`

Its stated solution scheme is two independent full-mixture Rachford-Rice flashes: an AQ-BIP pass
from which the aqueous result is taken, and an NA-BIP pass from which the non-aqueous result is
taken. Its `sw_original` mode retains the original Søreide-Whitson alpha/BIP family for the
relevant species. This is contemporary engineering evidence, not retroactive proof that the 1992
paper defined a single unique full-flash algorithm in exactly this software form.

### 3. Alternative generalized theory: Xu, Haynes and Stadtherr (2005)

G. Xu, W. D. Haynes and M. A. Stadtherr, *Reliable Phase Stability Analysis for Asymmetric
Models*, Fluid Phase Equilibria 235 (2005) 152-165, DOI
`10.1016/j.fluid.2005.06.016`:

`https://academicweb.nd.edu/~markst/xu_fpe2005.pdf`

The paper treats phase stability when different thermodynamic models are used for different phase
types. A later author review makes the core construction explicit: the system Gibbs surface is
the lower envelope of the phase-model surfaces, and all model-specific TPD functions are formed
against the same feed tangent plane. Xu et al. introduce a binary model variable in a pseudo-TPD
to handle the asymmetric model choice.

Reference review:

`https://www3.nd.edu/~markst/stadtherr-ijrs2007.pdf`

This is a sound **general asymmetric-model framework**, but its deterministic/global guarantee
comes from interval analysis. MPMC_HNU's current multistart TPD search is deliberately finite and
must retain `global_stability_proven=false`; it cannot inherit Xu et al.'s global guarantee merely
by adopting the same thermodynamic formulation.

## Required algorithm identities

Two different algorithm profiles are reserved conceptually. These strings are documentation
contracts only in this increment; no public C++ enum is introduced yet.

| Algorithm profile | Intended role | May publish one accepted joint phase set? |
| --- | --- | --- |
| `SW92-equilibrium/whitson-dual-model-observables/v1` | Evidence-aligned engineering compatibility / mutual-solubility observables from separate AQ and NA model runs | **No** |
| `SW92-equilibrium/xu-asymmetric-gibbs/v1` | Generalized joint equilibrium with explicit AQ/NA phase models and common asymmetric stability | **Yes, only after joint balance/equilibrium/stability checks** |

Neither profile changes the thermodynamics identity
`SW92/corrected-original/PR76-base/NaCl-molality`.

## Profile A: `SW92-equilibrium/whitson-dual-model-observables/v1`

### Purpose

This is the primary compatibility route because it follows the phase-specific two-model pattern
visible in the original SW92 model design and in current Whitson engineering practice. It is an
**observable calculation composed of two internally separate equilibrium calculations**, not a
single thermodynamically unified AQ+NA phase split.

### Inputs

A calculation must use:

- one ordered `Sw92Phase` thermodynamic snapshot;
- one finite `p > 0` and `T > 0` within the declared data applicability;
- one finite nonnegative composition/feed under the selected application contract;
- one fixed NaCl molality under the existing `NaCl-molality` profile;
- explicit root/stability/split options for each model pass;
- the algorithm profile ID above.

The AQ and NA passes must use the **same dataset/revision/component order/pure-component data**.
Only the family-specific BIP parameterization changes between the passes.

### Two independent model passes

1. **AQ-model pass**: every SW92 water-pair property evaluation in this entire stability/flash run
   uses `SwPhaseFamily::aqueous`. The pass performs its own phase stability, material balance and
   fugacity-equilibrium checks under that one thermodynamic model.
2. **NA-model pass**: the complete calculation is repeated using
   `SwPhaseFamily::nonaqueous` for every SW92 water-pair property evaluation in that run.
3. Each pass retains its own accepted/unstable/indeterminate status, phase fractions,
   compositions, roots, residuals and diagnostics.
4. An application may then extract the aqueous-oriented observable from the AQ pass and the
   non-aqueous-oriented observable from the NA pass.

No phase fraction from the AQ pass may be combined with a phase composition or phase fraction
from the NA pass to fabricate a common material balance.

### Target-phase identification

The contract deliberately does not introduce a universal `x_water > constant` rule. Mature
software may use such engineering heuristics, but MPMC_HNU will not silently adopt them.

For narrowly defined binary compatibility regressions, ordering two distinct accepted phases by
water mole fraction (water-richer versus water-poorer) is permitted without a fixed cutoff; a
numerically unresolved tie is indeterminate. A reusable multicomponent physical phase-label
contract must be audited separately before production observable extraction is generalized.

### Output semantics

A future result for this profile must preserve both model-run results separately. It may expose
application observables derived from the two runs, but:

- it must **not** return the two extracted phases through `PtPhaseSetResult::accepted_phase_set()`;
- it must **not** claim `z_i = sum(beta_alpha*x_i_alpha)` across one AQ-extracted and one
  NA-extracted phase;
- it must **not** claim one common Gibbs minimum or one common tangent plane across the two runs;
- any ratio formed between NA-pass and AQ-pass target compositions must be named as an
  algorithm-specific **cross-model equilibrium ratio**, not silently reused as an ordinary
  Rachford-Rice `K` belonging to one joint material balance.

Success means that the requested compatibility observables are available from individually
accepted model runs. It is not a successful `mutual_solubility` joint phase-set decision under the
project's full flash contract.

### Relation to the project's mutual-solubility goal

This profile is valuable for reproducing documented SW92/Whitson engineering behavior and for
building traceable regression evidence. By itself it does not satisfy the project's stronger
future requirement that all actually coexisting phases participate in one common component
inventory and one set of equilibrium constraints. It may provide comparison data and numerical
initial guesses for a later joint solver, never acceptance evidence for that solver.

## Profile B: `SW92-equilibrium/xu-asymmetric-gibbs/v1`

### Purpose

This is a generalized research route for a **single joint equilibrium problem** whose possible
phase models include SW92 aqueous and non-aqueous families. It is motivated by Xu et al.'s
asymmetric-model stability theory, not claimed to be the original SW92 algorithm.

### Mandatory scientific preconditions

Before this profile can be implemented as an accepted production equilibrium algorithm, an audit
must establish that:

1. AQ and NA reduced chemical potentials/fugacities use compatible component reference states at
   the same `p,T` and fixed molality;
2. both family Gibbs functions are defined over the composition domain traversed by the chosen
   stability/phase-split algorithm, with empirical-validity limitations documented separately;
3. the fixed-molality state variable is intentionally external and no conserved salt inventory is
   being claimed;
4. the known SW92 thermodynamic-consistency caveat from the 1992 paper is retained in the model
   identity and validation report rather than hidden by the solver.

Common formula shape for `ln(phi)` is necessary but not by itself sufficient evidence for these
preconditions.

### Feed reference and asymmetric stability

For each family, the existing fixed-family evaluator first resolves the mechanically admissible,
same-family minimum-Gibbs cubic root at the feed composition. A generalized asymmetric test then
uses the lower family Gibbs surface at the feed as the reference surface.

If the feed reference is resolved, define the common reduced activity/tangent data for each active
component from that one reference. Every family search must then evaluate its TPD against **the
same** feed tangent plane. Conceptually,

```text
D(w) = min(D_AQ(w), D_NA(w))
```

with both terms based on the same feed `g0` and feed slope/reference activities.

The Xu binary-model-variable pseudo-TPD is the theoretical reference treatment for model-switch
nonsmoothness. A first MPMC_HNU implementation that merely runs two finite family-specific
multistart searches against the same imposed reference remains a finite search and must report
`global_stability_proven=false`.

A feed AQ/NA Gibbs near-tie is a model-envelope nonsmooth point. Until explicit pseudo-TPD or an
otherwise validated tie treatment is implemented, such a state must remain `indeterminate`; no
family is selected by index, water-fraction cutoff or execution order.

### Joint phase split

A future accepted phase set under this profile must solve one common problem. Each actual phase
carries an explicit family identity, while cubic root index remains only a local branch
diagnostic. At minimum the accepted set must satisfy, within audited tolerances:

```text
sum(alpha, beta_alpha) = 1
sum(i, x_i_alpha) = 1
z_i = sum(alpha, beta_alpha*x_i_alpha)
ln(x_i_alpha) + ln(phi_i_alpha) = common_i
```

for every active transferable component and every actually coexisting phase to which that
component is allowed to transfer. The phase set must then pass final asymmetric stability search
against the common tangent/reference activities.

Only after these checks may this profile publish an accepted `PtPhaseSetResult`. A converged set of
phase equations without final asymmetric stability is not an accepted result.

### Failure and resource semantics

AQ and NA property failures remain family-tagged. Resource budgets must not acquire an AQ-first or
NA-first bias: per-family search budgets should be symmetric/explicit and total evaluation counts
reported separately. If no robust negative witness exists but one required family search is
indeterminate, the combined stability result is indeterminate.

The current generic stability status meanings are retained; this profile does not redefine
`no_instability_found` as a global proof.

## Relationship to current MPMC_HNU interfaces

- `Sw92FamilyStabilityEvaluator` is valid for both profiles and remains fixed-family.
- `pt_stability.hpp` does not need to know SW92 formulas or phase-family policy.
- The generic `PtPhaseSetResult` remains model-independent.
- Profile A must not project its two independent model-pass outputs into one accepted generic phase
  set.
- Profile B may use the generic phase-set contract only after joint material balance, common
  chemical-potential equality and asymmetric final stability have all been checked.
- Thermodynamics remains below flash; no equilibrium-algorithm selection is moved into
  `Sw92Phase`.

## Mature-software architecture observations

These are architecture references only; no source code or numerical constants are copied.

- **ThermoPack** exposes liquid, vapor, minimum-Gibbs and unidentified-single-phase flags as
  separate concepts. This supports keeping phase identity separate from cubic root index:
  `https://thermotools.github.io/thermopack/vcurrent/phase_flags.html`.
- **Clapeyron.jl** exposes distinct TP-flash method objects (`RRTPFlash`, `MichelsenTPFlash`,
  `MultiPhaseTPFlash`, etc.) independently of the thermodynamic model. This supports explicit
  equilibrium-algorithm profile identity:
  `https://clapeyronthermo.github.io/Clapeyron.jl/dev/properties/flash/`.
- **DWSIM** separates property packages from `FlashAlgorithm` operations such as stability tests,
  phase identification and PT flash. This supports the existing MPMC_HNU dependency direction:
  `https://dwsim.org/api_help/html/Methods_T_DWSIM_Thermodynamics_PropertyPackages_Auxiliary_FlashAlgorithms_FlashAlgorithm.htm`.
- **NeqSim** separates `SystemSoreideWhitson`, `PhaseSoreideWhitson`, the Whitson-Soreide mixing
  rule and thermodynamic operations. Some current mixing-rule paths classify aqueous behavior by
  a fixed water-fraction threshold; MPMC_HNU explicitly treats that as an implementation heuristic
  to avoid, not a scientific definition to copy:
  `https://github.com/equinor/neqsim`.

## Implementation gates and recommended order

The contract changes the recommended next numerical order.

### Gate 1: one-family symmetric SW92 VLE

Implement and validate a **single-family** SW92 phase-split provider/bridge first. One complete run
must use only AQ or only NA property evaluations, reusing the existing generic VLE/stability
machinery where scientifically compatible. This is the missing numerical primitive required by
Profile A.

Acceptance evidence should include fixed-state binary CO2/H2O and CH4/H2O cases with traceable
source/model values, root/phase disappearance behavior, material balance, fugacity equality and
final same-family stability. Existing PR76 tolerances and references must not be silently changed.

### Gate 2: dual-model compatibility orchestration

After both family runs are individually reliable, implement
`SW92-equilibrium/whitson-dual-model-observables/v1` as an orchestration/result layer that retains
the two model runs independently. Reproduce selected original-1992 and current-Whitson
`sw_original` observables without presenting the cross-model output as one material-balanced phase
set.

### Gate 3: asymmetric-Gibbs research route, only if required

If the project requires a thermodynamically joint SW92 mutual-solubility phase set, first perform
the reference-state/domain audit above and add independent cross-family feed-Gibbs and
common-reference TPD anchors. Then implement the Xu-style generalized profile under its distinct
algorithm ID. Do not call it `SW92 original`.

### Gate 4: three-phase

Only after a joint two-family equilibrium formulation has passed independent binary/two-phase
validation should it be generalized to three-phase phase-set search/split. Two Profile-A model
runs must never be stitched together and labeled a three-phase equilibrium.

## Verification policy for this contract

This document is a scientific/architecture contract, not a numerical implementation. The present
increment therefore requires source/provenance review, internal-link review and final diff review,
not C++/CTest/Decimal execution. Numerical CI becomes mandatory when any evaluator, split solver,
status/result type, tolerance, build target or executable reference calculation is changed.
