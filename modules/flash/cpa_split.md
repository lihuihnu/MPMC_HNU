# CPA two-phase PT split baseline

## Scope

This gate connects the validated CPA stability provider to the existing model-neutral two-phase PT split equations and phase-set review. It does not introduce a second Rachford-Rice, fugacity-equilibrium, Gibbs, or common-tangent algorithm.

CPA-specific convention:

`CPA/PT/two-phase/density-side-logK-RR-common-tangent/v1`.

The full route is:

```text
feed CPA stability
    -> robust negative-TPD witness
    -> material-balanced logK/RR seed
    -> two-phase fugacity-equality iteration
    -> feed-vs-pair Gibbs check
    -> final common-tangent CPA stability review
```

Only `two_phase_no_instability_found` is an accepted two-phase result. A converged pair whose final stability is unstable or indeterminate is not published as an accepted phase set.

## Split phase-property policy

CPA stability and CPA split need different root policies.

For stability, every composition is evaluated on the minimum-Gibbs mechanically admissible root because TPD requires the lower allowed phase-property envelope.

For the fixed two-phase iteration, two numerical candidate sides are required. `CpaVleEvaluator` therefore selects among resolved simple roots with `dP/drho>0`:

- `liquid_candidate`: highest molar density admissible root;
- `vapor_candidate`: lowest molar density admissible root.

These names are inherited from the generic `PtPhaseRole` interface and are candidate-side labels only. They are not an authoritative physical morphology classifier. If a composition has only one admissible root, both requested roles may return that root; the generic split must still establish distinct compositions/Z and all equilibrium/acceptance conditions.

A CPA near-multiple/tangent root topology, root-count quota, incomplete bounded root search, no admissible root, or invalid property is propagated through the existing `StabilityPropertyError` contract. It is never converted to a default phase property.

## Inner/outer numerical coordination

The generic split acceptance remains at the existing absolute log-fugacity tolerance (`1e-11` by default). During the first hosted regression of this gate, the standalone CPA PT root default (`1e-10` relative pressure closure) produced an approximately `4e-11` fugacity residual floor at the independent equilibrium anchor. The correct fix is to make the **inner CPA density root more accurate**, not to weaken the outer flash criterion.

`CpaVleEvaluator` therefore uses a split-specific CPA root default:

```text
pressure_relative_tolerance = 1e-12
pressure_absolute_tolerance = 1e-7 Pa
```

The standalone `CpaPtPhase` default remains unchanged. Callers can still pass explicit CPA PT options, but production flash validation must ensure the chosen inner accuracy is compatible with the requested outer fugacity tolerance.

## Attempt-budget fairness

The same hosted diagnosis exposed a generic orchestration issue: one negative-TPD witness is deliberately tried in both material-balanced role assignments. A wrong first assignment could consume the entire global split property budget, preventing the alternate assignment from running.

`PtSplitOptions` therefore gains a bounded `max_evaluations_per_attempt` in addition to the existing global `iteration.max_evaluations`. The underlying RR, fugacity residuals, line search, Gibbs selection and final stability equations are unchanged. The historical repository default remains 20000 provider calls per attempt, while the CPA VLE default sets:

```text
max_evaluations_per_attempt = 8192
```

inside the unchanged 20000 global split budget. Thus an unresolved first density-side assignment cannot starve the second assignment, while all work remains bounded. This option is an orchestration resource control, not a thermodynamic tolerance.

## Reused generic equations and acceptance

`iterate_pt_split(...)` still owns the same equations. For each log-K state it:

1. solves Rachford-Rice for an interior phase fraction;
2. reconstructs liquid/vapor candidate compositions;
3. evaluates the two requested CPA density sides;
4. drives all component fugacity residuals

```text
r_i = ln(x_i) - ln(y_i) + ln(phi_i^L) - ln(phi_i^V)
```

with the existing guarded log-K/Gibbs descent;
5. checks absolute and positive-feed relative material balance;
6. requires non-disappearing, composition/Z-distinct phases.

`solve_pt_vle(...)` continues to own initial finite stability, negative-witness seeding, candidate Gibbs selection, feed-vs-pair Gibbs comparison and final common-tangent stability review using the all-root CPA stability provider. Its only behavior extension in this gate is the optional per-attempt resource cap described above.

As for PR76, finite root search + finite TPD multistart means `global_stability_proven=false`.

## Structural validation

Current tests use an explicitly synthetic non-associating CPA limit, so the PT equations can be independently checked against SRK while still exercising the exact CPA production adapters.

At `P=3 MPa`, `T=180 K`, `z=(0.5,0.5)`, an independent standard-SRK solve gives the structural anchor approximately

```text
x_A(liquid-side) = 0.7246279570576378
y_A(vapor-side)  = 0.4605075426594610
beta_vapor       = 0.8504755589206291
```

and the same feed has a robust negative TPD at `w=(0.78,0.22)`.

The focused suite checks:

- three-root pure SRK limit selects distinct highest-density/lowest-density candidate sides;
- direct fixed-seed CPA split reaches the original `1e-11` fugacity tolerance and the independent equilibrium anchor;
- the complete `stability -> split -> final review` route actually attempts both witness role assignments, prevents the unresolved first assignment from exhausting the global budget, accepts the converged alternate assignment, lowers Gibbs relative to the feed, and passes final common-tangent review;
- runtime component permutation preserves phase fraction and remaps both phase compositions by component identity;
- public-header self containment.

These fixtures are software/algorithm regression only, not CPA experimental validation.

## Boundary

After this gate CPA has an accepted maximum-two-phase PT baseline under the declared finite-search contract. It still does not provide:

- a generic authoritative phase-set publication adapter;
- two-to-three additional-phase orchestration;
- fixed-three-phase CPA phase-side adapter;
- phase-disappearance neighbor orchestration beyond the generic two-phase evidence already retained;
- physical CPA parameter datasets/experimental validation;
- flash sensitivities.

The next gate should first publish the accepted CPA 1/2-phase result through the existing generic phase-set/backend contract, then extend the same numerical density-side semantics to the existing model-neutral `PtThreePhase` primitive and max3 orchestration without inventing liquid/L1/L2 morphology labels.
