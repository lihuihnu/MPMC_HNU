# CPA two-phase PT split baseline

## Scope

This gate connects the validated CPA stability provider to the existing model-neutral two-phase PT split orchestration. It does not introduce a second Rachford-Rice or common-tangent implementation.

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

## Reused generic equations and acceptance

`iterate_pt_split(...)` remains unchanged. For each log-K state it:

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

`solve_pt_vle(...)` also remains unchanged. It owns initial finite stability, negative-witness seeding, candidate Gibbs selection, feed-vs-pair Gibbs comparison and final common-tangent stability review using the all-root CPA stability provider.

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
- direct fixed-seed CPA split converges to the independent equilibrium anchor while satisfying material balance and fugacity equality;
- the complete `stability -> split -> final review` route accepts the two-phase state and lowers Gibbs relative to the original feed;
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

The next gate should first publish the accepted CPA 1/2-phase result through the existing generic phase-set/backend contract, then extend the same density-side semantics to the existing model-neutral `PtThreePhase` primitive and max3 orchestration without inventing liquid/L1/L2 morphology labels.
