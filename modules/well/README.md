# Well module

The well module owns well-model physics independently of the reservoir spatial
discretization and solver backend.

Current layers:

- `mpmc::well`: model-specific cell-connection geometry/property contracts.
  The current baseline is the axis-aligned Cartesian/diagonal-K Peaceman well
  index.
- `mpmc::well_discretization`: coupling of a frozen well connection to
  reservoir-side constitutive quantities. The current contract multiplies
  Peaceman `WI [m3]` by the existing local phase mobility
  `lambda [1/(Pa s)]` and propagates its natural-variable derivative.

The hydraulic-conductance layer defines

`C_alpha = WI * lambda_alpha [m3/(Pa s)]`.

The next connection-local layer defines pressure drawdown and one-phase
reservoir-volume rate with production-positive sign:

`Delta p_alpha = p_alpha,cell - p_bhp`

`q_alpha = C_alpha * Delta p_alpha [m3/s]`.

It carries the full reservoir natural-variable product-rule derivative and the
explicit derivative `dq_alpha/dp_bhp = -C_alpha`. In this layer `p_bhp` is
still a supplied scalar, not a globally numbered unknown.

The connection component layer now converts the three phase rates to canonical
component molar rates:

`n_dot_i = sum_alpha q_alpha c_alpha x_alpha,i [mol/s]`.

It consumes the existing phase molar-density primal/Jacobian and reconstructs
composition derivatives from the frozen natural-variable composition pivot. It
also propagates the explicit BHP derivative and verifies component-to-total
molar closure for primal values and derivatives.

The connection energy layer now defines directional advective energy rate:

`E_dot_alpha = q_alpha rho_alpha h_selected`, summed over phases.

For production and exact-zero ties, `h_selected` is reservoir phase enthalpy.
For injection, `h_selected` is explicit well-side injection enthalpy supplied
by the caller; reservoir enthalpy is deliberately not reused for injected
fluid. Reservoir mass density remains the volume-to-mass conversion because
the phase rate is reservoir-volume based. The reservoir and BHP derivatives
follow the same directional enthalpy selection.

The connection-to-cell-source adapter now bridges the well sign convention to
the model-neutral finite-volume source convention. Well component/energy rates
are production-positive (cell -> well), while `CellSourceLinearization3D` is
injection-positive (into cell), so the adapter applies exactly one sign
reversal to rates and reservoir Jacobians. Explicit BHP derivatives are kept
as a sidecar and are not inserted into the reservoir natural-variable block.

The adapter is intentionally compatible with the existing
`normalize_cell_source_by_bulk_volume()` contract: after its second sign
change, a production-positive well rate becomes a positive outward
conservation residual. This prevents a hidden double-negation convention.

The PETSc bridge lives under `well/discretization/petsc`. It binds one
explicit Peaceman connection with frozen BHP to the existing owner-only
mixed-cardinality cell-source callback. The target may now be a frozen 1P, 2P
or 3P natural-variable cell. The bridge consumes exactly the active phases in
that chart; it does not pad inactive phases with fictitious mobility, density,
composition or enthalpy.

The cardinality-neutral well-discretization layer evaluates

`q_alpha = WI * lambda_alpha * (p_alpha - p_bhp)`

over the current active phase set, then forms component molar and directional
advective-energy rates with the same production-positive convention used by
the original 3P chain. The well-side injection-enthalpy payload must contain
exactly one value per active phase. Reservoir natural-variable derivatives and
explicit BHP derivatives are retained for every supported cardinality.

The bridge reuses the validated `CellSourceLinearization3D` sign/normalization
contract, so component/energy well terms enter the existing
bulk-volume-normalized fully implicit residual and diagonal reservoir Jacobian.
BHP remains a frozen parameter and is not globally numbered. The former
three-phase evaluator names remain compatibility entry points, while the
production path uses the variable-cardinality evaluator.

Phase-transition-aware rebinding is now explicit rather than slot based. A
rebindable fixed-BHP context retains a stable physical-phase ->
well-side-injection-enthalpy registry. After an accepted outer phase-set
rebuild, orchestration locates the same completion by stable cell
`GlobalEntityId` and resolves the rebuilt `FrozenActivePhaseIdentityMap`
back to the new active slots. Peaceman connection data, frozen BHP and source
provenance are preserved; the reservoir Jacobian cardinality is rebuilt from
the new natural-variable chart. The first regression covers a stable-cell
`2P -> 3P` restart followed by a separate accepted physical timestep.

The well module still does not define conductive well/reservoir heat exchange,
wellbore heat loss, rate-control equations, global well unknowns, control
switching, completion creation/deletion, cross-cell well migration or
scheduling.


## Single-well multi-connection fixed-BHP aggregation

One logical fixed-BHP well may now own multiple Peaceman connections targeting
unique stable reservoir cell `GlobalEntityId` values. Every connection shares
exactly one frozen BHP and remains an ordinary variable-cardinality cell source;
there is still no globally numbered well pressure unknown.

The PETSc multi-connection context sorts connections by stable cell identity,
rejects duplicate stable-cell completions and rejects inconsistent BHP values.
The existing mixed-cardinality assembly invokes the source callback only for
locally owned cells, so a ghost copy of a completion never contributes a second
source row. The callback itself is stateless and only dispatches the unique
matching connection; it never accumulates well totals during SNES evaluation.

Whole-well component and energy rates are aggregated separately from already
authoritative connection evaluations. The aggregation is production-positive
and contains no synthetic cross-cell Jacobian: each connection Jacobian remains
in its own reservoir cell block, while distributed callers sum owner-side
connection primals and then perform their MPI reduction. This keeps Newton or
line-search callback re-entry from multiplying reported well totals.

The current regression spans one 2P connection owned by rank 0 and one 3P
connection owned by rank 1, verifies exactly two authoritative contributions,
and closes global component/energy conservation against their summed well
rate.

A multi-connection context can also rebind exactly one completion by stable
cell identity. The selected connection delegates to the same stable
physical-phase rebinding contract used by a single completion, while all other
connection objects are copied unchanged and the complete logical well is
revalidated for unique stable cells and one shared frozen BHP. The local
transition regression drives only cell30 through `2P -> 3P -> 2P`; cell60
remains frozen 3P. After each restart, the well still has exactly two
authoritative owner-side connection rates and the post-rebuild accepted
physical timestep closes component/energy conservation against their sum.

Rate control, global well unknowns, control switching, wellbore pressure drop,
crossflow control, multi-well networks and scheduling remain outside this
contract.


## Fixed-total-molar-rate control

The first monolithic rate-control bridge adds exactly one well unknown to the
frozen multi-connection reservoir system: bottom-hole pressure. Reservoir
scalar indices remain unchanged and the BHP scalar is appended as the final
PETSc global scalar, with one authoritative owner. The control residual is
production-positive total molar rate:

`R_w = sum_connections sum_components n_dot_i - n_dot_target`.

Each Newton evaluation injects the current global BHP into the existing
connection-local Peaceman source calculation. The reservoir block therefore
continues to use the validated variable-cardinality source physics, while the
augmented Jacobian contains all four analytic blocks: the existing reservoir
`J_rr`, reservoir-source derivatives with respect to BHP `J_rw`, the
whole-well total-rate derivatives with respect to every authoritative
completion's reservoir variables `J_wr`, and the summed explicit BHP
derivative `J_ww`. Component/energy source terms remain owner-only; the well
row is assembled from the same authoritative connection evaluations, so ghost
connections do not double count.

This baseline uses PETSc SNES Newton line search with GMRES + restricted ASM,
with a frozen row-equilibration vector built from the initial augmented
analytic Jacobian. The controlled solve freezes every cell phase set. BHP/rate
switching, BHP limits, surface-rate conversions, phase transitions during the
controlled solve, multi-well networks and wellbore pressure-drop models remain
outside this contract.


### Fixed-total-molar-rate controlled physical timestep

The fixed-total-molar-rate augmented solve now has a production physical-timestep
bridge. It deliberately reuses the existing adaptive timestep controller and
accepted-history clock rather than introducing a second retry/commit state
machine.

For every nonlinear trial the bridge changes only the reservoir trial
backward-Euler timestep, then recreates the one-scalar augmented
`[q_reservoir, p_bhp]` system from the same accepted reservoir state/history
and the same entry accepted BHP initial guess. A nonconverged augmented solve is
therefore disposable and may be cut back without rebasing reservoir history or
carrying a failed-trial BHP into the next retry.

On acceptance, the augmented solution is split at the existing numbering
boundary. Only the reservoir state is passed to
`commit_accepted_physical_timestep_3d()`, which owns accumulation-history and
physical-time advancement. The converged BHP is committed separately only as
the next physical timestep's nonlinear initial guess. BHP never enters the
backward-Euler reservoir history.

Terminal adaptive rejection restores the entry reservoir trial timestep and
the entry BHP evaluator state. Phase cardinality and rate-control mode remain
frozen throughout this bridge. BHP/rate switching, BHP limits, phase
transitions during the controlled solve, surface/phase-rate controls and
multi-well control remain outside this contract.


### Minimum-BHP rate-to-BHP control switching

The fixed-total-molar-rate physical-timestep driver now accepts an optional
producer minimum-BHP constraint. Fixed-rate control remains the entry mode. If
its augmented nonlinear solve converges with

`p_bhp < p_min`,

that 43-scalar candidate is diagnostic evidence only and is destroyed without
advancing reservoir history or physical time. The source evaluator is then
fixed at exactly `p_min`, and the existing 42-scalar reservoir-only nonlinear
system is re-solved at the same physical timestep from the original accepted
reservoir state/history.

Once the minimum-BHP constraint is triggered it is sticky for all nonlinear
retries of that physical timestep. The source context explicitly enters a
reservoir-only fixed-BHP evaluation mode: Peaceman source physics and owner-only
cell insertion are unchanged, while the authoritative connection sidecar used
only by the augmented well row is disabled. This prevents repeated reservoir
SNES function/Jacobian evaluations from being misclassified as duplicate
completions. If the fixed-BHP re-solve diverges, normal adaptive timestep
cutback applies, but the retry stays in minimum-BHP mode and still starts from
the accepted reservoir baseline. A terminal timestep rejection discards the
trial switch and restores the entry BHP/dt; there is no accepted control-mode
mutation without an accepted physical step.

On acceptance, only the final fixed-BHP reservoir state is committed. The
accepted BHP is exactly `p_min`; the discarded rate candidate is retained in
the driver report for audit but never enters backward-Euler history. Because the
BHP limit is binding, the final well rate is allowed to differ from the original
rate target. BHP-to-rate switching, maximum-BHP constraints, surface/phase-rate
controls, multi-well control priority and schedule logic remain outside this
contract.


### Accepted well-control state across timesteps

Well control mode is now part of explicit accepted state rather than a
single-driver-call flag. `AcceptedFixedTotalMolarRateWellControlState3D`
stores both the accepted control mode and its BHP value.

A fixed-rate accepted step retains `fixed_total_molar_rate` and stores the
converged BHP only as the next rate-control nonlinear initial guess. When a
minimum-BHP switch is accepted, the state is committed atomically as
`minimum_bottom_hole_pressure` with BHP exactly equal to `p_min`. The next
physical timestep reads that accepted state and enters reservoir-only fixed-BHP
mode immediately; it does not create or solve the augmented rate-control system
first.

Control state follows the same transactional boundary as reservoir history and
the accepted physical clock. Nonlinear retry and timestep cutback never mutate
the accepted control state. Terminal rejection restores the entry source mode
and leaves accepted mode/BHP unchanged. No automatic minimum-BHP-to-rate
switchback or hysteresis is introduced in this contract.

The former `double* accepted_bottom_hole_pressure_pa` driver overload remains
as a single-step compatibility entry and reconstructs fixed-rate mode on each
call. Multi-timestep callers that require persistent control ownership use the
accepted well-control state overload.


### Guarded minimum-BHP-to-rate reactivation

An accepted minimum-BHP well can now reactivate fixed-total-molar-rate control,
but only through two explicit positive hysteresis margins:
`minimum_bhp_release_rate_margin_mol_per_s` and
`minimum_bhp_release_pressure_margin_pa`. The margins are configured as a
pair; omitting either disables reactivation rather than silently guessing a
scale.

For an entry minimum-BHP timestep, the reservoir-only `p_bhp=p_min` solve is
always performed first. Its converged state is a feasibility probe. The driver
re-evaluates the owner-only connection sources once at that probe state and
MPI-reduces whole-well molar production. No augmented rate system is created
unless

`q_bhp >= q_target + delta_q_release`.

If the capacity guard passes, the driver solves the 43-scalar augmented
fixed-total-molar-rate system from the same accepted reservoir history and the
same physical timestep. The rate candidate is accepted only when

`p_bhp_rate >= p_min + delta_p_release`.

When both guards pass, the fixed-BHP probe is discarded and only the augmented
rate candidate is committed; accepted control state becomes
`fixed_total_molar_rate`. If the pressure guard fails, or if the reactivation
rate solve does not converge, the already converged fixed-BHP probe remains the
final candidate and accepted control stays minimum-BHP. Thus reactivation
cannot invalidate an otherwise usable physical timestep.

A rate-to-BHP switch within the same timestep never immediately attempts the
reverse transition: reactivation is considered only when the *entry accepted*
control mode was already minimum-BHP. This prevents same-step control
oscillation. The accepted rate mode then persists normally into subsequent
physical timesteps.


### PETSc-independent single-well control arbitration kernel

The `rate <-> minimum-BHP` control policy is now owned by the core
`mpmc::well` module in `single_well_control_policy.hpp`. The kernel contains
no PETSc, MPI, reservoir-assembly or Peaceman-evaluation types. It owns control
mode, accepted control state, policy validation, and the three pure arbitration
decisions used by the production driver.

The comparison semantics are frozen and unit tested: `p_rate < p_min` switches
to BHP while equality remains rate controlled; `q_bhp >= q_target + delta_q`
permits a rate reactivation probe including equality; and
`p_rate >= p_min + delta_p` accepts reactivated rate control including
equality.

The PETSc physical-timestep driver retains execution responsibilities only:
constructing and solving candidates, owner-only source evaluation, MPI
reduction, adaptive retry/cutback and transactional reservoir/time commit. It
passes scalar candidate summaries to the core arbitration kernel and consumes
the returned decisions. Existing PETSc-facing control-mode and accepted-state
type names are compatibility aliases to the core types.

A separate well-core target tests hold/switch decisions, capacity and pressure
deadbands, all equality boundaries, disabled reactivation, accepted-state
construction and invalid policy configuration without PETSc or MPI.
