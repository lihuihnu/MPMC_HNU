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
wellbore heat loss, multi-connection aggregation, rate-control equations,
global well unknowns, control switching, completion creation/deletion,
cross-cell well migration or scheduling.


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
rate. Rate control, global well unknowns, control switching, wellbore pressure
drop, crossflow control, multi-well networks and scheduling remain outside this
contract.
