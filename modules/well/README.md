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

The well module still does not connect this adapter to the production source
evaluator/PETSc assembly, and still does not define conductive well/reservoir
heat exchange, wellbore heat loss, multi-connection aggregation, rate-control
equations, global well unknowns or scheduling.
