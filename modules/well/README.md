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

The well module still does not define rate-control equations, component
splitting, mass/molar conversion, multi-connection aggregation, wellbore
hydrostatics/friction, source assembly, PETSc objects or scheduling.
