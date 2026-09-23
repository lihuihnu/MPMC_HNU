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

The hydraulic-conductance layer stops at

`C_alpha = WI * lambda_alpha [m3/(Pa s)]`.

It does not define pressure drawdown, phase/component well rates, BHP/rate
controls, well unknowns, wellbore hydraulics, hydrostatic correction, source
assembly, PETSc objects or scheduling. Those remain separate future layers.
