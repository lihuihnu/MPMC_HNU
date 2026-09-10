# Physics thermodynamic-closure consumption contract

## Scope

This is the first concrete `physics` layer increment. It consumes an already solved PR76 VLE PT flash result and exposes a small thermodynamic closure snapshot for future conservation-law and discretization code. It does **not** create mesh, flux, saturation, Darcy, time stepping, global residual assembly, a nonlinear solver, viscosity, enthalpy or an energy equation.

The public generic header is `mpmc/physics/thermodynamic_closure.hpp`; the current PR76 adapter is `mpmc/physics/pr76_thermodynamic_closure.hpp`. The CMake target is `mpmc::physics`, which depends on the opt-in `mpmc::flash_sensitivity` target. Dependency direction is therefore `physics -> flash sensitivity -> flash/thermodynamics + ad`; no lower layer depends on physics.

## Atomic snapshot and two-axis availability

`ThermodynamicClosureSnapshot` separates two decisions:

- `primal_status`: whether this snapshot owns a thermodynamic primal that may be used by a future conservation residual or line-search merit evaluation;
- `linearization_status`: whether this **same** snapshot owns a local derivative that may seed a Newton Jacobian.

The primal and derivative are `std::optional` values built from scratch on every call. A derivative from an older state is never copied into a new closure result by the adapter. A caller replacing its snapshot therefore cannot obtain `new primal + old derivative` from this API.

For the current adapter, an accepted two-phase flash can remain a valid primal when sensitivity is unavailable because of a phase-boundary guard, local ill-conditioning, unsupported feed support, derivative property failure or derivative arithmetic failure. In those cases `residual_available()==true`, `can_seed_newton()==false`, and `linearization` is empty. If the sensitivity path reports `solution_not_accepted`, or the flash itself is `phase_set_unstable`/`indeterminate`, the adapter does not publish a usable primal.

A stable single-phase flash is a valid thermodynamic regime in general, but the v1 **two-phase** consumption adapter intentionally does not implement its closure yet and reports `phase_regime_not_implemented`; it does not manufacture a second phase.

## Coordinates and quantities

The derivative columns are inherited exactly from the accepted flash sensitivity contract:

```text
q = (p_Pa, T_K, z_0, ..., z_{N-2})
z_{N-1} = 1 - sum(z_0, ..., z_{N-2})
```

The primal contains the liquid/vapor mole phase fractions, phase compositions, compressibility factors and phase molar densities. The flash phase fraction is explicitly a **mole phase fraction, not pore-volume saturation**.

The existing PR76 phase contract defines

```text
Z = p * v_m / (R*T)
```

so the adapter introduces no new empirical model when it computes

```text
c = 1/v_m = p/(Z*R*T)  [mol/m^3].
```

For an available linearization, the adapter obtains the local partial derivative of each selected PR76 `Z` branch with respect to `(p,T,w_reduced)` from the existing AD/simple-root IFT phase kernel. It then combines that with the converged flash `dx/dq` or `dy/dq`:

```text
dZ_alpha/dq = partial Z_alpha/partial(p,T,w_alpha) * d(p,T,w_alpha)/dq.
```

Molar-density derivatives use the exact chain rule

```text
dc_alpha/c_alpha = dp/p - dZ_alpha/Z_alpha - dT/T.
```

No flash iteration, RR bisection, line search, TPD search or root iteration is differentiated here.

## Intended future solver policy

This module does not implement Newton globalization, but it makes the required distinction machine-readable:

- future residual/merit evaluation may consume a snapshot only when `residual_available()` is true;
- a future Newton base/Jacobian assembly may consume the derivative only when `can_seed_newton()` is true;
- an unavailable derivative must not be replaced by zeros, a previous-cell/previous-iterate derivative, clipping, or a hidden finite-difference fallback;
- a future solver may deliberately implement a documented global inexact/quasi-Newton strategy, but that is a separate numerical algorithm and must not arise as an implicit per-cell fallback.

## Current non-capabilities

The adapter does not provide single-phase closure values, phase switching, bubble/dew or phase-disappearance derivatives, critical/root-switch derivatives, zero-support changes, LLE/three-phase states, saturation, mass density, mobility, transport properties, enthalpy/internal energy, SW or CPA closure. Those require separate model and validation increments.
