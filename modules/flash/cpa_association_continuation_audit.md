# CPA association-state continuity audit

## Scope

This is a **test-only audit before optimization**. It does not modify the
production CPA association solver, density-root search, stability search, flash
iteration, tolerances, root topology rules, or failure semantics.

The audit answers one question only:

> Within one bounded CPA density-root search, can a converged association site-
> fraction state at a nearby density be reused as the initial fixed-point state
> at the next density evaluation without changing the validated thermodynamic
> result, and how many fixed-point sweeps could this avoid after conservative
> fallback accounting?

No continuation/reuse path is enabled in production by this slice.

## Fixed physical workload

The audit uses the same validated 333.15 K methanol-water five-state VLE fixture
as the CPA performance baseline:

- existing Kontogeorgis/Folas literature parameter snapshot;
- five Kurihara pressure states;
- existing Kurihara liquid/vapor compositions as initialization only;
- `automatic_starts=false` for initial and final stability, exactly as in the
  established physical-validation path;
- unchanged production `CpaPtPhase`, root options, stability and PT split logic.

The current production solve is run first. A recording provider then repeats the
same flash path and records every production `CpaPtPhase::roots` request. Its
outer status, evaluation counts, phase fractions and phase compositions must
match the production solve before any continuation experiment is considered.

## Root-search replay contract

Every recorded root request is replayed test-side using the current production
root-search ordering and **cold production phase pressure** to drive all scan
signs, brackets, bisections and root candidates.

The continuation shadow is evaluated beside that cold path only. It is never
allowed to select a bracket, alter a pressure sign, change a density candidate,
or change root topology.

The replay must reproduce for every root request:

- `CpaPtRootStatus`;
- production `roots.evaluations` count;
- root count;
- root density and slope sign;
- cold `ln(phi)` to a strict replay-only numerical guard.

This validates that the shadow experiment covers the same density evaluations
as production rather than a simplified surrogate scan.

## Shadow continuation policy

The association cache is local to **one root search** and is cleared before the
next root search. For each new reduced density `u`:

1. the normal production density-state phase evaluation is performed cold and is
   the scientific reference;
2. the nearest previously accepted density state in that root search is chosen;
3. its converged site fractions are used only as the initial state of a test-side
   reproduction of the existing fixed-point equations;
4. if no previous state exists, the first density uses the normal all-ones cold
   initialization;
5. if the warm shadow does not converge or fails an equivalence guard, the audit
   counts a deterministic fallback to the cold production result;
6. the accepted warm or fallback-cold site fractions are cached for subsequent
   nearby densities.

Nearest-density reuse is used rather than merely the previous execution in time,
because the bisection stage can jump back from the end of the monotone scan into
a local bracket. The audit is about thermodynamic continuity in density space.

## Predeclared equivalence guards

These guards are fixed before inspecting continuation results:

- site fractions: `max |Delta X_A| <= 1e-10` (audit-only internal guard);
- total phase pressure: reuse the frozen ThermoPack parity v1
  `max |Delta P| <= 5e-6 Pa`;
- accepted-root fugacity coefficients: reuse the frozen ThermoPack parity v1
  `max |Delta ln(phi)| <= 1e-10`.

The site-fraction guard is not a new production tolerance. The production
fixed-point stopping tolerance remains unchanged at `1e-12`.

A warm attempt that violates a guard is not rescued by loosening a tolerance; it
is counted as a fallback. Root-level `ln(phi)` violation also adds the cold solve
cost to the effective sweep accounting.

## Sweep accounting

The audit records:

- total production cold fixed-point sweeps;
- seeded warm-attempt sweeps;
- warm-better / equal / worse counts;
- non-convergence, site-fraction, pressure and root-`ln(phi)` fallback counts;
- exact-density cache hits;
- effective sweeps after charging every fallback for both the attempted warm
  solve and the required cold solve;
- the resulting effective sweep-reduction fraction.

The production cold reference computations performed by the audit itself are
not counted as hypothetical future continuation overhead. The effective counter
models a conservative `warm attempt -> cold fallback when needed` strategy.

## Interpretation boundary

A large sweep reduction with low fallback rate would only establish that
association-state continuation is a plausible next production optimization
candidate for this fixed two-phase workload. It would **not** authorize changing
root scan density, pressure/root tolerances, association convergence tolerance,
stability checks, or flash acceptance logic.

Before any production change, the continuation mechanism would still need a
separate design/review that preserves deterministic fallback and the already
frozen CPA/ThermoPack parity and physical-validation gates.
