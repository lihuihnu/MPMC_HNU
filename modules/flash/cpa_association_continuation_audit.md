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

These guards were fixed before inspecting continuation results:

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

## Hosted audit result

Workflow:

```text
CPA association continuity audit
run id: 35217384468
attempts: 1 and 2
runner: official ubuntu-24.04
compiler: Clang 18.1.3 Release
head: a1f70012ce3dd55f94f1a3889c759ed2e69cd638
```

The second attempt was an explicit rerun of the same SHA on a different hosted
Azure region. **Every printed count and floating-point result below was identical
between the two attempts.** The result is therefore treated as deterministic
algorithmic evidence rather than hosted wall-time evidence.

### Production-path coverage

The production and recording providers agree for all five flashes. Recorded root
calls by pressure are:

| P / Pa | root calls |
| ---: | ---: |
| 39,223 | 93 |
| 48,852 | 98 |
| 56,652 | 108 |
| 63,998 | 110 |
| 72,832 | 88 |
| **total** | **497** |

The root-search replay then covers exactly:

```text
root searches       = 497
density evaluations = 305,928
```

These counts match the previously established CPA hot-path performance baseline.
The replay reproduced the production root-search status, evaluation count, root
count, root density/slope and cold `ln(phi)` for every recorded request.

### Continuation/fallback result

```text
seeded density evaluations       = 305,431
exact-density seed hits           = 1,491
accepted continuations            = 297,415
fallback evaluations              = 8,016
non-convergence fallbacks         = 0
site-fraction-guard fallbacks     = 0
pressure-guard fallbacks          = 8,016
root-ln(phi)-guard fallbacks      = 0
fallback rate                     = 0.02624488018570479
```

Thus **97.3755%** of seeded density evaluations pass the predeclared continuation
acceptance guards directly. Every fallback is caused by the deliberately strict
`5e-6 Pa` total-pressure guard; none is caused by failure to converge, the
`1e-10` site-fraction guard, or the `1e-10` accepted-root `ln(phi)` guard.

The raw maximum differences before fallback are:

```text
max |Delta X_A|      = 5.9177662770082406e-13
max |Delta P|        = 2.44140625e-4 Pa
max |Delta ln(phi)|  = 1.8989254613188677e-12
```

The site-fraction and accepted-root fugacity differences are far inside their
predeclared guards. The pressure observable is more sensitive: extremely small
site-fraction differences can be amplified by the large association-pressure
term and physical/association cancellation at dense states. This is precisely
why a future production continuation design must retain an observable-based
fallback rather than assuming that convergence of `X_A` alone is sufficient.

### Fixed-point work

Production cold baseline:

```text
cold fixed-point sweeps = 12,548,901
mean cold sweeps / association solve = 41.019131952616299
max cold sweeps = 43
```

Nearest-density shadow continuation:

```text
mean seeded warm-attempt sweeps = 21.989074455441653
max warm sweeps                 = 34
warm better                     = 297,071
warm equal                      = 344
warm worse                      = 0
```

After conservatively charging every rejected warm attempt for both its warm
work and a complete cold fallback:

```text
effective sweeps with fallback = 7,068,639
sweep reduction fraction       = 0.43671250574054254
```

So this fixed physical workload shows an **effective 43.67% reduction in
association fixed-point sweeps even with the strict pressure-based fallback**.
The raw seeded warm solve reduces the mean sweep count from about `41.02` to
`21.99`, and no accepted seeded state requires more sweeps than its cold solve.

## Scientific conclusion

For this validated methanol-water two-phase workload, association site fractions
are strongly continuous across neighboring density states inside one root search.
Nearest-density continuation is therefore a credible production optimization
candidate, not merely a microbenchmark artifact.

The audit also identifies the required safety property: **continuation must be
advisory, with deterministic cold fallback based on thermodynamic observables.**
The pressure guard catches 2.62% of seeded evaluations even though `X_A` and root
`ln(phi)` remain extremely close to the cold solution.

This result does not justify weakening the association tolerance, root pressure
tolerance, density scan, stability review, or flash acceptance gates.

## Interpretation boundary

A large sweep reduction with low fallback rate only establishes that
association-state continuation is a plausible next production optimization
candidate for this fixed two-phase workload. It does **not** authorize changing
root scan density, pressure/root tolerances, association convergence tolerance,
stability checks, or flash acceptance logic.

Before any production change, the continuation mechanism still needs a separate
design/review that preserves:

- per-root-search state ownership (no cross-composition/cross-state leakage);
- nearest-density or otherwise explicitly bounded seed selection;
- unchanged fixed-point equations and stopping tolerance;
- deterministic cold fallback;
- unchanged root topology and pressure signs;
- frozen CPA/ThermoPack parity gates and existing physical validation.
