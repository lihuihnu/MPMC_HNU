# CPA flash hot-path performance baseline

## Scope

This baseline is an audit before optimization. It does **not** change CPA
thermodynamics, root selection, association iteration, stability search, flash
iteration, tolerances, or failure semantics, and it does not create a hosted
wall-time performance gate.

The purpose is to answer two questions before any production optimization:

1. does the current PR show an observable performance regression relative to the
   current `main` revision on the same hosted runner; and
2. which nested part of the validated CPA VLE workload is structurally
   responsible for most repeated work?

## Fixed workload

The workload is the existing independent 333.15 K methanol-water physical
validation fixture:

- five Kurihara pressure states;
- the existing Kontogeorgis/Folas literature CPA parameter snapshot;
- the existing experimental liquid/vapor compositions as initialization only;
- normalized feed as already defined by the physical-validation regression;
- `automatic_starts=false` for initial and final stability, exactly matching the
  established physical regression path;
- current production `CpaPtPhase`, `CpaVleEvaluator`, stability and PT split
  implementations and their unchanged numerical defaults.

This baseline therefore measures the validated warm-start two-phase path. It is
not a cold-start benchmark and does not claim to characterize every CPA state or
max3/VLLE orchestration.

## Same-runner main/head comparison

Workflow:

```text
.github/workflows/cpa_performance_audit.yml
```

The workflow runs on the official `ubuntu-24.04` image and records exact
`main`/PR-head SHAs. Both production trees are built inside the **same job** from
the same benchmark source. Normal timing binaries use Clang Release.

Measurements are interleaved rather than running all `main` samples followed by
all head samples:

- five micro-hot-path samples for each revision;
- three complete five-state flash samples for each revision;
- odd samples run `main -> head`, even samples run `head -> main`.

The report records medians, ranges and median absolute deviations. Hosted wall
time remains descriptive evidence only. No single run, median ratio, or runner
sample is promoted to a hard performance acceptance threshold by this audit.

The performance-audit workflow deliberately does not cancel an already-started
paired run when an unrelated commit lands on the same PR branch. Every run pins
and prints its own exact head SHA, so a completed main/head pair remains an
internally consistent observation rather than mixing revisions.

## Timed hot-path layers

The uninstrumented Release micro benchmark times the same fixed phase states at
five pressures for:

- direct association solve;
- density-state CPA phase-property evaluation;
- complete density-root search;
- CPA stability adapter phase evaluation;
- CPA split adapter phase evaluation;
- complete five-state PT VLE flash workload.

These layers are nested. Their durations must **not** be added together. The
micro measurements are used to understand per-call cost and nesting, while the
full-flash timing is the end-to-end observation.

## Deterministic solver counters

The normal full workload records counters already exposed by production result
objects, including:

- initial and final stability property evaluations;
- stability trial count and trial iterations;
- split property evaluations;
- split attempt count, attempt evaluations, iterations and backtracks.

These are the authoritative stability/flash evaluation counts for the fixed
workload. They are deterministic algorithmic counts, not clock measurements.

## Exact internal structural counts without production profiling hooks

The production API does not expose an accumulated count of every density-state
property evaluation or every association fixed-point sweep performed inside all
root searches. Adding runtime profiling hooks to production solely for this
audit would itself change the hot path.

Instead, a separate GCC Release build enables gcov **only for structural
counting**. The same five-state workload is run once and line execution counts
are extracted from existing production headers for:

- `CpaPtPhase::roots` calls;
- root-search density/property evaluations, including both existing
  `++result.evaluations` sites;
- `evaluate_cpa_phase_at_density` executions;
- `solve_cpa_association` calls;
- association fixed-point iteration sweeps.

Stability-adapter and split-adapter counts are intentionally **not** inferred
from gcov because those thin inline wrappers can lose independent line records
under Release optimization. Their counts are already available directly from
`StabilityResult.evaluations` and `PtSplitResult.split_evaluations`.

The collector verifies the current thermodynamic nesting invariant:

```text
phase evaluations == association solve calls
```

Coverage-instrumented elapsed time is never used as performance evidence.

## Frozen baseline observation: official hosted run #12

Successful workflow run:

```text
CPA flash performance audit / run 12
run id: 35215794105
runner: ubuntu-24.04
Clang: 18.1.3
GCC coverage counter build: 13.3.0
main: f9d65c9e5de03a6ac64fe96f6458faf311ae6fea
PR head: 510606debf432167782a6d4c12d877023922442f
```

### Repeated uninstrumented Release timing

| Metric | `main` median | PR-head median | head/main |
| --- | ---: | ---: | ---: |
| five-state full flash | 22.702 s | 22.619 s | 0.9963 |
| association solve / call | 67.866 us | 67.830 us | 0.9995 |
| density-state phase evaluation / call | 67.924 us | 68.090 us | 1.0024 |
| density-root search / call | 45.263 ms | 45.273 ms | 1.0002 |
| stability adapter / call | 44.853 ms | 44.880 ms | 1.0006 |
| split adapter / call | 45.423 ms | 45.636 ms | 1.0047 |

The generated report classifies this as:

```text
no_clear_hosted_runner_regression_signal
```

The five-state full-flash median differs by about `-0.37%`; individual micro
ratios straddle 1.0. Earlier paired observations on other hosted runners also
changed direction at roughly percent scale. Therefore this audit does **not**
claim a speedup and finds no reproducible regression signal.

### Deterministic outer solver work

`main` and the measured PR head are identical for every recorded outer counter:

```text
states                         = 5
initial_stability_evals        = 81
final_stability_evals          = 80
stability_trials               = 25
stability_trial_iterations     = 131
split_evals                    = 336
split_attempts                 = 10
split_attempt_evals            = 336
split_attempt_iterations       = 152
split_backtracks               = 650
```

Thus the full five-state workload makes `81 + 80 + 336 = 497` CPA phase-provider
calls/root searches through the stability and split paths.

### Exact nested production work

The Release gcov structural audit is also identical for `main` and PR head:

| Production hot-path count | `main` | PR head |
| --- | ---: | ---: |
| root-search calls | 497 | 497 |
| root density/property evaluations | 305,928 | 305,928 |
| density-state phase evaluations | 305,928 | 305,928 |
| association solves | 305,928 | 305,928 |
| association fixed-point sweeps | 12,548,901 | 12,548,901 |

Derived multipliers:

```text
phase evaluations / root search        = 615.5493
association sweeps / association solve = 41.0191
```

The independent micro measurements give:

```text
association solve time / density-state phase-evaluation time = 0.996
```

So the measured hot path has the following nested shape:

```text
497 root searches
  -> 305,928 density-state phase evaluations
     -> 305,928 association solves
        -> 12,548,901 association fixed-point sweeps
```

The median direct root call is about `45.27 ms`; multiplying the structural
`497` root calls by that micro scale is consistent with essentially the whole
`~22.6 s` full workload. This is supporting triangulation, not an exclusive
profiler attribution: the timed layers are nested and the micro states are not
summed as independent costs.

## Bottleneck conclusion

The first established CPA flash hot-path bottleneck is:

> **the density-root scan multiplying repeated association fixed-point solves.**

Two independent observations support this conclusion:

1. each root search performs about `615.55` density-state phase evaluations in
   the fixed full workload; and
2. a direct association solve consumes about `99.6%` of the direct
   density-state phase-evaluation time in the micro benchmark, with about
   `41.02` fixed-point sweeps per association solve in the full workload.

The audit does **not** conclude that the validated `scan_intervals=512` should be
reduced, that density/root tolerances should be relaxed, or that association
convergence criteria should be weakened. Those changes could alter root topology
or scientific acceptance and are outside this baseline.

## First production hotspot to investigate next

Before changing any production implementation, the first optimization-focused
audit should target **repeated cold association fixed-point work across adjacent
density evaluations inside one bounded root search**. In particular, evaluate
whether a strictly local association-state continuation/reuse strategy can
reduce repeated fixed-point work while preserving all of the following:

- the same density scan and root-topology coverage;
- the same pressure/root tolerances;
- the same converged association solution and failure semantics;
- deterministic fallback to the current cold solve when continuation is not
  acceptable;
- the frozen ThermoPack parity envelopes and existing physical validation.

That is a recommendation for the next audit only. No such optimization is
implemented by this baseline slice.

## Interpretation rules

The generated report may describe an observed hosted-runner signal, but it does
not fail CI based on timing. In particular:

- small main/head timing differences inside observed runner dispersion are not
  called a speedup or regression;
- a possible timing regression is reported as a signal requiring repetition,
  not as proof;
- deterministic structural-count changes are reported separately from timing;
- a bottleneck is identified from nesting plus repeated-work counts, not merely
  from whichever standalone microbenchmark has the largest wall time.

Any later optimization must preserve the already frozen CPA/ThermoPack parity
contract and all existing material-balance, fugacity, stability and failure
semantics. Performance changes should be compared against this fixed workload
before a production optimization is accepted.
