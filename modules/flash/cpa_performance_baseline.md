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

The normal full workload also records counters already exposed by production
result objects, including:

- initial and final stability property evaluations;
- stability trial count and trial iterations;
- split property evaluations;
- split attempt count, attempt evaluations, iterations and backtracks.

These are deterministic algorithmic counts for the fixed state/workload, not
clock measurements.

## Exact internal structural counts without production profiling hooks

The production API does not expose an accumulated count of every density-state
property evaluation or every association fixed-point sweep performed inside all
root searches. Adding runtime profiling hooks to production solely for this
audit would itself change the hot path.

Instead, a separate GCC Release build enables gcov **only for structural
counting**. The same five-state workload is run once and line execution counts
are extracted from the existing production headers for:

- `evaluate_cpa_phase_at_density` executions;
- `CpaPtPhase::roots` calls;
- root-search density evaluations (`++result.evaluations`);
- `solve_cpa_association` calls;
- association fixed-point iteration sweeps;
- CPA stability-adapter evaluations;
- CPA split-adapter evaluations.

The collector also verifies two invariants of the current implementation:

```text
phase evaluations == association solve calls
root searches == stability-adapter evaluations + split-adapter evaluations
```

Coverage-instrumented elapsed time is never used as performance evidence.

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
