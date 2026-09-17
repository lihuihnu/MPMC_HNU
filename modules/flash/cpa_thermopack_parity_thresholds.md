# CPA ThermoPack numerical-parity acceptance envelope

## Scope

This document freezes the first numerical acceptance envelope for the named CPA
parity snapshot:

```text
dataset_id:
  thermopack-d68c794__meoh-h2o__mpmc-literature-folas-cr1-aligned
revision:
  v1__Tc-thermopack-component-db__other-cpa-fields-preserved-from-literature-fixture
external software:
  thermotools/thermopack@d68c794c7342bfc6938eb424a1fbb88b7780b738
threshold contract:
  MPMC_HNU/CPA/ThermoPack-parity-thresholds/v1
```

The envelope applies only to the five frozen methanol-water TP flash states at
333.15 K and the corresponding ten liquid/vapor phase-property states. It is a
**cross-implementation numerical parity contract**, not an experimental-accuracy
claim, not a general CPA error bound, and not a proof of global Gibbs stability.

No production CPA parameter, EOS/association equation, density-root solver,
stability algorithm, flash orchestration or production tolerance is changed by
this contract.

## Evidence used before freezing the thresholds

The thresholds were selected only after the formulation and parameter identity
were locked and before any production CPA implementation change.

Relevant official GitHub-hosted evidence includes:

- CPA physical-validation run `35208863756`, Ubuntu 24.04 / Clang Release;
- CPA physical-validation run `35210263381`, Ubuntu 24.04 / Clang Release;
- CPA physical-validation run `35210263381`, Windows Server 2022 / MSVC Release;
- pinned ThermoPack external-oracle run `35210263457`, which rebuilt the exact
  external revision and reproduced both frozen external references.

The two independent Clang hosted runs reproduced the parity summary at printed
double precision. The Windows/MSVC run showed only last-digit platform variation
in flash observables; phase-kernel maxima were effectively unchanged.

Observed maxima before the gate was introduced were:

| Observable | Ubuntu/Clang observed max | Windows/MSVC observed max |
| --- | ---: | ---: |
| `|Δ beta_vapor|` | `1.1870388436e-7` | `1.1870390812e-7` |
| `|Δ x_MeOH|` | `4.0481111463e-8` | `4.0481118957e-8` |
| `|Δ y_MeOH|` | `1.5360046568e-11` | `1.5361600880e-11` |
| relative density error | `4.0089206520e-11` | `4.0089206520e-11` |
| `|ΔZ|` | `3.7996050750e-11` | `3.7996050750e-11` |
| `|Δ ln(phi)|` | `5.8310710427e-11` | `5.8310710427e-11` |
| common-TV `|ΔP_physical|` / Pa | `2.0861625671e-7` | `2.0861625671e-7` |
| common-TV `|ΔP_association|` / Pa | `3.8444995880e-6` | `3.8444995880e-6` |
| common-TV `|ΔP_total|` / Pa | `2.4437904358e-6` | `2.4437904358e-6` |
| common-TV `|Δ(mu_cubic/RT)|` | `1.7763568394e-15` | `1.7763568394e-15` |
| common-TV `|Δ(mu_assoc/RT)|` | `3.3396049814e-11` | `3.3396049814e-11` |

The envelope is intentionally expressed as simple fixed values with finite
cross-platform margin. It is not defined as `current_result * factor`, so a
future implementation cannot silently move the gate by changing the baseline.

## Frozen v1 acceptance envelope

### Five-state TP flash

For every frozen state, existing phase-count, material-balance, fugacity and
final-stability requirements remain mandatory. In addition, the maximum absolute
MPMC_HNU minus ThermoPack error over the five states must satisfy:

```text
max |Δ beta_vapor| <= 2.0e-7
max |Δ x_MeOH|     <= 6.0e-8
max |Δ y_MeOH|     <= 3.0e-11
```

These gates apply to the **named parity snapshot only**. The separate literature
snapshot remains a physical-validation dataset and is not required to satisfy
this external-software parity envelope.

### Ten phase-property states

For the lower-/upper-density mechanically admissible roots selected by the
existing CPA test convention at the frozen ThermoPack phase compositions:

```text
max relative |Δ rho| <= 1.0e-10
max |Δ Z|            <= 1.0e-10
max |Δ ln(phi_i)|    <= 1.0e-10
```

For the common-`(T,V,x)` decomposition audit:

```text
max |Δ P_physical|       <= 5.0e-7 Pa
max |Δ P_association|    <= 5.0e-6 Pa
max |Δ P_total|          <= 5.0e-6 Pa
max |Δ (mu_cubic/RT)|    <= 5.0e-15
max |Δ (mu_association/RT)| <= 5.0e-11
```

The pressure and chemical-potential decomposition checks are diagnostic parity
checks for this pinned formulation. They do not redefine the production solver's
pressure/root tolerances.

## Change-control rule

This is a **pre-production-change gate**. Once production CPA work begins, a
failed parity test must not be made green by simply increasing these constants in
the same corrective change.

Changing the envelope requires all of the following:

1. explicit evidence that the pinned external reference, parameter snapshot or
   supported numerical environment changed in a scientifically relevant way;
2. an updated evidence table from independent hosted runs;
3. an explicit new threshold-contract version rather than silently editing `v1`;
4. review of whether the change represents numerical portability, a parameter
   revision, or an actual loss of thermodynamic parity.

The machine-readable constants live in:

```text
tests/flash/cpa_physical_validation/cpa_thermopack_parity_thresholds.hpp
```

The flash and phase-kernel parity tests consume those constants directly.
