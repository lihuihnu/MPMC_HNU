# CPA ThermoPack formulation mapping and first external oracle

## Scope

This audit establishes whether the repository's existing 333.15 K
methanol(2B)-water(4C) CPA validation fixture can be compared meaningfully with
ThermoPack before any production CPA implementation is changed.

This increment is deliberately external-reference only:

- no `modules/thermodynamics` CPA implementation is modified;
- no `modules/flash/include` CPA solver/backend implementation is modified;
- no production tolerance is changed;
- no ThermoPack output is used to fit CPA parameters;
- no three-phase claim is made.

The external program is pinned to:

```text
repository: thermotools/thermopack
commit:     d68c794c7342bfc6938eb424a1fbb88b7780b738
```

The corresponding upstream hosted `unittests` and `cibuildwheel` runs completed
successfully at that commit. MPMC_HNU nevertheless builds the pinned source
itself on an official GitHub-hosted runner before generating the oracle.

## Formulation mapping

### Physical contribution and alpha

MPMC_HNU profile:

```text
CPA/SRK-physical/simplified-rdf-1.9eta/explicit-site-pairs/v1
```

ThermoPack is initialized through `SRK_CPA("MEOH,H2O", mixing="vdW",
alpha="Classic")`. ThermoPack's CPA memo states that the common SRK form uses

```text
A^CPA = A^ideal + A^SRK + A^assoc
```

and conventional quadratic-`a` / linear-`b` cubic mixing:

```text
a_ij = sqrt(a_i a_j) (1-k_ij)
b_mix = sum_i x_i b_i
```

This matches the physical contribution and mixing convention used by the current
MPMC_HNU CPA profile for this binary fixture.

### Simplified radial distribution function

ThermoPack pinned source (`src/saft_rdf.f90`) selects, when simplified CPA is
active,

```text
g = 1 / (1 - 1.9 eta)
eta = b_mix rho / 4
```

which is the same radial-distribution expression frozen by the MPMC_HNU profile.
The oracle generator explicitly calls `set_cpa_formulation(True, False)`; it does
not rely on an implicit formulation default.

### Association equations and CR-1

ThermoPack's CPA memo gives the same site-fraction structure used by MPMC_HNU,

```text
X_Ai = 1 / [1 + rho sum_j x_j sum_B X_Bj Delta_AiBj]
```

with

```text
Delta_AiBj = g [exp(epsilon_AiBj / RT) - 1] b_ij beta_AiBj
b_ij = (b_i + b_j)/2
```

and documents the CR-1 rules

```text
epsilon_cross = (epsilon_i + epsilon_j)/2
beta_cross    = sqrt(beta_i beta_j)
```

The pinned `MEOH/H2O` CPA binary record uses `ARITHMETIC` epsilon and
`GEOMETRIC` beta combining rules. Setting the optional Elliott formulation flag
to `False` retains ThermoPack's STANDARD association-Delta rule, i.e. this CR-1
route rather than Elliott's rule.

## Pure-parameter mapping

ThermoPack's public CPA tuning API defines the five-value vector as

```text
[a0, b, epsilon, beta, c1]
a0      [Pa L^2 mol^-2]
b       [L mol^-1]
epsilon [J mol^-1]
beta    [-]
c1      [-]
```

The MPMC_HNU fixture stores the physical parameters in SI. Converting only the
volume unit gives an exact numerical mapping to ThermoPack's public API:

| Component | Site scheme | Quantity | MPMC_HNU | ThermoPack API value |
| --- | --- | --- | ---: | ---: |
| Methanol | 2B | `a0` | `0.40531 Pa m^6 mol^-2` | `4.0531e5 Pa L^2 mol^-2` |
| Methanol | 2B | `b` | `3.0978e-5 m^3 mol^-1` | `0.030978 L mol^-1` |
| Methanol | 2B | `epsilon` | `24591 J mol^-1` | `24591 J mol^-1` |
| Methanol | 2B | `beta` | `0.0161` | `0.0161` |
| Methanol | 2B | `c1` | `0.43102` | `0.43102` |
| Water | 4C | `a0` | `0.12277 Pa m^6 mol^-2` | `1.2277e5 Pa L^2 mol^-2` |
| Water | 4C | `b` | `1.4515e-5 m^3 mol^-1` | `0.014515 L mol^-1` |
| Water | 4C | `epsilon` | `16655 J mol^-1` | `16655 J mol^-1` |
| Water | 4C | `beta` | `0.0692` | `0.0692` |
| Water | 4C | `c1` | `0.67359` | `0.67359` |

The pinned ThermoPack `fluids/Methanol.json` and `fluids/Water.json` contain the
same numerical CPA pure parameters and 2B/4C association schemes. The generator
still calls `set_pure_params(...)` explicitly and verifies `get_pure_params(...)`
read-back, so future database-default drift cannot silently change the oracle.

## Binary-parameter mismatch that must be corrected before comparison

The MPMC_HNU physical fixture uses the Folas 333.15 K CR-1 value

```text
kij_a(MEOH,H2O) = -0.055
```

The pinned ThermoPack default CPA binary record instead contains

```text
kij_a(MEOH,H2O) = -0.09
kij_eps          = 0
beta rule        = GEOMETRIC
```

Therefore **ThermoPack default initialization is not a valid oracle for the
MPMC_HNU fixture**. The independent generator explicitly executes

```text
set_kij(1, 2, -0.055, 0.0)
```

and requires `get_kij(1,2)` to return that exact matched pair before any flash is
run. This is parameter alignment, not fitting: `-0.055` is the existing
literature parameter already frozen by the repository before this comparison.

## Oracle state construction

The first oracle uses exactly the five traceable Kurihara et al. 333.15 K
pressure points already owned by `cpa_physical_validation`. For each state, the
overall methanol fraction is constructed as

```text
z_MEOH = 0.9 x_MEOH,experimental + 0.1 y_MEOH,experimental
```

matching the existing repository carrier-feed convention. The experimental
phase compositions are **not** passed to ThermoPack as equilibrium starts or
expected results; only `T`, `P`, and `z` are supplied to ThermoPack's own
`two_phase_tpflash` routine.

The generated oracle records:

- ThermoPack revision and model switches;
- pure-parameter and `kij` read-back;
- `T`, `P`, and overall feed;
- ThermoPack phase code;
- vapor and liquid mole phase fractions;
- ThermoPack liquid and vapor component mole fractions;
- a material-balance reconstruction residual used only to validate extraction.

The generator is
`tests/flash/cpa_thermopack_oracle/generate_oracle.py`. It imports ThermoPack and
Python standard-library modules only; it does not import MPMC_HNU production or
test code.

## Reproduction gate

`.github/workflows/cpa_thermopack_oracle.yml` runs on an official
`ubuntu-24.04` hosted runner and:

1. checks out the MPMC_HNU PR;
2. clones `thermotools/thermopack`;
3. detaches exactly at `d68c794c7342bfc6938eb424a1fbb88b7780b738`;
4. builds the pinned ThermoPack source in Release mode with its documented CMake
   path;
5. configures the upstream Python wrapper;
6. runs the independent generator.

The numerical oracle is frozen into the repository only after this hosted run
succeeds. A failed external build, parameter read-back, two-phase classification,
or material-balance extraction is a blocker; no MPMC_HNU production change is
permitted in this slice to work around it.

## Source chain

ThermoPack pinned sources used by the mapping audit:

- `addon/pycThermopack/thermopack/cpa.py` — SRK-CPA initialization,
  `set_pure_params`, `get_pure_params`, `set_kij`, `get_kij`, formulation switch;
- `src/saft_rdf.f90` — simplified `1/(1-1.9 eta)` RDF;
- `docs/memo/CPA/cpa.tex` — CPA association equation and CR-1 rules;
- `fluids/Methanol.json` — methanol CPA pure parameters / 2B scheme;
- `fluids/Water.json` — water CPA pure parameters / 4C scheme;
- `binaries/CPA.json` — MEOH/H2O default binary record and combining rules;
- `addon/pycThermopack/thermopack/utils.py` — TP `FlashResult` fields.

MPMC_HNU source chain remains the already documented Kontogeorgis/Folas/Kurihara
chain in `cpa_physical_validation.md` and
`tests/flash/cpa_physical_validation/test_support.hpp`.

## Explicit non-claims

This mapping proves only that a formulation-matched **two-phase** external oracle
can be generated without changing the production CPA implementation. It does not
yet prove MPMC_HNU/ThermoPack numerical parity, does not establish CPA VLLE
accuracy, does not promote a numerical root side to physical morphology, and does
not alter `global_stability_proven=false`.
