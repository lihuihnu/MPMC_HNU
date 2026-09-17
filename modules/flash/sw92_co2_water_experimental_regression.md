# SW92 CO2-H2O experimental absolute-error regression

## Purpose

This regression provides an **independent experimental magnitude check** for the corrected-original SW92 freshwater CO2/H2O PT-flash path. It is separate from formula goldens, Decimal(80) implementation references, PR/SW common-limit checks, and structural backend-contract tests.

The test does **not** refit SW92 coefficients or solver tolerances. It freezes the current corrected-original model against traceable experimental composition data and fails if future code changes materially degrade that agreement.

## Experimental source

Original measurements:

- A. Zawisza and B. Malesinska, **“Solubility of carbon dioxide in liquid water and of water in gaseous carbon dioxide in the range 0.2-5 MPa and at temperatures up to 473 K,”** *Journal of Chemical & Engineering Data* 26 (1981) 388–391.
- DOI: <https://doi.org/10.1021/je00026a012>

Machine-readable compilation used for the exact regression rows:

- IUPAC-NIST Solubility Database, SRD 106, CO2 + H2O, system `62_57`.
- <https://srdata.nist.gov/solubility/sol_detail.aspx?sysID=62_57>

NIST identifies component 1 as CO2 and reports the `x1` values used below. The NIST page states that these mole fractions were calculated by the compiler from the original measurements. Consequently, the repository acceptance budgets are **model regression envelopes**, not claims about the experimental uncertainty of `x_CO2`.

## Frozen freshwater points

The SW backend uses `nacl_molality_mol_per_kg_water = 0`. Nine low-pressure CO2/H2O points are frozen across three isotherms:

| T / K | P / bar | NIST `x_CO2` | SW92 `x_CO2` at audit | absolute error |
| ---: | ---: | ---: | ---: | ---: |
| 323.15 | 8.92 | 0.0030 | 0.0023452922 | 0.0006547078 |
| 323.15 | 14.41 | 0.0046 | 0.0037483466 | 0.0008516534 |
| 323.15 | 25.05 | 0.0076 | 0.0063259950 | 0.0012740050 |
| 348.15 | 2.37 | 0.0005 | 0.0003936212 | 0.0001063788 |
| 348.15 | 21.00 | 0.0046 | 0.0038930188 | 0.0007069812 |
| 348.15 | 35.91 | 0.0076 | 0.0064297038 | 0.0011702962 |
| 373.15 | 3.56 | 0.0005 | 0.0004156508 | 0.0000843492 |
| 373.15 | 27.10 | 0.0046 | 0.0040457373 | 0.0005542627 |
| 373.15 | 45.60 | 0.0076 | 0.0066145951 | 0.0009854049 |

Audit aggregate:

- maximum absolute error: `0.0012740050` mole fraction;
- mean absolute error: `0.0007097821` mole fraction.

## Acceptance gate

`tests/runtime/pt_service/sw92_co2_water_experimental_regression_test.cpp` requires:

- all nine states are structurally valid and scientifically accepted by the SW92 PT service path;
- each state closes as a two-phase CO2/H2O equilibrium;
- the water-richer composition coordinate is selected numerically by the largest water mole fraction, without publishing an aqueous/water phase identity;
- **maximum absolute `x_CO2` error <= 0.0015**;
- **mean absolute `x_CO2` error <= 0.0010**.

The two error thresholds correspond to 0.15 mol% maximum and 0.10 mol% mean absolute composition error. They are intentionally round regression budgets with margin above the audited current values. They must not be relaxed merely to keep CI green; any exceedance requires a scientific audit of equations, parameter routing, units, phase selection, and solver closure.

## Scope and interpretation

This is a validation of the repository's **SW92/corrected-original** implementation against an independent CO2/H2O experimental dataset. It does not claim that SW92 is the most accurate modern model, and it does not compare different SW refresh/refit parameterizations.

The overall service request uses an equimolar CO2/H2O feed only as a carrier that places the state inside the two-phase region. The compared quantity is the equilibrium CO2 mole fraction on the water-richer composition coordinate; it is an intensive phase-equilibrium result. The public role-neutral phase contract remains unchanged.
