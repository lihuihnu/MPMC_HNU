# CPA parameter snapshot identities for physical validation and ThermoPack parity

## Purpose

The methanol-water CPA work now uses two intentionally distinct parameter
snapshots. They share the same CPA formulation and almost all numerical fields,
but they are **not the same dataset** because Classic alpha consumes different
critical-temperature records.

This separation prevents an external-software comparison from silently changing
the repository's literature fixture, and prevents a literature validation result
from being mislabeled as ThermoPack-default parity.

No production CPA equation, solver, tolerance, phase selection or flash
orchestration is changed by these snapshots.

## Snapshot A: literature physical-validation dataset

Identity:

```text
dataset_id:
  literature-cpa-water-methanol-cr1-333.15K
revision:
  Kontogeorgis-2008-pure__Folas-CR1-kij__Kurihara-1995-VLE
```

This remains the repository's independent physical-validation dataset. Its key
sources are:

- Kontogeorgis et al. (2008), DOI `10.2516/ogst:2008025` for the CPA pure
  parameters and association model records;
- Folas thesis Table 2.1 for the 333.15 K water-methanol CR-1 cubic
  `kij=-0.055`;
- Kurihara et al. (1995), DOI `10.1021/je00019a033` for the five experimental
  VLE points.

The Classic-alpha critical temperatures remain exactly the literature values
already frozen before the ThermoPack comparison:

```text
METHANOL Tc = 512.64 K
WATER    Tc = 647.29 K
```

These values are **not** changed merely to match an external program.

## Snapshot B: pinned ThermoPack parity dataset

Identity:

```text
dataset_id:
  thermopack-d68c794__meoh-h2o__mpmc-literature-folas-cr1-aligned
revision:
  v1__Tc-thermopack-component-db__other-cpa-fields-preserved-from-literature-fixture
```

This is deliberately not called a `ThermoPack default` dataset. The independent
oracle does not use ThermoPack's complete default parameter selection. It builds
a formulation-matched comparison state with a mixed provenance chain:

| Field group | Value origin in parity snapshot | Required provenance |
| --- | --- | --- |
| `Tc`, consumed by Classic alpha | pinned ThermoPack component database | `thermotools/thermopack@d68c794c7342bfc6938eb424a1fbb88b7780b738`, `fluids/Methanol.json` / `fluids/Water.json` |
| `a0`, `b`, `c1` | unchanged repository literature fixture; explicitly written to ThermoPack and read back | original Kontogeorgis literature provenance retained verbatim |
| association site topology | unchanged repository literature fixture | original literature provenance / declared site records retained |
| self `epsilon`, `beta` | unchanged repository literature fixture | original Kontogeorgis literature provenance retained verbatim |
| CR-1 cross association | unchanged explicit repository CR-1 records | original CR-1 literature provenance retained verbatim |
| cubic `kij` | unchanged Folas `-0.055` record, explicitly written to ThermoPack | Folas literature provenance retained verbatim |

The only scientific field values changed relative to Snapshot A are therefore:

```text
METHANOL Tc: 512.64 K -> 512.6 K
WATER    Tc: 647.29 K -> 647.3 K
```

The pinned ThermoPack source uses those component-database temperatures in
`cbeos%single(i)%tc` for the Classic-alpha reduced temperature. Its CPA tuning
API changes `a0`, `b`, association parameters and `c1`, but does not replace this
component `Tc`.

## Important non-default binary parameter

Pinned ThermoPack's default MEOH/H2O cubic interaction record is `kij=-0.09`.
The parity oracle instead explicitly sets `kij=-0.055`, because the purpose is
to compare **the same scientific model and parameter set**, not to compare two
packages' unrelated defaults.

Therefore the parity snapshot must retain:

```text
kij(METHANOL,WATER) = -0.055
source kind          = literature
source identity      = Folas thesis Table 2.1
```

Changing that provenance to `ThermoPack database`, or silently replacing the
value with `-0.09`, is a snapshot-identity failure.

## Machine-readable construction and gate

The test-only constructor is:

```text
tests/flash/cpa_physical_validation/cpa_thermopack_parameter_snapshot.hpp
```

It starts from the existing literature `CpaParameterSet`, copies all records and
provenance, assigns a distinct dataset identity, and replaces only the two `Tc`
records with the values/provenance frozen by the pinned ThermoPack phase-kernel
reference.

The focused identity regression is:

```text
flash.cpa_physical_validation.thermopack_parameter_identity
```

It requires, for both component orders:

- literature and parity `dataset_id` / `revision` remain distinct;
- component identities remain aligned with runtime ordering;
- only `critical_temperature_k` values/provenance differ;
- parity `Tc` is exactly `512.6 K / 647.3 K` and points to the pinned component
  JSON records;
- `a0`, `b`, `c1`, site topology, all self/cross association records and their
  provenance remain byte-for-byte/equality-equivalent to the literature
  snapshot;
- cubic `kij` remains exactly `-0.055` with Folas literature provenance;
- the parity snapshot cannot silently collapse back into the literature
  snapshot or into ThermoPack's default `kij=-0.09` dataset.

## Canonical parity consumers

The named snapshot above is now the **single parameter-construction source** for
both external parity paths:

```text
tests/flash/cpa_physical_validation/cpa_thermopack_delta_test.cpp
tests/flash/cpa_physical_validation/cpa_thermopack_phase_kernel_audit_test.cpp
```

Both call:

```text
cpa_thermopack_snapshot::parameters(false)
```

and require the returned `dataset_id` and `revision` to match the named parity
snapshot before evaluating ThermoPack deltas. Neither test independently copies
CPA records, manufactures a temporary `Tc-only` parameter set, nor carries a
second copy of ThermoPack `Tc` provenance.

The literature snapshot may still be evaluated alongside parity as a diagnostic
to show the consequence of using a scientifically distinct dataset. Such output
is explicitly labeled `LITERATURE_*`; it is not part of the formulation-matched
parity path.

This single-source rule is important before numerical parity thresholds are
frozen: future tolerance evidence must refer to this named snapshot identity,
not to an anonymous counterfactual whose fields could drift independently.

## Scientific interpretation

The two snapshots answer different questions:

- **literature snapshot:** does the repository's documented CPA parameterization
  reproduce the independent Kurihara VLE data within its declared validation
  envelope?
- **ThermoPack parity snapshot:** when both programs are given the same CPA
  formulation and the same effective parameter values, do their numerical
  phase-property and flash results agree?

A future production parameter catalog may expose both datasets explicitly, but
this audit does not promote the parity snapshot into production and does not
replace the literature dataset. Any such promotion requires a separate task and
source review.
