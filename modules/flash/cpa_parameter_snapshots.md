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

### Dataset semantics of the literature `Tc`

The 2008 paper does not present these critical temperatures as detached generic
fluid metadata. Table 1a is explicitly a table of **CPA pure-compound
parameters for the associating compounds used in that work**. In the same row it
lists `Tc`, `b`, `a0`, `c1`, association energy and association volume. The paper
also states that the CPA pure-component parameters are fitted against vapor
pressure and liquid-density data.

Therefore, for repository provenance purposes, `512.64 K` and `647.29 K` belong
to the same literature CPA parameterization identity as the corresponding
`a0/b/c1/epsilon/beta` records. They must not be silently replaced by a generic
critical-property database value while retaining the original literature dataset
identifier.

Primary source:

```text
Kontogeorgis, Folas, Muro Sunè, Roca Leon & Michelsen (2008)
Oil & Gas Science and Technology 63(3), 305-319
DOI: 10.2516/ogst:2008025
Table 1a: CPA pure compound parameters for associating compounds
```

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

### Dataset semantics of the ThermoPack DB `Tc`

At the pinned ThermoPack revision, methanol and water store these values in the
component-level `critical` object:

```text
fluids/Methanol.json: critical.temperature = 512.6 K
fluids/Water.json:    critical.temperature = 647.3 K
```

Those records are marked `ref: "Default"` and contain no `bib_reference`. The CPA
records are separate nested `SRK -> CPA-*` objects with their own CPA parameter
references. In particular, methanol `CPA-1` references Kontogeorgis 2008 while
water `CPA-1` references the water CPA literature chain; the component-level
critical-temperature record is not embedded in those CPA parameter blocks.

The pinned cubic initialization path assigns `cbeos%single(i)%Tc` from the
component record (`comp(i)%p_comp%tc`) unless an explicit critical temperature is
supplied. The CPA pure-parameter setter separately replaces the CPA/cubic `a0`,
`b`, Classic-alpha coefficient and association parameters; it does not redefine
the component critical-temperature record. Consequently ThermoPack's effective
Classic-alpha reduced temperature uses the component-database `Tc`.

Pinned source locators:

```text
thermotools/thermopack@d68c794c7342bfc6938eb424a1fbb88b7780b738
  fluids/Methanol.json
  fluids/Water.json
  src/cbselect.f90
  src/saft_interface.f90
```

This means the ThermoPack values are best described here as **pinned
component-database critical-property metadata consumed by the ThermoPack
Classic-alpha implementation**, not as the `Tc` entries of the repository's
Kontogeorgis-2008 literature parameter snapshot.

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

## Production parameter policy decision

**Decision for the current CPA backend: keep the literature snapshot as the
production/physical-validation parameter identity. Keep the ThermoPack-aligned
`Tc` snapshot test-only as a named parity fixture.**

Rationale:

1. The literature snapshot has a complete, traceable scientific chain and has
   already been validated against independent Kurihara methanol-water VLE data.
2. Its `Tc` values are part of the same published CPA pure-parameter table as
   `a0/b/c1/epsilon/beta`; replacing only `Tc` while retaining the literature
   identity would create an undeclared mixed-provenance parameterization.
3. The ThermoPack component-database `Tc` records serve a different purpose:
   they are package-level critical-property metadata used by ThermoPack's cubic
   initialization and Classic alpha. Their pinned JSON records do not provide a
   source trail demonstrating that `512.6/647.3 K` should replace the published
   literature values in MPMC_HNU's physical-validation dataset.
4. The single-factor audit shows that changing only these two `Tc` values nearly
   eliminates the observed MPMC_HNU/ThermoPack numerical delta. That is strong
   evidence about **software parity provenance**, but it is not evidence that the
   ThermoPack DB values are a scientifically superior production parameter set.
5. External-software parity and independent physical validation answer different
   questions and must retain separate dataset identities.

Repository policy until a separate parameter-catalog review explicitly changes
it:

```text
production / independent physical validation
  -> literature-cpa-water-methanol-cr1-333.15K
  -> Tc(MeOH)=512.64 K, Tc(H2O)=647.29 K

ThermoPack numerical parity only
  -> thermopack-d68c794__meoh-h2o__mpmc-literature-folas-cr1-aligned
  -> Tc(MeOH)=512.6 K, Tc(H2O)=647.3 K
  -> test-side fixture only
```

The parity snapshot must **not** be promoted into the production CPA parameter
catalog merely because it reduces external-software deltas. Such a promotion
would require a separate scientific task that establishes a source-complete
parameter-set identity, validates it against independent physical data, and
states whether the entire CPA parameter set was fitted consistently with the
selected critical properties.

## Scientific interpretation

The two snapshots answer different questions:

- **literature snapshot:** does the repository's documented CPA parameterization
  reproduce the independent Kurihara VLE data within its declared validation
  envelope?
- **ThermoPack parity snapshot:** when both programs are given the same CPA
  formulation and the same effective parameter values, do their numerical
  phase-property and flash results agree?

The production decision above intentionally keeps those questions separate.
No production parameter has been changed by this audit.
