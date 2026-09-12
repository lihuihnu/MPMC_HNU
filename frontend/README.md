# MPMC_HNU frontend — Profile-C PT shell

This directory contains the first web frontend increment for MPMC_HNU. It is a
React + TypeScript + Vite application focused on the validated SW92 Profile-C
PT topology result contract.

## Scientific boundary

The browser is **not** an EOS or flash implementation. It must not duplicate,
approximate or reinterpret the C++ thermodynamics and phase-equilibrium logic.

The authoritative backend entry is conceptually:

```text
solve_sw92_phase_assigned_pt(p, T, z, model, molality, options)
```

The frontend contract mirrors its topology statuses and phase-instance fields.
In particular:

- aqueous phase instances map to the SW92 AQ family;
- nonaqueous instances remain `nonaqueous_unclassified`;
- the UI must not relabel NA as liquid or vapor;
- locally-closed finite-search candidates are not displayed as globally proven
  phase sets;
- `globalStabilityProven`, `acceptedPhaseSetPublished` and
  `morphologyResolved` are explicitly false in this v1 contract;
- input compositions are validated but never normalized, clipped or repaired by
  the browser.

## Backend connection

No compute transport is fabricated in this increment. `FlashClient` is the
single injectable boundary. The production entry currently uses
`unconfiguredFlashClient`, which disables calculation submission and reports
that the compute service is not configured.

A later service increment should implement this interface against an explicit
versioned MPMC_HNU transport contract. The service remains responsible for
server-side validation, model/parameter capability checks and all C++
thermodynamic calculations.

Tests may inject deterministic response fixtures to validate presentation, but
fixtures must never be shipped as production calculation results.

## Current UI

The initial workspace provides:

- pressure input in MPa, converted to Pa only at the request boundary;
- temperature in K;
- NaCl molality in mol/kg H2O;
- runtime ordered component IDs and overall mole fractions;
- exact backend-aligned topology status display;
- one-, two- and three-phase candidate cards;
- phase fraction, ordered composition, AQ/NA family and optional Z;
- explicit diagnostic and scientific-validity messaging;
- responsive desktop/mobile layout.

## Development

Requires Node 24.x.

```bash
cd frontend
npm install --ignore-scripts --no-audit --no-fund
npm run typecheck
npm test
npm run build
npm run dev
```

Dependencies are pinned to exact versions in `package.json`. The frontend is an
independent web build and is not a prerequisite for building or testing the C++
scientific core.

## Next integration gate

The transport-neutral C++ [model-neutral PT service boundary](../modules/runtime/README.md)
now defines configured-backend discovery, runtime component inventories,
ID-keyed requests, variable accepted phase counts, provenance, and explicit
indeterminate/error outcomes across PR76, SW92 Profile-C, and CPA. This frontend
has not yet been migrated to or connected through that contract.

The next functional increment is a versioned Protobuf/gRPC-Web adapter and a
frontend migration that discovers the selected backend before building the
component form. The adapter must call `mpmc::runtime::PtService`, not the
model-specific Profile-C API, and must preserve non-accepted outcomes without
fabricating phase data. Until that transport exists, the disabled compute state
remains intentional.
