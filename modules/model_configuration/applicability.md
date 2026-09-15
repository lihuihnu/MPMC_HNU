# Model applicability endpoint semantics

`ThermodynamicModelDefinition::applicability` is a provenance-bearing declaration of the state range supported by the supplied parameter dataset. It is not an algebraic domain extension and it is not a promise that states with undeclared bounds are validated.

## Public contract

Temperature and pressure each have an optional lower endpoint and optional upper endpoint. Every declared endpoint is finite, strictly positive SI (`K` or `Pa`). A missing endpoint means **unknown**, not `-infinity` or `+infinity`.

Each endpoint is closed by default. The corresponding `*_exclusive` flag changes only that declared endpoint to open. An exclusive flag without its endpoint is invalid. When both endpoints are present, lower must not exceed upper; equal endpoints are legal only when both are closed, representing the single declared state.

`ModelApplicability::assess_*` has three outcomes:

- `outside_declared_bounds`: any declared endpoint is definitely violated, including equality at an open endpoint;
- `inside_declared_bounds`: both endpoints of that axis are declared and the value satisfies both;
- `unknown`: no declared endpoint is violated, but at least one endpoint is missing.

The combined `assess(T, p)` is outside if either axis is definitely outside, inside only if both axes are inside, and otherwise unknown. Non-finite query values assess as unknown; solve-request validation rejects non-finite or non-positive P/T separately.

## PR76 mapping and solve enforcement

The existing thermodynamics `Applicability` type stores only complete closed intervals. The model-configuration adapter therefore uses a conservative bridge:

- a complete public interval maps to the same native closed numeric envelope;
- one-sided public bounds do **not** fabricate the missing endpoint and map to an absent native interval;
- open endpoints are not moved with `nextafter`, rounded, clipped, or replaced by invented limits.

`Pr76ExecutableModel` supplements the native kernel at the public boundary. It rejects a definite violation that the native closed-interval representation cannot express: one-sided violations and exact equality at an open endpoint. Complete closed-interval violations continue to come from the unchanged native PR76 property path. The rejection remains a `std::domain_error` plus `Pr76SolveRequestError` field location (`temperature_k` or `pressure_pa`).

This layer does not alter PR76 formulas, roots, stability, split, three-phase logic, tolerances, or phase acceptance. States that are merely `unknown` are not relabeled as scientifically validated; they remain executable algebraic states unless another model/domain rule rejects them.

## Versioned wire contract

`ModelConfigurationService` carries the four endpoint scalars with presence and four additive boolean exclusivity flags. `false` preserves the historical closed-endpoint default. Create/describe snapshots round-trip one-sided presence and open/closed flags exactly. The old `PtFlashService v1` is unchanged.

The UI remains deferred. Endpoint editing and user-facing visualization must consume this contract rather than reinterpret missing bounds as unlimited ranges.
