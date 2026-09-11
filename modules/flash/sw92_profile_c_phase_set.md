# SW92 Profile-C authoritative PT phase-set publication

## Purpose and audited boundary

This increment is the publication layer for the already validated SW92
Profile-C PT topology graph. It does **not** add another flash algorithm.
Production calculation is performed only by
`solve_sw92_phase_assigned_pt_boundary_aware(...)`; the new public entry point
runs that driver once and then projects its owned, already-converged evidence
into the model-independent `PtPhaseSetResult` contract.

No EOS equation, binary interaction parameter, root-selection rule, equilibrium
equation, tolerance, finite-stability rule, frontend type or service contract is
changed here.

The authoritative word is deliberately scoped: `accepted` means accepted under
the declared Profile-C finite-search/topology contract. It is not a mathematical
proof of the global Gibbs minimum. Therefore
`global_stability_proven == false` remains mandatory.

## Production entry point

```cpp
#include <mpmc/flash/sw92_profile_c_phase_set.hpp>

const auto result = mpmc::flash::solve_sw92_profile_c_pt_phase_set(
    pressure_pa, temperature_k, feed, sw92_model, nacl_molality);

if (const auto* phases = result.solution.accepted_phase_set()) {
    // phases->phases contains the authoritative 1/2/3-phase PT result accepted
    // by the current Profile-C finite topology contract.
}
```

The solver path is exactly:

```text
solve_sw92_profile_c_pt_phase_set
  -> solve_sw92_phase_assigned_pt_boundary_aware
  -> project_sw92_profile_c_pt_phase_set
```

The projector is separately public so tests/downstream adapters can project an
already-owned boundary-aware result without running the flash twice. It performs
no property call or numerical solve.

## Pure projection rule

The projection copies only phase properties already retained by the validated
source chain:

- `W+H` C1 phases: composition, mole phase fraction, `StabilityPhase activity`
  and Z from the converged joint C1 state;
- `W+H0+H1` C2b.1 phases: the same fields from the converged three-phase state;
- no-W two-H state: the accepted fixed-NA `PtSplitState` liquid/vapor numerical
  slots, without interpreting those slots as physical liquid/vapor morphology;
- no-W one-H state: the accepted fixed-NA reference activity and feed
  composition. That source does not retain a selected Z, so the generic
  `compressibility_factor` is intentionally `nullopt` rather than recomputed.

The adapter performs structural/provenance checks only. It does not call the EOS,
select a root, rerun TPD, rerun material balance, normalize composition, change a
tolerance or invent a missing property.

## Publication gate

`PtPhaseSetStatus::accepted` is emitted only for the four boundary-aware
Profile-C locally-closed states:

| Profile-C status | Published count | Physical-role metadata |
| --- | ---: | --- |
| `no_w_single_h_locally_closed` | 1 | `nonaqueous_unclassified` / NA |
| `no_w_two_h_locally_closed` | 2 | both `nonaqueous_unclassified` / NA |
| `w_h_locally_closed` | 2 | W = aqueous/AQ, H = `nonaqueous_unclassified`/NA |
| `w_h0_h1_locally_closed` | 3 | W = aqueous/AQ, H0/H1 = `nonaqueous_unclassified`/NA |

If a three-phase disappearance route is present, publication additionally
requires the complete C1 -> C2a1 -> C2b.1 -> C2b.2 provenance chain and the
fresh boundary re-solve result. A route status alone never authorizes phase
deletion.

`higher_phase_count_or_wrong_candidate`, unresolved topology, numerical failure,
source-chain inconsistency, malformed phase data or missing property provenance
are conservatively published as `PtPhaseSetStatus::indeterminate` with no
accepted phase set. This first authoritative adapter intentionally does not widen
`phase_set_unstable` semantics.

## Model-independent payload versus SW92 identity

`PtCandidatePhase` remains model-independent and does not gain aqueous/liquid/
vapor labels. `Sw92ProfileCPtPhaseSetResult` therefore carries a parallel
`phase_metadata` array containing:

- the Profile-C physical role (`aqueous` or `nonaqueous_unclassified`);
- the SW92 thermodynamic family (`AQ` or `NA`).

This preserves the separation established by the earlier C2a2 audit:
thermodynamic family, cubic-root branch, phase instance and physical morphology
are different concepts. In particular, H0/H1 are not labelled liquid/vapor from
Z ordering or branch index. Their ordering remains representation-only.

`morphology_resolved == false` is retained explicitly.

## Provenance

The SW92-specific wrapper preserves:

- input feed sum and normalized ordered feed;
- pressure and temperature;
- prescribed NaCl molality;
- dataset ID and revision;
- ordered component IDs;
- SW92 model/phase/equilibrium/orchestration conventions;
- boundary re-solve convention;
- publication convention.

The projector checks that acceptance-critical source objects share the same
p/T/feed/molality/dataset/revision/component snapshot and expected Profile-C
conventions. A mismatched source chain cannot become an accepted generic phase
set merely because a top-level status was manually changed.

## Validation

The focused regression covers:

1. traceable dry CO2/H2O binary -> authoritative one-phase no-W H;
2. traceable wet CO2/H2O binary -> authoritative W+H two-phase result;
3. Mortezazadeh-Rasaei Sample-6 -> authoritative W+H0+H1 three-phase result,
   including exact projection of phase fraction/composition/Z/activity and common
   reduced chemical potentials;
4. fresh physical W+H disappearance-neighbor re-solve remains closed;
5. the physical H0+H1 edge retains its incipient-W witness and is not published
   as an authoritative no-W pair;
6. mismatched C2/fixed-family provenance is downgraded to indeterminate;
7. runtime component-order permutation invariance;
8. H0/H1 slot exchange changes representation order only; both remain
   `nonaqueous_unclassified` NA phase instances;
9. public-header self containment.

GitHub-hosted GCC Debug + ASan/UBSan, Clang Release and MSVC Release are selected
for the new target. Because production equations are unchanged, the affected
regression set is limited to the new publication suite, existing boundary-aware
Profile-C PT tests, existing top-level Profile-C PT tests, and the generic
phase-set representation tests. The independent Sample-6 Decimal(80) oracle is
regenerated before those C++ checks.

## Remaining boundary

This publication gate completes the current fixed-molality SW92 Profile-C PT
1/2/3-phase result contract under finite search. It still does not provide:

- a mathematical global-stability proof;
- a validated H0/H1 liquid/vapor or LV/LL morphology classifier;
- salt-inventory conservation (NaCl molality remains a prescribed model
  parameter in this profile);
- SW92 flash sensitivities or a SW92 physics-closure adapter;
- CPA or later mesh/flow/discretization work.

Those are separate increments and must not be inferred from an accepted
Profile-C `PtPhaseSetResult`.
