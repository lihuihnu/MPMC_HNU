# SW92 PT equilibrium, transition, and flash correctness contract

## Goal lock

**Current goal:** make the Søreide–Whitson (SW92) PT phase-equilibrium publication contract complete and scientifically honest across one-, two-, and three-phase results, including phase-transition evidence, without assigning physical identities to flash phases.

**Completion conditions:**

- the SW92 PT backend continues to solve the existing corrected-original SW92 model with explicit AQ/NA parameter families;
- accepted one-, two-, and three-phase candidates continue to pass the established material-balance, fugacity-equilibrium, finite stability-review, and boundary-re-solve gates already owned by the Profile-C solver chain;
- all supported phase-count transitions are declared with their real capability level, including an explicit detection-only `2 -> 1` edge rather than a fabricated closed one-phase solve;
- the authoritative PT publication exposes only numerical phase fractions, compositions, phase activities/properties, phase count, transition evidence, and provenance/capability information;
- no authoritative flash phase is labelled water, aqueous, non-aqueous, oil, gas, liquid, vapor, `W`, `H`, AQ, NA, or by a cubic-root identity;
- no frontend file is changed.

**This change does not:** alter the SW92 EOS equations, BIP correlations, root solver, TPD/stability equations, split iterations, convergence tolerances, salinity convention, component data, PR/CPA backends, or frontend behavior.

## Literature-backed meaning of the SW92 water branch

The repository uses one precise solver-side meaning for the phrase **SW92 water branch**:

> the equilibrium branch evaluated with the Søreide–Whitson **aqueous (`AQ`) BIP parameter family** which, when it coexists with one or more `NA` branches in the Profile-C orchestration, must also satisfy the existing water-enrichment admissibility guard relative to those coexisting `NA` branches.

This definition is intentionally narrower than a general physical phase label.

Søreide and Whitson introduced the modified Peng–Robinson formulation to predict mutual solubility of brine/hydrocarbon mixtures and fitted **two sets of binary interaction parameters** for the water/brine and hydrocarbon-side solubility problems. The 2026 Burgoyne–Nielsen refresh describes the same architecture as a **dual-flash scheme**: aqueous-phase BIPs `k^AQ` control dissolved-gas content, while non-aqueous BIPs `k^NA` control water content on the non-aqueous side.

Therefore:

1. `SwPhaseFamily::aqueous` and `SwPhaseFamily::nonaqueous` are **EOS parameter-family assignments** used by the solver.
2. An AQ family assignment is not, by itself, a portable phase identity exposed to callers.
3. An NA family assignment is not a gas/oil/liquid/vapor identity. Existing repository audits have already shown that cubic-root morphology and water-richness cannot be inferred safely from `NA` alone.
4. Internal Profile-C `W/H` slots remain orchestration variables only. They are not part of the public phase semantics.
5. The public phase vector is anonymous. Its order is numerical/provenance output, not an identity contract.

### Sources

- I. Søreide and C. H. Whitson, **“Peng-Robinson predictions for hydrocarbons, CO2, N2, and H2S with pure water and NaCl brine,”** *Fluid Phase Equilibria* 77 (1992) 217–240. DOI: <https://doi.org/10.1016/0378-3812(92)85105-H>.
- M. Burgoyne and M. H. Nielsen, **“Refreshed Søreide–Whitson framework for gas solubility in water and brine with extension to hydrogen,”** *Fluid Phase Equilibria* (2026), article 114824. DOI: <https://doi.org/10.1016/j.fluid.2026.114824>.
- Reproducible code/data accompanying the 2026 refresh: <https://github.com/mwburgoyne/SW_Framework_Refresh>.

## Solver/publication boundary

The existing Profile-C implementation is retained as the numerical and provenance gate:

1. fixed-family SW92 phase properties and cubic-root selection;
2. fixed-family stability searches;
3. NA-only one/two-phase route;
4. targeted AQ-family appearance search using a water-enriched start/admissibility guard;
5. AQ/NA two-phase joint equilibrium;
6. additional NA-family witness search;
7. AQ/NA/NA three-phase equilibrium and final multiplicity review;
8. boundary-aware re-solves for supported three-phase disappearance routes;
9. generic owned `PtPhaseSetResult` publication.

`sw92_pt_flash.hpp` is now the authoritative role-neutral publication boundary. It reuses the existing Profile-C integrity projection and returns only the generic `PtPhaseSetResult`. Provider-specific role/family metadata are not propagated through this boundary.

The runtime `Sw92ProfileCPtFlashBackend` follows the same rule:

- `phase_metadata` is empty;
- `phase_metadata_namespace` is empty;
- `morphology_resolved` is false;
- `publication_profile` is `SW92/PT/role-neutral-phase-set/v1`;
- transition evidence reports only phase counts, triggers, capability/resolution state, and solver diagnostics.

The older provider-side Profile-C metadata types remain implementation/provenance machinery for now so the numerical solver chain is not destabilized in the same PR. They are not the authoritative flash output contract.

## Phase-equilibrium acceptance

A published accepted phase set is not accepted from phase count alone. The upstream solver chain remains responsible for the established gates:

| Gate | Required behavior |
| --- | --- |
| Snapshot/provenance | Pressure, temperature, feed, component order, dataset/revision, salinity, model profile, and orchestration profile must remain consistent through projection. |
| Material balance | Phase fractions/compositions must reconstruct the feed within the solver’s declared numerical tolerances. |
| Fugacity equilibrium | The appropriate SW92 family is used internally for each solver branch and the joint equilibrium residual must satisfy the existing acceptance threshold. |
| Phase multiplicity | One/two/three-phase topology must be closed by the corresponding finite solver/review route; a witness alone is not an accepted phase set. |
| Stability review | The existing finite multistart TPD/stability review is preserved. `global_stability_proven` remains `false`. |
| Boundary behavior | A disappearing phase is never silently dropped. Fresh neighbor re-solves are required where the backend advertises them. |
| Publication integrity | Malformed or provenance-inconsistent sources are downgraded to `indeterminate`. |

Finite multistart stability review is useful numerical evidence, not a mathematical proof of global Gibbs-energy optimality. This PR does not upgrade that claim.

## Phase-transition capability

The public capability contract remains explicit rather than optimistic:

| Transition | Support | Meaning |
| --- | --- | --- |
| `1 -> 2` | `fresh_target_resolve` | instability/topology witness is followed by a fresh two-phase solve and review before acceptance |
| `2 -> 1` | `detection_only` | phase disappearance can be detected, but there is no independent closed one-phase target re-solve in this solver path; no one-phase state is fabricated |
| `2 -> 3` | `fresh_target_resolve` | additional-family witness is followed by a fresh three-phase solve and final multiplicity review |
| `3 -> 2` | `fresh_target_resolve` | supported disappearance routes are re-solved in the lower topology and accepted only if locally closed |
| `3 -> 1` | `fresh_target_resolve` | supported direct lower-topology boundary route is re-solved and accepted only if locally closed |

A failed or unavailable target re-solve remains `target_resolve_failed`, `target_resolve_required`, `broader_topology_required`, or `indeterminate` as appropriate. It is not converted into a lower phase count merely to make the transition table look complete.

## Regression scope

This PR intentionally reuses the existing scientific fixtures rather than introducing fabricated data:

- SW92 binary water/gas backend equivalence and capability contract;
- the existing physical Sample-6 three-phase fixture at 10 MPa and 350 K;
- existing three-phase disappearance/boundary neighbor regressions;
- existing SW92 family, stability, phase-assigned two-phase, three-phase closure, sensitivity, and experimental-line regressions selected by the repository CI rules;
- common PR76/CPA/SW92 max-three-phase backend freeze, to ensure the role-neutral generic contract does not regress the other EOS backends.

No frontend test is required because no frontend source is modified.
