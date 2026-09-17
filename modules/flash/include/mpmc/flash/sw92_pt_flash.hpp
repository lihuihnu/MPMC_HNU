#ifndef MPMC_FLASH_SW92_PT_FLASH_HPP
#define MPMC_FLASH_SW92_PT_FLASH_HPP

#include <mpmc/flash/sw92_profile_c_phase_set.hpp>

#include <span>
#include <string_view>
#include <utility>

namespace mpmc::flash {

// Public SW92 PT publication contract.  The solver may use phase-family
// assignments internally, but the published phase vector is deliberately
// anonymous: vector position, cubic root branch, AQ/NA parameter family, and
// internal W/H topology slots are not phase identities.
inline constexpr std::string_view sw92_pt_flash_publication_convention =
    "SW92/PT/role-neutral-phase-set/v1";

// Repository definition used by the SW92 orchestration when it needs the
// model's water-bearing branch internally.  "AQ" denotes the Soreide-Whitson
// aqueous BIP parameter family.  A mixed AQ/NA candidate may retain that AQ
// branch only when it is water-enriched relative to the coexisting NA branch
// under the declared numerical guard.  This is a solver-side admissibility
// rule, never an identity attached to a published flash phase.
inline constexpr std::string_view sw92_aq_water_branch_semantics =
    "AQ-parameterized water-enriched equilibrium branch; solver-side only";

// Pure role-neutral projection.  The existing Profile-C projector remains the
// provenance/integrity gate; only its generic PtPhaseSetResult is released.
// No provider role/family sidecar is propagated across this public boundary.
[[nodiscard]] inline PtPhaseSetResult project_sw92_pt_flash_phase_set(
    const Sw92PhaseAssignedBoundaryAwareResult& source) {
    auto provider_projection = project_sw92_profile_c_pt_phase_set(source);
    return std::move(provider_projection.solution);
}

// Production role-neutral PT entry point.  This does not alter EOS evaluation,
// stability searches, split iteration, boundary resolution, or tolerances; it
// only makes the authoritative publication contract explicit.
[[nodiscard]] inline PtPhaseSetResult solve_sw92_pt_flash(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    const thermodynamics::Sw92Phase<double>& model,
    double nacl_molality_mol_per_kg_water,
    Sw92PhaseAssignedPtOptions options = {}) {
    const auto source = solve_sw92_phase_assigned_pt_boundary_aware(
        pressure_pa, temperature_k, feed, model,
        nacl_molality_mol_per_kg_water, std::move(options));
    return project_sw92_pt_flash_phase_set(source);
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_PT_FLASH_HPP
