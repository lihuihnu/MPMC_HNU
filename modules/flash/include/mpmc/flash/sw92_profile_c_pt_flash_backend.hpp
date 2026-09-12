#ifndef MPMC_FLASH_SW92_PROFILE_C_PT_FLASH_BACKEND_HPP
#define MPMC_FLASH_SW92_PROFILE_C_PT_FLASH_BACKEND_HPP

#include <mpmc/flash/pt_flash_backend.hpp>
#include <mpmc/flash/sw92_profile_c_phase_set.hpp>

#include <optional>
#include <string>
#include <string_view>
#include <utility>

namespace mpmc::flash {

inline constexpr std::string_view sw92_profile_c_pt_flash_backend_id =
    "SW92/Profile-C/PT/phase-set-backend/v1";
inline constexpr std::string_view sw92_profile_c_pt_flash_backend_configuration_profile =
    "SW92/Profile-C/fixed-molality/backend-configuration/v1";
inline constexpr std::string_view sw92_profile_c_phase_metadata_namespace =
    "SW92/Profile-C/phase-metadata/v1";
inline constexpr std::string_view sw92_profile_c_transition_evidence_profile =
    "SW92/Profile-C/phase-transition-evidence/v1";

struct Sw92ProfileCPtFlashBackendOptions {
    double nacl_molality_mol_per_kg_water{};
    Sw92PhaseAssignedPtOptions flash;
};

class Sw92ProfileCPtFlashBackend final : public PtFlashBackend {
public:
    Sw92ProfileCPtFlashBackend(
        const thermodynamics::Sw92Phase<double>& model,
        Sw92ProfileCPtFlashBackendOptions options = {})
        : model_(model), options_(std::move(options)),
          capability_(build_capability(model_, options_)) {}

    [[nodiscard]] const PtFlashBackendCapability& capability()
        const noexcept override {
        return capability_;
    }

    [[nodiscard]] PtFlashBackendResult solve(
        const PtFlashRequest& request) override {
        // Run the established boundary-aware Profile-C driver exactly once and
        // retain its owned topology chain for transition evidence. Publication is
        // still the existing pure projection used by solve_sw92_profile_c_pt_phase_set.
        const auto source = solve_sw92_phase_assigned_pt_boundary_aware(
            request.pressure_pa, request.temperature_k, request.feed, model_,
            options_.nacl_molality_mol_per_kg_water, options_.flash);
        const auto published = project_sw92_profile_c_pt_phase_set(source);

        PtFlashBackendResult result;
        result.capability = capability_;
        result.solution = published.solution;
        result.transition_report = build_transition_report(source);
        result.provider_result_convention = published.publication_convention;
        result.morphology_resolved =
            Sw92ProfileCPtPhaseSetResult::morphology_resolved;
        result.phase_metadata.reserve(published.phase_metadata.size());
        for (const auto& metadata : published.phase_metadata) {
            result.phase_metadata.push_back({
                physical_role_id(metadata.physical_role),
                thermodynamic_family_id(metadata.thermodynamic_family)});
        }

        if (published.dataset_id != capability_.dataset_id ||
            published.revision != capability_.revision ||
            published.component_ids != capability_.component_ids ||
            published.model_profile != capability_.model_profile ||
            published.orchestration_convention != capability_.algorithm_profile ||
            published.publication_convention != capability_.publication_profile ||
            published.nacl_molality_mol_per_kg_water !=
                options_.nacl_molality_mol_per_kg_water) {
            reject_adapter_result(
                result, "SW92 Profile-C backend: provider/model provenance mismatch");
            return result;
        }

        if (!result.structurally_valid()) {
            reject_adapter_result(
                result,
                "SW92 Profile-C backend: generic result integrity guard failed");
        }
        return result;
    }

private:
    [[nodiscard]] static PtFlashBackendCapability build_capability(
        const thermodynamics::Sw92Phase<double>& model,
        const Sw92ProfileCPtFlashBackendOptions& options) {
        PtFlashBackendCapability capability;
        capability.backend_id =
            std::string(sw92_profile_c_pt_flash_backend_id);
        capability.model_profile =
            std::string(thermodynamics::sw92_corrected_profile);
        capability.algorithm_profile =
            std::string(sw92_phase_assigned_pt_convention);
        capability.publication_profile =
            std::string(sw92_profile_c_phase_set_publication_convention);
        capability.configuration_profile =
            std::string(sw92_profile_c_pt_flash_backend_configuration_profile);
        capability.scalar_settings.push_back({
            "nacl_molality_mol_per_kg_water",
            options.nacl_molality_mol_per_kg_water,
            "mol/kg_H2O"});
        const auto& parameters = model.parameters();
        capability.dataset_id = parameters.dataset_id();
        capability.revision = parameters.revision();
        for (const auto& component : parameters.components().items()) {
            capability.component_ids.push_back(component.id);
        }
        capability.supported_phase_counts = {1U, 2U, 3U};
        capability.transition_capability.edges = {
            {1U, 2U, PtPhaseTransitionSupport::fresh_target_resolve, true},
            {2U, 1U, PtPhaseTransitionSupport::detection_only, true},
            {2U, 3U, PtPhaseTransitionSupport::fresh_target_resolve, true},
            {3U, 2U, PtPhaseTransitionSupport::fresh_target_resolve, true},
            {3U, 1U, PtPhaseTransitionSupport::fresh_target_resolve, true}};
        capability.performs_initial_stability_search = true;
        capability.performs_final_phase_set_review = true;
        capability.performs_boundary_neighbor_resolve = true;
        capability.global_stability_proven = false;
        capability.phase_metadata_namespace =
            std::string(sw92_profile_c_phase_metadata_namespace);
        return capability;
    }

    [[nodiscard]] static PtPhaseTransitionReport build_transition_report(
        const Sw92PhaseAssignedBoundaryAwareResult& source) {
        PtPhaseTransitionReport report;
        const auto add = [&](
            std::size_t from, std::optional<std::size_t> to,
            PtPhaseTransitionTrigger trigger,
            PtPhaseTransitionResolution resolution,
            bool fresh_attempted, bool target_closed,
            std::string diagnostic) {
            report.evidence.push_back({
                from, to, trigger, resolution, fresh_attempted, target_closed,
                std::string(sw92_profile_c_transition_evidence_profile),
                std::move(diagnostic)});
        };

        if (source.boundary) {
            switch (source.boundary->status) {
            case Sw92PhaseAssignedBoundaryStatus::resolved_to_w_h:
            case Sw92PhaseAssignedBoundaryStatus::resolved_to_no_w_two_h:
                add(3U, 2U, PtPhaseTransitionTrigger::provider_boundary_route,
                    PtPhaseTransitionResolution::accepted_target,
                    source.boundary->re_solve_attempted, true,
                    source.boundary->diagnostic);
                return report;
            case Sw92PhaseAssignedBoundaryStatus::resolved_to_no_w_single_h:
                add(3U, 1U, PtPhaseTransitionTrigger::provider_boundary_route,
                    PtPhaseTransitionResolution::accepted_target,
                    source.boundary->re_solve_attempted, true,
                    source.boundary->diagnostic);
                return report;
            case Sw92PhaseAssignedBoundaryStatus::neighbor_topology_not_closed:
            case Sw92PhaseAssignedBoundaryStatus::higher_phase_count_or_wrong_candidate:
                add(3U, std::nullopt,
                    PtPhaseTransitionTrigger::phase_disappearance,
                    PtPhaseTransitionResolution::target_resolve_failed,
                    source.boundary->re_solve_attempted, false,
                    source.boundary->diagnostic);
                return report;
            case Sw92PhaseAssignedBoundaryStatus::source_chain_inconsistent:
            case Sw92PhaseAssignedBoundaryStatus::numerical_indeterminate:
                add(3U, std::nullopt,
                    PtPhaseTransitionTrigger::phase_disappearance,
                    PtPhaseTransitionResolution::indeterminate,
                    source.boundary->re_solve_attempted, false,
                    source.boundary->diagnostic);
                return report;
            case Sw92PhaseAssignedBoundaryStatus::no_boundary_route:
                break;
            }
        }

        switch (source.status) {
        case Sw92PhaseAssignedPtStatus::no_w_single_h_locally_closed:
            break;
        case Sw92PhaseAssignedPtStatus::no_w_two_h_locally_closed:
            add(1U, 2U, PtPhaseTransitionTrigger::initial_stability_witness,
                PtPhaseTransitionResolution::accepted_target,
                true, true,
                "fixed-NA feed instability was followed by a fresh two-phase solve and final review");
            break;
        case Sw92PhaseAssignedPtStatus::w_h_locally_closed:
            add(1U, 2U, PtPhaseTransitionTrigger::provider_topology_witness,
                PtPhaseTransitionResolution::accepted_target,
                true, true,
                "an aqueous-appearance witness was followed by a fresh W(AQ)+H(NA) joint solve and H-side review");
            break;
        case Sw92PhaseAssignedPtStatus::w_h0_h1_locally_closed:
            add(2U, 3U, PtPhaseTransitionTrigger::provider_topology_witness,
                PtPhaseTransitionResolution::accepted_target,
                true, true,
                "an additional-H witness was followed by a fresh W(AQ)+H0(NA)+H1(NA) solve and final multiplicity review");
            break;
        case Sw92PhaseAssignedPtStatus::higher_phase_count_or_wrong_candidate:
            add(source.phases.empty() ? 2U : source.phases.size(), std::nullopt,
                PtPhaseTransitionTrigger::final_phase_set_instability,
                PtPhaseTransitionResolution::broader_topology_required,
                false, false,
                source.diagnostic);
            break;
        case Sw92PhaseAssignedPtStatus::topology_unresolved:
            if (source.base.no_w.status ==
                Sw92PhaseAssignedNoWStatus::phase_disappearance_unresolved) {
                add(2U, 1U, PtPhaseTransitionTrigger::phase_disappearance,
                    PtPhaseTransitionResolution::target_resolve_required,
                    false, false,
                    source.base.no_w.diagnostic);
            } else if (source.base.c2a1 &&
                source.base.c2a1->status ==
                    Sw92PhaseAssignedHSideWitnessStatus::
                        additional_nonaqueous_phase_witness_found) {
                add(2U, 3U, PtPhaseTransitionTrigger::provider_topology_witness,
                    PtPhaseTransitionResolution::candidate_not_accepted,
                    source.base.c2b1.has_value(), false,
                    source.diagnostic);
            }
            break;
        case Sw92PhaseAssignedPtStatus::numerical_indeterminate:
            if (source.base.c2a1 &&
                source.base.c2a1->status ==
                    Sw92PhaseAssignedHSideWitnessStatus::
                        additional_nonaqueous_phase_witness_found) {
                add(2U, 3U, PtPhaseTransitionTrigger::provider_topology_witness,
                    PtPhaseTransitionResolution::indeterminate,
                    source.base.c2b1.has_value(), false,
                    source.diagnostic);
            }
            break;
        }
        return report;
    }

    [[nodiscard]] static std::string physical_role_id(
        Sw92PhaseAssignedPtPhysicalRole role) {
        switch (role) {
        case Sw92PhaseAssignedPtPhysicalRole::aqueous:
            return "aqueous";
        case Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified:
            return "nonaqueous_unclassified";
        }
        return "unknown";
    }

    [[nodiscard]] static std::string thermodynamic_family_id(
        thermodynamics::SwPhaseFamily family) {
        switch (family) {
        case thermodynamics::SwPhaseFamily::aqueous:
            return "aqueous";
        case thermodynamics::SwPhaseFamily::nonaqueous:
            return "nonaqueous";
        }
        return "unknown";
    }

    static void reject_adapter_result(
        PtFlashBackendResult& result, std::string diagnostic) {
        result.solution.status = PtPhaseSetStatus::indeterminate;
        result.solution.candidate_phase_set.reset();
        result.solution.global_stability_proven = false;
        result.phase_metadata.clear();
        result.transition_report.evidence.clear();
        result.solution.diagnostic = std::move(diagnostic);
    }

    const thermodynamics::Sw92Phase<double>& model_;
    Sw92ProfileCPtFlashBackendOptions options_;
    PtFlashBackendCapability capability_;
};

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_PROFILE_C_PT_FLASH_BACKEND_HPP
