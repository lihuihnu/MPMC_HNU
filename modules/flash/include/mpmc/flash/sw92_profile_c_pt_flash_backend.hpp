#ifndef MPMC_FLASH_SW92_PROFILE_C_PT_FLASH_BACKEND_HPP
#define MPMC_FLASH_SW92_PROFILE_C_PT_FLASH_BACKEND_HPP

#include <mpmc/flash/pt_flash_backend.hpp>
#include <mpmc/flash/sw92_profile_c_phase_set.hpp>

#include <string>
#include <string_view>
#include <utility>

namespace mpmc::flash {

inline constexpr std::string_view sw92_profile_c_pt_flash_backend_id =
    "SW92/Profile-C/PT/phase-set-backend/v1";
inline constexpr std::string_view sw92_profile_c_phase_metadata_namespace =
    "SW92/Profile-C/phase-metadata/v1";

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
          capability_(build_capability(model_)) {}

    [[nodiscard]] const PtFlashBackendCapability& capability()
        const noexcept override {
        return capability_;
    }

    [[nodiscard]] PtFlashBackendResult solve(
        const PtFlashRequest& request) override {
        const auto source = solve_sw92_profile_c_pt_phase_set(
            request.pressure_pa, request.temperature_k, request.feed, model_,
            options_.nacl_molality_mol_per_kg_water, options_.flash);

        PtFlashBackendResult result;
        result.capability = capability_;
        result.solution = source.solution;
        result.provider_result_convention = source.publication_convention;
        result.morphology_resolved =
            Sw92ProfileCPtPhaseSetResult::morphology_resolved;
        result.phase_metadata.reserve(source.phase_metadata.size());
        for (const auto& metadata : source.phase_metadata) {
            result.phase_metadata.push_back({
                physical_role_id(metadata.physical_role),
                thermodynamic_family_id(metadata.thermodynamic_family)});
        }

        if (source.dataset_id != capability_.dataset_id ||
            source.revision != capability_.revision ||
            source.component_ids != capability_.component_ids ||
            source.model_profile != capability_.model_profile ||
            source.orchestration_convention != capability_.algorithm_profile ||
            source.publication_convention != capability_.publication_profile ||
            source.nacl_molality_mol_per_kg_water !=
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
        const thermodynamics::Sw92Phase<double>& model) {
        PtFlashBackendCapability capability;
        capability.backend_id =
            std::string(sw92_profile_c_pt_flash_backend_id);
        capability.model_profile =
            std::string(thermodynamics::sw92_corrected_profile);
        capability.algorithm_profile =
            std::string(sw92_phase_assigned_pt_convention);
        capability.publication_profile =
            std::string(sw92_profile_c_phase_set_publication_convention);
        const auto& parameters = model.parameters();
        capability.dataset_id = parameters.dataset_id();
        capability.revision = parameters.revision();
        for (const auto& component : parameters.components().items()) {
            capability.component_ids.push_back(component.id);
        }
        capability.supported_phase_counts = {1U, 2U, 3U};
        capability.performs_initial_stability_search = true;
        capability.performs_final_phase_set_review = true;
        capability.performs_boundary_neighbor_resolve = true;
        capability.global_stability_proven = false;
        capability.phase_metadata_namespace =
            std::string(sw92_profile_c_phase_metadata_namespace);
        return capability;
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
        result.solution.diagnostic = std::move(diagnostic);
    }

    const thermodynamics::Sw92Phase<double>& model_;
    Sw92ProfileCPtFlashBackendOptions options_;
    PtFlashBackendCapability capability_;
};

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_PROFILE_C_PT_FLASH_BACKEND_HPP
