#include <mpmc/flow_discretization_petsc/post_snes_pt_flash_phase_transition_scanner.hpp>
#include <mpmc/flash/sw92_profile_c_pt_flash_backend.hpp>

#include "../../flash/sw92_phase_assigned_pt/physical_sample6.hpp"

#include <petscsys.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace fdp = mpmc::flow_discretization_petsc;
namespace fl = mpmc::flash;
namespace mesh = mpmc::mesh;
namespace sample6 = sw92_profile_c_sample6;
namespace th = mpmc::thermodynamics;

PetscErrorCode resolve_density_from_published_z(
    const fl::PtFlashBackendResult& result,
    void*,
    std::optional<std::vector<double>>* output) {
    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    const auto* accepted = result.accepted_phase_set();
    if (accepted == nullptr) {
        return PETSC_SUCCESS;
    }

    std::vector<double> densities;
    densities.reserve(accepted->phases.size());
    for (const auto& phase : accepted->phases) {
        if (!phase.compressibility_factor ||
            !std::isfinite(*phase.compressibility_factor) ||
            !(*phase.compressibility_factor > 0.0)) {
            return PETSC_SUCCESS;
        }
        const double denominator =
            *phase.compressibility_factor *
            th::Sw92Pure<double>::gas_constant() *
            result.solution.temperature_k;
        const double density =
            result.solution.pressure_pa / denominator;
        if (!std::isfinite(density) || !(density > 0.0)) {
            return PETSC_SUCCESS;
        }
        densities.push_back(density);
    }
    output->emplace(std::move(densities));
    return PETSC_SUCCESS;
}

class IndeterminateBackend final : public fl::PtFlashBackend {
public:
    explicit IndeterminateBackend(std::vector<std::string> component_ids) {
        capability_.backend_id = "test/indeterminate";
        capability_.model_profile = "test/model";
        capability_.algorithm_profile = "test/algorithm";
        capability_.publication_profile = "test/publication";
        capability_.configuration_profile = "test/configuration";
        capability_.dataset_id = "test/dataset";
        capability_.revision = "test/revision";
        capability_.component_ids = std::move(component_ids);
        capability_.supported_phase_counts = {1U, 2U, 3U};
        capability_.transition_capability.edges = {
            {1U, 2U, fl::PtPhaseTransitionSupport::fresh_target_resolve, true},
            {2U, 1U, fl::PtPhaseTransitionSupport::detection_only, true},
            {2U, 3U, fl::PtPhaseTransitionSupport::fresh_target_resolve, true},
            {3U, 2U, fl::PtPhaseTransitionSupport::fresh_target_resolve, true}};
        capability_.performs_initial_stability_search = true;
        capability_.performs_final_phase_set_review = true;
    }

    [[nodiscard]] const fl::PtFlashBackendCapability&
    capability() const noexcept override {
        return capability_;
    }

    [[nodiscard]] fl::PtFlashBackendResult
    solve(const fl::PtFlashRequest& request) override {
        fl::PtFlashBackendResult result;
        result.capability = capability_;
        result.solution.capability.maximum_phase_count =
            capability_.maximum_phase_count();
        result.solution.status = fl::PtPhaseSetStatus::indeterminate;
        result.solution.pressure_pa = request.pressure_pa;
        result.solution.temperature_k = request.temperature_k;
        result.solution.feed = request.feed;
        result.solution.diagnostic =
            "intentional scanner indeterminate fixture";
        result.provider_result_convention = "test/result/v1";
        return result;
    }

private:
    fl::PtFlashBackendCapability capability_;
};

void run_scanner_regression() {
    const auto model = sample6::model();
    fl::Sw92ProfileCPtFlashBackend backend(model);

    fdp::PostSnesPtFlashSourceCellSnapshot3D source;
    source.cell_global = mesh::GlobalEntityId{UINT64_C(7001)};
    source.source_phase_count = 2U;
    source.pressure_pa = 1.0e7;
    source.temperature_k = 350.0;
    source.component_ids = sample6::normal_order();
    source.overall_composition = sample6::feed();

    std::optional<fdp::PostSnesPhaseTransitionProposal3D> proposal;
    fdp::PostSnesPhaseTransitionScanStatus3D status =
        fdp::PostSnesPhaseTransitionScanStatus3D::indeterminate;

    PetscErrorCode error =
        fdp::scan_post_snes_pt_flash_source_cell_3d(
            source,
            backend,
            {&resolve_density_from_published_z, nullptr},
            &proposal,
            &status);
    if (error != PETSC_SUCCESS ||
        status != fdp::PostSnesPhaseTransitionScanStatus3D::complete ||
        !proposal.has_value() ||
        proposal->cell_global != source.cell_global ||
        proposal->candidate.source_phase_count != 2U ||
        proposal->candidate.target_phase_count != 3U ||
        proposal->candidate.target_phases.size() != 3U ||
        proposal->candidate.component_ids != source.component_ids ||
        proposal->candidate.evidence_profile.empty()) {
        throw std::runtime_error(
            "real SW92 Sample-6 2->3 PT scanner regression failed");
    }
    for (const auto& phase : proposal->candidate.target_phases) {
        if (!std::isfinite(phase.molar_density_mol_per_m3) ||
            !(phase.molar_density_mol_per_m3 > 0.0)) {
            throw std::runtime_error(
                "PT scanner target density resolver failed");
        }
    }

    source.source_phase_count = 3U;
    proposal.reset();
    status = fdp::PostSnesPhaseTransitionScanStatus3D::indeterminate;
    error =
        fdp::scan_post_snes_pt_flash_source_cell_3d(
            source,
            backend,
            {&resolve_density_from_published_z, nullptr},
            &proposal,
            &status);
    if (error != PETSC_SUCCESS ||
        status != fdp::PostSnesPhaseTransitionScanStatus3D::complete ||
        proposal.has_value()) {
        throw std::runtime_error(
            "same-cardinality accepted PT result generated a false transition");
    }

    source.source_phase_count = 2U;
    IndeterminateBackend indeterminate_backend{source.component_ids};
    proposal.reset();
    status = fdp::PostSnesPhaseTransitionScanStatus3D::complete;
    error =
        fdp::scan_post_snes_pt_flash_source_cell_3d(
            source,
            indeterminate_backend,
            {&resolve_density_from_published_z, nullptr},
            &proposal,
            &status);
    if (error != PETSC_SUCCESS ||
        status != fdp::PostSnesPhaseTransitionScanStatus3D::indeterminate ||
        proposal.has_value()) {
        throw std::runtime_error(
            "indeterminate PT result was incorrectly treated as stable");
    }
}

} // namespace

void post_snes_pt_flash_phase_transition_scanner_test() {
    int rank = -1;
    if (MPI_Comm_rank(PETSC_COMM_WORLD, &rank) != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI_Comm_rank failed in PT scanner regression");
    }

    int local_success = 1;
    if (rank == 0) {
        try {
            run_scanner_regression();
        } catch (...) {
            local_success = 0;
        }
    }

    int global_success = 0;
    if (MPI_Allreduce(
            &local_success,
            &global_success,
            1,
            MPI_INT,
            MPI_MIN,
            PETSC_COMM_WORLD) != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI_Allreduce failed in PT scanner regression");
    }
    if (global_success == 0) {
        throw std::runtime_error(
            "post-SNES PT flash phase-transition scanner regression failed");
    }
}
