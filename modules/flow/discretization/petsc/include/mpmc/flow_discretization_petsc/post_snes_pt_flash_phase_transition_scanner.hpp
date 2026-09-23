#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_POST_SNES_PT_FLASH_PHASE_TRANSITION_SCANNER_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_POST_SNES_PT_FLASH_PHASE_TRANSITION_SCANNER_HPP

#include <mpmc/flow/phase_set_transition_flash_adapter.hpp>
#include <mpmc/flow_discretization_petsc/post_snes_phase_transition_controller.hpp>

#include <petscvec.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    post_snes_pt_flash_phase_transition_scanner_convention =
        "flow_discretization_petsc/post-snes-pt-flash-phase-transition-scanner/v1";

struct PostSnesPtFlashSourceCellSnapshot3D {
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    std::size_t source_phase_count{};
    double pressure_pa{};
    double temperature_k{};
    std::vector<std::string> component_ids;
    std::vector<double> overall_composition;
};

using PostSnesPtFlashTargetPhaseMolarDensityResolver3D =
    PetscErrorCode (*)(
        const mpmc::flash::PtFlashBackendResult& result,
        void* user_context,
        std::optional<std::vector<double>>* densities);

struct PostSnesPtFlashTargetPhaseMolarDensityBinding3D {
    PostSnesPtFlashTargetPhaseMolarDensityResolver3D resolver{};
    void* user_context{};
};

struct PostSnesPtFlashPhaseTransitionScannerContext3D {
    mpmc::flash::PtFlashBackend* backend{};
    PostSnesPtFlashTargetPhaseMolarDensityBinding3D target_density;
};

namespace post_snes_pt_flash_scanner_detail {

[[nodiscard]] inline bool near_roundoff(double first, double second) {
    if (!std::isfinite(first) || !std::isfinite(second)) {
        return false;
    }
    const double scale = std::max({1.0, std::abs(first), std::abs(second)});
    return std::abs(first - second) <=
        8192.0 * std::numeric_limits<double>::epsilon() * scale;
}

[[nodiscard]] inline bool same_feed(
    std::span<const double> first,
    std::span<const double> second) {
    if (first.size() != second.size()) {
        return false;
    }
    for (std::size_t component = 0U; component < first.size(); ++component) {
        if (!near_roundoff(first[component], second[component])) {
            return false;
        }
    }
    return true;
}

inline void validate_source_snapshot(
    const PostSnesPtFlashSourceCellSnapshot3D& source) {
    if (source.source_phase_count == 0U ||
        source.source_phase_count > mpmc::flow::fixed_three_phase_count ||
        !std::isfinite(source.pressure_pa) || !(source.pressure_pa > 0.0) ||
        !std::isfinite(source.temperature_k) || !(source.temperature_k > 0.0) ||
        source.component_ids.size() < 2U ||
        source.overall_composition.size() != source.component_ids.size()) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: malformed post-SNES PT flash source cell");
    }

    long double sum = 0.0L;
    for (std::size_t component = 0U;
         component < source.component_ids.size();
         ++component) {
        if (source.component_ids[component].empty() ||
            !std::isfinite(source.overall_composition[component]) ||
            !(source.overall_composition[component] > 0.0)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: PT flash source composition must remain on strict positive support");
        }
        for (std::size_t previous = 0U; previous < component; ++previous) {
            if (source.component_ids[previous] == source.component_ids[component]) {
                throw std::invalid_argument(
                    "mpmc::flow_discretization_petsc: duplicate PT flash source component identity");
            }
        }
        sum += static_cast<long double>(source.overall_composition[component]);
    }
    const long double tolerance =
        4096.0L *
        static_cast<long double>(std::numeric_limits<double>::epsilon()) *
        static_cast<long double>(source.component_ids.size());
    if (!std::isfinite(sum) || std::abs(sum - 1.0L) > tolerance) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: PT flash source composition is not normalized");
    }
}

[[nodiscard]] inline mpmc::flow::PoreVolumeComponentAccumulationSnapshot3P
current_component_accumulation(
    const MixedCardinalityPhysicalCurrentCellLinearization3D& current,
    double porosity) {
    return std::visit(
        [porosity](const auto& typed)
            -> mpmc::flow::PoreVolumeComponentAccumulationSnapshot3P {
            using Typed = std::decay_t<decltype(typed)>;
            if constexpr (
                std::is_same_v<Typed, SinglePhaseCurrentCellLinearization3D>) {
                return mpmc::flow::build_single_phase_component_accumulation(
                    typed.state, porosity);
            } else if constexpr (
                std::is_same_v<Typed, TwoPhaseCurrentCellLinearization3D>) {
                return mpmc::flow::build_two_phase_component_accumulation(
                    typed.state, porosity);
            } else {
                return mpmc::flow::build_pore_volume_component_accumulation(
                    typed.state, porosity);
            }
        },
        current);
}

[[nodiscard]] inline bool has_common_pt_pressure(
    const MixedCardinalityPhysicalCurrentCellLinearization3D& current) {
    return std::visit(
        [](const auto& typed) {
            using Typed = std::decay_t<decltype(typed)>;
            if constexpr (
                std::is_same_v<Typed, FixedThreePhaseCurrentCellLinearization3D>) {
                const double reference = typed.state.reference_pressure_pa();
                for (double pressure :
                     typed.saturation_constitutive.phase_pressure_pa) {
                    if (!near_roundoff(pressure, reference)) {
                        return false;
                    }
                }
            }
            return true;
        },
        current);
}

[[nodiscard]] inline PostSnesPtFlashSourceCellSnapshot3D make_source_snapshot(
    mpmc::mesh::GlobalEntityId cell_global,
    const MixedCardinalityPhysicalCurrentCellLinearization3D& current,
    double porosity) {
    const auto accumulation = current_component_accumulation(current, porosity);

    long double total = 0.0L;
    for (double amount :
         accumulation.component_accumulation_mol_per_bulk_m3) {
        if (!std::isfinite(amount) || !(amount > 0.0)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: current component inventory cannot define a PT feed");
        }
        total += static_cast<long double>(amount);
    }
    if (!std::isfinite(total) || !(total > 0.0L)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: invalid current total component inventory");
    }

    PostSnesPtFlashSourceCellSnapshot3D result;
    result.cell_global = cell_global;
    result.component_ids = accumulation.component_ids;
    result.overall_composition.resize(result.component_ids.size());
    for (std::size_t component = 0U;
         component < result.overall_composition.size();
         ++component) {
        result.overall_composition[component] =
            static_cast<double>(
                static_cast<long double>(
                    accumulation.component_accumulation_mol_per_bulk_m3[component]) /
                total);
    }

    std::visit(
        [&](const auto& typed) {
            result.source_phase_count = typed.state.layout().phase_count();
            result.pressure_pa = typed.state.reference_pressure_pa();
            result.temperature_k = typed.state.temperature_k();
        },
        current);
    validate_source_snapshot(result);
    return result;
}

} // namespace post_snes_pt_flash_scanner_detail

[[nodiscard]] inline PetscErrorCode scan_post_snes_pt_flash_source_cell_3d(
    const PostSnesPtFlashSourceCellSnapshot3D& source,
    mpmc::flash::PtFlashBackend& backend,
    PostSnesPtFlashTargetPhaseMolarDensityBinding3D target_density,
    std::optional<PostSnesPhaseTransitionProposal3D>* proposal,
    PostSnesPhaseTransitionScanStatus3D* scan_status) {
    using namespace post_snes_pt_flash_scanner_detail;

    if (proposal == nullptr || scan_status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    proposal->reset();
    *scan_status = PostSnesPhaseTransitionScanStatus3D::complete;

    try {
        validate_source_snapshot(source);
        const auto& capability = backend.capability();
        if (!capability.structurally_valid() ||
            capability.component_ids != source.component_ids ||
            !capability.supports_phase_count(source.source_phase_count)) {
            return PETSC_ERR_ARG_INCOMP;
        }

        mpmc::flash::PtFlashRequest request{
            source.pressure_pa,
            source.temperature_k,
            source.overall_composition};
        auto result = backend.solve(request);

        if (!result.structurally_valid() ||
            result.capability.backend_id != capability.backend_id ||
            result.capability.model_profile != capability.model_profile ||
            result.capability.configuration_profile !=
                capability.configuration_profile ||
            result.capability.dataset_id != capability.dataset_id ||
            result.capability.revision != capability.revision ||
            result.capability.component_ids != capability.component_ids ||
            !near_roundoff(result.solution.pressure_pa, request.pressure_pa) ||
            !near_roundoff(result.solution.temperature_k, request.temperature_k) ||
            !same_feed(result.solution.feed, request.feed)) {
            *scan_status = PostSnesPhaseTransitionScanStatus3D::indeterminate;
            return PETSC_SUCCESS;
        }

        const auto* accepted = result.accepted_phase_set();
        if (accepted == nullptr) {
            *scan_status = PostSnesPhaseTransitionScanStatus3D::indeterminate;
            return PETSC_SUCCESS;
        }
        if (accepted->phases.size() == source.source_phase_count) {
            return PETSC_SUCCESS;
        }

        if (target_density.resolver == nullptr) {
            *scan_status = PostSnesPhaseTransitionScanStatus3D::indeterminate;
            return PETSC_SUCCESS;
        }

        std::optional<std::vector<double>> densities;
        const PetscErrorCode density_error =
            target_density.resolver(
                result,
                target_density.user_context,
                &densities);
        if (density_error != PETSC_SUCCESS) {
            return density_error;
        }
        if (!densities.has_value()) {
            *scan_status = PostSnesPhaseTransitionScanStatus3D::indeterminate;
            return PETSC_SUCCESS;
        }

        auto candidate =
            mpmc::flow::make_phase_set_transition_candidate_from_flash(
                source.source_phase_count,
                result,
                *densities);
        if (!candidate.has_value()) {
            *scan_status = PostSnesPhaseTransitionScanStatus3D::indeterminate;
            return PETSC_SUCCESS;
        }
        proposal->emplace(
            PostSnesPhaseTransitionProposal3D{
                source.cell_global,
                std::move(*candidate)});
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        return PETSC_ERR_ARG_INCOMP;
    } catch (const std::out_of_range&) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    } catch (const std::exception&) {
        *scan_status = PostSnesPhaseTransitionScanStatus3D::indeterminate;
        return PETSC_SUCCESS;
    }
}

[[nodiscard]] inline PetscErrorCode
scan_post_snes_pt_flash_phase_transitions_3d(
    const PhaseTransitionRebuiltNaturalVariableSystem3D& system,
    Vec converged_state,
    const VariableCardinalityNaturalVariableSnesSolveReport3D& solve_report,
    void* raw_context,
    PostSnesPhaseTransitionScanStatus3D* scan_status,
    std::vector<PostSnesPhaseTransitionProposal3D>* local_owned_proposals) {
    using namespace post_snes_pt_flash_scanner_detail;

    if (raw_context == nullptr || scan_status == nullptr ||
        local_owned_proposals == nullptr || converged_state == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    local_owned_proposals->clear();
    *scan_status = PostSnesPhaseTransitionScanStatus3D::complete;

    if (static_cast<int>(solve_report.converged_reason) <= 0) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }

    auto* context =
        static_cast<PostSnesPtFlashPhaseTransitionScannerContext3D*>(raw_context);
    if (context->backend == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    const auto& numbering = system.numbering();
    for (const auto& record : numbering.cells()) {
        if (record.owner_rank != numbering.local_rank()) {
            continue;
        }

        std::vector<PetscInt> indices(record.scalar_count);
        std::vector<PetscScalar> petsc_values(record.scalar_count);
        for (std::size_t slot = 0U; slot < record.scalar_count; ++slot) {
            indices[slot] =
                record.petsc_global_scalar_start + static_cast<PetscInt>(slot);
        }
        PetscErrorCode error =
            VecGetValues(
                converged_state,
                static_cast<PetscInt>(indices.size()),
                indices.data(),
                petsc_values.data());
        if (error != PETSC_SUCCESS) {
            return error;
        }

        std::vector<double> natural_variables(record.scalar_count);
        for (std::size_t slot = 0U; slot < record.scalar_count; ++slot) {
            natural_variables[slot] =
                static_cast<double>(petsc_values[slot]);
        }

        std::optional<MixedCardinalityPhysicalCurrentCellLinearization3D> current;
        double porosity = 0.0;
        NaturalVariableSnesEvaluationStatus3D cell_status =
            NaturalVariableSnesEvaluationStatus3D::success;
        error =
            system.evaluate_current_cell_for_phase_transition(
                record.cell,
                natural_variables,
                &current,
                &porosity,
                &cell_status);
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (cell_status != NaturalVariableSnesEvaluationStatus3D::success ||
            !current.has_value()) {
            *scan_status = PostSnesPhaseTransitionScanStatus3D::indeterminate;
            local_owned_proposals->clear();
            return PETSC_SUCCESS;
        }

        if (!has_common_pt_pressure(*current)) {
            *scan_status = PostSnesPhaseTransitionScanStatus3D::indeterminate;
            local_owned_proposals->clear();
            return PETSC_SUCCESS;
        }

        PostSnesPtFlashSourceCellSnapshot3D source;
        try {
            source = make_source_snapshot(
                record.cell_global,
                *current,
                porosity);
        } catch (const std::exception&) {
            *scan_status = PostSnesPhaseTransitionScanStatus3D::indeterminate;
            local_owned_proposals->clear();
            return PETSC_SUCCESS;
        }

        std::optional<PostSnesPhaseTransitionProposal3D> proposal;
        PostSnesPhaseTransitionScanStatus3D cell_scan_status =
            PostSnesPhaseTransitionScanStatus3D::complete;
        error =
            scan_post_snes_pt_flash_source_cell_3d(
                source,
                *context->backend,
                context->target_density,
                &proposal,
                &cell_scan_status);
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (cell_scan_status == PostSnesPhaseTransitionScanStatus3D::indeterminate) {
            *scan_status = PostSnesPhaseTransitionScanStatus3D::indeterminate;
            local_owned_proposals->clear();
            return PETSC_SUCCESS;
        }
        if (proposal.has_value()) {
            local_owned_proposals->push_back(std::move(*proposal));
        }
    }

    return PETSC_SUCCESS;
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_POST_SNES_PT_FLASH_PHASE_TRANSITION_SCANNER_HPP
