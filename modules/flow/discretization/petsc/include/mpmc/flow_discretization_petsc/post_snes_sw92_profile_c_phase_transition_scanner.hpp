#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_POST_SNES_SW92_PROFILE_C_PHASE_TRANSITION_SCANNER_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_POST_SNES_SW92_PROFILE_C_PHASE_TRANSITION_SCANNER_HPP

#include <mpmc/flash/sw92_profile_c_phase_set.hpp>
#include <mpmc/flash/sw92_profile_c_pt_flash_backend.hpp>
#include <mpmc/flow/phase_set_transition_flash_adapter.hpp>
#include <mpmc/flow_discretization_petsc/post_snes_pt_flash_phase_transition_scanner.hpp>
#include <mpmc/thermodynamics/selected_phase_density.hpp>

#include <petscvec.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    post_snes_sw92_profile_c_phase_transition_scanner_convention =
        "flow_discretization_petsc/post-snes-sw92-profile-c-authoritative-target/v1";

struct Sw92AuthoritativeTargetPhaseMaterialization3D {
    double mole_phase_fraction{};
    std::vector<double> composition;
    mpmc::flash::Sw92ProfileCPhaseMetadata metadata;
    mpmc::thermodynamics::Sw92SelectedPhase<double>
        selection;
    double molar_density_mol_per_m3{};
};

struct Sw92AuthoritativeTargetMaterialization3D {
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    double nacl_molality_mol_per_kg_water{};
    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string publication_convention;
    std::string transition_evidence_profile;
    std::string transition_diagnostic;
    std::vector<
        Sw92AuthoritativeTargetPhaseMaterialization3D>
        phases;
    mpmc::flow::PhaseSetTransitionProjection
        projection;
};

struct PostSnesSw92ProfileCTransitionScanCellResult3D {
    PostSnesPhaseTransitionProposal3D proposal;
    Sw92AuthoritativeTargetMaterialization3D target;
};

struct PostSnesSw92ProfileCPhaseTransitionScannerContext3D {
    const mpmc::thermodynamics::Sw92Phase<double>*
        model{};
    mpmc::flash::Sw92ProfileCPtFlashBackendOptions
        options;
    std::vector<
        Sw92AuthoritativeTargetMaterialization3D>
        local_owned_targets;

    [[nodiscard]]
    const Sw92AuthoritativeTargetMaterialization3D*
    find_target(
        mpmc::mesh::GlobalEntityId cell_global) const noexcept {
        const auto found =
            std::find_if(
                local_owned_targets.begin(),
                local_owned_targets.end(),
                [&](const auto& target) {
                    return target.cell_global ==
                        cell_global;
                });
        return found ==
                local_owned_targets.end()
            ? nullptr
            : &*found;
    }
};

namespace post_snes_sw92_profile_c_scanner_detail {

[[nodiscard]] inline bool near_roundoff(
    double first,
    double second) {
    return post_snes_pt_flash_scanner_detail::
        near_roundoff(first, second);
}

[[nodiscard]] inline bool same_vector(
    std::span<const double> first,
    std::span<const double> second) {
    if (first.size() != second.size()) {
        return false;
    }
    for (std::size_t i = 0U;
         i < first.size();
         ++i) {
        if (!near_roundoff(
                first[i],
                second[i])) {
            return false;
        }
    }
    return true;
}

inline void validate_model_source_identity(
    const PostSnesPtFlashSourceCellSnapshot3D& source,
    const mpmc::thermodynamics::Sw92Phase<double>&
        model,
    const mpmc::flash::
        Sw92ProfileCPtFlashBackendOptions& options) {
    post_snes_pt_flash_scanner_detail::
        validate_source_snapshot(source);
    if (!std::isfinite(
            options
                .nacl_molality_mol_per_kg_water) ||
        options
                .nacl_molality_mol_per_kg_water <
            0.0 ||
        model.size() !=
            source.component_ids.size()) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: SW92 scanner model/source cardinality or molality mismatch");
    }

    const auto& parameters =
        model.parameters();
    for (std::size_t component = 0U;
         component < source.component_ids.size();
         ++component) {
        if (source.component_ids[component] !=
            parameters.components()
                .at(component)
                .id) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: SW92 scanner component identity/order mismatch");
        }
    }
}

[[nodiscard]] inline
mpmc::flash::PtFlashBackendResult
make_generic_result_from_owned_source(
    const mpmc::flash::
        Sw92PhaseAssignedBoundaryAwareResult& source,
    const mpmc::thermodynamics::Sw92Phase<double>&
        model,
    const mpmc::flash::
        Sw92ProfileCPtFlashBackendOptions& options) {
    mpmc::flash::Sw92ProfileCPtFlashBackend
        capability_source{
            model,
            options};

    mpmc::flash::PtFlashBackendResult result;
    result.capability =
        capability_source.capability();
    result.solution =
        mpmc::flash::
            project_sw92_pt_flash_phase_set(
                source);
    result.transition_report =
        mpmc::flash::
            project_sw92_profile_c_transition_report(
                source);
    result.provider_result_convention =
        std::string{
            mpmc::flash::
                sw92_pt_flash_publication_convention};
    result.morphology_resolved = false;
    result.phase_metadata.clear();
    return result;
}

inline void validate_authoritative_projection_pair(
    const mpmc::flash::PtFlashBackendResult&
        generic,
    const mpmc::flash::
        Sw92ProfileCPtPhaseSetResult&
            authoritative,
    const PostSnesPtFlashSourceCellSnapshot3D&
        source,
    const mpmc::thermodynamics::Sw92Phase<double>&
        model,
    const mpmc::flash::
        Sw92ProfileCPtFlashBackendOptions& options) {
    if (!generic.structurally_valid() ||
        authoritative.solution.status !=
            mpmc::flash::PtPhaseSetStatus::
                accepted ||
        !authoritative
             .accepted_phase_set_published() ||
        authoritative.dataset_id !=
            model.parameters().dataset_id() ||
        authoritative.revision !=
            model.parameters().revision() ||
        authoritative.component_ids !=
            source.component_ids ||
        authoritative
                .nacl_molality_mol_per_kg_water !=
            options
                .nacl_molality_mol_per_kg_water ||
        !near_roundoff(
            authoritative.solution.pressure_pa,
            source.pressure_pa) ||
        !near_roundoff(
            authoritative.solution.temperature_k,
            source.temperature_k) ||
        !same_vector(
            authoritative.solution.feed,
            source.overall_composition)) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: SW92 authoritative/generic phase-set provenance mismatch");
    }

    const auto* generic_set =
        generic.accepted_phase_set();
    const auto* authoritative_set =
        authoritative.solution
            .accepted_phase_set();
    if (generic_set == nullptr ||
        authoritative_set == nullptr ||
        generic_set->phases.size() !=
            authoritative_set->phases.size() ||
        authoritative.phase_metadata.size() !=
            authoritative_set
                ->phases.size()) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: SW92 authoritative/generic accepted phase-set shape mismatch");
    }

    for (std::size_t phase = 0U;
         phase < generic_set->phases.size();
         ++phase) {
        const auto& first =
            generic_set->phases[phase];
        const auto& second =
            authoritative_set->phases[phase];
        if (!near_roundoff(
                first.mole_phase_fraction,
                second.mole_phase_fraction) ||
            !same_vector(
                first.composition,
                second.composition) ||
            first.activity.branch !=
                second.activity.branch ||
            first.activity.smooth !=
                second.activity.smooth) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: SW92 generic publication changed authoritative phase slot/branch provenance");
        }
    }
}

[[nodiscard]] inline
mpmc::flow::PoreVolumeComponentAccumulationSnapshot3P
unit_inventory_from_source(
    const PostSnesPtFlashSourceCellSnapshot3D&
        source) {
    return {
        0.5,
        source.component_ids,
        source.overall_composition,
        1.0};
}

[[nodiscard]] inline
Sw92AuthoritativeTargetMaterialization3D
materialize_authoritative_target(
    mpmc::mesh::GlobalEntityId cell_global,
    const PostSnesPtFlashSourceCellSnapshot3D&
        source,
    const mpmc::flash::
        Sw92ProfileCPtPhaseSetResult&
            authoritative,
    const mpmc::thermodynamics::Sw92Phase<double>&
        model,
    const mpmc::flow::
        PhaseSetTransitionCandidate& candidate) {
    const auto* accepted =
        authoritative.solution
            .accepted_phase_set();
    if (accepted == nullptr ||
        accepted->phases.size() !=
            candidate.target_phase_count ||
        authoritative.phase_metadata.size() !=
            candidate.target_phase_count) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: SW92 target materialization cardinality mismatch");
    }

    std::vector<
        Sw92AuthoritativeTargetPhaseMaterialization3D>
        phases;
    phases.reserve(
        candidate.target_phase_count);
    for (std::size_t phase = 0U;
         phase < candidate.target_phase_count;
         ++phase) {
        const auto& published =
            accepted->phases[phase];
        const auto& metadata =
            authoritative.phase_metadata[phase];

        mpmc::thermodynamics::
            Sw92SelectedPhase<double>
                selection{
                    authoritative
                        .nacl_molality_mol_per_kg_water,
                    metadata.thermodynamic_family,
                    published.activity.branch,
                    {}};

        mpmc::thermodynamics::
            Sw92PhaseWorkspace<double>
                workspace;
        const auto density =
            mpmc::thermodynamics::
                evaluate_selected_phase_molar_density(
                    model,
                    authoritative.solution.pressure_pa,
                    authoritative.solution.temperature_k,
                    std::span<const double>{
                        published.composition},
                    selection,
                    workspace);

        if (!std::isfinite(
                density
                    .molar_density_mol_per_m3) ||
            !(density
                  .molar_density_mol_per_m3 >
              0.0) ||
            !near_roundoff(
                density.molar_density_mol_per_m3,
                candidate.target_phases[phase]
                    .molar_density_mol_per_m3) ||
            candidate.target_phases[phase]
                    .provider_activity_branch !=
                std::optional<std::size_t>{
                    published.activity.branch} ||
            candidate.target_phases[phase]
                    .provider_activity_smooth !=
                published.activity.smooth) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: SW92 target candidate lost authoritative density/root provenance");
        }

        phases.push_back({
            published.mole_phase_fraction,
            published.composition,
            metadata,
            selection,
            density.molar_density_mol_per_m3});
    }

    auto projection =
        mpmc::flow::
            project_phase_set_transition_candidate(
                candidate,
                unit_inventory_from_source(
                    source));

    return {
        cell_global,
        authoritative
            .nacl_molality_mol_per_kg_water,
        authoritative.dataset_id,
        authoritative.revision,
        authoritative.component_ids,
        authoritative.publication_convention,
        candidate.evidence_profile,
        candidate.diagnostic,
        std::move(phases),
        std::move(projection)};
}

} // namespace post_snes_sw92_profile_c_scanner_detail

[[nodiscard]] inline PetscErrorCode
scan_post_snes_sw92_profile_c_source_cell_3d(
    const PostSnesPtFlashSourceCellSnapshot3D& source,
    const mpmc::thermodynamics::Sw92Phase<double>& model,
    mpmc::flash::Sw92ProfileCPtFlashBackendOptions
        options,
    std::optional<
        PostSnesSw92ProfileCTransitionScanCellResult3D>*
            transition,
    PostSnesPhaseTransitionScanStatus3D*
        scan_status) {
    using namespace
        post_snes_sw92_profile_c_scanner_detail;

    if (transition == nullptr ||
        scan_status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    transition->reset();
    *scan_status =
        PostSnesPhaseTransitionScanStatus3D::
            complete;

    try {
        validate_model_source_identity(
            source,
            model,
            options);

        const auto owned_source =
            mpmc::flash::
                solve_sw92_phase_assigned_pt_boundary_aware(
                    source.pressure_pa,
                    source.temperature_k,
                    source.overall_composition,
                    model,
                    options
                        .nacl_molality_mol_per_kg_water,
                    options.flash);
        const auto authoritative =
            mpmc::flash::
                project_sw92_profile_c_pt_phase_set(
                    owned_source);
        const auto generic =
            make_generic_result_from_owned_source(
                owned_source,
                model,
                options);

        if (authoritative.solution.status !=
                mpmc::flash::
                    PtPhaseSetStatus::accepted ||
            !authoritative
                 .accepted_phase_set_published() ||
            generic.accepted_phase_set() ==
                nullptr) {
            *scan_status =
                PostSnesPhaseTransitionScanStatus3D::
                    indeterminate;
            return PETSC_SUCCESS;
        }

        validate_authoritative_projection_pair(
            generic,
            authoritative,
            source,
            model,
            options);

        const auto* accepted =
            authoritative.solution
                .accepted_phase_set();
        if (accepted->phases.size() ==
            source.source_phase_count) {
            return PETSC_SUCCESS;
        }

        std::vector<double> densities;
        densities.reserve(
            accepted->phases.size());
        for (std::size_t phase = 0U;
             phase < accepted->phases.size();
             ++phase) {
            const auto& published =
                accepted->phases[phase];
            const auto& metadata =
                authoritative.phase_metadata[phase];
            mpmc::thermodynamics::
                Sw92SelectedPhase<double>
                    selection{
                        authoritative
                            .nacl_molality_mol_per_kg_water,
                        metadata.thermodynamic_family,
                        published.activity.branch,
                        {}};
            mpmc::thermodynamics::
                Sw92PhaseWorkspace<double>
                    workspace;
            const auto density =
                mpmc::thermodynamics::
                    evaluate_selected_phase_molar_density(
                        model,
                        authoritative
                            .solution.pressure_pa,
                        authoritative
                            .solution.temperature_k,
                        std::span<const double>{
                            published.composition},
                        selection,
                        workspace);
            if (!std::isfinite(
                    density
                        .molar_density_mol_per_m3) ||
                !(density
                      .molar_density_mol_per_m3 >
                  0.0)) {
                *scan_status =
                    PostSnesPhaseTransitionScanStatus3D::
                        indeterminate;
                return PETSC_SUCCESS;
            }
            densities.push_back(
                density.molar_density_mol_per_m3);
        }

        auto candidate =
            mpmc::flow::
                make_phase_set_transition_candidate_from_flash(
                    source.source_phase_count,
                    generic,
                    densities);
        if (!candidate.has_value()) {
            *scan_status =
                PostSnesPhaseTransitionScanStatus3D::
                    indeterminate;
            return PETSC_SUCCESS;
        }

        auto target =
            materialize_authoritative_target(
                source.cell_global,
                source,
                authoritative,
                model,
                *candidate);
        PostSnesPhaseTransitionProposal3D
            proposal{
                source.cell_global,
                std::move(*candidate)};
        transition->emplace(
            PostSnesSw92ProfileCTransitionScanCellResult3D{
                std::move(proposal),
                std::move(target)});
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        return PETSC_ERR_ARG_INCOMP;
    } catch (const std::out_of_range&) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    } catch (const std::exception&) {
        *scan_status =
            PostSnesPhaseTransitionScanStatus3D::
                indeterminate;
        transition->reset();
        return PETSC_SUCCESS;
    }
}

[[nodiscard]] inline PetscErrorCode
scan_post_snes_sw92_profile_c_phase_transitions_3d(
    const PhaseTransitionRebuiltNaturalVariableSystem3D& system,
    Vec converged_state,
    const VariableCardinalityNaturalVariableSnesSolveReport3D&
        solve_report,
    void* raw_context,
    PostSnesPhaseTransitionScanStatus3D*
        scan_status,
    std::vector<
        PostSnesPhaseTransitionProposal3D>*
            local_owned_proposals) {
    using namespace
        post_snes_pt_flash_scanner_detail;

    if (raw_context == nullptr ||
        scan_status == nullptr ||
        local_owned_proposals == nullptr ||
        converged_state == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    local_owned_proposals->clear();
    *scan_status =
        PostSnesPhaseTransitionScanStatus3D::
            complete;

    if (static_cast<int>(
            solve_report.converged_reason) <=
        0) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }

    auto* context =
        static_cast<
            PostSnesSw92ProfileCPhaseTransitionScannerContext3D*>(
                raw_context);
    context->local_owned_targets.clear();
    if (context->model == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    const auto& numbering =
        system.numbering();
    for (const auto& record :
         numbering.cells()) {
        if (record.owner_rank !=
            numbering.local_rank()) {
            continue;
        }

        std::vector<PetscInt>
            indices(record.scalar_count);
        std::vector<PetscScalar>
            petsc_values(record.scalar_count);
        for (std::size_t slot = 0U;
             slot < record.scalar_count;
             ++slot) {
            indices[slot] =
                record.petsc_global_scalar_start +
                static_cast<PetscInt>(slot);
        }
        PetscErrorCode error =
            VecGetValues(
                converged_state,
                static_cast<PetscInt>(
                    indices.size()),
                indices.data(),
                petsc_values.data());
        if (error != PETSC_SUCCESS) {
            return error;
        }

        std::vector<double>
            natural_variables(
                record.scalar_count);
        for (std::size_t slot = 0U;
             slot < record.scalar_count;
             ++slot) {
            natural_variables[slot] =
                static_cast<double>(
                    PetscRealPart(
                        petsc_values[slot]));
        }

        std::optional<
            MixedCardinalityPhysicalCurrentCellLinearization3D>
            current;
        double porosity = 0.0;
        NaturalVariableSnesEvaluationStatus3D
            cell_status =
                NaturalVariableSnesEvaluationStatus3D::
                    success;
        error =
            system
                .evaluate_current_cell_for_phase_transition(
                    record.cell,
                    natural_variables,
                    &current,
                    &porosity,
                    &cell_status);
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (cell_status !=
                NaturalVariableSnesEvaluationStatus3D::
                    success ||
            !current.has_value() ||
            !has_common_pt_pressure(*current)) {
            *scan_status =
                PostSnesPhaseTransitionScanStatus3D::
                    indeterminate;
            local_owned_proposals->clear();
            context->local_owned_targets.clear();
            return PETSC_SUCCESS;
        }

        PostSnesPtFlashSourceCellSnapshot3D
            source;
        try {
            source =
                make_source_snapshot(
                    record.cell_global,
                    *current,
                    porosity);
        } catch (const std::exception&) {
            *scan_status =
                PostSnesPhaseTransitionScanStatus3D::
                    indeterminate;
            local_owned_proposals->clear();
            context->local_owned_targets.clear();
            return PETSC_SUCCESS;
        }

        std::optional<
            PostSnesSw92ProfileCTransitionScanCellResult3D>
            transition;
        PostSnesPhaseTransitionScanStatus3D
            cell_status_scan =
                PostSnesPhaseTransitionScanStatus3D::
                    complete;
        error =
            scan_post_snes_sw92_profile_c_source_cell_3d(
                source,
                *context->model,
                context->options,
                &transition,
                &cell_status_scan);
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (cell_status_scan ==
            PostSnesPhaseTransitionScanStatus3D::
                indeterminate) {
            *scan_status =
                PostSnesPhaseTransitionScanStatus3D::
                    indeterminate;
            local_owned_proposals->clear();
            context->local_owned_targets.clear();
            return PETSC_SUCCESS;
        }
        if (transition.has_value()) {
            if (context->find_target(
                    record.cell_global) !=
                nullptr) {
                return PETSC_ERR_PLIB;
            }
            local_owned_proposals->push_back(
                std::move(
                    transition->proposal));
            context->local_owned_targets.push_back(
                std::move(
                    transition->target));
        }
    }

    return PETSC_SUCCESS;
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_POST_SNES_SW92_PROFILE_C_PHASE_TRANSITION_SCANNER_HPP
