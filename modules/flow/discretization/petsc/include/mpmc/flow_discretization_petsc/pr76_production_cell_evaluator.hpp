#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_PR76_PRODUCTION_CELL_EVALUATOR_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_PR76_PRODUCTION_CELL_EVALUATOR_HPP

#include <mpmc/flow/pr76_selected_phase_property_closure.hpp>
#include <mpmc/flow_discretization_petsc/selected_phase_production_cell_evaluator.hpp>

#include <petscsys.h>

#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    pr76_production_cell_evaluator_convention =
        "flow_discretization_petsc/pr76-production-cell-evaluator/v1";

using Pr76RockThermalStorageLinearization3D =
    SelectedPhaseRockThermalStorageLinearization3D;
using Pr76RockThermalStorageEvaluator3D =
    SelectedPhaseRockThermalStorageEvaluator3D;
using Pr76RockThermalStorageEvaluatorBinding3D =
    SelectedPhaseRockThermalStorageEvaluatorBinding3D;
using Pr76TwoPhaseRelativePermeabilityLinearization3D =
    SelectedPhaseTwoPhaseRelativePermeabilityLinearization3D;
using Pr76TwoPhaseRelativePermeabilityEvaluator3D =
    SelectedPhaseTwoPhaseRelativePermeabilityEvaluator3D;
using Pr76TwoPhaseRelativePermeabilityEvaluatorBinding3D =
    SelectedPhaseTwoPhaseRelativePermeabilityEvaluatorBinding3D;
using Pr76ThreePhaseSaturationConstitutiveEvaluator3D =
    SelectedPhaseThreePhaseSaturationConstitutiveEvaluator3D;
using Pr76ThreePhaseSaturationConstitutiveEvaluatorBinding3D =
    SelectedPhaseThreePhaseSaturationConstitutiveEvaluatorBinding3D;

struct Pr76ProductionCellEvaluatorTraits {
    template <typename Function>
    [[nodiscard]] static PetscErrorCode
    translate_exceptions(
        Function&& function,
        NaturalVariableSnesEvaluationStatus3D* status) {
        try {
            return std::forward<Function>(function)();
        } catch (const mpmc::flow::
                     Pr76SelectedPhasePcNoneCapabilityError&) {
            return PETSC_ERR_SUP;
        } catch (const std::domain_error&) {
            if (status != nullptr) {
                *status =
                    NaturalVariableSnesEvaluationStatus3D::
                        domain_error;
            }
            return PETSC_SUCCESS;
        } catch (const std::invalid_argument&) {
            if (status != nullptr) {
                *status =
                    NaturalVariableSnesEvaluationStatus3D::
                        domain_error;
            }
            return PETSC_SUCCESS;
        } catch (const std::range_error&) {
            if (status != nullptr) {
                *status =
                    NaturalVariableSnesEvaluationStatus3D::
                        domain_error;
            }
            return PETSC_SUCCESS;
        } catch (const mpmc::thermodynamics::
                     Pr76PhaseError&) {
            if (status != nullptr) {
                *status =
                    NaturalVariableSnesEvaluationStatus3D::
                        domain_error;
            }
            return PETSC_SUCCESS;
        } catch (const std::logic_error&) {
            return PETSC_ERR_PLIB;
        } catch (const std::exception&) {
            return PETSC_ERR_LIB;
        }
    }

    template <typename PropertyClosure>
    [[nodiscard]] static auto
    evaluate_property_chart(
        const mpmc::flow::NaturalVariableLayoutDescriptor& layout,
        std::span<const std::string> component_ids,
        std::span<const double> natural_variables,
        const PropertyClosure& closure,
        mpmc::ad::RuntimeJacobianWorkspace<double, 4U>& workspace) {
        return mpmc::flow::
            evaluate_pr76_selected_phase_property_chart(
                layout,
                component_ids,
                natural_variables,
                closure,
                workspace);
    }

    [[nodiscard]] static auto
    make_flow_1p(
        const mpmc::flow::
            Pr76SelectedPhasePropertyChartLinearization& properties,
        std::span<const double> natural_variables,
        double relative_permeability,
        std::vector<double> relative_permeability_gradient) {
        return mpmc::flow::
            make_pr76_selected_phase_flow_linearization_1p(
                properties,
                natural_variables,
                relative_permeability,
                std::move(relative_permeability_gradient));
    }

    [[nodiscard]] static auto
    make_flow_2p(
        const mpmc::flow::
            Pr76SelectedPhasePropertyChartLinearization& properties,
        std::span<const double> natural_variables,
        std::array<double, 2> relative_permeability,
        std::array<std::vector<double>, 2>
            relative_permeability_gradient) {
        return mpmc::flow::
            make_pr76_selected_phase_flow_linearization_2p(
                properties,
                natural_variables,
                relative_permeability,
                std::move(relative_permeability_gradient));
    }

    [[nodiscard]] static auto
    make_flow_3p(
        const mpmc::flow::
            Pr76SelectedPhasePropertyChartLinearization& properties,
        std::span<const double> natural_variables) {
        return mpmc::flow::
            make_pr76_selected_phase_flow_linearization_3p(
                properties,
                natural_variables);
    }

    static void require_pc_none(
        const mpmc::flow::NaturalVariableCellState3P& state,
        const mpmc::flow::
            ThreePhaseSaturationConstitutiveNaturalVariableLinearization3P&
                constitutive) {
        mpmc::flow::
            require_pr76_selected_phase_pc_none_capability(
                state,
                constitutive);
    }
};

template <typename PropertyClosure>
using Pr76SinglePhaseProductionCellEvaluatorContext3D =
    SelectedPhaseSinglePhaseProductionCellEvaluatorContext3D<
        PropertyClosure,
        Pr76ProductionCellEvaluatorTraits>;
template <typename PropertyClosure>
using Pr76TwoPhaseProductionCellEvaluatorContext3D =
    SelectedPhaseTwoPhaseProductionCellEvaluatorContext3D<
        PropertyClosure,
        Pr76ProductionCellEvaluatorTraits>;
template <typename PropertyClosure>
using Pr76ThreePhaseProductionCellEvaluatorContext3D =
    SelectedPhaseThreePhaseProductionCellEvaluatorContext3D<
        PropertyClosure,
        Pr76ProductionCellEvaluatorTraits>;

template <typename PropertyClosure>
[[nodiscard]] inline PetscErrorCode
evaluate_pr76_single_phase_production_cell_3d(
    mpmc::mesh::LocalIndex cell,
    mpmc::mesh::GlobalEntityId cell_global,
    std::span<const double> natural_variables,
    const mpmc::flow::NaturalVariableLayout1P& frozen_layout,
    std::span<const std::string> component_ids,
    void* raw_context,
    std::optional<SinglePhaseCurrentCellLinearization3D>* output,
    NaturalVariableSnesEvaluationStatus3D* status) {
    return evaluate_selected_phase_single_phase_production_cell_3d<
        PropertyClosure,
        Pr76ProductionCellEvaluatorTraits>(
            cell,
            cell_global,
            natural_variables,
            frozen_layout,
            component_ids,
            raw_context,
            output,
            status);
}

template <typename PropertyClosure>
[[nodiscard]] inline PetscErrorCode
evaluate_pr76_two_phase_production_cell_3d(
    mpmc::mesh::LocalIndex cell,
    mpmc::mesh::GlobalEntityId cell_global,
    std::span<const double> natural_variables,
    const mpmc::flow::NaturalVariableLayout2P& frozen_layout,
    std::span<const std::string> component_ids,
    void* raw_context,
    std::optional<TwoPhaseCurrentCellLinearization3D>* output,
    NaturalVariableSnesEvaluationStatus3D* status) {
    return evaluate_selected_phase_two_phase_production_cell_3d<
        PropertyClosure,
        Pr76ProductionCellEvaluatorTraits>(
            cell,
            cell_global,
            natural_variables,
            frozen_layout,
            component_ids,
            raw_context,
            output,
            status);
}

template <typename PropertyClosure>
[[nodiscard]] inline PetscErrorCode
evaluate_pr76_three_phase_production_cell_3d(
    mpmc::mesh::LocalIndex cell,
    mpmc::mesh::GlobalEntityId cell_global,
    std::span<const double> natural_variables,
    const mpmc::flow::NaturalVariableLayout3P& frozen_layout,
    std::span<const std::string> component_ids,
    void* raw_context,
    std::optional<FixedThreePhaseCurrentCellLinearization3D>* output,
    NaturalVariableSnesEvaluationStatus3D* status) {
    return evaluate_selected_phase_three_phase_production_cell_3d<
        PropertyClosure,
        Pr76ProductionCellEvaluatorTraits>(
            cell,
            cell_global,
            natural_variables,
            frozen_layout,
            component_ids,
            raw_context,
            output,
            status);
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_PR76_PRODUCTION_CELL_EVALUATOR_HPP
