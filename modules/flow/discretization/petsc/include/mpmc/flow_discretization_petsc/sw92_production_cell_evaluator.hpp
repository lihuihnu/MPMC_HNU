#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_SW92_PRODUCTION_CELL_EVALUATOR_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_SW92_PRODUCTION_CELL_EVALUATOR_HPP

#include <mpmc/flash/sw92_profile_c_phase_set.hpp>
#include <mpmc/flow/sw92_co2_water_properties.hpp>
#include <mpmc/flow_discretization_petsc/selected_phase_production_cell_evaluator.hpp>

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    sw92_production_cell_evaluator_convention =
        "flow_discretization_petsc/sw92-production-cell-evaluator/zero-salinity-co2-water/v1";

using Sw92RockThermalStorageLinearization3D =
    SelectedPhaseRockThermalStorageLinearization3D;
using Sw92RockThermalStorageEvaluator3D =
    SelectedPhaseRockThermalStorageEvaluator3D;
using Sw92RockThermalStorageEvaluatorBinding3D =
    SelectedPhaseRockThermalStorageEvaluatorBinding3D;
using Sw92TwoPhaseRelativePermeabilityLinearization3D =
    SelectedPhaseTwoPhaseRelativePermeabilityLinearization3D;
using Sw92TwoPhaseRelativePermeabilityEvaluator3D =
    SelectedPhaseTwoPhaseRelativePermeabilityEvaluator3D;
using Sw92TwoPhaseRelativePermeabilityEvaluatorBinding3D =
    SelectedPhaseTwoPhaseRelativePermeabilityEvaluatorBinding3D;
using Sw92ThreePhaseSaturationConstitutiveEvaluator3D =
    SelectedPhaseThreePhaseSaturationConstitutiveEvaluator3D;
using Sw92ThreePhaseSaturationConstitutiveEvaluatorBinding3D =
    SelectedPhaseThreePhaseSaturationConstitutiveEvaluatorBinding3D;

struct Sw92ProductionCellEvaluatorTraits {
    template <typename Function>
    [[nodiscard]] static PetscErrorCode
    translate_exceptions(
        Function&& function,
        NaturalVariableSnesEvaluationStatus3D* status) {
        try {
            return std::forward<Function>(function)();
        } catch (const mpmc::flow::
                     Sw92SelectedPhasePcNoneCapabilityError&) {
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
                     Sw92PhaseError&) {
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
            evaluate_sw92_selected_phase_property_chart(
                layout,
                component_ids,
                natural_variables,
                closure,
                workspace);
    }

    [[nodiscard]] static auto
    make_flow_1p(
        const mpmc::flow::
            Sw92SelectedPhasePropertyChartLinearization& properties,
        std::span<const double> natural_variables,
        double relative_permeability,
        std::vector<double> relative_permeability_gradient) {
        return mpmc::flow::
            make_sw92_selected_phase_flow_linearization_1p(
                properties,
                natural_variables,
                relative_permeability,
                std::move(relative_permeability_gradient));
    }

    [[nodiscard]] static auto
    make_flow_2p(
        const mpmc::flow::
            Sw92SelectedPhasePropertyChartLinearization& properties,
        std::span<const double> natural_variables,
        std::array<double, 2> relative_permeability,
        std::array<std::vector<double>, 2>
            relative_permeability_gradient) {
        return mpmc::flow::
            make_sw92_selected_phase_flow_linearization_2p(
                properties,
                natural_variables,
                relative_permeability,
                std::move(relative_permeability_gradient));
    }

    [[nodiscard]] static auto
    make_flow_3p(
        const mpmc::flow::
            Sw92SelectedPhasePropertyChartLinearization& properties,
        std::span<const double> natural_variables) {
        return mpmc::flow::
            make_sw92_selected_phase_flow_linearization_3p(
                properties,
                natural_variables);
    }

    static void require_pc_none(
        const mpmc::flow::NaturalVariableCellState3P& state,
        const mpmc::flow::
            ThreePhaseSaturationConstitutiveNaturalVariableLinearization3P&
                constitutive) {
        mpmc::flow::
            require_sw92_selected_phase_pc_none_capability(
                state,
                constitutive);
    }
};

template <typename PropertyClosure>
using Sw92SinglePhaseProductionCellEvaluatorContext3D =
    SelectedPhaseSinglePhaseProductionCellEvaluatorContext3D<
        PropertyClosure,
        Sw92ProductionCellEvaluatorTraits>;
template <typename PropertyClosure>
using Sw92TwoPhaseProductionCellEvaluatorContext3D =
    SelectedPhaseTwoPhaseProductionCellEvaluatorContext3D<
        PropertyClosure,
        Sw92ProductionCellEvaluatorTraits>;
template <typename PropertyClosure>
using Sw92ThreePhaseProductionCellEvaluatorContext3D =
    SelectedPhaseThreePhaseProductionCellEvaluatorContext3D<
        PropertyClosure,
        Sw92ProductionCellEvaluatorTraits>;

template <typename PropertyClosure>
[[nodiscard]] inline PetscErrorCode
evaluate_sw92_single_phase_production_cell_3d(
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
        Sw92ProductionCellEvaluatorTraits>(
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
evaluate_sw92_two_phase_production_cell_3d(
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
        Sw92ProductionCellEvaluatorTraits>(
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
evaluate_sw92_three_phase_production_cell_3d(
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
        Sw92ProductionCellEvaluatorTraits>(
            cell,
            cell_global,
            natural_variables,
            frozen_layout,
            component_ids,
            raw_context,
            output,
            status);
}

using Sw92FrozenNaturalVariableLayout3D =
    std::variant<
        mpmc::flow::NaturalVariableLayout1P,
        mpmc::flow::NaturalVariableLayout2P,
        mpmc::flow::NaturalVariableLayout3P>;

using Sw92Co2WaterSelectedPhasePropertyClosure3D =
    mpmc::flow::Sw92SelectedPhasePropertyClosure<
        double,
        mpmc::flow::Sw92Co2WaterPropertyProvider>;

struct Sw92Co2WaterFrozenCellMaterialization3D {
    Sw92FrozenNaturalVariableLayout3D frozen_layout;
    std::vector<double> natural_variables;
    std::vector<std::string> component_ids;
    Sw92Co2WaterSelectedPhasePropertyClosure3D
        property_closure;
};

namespace sw92_production_cell_evaluator_detail {

[[nodiscard]] inline std::size_t
dependent_component(
    std::span<const double> composition) {
    if (composition.size() < 2U) {
        throw std::invalid_argument(
            "SW92 production materializer requires at least two components");
    }
    const auto found =
        std::max_element(
            composition.begin(),
            composition.end());
    if (found == composition.end() ||
        !std::isfinite(*found) ||
        !(*found > 0.0)) {
        throw std::invalid_argument(
            "SW92 production materializer requires positive finite phase composition");
    }
    for (double value : composition) {
        if (!std::isfinite(value) ||
            !(value > 0.0)) {
            throw std::invalid_argument(
                "SW92 production materializer requires strict-positive phase composition");
        }
    }
    return static_cast<std::size_t>(
        std::distance(
            composition.begin(),
            found));
}

inline void validate_source_snapshot(
    const mpmc::flash::Sw92ProfileCPtPhaseSetResult& source,
    const mpmc::thermodynamics::Sw92Phase<double>& model) {
    if (source.solution.status !=
            mpmc::flash::PtPhaseSetStatus::accepted ||
        !source.accepted_phase_set_published() ||
        source.nacl_molality_mol_per_kg_water != 0.0 ||
        source.dataset_id !=
            model.parameters().dataset_id() ||
        source.revision !=
            model.parameters().revision() ||
        source.component_ids.size() !=
            model.size() ||
        source.solution.feed.size() !=
            model.size() ||
        source.solution.temperature_k <
            mpmc::flow::sw92_co2_water_property_detail::
                minimum_temperature_k ||
        source.solution.temperature_k >
            mpmc::flow::sw92_co2_water_property_detail::
                maximum_temperature_k) {
        throw std::invalid_argument(
            "SW92 production materializer requires an accepted zero-salinity sourced CO2/H2O Profile-C snapshot in the provider temperature interval");
    }
    for (std::size_t component = 0U;
         component < model.size();
         ++component) {
        if (source.component_ids[component] !=
            model.parameters()
                .components()
                .at(component)
                .id) {
            throw std::invalid_argument(
                "SW92 production materializer component identity/order mismatch");
        }
    }
}

} // namespace sw92_production_cell_evaluator_detail

[[nodiscard]] inline
Sw92Co2WaterFrozenCellMaterialization3D
materialize_sw92_co2_water_profile_c_frozen_cell_3d(
    const mpmc::flash::Sw92ProfileCPtPhaseSetResult& source,
    const mpmc::thermodynamics::Sw92Phase<double>& model) {
    using namespace sw92_production_cell_evaluator_detail;
    validate_source_snapshot(source, model);

    const auto* accepted =
        source.solution.accepted_phase_set();
    if (accepted == nullptr ||
        accepted->phases.empty() ||
        accepted->phases.size() > 3U ||
        accepted->phases.size() !=
            source.phase_metadata.size()) {
        throw std::invalid_argument(
            "SW92 production materializer accepted phase-set shape mismatch");
    }

    const std::size_t n = model.size();
    const std::size_t p =
        accepted->phases.size();
    std::vector<std::size_t> dependent;
    std::vector<
        mpmc::thermodynamics::Sw92SelectedPhase<double>>
        selections;
    std::vector<double> volume_weights;
    dependent.reserve(p);
    selections.reserve(p);
    volume_weights.reserve(p);

    for (std::size_t phase = 0U;
         phase < p;
         ++phase) {
        const auto& published =
            accepted->phases[phase];
        const auto& metadata =
            source.phase_metadata[phase];
        if (!std::isfinite(
                published.mole_phase_fraction) ||
            !(published.mole_phase_fraction > 0.0) ||
            published.composition.size() != n ||
            !published.activity.smooth ||
            published.activity.ln_phi.size() != n) {
            throw std::invalid_argument(
                "SW92 production materializer received malformed accepted phase data");
        }
        dependent.push_back(
            dependent_component(
                published.composition));
        selections.push_back({
            source.nacl_molality_mol_per_kg_water,
            metadata.thermodynamic_family,
            published.activity.branch,
            {}});

        mpmc::thermodynamics::
            Sw92PhaseWorkspace<double>
                workspace;
        const auto density =
            mpmc::thermodynamics::
                evaluate_selected_phase_molar_density(
                    model,
                    source.solution.pressure_pa,
                    source.solution.temperature_k,
                    std::span<const double>{
                        published.composition},
                    selections.back(),
                    workspace);
        const double weight =
            published.mole_phase_fraction /
            density.molar_density_mol_per_m3;
        if (!std::isfinite(weight) ||
            !(weight > 0.0)) {
            throw std::range_error(
                "SW92 production materializer cannot construct positive phase volume");
        }
        volume_weights.push_back(weight);
    }

    const double volume_sum =
        std::accumulate(
            volume_weights.begin(),
            volume_weights.end(),
            0.0);
    if (!std::isfinite(volume_sum) ||
        !(volume_sum > 0.0)) {
        throw std::range_error(
            "SW92 production materializer phase-volume sum is invalid");
    }

    mpmc::flow::NaturalVariableLayoutDescriptor descriptor{
        n,
        p,
        dependent};
    std::vector<double> q(
        descriptor.unknown_count(),
        0.0);
    q[descriptor.pressure_unknown_index()] =
        source.solution.pressure_pa;
    q[descriptor.temperature_unknown_index()] =
        source.solution.temperature_k;

    for (std::size_t phase = 0U;
         phase + 1U < p;
         ++phase) {
        const auto column =
            descriptor.independent_saturation_unknown_index(
                static_cast<mpmc::flow::PhaseSlot3>(
                    phase));
        if (!column) {
            throw std::logic_error(
                "SW92 production materializer saturation column missing");
        }
        q[*column] =
            volume_weights[phase] /
            volume_sum;
    }

    for (std::size_t phase = 0U;
         phase < p;
         ++phase) {
        const auto slot =
            static_cast<mpmc::flow::PhaseSlot3>(
                phase);
        const auto& composition =
            accepted->phases[phase].composition;
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const auto column =
                descriptor.independent_composition_unknown_index(
                    slot,
                    component);
            if (column) {
                q[*column] =
                    composition[component];
            }
        }
    }

    Sw92FrozenNaturalVariableLayout3D layout =
        p == 1U
            ? Sw92FrozenNaturalVariableLayout3D{
                  mpmc::flow::NaturalVariableLayout1P{
                      mpmc::flow::
                          NaturalVariableCompositionPivot1P::
                              from_dependent_component(
                                  n,
                                  dependent[0])}}
            : (p == 2U
                   ? Sw92FrozenNaturalVariableLayout3D{
                         mpmc::flow::NaturalVariableLayout2P{
                             mpmc::flow::
                                 NaturalVariableCompositionPivot2P::
                                     from_dependent_components(
                                         n,
                                         std::array<std::size_t, 2>{
                                             dependent[0],
                                             dependent[1]})}}
                   : Sw92FrozenNaturalVariableLayout3D{
                         mpmc::flow::NaturalVariableLayout3P{
                             mpmc::flow::
                                 NaturalVariableCompositionPivot3P::
                                     from_dependent_components(
                                         n,
                                         std::array<std::size_t, 3>{
                                             dependent[0],
                                             dependent[1],
                                             dependent[2]})}});

    auto closure =
        mpmc::flow::
            make_sw92_co2_water_selected_phase_property_closure(
                model,
                std::move(selections));

    return {
        std::move(layout),
        std::move(q),
        source.component_ids,
        std::move(closure)};
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_SW92_PRODUCTION_CELL_EVALUATOR_HPP
