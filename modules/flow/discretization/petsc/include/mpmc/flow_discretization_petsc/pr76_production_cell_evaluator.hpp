#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_PR76_PRODUCTION_CELL_EVALUATOR_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_PR76_PRODUCTION_CELL_EVALUATOR_HPP

#include <mpmc/flow/pr76_selected_phase_property_closure.hpp>
#include <mpmc/flow_discretization_petsc/fixed_three_phase_snes_assembly.hpp>
#include <mpmc/flow_discretization_petsc/single_phase_snes_assembly.hpp>
#include <mpmc/flow_discretization_petsc/two_phase_snes_assembly.hpp>

#include <petscsys.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    pr76_production_cell_evaluator_convention =
        "flow_discretization_petsc/pr76-production-cell-evaluator/v1";

struct Pr76RockThermalStorageLinearization3D {
    double volumetric_internal_energy_j_per_rock_m3{};
    std::vector<double> gradient;
    mpmc::flow::TransportPropertyProvenance provenance;
};

using Pr76RockThermalStorageEvaluator3D =
    PetscErrorCode (*)(
        const mpmc::flow::NaturalVariableLayoutDescriptor& layout,
        std::span<const double> natural_variables,
        void* user_context,
        std::optional<
            Pr76RockThermalStorageLinearization3D>* output,
        NaturalVariableSnesEvaluationStatus3D* status);

struct Pr76RockThermalStorageEvaluatorBinding3D {
    Pr76RockThermalStorageEvaluator3D evaluator{};
    void* user_context{};
};

struct Pr76TwoPhaseRelativePermeabilityLinearization3D {
    std::array<double, 2> relative_permeability{};
    std::array<std::vector<double>, 2> gradient;
};

using Pr76TwoPhaseRelativePermeabilityEvaluator3D =
    PetscErrorCode (*)(
        const mpmc::flow::NaturalVariableCellState2P& state,
        void* user_context,
        std::optional<
            Pr76TwoPhaseRelativePermeabilityLinearization3D>* output,
        NaturalVariableSnesEvaluationStatus3D* status);

struct Pr76TwoPhaseRelativePermeabilityEvaluatorBinding3D {
    Pr76TwoPhaseRelativePermeabilityEvaluator3D evaluator{};
    void* user_context{};
};

using Pr76ThreePhaseSaturationConstitutiveEvaluator3D =
    PetscErrorCode (*)(
        const mpmc::flow::NaturalVariableCellState3P& state,
        void* user_context,
        std::optional<
            mpmc::flow::
                ThreePhaseSaturationConstitutiveNaturalVariableLinearization3P>*
                    output,
        NaturalVariableSnesEvaluationStatus3D* status);

struct Pr76ThreePhaseSaturationConstitutiveEvaluatorBinding3D {
    Pr76ThreePhaseSaturationConstitutiveEvaluator3D evaluator{};
    void* user_context{};
};

namespace pr76_production_cell_evaluator_detail {

inline void validate_rock(
    const Pr76RockThermalStorageLinearization3D& rock,
    std::size_t q) {
    if (!std::isfinite(
            rock.volumetric_internal_energy_j_per_rock_m3) ||
        rock.gradient.size() != q) {
        throw std::invalid_argument(
            "PR76 production rock thermal-storage shape/value mismatch");
    }
    for (const double value : rock.gradient) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "PR76 production rock thermal-storage gradient is non-finite");
        }
    }
    mpmc::flow::phase_transport_detail::
        validate_provenance(
            rock.provenance,
            "PR76 production rock thermal storage");
}

inline void require_binding(
    const Pr76RockThermalStorageEvaluatorBinding3D& binding) {
    if (binding.evaluator == nullptr) {
        throw std::invalid_argument(
            "PR76 production rock thermal-storage evaluator is required");
    }
}

inline void require_binding(
    const Pr76TwoPhaseRelativePermeabilityEvaluatorBinding3D&
        binding) {
    if (binding.evaluator == nullptr) {
        throw std::invalid_argument(
            "PR76 production two-phase relative-permeability evaluator is required");
    }
}

inline void require_binding(
    const Pr76ThreePhaseSaturationConstitutiveEvaluatorBinding3D&
        binding) {
    if (binding.evaluator == nullptr) {
        throw std::invalid_argument(
            "PR76 production three-phase saturation evaluator is required");
    }
}

template <typename Function>
[[nodiscard]] PetscErrorCode translate_model_exceptions(
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
    } catch (const mpmc::thermodynamics::Pr76PhaseError&) {
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

inline std::optional<
    Pr76RockThermalStorageLinearization3D>
evaluate_rock(
    const Pr76RockThermalStorageEvaluatorBinding3D& binding,
    const mpmc::flow::NaturalVariableLayoutDescriptor& layout,
    std::span<const double> q,
    NaturalVariableSnesEvaluationStatus3D* status,
    PetscErrorCode* error) {
    require_binding(binding);
    std::optional<
        Pr76RockThermalStorageLinearization3D>
        rock;
    *error =
        binding.evaluator(
            layout,
            q,
            binding.user_context,
            &rock,
            status);
    if (*error != PETSC_SUCCESS ||
        *status !=
            NaturalVariableSnesEvaluationStatus3D::
                success) {
        return std::nullopt;
    }
    if (!rock) {
        *error = PETSC_ERR_ARG_WRONGSTATE;
        return std::nullopt;
    }
    validate_rock(
        *rock,
        layout.unknown_count());
    return rock;
}

} // namespace pr76_production_cell_evaluator_detail

template <typename PropertyClosure>
struct Pr76SinglePhaseProductionCellEvaluatorContext3D {
    const PropertyClosure* property_closure{};
    double relative_permeability{1.0};
    Pr76RockThermalStorageEvaluatorBinding3D rock;
    mpmc::ad::RuntimeJacobianWorkspace<double, 4U>
        workspace;
};

template <typename PropertyClosure>
[[nodiscard]] PetscErrorCode
evaluate_pr76_single_phase_production_cell_3d(
    mpmc::mesh::LocalIndex,
    mpmc::mesh::GlobalEntityId,
    std::span<const double> natural_variables,
    const mpmc::flow::NaturalVariableLayout1P&
        frozen_layout,
    std::span<const std::string> component_ids,
    void* raw_context,
    std::optional<
        SinglePhaseCurrentCellLinearization3D>* output,
    NaturalVariableSnesEvaluationStatus3D* status) {
    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    *status =
        NaturalVariableSnesEvaluationStatus3D::
            success;

    auto* context =
        static_cast<
            Pr76SinglePhaseProductionCellEvaluatorContext3D<
                PropertyClosure>*>(raw_context);

    return pr76_production_cell_evaluator_detail::
        translate_model_exceptions(
            [&]() -> PetscErrorCode {
                if (context->property_closure == nullptr ||
                    context->property_closure->phase_count() != 1U ||
                    natural_variables.size() !=
                        frozen_layout.unknown_count() ||
                    component_ids.size() !=
                        frozen_layout.component_count() ||
                    !std::isfinite(
                        context->relative_permeability) ||
                    context->relative_permeability < 0.0) {
                    return PETSC_ERR_ARG_INCOMP;
                }

                const auto descriptor =
                    frozen_layout.descriptor();
                auto properties =
                    mpmc::flow::
                        evaluate_pr76_selected_phase_property_chart(
                            descriptor,
                            component_ids,
                            natural_variables,
                            *context->property_closure,
                            context->workspace);

                std::vector<double> dkr(
                    descriptor.unknown_count(),
                    0.0);
                auto flow =
                    mpmc::flow::
                        make_pr76_selected_phase_flow_linearization_1p(
                            properties,
                            natural_variables,
                            context->relative_permeability,
                            std::move(dkr));

                PetscErrorCode error = PETSC_SUCCESS;
                auto rock =
                    pr76_production_cell_evaluator_detail::
                        evaluate_rock(
                            context->rock,
                            descriptor,
                            natural_variables,
                            status,
                            &error);
                if (error != PETSC_SUCCESS ||
                    *status !=
                        NaturalVariableSnesEvaluationStatus3D::
                            success) {
                    return error;
                }

                auto rock_linearization =
                    mpmc::flow::
                        make_single_phase_rock_thermal_storage_linearization(
                            flow.state,
                            rock->volumetric_internal_energy_j_per_rock_m3,
                            std::move(rock->gradient),
                            std::move(rock->provenance));

                output->emplace(
                    SinglePhaseCurrentCellLinearization3D{
                        std::move(flow.state),
                        std::move(flow.molar_density),
                        std::move(flow.transport),
                        std::move(flow.caloric),
                        std::move(rock_linearization)});
                return PETSC_SUCCESS;
            },
            status);
}

template <typename PropertyClosure>
struct Pr76TwoPhaseProductionCellEvaluatorContext3D {
    const PropertyClosure* property_closure{};
    Pr76TwoPhaseRelativePermeabilityEvaluatorBinding3D
        relative_permeability;
    Pr76RockThermalStorageEvaluatorBinding3D rock;
    mpmc::ad::RuntimeJacobianWorkspace<double, 4U>
        workspace;
};

template <typename PropertyClosure>
[[nodiscard]] PetscErrorCode
evaluate_pr76_two_phase_production_cell_3d(
    mpmc::mesh::LocalIndex,
    mpmc::mesh::GlobalEntityId,
    std::span<const double> natural_variables,
    const mpmc::flow::NaturalVariableLayout2P&
        frozen_layout,
    std::span<const std::string> component_ids,
    void* raw_context,
    std::optional<
        TwoPhaseCurrentCellLinearization3D>* output,
    NaturalVariableSnesEvaluationStatus3D* status) {
    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    *status =
        NaturalVariableSnesEvaluationStatus3D::
            success;
    auto* context =
        static_cast<
            Pr76TwoPhaseProductionCellEvaluatorContext3D<
                PropertyClosure>*>(raw_context);

    return pr76_production_cell_evaluator_detail::
        translate_model_exceptions(
            [&]() -> PetscErrorCode {
                pr76_production_cell_evaluator_detail::
                    require_binding(
                        context->relative_permeability);
                if (context->property_closure == nullptr ||
                    context->property_closure->phase_count() != 2U ||
                    natural_variables.size() !=
                        frozen_layout.unknown_count() ||
                    component_ids.size() !=
                        frozen_layout.component_count()) {
                    return PETSC_ERR_ARG_INCOMP;
                }

                const auto descriptor =
                    frozen_layout.descriptor();
                auto properties =
                    mpmc::flow::
                        evaluate_pr76_selected_phase_property_chart(
                            descriptor,
                            component_ids,
                            natural_variables,
                            *context->property_closure,
                            context->workspace);

                // Build the exact state once without inventing kr so the
                // configured constitutive callback can evaluate on it.
                std::array<double, 2> placeholder_kr{1.0, 1.0};
                std::array<std::vector<double>, 2>
                    placeholder_gradient{
                        std::vector<double>(
                            descriptor.unknown_count(),
                            0.0),
                        std::vector<double>(
                            descriptor.unknown_count(),
                            0.0)};
                auto provisional =
                    mpmc::flow::
                        make_pr76_selected_phase_flow_linearization_2p(
                            properties,
                            natural_variables,
                            placeholder_kr,
                            std::move(
                                placeholder_gradient));

                std::optional<
                    Pr76TwoPhaseRelativePermeabilityLinearization3D>
                    relative;
                PetscErrorCode error =
                    context->relative_permeability.evaluator(
                        provisional.state,
                        context->relative_permeability.user_context,
                        &relative,
                        status);
                if (error != PETSC_SUCCESS ||
                    *status !=
                        NaturalVariableSnesEvaluationStatus3D::
                            success) {
                    return error;
                }
                if (!relative) {
                    return PETSC_ERR_ARG_WRONGSTATE;
                }

                auto flow =
                    mpmc::flow::
                        make_pr76_selected_phase_flow_linearization_2p(
                            properties,
                            natural_variables,
                            relative->relative_permeability,
                            std::move(relative->gradient));

                auto rock =
                    pr76_production_cell_evaluator_detail::
                        evaluate_rock(
                            context->rock,
                            descriptor,
                            natural_variables,
                            status,
                            &error);
                if (error != PETSC_SUCCESS ||
                    *status !=
                        NaturalVariableSnesEvaluationStatus3D::
                            success) {
                    return error;
                }

                auto rock_linearization =
                    mpmc::flow::
                        make_two_phase_rock_thermal_storage_linearization(
                            flow.state,
                            rock->volumetric_internal_energy_j_per_rock_m3,
                            std::move(rock->gradient),
                            std::move(rock->provenance));

                output->emplace(
                    TwoPhaseCurrentCellLinearization3D{
                        std::move(flow.state),
                        std::move(flow.molar_density),
                        std::move(flow.transport),
                        std::move(flow.caloric),
                        std::move(rock_linearization),
                        std::move(flow.fugacity)});
                return PETSC_SUCCESS;
            },
            status);
}

template <typename PropertyClosure>
struct Pr76ThreePhaseProductionCellEvaluatorContext3D {
    const PropertyClosure* property_closure{};
    Pr76ThreePhaseSaturationConstitutiveEvaluatorBinding3D
        saturation_constitutive;
    Pr76RockThermalStorageEvaluatorBinding3D rock;
    mpmc::ad::RuntimeJacobianWorkspace<double, 4U>
        workspace;
};

template <typename PropertyClosure>
[[nodiscard]] PetscErrorCode
evaluate_pr76_three_phase_production_cell_3d(
    mpmc::mesh::LocalIndex,
    mpmc::mesh::GlobalEntityId,
    std::span<const double> natural_variables,
    const mpmc::flow::NaturalVariableLayout3P&
        frozen_layout,
    std::span<const std::string> component_ids,
    void* raw_context,
    std::optional<
        FixedThreePhaseCurrentCellLinearization3D>* output,
    NaturalVariableSnesEvaluationStatus3D* status) {
    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    *status =
        NaturalVariableSnesEvaluationStatus3D::
            success;
    auto* context =
        static_cast<
            Pr76ThreePhaseProductionCellEvaluatorContext3D<
                PropertyClosure>*>(raw_context);

    return pr76_production_cell_evaluator_detail::
        translate_model_exceptions(
            [&]() -> PetscErrorCode {
                pr76_production_cell_evaluator_detail::
                    require_binding(
                        context->saturation_constitutive);
                if (context->property_closure == nullptr ||
                    context->property_closure->phase_count() != 3U ||
                    natural_variables.size() !=
                        frozen_layout.unknown_count() ||
                    component_ids.size() !=
                        frozen_layout.component_count()) {
                    return PETSC_ERR_ARG_INCOMP;
                }

                const mpmc::flow::NaturalVariableLayoutDescriptor
                    descriptor{frozen_layout};
                auto properties =
                    mpmc::flow::
                        evaluate_pr76_selected_phase_property_chart(
                            descriptor,
                            component_ids,
                            natural_variables,
                            *context->property_closure,
                            context->workspace);
                auto flow =
                    mpmc::flow::
                        make_pr76_selected_phase_flow_linearization_3p(
                            properties,
                            natural_variables);

                std::optional<
                    mpmc::flow::
                        ThreePhaseSaturationConstitutiveNaturalVariableLinearization3P>
                    saturation;
                PetscErrorCode error =
                    context->saturation_constitutive.evaluator(
                        flow.state,
                        context->saturation_constitutive.user_context,
                        &saturation,
                        status);
                if (error != PETSC_SUCCESS ||
                    *status !=
                        NaturalVariableSnesEvaluationStatus3D::
                            success) {
                    return error;
                }
                if (!saturation) {
                    return PETSC_ERR_ARG_WRONGSTATE;
                }

                mpmc::flow::
                    require_pr76_selected_phase_pc_none_capability(
                        flow.state,
                        *saturation);

                auto rock =
                    pr76_production_cell_evaluator_detail::
                        evaluate_rock(
                            context->rock,
                            descriptor,
                            natural_variables,
                            status,
                            &error);
                if (error != PETSC_SUCCESS ||
                    *status !=
                        NaturalVariableSnesEvaluationStatus3D::
                            success) {
                    return error;
                }

                auto rock_linearization =
                    mpmc::flow::
                        make_stationary_rock_thermal_storage_linearization(
                            flow.state,
                            rock->volumetric_internal_energy_j_per_rock_m3,
                            std::move(rock->gradient),
                            std::move(rock->provenance));

                output->emplace(
                    FixedThreePhaseCurrentCellLinearization3D{
                        std::move(flow.state),
                        std::move(flow.molar_density),
                        std::move(flow.transport),
                        std::move(flow.caloric),
                        std::move(rock_linearization),
                        std::move(*saturation),
                        std::move(flow.fugacity)});
                return PETSC_SUCCESS;
            },
            status);
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_PR76_PRODUCTION_CELL_EVALUATOR_HPP
