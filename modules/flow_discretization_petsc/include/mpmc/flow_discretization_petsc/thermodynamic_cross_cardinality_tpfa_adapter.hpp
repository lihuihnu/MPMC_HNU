#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_THERMODYNAMIC_CROSS_CARDINALITY_TPFA_ADAPTER_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_THERMODYNAMIC_CROSS_CARDINALITY_TPFA_ADAPTER_HPP

#include <mpmc/flow/thermodynamics_absent_phase_extension.hpp>
#include <mpmc/flow_discretization_petsc/cross_cardinality_tpfa_bridge.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <variant>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    thermodynamic_cross_cardinality_tpfa_adapter_convention =
        "flow_discretization_petsc/thermodynamic-cross-cardinality-tpfa-adapter/v1";

using AbsentPhaseThermodynamicCoordinateResolver3D =
    PetscErrorCode (*)(
        const MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D&
            face_input,
        const mpmc::flow::
            CrossCardinalityFacePhaseIdentityBinding&
                phase_binding,
        CrossCardinalityAbsentPhaseSide3D absent_side,
        const MixedCardinalityPhysicalCurrentCellLinearization3D&
            owner,
        const MixedCardinalityPhysicalCurrentCellLinearization3D&
            neighbour,
        void* user_context,
        std::optional<
            mpmc::flow::
                AbsentPhaseThermodynamicCoordinateExtension>*
                    output,
        NaturalVariableSnesEvaluationStatus3D*
            status);

using AbsentPhaseThermodynamicProviderEvaluator3D =
    PetscErrorCode (*)(
        const mpmc::flow::
            AbsentPhaseThermodynamicCoordinateExtension&
                coordinates,
        void* user_context,
        std::optional<
            mpmc::flow::
                AbsentPhasePotentialExtensionLinearization>*
                    output,
        NaturalVariableSnesEvaluationStatus3D*
            status);

struct ThermodynamicAbsentPhaseExtensionAdapterBinding3D {
    AbsentPhaseThermodynamicCoordinateResolver3D
        coordinate_resolver{};
    void* coordinate_resolver_context{};
    AbsentPhaseThermodynamicProviderEvaluator3D
        provider_evaluator{};
    void* provider_context{};
};

namespace thermodynamic_cross_cardinality_detail {

[[nodiscard]] inline const
mpmc::flow::NaturalVariableStateIdentity3P&
state_identity(
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        cell) {
    return std::visit(
        [](const auto& typed)
            -> const mpmc::flow::
                NaturalVariableStateIdentity3P& {
            return typed.transport.state_identity;
        },
        cell);
}

[[nodiscard]] inline bool
near_roundoff(
    double first,
    double second) {
    if (!std::isfinite(first) ||
        !std::isfinite(second)) {
        return false;
    }
    const double scale =
        std::max(
            {1.0,
             std::abs(first),
             std::abs(second)});
    return std::abs(first - second) <=
        8192.0 *
            std::numeric_limits<double>::
                epsilon() *
            scale;
}

[[nodiscard]] inline bool
same_layout(
    const mpmc::flow::
        NaturalVariableLayoutDescriptor& first,
    const mpmc::flow::
        NaturalVariableLayoutDescriptor& second) {
    if (first.component_count() !=
            second.component_count() ||
        first.phase_count() !=
            second.phase_count() ||
        first.unknown_count() !=
            second.unknown_count()) {
        return false;
    }
    for (std::size_t phase = 0U;
         phase < first.phase_count();
         ++phase) {
        const auto slot =
            static_cast<
                mpmc::flow::PhaseSlot3>(
                    phase);
        if (first
                .dependent_composition_component(
                    slot) !=
            second
                .dependent_composition_component(
                    slot)) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] inline bool
same_vector(
    std::span<const double> first,
    std::span<const double> second) {
    if (first.size() != second.size()) {
        return false;
    }
    for (std::size_t index = 0U;
         index < first.size();
         ++index) {
        if (!near_roundoff(
                first[index],
                second[index])) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] inline bool
same_state_identity(
    const mpmc::flow::
        NaturalVariableStateIdentity3P& first,
    const mpmc::flow::
        NaturalVariableStateIdentity3P& second) {
    if (!same_layout(
            first.layout,
            second.layout) ||
        first.component_ids !=
            second.component_ids ||
        !near_roundoff(
            first.reference_pressure_pa,
            second.reference_pressure_pa) ||
        !near_roundoff(
            first.temperature_k,
            second.temperature_k)) {
        return false;
    }

    for (std::size_t phase = 0U;
         phase <
             mpmc::flow::
                 fixed_three_phase_count;
         ++phase) {
        if (!near_roundoff(
                first.saturation[phase],
                second.saturation[phase]) ||
            !same_vector(
                first.phase_composition[phase],
                second.phase_composition[phase])) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] inline bool
requested_side_is_absent(
    const mpmc::flow::
        CrossCardinalityFacePhaseIdentityBinding&
            phase_binding,
    CrossCardinalityAbsentPhaseSide3D side) {
    switch (side) {
    case CrossCardinalityAbsentPhaseSide3D::
        owner:
        return !phase_binding
                    .owner_active_phase_index
                    .has_value() &&
            phase_binding
                .neighbour_active_phase_index
                .has_value();
    case CrossCardinalityAbsentPhaseSide3D::
        neighbour:
        return phase_binding
                   .owner_active_phase_index
                   .has_value() &&
            !phase_binding
                 .neighbour_active_phase_index
                 .has_value();
    }
    return false;
}

[[nodiscard]] inline const
MixedCardinalityPhysicalCurrentCellLinearization3D&
absent_cell(
    CrossCardinalityAbsentPhaseSide3D side,
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        owner,
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        neighbour) {
    return side ==
            CrossCardinalityAbsentPhaseSide3D::
                owner
        ? owner
        : neighbour;
}

} // namespace thermodynamic_cross_cardinality_detail

template <class Provider>
[[nodiscard]] inline PetscErrorCode
evaluate_absent_phase_thermodynamic_provider_3d(
    const mpmc::flow::
        AbsentPhaseThermodynamicCoordinateExtension&
            coordinates,
    void* raw_context,
    std::optional<
        mpmc::flow::
            AbsentPhasePotentialExtensionLinearization>*
                output,
    NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    *status =
        NaturalVariableSnesEvaluationStatus3D::
            success;

    try {
        const auto* provider =
            static_cast<const Provider*>(
                raw_context);
        output->emplace(
            provider->evaluate(
                coordinates));
        return PETSC_SUCCESS;
    } catch (const std::domain_error&) {
        *status =
            NaturalVariableSnesEvaluationStatus3D::
                domain_error;
        return PETSC_SUCCESS;
    } catch (const std::range_error&) {
        *status =
            NaturalVariableSnesEvaluationStatus3D::
                domain_error;
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        return PETSC_ERR_ARG_INCOMP;
    } catch (const std::out_of_range&) {
        return PETSC_ERR_ARG_INCOMP;
    } catch (...) {
        return PETSC_ERR_LIB;
    }
}

[[nodiscard]] inline PetscErrorCode
evaluate_thermodynamic_absent_phase_extension_3d(
    const MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D&
        face_input,
    const mpmc::flow::
        CrossCardinalityFacePhaseIdentityBinding&
            phase_binding,
    CrossCardinalityAbsentPhaseSide3D absent_side,
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        owner,
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        neighbour,
    void* raw_context,
    std::optional<
        mpmc::flow::
            AbsentPhasePotentialExtensionLinearization>*
                output,
    NaturalVariableSnesEvaluationStatus3D*
        status) {
    using namespace
        thermodynamic_cross_cardinality_detail;

    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    *status =
        NaturalVariableSnesEvaluationStatus3D::
            success;

    auto* binding =
        static_cast<
            ThermodynamicAbsentPhaseExtensionAdapterBinding3D*>(
                raw_context);
    if (binding->coordinate_resolver ==
            nullptr ||
        binding->provider_evaluator ==
            nullptr ||
        !requested_side_is_absent(
            phase_binding,
            absent_side)) {
        return PETSC_ERR_ARG_INCOMP;
    }

    std::optional<
        mpmc::flow::
            AbsentPhaseThermodynamicCoordinateExtension>
        coordinates;
    NaturalVariableSnesEvaluationStatus3D
        coordinate_status =
            NaturalVariableSnesEvaluationStatus3D::
                success;
    PetscErrorCode error =
        binding->coordinate_resolver(
            face_input,
            phase_binding,
            absent_side,
            owner,
            neighbour,
            binding
                ->coordinate_resolver_context,
            &coordinates,
            &coordinate_status);
    if (error != PETSC_SUCCESS) {
        return error;
    }
    if (coordinate_status ==
        NaturalVariableSnesEvaluationStatus3D::
            domain_error) {
        *status = coordinate_status;
        return PETSC_SUCCESS;
    }
    if (coordinate_status !=
            NaturalVariableSnesEvaluationStatus3D::
                success ||
        !coordinates.has_value()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    try {
        coordinates->validate();
    } catch (const std::invalid_argument&) {
        return PETSC_ERR_ARG_INCOMP;
    }

    const auto& current_identity =
        state_identity(
            absent_cell(
                absent_side,
                owner,
                neighbour));
    if (coordinates->identity !=
            phase_binding.identity ||
        !same_state_identity(
            coordinates
                ->host_state_identity,
            current_identity)) {
        return PETSC_ERR_ARG_INCOMP;
    }

    std::optional<
        mpmc::flow::
            AbsentPhasePotentialExtensionLinearization>
        extension;
    NaturalVariableSnesEvaluationStatus3D
        provider_status =
            NaturalVariableSnesEvaluationStatus3D::
                success;
    error =
        binding->provider_evaluator(
            *coordinates,
            binding->provider_context,
            &extension,
            &provider_status);
    if (error != PETSC_SUCCESS) {
        return error;
    }
    if (provider_status ==
        NaturalVariableSnesEvaluationStatus3D::
            domain_error) {
        *status = provider_status;
        return PETSC_SUCCESS;
    }
    if (provider_status !=
            NaturalVariableSnesEvaluationStatus3D::
                success ||
        !extension.has_value() ||
        extension->identity() !=
            phase_binding.identity ||
        extension->input_count() !=
            current_identity.layout
                .unknown_count() ||
        !near_roundoff(
            extension->phase_pressure_pa(),
            coordinates
                ->phase_pressure_pa) ||
        !same_vector(
            extension
                ->phase_pressure_gradient(),
            coordinates
                ->phase_pressure_gradient)) {
        return PETSC_ERR_ARG_INCOMP;
    }

    output->emplace(
        std::move(*extension));
    return PETSC_SUCCESS;
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_THERMODYNAMIC_CROSS_CARDINALITY_TPFA_ADAPTER_HPP
