#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_CROSS_CARDINALITY_TPFA_BRIDGE_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_CROSS_CARDINALITY_TPFA_BRIDGE_HPP

#include <mpmc/flow_discretization_petsc/mixed_cardinality_physical_snes_assembly.hpp>

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
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    cross_cardinality_tpfa_bridge_convention =
        "flow_discretization_petsc/cross-cardinality-tpfa-bridge/v1";

enum class CrossCardinalityAbsentPhaseSide3D {
    owner,
    neighbour
};

using CrossCardinalityAbsentPhasePotentialExtensionEvaluator3D =
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
                AbsentPhasePotentialExtensionLinearization>*
                    output,
        NaturalVariableSnesEvaluationStatus3D*
            status);

struct CrossCardinalityTpfaBridgeBinding3D {
    CrossCardinalityAbsentPhasePotentialExtensionEvaluator3D
        absent_phase_extension_evaluator{};
    void* user_context{};
};

namespace cross_cardinality_tpfa_detail {

struct ActivePhaseAdvectivePayload3D {
    double phase_pressure_pa{};
    double mass_density_kg_per_m3{};
    double mobility_per_pa_s{};
    std::vector<double>
        phase_pressure_gradient;
    std::vector<double>
        mass_density_gradient;
    std::vector<double>
        mobility_gradient;

    std::vector<double>
        component_molar_content_mol_per_m3;
    std::vector<double>
        component_molar_content_jacobian;

    double advective_energy_density_j_per_m3{};
    std::vector<double>
        advective_energy_density_gradient;
};

struct PotentialSidePayload3D {
    bool active{};
    double phase_pressure_pa{};
    double mass_density_kg_per_m3{};
    double mobility_per_pa_s{};
    std::vector<double>
        phase_pressure_gradient;
    std::vector<double>
        mass_density_gradient;
    std::vector<double>
        mobility_gradient;
    std::optional<
        ActivePhaseAdvectivePayload3D>
        active_advective;
};

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

[[nodiscard]] inline std::size_t
phase_count(
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        cell) {
    return state_identity(cell)
        .layout.phase_count();
}

[[nodiscard]] inline std::size_t
component_count(
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        cell) {
    return state_identity(cell)
        .layout.component_count();
}

[[nodiscard]] inline bool
near_roundoff(
    double first,
    double second,
    double extra_scale = 0.0) {
    if (!std::isfinite(first) ||
        !std::isfinite(second) ||
        !std::isfinite(extra_scale) ||
        extra_scale < 0.0) {
        return false;
    }
    const double scale =
        std::max(
            {1.0,
             std::abs(first),
             std::abs(second),
             extra_scale});
    return std::abs(first - second) <=
        16384.0 *
            std::numeric_limits<double>::epsilon() *
            scale;
}

inline void require_gradient(
    std::span<const double> gradient,
    std::size_t q,
    const char* quantity) {
    if (gradient.size() != q) {
        throw std::invalid_argument(
            std::string{
                "mpmc::flow_discretization_petsc: cross-cardinality "} +
            quantity +
            " gradient shape mismatch");
    }
    for (const double value : gradient) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                std::string{
                    "mpmc::flow_discretization_petsc: cross-cardinality "} +
                quantity +
                " gradient contains non-finite derivative");
        }
    }
}

template <class CompositionDerivative>
[[nodiscard]] inline ActivePhaseAdvectivePayload3D
make_active_payload(
    double phase_pressure_pa,
    double mass_density_kg_per_m3,
    double mobility_per_pa_s,
    std::span<const double> phase_pressure_gradient,
    std::span<const double> mass_density_gradient,
    std::span<const double> mobility_gradient,
    double molar_density_mol_per_m3,
    std::span<const double> molar_density_gradient,
    std::span<const double> composition,
    double specific_enthalpy_j_per_kg,
    std::span<const double> specific_enthalpy_gradient,
    std::size_t q,
    CompositionDerivative&& d_composition) {
    if (!std::isfinite(phase_pressure_pa) ||
        !(phase_pressure_pa > 0.0) ||
        !std::isfinite(mass_density_kg_per_m3) ||
        !(mass_density_kg_per_m3 > 0.0) ||
        !std::isfinite(mobility_per_pa_s) ||
        mobility_per_pa_s < 0.0 ||
        !std::isfinite(molar_density_mol_per_m3) ||
        !(molar_density_mol_per_m3 > 0.0) ||
        !std::isfinite(specific_enthalpy_j_per_kg) ||
        composition.empty()) {
        throw std::invalid_argument(
            "mpmc::flow_discretization_petsc: invalid active cross-cardinality phase payload");
    }

    require_gradient(
        phase_pressure_gradient,
        q,
        "phase-pressure");
    require_gradient(
        mass_density_gradient,
        q,
        "mass-density");
    require_gradient(
        mobility_gradient,
        q,
        "mobility");
    require_gradient(
        molar_density_gradient,
        q,
        "molar-density");
    require_gradient(
        specific_enthalpy_gradient,
        q,
        "enthalpy");

    ActivePhaseAdvectivePayload3D
        result;
    result.phase_pressure_pa =
        phase_pressure_pa;
    result.mass_density_kg_per_m3 =
        mass_density_kg_per_m3;
    result.mobility_per_pa_s =
        mobility_per_pa_s;
    result.phase_pressure_gradient.assign(
        phase_pressure_gradient.begin(),
        phase_pressure_gradient.end());
    result.mass_density_gradient.assign(
        mass_density_gradient.begin(),
        mass_density_gradient.end());
    result.mobility_gradient.assign(
        mobility_gradient.begin(),
        mobility_gradient.end());

    const std::size_t n =
        composition.size();
    result.component_molar_content_mol_per_m3
        .resize(n);
    result.component_molar_content_jacobian
        .assign(
            n * q,
            0.0);

    for (std::size_t component = 0U;
         component < n;
         ++component) {
        const double x =
            composition[component];
        if (!std::isfinite(x) ||
            !(x > 0.0)) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: active phase composition must stay on positive support");
        }
        result.component_molar_content_mol_per_m3[
            component] =
            molar_density_mol_per_m3 *
            x;

        for (std::size_t column = 0U;
             column < q;
             ++column) {
            const double derivative =
                molar_density_gradient[column] *
                    x +
                molar_density_mol_per_m3 *
                    d_composition(
                        component,
                        column);
            if (!std::isfinite(derivative)) {
                throw std::range_error(
                    "mpmc::flow_discretization_petsc: component molar-content derivative is non-finite");
            }
            result.component_molar_content_jacobian[
                component * q +
                column] =
                derivative;
        }
    }

    result.advective_energy_density_j_per_m3 =
        mass_density_kg_per_m3 *
        specific_enthalpy_j_per_kg;
    if (!std::isfinite(
            result.advective_energy_density_j_per_m3)) {
        throw std::range_error(
            "mpmc::flow_discretization_petsc: active phase advective energy density is non-finite");
    }

    result.advective_energy_density_gradient
        .resize(q);
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        const double derivative =
            mass_density_gradient[column] *
                specific_enthalpy_j_per_kg +
            mass_density_kg_per_m3 *
                specific_enthalpy_gradient[
                    column];
        if (!std::isfinite(derivative)) {
            throw std::range_error(
                "mpmc::flow_discretization_petsc: active phase energy-density derivative is non-finite");
        }
        result.advective_energy_density_gradient[
            column] =
            derivative;
    }

    return result;
}

[[nodiscard]] inline ActivePhaseAdvectivePayload3D
extract_active_phase(
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        cell,
    std::size_t active_phase_index) {
    return std::visit(
        [active_phase_index](const auto& typed)
            -> ActivePhaseAdvectivePayload3D {
            using Typed =
                std::remove_cvref_t<
                    decltype(typed)>;

            if constexpr (
                std::is_same_v<
                    Typed,
                    SinglePhaseCurrentCellLinearization3D>) {
                if (active_phase_index != 0U) {
                    throw std::out_of_range(
                        "single-phase active phase index out of range");
                }
                const auto mobility =
                    mpmc::flow::
                        build_single_phase_mobility_linearization(
                            typed.state,
                            typed.transport);
                const auto composition =
                    typed.state.phase_composition();
                const std::size_t q =
                    typed.state.layout()
                        .unknown_count();
                return make_active_payload(
                    mobility.phase_pressure_pa,
                    mobility.mass_density_kg_per_m3,
                    mobility.mobility_per_pa_s,
                    mobility.phase_pressure_gradient,
                    mobility.mass_density_gradient,
                    mobility.mobility_gradient,
                    typed.molar_density
                        .molar_density_mol_per_m3,
                    typed.molar_density.gradient,
                    composition,
                    typed.caloric
                        .specific_enthalpy_j_per_kg,
                    typed.caloric
                        .specific_enthalpy_gradient,
                    q,
                    [&typed](
                        std::size_t component,
                        std::size_t column) {
                        return mpmc::flow::
                            single_phase_detail::
                                d_composition(
                                    typed.state.layout(),
                                    component,
                                    column);
                    });
            } else if constexpr (
                std::is_same_v<
                    Typed,
                    TwoPhaseCurrentCellLinearization3D>) {
                if (active_phase_index >= 2U) {
                    throw std::out_of_range(
                        "two-phase active phase index out of range");
                }
                const auto mobility =
                    mpmc::flow::
                        build_two_phase_mobility_linearization(
                            typed.state,
                            typed.transport);
                const auto composition =
                    typed.state.phase_composition(
                        active_phase_index);
                const std::size_t q =
                    typed.state.layout()
                        .unknown_count();
                return make_active_payload(
                    mobility.phase_pressure_pa[
                        active_phase_index],
                    mobility.mass_density_kg_per_m3[
                        active_phase_index],
                    mobility.mobility_per_pa_s[
                        active_phase_index],
                    mobility.phase_pressure_gradient[
                        active_phase_index],
                    mobility.mass_density_gradient[
                        active_phase_index],
                    mobility.mobility_gradient[
                        active_phase_index],
                    typed.molar_density
                        .molar_density_mol_per_m3[
                            active_phase_index],
                    typed.molar_density.gradient[
                        active_phase_index],
                    composition,
                    typed.caloric
                        .specific_enthalpy_j_per_kg[
                            active_phase_index],
                    typed.caloric
                        .specific_enthalpy_gradient[
                            active_phase_index],
                    q,
                    [&typed, active_phase_index](
                        std::size_t component,
                        std::size_t column) {
                        return mpmc::flow::
                            two_phase_detail::
                                d_composition(
                                    typed.state.layout(),
                                    active_phase_index,
                                    component,
                                    column);
                    });
            } else {
                static_assert(
                    std::is_same_v<
                        Typed,
                        FixedThreePhaseCurrentCellLinearization3D>);
                if (active_phase_index >=
                    mpmc::flow::
                        fixed_three_phase_count) {
                    throw std::out_of_range(
                        "three-phase active phase index out of range");
                }
                const auto mobility =
                    mpmc::flow::
                        build_local_phase_mobility_linearization(
                            typed.state,
                            typed.transport,
                            typed
                                .saturation_constitutive);
                const auto slot =
                    static_cast<
                        mpmc::flow::PhaseSlot3>(
                            active_phase_index);
                const auto composition =
                    typed.state.phase_composition(
                        slot);
                const std::size_t q =
                    typed.state.layout()
                        .unknown_count();
                return make_active_payload(
                    mobility.phase_pressure_pa[
                        active_phase_index],
                    mobility.mass_density_kg_per_m3[
                        active_phase_index],
                    mobility.mobility_per_pa_s[
                        active_phase_index],
                    mobility.phase_pressure_gradient[
                        active_phase_index],
                    mobility.mass_density_gradient[
                        active_phase_index],
                    mobility.mobility_gradient[
                        active_phase_index],
                    typed.molar_density
                        .molar_density_mol_per_m3[
                            active_phase_index],
                    typed.molar_density.gradient[
                        active_phase_index],
                    composition,
                    typed.caloric
                        .specific_enthalpy_j_per_kg[
                            active_phase_index],
                    typed.caloric
                        .specific_enthalpy_gradient[
                            active_phase_index],
                    q,
                    [&typed, active_phase_index](
                        std::size_t component,
                        std::size_t column) {
                        return mpmc::flow::
                            component_accumulation_detail::
                                d_composition(
                                    typed.state.layout(),
                                    active_phase_index,
                                    component,
                                    column);
                    });
            }
        },
        cell);
}

[[nodiscard]] inline PetscErrorCode
build_potential_side(
    const MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D&
        face_input,
    const mpmc::flow::
        CrossCardinalityFacePhaseIdentityBinding&
            phase_binding,
    CrossCardinalityAbsentPhaseSide3D side,
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        owner,
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        neighbour,
    const CrossCardinalityTpfaBridgeBinding3D&
        bridge,
    PotentialSidePayload3D* output,
    NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    const auto& cell =
        side ==
                CrossCardinalityAbsentPhaseSide3D::
                    owner
            ? owner
            : neighbour;
    const auto active_index =
        side ==
                CrossCardinalityAbsentPhaseSide3D::
                    owner
            ? phase_binding
                  .owner_active_phase_index
            : phase_binding
                  .neighbour_active_phase_index;
    const std::size_t q =
        state_identity(cell)
            .layout.unknown_count();

    if (active_index.has_value()) {
        auto active =
            extract_active_phase(
                cell,
                *active_index);
        output->active = true;
        output->phase_pressure_pa =
            active.phase_pressure_pa;
        output->mass_density_kg_per_m3 =
            active.mass_density_kg_per_m3;
        output->mobility_per_pa_s =
            active.mobility_per_pa_s;
        output->phase_pressure_gradient =
            active.phase_pressure_gradient;
        output->mass_density_gradient =
            active.mass_density_gradient;
        output->mobility_gradient =
            active.mobility_gradient;
        output->active_advective.emplace(
            std::move(active));
        return PETSC_SUCCESS;
    }

    if (bridge.absent_phase_extension_evaluator ==
        nullptr) {
        return PETSC_ERR_SUP;
    }

    std::optional<
        mpmc::flow::
            AbsentPhasePotentialExtensionLinearization>
        extension;
    NaturalVariableSnesEvaluationStatus3D
        extension_status =
            NaturalVariableSnesEvaluationStatus3D::
                success;
    const PetscErrorCode error =
        bridge.absent_phase_extension_evaluator(
            face_input,
            phase_binding,
            side,
            owner,
            neighbour,
            bridge.user_context,
            &extension,
            &extension_status);
    if (error != PETSC_SUCCESS) {
        return error;
    }
    if (extension_status ==
        NaturalVariableSnesEvaluationStatus3D::
            domain_error) {
        *status = extension_status;
        return PETSC_SUCCESS;
    }
    if (extension_status !=
            NaturalVariableSnesEvaluationStatus3D::
                success ||
        !extension.has_value() ||
        extension->identity() !=
            phase_binding.identity ||
        extension->input_count() !=
            q) {
        return PETSC_ERR_ARG_INCOMP;
    }

    output->active = false;
    output->phase_pressure_pa =
        extension->phase_pressure_pa();
    output->mass_density_kg_per_m3 =
        extension->mass_density_kg_per_m3();
    output->mobility_per_pa_s =
        mpmc::flow::
            CrossCardinalityAbsentPhaseSemantics::
                mobility_per_pa_s;
    output->phase_pressure_gradient.assign(
        extension
            ->phase_pressure_gradient()
            .begin(),
        extension
            ->phase_pressure_gradient()
            .end());
    output->mass_density_gradient.assign(
        extension
            ->mass_density_gradient()
            .begin(),
        extension
            ->mass_density_gradient()
            .end());
    output->mobility_gradient.assign(
        q,
        0.0);
    output->active_advective.reset();
    return PETSC_SUCCESS;
}

[[nodiscard]] inline bool
selects_owner(
    mpmc::flow::UpwindCellSelection3P
        selection) noexcept {
    return selection ==
               mpmc::flow::
                   UpwindCellSelection3P::
                       owner_negative_phase_potential ||
        selection ==
               mpmc::flow::
                   UpwindCellSelection3P::
                       owner_exact_zero_tie;
}

[[nodiscard]] inline
mpmc::flow::UpwindCellSelection3P
select_upwind(
    double phase_potential_difference_pa) {
    if (phase_potential_difference_pa <
        0.0) {
        return mpmc::flow::
            UpwindCellSelection3P::
                owner_negative_phase_potential;
    }
    if (phase_potential_difference_pa >
        0.0) {
        return mpmc::flow::
            UpwindCellSelection3P::
                neighbour_positive_phase_potential;
    }
    return mpmc::flow::
        UpwindCellSelection3P::
            owner_exact_zero_tie;
}

inline void normalize_component_result(
    mpmc::flow_discretization::
        NormalizedComponentFaceContributionLinearization3D*
            result,
    std::span<const double> rate,
    std::span<const double> owner_rate_jacobian,
    std::span<const double> neighbour_rate_jacobian,
    mpmc::flow_discretization::
        TwoCellBulkVolume3D bulk_volume,
    std::size_t owner_q,
    std::size_t neighbour_q) {
    if (result == nullptr) {
        throw std::invalid_argument(
            "null cross-cardinality component result");
    }
    const std::size_t n =
        rate.size();
    if (owner_rate_jacobian.size() !=
            n * owner_q ||
        neighbour_rate_jacobian.size() !=
            n * neighbour_q) {
        throw std::invalid_argument(
            "cross-cardinality component rate Jacobian shape mismatch");
    }

    result
        ->owner_component_contribution_mol_per_bulk_m3_s
        .resize(n);
    result
        ->neighbour_component_contribution_mol_per_bulk_m3_s
        .resize(n);
    result
        ->owner_row_owner_column_jacobian
        .resize(n * owner_q);
    result
        ->owner_row_neighbour_column_jacobian
        .resize(n * neighbour_q);
    result
        ->neighbour_row_owner_column_jacobian
        .resize(n * owner_q);
    result
        ->neighbour_row_neighbour_column_jacobian
        .resize(n * neighbour_q);

    result->owner_total_contribution_mol_per_bulk_m3_s =
        0.0;
    result->neighbour_total_contribution_mol_per_bulk_m3_s =
        0.0;
    result
        ->owner_total_row_owner_column_gradient
        .assign(owner_q, 0.0);
    result
        ->owner_total_row_neighbour_column_gradient
        .assign(neighbour_q, 0.0);
    result
        ->neighbour_total_row_owner_column_gradient
        .assign(owner_q, 0.0);
    result
        ->neighbour_total_row_neighbour_column_gradient
        .assign(neighbour_q, 0.0);

    for (std::size_t component = 0U;
         component < n;
         ++component) {
        const double value =
            rate[component];
        result
            ->owner_component_contribution_mol_per_bulk_m3_s[
                component] =
            value /
            bulk_volume.owner_bulk_volume_m3;
        result
            ->neighbour_component_contribution_mol_per_bulk_m3_s[
                component] =
            -value /
            bulk_volume.neighbour_bulk_volume_m3;

        result->owner_total_contribution_mol_per_bulk_m3_s +=
            value /
            bulk_volume.owner_bulk_volume_m3;
        result->neighbour_total_contribution_mol_per_bulk_m3_s -=
            value /
            bulk_volume.neighbour_bulk_volume_m3;

        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            const double derivative =
                owner_rate_jacobian[
                    component * owner_q +
                    column];
            result
                ->owner_row_owner_column_jacobian[
                    component * owner_q +
                    column] =
                derivative /
                bulk_volume.owner_bulk_volume_m3;
            result
                ->neighbour_row_owner_column_jacobian[
                    component * owner_q +
                    column] =
                -derivative /
                bulk_volume.neighbour_bulk_volume_m3;
            result
                ->owner_total_row_owner_column_gradient[
                    column] +=
                derivative /
                bulk_volume.owner_bulk_volume_m3;
            result
                ->neighbour_total_row_owner_column_gradient[
                    column] -=
                derivative /
                bulk_volume.neighbour_bulk_volume_m3;
        }
        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            const double derivative =
                neighbour_rate_jacobian[
                    component *
                        neighbour_q +
                    column];
            result
                ->owner_row_neighbour_column_jacobian[
                    component *
                        neighbour_q +
                    column] =
                derivative /
                bulk_volume.owner_bulk_volume_m3;
            result
                ->neighbour_row_neighbour_column_jacobian[
                    component *
                        neighbour_q +
                    column] =
                -derivative /
                bulk_volume.neighbour_bulk_volume_m3;
            result
                ->owner_total_row_neighbour_column_gradient[
                    column] +=
                derivative /
                bulk_volume.owner_bulk_volume_m3;
            result
                ->neighbour_total_row_neighbour_column_gradient[
                    column] -=
                derivative /
                bulk_volume.neighbour_bulk_volume_m3;
        }
    }
}

inline void normalize_energy_result(
    mpmc::flow_discretization::
        NormalizedEnergyFaceContributionLinearization3D*
            result,
    double energy_rate_w,
    std::span<const double>
        owner_energy_rate_gradient,
    std::span<const double>
        neighbour_energy_rate_gradient,
    mpmc::flow_discretization::
        TwoCellBulkVolume3D bulk_volume) {
    if (result == nullptr ||
        owner_energy_rate_gradient.size() !=
            result->owner_state_identity
                .layout.unknown_count() ||
        neighbour_energy_rate_gradient.size() !=
            result->neighbour_state_identity
                .layout.unknown_count()) {
        throw std::invalid_argument(
            "cross-cardinality energy rate Jacobian shape mismatch");
    }

    result->owner_contribution_w_per_bulk_m3 =
        energy_rate_w /
        bulk_volume.owner_bulk_volume_m3;
    result->neighbour_contribution_w_per_bulk_m3 =
        -energy_rate_w /
        bulk_volume.neighbour_bulk_volume_m3;

    result->owner_row_owner_column_gradient.assign(
        owner_energy_rate_gradient.begin(),
        owner_energy_rate_gradient.end());
    result->neighbour_row_owner_column_gradient =
        result->owner_row_owner_column_gradient;
    for (double& value :
         result->owner_row_owner_column_gradient) {
        value /=
            bulk_volume.owner_bulk_volume_m3;
    }
    for (double& value :
         result->neighbour_row_owner_column_gradient) {
        value *=
            -1.0 /
            bulk_volume.neighbour_bulk_volume_m3;
    }

    result->owner_row_neighbour_column_gradient.assign(
        neighbour_energy_rate_gradient.begin(),
        neighbour_energy_rate_gradient.end());
    result->neighbour_row_neighbour_column_gradient =
        result->owner_row_neighbour_column_gradient;
    for (double& value :
         result->owner_row_neighbour_column_gradient) {
        value /=
            bulk_volume.owner_bulk_volume_m3;
    }
    for (double& value :
         result->neighbour_row_neighbour_column_gradient) {
        value *=
            -1.0 /
            bulk_volume.neighbour_bulk_volume_m3;
    }
}

} // namespace cross_cardinality_tpfa_detail

/// Standard physical bridge for a face whose two cells do not share a
/// slot-aligned active phase set.
///
/// For every explicit physical phase identity in the pairing plan:
/// - active-side pressure/density/mobility and advective payload come from the
///   already-evaluated local physical closure;
/// - an inactive side contributes structural-zero mobility and must obtain only
///   hypothetical pressure/density (+ derivatives) from the explicit extension
///   evaluator;
/// - if the inactive side is selected as upwind, phase flux and its Jacobian are
///   exactly zero, so composition/enthalpy are never requested for that side;
/// - if the active side is selected as upwind, component and energy advection
///   use only that active side's real molar density/composition/enthalpy.
///
/// The face-density policy intentionally matches the existing fixed-cardinality
/// v1 kernels: arithmetic owner/neighbour mean. Static transmissibility and
/// thermal conductance remain frozen.
inline PetscErrorCode
evaluate_standard_cross_cardinality_tpfa_face_3d(
    const MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D&
        face_input,
    const mpmc::flow::
        CrossCardinalityFacePhaseIdentityPlan&
            phase_identity_plan,
    mpmc::flow_discretization::
        TwoCellBulkVolume3D bulk_volume,
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        owner,
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        neighbour,
    void* raw_context,
    std::optional<
        MixedCardinalityPhysicalFaceLinearization3D>*
            output,
    NaturalVariableSnesEvaluationStatus3D*
        status) {
    using namespace
        cross_cardinality_tpfa_detail;

    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    *status =
        NaturalVariableSnesEvaluationStatus3D::
            success;

    auto* bridge =
        static_cast<
            CrossCardinalityTpfaBridgeBinding3D*>(
                raw_context);

    try {
        mpmc::flow_discretization::detail::
            validate_materialized_entry(
                face_input.transmissibility);
        mpmc::flow::
            phase_potential_upwind_detail::
                require_finite_geometry(
                    face_input.gravity,
                    face_input
                        .owner_to_neighbour_displacement);

        const auto& owner_identity =
            state_identity(owner);
        const auto& neighbour_identity =
            state_identity(neighbour);
        const std::size_t owner_q =
            owner_identity.layout
                .unknown_count();
        const std::size_t neighbour_q =
            neighbour_identity.layout
                .unknown_count();
        const std::size_t n =
            owner_identity.layout
                .component_count();

        if (owner_identity.component_ids !=
                neighbour_identity.component_ids ||
            n == 0U ||
            phase_identity_plan
                    .owner_phase_count() !=
                owner_identity.layout
                    .phase_count() ||
            phase_identity_plan
                    .neighbour_phase_count() !=
                neighbour_identity.layout
                    .phase_count() ||
            !std::isfinite(
                bulk_volume.owner_bulk_volume_m3) ||
            !(bulk_volume.owner_bulk_volume_m3 >
              0.0) ||
            !std::isfinite(
                bulk_volume
                    .neighbour_bulk_volume_m3) ||
            !(bulk_volume
                  .neighbour_bulk_volume_m3 >
              0.0) ||
            !std::isfinite(
                face_input.thermal_conductance
                    .conductance_w_per_k) ||
            face_input.thermal_conductance
                    .conductance_w_per_k <
                0.0) {
            throw std::invalid_argument(
                "mpmc::flow_discretization_petsc: malformed cross-cardinality TPFA face metadata");
        }

        const double tf =
            face_input.transmissibility
                .static_transmissibility
                ->face_transmissibility_m3;
        const double gravity_projection =
            mpmc::flow::
                phase_potential_upwind_detail::
                    dot(
                        face_input.gravity,
                        face_input
                            .owner_to_neighbour_displacement);

        std::vector<double>
            component_rate(
                n,
                0.0);
        std::vector<double>
            owner_component_rate_jacobian(
                n * owner_q,
                0.0);
        std::vector<double>
            neighbour_component_rate_jacobian(
                n * neighbour_q,
                0.0);
        double advective_energy_rate = 0.0;
        std::vector<double>
            owner_energy_rate_gradient(
                owner_q,
                0.0);
        std::vector<double>
            neighbour_energy_rate_gradient(
                neighbour_q,
                0.0);

        for (const auto& binding :
             phase_identity_plan.bindings()) {
            PotentialSidePayload3D
                owner_side;
            PotentialSidePayload3D
                neighbour_side;

            PetscErrorCode error =
                build_potential_side(
                    face_input,
                    binding,
                    CrossCardinalityAbsentPhaseSide3D::
                        owner,
                    owner,
                    neighbour,
                    *bridge,
                    &owner_side,
                    status);
            if (error != PETSC_SUCCESS ||
                *status ==
                    NaturalVariableSnesEvaluationStatus3D::
                        domain_error) {
                return error;
            }
            error =
                build_potential_side(
                    face_input,
                    binding,
                    CrossCardinalityAbsentPhaseSide3D::
                        neighbour,
                    owner,
                    neighbour,
                    *bridge,
                    &neighbour_side,
                    status);
            if (error != PETSC_SUCCESS ||
                *status ==
                    NaturalVariableSnesEvaluationStatus3D::
                        domain_error) {
                return error;
            }

            const double face_density =
                0.5 *
                (owner_side
                     .mass_density_kg_per_m3 +
                 neighbour_side
                     .mass_density_kg_per_m3);
            const double delta =
                (neighbour_side.phase_pressure_pa -
                 owner_side.phase_pressure_pa) -
                face_density *
                    gravity_projection;
            if (!std::isfinite(face_density) ||
                !(face_density > 0.0) ||
                !std::isfinite(delta)) {
                throw std::range_error(
                    "mpmc::flow_discretization_petsc: cross-cardinality phase potential is invalid");
            }

            std::vector<double>
                owner_delta_gradient(
                    owner_q,
                    0.0);
            std::vector<double>
                neighbour_delta_gradient(
                    neighbour_q,
                    0.0);
            for (std::size_t column = 0U;
                 column < owner_q;
                 ++column) {
                owner_delta_gradient[column] =
                    -owner_side
                         .phase_pressure_gradient[
                             column] -
                    0.5 *
                        owner_side
                            .mass_density_gradient[
                                column] *
                        gravity_projection;
            }
            for (std::size_t column = 0U;
                 column < neighbour_q;
                 ++column) {
                neighbour_delta_gradient[column] =
                    neighbour_side
                        .phase_pressure_gradient[
                            column] -
                    0.5 *
                        neighbour_side
                            .mass_density_gradient[
                                column] *
                        gravity_projection;
            }

            const auto selection =
                select_upwind(delta);
            const bool owner_upwind =
                selects_owner(selection);
            const auto& upwind_side =
                owner_upwind
                    ? owner_side
                    : neighbour_side;

            const double mobility =
                upwind_side
                    .mobility_per_pa_s;
            const double flux =
                -tf *
                mobility *
                delta;
            if (!std::isfinite(flux)) {
                throw std::range_error(
                    "mpmc::flow_discretization_petsc: cross-cardinality Darcy phase flux is non-finite");
            }

            std::vector<double>
                owner_flux_gradient(
                    owner_q,
                    0.0);
            std::vector<double>
                neighbour_flux_gradient(
                    neighbour_q,
                    0.0);
            for (std::size_t column = 0U;
                 column < owner_q;
                 ++column) {
                const double dlambda =
                    owner_upwind
                        ? owner_side
                              .mobility_gradient[
                                  column]
                        : 0.0;
                owner_flux_gradient[column] =
                    -tf *
                    (mobility *
                         owner_delta_gradient[
                             column] +
                     delta * dlambda);
            }
            for (std::size_t column = 0U;
                 column < neighbour_q;
                 ++column) {
                const double dlambda =
                    owner_upwind
                        ? 0.0
                        : neighbour_side
                              .mobility_gradient[
                                  column];
                neighbour_flux_gradient[column] =
                    -tf *
                    (mobility *
                         neighbour_delta_gradient[
                             column] +
                     delta * dlambda);
            }

            if (!upwind_side.active) {
                if (mobility != 0.0 ||
                    flux != 0.0 ||
                    std::any_of(
                        owner_flux_gradient.begin(),
                        owner_flux_gradient.end(),
                        [](double value) {
                            return value != 0.0;
                        }) ||
                    std::any_of(
                        neighbour_flux_gradient.begin(),
                        neighbour_flux_gradient.end(),
                        [](double value) {
                            return value != 0.0;
                        })) {
                    throw std::runtime_error(
                        "mpmc::flow_discretization_petsc: inactive upwind phase violated structural-zero flux contract");
                }
                continue;
            }

            if (!upwind_side
                     .active_advective
                     .has_value()) {
                throw std::runtime_error(
                    "mpmc::flow_discretization_petsc: active upwind phase lacks advective payload");
            }
            const auto& advective =
                *upwind_side.active_advective;

            for (std::size_t component = 0U;
                 component < n;
                 ++component) {
                const double content =
                    advective
                        .component_molar_content_mol_per_m3[
                            component];
                component_rate[component] +=
                    content * flux;

                for (std::size_t column = 0U;
                     column < owner_q;
                     ++column) {
                    const double dcontent =
                        owner_upwind
                            ? advective
                                  .component_molar_content_jacobian[
                                      component *
                                          owner_q +
                                      column]
                            : 0.0;
                    owner_component_rate_jacobian[
                        component * owner_q +
                        column] +=
                        dcontent * flux +
                        content *
                            owner_flux_gradient[
                                column];
                }
                for (std::size_t column = 0U;
                     column < neighbour_q;
                     ++column) {
                    const double dcontent =
                        owner_upwind
                            ? 0.0
                            : advective
                                  .component_molar_content_jacobian[
                                      component *
                                          neighbour_q +
                                      column];
                    neighbour_component_rate_jacobian[
                        component *
                            neighbour_q +
                        column] +=
                        dcontent * flux +
                        content *
                            neighbour_flux_gradient[
                                column];
                }
            }

            const double energy_density =
                advective
                    .advective_energy_density_j_per_m3;
            advective_energy_rate +=
                energy_density *
                flux;
            for (std::size_t column = 0U;
                 column < owner_q;
                 ++column) {
                const double d_energy_density =
                    owner_upwind
                        ? advective
                              .advective_energy_density_gradient[
                                  column]
                        : 0.0;
                owner_energy_rate_gradient[
                    column] +=
                    d_energy_density *
                        flux +
                    energy_density *
                        owner_flux_gradient[
                            column];
            }
            for (std::size_t column = 0U;
                 column < neighbour_q;
                 ++column) {
                const double d_energy_density =
                    owner_upwind
                        ? 0.0
                        : advective
                              .advective_energy_density_gradient[
                                  column];
                neighbour_energy_rate_gradient[
                    column] +=
                    d_energy_density *
                        flux +
                    energy_density *
                        neighbour_flux_gradient[
                            column];
            }
        }

        const double thermal_conductance =
            face_input.thermal_conductance
                .conductance_w_per_k;
        const double conductive_energy_rate =
            thermal_conductance *
            (owner_identity.temperature_k -
             neighbour_identity.temperature_k);
        const double total_energy_rate =
            advective_energy_rate +
            conductive_energy_rate;
        owner_energy_rate_gradient[
            owner_identity.layout
                .temperature_unknown_index()] +=
            thermal_conductance;
        neighbour_energy_rate_gradient[
            neighbour_identity.layout
                .temperature_unknown_index()] -=
            thermal_conductance;

        MixedCardinalityPhysicalFaceLinearization3D
            result;
        result.component.face =
            face_input.face;
        result.component.bulk_volume =
            bulk_volume;
        result.component.component_ids =
            owner_identity.component_ids;
        result.component.owner_state_identity =
            owner_identity;
        result.component.neighbour_state_identity =
            neighbour_identity;

        normalize_component_result(
            &result.component,
            component_rate,
            owner_component_rate_jacobian,
            neighbour_component_rate_jacobian,
            bulk_volume,
            owner_q,
            neighbour_q);

        result.energy.face =
            face_input.face;
        result.energy.bulk_volume = {
            bulk_volume.owner_bulk_volume_m3,
            bulk_volume.neighbour_bulk_volume_m3};
        result.energy.owner_state_identity =
            owner_identity;
        result.energy.neighbour_state_identity =
            neighbour_identity;

        normalize_energy_result(
            &result.energy,
            total_energy_rate,
            owner_energy_rate_gradient,
            neighbour_energy_rate_gradient,
            bulk_volume);

        output->emplace(
            std::move(result));
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        *status =
            NaturalVariableSnesEvaluationStatus3D::
                domain_error;
        return PETSC_SUCCESS;
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_CROSS_CARDINALITY_TPFA_BRIDGE_HPP
