#ifndef MPMC_FLOW_DISCRETIZATION_PETSC_SW92_TRANSACTIONAL_PHASE_TRANSITION_RESTART_HPP
#define MPMC_FLOW_DISCRETIZATION_PETSC_SW92_TRANSACTIONAL_PHASE_TRANSITION_RESTART_HPP

#include <mpmc/flow_discretization_petsc/cell_scoped_mixed_cardinality_evaluator_dispatcher.hpp>
#include <mpmc/flow_discretization_petsc/post_snes_sw92_profile_c_phase_transition_scanner.hpp>
#include <mpmc/flow_discretization_petsc/sw92_production_cell_evaluator.hpp>

#include <petscvec.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <type_traits>
#include <variant>
#include <vector>

namespace mpmc::flow_discretization_petsc {

inline constexpr std::string_view
    sw92_transactional_phase_transition_restart_convention =
        "flow_discretization_petsc/sw92-transactional-phase-transition-restart/no-cross-cardinality-face/v1";

using Sw92TransitionTargetIdentityResolver3D =
    PetscErrorCode (*)(
        mpmc::mesh::GlobalEntityId cell_global,
        const Sw92AuthoritativeTargetMaterialization3D&
            target,
        const mpmc::flow::FrozenActivePhaseIdentityMap&
            source_active_phases,
        void* user_context,
        std::optional<
            mpmc::flow::FrozenActivePhaseIdentityMap>*
                target_active_phases);

struct Sw92TransitionTargetIdentityBinding3D {
    Sw92TransitionTargetIdentityResolver3D resolver{};
    void* user_context{};
};

struct Sw92TransactionalRebuildBaselineCell3D {
    mpmc::mesh::LocalIndex cell{
        mpmc::mesh::LocalIndex::value_type{0}};
    mpmc::mesh::GlobalEntityId cell_global{
        mpmc::mesh::GlobalEntityId::value_type{0}};
    double bulk_volume_m3{};
    double porosity{};
    std::vector<std::string> component_ids;
    std::optional<
        mpmc::flow::
            PoreVolumeComponentAccumulationSnapshot3P>
        previous_component_accumulation;
    std::optional<
        mpmc::flow::
            PoreVolumeEnergyAccumulationSnapshot3P>
        previous_energy_accumulation;
};

namespace sw92_transactional_restart_detail {

[[nodiscard]] inline bool near_roundoff(
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
            std::numeric_limits<double>::epsilon() *
            scale;
}

[[nodiscard]] inline const
mpmc::flow::NaturalVariableStateIdentity3P&
state_identity(
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        current) {
    return std::visit(
        [](const auto& typed)
            -> const mpmc::flow::
                NaturalVariableStateIdentity3P& {
            return typed.transport.state_identity;
        },
        current);
}

[[nodiscard]] inline
mpmc::flow::PoreVolumeComponentAccumulationSnapshot3P
current_component_inventory(
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        current,
    double porosity) {
    return std::visit(
        [porosity](const auto& typed)
            -> mpmc::flow::
                PoreVolumeComponentAccumulationSnapshot3P {
            using Typed =
                std::decay_t<decltype(typed)>;
            if constexpr (
                std::is_same_v<
                    Typed,
                    SinglePhaseCurrentCellLinearization3D>) {
                return mpmc::flow::
                    build_single_phase_component_accumulation(
                        typed.state,
                        porosity);
            } else if constexpr (
                std::is_same_v<
                    Typed,
                    TwoPhaseCurrentCellLinearization3D>) {
                return mpmc::flow::
                    build_two_phase_component_accumulation(
                        typed.state,
                        porosity);
            } else {
                return mpmc::flow::
                    build_pore_volume_component_accumulation(
                        typed.state,
                        porosity);
            }
        },
        current);
}


[[nodiscard]] inline
mpmc::flow::PoreVolumeEnergyAccumulationSnapshot3P
current_energy_inventory(
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        current,
    double porosity) {
    return std::visit(
        [porosity](const auto& typed)
            -> mpmc::flow::
                PoreVolumeEnergyAccumulationSnapshot3P {
            using Typed =
                std::decay_t<decltype(typed)>;
            if constexpr (
                std::is_same_v<
                    Typed,
                    SinglePhaseCurrentCellLinearization3D>) {
                return mpmc::flow::
                    build_single_phase_energy_accumulation_snapshot(
                        typed.state,
                        porosity,
                        typed.transport,
                        typed.caloric,
                        typed.rock);
            } else if constexpr (
                std::is_same_v<
                    Typed,
                    TwoPhaseCurrentCellLinearization3D>) {
                return mpmc::flow::
                    build_two_phase_energy_accumulation_snapshot(
                        typed.state,
                        porosity,
                        typed.transport,
                        typed.caloric,
                        typed.rock);
            } else {
                return mpmc::flow::
                    build_pore_volume_energy_accumulation_snapshot(
                        typed.state,
                        porosity,
                        typed.transport,
                        typed.caloric,
                        typed.rock);
            }
        },
        current);
}

struct AggregateStorageLinearization3D {
    double total_molar_accumulation{};
    double total_internal_energy{};
    double d_molar_dp{};
    double d_molar_dt{};
    double d_energy_dp{};
    double d_energy_dt{};
};

[[nodiscard]] inline
AggregateStorageLinearization3D
aggregate_storage_linearization(
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        current,
    double porosity) {
    return std::visit(
        [porosity](const auto& typed)
            -> AggregateStorageLinearization3D {
            using Typed =
                std::decay_t<decltype(typed)>;

            mpmc::flow::
                PoreVolumeComponentAccumulationLinearization3P
                component;
            mpmc::flow::
                PoreVolumeEnergyAccumulationLinearization3P
                energy;

            if constexpr (
                std::is_same_v<
                    Typed,
                    SinglePhaseCurrentCellLinearization3D>) {
                component =
                    mpmc::flow::
                        build_single_phase_component_accumulation_linearization(
                            typed.state,
                            porosity,
                            typed.molar_density);
                energy =
                    mpmc::flow::
                        build_single_phase_energy_accumulation_linearization(
                            typed.state,
                            porosity,
                            typed.transport,
                            typed.caloric,
                            typed.rock);
            } else if constexpr (
                std::is_same_v<
                    Typed,
                    TwoPhaseCurrentCellLinearization3D>) {
                component =
                    mpmc::flow::
                        build_two_phase_component_accumulation_linearization(
                            typed.state,
                            porosity,
                            typed.molar_density);
                energy =
                    mpmc::flow::
                        build_two_phase_energy_accumulation_linearization(
                            typed.state,
                            porosity,
                            typed.transport,
                            typed.caloric,
                            typed.rock);
            } else {
                component =
                    mpmc::flow::
                        build_pore_volume_component_accumulation_linearization(
                            typed.state,
                            porosity,
                            typed.molar_density);
                energy =
                    mpmc::flow::
                        build_pore_volume_energy_accumulation_linearization(
                            typed.state,
                            porosity,
                            typed.transport,
                            typed.caloric,
                            typed.rock);
            }

            if (component.total_accumulation_gradient.size() < 2U ||
                energy.total_internal_energy_gradient.size() < 2U) {
                throw std::invalid_argument(
                    "SW92 transactional aggregate storage lacks p/T derivatives");
            }

            const auto component_snapshot =
                current_component_inventory(
                    MixedCardinalityPhysicalCurrentCellLinearization3D{
                        typed},
                    porosity);
            return {
                component_snapshot
                    .total_accumulation_mol_per_bulk_m3,
                energy
                    .total_internal_energy_j_per_bulk_m3,
                component
                    .total_accumulation_gradient[0],
                component
                    .total_accumulation_gradient[1],
                energy
                    .total_internal_energy_gradient[0],
                energy
                    .total_internal_energy_gradient[1]};
        },
        current);
}

struct ConservationStorageLinearization3D {
    mpmc::flow::PoreVolumeComponentAccumulationSnapshot3P
        component_snapshot;
    mpmc::flow::PoreVolumeEnergyAccumulationSnapshot3P
        energy_snapshot;
    mpmc::flow::PoreVolumeComponentAccumulationLinearization3P
        component_linearization;
    mpmc::flow::PoreVolumeEnergyAccumulationLinearization3P
        energy_linearization;
};

[[nodiscard]] inline
ConservationStorageLinearization3D
conservation_storage_linearization(
    const MixedCardinalityPhysicalCurrentCellLinearization3D&
        current,
    double porosity) {
    return std::visit(
        [porosity](const auto& typed)
            -> ConservationStorageLinearization3D {
            using Typed =
                std::decay_t<decltype(typed)>;
            if constexpr (
                std::is_same_v<
                    Typed,
                    SinglePhaseCurrentCellLinearization3D>) {
                return {
                    mpmc::flow::
                        build_single_phase_component_accumulation(
                            typed.state,
                            porosity),
                    mpmc::flow::
                        build_single_phase_energy_accumulation_snapshot(
                            typed.state,
                            porosity,
                            typed.transport,
                            typed.caloric,
                            typed.rock),
                    mpmc::flow::
                        build_single_phase_component_accumulation_linearization(
                            typed.state,
                            porosity,
                            typed.molar_density),
                    mpmc::flow::
                        build_single_phase_energy_accumulation_linearization(
                            typed.state,
                            porosity,
                            typed.transport,
                            typed.caloric,
                            typed.rock)};
            } else if constexpr (
                std::is_same_v<
                    Typed,
                    TwoPhaseCurrentCellLinearization3D>) {
                return {
                    mpmc::flow::
                        build_two_phase_component_accumulation(
                            typed.state,
                            porosity),
                    mpmc::flow::
                        build_two_phase_energy_accumulation_snapshot(
                            typed.state,
                            porosity,
                            typed.transport,
                            typed.caloric,
                            typed.rock),
                    mpmc::flow::
                        build_two_phase_component_accumulation_linearization(
                            typed.state,
                            porosity,
                            typed.molar_density),
                    mpmc::flow::
                        build_two_phase_energy_accumulation_linearization(
                            typed.state,
                            porosity,
                            typed.transport,
                            typed.caloric,
                            typed.rock)};
            } else {
                return {
                    mpmc::flow::
                        build_pore_volume_component_accumulation(
                            typed.state,
                            porosity),
                    mpmc::flow::
                        build_pore_volume_energy_accumulation_snapshot(
                            typed.state,
                            porosity,
                            typed.transport,
                            typed.caloric,
                            typed.rock),
                    mpmc::flow::
                        build_pore_volume_component_accumulation_linearization(
                            typed.state,
                            porosity,
                            typed.molar_density),
                    mpmc::flow::
                        build_pore_volume_energy_accumulation_linearization(
                            typed.state,
                            porosity,
                            typed.transport,
                            typed.caloric,
                            typed.rock)};
            }
        },
        current);
}

inline bool solve_dense_partial_pivot(
    std::vector<double> matrix,
    std::vector<double> rhs,
    std::vector<double>* solution) {
    if (solution == nullptr ||
        rhs.empty() ||
        matrix.size() !=
            rhs.size() * rhs.size()) {
        return false;
    }
    const std::size_t n = rhs.size();
    double matrix_scale = 0.0;
    for (double value : matrix) {
        if (!std::isfinite(value)) {
            return false;
        }
        matrix_scale =
            std::max(
                matrix_scale,
                std::abs(value));
    }
    if (!(matrix_scale > 0.0) ||
        !std::isfinite(matrix_scale)) {
        return false;
    }

    for (std::size_t column = 0U;
         column < n;
         ++column) {
        std::size_t pivot = column;
        double pivot_value =
            std::abs(
                matrix[
                    column * n +
                    column]);
        for (std::size_t row =
                 column + 1U;
             row < n;
             ++row) {
            const double value =
                std::abs(
                    matrix[
                        row * n +
                        column]);
            if (value > pivot_value) {
                pivot = row;
                pivot_value = value;
            }
        }
        if (!std::isfinite(pivot_value) ||
            pivot_value <=
                16384.0 *
                    std::numeric_limits<double>::epsilon() *
                    matrix_scale) {
            return false;
        }
        if (pivot != column) {
            for (std::size_t entry = 0U;
                 entry < n;
                 ++entry) {
                std::swap(
                    matrix[
                        column * n +
                        entry],
                    matrix[
                        pivot * n +
                        entry]);
            }
            std::swap(
                rhs[column],
                rhs[pivot]);
        }

        const double diagonal =
            matrix[
                column * n +
                column];
        for (std::size_t row =
                 column + 1U;
             row < n;
             ++row) {
            const double factor =
                matrix[
                    row * n +
                    column] /
                diagonal;
            matrix[
                row * n +
                column] = 0.0;
            for (std::size_t entry =
                     column + 1U;
                 entry < n;
                 ++entry) {
                matrix[
                    row * n +
                    entry] -=
                    factor *
                    matrix[
                        column * n +
                        entry];
            }
            rhs[row] -=
                factor * rhs[column];
        }
    }

    solution->assign(
        n,
        0.0);
    for (std::size_t reverse = 0U;
         reverse < n;
         ++reverse) {
        const std::size_t row =
            n - 1U - reverse;
        double value = rhs[row];
        for (std::size_t column =
                 row + 1U;
             column < n;
             ++column) {
            value -=
                matrix[
                    row * n +
                    column] *
                (*solution)[column];
        }
        const double diagonal =
            matrix[
                row * n +
                row];
        if (!std::isfinite(diagonal) ||
            std::abs(diagonal) <=
                16384.0 *
                    std::numeric_limits<double>::epsilon() *
                    matrix_scale) {
            return false;
        }
        (*solution)[row] =
            value / diagonal;
        if (!std::isfinite(
                (*solution)[row])) {
            return false;
        }
    }
    return true;
}

template <typename Evaluate>
inline void conservative_storage_anchor(
    std::vector<double>* q,
    double porosity,
    const mpmc::flow::
        PoreVolumeComponentAccumulationSnapshot3P&
            target_component,
    const mpmc::flow::
        PoreVolumeEnergyAccumulationSnapshot3P&
            target_energy,
    Evaluate&& evaluate) {
    if (q == nullptr ||
        q->size() < 2U ||
        target_component.component_ids.size() < 2U ||
        target_component
                .component_accumulation_mol_per_bulk_m3
                .size() !=
            target_component.component_ids.size() ||
        !std::isfinite(
            target_energy
                .total_internal_energy_j_per_bulk_m3)) {
        throw std::invalid_argument(
            "SW92 transactional conservation anchor inputs are invalid");
    }

    const std::size_t n =
        target_component.component_ids.size();
    const std::size_t row_count =
        n + 1U;
    const std::size_t column_count =
        q->size();

    auto residual_and_jacobian =
        [&](std::span<const double> state,
            std::vector<double>* residual,
            std::vector<double>* jacobian,
            double* norm)
            -> bool {
        const auto current =
            evaluate(state);
        if (!current.has_value()) {
            return false;
        }
        const auto storage =
            conservation_storage_linearization(
                *current,
                porosity);
        if (storage.component_snapshot.component_ids !=
                target_component.component_ids ||
            storage.component_linearization.input_count !=
                column_count ||
            storage.energy_linearization
                    .total_internal_energy_gradient
                    .size() !=
                column_count) {
            return false;
        }

        residual->assign(
            row_count,
            0.0);
        jacobian->assign(
            row_count * column_count,
            0.0);
        double sum_square = 0.0;

        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const double scale =
                std::max(
                    1.0,
                    std::abs(
                        target_component
                            .component_accumulation_mol_per_bulk_m3[
                                component]));
            const double value =
                (storage.component_snapshot
                     .component_accumulation_mol_per_bulk_m3[
                         component] -
                 target_component
                     .component_accumulation_mol_per_bulk_m3[
                         component]) /
                scale;
            if (!std::isfinite(value)) {
                return false;
            }
            (*residual)[component] =
                value;
            sum_square +=
                value * value;
            for (std::size_t column = 0U;
                 column < column_count;
                 ++column) {
                (*jacobian)[
                    component * column_count +
                    column] =
                    storage.component_linearization
                        .d_component(
                            component,
                            column) /
                    scale;
            }
        }

        const double energy_scale =
            std::max(
                1.0,
                std::abs(
                    target_energy
                        .total_internal_energy_j_per_bulk_m3));
        const double energy_residual =
            (storage.energy_snapshot
                 .total_internal_energy_j_per_bulk_m3 -
             target_energy
                 .total_internal_energy_j_per_bulk_m3) /
            energy_scale;
        if (!std::isfinite(
                energy_residual)) {
            return false;
        }
        (*residual)[n] =
            energy_residual;
        sum_square +=
            energy_residual *
            energy_residual;
        for (std::size_t column = 0U;
             column < column_count;
             ++column) {
            (*jacobian)[
                n * column_count +
                column] =
                storage.energy_linearization
                    .total_internal_energy_gradient[
                        column] /
                energy_scale;
        }

        *norm =
            std::sqrt(
                sum_square);
        return std::isfinite(*norm);
    };

    constexpr std::size_t maximum_iterations = 8U;
    constexpr std::size_t maximum_backtracks = 12U;
    constexpr double target_norm = 1.0e-8;

    for (std::size_t iteration = 0U;
         iteration < maximum_iterations;
         ++iteration) {
        std::vector<double> residual;
        std::vector<double> jacobian;
        double norm = 0.0;
        if (!residual_and_jacobian(
                *q,
                &residual,
                &jacobian,
                &norm)) {
            throw std::range_error(
                "SW92 transactional conservation anchor cannot evaluate target state");
        }
        if (norm <= target_norm) {
            return;
        }

        std::vector<double> column_scale(
            column_count,
            1.0);
        column_scale[0] =
            std::max(
                1.0e5,
                0.1 *
                    std::abs(
                        (*q)[0]));
        column_scale[1] =
            std::max(
                10.0,
                0.1 *
                    std::abs(
                        (*q)[1]));

        std::vector<double> scaled_jacobian(
            jacobian.size(),
            0.0);
        for (std::size_t row = 0U;
             row < row_count;
             ++row) {
            for (std::size_t column = 0U;
                 column < column_count;
                 ++column) {
                scaled_jacobian[
                    row * column_count +
                    column] =
                    jacobian[
                        row * column_count +
                        column] *
                    column_scale[column];
            }
        }

        std::vector<double> normal(
            row_count * row_count,
            0.0);
        for (std::size_t row = 0U;
             row < row_count;
             ++row) {
            for (std::size_t other = 0U;
                 other < row_count;
                 ++other) {
                double value = 0.0;
                for (std::size_t column = 0U;
                     column < column_count;
                     ++column) {
                    value +=
                        scaled_jacobian[
                            row * column_count +
                            column] *
                        scaled_jacobian[
                            other * column_count +
                            column];
                }
                normal[
                    row * row_count +
                    other] =
                    value;
            }
        }

        double diagonal_scale = 0.0;
        for (std::size_t row = 0U;
             row < row_count;
             ++row) {
            diagonal_scale =
                std::max(
                    diagonal_scale,
                    std::abs(
                        normal[
                            row * row_count +
                            row]));
        }
        if (!(diagonal_scale > 0.0) ||
            !std::isfinite(
                diagonal_scale)) {
            throw std::range_error(
                "SW92 transactional conservation anchor normal matrix is empty");
        }
        const double regularization =
            1024.0 *
            std::numeric_limits<double>::epsilon() *
            diagonal_scale;
        for (std::size_t row = 0U;
             row < row_count;
             ++row) {
            normal[
                row * row_count +
                row] +=
                regularization;
        }

        std::vector<double> rhs(
            row_count,
            0.0);
        for (std::size_t row = 0U;
             row < row_count;
             ++row) {
            rhs[row] =
                -residual[row];
        }
        std::vector<double> dual;
        if (!solve_dense_partial_pivot(
                std::move(normal),
                std::move(rhs),
                &dual)) {
            throw std::range_error(
                "SW92 transactional conservation anchor normal solve failed");
        }

        std::vector<double> step(
            column_count,
            0.0);
        for (std::size_t column = 0U;
             column < column_count;
             ++column) {
            double normalized_step = 0.0;
            for (std::size_t row = 0U;
                 row < row_count;
                 ++row) {
                normalized_step +=
                    scaled_jacobian[
                        row * column_count +
                        column] *
                    dual[row];
            }
            step[column] =
                column_scale[column] *
                normalized_step;
            if (!std::isfinite(
                    step[column])) {
                throw std::range_error(
                    "SW92 transactional conservation anchor step is non-finite");
            }
        }

        bool accepted = false;
        double damping = 1.0;
        for (std::size_t backtrack = 0U;
             backtrack < maximum_backtracks;
             ++backtrack) {
            auto trial = *q;
            for (std::size_t column = 0U;
                 column < column_count;
                 ++column) {
                trial[column] +=
                    damping *
                    step[column];
            }
            if (!(trial[0] > 0.0) ||
                !(trial[1] > 0.0)) {
                damping *= 0.5;
                continue;
            }

            std::vector<double> trial_residual;
            std::vector<double> trial_jacobian;
            double trial_norm = 0.0;
            if (residual_and_jacobian(
                    trial,
                    &trial_residual,
                    &trial_jacobian,
                    &trial_norm) &&
                trial_norm < norm) {
                *q =
                    std::move(trial);
                accepted = true;
                break;
            }
            damping *= 0.5;
        }
        if (!accepted) {
            throw std::range_error(
                "SW92 transactional conservation anchor line search failed");
        }
    }

    std::vector<double> final_residual;
    std::vector<double> final_jacobian;
    double final_norm = 0.0;
    if (!residual_and_jacobian(
            *q,
            &final_residual,
            &final_jacobian,
            &final_norm) ||
        final_norm > 1.0e-6) {
        throw std::range_error(
            "SW92 transactional conservation anchor did not close component/energy storage");
    }
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
        const double scale =
            std::max(
                {1.0,
                 std::abs(first[i]),
                 std::abs(second[i])});
        if (!std::isfinite(first[i]) ||
            !std::isfinite(second[i]) ||
            std::abs(first[i] - second[i]) >
                8192.0 *
                    std::numeric_limits<double>::
                        epsilon() *
                    scale) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] inline PetscErrorCode
reject_absent_phase_provider(
    const mpmc::flow::
        AbsentPhaseThermodynamicCoordinateExtension&,
    void* user_context,
    std::optional<
        mpmc::flow::
            AbsentPhasePotentialExtensionLinearization>*,
    NaturalVariableSnesEvaluationStatus3D*) {
    return user_context == nullptr
        ? PETSC_ERR_ARG_NULL
        : PETSC_ERR_SUP;
}

[[nodiscard]] inline
std::vector<double>
owned_q(
    Vec state,
    const VariableCardinalityNaturalVariableCellDof3D&
        record) {
    std::vector<PetscInt>
        indices(record.scalar_count);
    std::vector<PetscScalar>
        values(record.scalar_count);
    for (std::size_t slot = 0U;
         slot < record.scalar_count;
         ++slot) {
        indices[slot] =
            record.petsc_global_scalar_start +
            static_cast<PetscInt>(slot);
    }
    const PetscErrorCode error =
        VecGetValues(
            state,
            static_cast<PetscInt>(
                indices.size()),
            indices.data(),
            values.data());
    if (error != PETSC_SUCCESS) {
        throw std::runtime_error(
            "SW92 transactional rebuild cannot read owned converged q");
    }
    std::vector<double>
        q(record.scalar_count);
    for (std::size_t slot = 0U;
         slot < record.scalar_count;
         ++slot) {
        q[slot] =
            static_cast<double>(
                PetscRealPart(values[slot]));
    }
    return q;
}

} // namespace sw92_transactional_restart_detail

template <typename Provider>
    requires std::copy_constructible<Provider>
struct Sw92TransactionalPhaseTransitionRuntime3D {
    using Closure =
        mpmc::flow::
            Sw92SelectedPhasePropertyClosure<
                double,
                Provider>;
    using OneContext =
        Sw92SinglePhaseProductionCellEvaluatorContext3D<
            Closure>;
    using TwoContext =
        Sw92TwoPhaseProductionCellEvaluatorContext3D<
            Closure>;
    using ThreeContext =
        Sw92ThreePhaseProductionCellEvaluatorContext3D<
            Closure>;

    std::vector<std::unique_ptr<Closure>>
        closures;
    std::vector<std::unique_ptr<OneContext>>
        one_phase_contexts;
    std::vector<std::unique_ptr<TwoContext>>
        two_phase_contexts;
    std::vector<std::unique_ptr<ThreeContext>>
        three_phase_contexts;
    std::unique_ptr<
        CellScopedMixedCardinalityEvaluatorDispatcher3D>
        dispatcher;
};

template <typename Provider>
    requires std::copy_constructible<Provider>
struct Sw92TransactionalPhaseTransitionRebuildContext3D {
    using Runtime =
        Sw92TransactionalPhaseTransitionRuntime3D<
            Provider>;

    MPI_Comm comm{MPI_COMM_NULL};
    const mpmc::discretization_petsc::
        ParallelOwnedConnectionSchedule3D*
        schedule{};
    const mpmc::mesh::PartitionSnapshot*
        partition{};
    const mpmc::discretization_petsc::
        PetscMpiAijSymbolicPreallocation3D*
        cell_bridge{};
    const mpmc::discretization_petsc::
        OwnedCellStructuralColumnPatternSnapshot3D*
        cell_pattern{};
    const mpmc::thermodynamics::Sw92Phase<double>*
        model{};
    Provider provider;
    mpmc::flow::SelectedPhasePropertyProvenance
        provenance;
    PostSnesSw92ProfileCPhaseTransitionScannerContext3D*
        scanner_context{};
    Sw92TransitionTargetIdentityBinding3D
        target_identity;
    MixedCardinalityPhysicalCellEvaluatorBindings3D
        initial_evaluators;
    double single_phase_relative_permeability{1.0};
    Sw92TwoPhaseRelativePermeabilityEvaluatorBinding3D
        two_phase_relative_permeability;
    Sw92ThreePhaseSaturationConstitutiveEvaluatorBinding3D
        three_phase_saturation;
    Sw92RockThermalStorageEvaluatorBinding3D
        rock_storage;
    MixedCardinalityPhysicalCellSourceEvaluatorBinding3D
        cell_source_evaluator;
    // Optional projection used only when the caller intentionally wants to
    // adjust the authoritative target q onto the frozen component/energy
    // storage manifold before the rebuilt SNES solve.  Source-complete
    // production property models should normally leave this false and let
    // SNES satisfy equilibrium and conservation simultaneously.
    bool project_target_to_conservation_storage{false};
    std::vector<
        Sw92TransactionalRebuildBaselineCell3D>
        baseline_cells;
    std::vector<std::unique_ptr<Runtime>>
        runtimes;
    int absent_provider_token{1};

    [[nodiscard]]
    MixedCardinalityPhysicalCellEvaluatorBindings3D
    active_evaluators() noexcept {
        return runtimes.empty()
            ? initial_evaluators
            : runtimes.back()
                  ->dispatcher
                  ->bindings();
    }
};

template <typename Provider>
    requires std::copy_constructible<Provider>
[[nodiscard]] inline PetscErrorCode
rebuild_sw92_transactional_phase_transition_system_3d(
    const PhaseTransitionRebuiltNaturalVariableSystem3D&
        current_system,
    Vec converged_state,
    const VariableCardinalityNaturalVariableSnesSolveReport3D&
        solve_report,
    std::span<
        const PostSnesPhaseTransitionProposal3D>
        local_owned_proposals,
    std::span<
        const AcceptedPhaseTransitionSummary3D>
        accepted_global_batch,
    void* raw_context,
    std::unique_ptr<
        PhaseTransitionRebuiltNaturalVariableSystem3D>*
            rebuilt_system) {
    using namespace
        sw92_transactional_restart_detail;
    using Context =
        Sw92TransactionalPhaseTransitionRebuildContext3D<
            Provider>;
    using Runtime =
        Sw92TransactionalPhaseTransitionRuntime3D<
            Provider>;
    using Closure =
        typename Runtime::Closure;

    if (raw_context == nullptr ||
        rebuilt_system == nullptr ||
        converged_state == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    rebuilt_system->reset();

    auto* context =
        static_cast<Context*>(
            raw_context);
    if (context->comm == MPI_COMM_NULL ||
        context->schedule == nullptr ||
        context->partition == nullptr ||
        context->cell_bridge == nullptr ||
        context->cell_pattern == nullptr ||
        context->model == nullptr ||
        context->scanner_context == nullptr ||
        context->target_identity.resolver == nullptr ||
        context->rock_storage.evaluator == nullptr ||
        context->two_phase_relative_permeability.evaluator ==
            nullptr ||
        context->three_phase_saturation.evaluator ==
            nullptr ||
        !context->schedule
             ->assembly_rows()
             .empty() ||
        static_cast<int>(
            solve_report.converged_reason) <=
            0 ||
        current_system.time_step_seconds() <=
            0.0 ||
        context->baseline_cells.size() !=
            current_system
                .numbering()
                .local_cell_count()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    const auto& numbering =
        current_system.numbering();
    for (const auto& record :
         numbering.cells()) {
        if (record.owner_rank !=
            numbering.local_rank()) {
            // SW92 has no family/root-aware absent-phase provider yet.  The
            // no-face transactional contract therefore refuses ghost overlap
            // rather than pretending cross-cardinality transport is available.
            return PETSC_ERR_SUP;
        }
    }

    std::vector<std::optional<
        MixedCardinalityPhysicalCurrentCellLinearization3D>>
        current;
    std::vector<double> porosities;
    NaturalVariableSnesEvaluationStatus3D
        evaluation_status =
            NaturalVariableSnesEvaluationStatus3D::
                success;
    PetscErrorCode error =
        current_system
            .evaluate_local_cells_for_phase_transition(
                converged_state,
                &current,
                &porosities,
                &evaluation_status);
    if (error != PETSC_SUCCESS) {
        return error;
    }
    if (evaluation_status !=
            NaturalVariableSnesEvaluationStatus3D::
                success ||
        current.size() !=
            context->baseline_cells.size() ||
        porosities.size() !=
            current.size()) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }

    const auto fallback =
        context->active_evaluators();
    auto runtime =
        std::make_unique<Runtime>();
    runtime->closures.reserve(
        context->baseline_cells.size());
    runtime->one_phase_contexts.reserve(
        context->baseline_cells.size());
    runtime->two_phase_contexts.reserve(
        context->baseline_cells.size());
    runtime->three_phase_contexts.reserve(
        context->baseline_cells.size());

    std::vector<
        CellScopedMixedCardinalityEvaluatorDispatcher3D::
            SingleEntry>
        one_entries;
    std::vector<
        CellScopedMixedCardinalityEvaluatorDispatcher3D::
            TwoEntry>
        two_entries;
    std::vector<
        CellScopedMixedCardinalityEvaluatorDispatcher3D::
            ThreeEntry>
        three_entries;
    std::vector<
        FrozenPhaseTransitionRebuildCell3D>
        cells;
    cells.reserve(
        context->baseline_cells.size());

    try {
        for (std::size_t local = 0U;
             local <
                 context->baseline_cells.size();
             ++local) {
            const auto& baseline =
                context->baseline_cells[local];
            const auto cell =
                mpmc::mesh::LocalIndex{
                    static_cast<
                        mpmc::mesh::LocalIndex::value_type>(
                            local)};
            const auto& record =
                numbering.cell(cell);
            if (baseline.cell != cell ||
                baseline.cell_global !=
                    record.cell_global ||
                baseline.component_ids.size() !=
                    numbering.component_count() ||
                !baseline
                     .previous_component_accumulation
                     .has_value() ||
                !baseline
                     .previous_energy_accumulation
                     .has_value() ||
                !current[local].has_value() ||
                !std::isfinite(
                    baseline.bulk_volume_m3) ||
                !(baseline.bulk_volume_m3 > 0.0) ||
                !std::isfinite(
                    baseline.porosity) ||
                !(baseline.porosity > 0.0) ||
                !(baseline.porosity < 1.0) ||
                !near_roundoff(
                    baseline.porosity,
                    porosities[local])) {
                throw std::invalid_argument(
                    "SW92 transactional rebuild baseline/current mismatch");
            }

            const auto& source_active =
                current_system
                    .coordinate_registry()
                    .cell(
                        record.cell_global)
                    .active_phases;
            const auto proposal =
                std::find_if(
                    local_owned_proposals.begin(),
                    local_owned_proposals.end(),
                    [&](const auto& item) {
                        return item.cell_global ==
                            record.cell_global;
                    });
            const auto* target =
                context->scanner_context
                    ->find_target(
                        record.cell_global);

            if ((proposal ==
                     local_owned_proposals.end()) !=
                (target == nullptr)) {
                throw std::invalid_argument(
                    "SW92 transactional proposal/sidecar mismatch");
            }

            if (proposal ==
                local_owned_proposals.end()) {
                const auto& identity =
                    state_identity(
                        *current[local]);
                auto q =
                    owned_q(
                        converged_state,
                        record);
                cells.push_back({
                    baseline.cell,
                    baseline.cell_global,
                    baseline.bulk_volume_m3,
                    baseline.porosity,
                    baseline.component_ids,
                    identity.layout,
                    std::move(q),
                    source_active,
                    {},
                    baseline
                        .previous_component_accumulation,
                    baseline
                        .previous_energy_accumulation,
                    "SW92/transactional/no-transition/v1"});

                if (record.phase_count == 1U) {
                    one_entries.push_back({
                        record.cell_global,
                        fallback.single_phase});
                } else if (
                    record.phase_count == 2U) {
                    two_entries.push_back({
                        record.cell_global,
                        fallback.two_phase});
                } else if (
                    record.phase_count == 3U) {
                    three_entries.push_back({
                        record.cell_global,
                        fallback.three_phase});
                } else {
                    throw std::invalid_argument(
                        "SW92 transactional current cardinality is not 1/2/3");
                }
                continue;
            }

            const auto accepted =
                std::find_if(
                    accepted_global_batch.begin(),
                    accepted_global_batch.end(),
                    [&](const auto& item) {
                        return item.cell_global ==
                            record.cell_global;
                    });
            if (accepted ==
                    accepted_global_batch.end() ||
                accepted->source_phase_count !=
                    record.phase_count ||
                accepted->target_phase_count !=
                    target->projection
                        .target_phase_count() ||
                target->projection
                        .source_phase_count() !=
                    record.phase_count ||
                proposal->candidate
                        .target_phase_count !=
                    target->phases.size() ||
                target->component_ids !=
                    baseline.component_ids) {
                throw std::invalid_argument(
                    "SW92 transactional accepted batch/sidecar mismatch");
            }

            std::optional<
                mpmc::flow::
                    FrozenActivePhaseIdentityMap>
                target_active;
            error =
                context->target_identity
                    .resolver(
                        record.cell_global,
                        *target,
                        source_active,
                        context->target_identity
                            .user_context,
                        &target_active);
            if (error != PETSC_SUCCESS) {
                return error;
            }
            if (!target_active.has_value() ||
                target_active->phase_count() !=
                    proposal->candidate
                        .target_phase_count) {
                throw std::invalid_argument(
                    "SW92 transactional target identity resolution is indeterminate");
            }

            const auto current_inventory =
                current_component_inventory(
                    *current[local],
                    baseline.porosity);
            auto rebuilt =
                make_accepted_phase_transition_rebuild_cell_3d(
                    baseline.cell,
                    baseline.cell_global,
                    baseline.bulk_volume_m3,
                    baseline.porosity,
                    proposal->candidate,
                    current_inventory,
                    *baseline
                         .previous_component_accumulation,
                    *baseline
                         .previous_energy_accumulation,
                    source_active,
                    *target_active,
                    {});

            if (rebuilt.target_layout.phase_count() !=
                    target->projection
                        .target_phase_count() ||
                !same_vector(
                    rebuilt
                        .target_natural_variables,
                    target->projection
                        .natural_variables())) {
                throw std::invalid_argument(
                    "SW92 transactional outer-rebuild projection differs from authoritative sidecar");
            }

            std::vector<
                mpmc::thermodynamics::
                    Sw92SelectedPhase<double>>
                selections;
            selections.reserve(
                target->phases.size());
            for (const auto& phase :
                 target->phases) {
                selections.push_back(
                    phase.selection);
            }
            auto closure =
                std::make_unique<Closure>(
                    *context->model,
                    std::move(selections),
                    context->provider,
                    context->provenance);
            Closure* closure_ptr =
                closure.get();
            runtime->closures.push_back(
                std::move(closure));

            switch (
                proposal->candidate
                    .target_phase_count) {
            case 1U: {
                auto evaluator =
                    std::make_unique<
                        typename Runtime::OneContext>();
                evaluator->property_closure =
                    closure_ptr;
                evaluator->relative_permeability =
                    context
                        ->single_phase_relative_permeability;
                evaluator->rock =
                    context->rock_storage;
                auto* evaluator_ptr =
                    evaluator.get();
                if (context
                        ->project_target_to_conservation_storage) {
                conservative_storage_anchor(
                    &rebuilt.target_natural_variables,
                    baseline.porosity,
                    current_inventory,
                    current_energy_inventory(
                        *current[local],
                        baseline.porosity),
                    [&](std::span<const double> q)
                        -> std::optional<
                            MixedCardinalityPhysicalCurrentCellLinearization3D> {
                        std::optional<
                            SinglePhaseCurrentCellLinearization3D>
                            value;
                        NaturalVariableSnesEvaluationStatus3D
                            status =
                                NaturalVariableSnesEvaluationStatus3D::
                                    success;
                        const auto typed_layout =
                            phase_transition_outer_rebuild_detail::
                                layout_1p(
                                    rebuilt.target_layout);
                        const PetscErrorCode evaluate_error =
                            evaluate_sw92_single_phase_production_cell_3d<
                                Closure>(
                                    baseline.cell,
                                    baseline.cell_global,
                                    q,
                                    typed_layout,
                                    baseline.component_ids,
                                    evaluator_ptr,
                                    &value,
                                    &status);
                        if (evaluate_error != PETSC_SUCCESS ||
                            status !=
                                NaturalVariableSnesEvaluationStatus3D::
                                    success ||
                            !value.has_value()) {
                            return std::nullopt;
                        }
                        return MixedCardinalityPhysicalCurrentCellLinearization3D{
                            std::move(*value)};
                    });
                }
                runtime->one_phase_contexts
                    .push_back(
                        std::move(evaluator));
                one_entries.push_back({
                    record.cell_global,
                    {
                        &evaluate_sw92_single_phase_production_cell_3d<
                            Closure>,
                        evaluator_ptr}});
                break;
            }
            case 2U: {
                auto evaluator =
                    std::make_unique<
                        typename Runtime::TwoContext>();
                evaluator->property_closure =
                    closure_ptr;
                evaluator->relative_permeability =
                    context
                        ->two_phase_relative_permeability;
                evaluator->rock =
                    context->rock_storage;
                auto* evaluator_ptr =
                    evaluator.get();
                if (context
                        ->project_target_to_conservation_storage) {
                conservative_storage_anchor(
                    &rebuilt.target_natural_variables,
                    baseline.porosity,
                    current_inventory,
                    current_energy_inventory(
                        *current[local],
                        baseline.porosity),
                    [&](std::span<const double> q)
                        -> std::optional<
                            MixedCardinalityPhysicalCurrentCellLinearization3D> {
                        std::optional<
                            TwoPhaseCurrentCellLinearization3D>
                            value;
                        NaturalVariableSnesEvaluationStatus3D
                            status =
                                NaturalVariableSnesEvaluationStatus3D::
                                    success;
                        const auto typed_layout =
                            phase_transition_outer_rebuild_detail::
                                layout_2p(
                                    rebuilt.target_layout);
                        const PetscErrorCode evaluate_error =
                            evaluate_sw92_two_phase_production_cell_3d<
                                Closure>(
                                    baseline.cell,
                                    baseline.cell_global,
                                    q,
                                    typed_layout,
                                    baseline.component_ids,
                                    evaluator_ptr,
                                    &value,
                                    &status);
                        if (evaluate_error != PETSC_SUCCESS ||
                            status !=
                                NaturalVariableSnesEvaluationStatus3D::
                                    success ||
                            !value.has_value()) {
                            return std::nullopt;
                        }
                        return MixedCardinalityPhysicalCurrentCellLinearization3D{
                            std::move(*value)};
                    });
                }
                runtime->two_phase_contexts
                    .push_back(
                        std::move(evaluator));
                two_entries.push_back({
                    record.cell_global,
                    {
                        &evaluate_sw92_two_phase_production_cell_3d<
                            Closure>,
                        evaluator_ptr}});
                break;
            }
            case 3U: {
                auto evaluator =
                    std::make_unique<
                        typename Runtime::ThreeContext>();
                evaluator->property_closure =
                    closure_ptr;
                evaluator->saturation_constitutive =
                    context
                        ->three_phase_saturation;
                evaluator->rock =
                    context->rock_storage;
                auto* evaluator_ptr =
                    evaluator.get();
                if (context
                        ->project_target_to_conservation_storage) {
                conservative_storage_anchor(
                    &rebuilt.target_natural_variables,
                    baseline.porosity,
                    current_inventory,
                    current_energy_inventory(
                        *current[local],
                        baseline.porosity),
                    [&](std::span<const double> q)
                        -> std::optional<
                            MixedCardinalityPhysicalCurrentCellLinearization3D> {
                        std::optional<
                            FixedThreePhaseCurrentCellLinearization3D>
                            value;
                        NaturalVariableSnesEvaluationStatus3D
                            status =
                                NaturalVariableSnesEvaluationStatus3D::
                                    success;
                        const auto typed_layout =
                            static_cast<
                                mpmc::flow::NaturalVariableLayout3P>(
                                    rebuilt.target_layout);
                        const PetscErrorCode evaluate_error =
                            evaluate_sw92_three_phase_production_cell_3d<
                                Closure>(
                                    baseline.cell,
                                    baseline.cell_global,
                                    q,
                                    typed_layout,
                                    baseline.component_ids,
                                    evaluator_ptr,
                                    &value,
                                    &status);
                        if (evaluate_error != PETSC_SUCCESS ||
                            status !=
                                NaturalVariableSnesEvaluationStatus3D::
                                    success ||
                            !value.has_value()) {
                            return std::nullopt;
                        }
                        return MixedCardinalityPhysicalCurrentCellLinearization3D{
                            std::move(*value)};
                    });
                }
                runtime->three_phase_contexts
                    .push_back(
                        std::move(evaluator));
                three_entries.push_back({
                    record.cell_global,
                    {
                        &evaluate_sw92_three_phase_production_cell_3d<
                            Closure>,
                        evaluator_ptr}});
                break;
            }
            default:
                throw std::invalid_argument(
                    "SW92 transactional target cardinality is not 1/2/3");
            }

            cells.push_back(
                std::move(rebuilt));
        }

        for (const auto& accepted :
             accepted_global_batch) {
            if (context->partition
                    ->contains_global(
                        mpmc::mesh::EntityKind::cell,
                        accepted.cell_global)) {
                const auto local =
                    context->partition
                        ->local_index(
                            mpmc::mesh::EntityKind::cell,
                            accepted.cell_global);
                if (!context->partition
                         ->is_owned(
                             mpmc::mesh::EntityKind::cell,
                             local)) {
                    return PETSC_ERR_SUP;
                }
            }
        }
    } catch (const std::exception&) {
        return PETSC_ERR_ARG_INCOMP;
    }

    try {
        runtime->dispatcher =
            std::make_unique<
                CellScopedMixedCardinalityEvaluatorDispatcher3D>(
                    std::move(one_entries),
                    std::move(two_entries),
                    std::move(three_entries));
    } catch (const std::exception&) {
        return PETSC_ERR_ARG_INCOMP;
    }

    std::unique_ptr<
        PhaseTransitionRebuiltNaturalVariableSystem3D>
        next;
    error =
        rebuild_phase_transition_natural_variable_system_3d(
            context->comm,
            *context->schedule,
            *context->partition,
            *context->cell_bridge,
            *context->cell_pattern,
            current_system
                .time_step_seconds(),
            std::move(cells),
            {},
            runtime->dispatcher
                ->bindings(),
            context->cell_source_evaluator,
            &reject_absent_phase_provider,
            &context->absent_provider_token,
            &next);
    if (error != PETSC_SUCCESS ||
        next == nullptr) {
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_PLIB;
    }

    context->runtimes.push_back(
        std::move(runtime));
    *rebuilt_system =
        std::move(next);
    return PETSC_SUCCESS;
}

} // namespace mpmc::flow_discretization_petsc

#endif // MPMC_FLOW_DISCRETIZATION_PETSC_SW92_TRANSACTIONAL_PHASE_TRANSITION_RESTART_HPP
