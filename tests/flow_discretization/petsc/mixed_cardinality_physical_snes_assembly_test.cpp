#include <mpmc/flow_discretization_petsc/physical_timestep_driver.hpp>
#include <mpmc/flow_discretization_petsc/post_snes_phase_transition_controller.hpp>
#include <mpmc/flow_discretization_petsc/post_snes_phase_transition_handoff_scanner.hpp>
#include <mpmc/well_discretization_petsc/fixed_bhp_well_source_evaluator.hpp>
#include <mpmc/well_discretization_petsc/fixed_total_molar_rate_control.hpp>
#include <mpmc/well_discretization_petsc/fixed_total_molar_rate_timestep_driver.hpp>
#include <mpmc/thermodynamics/pr_parameters.hpp>

#include <petscmat.h>
#include <petscvec.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace {

namespace mesh = mpmc::mesh;
namespace flow = mpmc::flow;
namespace disc = mpmc::discretization;
namespace dp = mpmc::discretization_petsc;
namespace fd = mpmc::flow_discretization;
namespace fdp = mpmc::flow_discretization_petsc;
namespace wd = mpmc::well_discretization;
namespace wdp = mpmc::well_discretization_petsc;
namespace well = mpmc::well;
namespace th = mpmc::thermodynamics;

void require_collective(
    bool condition,
    std::string_view message) {
    int local = condition ? 1 : 0;
    int global = 0;
    if (MPI_Allreduce(
            &local,
            &global,
            1,
            MPI_INT,
            MPI_MIN,
            PETSC_COMM_WORLD) !=
        MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI_Allreduce failed in mixed physical dispatcher regression");
    }
    if (global == 0) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_collective(
    double actual,
    double expected,
    double relative,
    double absolute) {
    require_collective(
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                absolute +
                    relative *
                        std::max(
                            std::abs(actual),
                            std::abs(expected)),
        "numeric mismatch in mixed physical dispatcher regression");
}

[[nodiscard]] std::size_t
phase_count_for(
    std::uint64_t stable) {
    if (stable == UINT64_C(10) ||
        stable == UINT64_C(20)) {
        return 1U;
    }
    if (stable == UINT64_C(30) ||
        stable == UINT64_C(40)) {
        return 2U;
    }
    if (stable == UINT64_C(50) ||
        stable == UINT64_C(60)) {
        return 3U;
    }
    throw std::invalid_argument(
        "unknown mixed physical stable cell");
}

[[nodiscard]] flow::NaturalVariableLayout1P
layout_1p() {
    return flow::NaturalVariableLayout1P{
        flow::NaturalVariableCompositionPivot1P::
            fixed_last(3U)};
}

[[nodiscard]] flow::NaturalVariableLayout2P
layout_2p() {
    return flow::NaturalVariableLayout2P{
        flow::NaturalVariableCompositionPivot2P::
            fixed_last(3U)};
}

[[nodiscard]] flow::NaturalVariableLayout3P
layout_3p() {
    return flow::NaturalVariableLayout3P{
        flow::NaturalVariableCompositionPivot3P::
            fixed_last(3U)};
}

[[nodiscard]] std::vector<double>
target_1p() {
    auto layout = layout_1p();
    std::vector<double> q(
        layout.unknown_count(),
        0.0);
    q[layout.pressure_unknown_index()] =
        10.0;
    q[layout.temperature_unknown_index()] =
        8.0;
    q[*layout.independent_composition_unknown_index(
        0U)] = 0.20;
    q[*layout.independent_composition_unknown_index(
        1U)] = 0.30;
    return q;
}

[[nodiscard]] std::vector<double>
target_2p() {
    auto layout = layout_2p();
    std::vector<double> q(
        layout.unknown_count(),
        0.0);
    q[layout.pressure_unknown_index()] =
        20.0;
    q[layout.temperature_unknown_index()] =
        9.0;
    q[layout.independent_saturation_unknown_index()] =
        0.55;
    const std::array<
        std::array<double, 3>,
        2>
        x{{
            {{0.20, 0.30, 0.50}},
            {{0.40, 0.20, 0.40}}
        }};
    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        phase,
                        component);
            if (column) {
                q[*column] =
                    x[phase][component];
            }
        }
    }
    return q;
}

[[nodiscard]] std::vector<double>
target_3p() {
    auto layout = layout_3p();
    std::vector<double> q(
        layout.unknown_count(),
        0.0);
    q[layout.pressure_unknown_index()] =
        30.0;
    q[layout.temperature_unknown_index()] =
        10.0;
    const auto s0 =
        layout.independent_saturation_unknown_index(
            flow::PhaseSlot3::phase0);
    const auto s1 =
        layout.independent_saturation_unknown_index(
            flow::PhaseSlot3::phase1);
    if (!s0 || !s1) {
        throw std::logic_error(
            "three-phase target layout lacks saturation coordinates");
    }
    q[*s0] = 0.20;
    q[*s1] = 0.30;

    const std::array<
        std::array<double, 3>,
        3>
        x{{
            {{0.20, 0.30, 0.50}},
            {{0.40, 0.20, 0.40}},
            {{0.10, 0.50, 0.40}}
        }};
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<flow::PhaseSlot3>(
                phase);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        slot,
                        component);
            if (column) {
                q[*column] =
                    x[phase][component];
            }
        }
    }
    return q;
}

[[nodiscard]] std::vector<double>
target_state(
    std::uint64_t stable) {
    switch (phase_count_for(stable)) {
    case 1U:
        return target_1p();
    case 2U:
        return target_2p();
    case 3U:
        return target_3p();
    default:
        throw std::logic_error(
            "invalid target phase count");
    }
}

struct DispatchAudit {
    std::uint64_t single_calls{};
    std::uint64_t two_calls{};
    std::uint64_t three_calls{};
    bool source_enabled{};
    bool well_timestep_compressibility{};
    std::uint64_t source_evaluator_calls{};
    std::uint64_t source_cell20_calls{};
};

PetscErrorCode evaluate_explicit_cell_source(
    mesh::LocalIndex,
    mesh::GlobalEntityId cell_global,
    std::span<const double> natural_variables,
    const fdp::MixedCardinalityPhysicalCurrentCellLinearization3D&
        current,
    void* raw_context,
    std::optional<fd::CellSourceLinearization3D>* output,
    fdp::NaturalVariableSnesEvaluationStatus3D* status) {
    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    auto* audit =
        static_cast<DispatchAudit*>(raw_context);
    ++audit->source_evaluator_calls;
    *status =
        fdp::NaturalVariableSnesEvaluationStatus3D::
            success;
    if (!audit->source_enabled ||
        cell_global !=
            mesh::GlobalEntityId{UINT64_C(20)}) {
        return PETSC_SUCCESS;
    }
    ++audit->source_cell20_calls;

    const auto identity =
        std::visit(
            [](const auto& typed) {
                return typed.transport.state_identity;
            },
            current);
    const std::size_t q =
        natural_variables.size();
    if (identity.component_ids.size() != 3U ||
        q != identity.layout.unknown_count()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    fd::CellSourceLinearization3D source;
    source.provenance =
        "test/mixed-cell-source/v1";
    source.component_ids =
        identity.component_ids;
    source.input_count = q;
    source.component_molar_rate_mol_per_s =
        {2.0, -1.0, 0.5};
    source.component_molar_rate_jacobian_mol_per_s.assign(
        3U * q,
        0.0);
    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        for (std::size_t column = 0U;
             column < q;
             ++column) {
            source.component_molar_rate_jacobian_mol_per_s[
                component * q + column] =
                0.1 *
                static_cast<double>(
                    (component + 1U) *
                    (column + 1U));
        }
    }
    source.energy_rate_w = 20.0;
    source.energy_rate_gradient_w.resize(q);
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        source.energy_rate_gradient_w[column] =
            0.2 *
            static_cast<double>(
                column + 1U);
    }
    output->emplace(std::move(source));
    return PETSC_SUCCESS;
}

struct FixedBhpWellSourceAudit {
    wdp::
        FixedBhpPeacemanWellSourceEvaluatorContext3D*
            context{};
    std::uint64_t evaluator_calls{};
    std::uint64_t target_calls{};
};

PetscErrorCode
evaluate_audited_fixed_bhp_well_source(
    mesh::LocalIndex cell,
    mesh::GlobalEntityId cell_global,
    std::span<const double> natural_variables,
    const fdp::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            current,
    void* raw_context,
    std::optional<fd::CellSourceLinearization3D>* output,
    fdp::NaturalVariableSnesEvaluationStatus3D* status) {
    if (raw_context == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* audit =
        static_cast<FixedBhpWellSourceAudit*>(
            raw_context);
    if (audit->context == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    ++audit->evaluator_calls;
    if (cell_global ==
        audit->context->target_cell_global()) {
        ++audit->target_calls;
    }
    return wdp::
        evaluate_fixed_bhp_peaceman_well_source_3d(
            cell,
            cell_global,
            natural_variables,
            current,
            audit->context,
            output,
            status);
}

struct FixedBhpMultiConnectionWellSourceAudit {
    wdp::
        FixedBhpMultiConnectionWellSourceEvaluatorContext3D*
            context{};
    std::uint64_t evaluator_calls{};
    std::uint64_t cell30_calls{};
    std::uint64_t cell60_calls{};
};

PetscErrorCode
evaluate_audited_fixed_bhp_multi_connection_well_source(
    mesh::LocalIndex cell,
    mesh::GlobalEntityId cell_global,
    std::span<const double> natural_variables,
    const fdp::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            current,
    void* raw_context,
    std::optional<fd::CellSourceLinearization3D>* output,
    fdp::NaturalVariableSnesEvaluationStatus3D* status) {
    if (raw_context == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* audit =
        static_cast<
            FixedBhpMultiConnectionWellSourceAudit*>(
                raw_context);
    if (audit->context == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    ++audit->evaluator_calls;
    if (audit->context
            ->find_connection(cell_global) !=
        nullptr) {
        if (cell_global ==
            mesh::GlobalEntityId{
                UINT64_C(30)}) {
            ++audit->cell30_calls;
        } else if (
            cell_global ==
            mesh::GlobalEntityId{
                UINT64_C(60)}) {
            ++audit->cell60_calls;
        }
    }

    return wdp::
        evaluate_fixed_bhp_multi_connection_well_source_3d(
            cell,
            cell_global,
            natural_variables,
            current,
            audit->context,
            output,
            status);
}

PetscErrorCode
evaluate_1p(
    mesh::LocalIndex,
    mesh::GlobalEntityId,
    std::span<const double> natural_variables,
    const flow::NaturalVariableLayout1P&
        frozen_layout,
    std::span<const std::string> component_ids,
    void* raw_context,
    std::optional<
        fdp::SinglePhaseCurrentCellLinearization3D>*
        output,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    auto* audit =
        static_cast<DispatchAudit*>(
            raw_context);
    ++audit->single_calls;

    try {
        if (natural_variables.size() !=
                frozen_layout.unknown_count() ||
            component_ids.size() !=
                frozen_layout.component_count()) {
            return PETSC_ERR_ARG_INCOMP;
        }
        const double p =
            natural_variables[
                frozen_layout
                    .pressure_unknown_index()];
        const double temperature =
            natural_variables[
                frozen_layout
                    .temperature_unknown_index()];
        if (!std::isfinite(p) ||
            !(p > 0.0) ||
            !std::isfinite(temperature) ||
            !(temperature > 0.0)) {
            *status =
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
            return PETSC_SUCCESS;
        }

        flow::NaturalVariableCellStateInput1P
            input;
        input.component_ids.assign(
            component_ids.begin(),
            component_ids.end());
        input.reference_pressure_pa = p;
        input.temperature_k = temperature;
        input.composition_pivot =
            frozen_layout.composition_pivot();
        for (std::size_t rank = 0U;
             rank < 2U;
             ++rank) {
            const auto component =
                frozen_layout
                    .independent_composition_component(
                        rank);
            const auto column =
                frozen_layout
                    .independent_composition_unknown_index(
                        component);
            if (!column) {
                return PETSC_ERR_PLIB;
            }
            input.independent_composition
                .push_back(
                    natural_variables[*column]);
        }
        const bool synthetic_compressibility =
            audit
                ->well_timestep_compressibility;
        const double pressure_delta =
            p - 10.0;
        const double molar_density =
            5.0 +
            (synthetic_compressibility
                 ? 0.05 * pressure_delta
                 : 0.0);
        const double mass_density =
            2.0 +
            (synthetic_compressibility
                 ? 0.02 * pressure_delta
                 : 0.0);
        if (!(molar_density > 0.0) ||
            !(mass_density > 0.0)) {
            *status =
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
            return PETSC_SUCCESS;
        }
        input.phase_properties =
            flow::PhasePropertyPrerequisiteInput{
                molar_density,
                mass_density,
                1.0,
                temperature,
                0.5 * temperature};

        auto state =
            flow::NaturalVariableCellState1P::
                create(std::move(input));
        const std::size_t q =
            frozen_layout.unknown_count();
        std::vector<double> zero(
            q,
            0.0);
        std::vector<double> molar_density_gradient(
            q,
            0.0);
        std::vector<double> mass_density_gradient(
            q,
            0.0);
        if (synthetic_compressibility) {
            molar_density_gradient[
                frozen_layout
                    .pressure_unknown_index()] =
                0.05;
            mass_density_gradient[
                frozen_layout
                    .pressure_unknown_index()] =
                0.02;
        }
        auto molar =
            flow::
                make_single_phase_molar_density_linearization(
                    state,
                    molar_density_gradient);
        const flow::TransportPropertyProvenance
            provenance{
                "mixed-physical-dispatch",
                "controlled-regression",
                "v1"};
        auto transport =
            flow::
                make_single_phase_transport_linearization(
                    state,
                    mass_density_gradient,
                    zero,
                    1.0,
                    zero,
                    provenance,
                    provenance);
        std::vector<double> dh(
            q,
            0.0);
        std::vector<double> du(
            q,
            0.0);
        dh[frozen_layout
               .temperature_unknown_index()] =
            1.0;
        du[frozen_layout
               .temperature_unknown_index()] =
            0.5;
        auto caloric =
            flow::
                make_single_phase_caloric_linearization(
                    state,
                    dh,
                    du,
                    provenance,
                    provenance);
        std::vector<double> rock_gradient(
            q,
            0.0);
        rock_gradient[
            frozen_layout
                .temperature_unknown_index()] =
            2.0;
        auto rock =
            flow::
                make_single_phase_rock_thermal_storage_linearization(
                    state,
                    2.0 * temperature,
                    rock_gradient,
                    provenance);

        output->emplace(
            fdp::
                SinglePhaseCurrentCellLinearization3D{
                    std::move(state),
                    std::move(molar),
                    std::move(transport),
                    std::move(caloric),
                    std::move(rock)});
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                domain_error;
        return PETSC_SUCCESS;
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }
}

PetscErrorCode
evaluate_2p(
    mesh::LocalIndex,
    mesh::GlobalEntityId,
    std::span<const double> natural_variables,
    const flow::NaturalVariableLayout2P&
        frozen_layout,
    std::span<const std::string> component_ids,
    void* raw_context,
    std::optional<
        fdp::TwoPhaseCurrentCellLinearization3D>*
        output,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    auto* audit =
        static_cast<DispatchAudit*>(
            raw_context);
    ++audit->two_calls;

    try {
        if (natural_variables.size() !=
                frozen_layout.unknown_count() ||
            component_ids.size() !=
                frozen_layout.component_count()) {
            return PETSC_ERR_ARG_INCOMP;
        }

        const double p =
            natural_variables[
                frozen_layout
                    .pressure_unknown_index()];
        const double temperature =
            natural_variables[
                frozen_layout
                    .temperature_unknown_index()];
        const double s0 =
            natural_variables[
                frozen_layout
                    .independent_saturation_unknown_index()];
        if (!std::isfinite(p) ||
            !(p > 0.0) ||
            !std::isfinite(temperature) ||
            !(temperature > 0.0) ||
            !std::isfinite(s0) ||
            !(s0 > 0.0) ||
            !(s0 < 1.0)) {
            *status =
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
            return PETSC_SUCCESS;
        }

        flow::NaturalVariableCellStateInput2P
            input;
        input.component_ids.assign(
            component_ids.begin(),
            component_ids.end());
        input.reference_pressure_pa = p;
        input.temperature_k = temperature;
        input.independent_saturation = s0;
        input.composition_pivot =
            frozen_layout.composition_pivot();

        for (std::size_t phase = 0U;
             phase < 2U;
             ++phase) {
            for (std::size_t rank = 0U;
                 rank < 2U;
                 ++rank) {
                const auto component =
                    frozen_layout
                        .independent_composition_component(
                            phase,
                            rank);
                const auto column =
                    frozen_layout
                        .independent_composition_unknown_index(
                            phase,
                            component);
                if (!column) {
                    return PETSC_ERR_PLIB;
                }
                input.independent_phase_compositions[
                    phase]
                    .push_back(
                        natural_variables[
                            *column]);
            }
        }

        const bool synthetic_compressibility =
            audit->well_timestep_compressibility;
        const double pressure_delta =
            p - 20.0;
        const std::array<double, 2>
            molar_density{
                6.0 +
                    (synthetic_compressibility
                         ? 0.06 * pressure_delta
                         : 0.0),
                3.0 +
                    (synthetic_compressibility
                         ? 0.03 * pressure_delta
                         : 0.0)};
        const std::array<double, 2>
            mass_density{
                2.0 +
                    (synthetic_compressibility
                         ? 0.02 * pressure_delta
                         : 0.0),
                1.5 +
                    (synthetic_compressibility
                         ? 0.015 * pressure_delta
                         : 0.0)};
        for (std::size_t phase = 0U;
             phase < 2U;
             ++phase) {
            if (!(molar_density[phase] > 0.0) ||
                !(mass_density[phase] > 0.0)) {
                *status =
                    fdp::NaturalVariableSnesEvaluationStatus3D::
                        domain_error;
                return PETSC_SUCCESS;
            }
        }

        input.phase_properties[0] =
            flow::PhasePropertyPrerequisiteInput{
                molar_density[0],
                mass_density[0],
                1.0,
                temperature,
                0.5 * temperature};
        input.phase_properties[1] =
            flow::PhasePropertyPrerequisiteInput{
                molar_density[1],
                mass_density[1],
                1.0,
                1.2 * temperature,
                0.8 * temperature};

        auto state =
            flow::NaturalVariableCellState2P::
                create(std::move(input));
        const std::size_t q =
            frozen_layout.unknown_count();
        std::array<std::vector<double>, 2>
            zero{
                std::vector<double>(q, 0.0),
                std::vector<double>(q, 0.0)};
        std::array<std::vector<double>, 2>
            molar_density_gradient = zero;
        std::array<std::vector<double>, 2>
            mass_density_gradient = zero;
        if (synthetic_compressibility) {
            const auto pressure_column =
                frozen_layout
                    .pressure_unknown_index();
            molar_density_gradient[0][
                pressure_column] =
                0.06;
            molar_density_gradient[1][
                pressure_column] =
                0.03;
            mass_density_gradient[0][
                pressure_column] =
                0.02;
            mass_density_gradient[1][
                pressure_column] =
                0.015;
        }
        auto molar =
            flow::
                make_two_phase_molar_density_linearization(
                    state,
                    molar_density_gradient);
        const flow::TransportPropertyProvenance
            provenance{
                "mixed-physical-dispatch",
                "controlled-regression",
                "v1"};
        auto transport =
            flow::
                make_two_phase_transport_linearization(
                    state,
                    mass_density_gradient,
                    zero,
                    {0.60, 0.40},
                    zero,
                    provenance,
                    provenance);
        std::array<std::vector<double>, 2>
            dh{
                std::vector<double>(q, 0.0),
                std::vector<double>(q, 0.0)};
        std::array<std::vector<double>, 2>
            du{
                std::vector<double>(q, 0.0),
                std::vector<double>(q, 0.0)};
        dh[0][frozen_layout
                  .temperature_unknown_index()] =
            1.0;
        dh[1][frozen_layout
                  .temperature_unknown_index()] =
            1.2;
        du[0][frozen_layout
                  .temperature_unknown_index()] =
            0.5;
        du[1][frozen_layout
                  .temperature_unknown_index()] =
            0.8;
        auto caloric =
            flow::
                make_two_phase_caloric_linearization(
                    state,
                    dh,
                    du,
                    provenance,
                    provenance);
        std::vector<double> rock_gradient(
            q,
            0.0);
        rock_gradient[
            frozen_layout
                .temperature_unknown_index()] =
            2.0;
        auto rock =
            flow::
                make_two_phase_rock_thermal_storage_linearization(
                    state,
                    2.0 * temperature,
                    rock_gradient,
                    provenance);

        flow::TwoPhaseFugacityEquilibriumLinearization
            fugacity;
        fugacity.layout =
            frozen_layout.descriptor();
        fugacity.component_ids.assign(
            component_ids.begin(),
            component_ids.end());
        fugacity.residual.assign(
            3U,
            0.0);
        fugacity.input_count = q;
        fugacity.jacobian.assign(
            3U * q,
            0.0);
        for (std::size_t row = 0U;
             row < 3U;
             ++row) {
            fugacity.jacobian[
                row * q +
                (3U + row)] =
                1.0;
        }

        output->emplace(
            fdp::
                TwoPhaseCurrentCellLinearization3D{
                    std::move(state),
                    std::move(molar),
                    std::move(transport),
                    std::move(caloric),
                    std::move(rock),
                    std::move(fugacity)});
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                domain_error;
        return PETSC_SUCCESS;
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }
}

struct ConstantRelativePermeability3P {
    template <typename Number>
    flow::RelativePermeabilityEvaluation3P<Number>
    operator()(
        const flow::ThreePhaseSaturationState3P<Number>&)
        const {
        return {
            std::array<Number, 3>{
                Number{0.50},
                Number{0.30},
                Number{0.20}}};
    }
};

PetscErrorCode
evaluate_3p(
    mesh::LocalIndex,
    mesh::GlobalEntityId,
    std::span<const double> natural_variables,
    const flow::NaturalVariableLayout3P&
        frozen_layout,
    std::span<const std::string> component_ids,
    void* raw_context,
    std::optional<
        fdp::FixedThreePhaseCurrentCellLinearization3D>*
        output,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    auto* audit =
        static_cast<DispatchAudit*>(
            raw_context);
    ++audit->three_calls;

    try {
        if (natural_variables.size() !=
                frozen_layout.unknown_count() ||
            component_ids.size() !=
                frozen_layout.component_count()) {
            return PETSC_ERR_ARG_INCOMP;
        }

        const std::size_t q =
            frozen_layout.unknown_count();
        const double p =
            natural_variables[
                frozen_layout
                    .pressure_unknown_index()];
        const double temperature =
            natural_variables[
                frozen_layout
                    .temperature_unknown_index()];
        const auto s0 =
            frozen_layout
                .independent_saturation_unknown_index(
                    flow::PhaseSlot3::phase0);
        const auto s1 =
            frozen_layout
                .independent_saturation_unknown_index(
                    flow::PhaseSlot3::phase1);
        if (!s0 || !s1 ||
            !std::isfinite(p) ||
            !(p > 0.0) ||
            !std::isfinite(temperature) ||
            !(temperature > 0.0)) {
            *status =
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
            return PETSC_SUCCESS;
        }

        flow::NaturalVariableCellStateInput3P
            input;
        input.component_ids.assign(
            component_ids.begin(),
            component_ids.end());
        input.reference_pressure_pa = p;
        input.temperature_k = temperature;
        input.independent_saturations = {
            natural_variables[*s0],
            natural_variables[*s1]};
        input.composition_pivot =
            frozen_layout.composition_pivot();

        for (std::size_t phase = 0U;
             phase < 3U;
             ++phase) {
            const auto slot =
                static_cast<flow::PhaseSlot3>(
                    phase);
            for (std::size_t rank = 0U;
                 rank < 2U;
                 ++rank) {
                const auto component =
                    frozen_layout
                        .independent_composition_component(
                            slot,
                            rank);
                const auto column =
                    frozen_layout
                        .independent_composition_unknown_index(
                            slot,
                            component);
                if (!column) {
                    return PETSC_ERR_PLIB;
                }
                input.independent_phase_compositions[
                    phase]
                    .push_back(
                        natural_variables[
                            *column]);
            }
        }

        const bool synthetic_compressibility =
            audit->well_timestep_compressibility;
        const double pressure_delta =
            p - 30.0;
        const std::array<double, 3>
            molar_density{
                2.0 +
                    (synthetic_compressibility
                         ? 0.02 * pressure_delta
                         : 0.0),
                4.0 +
                    (synthetic_compressibility
                         ? 0.04 * pressure_delta
                         : 0.0),
                7.0 +
                    (synthetic_compressibility
                         ? 0.07 * pressure_delta
                         : 0.0)};
        const std::array<double, 3>
            mass_density{
                1.0 +
                    (synthetic_compressibility
                         ? 0.01 * pressure_delta
                         : 0.0),
                2.0 +
                    (synthetic_compressibility
                         ? 0.02 * pressure_delta
                         : 0.0),
                3.0 +
                    (synthetic_compressibility
                         ? 0.03 * pressure_delta
                         : 0.0)};
        const std::array<double, 3>
            viscosity{
                1.0, 1.2, 1.4};
        for (std::size_t phase = 0U;
             phase < 3U;
             ++phase) {
            if (!(molar_density[phase] > 0.0) ||
                !(mass_density[phase] > 0.0)) {
                *status =
                    fdp::NaturalVariableSnesEvaluationStatus3D::
                        domain_error;
                return PETSC_SUCCESS;
            }
        }
        for (std::size_t phase = 0U;
             phase < 3U;
             ++phase) {
            const double phase_scale =
                static_cast<double>(
                    phase + 1U);
            input.phase_properties[phase] =
                flow::PhasePropertyPrerequisiteInput{
                    molar_density[phase],
                    mass_density[phase],
                    viscosity[phase],
                    phase_scale *
                        temperature,
                    0.5 *
                        phase_scale *
                        temperature};
        }

        auto state =
            flow::NaturalVariableCellState3P::
                create(std::move(input));

        std::array<std::vector<double>, 3>
            zero;
        std::array<std::vector<double>, 3>
            molar_density_gradient;
        std::array<std::vector<double>, 3>
            mass_density_gradient;
        std::array<std::vector<double>, 3>
            dh;
        std::array<std::vector<double>, 3>
            du;
        for (std::size_t phase = 0U;
             phase < 3U;
             ++phase) {
            zero[phase].assign(
                q,
                0.0);
            molar_density_gradient[phase].assign(
                q,
                0.0);
            mass_density_gradient[phase].assign(
                q,
                0.0);
            if (synthetic_compressibility) {
                const auto pressure_column =
                    frozen_layout
                        .pressure_unknown_index();
                const double phase_scale =
                    static_cast<double>(
                        phase + 1U);
                molar_density_gradient[phase][
                    pressure_column] =
                    0.02 * phase_scale;
                mass_density_gradient[phase][
                    pressure_column] =
                    0.01 * phase_scale;
            }
            dh[phase].assign(
                q,
                0.0);
            du[phase].assign(
                q,
                0.0);
            const double phase_scale =
                static_cast<double>(
                    phase + 1U);
            dh[phase][
                frozen_layout
                    .temperature_unknown_index()] =
                phase_scale;
            du[phase][
                frozen_layout
                    .temperature_unknown_index()] =
                0.5 * phase_scale;
        }

        flow::
            PhaseMolarDensityNaturalVariableLinearization3P
            molar{
                frozen_layout,
                molar_density,
                molar_density_gradient};
        const flow::TransportPropertyProvenance
            provenance{
                "mixed-physical-dispatch",
                "controlled-regression",
                "v1"};
        auto transport =
            flow::
                make_phase_transport_property_linearization(
                    state,
                    mass_density_gradient,
                    zero,
                    provenance,
                    provenance);
        auto caloric =
            flow::
                make_phase_caloric_property_linearization(
                    state,
                    dh,
                    du,
                    provenance,
                    provenance);
        std::vector<double> rock_gradient(
            q,
            0.0);
        rock_gradient[
            frozen_layout
                .temperature_unknown_index()] =
            2.0;
        auto rock =
            flow::
                make_stationary_rock_thermal_storage_linearization(
                    state,
                    2.0 * temperature,
                    rock_gradient,
                    provenance);

        const auto saturation =
            flow::
                evaluate_three_phase_saturation_constitutive(
                    p,
                    natural_variables[*s0],
                    natural_variables[*s1],
                    ConstantRelativePermeability3P{},
                    flow::NoCapillaryPressure3P{});
        flow::
            ThreePhaseSaturationCoordinateDerivatives3P
            saturation_derivatives{};
        auto saturation_linearization =
            flow::
                make_saturation_constitutive_natural_variable_linearization(
                    state,
                    saturation,
                    saturation_derivatives);

        std::vector<double> fugacity_values(
            6U,
            0.0);
        std::vector<double> fugacity_jacobian(
            6U * q,
            0.0);
        for (std::size_t row = 0U;
             row < 6U;
             ++row) {
            fugacity_jacobian[
                row * q +
                (4U + row)] =
                1.0;
        }
        flow::
            FugacityEquilibriumResidualLinearization3P
            fugacity{
                frozen_layout,
                std::vector<std::string>{
                    component_ids.begin(),
                    component_ids.end()},
                flow::
                    FugacityEquilibriumResidual3P<double>{
                        3U,
                        std::move(
                            fugacity_values)},
                q,
                std::move(
                    fugacity_jacobian)};

        output->emplace(
            fdp::
                FixedThreePhaseCurrentCellLinearization3D{
                    std::move(state),
                    std::move(molar),
                    std::move(transport),
                    std::move(caloric),
                    std::move(rock),
                    std::move(
                        saturation_linearization),
                    std::move(fugacity)});
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        *status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                domain_error;
        return PETSC_SUCCESS;
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }
}

[[nodiscard]] mesh::Topology
make_topology() {
    mesh::Topology::EntityIds ids;
    ids.faces = {
        mesh::GlobalEntityId{UINT64_C(100)},
        mesh::GlobalEntityId{UINT64_C(200)},
        mesh::GlobalEntityId{UINT64_C(300)},
        mesh::GlobalEntityId{UINT64_C(400)}};
    ids.cells = {
        mesh::GlobalEntityId{UINT64_C(10)},
        mesh::GlobalEntityId{UINT64_C(20)},
        mesh::GlobalEntityId{UINT64_C(30)},
        mesh::GlobalEntityId{UINT64_C(40)},
        mesh::GlobalEntityId{UINT64_C(50)},
        mesh::GlobalEntityId{UINT64_C(60)}};
    return {
        std::move(ids),
        {}};
}

[[nodiscard]] mesh::PartitionSnapshot
make_partition(
    int rank) {
    auto topology =
        make_topology();
    mesh::EntityOwnerRanks owners;
    owners.faces = {
        mesh::PartitionRank{0U},
        mesh::PartitionRank{0U},
        mesh::PartitionRank{0U},
        mesh::PartitionRank{0U}};
    owners.cells = {
        mesh::PartitionRank{0U},
        mesh::PartitionRank{1U},
        mesh::PartitionRank{0U},
        mesh::PartitionRank{1U},
        mesh::PartitionRank{0U},
        mesh::PartitionRank{1U}};
    return mesh::PartitionSnapshot::create(
        topology,
        mesh::PartitionRank{
            static_cast<
                mesh::PartitionRank::value_type>(
                    rank)},
        2U,
        std::move(owners));
}

[[nodiscard]]
dp::AssemblyReadyInternalConnectionRow3D
connection(
    std::size_t face,
    std::size_t owner,
    std::size_t neighbour,
    std::uint64_t face_global,
    std::uint64_t owner_global,
    std::uint64_t neighbour_global) {
    return {
        mesh::LocalIndex{
            static_cast<
                mesh::LocalIndex::value_type>(
                    face)},
        mesh::GlobalEntityId{
            face_global},
        mesh::LocalIndex{
            static_cast<
                mesh::LocalIndex::value_type>(
                    owner)},
        mesh::GlobalEntityId{
            owner_global},
        mesh::LocalIndex{
            static_cast<
                mesh::LocalIndex::value_type>(
                    neighbour)},
        mesh::GlobalEntityId{
            neighbour_global},
        2.0e-12};
}

[[nodiscard]]
dp::ParallelOwnedConnectionSchedule3D
make_schedule(
    int rank,
    bool include_cross_cardinality) {
    std::vector<
        dp::AssemblyReadyInternalConnectionRow3D>
        rows{
            connection(
                0U, 0U, 1U,
                UINT64_C(100),
                UINT64_C(10),
                UINT64_C(20)),
            connection(
                1U, 2U, 3U,
                UINT64_C(200),
                UINT64_C(30),
                UINT64_C(40)),
            connection(
                2U, 4U, 5U,
                UINT64_C(300),
                UINT64_C(50),
                UINT64_C(60))};
    if (include_cross_cardinality) {
        rows.push_back(
            connection(
                3U, 1U, 2U,
                UINT64_C(400),
                UINT64_C(20),
                UINT64_C(30)));
    }

    if (rank == 0) {
        return {
            mesh::PartitionRank{0U},
            2U,
            std::move(rows),
            {}};
    }
    return {
        mesh::PartitionRank{1U},
        2U,
        {},
        std::move(rows)};
}

[[nodiscard]]
dp::PetscMpiAijSymbolicPreallocation3D
make_cell_bridge(
    int rank) {
    const std::vector<PetscInt>
        local_cell_rows{
            0, 3, 1, 4, 2, 5};
    if (rank == 0) {
        return {
            mesh::PartitionRank{0U},
            2U,
            0,
            3,
            6,
            {
                mesh::LocalIndex{0U},
                mesh::LocalIndex{2U},
                mesh::LocalIndex{4U}},
            {
                mesh::GlobalEntityId{UINT64_C(10)},
                mesh::GlobalEntityId{UINT64_C(30)},
                mesh::GlobalEntityId{UINT64_C(50)}},
            {0, 1, 2},
            {1, 1, 1},
            {1, 2, 1},
            local_cell_rows};
    }
    return {
        mesh::PartitionRank{1U},
        2U,
        3,
        6,
        6,
        {
            mesh::LocalIndex{1U},
            mesh::LocalIndex{3U},
            mesh::LocalIndex{5U}},
        {
            mesh::GlobalEntityId{UINT64_C(20)},
            mesh::GlobalEntityId{UINT64_C(40)},
            mesh::GlobalEntityId{UINT64_C(60)}},
        {3, 4, 5},
        {1, 1, 1},
        {2, 1, 1},
        local_cell_rows};
}

[[nodiscard]]
dp::OwnedCellStructuralColumnPatternSnapshot3D
make_cell_pattern(
    int rank) {
    if (rank == 0) {
        return {
            mesh::PartitionRank{0U},
            2U,
            6U,
            0,
            3,
            6,
            {
                {
                    mesh::LocalIndex{0U},
                    mesh::GlobalEntityId{UINT64_C(10)},
                    0, 0U, 1U, 0U, 1U},
                {
                    mesh::LocalIndex{2U},
                    mesh::GlobalEntityId{UINT64_C(30)},
                    1, 1U, 1U, 1U, 2U},
                {
                    mesh::LocalIndex{4U},
                    mesh::GlobalEntityId{UINT64_C(50)},
                    2, 2U, 1U, 3U, 1U}},
            {0, 1, 2},
            {3, 3, 4, 5}};
    }
    return {
        mesh::PartitionRank{1U},
        2U,
        6U,
        3,
        6,
        6,
        {
            {
                mesh::LocalIndex{1U},
                mesh::GlobalEntityId{UINT64_C(20)},
                3, 0U, 1U, 0U, 2U},
            {
                mesh::LocalIndex{3U},
                mesh::GlobalEntityId{UINT64_C(40)},
                4, 1U, 1U, 2U, 1U},
            {
                mesh::LocalIndex{5U},
                mesh::GlobalEntityId{UINT64_C(60)},
                5, 2U, 1U, 3U, 1U}},
        {3, 4, 5},
        {0, 1, 1, 2}};
}

[[nodiscard]]
disc::CombinedTransmissibilityAdmissibility3D
direct_admissibility() {
    return {
        disc::
            CombinedTransmissibilityAdmissibilityDisposition3D::
                direct_normal_projection_k_orthogonal_candidate,
        {
            disc::
                TransmissibilityGeometryDisposition3D::
                    direct_normal_projection_allowed,
            0.0,
            0.05},
        {
            disc::
                KOrthogonalityDisposition3D::
                    k_orthogonal_within_policy,
            std::nullopt,
            std::nullopt,
            {
                std::nullopt,
                std::nullopt},
            0.05}};
}

[[nodiscard]]
disc::TpfaInternalFaceTransmissibilityEntry3D
materialized_entry(
    std::size_t face) {
    return {
        mesh::LocalIndex{
            static_cast<
                mesh::LocalIndex::value_type>(
                    face)},
        disc::
            TpfaInternalFaceTransmissibilityDisposition3D::
                materialized,
        direct_admissibility(),
        disc::TpfaStaticFaceTransmissibility3D{
            disc::
                TpfaStaticFaceTransmissibilityDisposition3D::
                    positive_harmonic_combination,
            1.0,
            2.0e-12}};
}

[[nodiscard]]
fdp::MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D
make_face_input(
    std::size_t face,
    std::uint64_t face_global) {
    return {
        mesh::LocalIndex{
            static_cast<
                mesh::LocalIndex::value_type>(
                    face)},
        mesh::GlobalEntityId{
            face_global},
        materialized_entry(face),
        flow::GravityVector3D{
            0.0, 0.0, 0.0},
        flow::OwnerToNeighbourDisplacement3D{
            1.0, 0.0, 0.0},
        mpmc::flow_discretization::
            StaticThermalFaceConductance3D{
                0.0}};
}

[[nodiscard]]
std::vector<
    fdp::
        MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D>
make_face_inputs(
    int rank,
    bool include_cross_cardinality) {
    if (rank != 0) {
        return {};
    }
    std::vector<
        fdp::
            MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D>
        result{
            make_face_input(
                0U,
                UINT64_C(100)),
            make_face_input(
                1U,
                UINT64_C(200)),
            make_face_input(
                2U,
                UINT64_C(300))};
    if (include_cross_cardinality) {
        result.push_back(
            make_face_input(
                3U,
                UINT64_C(400)));
    }
    return result;
}

[[nodiscard]]
fdp::MixedCardinalityPhysicalCurrentCellLinearization3D
evaluate_target(
    std::uint64_t stable,
    DispatchAudit* audit) {
    const std::vector<std::string>
        ids{"A", "B", "C"};
    const auto q =
        target_state(stable);
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;

    if (phase_count_for(stable) == 1U) {
        std::optional<
            fdp::SinglePhaseCurrentCellLinearization3D>
            output;
        auto layout = layout_1p();
        const auto error =
            evaluate_1p(
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{stable},
                q,
                layout,
                ids,
                audit,
                &output,
                &status);
        if (error != PETSC_SUCCESS ||
            status !=
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success ||
            !output.has_value()) {
            throw std::runtime_error(
                "failed direct 1P target evaluation");
        }
        return std::move(*output);
    }

    if (phase_count_for(stable) == 2U) {
        std::optional<
            fdp::TwoPhaseCurrentCellLinearization3D>
            output;
        auto layout = layout_2p();
        const auto error =
            evaluate_2p(
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{stable},
                q,
                layout,
                ids,
                audit,
                &output,
                &status);
        if (error != PETSC_SUCCESS ||
            status !=
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success ||
            !output.has_value()) {
            throw std::runtime_error(
                "failed direct 2P target evaluation");
        }
        return std::move(*output);
    }

    std::optional<
        fdp::FixedThreePhaseCurrentCellLinearization3D>
        output;
    auto layout = layout_3p();
    const auto error =
        evaluate_3p(
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{stable},
            q,
            layout,
            ids,
            audit,
            &output,
            &status);
    if (error != PETSC_SUCCESS ||
        status !=
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success ||
        !output.has_value()) {
        throw std::runtime_error(
            "failed direct 3P target evaluation");
    }
    return std::move(*output);
}

[[nodiscard]]
std::vector<
    fdp::MixedCardinalityPhysicalSnesCellInput3D>
make_cell_inputs(
    int rank,
    DispatchAudit* audit) {
    std::vector<
        fdp::MixedCardinalityPhysicalSnesCellInput3D>
        result;
    result.reserve(6U);

    for (std::size_t local = 0U;
         local < 6U;
         ++local) {
        const std::uint64_t stable =
            static_cast<std::uint64_t>(
                (local + 1U) * 10U);
        const bool owned =
            (rank == 0 &&
             local % 2U == 0U) ||
            (rank == 1 &&
             local % 2U == 1U);
        const double bulk_volume =
            2.0 +
            static_cast<double>(local);
        const double porosity =
            0.25 +
            0.005 *
                static_cast<double>(
                    local);
        auto current =
            evaluate_target(
                stable,
                audit);

        if (phase_count_for(stable) == 1U) {
            fdp::SinglePhaseSnesCellInput3D
                input;
            input.cell =
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            local)};
            input.cell_global =
                mesh::GlobalEntityId{stable};
            input.bulk_volume_m3 =
                bulk_volume;
            input.porosity =
                porosity;
            input.frozen_layout =
                layout_1p();
            input.component_ids =
                {"A", "B", "C"};
            if (owned) {
                const auto& typed =
                    std::get<
                        fdp::
                            SinglePhaseCurrentCellLinearization3D>(
                                current);
                input.previous_component_accumulation =
                    flow::
                        build_single_phase_component_accumulation(
                            typed.state,
                            porosity);
                input.previous_energy_accumulation =
                    flow::
                        build_single_phase_energy_accumulation_snapshot(
                            typed.state,
                            porosity,
                            typed.transport,
                            typed.caloric,
                            typed.rock);
            }
            result.emplace_back(
                std::move(input));
            continue;
        }

        if (phase_count_for(stable) == 2U) {
            fdp::TwoPhaseSnesCellInput3D
                input;
            input.cell =
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            local)};
            input.cell_global =
                mesh::GlobalEntityId{stable};
            input.bulk_volume_m3 =
                bulk_volume;
            input.porosity =
                porosity;
            input.frozen_layout =
                layout_2p();
            input.component_ids =
                {"A", "B", "C"};
            if (owned) {
                const auto& typed =
                    std::get<
                        fdp::
                            TwoPhaseCurrentCellLinearization3D>(
                                current);
                input.previous_component_accumulation =
                    flow::
                        build_two_phase_component_accumulation(
                            typed.state,
                            porosity);
                input.previous_energy_accumulation =
                    flow::
                        build_two_phase_energy_accumulation_snapshot(
                            typed.state,
                            porosity,
                            typed.transport,
                            typed.caloric,
                            typed.rock);
            }
            result.emplace_back(
                std::move(input));
            continue;
        }

        fdp::FixedThreePhaseSnesCellInput3D
            input;
        input.cell =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)};
        input.cell_global =
            mesh::GlobalEntityId{stable};
        input.bulk_volume_m3 =
            bulk_volume;
        input.porosity =
            porosity;
        input.frozen_layout =
            layout_3p();
        input.component_ids =
            {"A", "B", "C"};
        if (owned) {
            const auto& typed =
                std::get<
                    fdp::
                        FixedThreePhaseCurrentCellLinearization3D>(
                            current);
            input.previous_component_accumulation =
                flow::
                    build_pore_volume_component_accumulation(
                        typed.state,
                        porosity);
            input.previous_energy_accumulation =
                flow::
                    build_pore_volume_energy_accumulation_snapshot(
                        typed.state,
                        porosity,
                        typed.transport,
                        typed.caloric,
                        typed.rock);
        }
        result.emplace_back(
            std::move(input));
    }
    return result;
}

[[nodiscard]]
std::vector<
    flow::FrozenActivePhaseIdentityMap>
make_phase_identity_maps() {
    const flow::FrozenPhysicalPhaseIdentity
        aqueous{
            "fixture/mixed-physical-dispatch",
            "aqueous"};
    const flow::FrozenPhysicalPhaseIdentity
        hydrocarbon0{
            "fixture/mixed-physical-dispatch",
            "hydrocarbon-0"};
    const flow::FrozenPhysicalPhaseIdentity
        hydrocarbon1{
            "fixture/mixed-physical-dispatch",
            "hydrocarbon-1"};

    return {
        flow::FrozenActivePhaseIdentityMap{
            {aqueous}},
        flow::FrozenActivePhaseIdentityMap{
            {aqueous}},
        flow::FrozenActivePhaseIdentityMap{
            {
                aqueous,
                hydrocarbon0}},
        flow::FrozenActivePhaseIdentityMap{
            {
                aqueous,
                hydrocarbon0}},
        flow::FrozenActivePhaseIdentityMap{
            {
                aqueous,
                hydrocarbon0,
                hydrocarbon1}},
        flow::FrozenActivePhaseIdentityMap{
            {
                aqueous,
                hydrocarbon0,
                hydrocarbon1}}};
}


[[nodiscard]]
flow::FrozenActivePhaseIdentityMap
phase_identity_map_for_stable(
    std::uint64_t stable) {
    const auto maps =
        make_phase_identity_maps();
    const std::size_t index =
        static_cast<std::size_t>(
            stable / UINT64_C(10) -
            UINT64_C(1));
    return maps.at(index);
}


th::Provenance thermodynamic_adapter_source(
    std::string locator) {
    return {
        th::SourceKind::synthetic_test,
        "MPMC_HNU thermodynamic cross-cardinality adapter regression",
        "v1",
        std::move(locator),
        "Manufactured values verify identity-to-selected-branch TPFA integration only",
        "tests/flow_discretization/petsc/mixed_cardinality_physical_snes_assembly_test.cpp",
        "Repository structural regression"};
}

th::SourcedScalar thermodynamic_adapter_scalar(
    double value,
    th::Unit unit,
    std::string locator) {
    return {
        value,
        unit,
        thermodynamic_adapter_source(
            std::move(locator)),
        "SI",
        "identity"};
}

th::Component thermodynamic_adapter_component(
    std::string id,
    double molar_mass_kg_per_mol) {
    const auto definition =
        thermodynamic_adapter_source(
            "component-" + id);
    return {
        id,
        id,
        th::ComponentKind::pure,
        definition,
        thermodynamic_adapter_scalar(
            molar_mass_kg_per_mol,
            th::Unit::kilogram_per_mole,
            "molar-mass-" + id)};
}

th::PrParameterSet
thermodynamic_adapter_pr_parameters() {
    std::vector<th::Component> catalog{
        thermodynamic_adapter_component(
            "A", 0.020),
        thermodynamic_adapter_component(
            "B", 0.030),
        thermodynamic_adapter_component(
            "C", 0.040)};
    const std::vector<std::string>
        order{"A", "B", "C"};

    th::PrParameterInput input;
    input.model_id =
        std::string{th::pr76_profile};
    input.dataset_id =
        "synthetic-cross-cardinality-pr76";
    input.revision = "v1";
    input.applicability = {
        std::nullopt,
        std::nullopt,
        thermodynamic_adapter_source(
            "applicability")};

    const std::array<double, 3>
        critical_temperature{
            20.0, 25.0, 30.0};
    const std::array<double, 3>
        critical_pressure{
            1.0e6, 1.2e6, 1.5e6};
    const std::array<double, 3>
        acentric_factor{
            0.0, 0.1, 0.2};
    for (std::size_t index = 0U;
         index < order.size();
         ++index) {
        input.pure.push_back({
            order[index],
            thermodynamic_adapter_scalar(
                critical_temperature[index],
                th::Unit::kelvin,
                order[index] + "-Tc"),
            thermodynamic_adapter_scalar(
                critical_pressure[index],
                th::Unit::pascal,
                order[index] + "-Pc"),
            thermodynamic_adapter_scalar(
                acentric_factor[index],
                th::Unit::dimensionless,
                order[index] + "-omega")});
    }
    for (std::size_t first = 0U;
         first < order.size();
         ++first) {
        for (std::size_t second =
                 first + 1U;
             second < order.size();
             ++second) {
            input.binary.push_back({
                order[first],
                order[second],
                thermodynamic_adapter_scalar(
                    0.0,
                    th::Unit::dimensionless,
                    order[first] + "-" +
                        order[second] +
                        "-kij")});
        }
    }
    return th::PrParameterSet::create(
        catalog,
        order,
        input,
        th::DataPolicy::
            allow_synthetic_tests);
}

flow::FrozenPhysicalPhaseIdentity
mixed_physical_phase_identity(
    std::string key) {
    return {
        "fixture/mixed-physical-dispatch",
        std::move(key)};
}

flow::PhaseIdentityContinuationSnapshot
thermodynamic_adapter_continuation() {
    flow::PhaseSetTransitionCandidate
        candidate;
    candidate.source_phase_count = 1U;
    candidate.target_phase_count = 2U;
    candidate.trigger =
        flow::
            PhaseSetTransitionTrigger::
                stability_witness;
    candidate.status =
        flow::
            PhaseSetTransitionCandidateStatus::
                target_resolved;
    candidate.evidence_profile =
        "fixture/cross-cardinality-fresh-target/v1";
    candidate.diagnostic =
        "explicit 1P->2P phase identity continuation";

    return flow::
        make_phase_identity_continuation_snapshot(
            candidate,
            flow::FrozenActivePhaseIdentityMap{
                {
                    mixed_physical_phase_identity(
                        "aqueous")}},
            flow::FrozenActivePhaseIdentityMap{
                {
                    mixed_physical_phase_identity(
                        "aqueous"),
                    mixed_physical_phase_identity(
                        "hydrocarbon-0")}});
}

flow::Pr76AbsentPhasePotentialExtensionProvider<double>
make_thermodynamic_adapter_pr_provider(
    const th::Pr76Phase<double>& model) {
    const auto continuation =
        thermodynamic_adapter_continuation();
    auto registry =
        flow::
            make_transition_selected_phase_branch_registry(
                continuation,
                std::vector<
                    flow::
                        FrozenSelectedPhaseBranchBinding<
                            th::Pr76SelectedPhase>>{
                    {
                        mixed_physical_phase_identity(
                            "aqueous"),
                        {0U, {}},
                        "fixture/aqueous-root0"},
                    {
                        mixed_physical_phase_identity(
                            "hydrocarbon-0"),
                        {0U, {}},
                        "fixture/hydrocarbon0-root0"}});
    return {
        model,
        std::move(registry)};
}

struct ThermodynamicCoordinateResolverAudit {
    std::uint64_t calls{};
};

PetscErrorCode
resolve_absent_phase_thermodynamic_coordinates(
    const fdp::
        MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D&,
    const flow::
        CrossCardinalityFacePhaseIdentityBinding&
            phase_binding,
    fdp::CrossCardinalityAbsentPhaseSide3D
        absent_side,
    const fdp::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            owner,
    const fdp::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            neighbour,
    void* raw_context,
    std::optional<
        flow::
            AbsentPhaseThermodynamicCoordinateExtension>*
                output,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (raw_context == nullptr ||
        output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    auto* audit =
        static_cast<
            ThermodynamicCoordinateResolverAudit*>(
                raw_context);
    ++audit->calls;

    try {
        if (phase_binding
                .identity.opaque_phase_key !=
                "hydrocarbon-0" ||
            absent_side !=
                fdp::
                    CrossCardinalityAbsentPhaseSide3D::
                        owner) {
            return PETSC_ERR_ARG_INCOMP;
        }

        const auto host_identity =
            std::visit(
                [](const auto& typed) {
                    return typed.transport
                        .state_identity;
                },
                owner);
        (void)neighbour;

        const std::size_t q =
            host_identity.layout
                .unknown_count();
        const std::size_t n =
            host_identity.layout
                .component_count();
        if (n != 3U) {
            return PETSC_ERR_ARG_INCOMP;
        }

        flow::
            AbsentPhaseThermodynamicCoordinateExtension
            coordinates;
        coordinates.identity =
            phase_binding.identity;
        coordinates.host_state_identity =
            host_identity;
        coordinates.phase_pressure_pa =
            host_identity
                .reference_pressure_pa;
        coordinates.phase_pressure_gradient
            .assign(q, 0.0);
        coordinates.phase_pressure_gradient[
            host_identity.layout
                .pressure_unknown_index()] =
            1.0;
        coordinates.hypothetical_composition =
            {0.40, 0.20, 0.40};
        coordinates
            .hypothetical_composition_jacobian
            .assign(n * q, 0.0);
        coordinates.provenance =
            "fixture/stability-hypothetical-composition/v1";
        coordinates.validate();

        output->emplace(
            std::move(coordinates));
        *status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        return PETSC_ERR_ARG_INCOMP;
    } catch (...) {
        return PETSC_ERR_LIB;
    }
}

PetscErrorCode
evaluate_absent_phase_extension(
    const fdp::
        MixedCardinalityPhysicalSnesAuthoritativeFaceInput3D&,
    const flow::
        CrossCardinalityFacePhaseIdentityBinding&
            phase_binding,
    fdp::CrossCardinalityAbsentPhaseSide3D
        absent_side,
    const fdp::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            owner,
    const fdp::
        MixedCardinalityPhysicalCurrentCellLinearization3D&
            neighbour,
    void*,
    std::optional<
        flow::
            AbsentPhasePotentialExtensionLinearization>*
                output,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (output == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();

    try {
        const auto& absent_cell =
            absent_side ==
                    fdp::
                        CrossCardinalityAbsentPhaseSide3D::
                            owner
                ? owner
                : neighbour;
        const auto identity =
            std::visit(
                [](const auto& typed) {
                    return typed.transport.state_identity;
                },
                absent_cell);
        const std::size_t q =
            identity.layout.unknown_count();

        std::vector<double>
            pressure_gradient(
                q,
                0.0);
        pressure_gradient[
            identity.layout
                .pressure_unknown_index()] =
            1.0;
        std::vector<double>
            density_gradient(
                q,
                0.0);

        double density =
            1.0;
        if (phase_binding
                .identity.opaque_phase_key ==
            "hydrocarbon-0") {
            density = 1.5;
        } else if (
            phase_binding
                .identity.opaque_phase_key ==
            "hydrocarbon-1") {
            density = 1.8;
        }

        output->emplace(
            phase_binding.identity,
            identity.reference_pressure_pa,
            density,
            q,
            std::move(
                pressure_gradient),
            std::move(
                density_gradient),
            "fixture/cross-cardinality-potential-extension/v1");
        *status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
        return PETSC_SUCCESS;
    } catch (const std::invalid_argument&) {
        *status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
        return PETSC_SUCCESS;
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }
}

void check_direct_cross_cardinality_bridge(
    std::uint64_t owner_stable,
    std::uint64_t neighbour_stable,
    std::size_t face_index,
    DispatchAudit* audit) {
    const auto owner =
        evaluate_target(
            owner_stable,
            audit);
    const auto neighbour =
        evaluate_target(
            neighbour_stable,
            audit);
    const auto plan =
        flow::
            make_cross_cardinality_face_phase_identity_plan(
                phase_identity_map_for_stable(
                    owner_stable),
                phase_identity_map_for_stable(
                    neighbour_stable));

    fdp::CrossCardinalityTpfaBridgeBinding3D
        bridge_binding{
            &evaluate_absent_phase_extension,
            nullptr};
    std::optional<
        fdp::
            MixedCardinalityPhysicalFaceLinearization3D>
        face;
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;

    const auto face_input =
        make_face_input(
            face_index,
            UINT64_C(1000) +
                static_cast<std::uint64_t>(
                    face_index));
    const PetscErrorCode error =
        fdp::
            evaluate_standard_cross_cardinality_tpfa_face_3d(
                face_input,
                plan,
                {2.0, 3.0},
                owner,
                neighbour,
                &bridge_binding,
                &face,
                &status);

    require_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            face.has_value(),
        "standard cross-cardinality TPFA bridge evaluation failed");

    const auto& component =
        face->component;
    const auto& energy =
        face->energy;
    require_collective(
        component.owner_state_identity
                    .layout.phase_count() ==
                phase_count_for(
                    owner_stable) &&
            component.neighbour_state_identity
                    .layout.phase_count() ==
                phase_count_for(
                    neighbour_stable) &&
            component
                    .owner_row_neighbour_column_jacobian
                    .size() ==
                component.component_count() *
                    component
                        .neighbour_state_identity
                        .layout.unknown_count() &&
            component
                    .neighbour_row_owner_column_jacobian
                    .size() ==
                component.component_count() *
                    component
                        .owner_state_identity
                        .layout.unknown_count(),
        "cross-cardinality TPFA rectangular Jacobian shape mismatch");

    double component_norm = 0.0;
    for (const double value :
         component
             .owner_component_contribution_mol_per_bulk_m3_s) {
        component_norm +=
            std::abs(value);
    }
    require_collective(
        component_norm > 0.0 &&
            std::abs(
                energy
                    .owner_contribution_w_per_bulk_m3) >
                0.0,
        "cross-cardinality bridge did not produce physical component/energy flux");

    for (std::size_t component_index = 0U;
         component_index <
            component.component_count();
         ++component_index) {
        const double owner_rate =
            component
                .owner_component_contribution_mol_per_bulk_m3_s[
                    component_index] *
            component.bulk_volume
                .owner_bulk_volume_m3;
        const double neighbour_rate =
            component
                .neighbour_component_contribution_mol_per_bulk_m3_s[
                    component_index] *
            component.bulk_volume
                .neighbour_bulk_volume_m3;
        require_collective(
            std::abs(
                owner_rate +
                neighbour_rate) <
                1.0e-12 *
                    std::max(
                        1.0,
                        std::abs(
                            owner_rate)),
            "cross-cardinality component face rate is not conservative");
    }

    const double owner_energy =
        energy
            .owner_contribution_w_per_bulk_m3 *
        energy.bulk_volume
            .owner_bulk_volume_m3;
    const double neighbour_energy =
        energy
            .neighbour_contribution_w_per_bulk_m3 *
        energy.bulk_volume
            .neighbour_bulk_volume_m3;
    require_collective(
        std::abs(
            owner_energy +
            neighbour_energy) <
            1.0e-12 *
                std::max(
                    1.0,
                    std::abs(
                        owner_energy)),
        "cross-cardinality energy face rate is not conservative");
}

void insert_target(
    Vec state,
    const fdp::
        VariableCardinalityNaturalVariableNumbering3D&
            numbering) {
    bool local_ok = true;
    for (const auto& cell :
         numbering.cells()) {
        if (cell.owner_rank !=
            numbering.local_rank()) {
            continue;
        }
        const auto target =
            target_state(
                cell.cell_global.value());
        local_ok =
            local_ok &&
            target.size() ==
                cell.scalar_count;
        for (std::size_t slot = 0U;
             slot < cell.scalar_count;
             ++slot) {
            const PetscInt index =
                cell.petsc_global_scalar_start +
                static_cast<PetscInt>(
                    slot);
            const PetscScalar value =
                static_cast<PetscScalar>(
                    target[slot]);
            local_ok =
                local_ok &&
                VecSetValues(
                    state,
                    1,
                    &index,
                    &value,
                    INSERT_VALUES) ==
                    PETSC_SUCCESS;
        }
    }
    require_collective(
        local_ok,
        "failed to insert mixed physical target state");
    require_collective(
        VecAssemblyBegin(state) ==
                PETSC_SUCCESS &&
            VecAssemblyEnd(state) ==
                PETSC_SUCCESS,
        "failed to assemble mixed physical target state");
}


fdp::FrozenAbsentPhaseCoordinateEntry3D
make_frozen_absent_coordinate_entry(
    std::uint64_t stable,
    flow::FrozenPhysicalPhaseIdentity identity,
    std::vector<double> composition,
    std::string branch_provenance,
    DispatchAudit* audit) {
    const auto current =
        evaluate_target(
            stable,
            audit);
    const auto host_identity =
        std::visit(
            [](const auto& typed) {
                return typed.transport
                    .state_identity;
            },
            current);
    const std::size_t q =
        host_identity.layout
            .unknown_count();
    const std::size_t n =
        host_identity.layout
            .component_count();
    if (composition.size() != n) {
        throw std::invalid_argument(
            "frozen absent coordinate composition size mismatch");
    }

    flow::
        AbsentPhaseThermodynamicCoordinateExtension
        coordinates;
    coordinates.identity =
        identity;
    coordinates.host_state_identity =
        host_identity;
    coordinates.phase_pressure_pa =
        host_identity
            .reference_pressure_pa;
    coordinates.phase_pressure_gradient
        .assign(q, 0.0);
    coordinates.phase_pressure_gradient[
        host_identity.layout
            .pressure_unknown_index()] =
        1.0;
    coordinates.hypothetical_composition =
        std::move(composition);
    coordinates
        .hypothetical_composition_jacobian
        .assign(n * q, 0.0);
    coordinates.provenance =
        "fixture/frozen-absent-coordinate-reference/v1";
    coordinates.validate();

    return {
        std::move(identity),
        std::move(coordinates),
        std::move(branch_provenance)};
}

flow::Pr76AbsentPhasePotentialExtensionProvider<double>
make_outer_rebuild_pr_provider(
    const th::Pr76Phase<double>& model) {
    flow::PhaseSetTransitionCandidate
        candidate;
    candidate.source_phase_count = 1U;
    candidate.target_phase_count = 3U;
    candidate.status =
        flow::
            PhaseSetTransitionCandidateStatus::
                target_resolved;
    candidate.trigger =
        flow::
            PhaseSetTransitionTrigger::
                stability_witness;
    candidate.evidence_profile =
        "fixture/outer-rebuild-phase-universe/v1";

    const auto continuation =
        flow::
            make_phase_identity_continuation_snapshot(
                candidate,
                flow::FrozenActivePhaseIdentityMap{
                    {
                        mixed_physical_phase_identity(
                            "aqueous")}},
                flow::FrozenActivePhaseIdentityMap{
                    {
                        mixed_physical_phase_identity(
                            "aqueous"),
                        mixed_physical_phase_identity(
                            "hydrocarbon-0"),
                        mixed_physical_phase_identity(
                            "hydrocarbon-1")}});

    auto registry =
        flow::
            make_transition_selected_phase_branch_registry(
                continuation,
                std::vector<
                    flow::
                        FrozenSelectedPhaseBranchBinding<
                            th::Pr76SelectedPhase>>{
                    {
                        mixed_physical_phase_identity(
                            "aqueous"),
                        {0U, {}},
                        "fixture/aqueous-root0"},
                    {
                        mixed_physical_phase_identity(
                            "hydrocarbon-0"),
                        {0U, {}},
                        "fixture/hydrocarbon0-root0"},
                    {
                        mixed_physical_phase_identity(
                            "hydrocarbon-1"),
                        {0U, {}},
                        "fixture/hydrocarbon1-root0"}});

    return {
        model,
        std::move(registry)};
}

std::vector<
    fdp::FrozenPhaseTransitionRebuildCell3D>
make_outer_rebuild_cells(
    int rank,
    DispatchAudit* audit) {
    auto base_inputs =
        make_cell_inputs(
            rank,
            audit);
    auto phase_maps =
        make_phase_identity_maps();

    std::vector<
        fdp::
            FrozenPhaseTransitionRebuildCell3D>
        result;
    result.reserve(
        base_inputs.size());

    for (std::size_t local = 0U;
         local < base_inputs.size();
         ++local) {
        const std::uint64_t stable =
            static_cast<std::uint64_t>(
                (local + 1U) * 10U);

        fdp::FrozenPhaseTransitionRebuildCell3D
            snapshot;
        std::visit(
            [&](const auto& typed) {
                snapshot.cell =
                    typed.cell;
                snapshot.cell_global =
                    typed.cell_global;
                snapshot.bulk_volume_m3 =
                    typed.bulk_volume_m3;
                snapshot.porosity =
                    typed.porosity;
                snapshot.component_ids =
                    typed.component_ids;
                using Typed =
                    std::remove_cvref_t<
                        decltype(typed)>;
                if constexpr (
                    std::is_same_v<
                        Typed,
                        fdp::
                            FixedThreePhaseSnesCellInput3D>) {
                    snapshot.target_layout =
                        flow::
                            NaturalVariableLayoutDescriptor{
                                typed.frozen_layout};
                } else {
                    snapshot.target_layout =
                        typed.frozen_layout
                            .descriptor();
                }
                snapshot.previous_component_accumulation =
                    typed.previous_component_accumulation;
                snapshot.previous_energy_accumulation =
                    typed.previous_energy_accumulation;
            },
            base_inputs[local]);

        snapshot.target_natural_variables =
            target_state(stable);
        snapshot.target_active_phases =
            phase_maps[local];
        snapshot.transition_evidence_profile =
            "fixture/frozen-active-set/no-transition/v1";

        if (stable == UINT64_C(20)) {
            snapshot.frozen_absent_phases.push_back(
                make_frozen_absent_coordinate_entry(
                    stable,
                    mixed_physical_phase_identity(
                        "hydrocarbon-0"),
                    {0.40, 0.20, 0.40},
                    "fixture/hydrocarbon0-root0",
                    audit));
            snapshot.frozen_absent_phases.push_back(
                make_frozen_absent_coordinate_entry(
                    stable,
                    mixed_physical_phase_identity(
                        "hydrocarbon-1"),
                    {0.10, 0.50, 0.40},
                    "fixture/hydrocarbon1-root0",
                    audit));
        } else if (
            stable == UINT64_C(40)) {
            snapshot.frozen_absent_phases.push_back(
                make_frozen_absent_coordinate_entry(
                    stable,
                    mixed_physical_phase_identity(
                        "hydrocarbon-1"),
                    {0.10, 0.50, 0.40},
                    "fixture/hydrocarbon1-root0",
                    audit));
        }

        result.push_back(
            std::move(snapshot));
    }

    // Accept a real 2P -> 3P transition on stable cell30.  The target phase
    // compositions are identical to the source overall inventory so the
    // transition projection passes the existing material-balance gate.
    constexpr std::size_t transitioned_local = 2U;
    constexpr double porosity = 0.26;
    const auto source_current =
        evaluate_target(
            UINT64_C(30),
            audit);
    const auto& source_two_phase =
        std::get<
            fdp::
                TwoPhaseCurrentCellLinearization3D>(
                    source_current);
    const auto source_component =
        flow::
            build_two_phase_component_accumulation(
                source_two_phase.state,
                porosity);
    const auto source_energy =
        flow::
            build_two_phase_energy_accumulation_snapshot(
                source_two_phase.state,
                porosity,
                source_two_phase.transport,
                source_two_phase.caloric,
                source_two_phase.rock);

    double inventory_total = 0.0;
    for (const double value :
         source_component
             .component_accumulation_mol_per_bulk_m3) {
        inventory_total += value;
    }
    std::vector<double>
        overall;
    overall.reserve(3U);
    for (const double value :
         source_component
             .component_accumulation_mol_per_bulk_m3) {
        overall.push_back(
            value / inventory_total);
    }

    flow::PhaseSetTransitionCandidate
        candidate;
    candidate.source_phase_count = 2U;
    candidate.target_phase_count = 3U;
    candidate.trigger =
        flow::
            PhaseSetTransitionTrigger::
                stability_witness;
    candidate.status =
        flow::
            PhaseSetTransitionCandidateStatus::
                target_resolved;
    candidate.pressure_pa = 20.0;
    candidate.temperature_k = 9.0;
    candidate.component_ids =
        {"A", "B", "C"};
    candidate.target_phases = {
        {0.20, overall, 2.0},
        {0.30, overall, 4.0},
        {0.50, overall, 7.0}};
    candidate.evidence_profile =
        "fixture/accepted-2p-to-3p/v1";
    candidate.diagnostic =
        "controlled outer rebuild transition";

    result[transitioned_local] =
        fdp::
            make_accepted_phase_transition_rebuild_cell_3d(
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            transitioned_local)},
                mesh::GlobalEntityId{
                    UINT64_C(30)},
                4.0,
                porosity,
                candidate,
                source_component,
                source_component,
                source_energy,
                phase_maps[transitioned_local],
                flow::FrozenActivePhaseIdentityMap{
                    {
                        mixed_physical_phase_identity(
                            "aqueous"),
                        mixed_physical_phase_identity(
                            "hydrocarbon-0"),
                        mixed_physical_phase_identity(
                            "hydrocarbon-1")}},
                {});

    return result;
}


std::vector<double>
controller_transitioned_3p_state() {
    auto q =
        target_3p();
    const auto layout =
        layout_3p();
    q[layout.pressure_unknown_index()] =
        20.0;
    q[layout.temperature_unknown_index()] =
        9.0;
    return q;
}

std::vector<
    fdp::FrozenPhaseTransitionRebuildCell3D>
make_controller_cells(
    int rank,
    bool cell30_is_three_phase,
    DispatchAudit* audit) {
    auto base_inputs =
        make_cell_inputs(
            rank,
            audit);
    auto phase_maps =
        make_phase_identity_maps();

    std::vector<
        fdp::
            FrozenPhaseTransitionRebuildCell3D>
        result;
    result.reserve(
        base_inputs.size());

    for (std::size_t local = 0U;
         local < base_inputs.size();
         ++local) {
        const std::uint64_t stable =
            static_cast<std::uint64_t>(
                (local + 1U) * 10U);
        fdp::FrozenPhaseTransitionRebuildCell3D
            snapshot;
        std::visit(
            [&](const auto& typed) {
                using Typed =
                    std::remove_cvref_t<
                        decltype(typed)>;
                snapshot.cell =
                    typed.cell;
                snapshot.cell_global =
                    typed.cell_global;
                snapshot.bulk_volume_m3 =
                    typed.bulk_volume_m3;
                snapshot.porosity =
                    typed.porosity;
                snapshot.component_ids =
                    typed.component_ids;
                if constexpr (
                    std::is_same_v<
                        Typed,
                        fdp::
                            FixedThreePhaseSnesCellInput3D>) {
                    snapshot.target_layout =
                        flow::
                            NaturalVariableLayoutDescriptor{
                                typed.frozen_layout};
                } else {
                    snapshot.target_layout =
                        typed.frozen_layout
                            .descriptor();
                }
                snapshot.previous_component_accumulation =
                    typed.previous_component_accumulation;
                snapshot.previous_energy_accumulation =
                    typed.previous_energy_accumulation;
            },
            base_inputs[local]);

        snapshot.target_natural_variables =
            target_state(stable);
        snapshot.target_active_phases =
            phase_maps[local];
        snapshot.transition_evidence_profile =
            "fixture/controller-frozen-active-set/v1";

        if (stable == UINT64_C(40)) {
            snapshot.frozen_absent_phases.push_back(
                make_frozen_absent_coordinate_entry(
                    stable,
                    mixed_physical_phase_identity(
                        "hydrocarbon-1"),
                    {0.10, 0.50, 0.40},
                    "fixture/hydrocarbon1-root0",
                    audit));
        }

        result.push_back(
            std::move(snapshot));
    }

    if (!cell30_is_three_phase) {
        return result;
    }

    constexpr std::size_t local = 2U;
    const auto q =
        controller_transitioned_3p_state();
    const auto layout =
        layout_3p();

    // This controller-policy fixture isolates outer orchestration from the
    // already-covered BE history-migration contract.  Seed the rebuilt 3P
    // previous state from the same controlled 3P closure so the restarted
    // SNES is an exact zero-residual solve; make_accepted_phase_transition_
    // rebuild_cell_3d above separately verifies preservation of real 2P
    // histories through an accepted transition.
    std::optional<
        fdp::
            FixedThreePhaseCurrentCellLinearization3D>
        controlled_three_phase;
    fdp::NaturalVariableSnesEvaluationStatus3D
        controlled_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    const std::vector<std::string>
        controlled_ids{
            "A", "B", "C"};
    const PetscErrorCode controlled_error =
        evaluate_3p(
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)},
            mesh::GlobalEntityId{
                UINT64_C(30)},
            q,
            layout,
            controlled_ids,
            audit,
            &controlled_three_phase,
            &controlled_status);
    if (controlled_error != PETSC_SUCCESS ||
        controlled_status !=
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success ||
        !controlled_three_phase.has_value()) {
        throw std::runtime_error(
            "failed to build controlled 3P controller state");
    }

    fdp::FrozenPhaseTransitionRebuildCell3D
        transitioned;
    transitioned.cell =
        mesh::LocalIndex{
            static_cast<
                mesh::LocalIndex::value_type>(
                    local)};
    transitioned.cell_global =
        mesh::GlobalEntityId{
            UINT64_C(30)};
    transitioned.bulk_volume_m3 =
        4.0;
    transitioned.porosity =
        0.26;
    transitioned.component_ids =
        {"A", "B", "C"};
    transitioned.target_layout =
        flow::NaturalVariableLayoutDescriptor{
            layout};
    transitioned.target_natural_variables =
        q;
    transitioned.target_active_phases =
        flow::FrozenActivePhaseIdentityMap{
            {
                mixed_physical_phase_identity(
                    "aqueous"),
                mixed_physical_phase_identity(
                    "hydrocarbon-0"),
                mixed_physical_phase_identity(
                    "hydrocarbon-1")}};
    transitioned.transition_evidence_profile =
        "fixture/controller-accepted-2p-to-3p/v1";
    if (rank == 0) {
        transitioned
            .previous_component_accumulation =
            flow::
                build_pore_volume_component_accumulation(
                    controlled_three_phase->state,
                    transitioned.porosity);
        transitioned
            .previous_energy_accumulation =
            flow::
                build_pore_volume_energy_accumulation_snapshot(
                    controlled_three_phase->state,
                    transitioned.porosity,
                    controlled_three_phase->transport,
                    controlled_three_phase->caloric,
                    controlled_three_phase->rock);
    }

    result[local] =
        std::move(
            transitioned);
    return result;
}

struct ControllerFixture {
    ControllerFixture(
        int input_rank,
        const dp::
            ParallelOwnedConnectionSchedule3D*
                input_schedule,
        const mesh::PartitionSnapshot*
            input_partition,
        const dp::
            PetscMpiAijSymbolicPreallocation3D*
                input_bridge,
        const dp::
            OwnedCellStructuralColumnPatternSnapshot3D*
                input_pattern,
        DispatchAudit* input_audit,
        const th::Pr76Phase<double>*
            input_pr_model,
        flow::
            Pr76AbsentPhasePotentialExtensionProvider<
                double>*
                input_provider,
        bool input_oscillate_after_restart,
        std::size_t input_rebuild_calls,
        bool input_scan_indeterminate =
            false,
        bool input_phase_disappearance_only =
            false)
        : rank(input_rank),
          schedule(input_schedule),
          partition(input_partition),
          bridge(input_bridge),
          pattern(input_pattern),
          audit(input_audit),
          pr_model(input_pr_model),
          provider(input_provider),
          oscillate_after_restart(
              input_oscillate_after_restart),
          rebuild_calls(
              input_rebuild_calls),
          scan_indeterminate(
              input_scan_indeterminate),
          phase_disappearance_only(
              input_phase_disappearance_only) {}

    int rank{};
    const dp::
        ParallelOwnedConnectionSchedule3D*
            schedule{};
    const mesh::PartitionSnapshot*
        partition{};
    const dp::
        PetscMpiAijSymbolicPreallocation3D*
            bridge{};
    const dp::
        OwnedCellStructuralColumnPatternSnapshot3D*
            pattern{};
    DispatchAudit* audit{};
    const th::Pr76Phase<double>* pr_model{};
    flow::
        Pr76AbsentPhasePotentialExtensionProvider<
            double>*
            provider{};
    bool oscillate_after_restart{};
    std::size_t rebuild_calls{};
    bool scan_indeterminate{};
    bool phase_disappearance_only{};

    wdp::
        FixedBhpPeacemanWellSourceEvaluatorContext3D*
            well_context{};
    std::optional<
        wdp::
            FixedBhpPeacemanWellSourceEvaluatorContext3D>
        rebound_well_context;
    FixedBhpWellSourceAudit
        well_source_audit{};
    bool well_rebound{};

    wdp::
        FixedBhpMultiConnectionWellSourceEvaluatorContext3D*
            multi_well_context{};
    std::optional<
        wdp::
            FixedBhpMultiConnectionWellSourceEvaluatorContext3D>
        rebound_multi_well_context;
    FixedBhpMultiConnectionWellSourceAudit
        multi_well_source_audit{};
    bool multi_well_rebound{};
};

fdp::PostSnesPhaseTransitionProposal3D
controller_two_to_three_proposal() {
    flow::PhaseSetTransitionCandidate
        candidate;
    candidate.status =
        flow::
            PhaseSetTransitionCandidateStatus::
                target_resolved;
    candidate.component_ids =
        {"A", "B", "C"};
    candidate.pressure_pa = 20.0;
    candidate.temperature_k = 9.0;
    candidate.evidence_profile =
        "fixture/post-snes-scan/v1";
    candidate.source_phase_count = 2U;
    candidate.target_phase_count = 3U;
    candidate.trigger =
        flow::
            PhaseSetTransitionTrigger::
                stability_witness;
    candidate.target_phases = {
        {0.3, {0.3, 0.3, 0.4}, 2.0},
        {0.3, {0.3, 0.3, 0.4}, 4.0},
        {0.4, {0.3, 0.3, 0.4}, 7.0}};
    return {
        mesh::GlobalEntityId{
            UINT64_C(30)},
        std::move(candidate)};
}

PetscErrorCode
controller_scan(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&
            system,
    Vec,
    const fdp::
        VariableCardinalityNaturalVariableSnesSolveReport3D&
            solve_report,
    void* raw_context,
    fdp::PostSnesPhaseTransitionScanStatus3D*
        scan_status,
    std::vector<
        fdp::PostSnesPhaseTransitionProposal3D>*
            output) {
    if (raw_context == nullptr ||
        scan_status == nullptr ||
        output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    *scan_status =
        fdp::PostSnesPhaseTransitionScanStatus3D::
            complete;
    output->clear();
    auto* fixture =
        static_cast<
            ControllerFixture*>(
                raw_context);

    if (static_cast<int>(
            solve_report.converged_reason) <=
        0) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }
    if (fixture->scan_indeterminate) {
        *scan_status =
            fdp::PostSnesPhaseTransitionScanStatus3D::
                indeterminate;
        return PETSC_SUCCESS;
    }

    const auto& record =
        system.numbering().cell(
            mesh::LocalIndex{2U});
    if (record.owner_rank !=
        system.numbering()
            .local_rank()) {
        return PETSC_SUCCESS;
    }

    if (record.phase_count == 2U) {
        if (!fixture
                 ->phase_disappearance_only) {
            output->push_back(
                controller_two_to_three_proposal());
        }
        return PETSC_SUCCESS;
    }

    if (record.phase_count == 3U &&
        (fixture
             ->oscillate_after_restart ||
         fixture
             ->phase_disappearance_only)) {
        flow::PhaseSetTransitionCandidate
            candidate;
        candidate.status =
            flow::
                PhaseSetTransitionCandidateStatus::
                    target_resolved;
        candidate.component_ids =
            {"A", "B", "C"};
        candidate.pressure_pa = 20.0;
        candidate.temperature_k = 9.0;
        candidate.evidence_profile =
            "fixture/post-snes-scan/v1";
        candidate.source_phase_count = 3U;
        candidate.target_phase_count = 2U;
        candidate.trigger =
            flow::
                PhaseSetTransitionTrigger::
                    phase_disappearance;
        candidate.target_phases = {
            {0.5, {0.3, 0.3, 0.4}, 4.0},
            {0.5, {0.3, 0.3, 0.4}, 4.0}};
        output->push_back(
            {
                record.cell_global,
                std::move(candidate)});
    }
    return PETSC_SUCCESS;
}

PetscErrorCode
controller_rebuild(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&,
    Vec,
    const fdp::
        VariableCardinalityNaturalVariableSnesSolveReport3D&,
    std::span<
        const fdp::
            PostSnesPhaseTransitionProposal3D>
        local_owned_proposals,
    std::span<
        const fdp::
            AcceptedPhaseTransitionSummary3D>
        accepted_global_batch,
    void* raw_context,
    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>*
                output) {
    if (raw_context == nullptr ||
        output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    auto* fixture =
        static_cast<
            ControllerFixture*>(
                raw_context);

    if (accepted_global_batch.size() !=
            1U ||
        accepted_global_batch[0]
                .cell_global !=
            mesh::GlobalEntityId{
                UINT64_C(30)}) {
        return PETSC_ERR_ARG_INCOMP;
    }
    const bool target_three_phase =
        accepted_global_batch[0]
                .source_phase_count ==
            2U &&
        accepted_global_batch[0]
                .target_phase_count ==
            3U;
    const bool target_two_phase =
        accepted_global_batch[0]
                .source_phase_count ==
            3U &&
        accepted_global_batch[0]
                .target_phase_count ==
            2U;
    if (!target_three_phase &&
        !target_two_phase) {
        return PETSC_ERR_ARG_INCOMP;
    }

    if (fixture->rank == 0) {
        if (local_owned_proposals.size() !=
                1U ||
            local_owned_proposals[0]
                    .cell_global !=
                mesh::GlobalEntityId{
                    UINT64_C(30)} ||
            local_owned_proposals[0]
                    .candidate
                    .source_phase_count !=
                accepted_global_batch[0]
                    .source_phase_count ||
            local_owned_proposals[0]
                    .candidate
                    .target_phase_count !=
                accepted_global_batch[0]
                    .target_phase_count) {
            return PETSC_ERR_ARG_INCOMP;
        }
    } else if (
        !local_owned_proposals.empty()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    ++fixture->rebuild_calls;
    auto cells =
        make_controller_cells(
            fixture->rank,
            target_three_phase,
            fixture->audit);
    auto faces =
        make_face_inputs(
            fixture->rank,
            false);

    fdp::
        MixedCardinalityPhysicalCellSourceEvaluatorBinding3D
        source_binding{};
    if (fixture->well_context != nullptr) {
        const auto target =
            std::find_if(
                cells.begin(),
                cells.end(),
                [&](const auto& cell) {
                    return cell.cell_global ==
                        fixture
                            ->well_context
                            ->target_cell_global();
                });
        if (target == cells.end()) {
            return PETSC_ERR_ARG_INCOMP;
        }
        try {
            fixture->rebound_well_context.emplace(
                wdp::
                    rebind_fixed_bhp_peaceman_well_source_context_3d(
                        *fixture->well_context,
                        target->cell_global,
                        target->target_active_phases));
        } catch (...) {
            return PETSC_ERR_ARG_INCOMP;
        }
        fixture->well_context =
            &*fixture->rebound_well_context;
        fixture->well_source_audit.context =
            fixture->well_context;
        fixture->well_rebound = true;
        source_binding = {
            &evaluate_audited_fixed_bhp_well_source,
            &fixture->well_source_audit};
    } else if (
        fixture->multi_well_context != nullptr) {
        const auto transitioned =
            accepted_global_batch.front()
                .cell_global;
        const auto target =
            std::find_if(
                cells.begin(),
                cells.end(),
                [&](const auto& cell) {
                    return cell.cell_global ==
                        transitioned;
                });
        if (target == cells.end()) {
            return PETSC_ERR_ARG_INCOMP;
        }
        try {
            fixture
                ->rebound_multi_well_context
                .emplace(
                    wdp::
                        rebind_fixed_bhp_multi_connection_well_source_context_3d(
                            *fixture
                                 ->multi_well_context,
                            transitioned,
                            target
                                ->target_active_phases));
        } catch (...) {
            return PETSC_ERR_ARG_INCOMP;
        }
        fixture->multi_well_context =
            &*fixture
                 ->rebound_multi_well_context;
        fixture
            ->multi_well_source_audit
            .context =
            fixture->multi_well_context;
        fixture->multi_well_rebound =
            true;
        source_binding = {
            &evaluate_audited_fixed_bhp_multi_connection_well_source,
            &fixture
                 ->multi_well_source_audit};
    }

    return fdp::
        rebuild_phase_transition_natural_variable_system_3d(
            PETSC_COMM_WORLD,
            *fixture->schedule,
            *fixture->partition,
            *fixture->bridge,
            *fixture->pattern,
            1.0,
            std::move(cells),
            std::move(faces),
            {
                {&evaluate_1p, fixture->audit},
                {&evaluate_2p, fixture->audit},
                {&evaluate_3p, fixture->audit}},
            source_binding,
            &fdp::
                evaluate_absent_phase_thermodynamic_provider_3d<
                    flow::
                        Pr76AbsentPhasePotentialExtensionProvider<
                            double>>,
            fixture->provider,
            output);
}

std::unique_ptr<
    fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
make_controller_initial_system(
    ControllerFixture* fixture) {
    auto cells =
        make_controller_cells(
            fixture->rank,
            false,
            fixture->audit);
    auto faces =
        make_face_inputs(
            fixture->rank,
            false);
    fdp::
        MixedCardinalityPhysicalCellSourceEvaluatorBinding3D
        source_binding{};
    if (fixture->well_context != nullptr) {
        fixture->well_source_audit.context =
            fixture->well_context;
        source_binding = {
            &evaluate_audited_fixed_bhp_well_source,
            &fixture->well_source_audit};
    } else if (
        fixture->multi_well_context != nullptr) {
        fixture
            ->multi_well_source_audit
            .context =
            fixture->multi_well_context;
        source_binding = {
            &evaluate_audited_fixed_bhp_multi_connection_well_source,
            &fixture
                 ->multi_well_source_audit};
    }

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        system;
    const PetscErrorCode error =
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_WORLD,
                *fixture->schedule,
                *fixture->partition,
                *fixture->bridge,
                *fixture->pattern,
                1.0,
                std::move(cells),
                std::move(faces),
                {
                    {&evaluate_1p, fixture->audit},
                    {&evaluate_2p, fixture->audit},
                    {&evaluate_3p, fixture->audit}},
                source_binding,
                &fdp::
                    evaluate_absent_phase_thermodynamic_provider_3d<
                        flow::
                            Pr76AbsentPhasePotentialExtensionProvider<
                                double>>,
                fixture->provider,
                &system);
    if (error != PETSC_SUCCESS ||
        system == nullptr) {
        throw std::runtime_error(
            "failed to build controller initial system");
    }
    return system;
}

void run_controller_case(
    ControllerFixture* fixture,
    std::size_t max_restarts,
    fdp::PostSnesPhaseTransitionOutcome3D
        expected_outcome,
    std::size_t expected_restarts) {
    fixture->rebuild_calls = 0U;
    auto initial =
        make_controller_initial_system(
            fixture);

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        final_system;
    Vec final_state = nullptr;
    std::optional<
        fdp::
            PostSnesPhaseTransitionControllerReport3D>
        controller_report;
    const PetscErrorCode error =
        fdp::
            solve_nonlinear_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                std::move(initial),
                {
                    &controller_scan,
                    fixture,
                    &controller_rebuild,
                    fixture},
                {max_restarts},
                &final_system,
                &final_state,
                &controller_report);

    require_collective(
        error == PETSC_SUCCESS &&
            final_system != nullptr &&
            final_state != nullptr &&
            controller_report.has_value() &&
            controller_report->outcome ==
                expected_outcome &&
            controller_report
                    ->transition_restarts ==
                expected_restarts,
        "post-SNES phase-transition controller outcome mismatch");

    const bool expected_accepted =
        expected_outcome ==
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                stable_phase_set;
    require_collective(
        controller_report
                ->timestep_accepted() ==
            expected_accepted,
        "post-SNES controller timestep acceptance semantics mismatch");

    if (expected_outcome ==
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                stable_phase_set) {
        require_collective(
            controller_report
                    ->generations.size() ==
                2U &&
                final_system
                        ->numbering()
                        .cell(
                            mesh::LocalIndex{2U})
                        .phase_count ==
                    3U &&
                final_system
                        ->numbering()
                        .petsc_global_scalar_count() ==
                    45,
            "stable controller path did not finish on rebuilt 3P phase set");
    } else if (
        expected_outcome ==
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                phase_set_cycle_detected) {
        require_collective(
            controller_report
                    ->generations.size() ==
                2U &&
                controller_report
                    ->generations[1]
                    .accepted_transition_batch
                    .size() ==
                1U &&
                fixture->rebuild_calls ==
                    2U,
            "physical phase-set cycle was not detected after the reverse rebuild");
    } else if (
        expected_outcome ==
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                phase_set_scan_indeterminate) {
        require_collective(
            controller_report
                    ->generations.size() ==
                1U &&
                controller_report
                    ->generations.front()
                    .phase_transition_scan_status ==
                    fdp::
                        PostSnesPhaseTransitionScanStatus3D::
                            indeterminate &&
                fixture->rebuild_calls ==
                    0U,
            "indeterminate phase-set scan was incorrectly accepted or rebuilt");
    } else if (
        expected_outcome ==
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                transition_restart_budget_exhausted) {
        require_collective(
            controller_report
                    ->generations.size() ==
                1U &&
                fixture->rebuild_calls ==
                    0U,
            "transition budget exhausted path rebuilt the system");
    }

    require_collective(
        VecDestroy(&final_state) ==
            PETSC_SUCCESS,
        "post-SNES controller final state cleanup failed");
}

[[nodiscard]] bool
local_solution_matches_target(
    const fdp::
        VariableCardinalityNaturalVariableSnesSolveReport3D&
            report) {
    for (const auto& entry :
         report.locally_owned_solution) {
        const auto target =
            target_state(
                entry.cell_global.value());
        if (entry.natural_variable_slot >=
                target.size() ||
            std::abs(
                entry.value -
                target[
                    entry.natural_variable_slot]) >
                1.0e-10 *
                    std::max(
                        1.0,
                        std::abs(
                            target[
                                entry.natural_variable_slot]))) {
            return false;
        }
    }
    return true;
}


struct FrozenWellTimestepControlAudit {
    mesh::LocalIndex target_cell{
        mesh::LocalIndex::value_type{0}};
    mesh::GlobalEntityId target_cell_global{
        mesh::GlobalEntityId::value_type{0}};
    std::size_t expected_phase_count{};
    std::size_t expected_scalar_count{};
    std::size_t scans{};
    std::size_t rebuild_calls{};
};

PetscErrorCode
frozen_well_timestep_scan(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&
            system,
    Vec,
    const fdp::
        VariableCardinalityNaturalVariableSnesSolveReport3D&
            solve_report,
    void* raw_context,
    fdp::PostSnesPhaseTransitionScanStatus3D*
        scan_status,
    std::vector<
        fdp::PostSnesPhaseTransitionProposal3D>*
            output) {
    if (raw_context == nullptr ||
        scan_status == nullptr ||
        output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* audit =
        static_cast<
            FrozenWellTimestepControlAudit*>(
                raw_context);
    *scan_status =
        fdp::PostSnesPhaseTransitionScanStatus3D::
            complete;
    output->clear();

    if (static_cast<int>(
            solve_report.converged_reason) <=
        0) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }

    const auto& well_cell =
        system.numbering().cell(
            audit->target_cell);
    if (well_cell.cell_global !=
            audit->target_cell_global ||
        well_cell.phase_count !=
            audit->expected_phase_count ||
        well_cell.scalar_count !=
            audit->expected_scalar_count) {
        return PETSC_ERR_ARG_INCOMP;
    }

    ++audit->scans;
    return PETSC_SUCCESS;
}

struct FrozenMultiWellTimestepControlAudit {
    std::size_t expected_cell30_phase_count{2U};
    std::size_t expected_cell30_scalar_count{7U};
    PetscInt expected_global_scalar_count{42};
    std::size_t scans{};
    std::size_t rebuild_calls{};
};

PetscErrorCode
frozen_multi_well_timestep_scan(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&
            system,
    Vec,
    const fdp::
        VariableCardinalityNaturalVariableSnesSolveReport3D&
            solve_report,
    void* raw_context,
    fdp::PostSnesPhaseTransitionScanStatus3D*
        scan_status,
    std::vector<
        fdp::PostSnesPhaseTransitionProposal3D>*
            output) {
    if (raw_context == nullptr ||
        scan_status == nullptr ||
        output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* audit =
        static_cast<
            FrozenMultiWellTimestepControlAudit*>(
                raw_context);
    *scan_status =
        fdp::PostSnesPhaseTransitionScanStatus3D::
            complete;
    output->clear();

    if (static_cast<int>(
            solve_report.converged_reason) <=
        0) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }

    const auto& cell30 =
        system.numbering().cell(
            mesh::LocalIndex{2U});
    const auto& cell60 =
        system.numbering().cell(
            mesh::LocalIndex{5U});
    if (system.numbering()
                .petsc_global_scalar_count() !=
            audit->expected_global_scalar_count ||
        cell30.cell_global !=
            mesh::GlobalEntityId{
                UINT64_C(30)} ||
        cell30.phase_count !=
            audit->expected_cell30_phase_count ||
        cell30.scalar_count !=
            audit->expected_cell30_scalar_count ||
        cell60.cell_global !=
            mesh::GlobalEntityId{
                UINT64_C(60)} ||
        cell60.phase_count != 3U ||
        cell60.scalar_count != 10U) {
        return PETSC_ERR_ARG_INCOMP;
    }

    ++audit->scans;
    return PETSC_SUCCESS;
}

PetscErrorCode
unexpected_frozen_multi_well_timestep_rebuild(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&,
    Vec,
    const fdp::
        VariableCardinalityNaturalVariableSnesSolveReport3D&,
    std::span<
        const fdp::
            PostSnesPhaseTransitionProposal3D>,
    std::span<
        const fdp::
            AcceptedPhaseTransitionSummary3D>,
    void* raw_context,
    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>*
                output) {
    if (raw_context == nullptr ||
        output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    auto* audit =
        static_cast<
            FrozenMultiWellTimestepControlAudit*>(
                raw_context);
    ++audit->rebuild_calls;
    return PETSC_ERR_PLIB;
}

PetscErrorCode
unexpected_frozen_well_timestep_rebuild(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&,
    Vec,
    const fdp::
        VariableCardinalityNaturalVariableSnesSolveReport3D&,
    std::span<
        const fdp::
            PostSnesPhaseTransitionProposal3D>,
    std::span<
        const fdp::
            AcceptedPhaseTransitionSummary3D>,
    void* raw_context,
    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>*
                output) {
    if (raw_context == nullptr ||
        output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    auto* audit =
        static_cast<
            FrozenWellTimestepControlAudit*>(
                raw_context);
    ++audit->rebuild_calls;
    return PETSC_ERR_PLIB;
}

std::unique_ptr<
    fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
make_frozen_fixed_bhp_timestep_system(
    int rank,
    const dp::
        ParallelOwnedConnectionSchedule3D&
            schedule,
    const mesh::PartitionSnapshot&
        partition,
    const dp::
        PetscMpiAijSymbolicPreallocation3D&
            bridge,
    const dp::
        OwnedCellStructuralColumnPatternSnapshot3D&
            pattern,
    DispatchAudit* dispatch_audit,
    FixedBhpWellSourceAudit* source_audit,
    flow::
        Pr76AbsentPhasePotentialExtensionProvider<
            double>* provider) {
    if (dispatch_audit == nullptr ||
        source_audit == nullptr ||
        source_audit->context == nullptr ||
        provider == nullptr) {
        throw std::invalid_argument(
            "invalid frozen fixed-BHP timestep fixture context");
    }

    auto cells =
        make_controller_cells(
            rank,
            false,
            dispatch_audit);
    auto faces =
        make_face_inputs(
            rank,
            false);

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        system;
    const PetscErrorCode error =
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_WORLD,
                schedule,
                partition,
                bridge,
                pattern,
                1.0,
                std::move(cells),
                std::move(faces),
                {
                    {&evaluate_1p, dispatch_audit},
                    {&evaluate_2p, dispatch_audit},
                    {&evaluate_3p, dispatch_audit}},
                {
                    &evaluate_audited_fixed_bhp_well_source,
                    source_audit},
                &fdp::
                    evaluate_absent_phase_thermodynamic_provider_3d<
                        flow::
                            Pr76AbsentPhasePotentialExtensionProvider<
                                double>>,
                provider,
                &system);
    if (error != PETSC_SUCCESS ||
        system == nullptr) {
        throw std::runtime_error(
            "failed to build frozen fixed-BHP physical timestep system");
    }
    return system;
}

std::unique_ptr<
    fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
make_frozen_multi_connection_fixed_bhp_timestep_system(
    int rank,
    const dp::
        ParallelOwnedConnectionSchedule3D&
            schedule,
    const mesh::PartitionSnapshot&
        partition,
    const dp::
        PetscMpiAijSymbolicPreallocation3D&
            bridge,
    const dp::
        OwnedCellStructuralColumnPatternSnapshot3D&
            pattern,
    DispatchAudit* dispatch_audit,
    FixedBhpMultiConnectionWellSourceAudit*
        source_audit,
    flow::
        Pr76AbsentPhasePotentialExtensionProvider<
            double>* provider) {
    if (dispatch_audit == nullptr ||
        source_audit == nullptr ||
        source_audit->context == nullptr ||
        provider == nullptr) {
        throw std::invalid_argument(
            "invalid multi-connection fixed-BHP timestep fixture context");
    }

    auto cells =
        make_controller_cells(
            rank,
            false,
            dispatch_audit);
    auto faces =
        make_face_inputs(
            rank,
            false);

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        system;
    const PetscErrorCode error =
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_WORLD,
                schedule,
                partition,
                bridge,
                pattern,
                1.0,
                std::move(cells),
                std::move(faces),
                {
                    {&evaluate_1p, dispatch_audit},
                    {&evaluate_2p, dispatch_audit},
                    {&evaluate_3p, dispatch_audit}},
                {
                    &evaluate_audited_fixed_bhp_multi_connection_well_source,
                    source_audit},
                &fdp::
                    evaluate_absent_phase_thermodynamic_provider_3d<
                        flow::
                            Pr76AbsentPhasePotentialExtensionProvider<
                                double>>,
                provider,
                &system);
    if (error != PETSC_SUCCESS ||
        system == nullptr) {
        throw std::runtime_error(
            "failed to build multi-connection fixed-BHP physical timestep system");
    }
    return system;
}

std::unique_ptr<
    fdp::PhaseTransitionRebuiltNaturalVariableSystem3D>
make_fixed_total_molar_rate_reservoir_system(
    int rank,
    const dp::
        ParallelOwnedConnectionSchedule3D&
            schedule,
    const mesh::PartitionSnapshot&
        partition,
    const dp::
        PetscMpiAijSymbolicPreallocation3D&
            bridge,
    const dp::
        OwnedCellStructuralColumnPatternSnapshot3D&
            pattern,
    DispatchAudit* dispatch_audit,
    wdp::
        FixedTotalMolarRateWellSourceEvaluatorContext3D*
            source_context,
    flow::
        Pr76AbsentPhasePotentialExtensionProvider<
            double>* provider,
    bool cell30_is_three_phase = false) {
    if (dispatch_audit == nullptr ||
        source_context == nullptr ||
        provider == nullptr) {
        throw std::invalid_argument(
            "invalid fixed-total-molar-rate reservoir fixture context");
    }

    auto cells =
        make_controller_cells(
            rank,
            cell30_is_three_phase,
            dispatch_audit);
    auto faces =
        make_face_inputs(
            rank,
            false);

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        system;
    const PetscErrorCode error =
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_WORLD,
                schedule,
                partition,
                bridge,
                pattern,
                1.0,
                std::move(cells),
                std::move(faces),
                {
                    {&evaluate_1p, dispatch_audit},
                    {&evaluate_2p, dispatch_audit},
                    {&evaluate_3p, dispatch_audit}},
                wdp::
                    fixed_total_molar_rate_well_source_binding_3d(
                        source_context),
                &fdp::
                    evaluate_absent_phase_thermodynamic_provider_3d<
                        flow::
                            Pr76AbsentPhasePotentialExtensionProvider<
                                double>>,
                provider,
                &system);
    if (error != PETSC_SUCCESS ||
        system == nullptr) {
        throw std::runtime_error(
            "failed to build fixed-total-molar-rate reservoir system");
    }
    return system;
}

std::array<double, 4>
owned_conserved_totals(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&
            system,
    Vec state) {
    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        current;
    std::vector<double>
        porosities;
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    const PetscErrorCode error =
        system.evaluate_local_cells_for_phase_transition(
            state,
            &current,
            &porosities,
            &status);
    if (error != PETSC_SUCCESS ||
        status !=
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success ||
        current.size() != porosities.size()) {
        throw std::runtime_error(
            "failed to evaluate conserved totals for frozen well timestep");
    }

    std::array<double, 4>
        totals{};
    for (std::size_t local = 0U;
         local < current.size();
         ++local) {
        const auto& record =
            system.numbering().cell(
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            local)});
        if (record.owner_rank !=
            system.numbering().local_rank()) {
            continue;
        }
        if (!current[local].has_value()) {
            throw std::runtime_error(
                "owned frozen well timestep cell has no current state");
        }

        const double bulk_volume_m3 =
            2.0 +
            static_cast<double>(local);
        std::visit(
            [&](const auto& typed) {
                using Typed =
                    std::remove_cvref_t<
                        decltype(typed)>;

                if constexpr (
                    std::is_same_v<
                        Typed,
                        fdp::
                            SinglePhaseCurrentCellLinearization3D>) {
                    const auto component =
                        flow::
                            build_single_phase_component_accumulation(
                                typed.state,
                                porosities[local]);
                    const auto energy =
                        flow::
                            build_single_phase_energy_accumulation_snapshot(
                                typed.state,
                                porosities[local],
                                typed.transport,
                                typed.caloric,
                                typed.rock);
                    for (std::size_t i = 0U;
                         i < 3U;
                         ++i) {
                        totals[i] +=
                            component
                                .component_accumulation_mol_per_bulk_m3[i] *
                            bulk_volume_m3;
                    }
                    totals[3] +=
                        energy
                            .total_internal_energy_j_per_bulk_m3 *
                        bulk_volume_m3;
                } else if constexpr (
                    std::is_same_v<
                        Typed,
                        fdp::
                            TwoPhaseCurrentCellLinearization3D>) {
                    const auto component =
                        flow::
                            build_two_phase_component_accumulation(
                                typed.state,
                                porosities[local]);
                    const auto energy =
                        flow::
                            build_two_phase_energy_accumulation_snapshot(
                                typed.state,
                                porosities[local],
                                typed.transport,
                                typed.caloric,
                                typed.rock);
                    for (std::size_t i = 0U;
                         i < 3U;
                         ++i) {
                        totals[i] +=
                            component
                                .component_accumulation_mol_per_bulk_m3[i] *
                            bulk_volume_m3;
                    }
                    totals[3] +=
                        energy
                            .total_internal_energy_j_per_bulk_m3 *
                        bulk_volume_m3;
                } else {
                    const auto component =
                        flow::
                            build_pore_volume_component_accumulation(
                                typed.state,
                                porosities[local]);
                    const auto energy =
                        flow::
                            build_pore_volume_energy_accumulation_snapshot(
                                typed.state,
                                porosities[local],
                                typed.transport,
                                typed.caloric,
                                typed.rock);
                    for (std::size_t i = 0U;
                         i < 3U;
                         ++i) {
                        totals[i] +=
                            component
                                .component_accumulation_mol_per_bulk_m3[i] *
                            bulk_volume_m3;
                    }
                    totals[3] +=
                        energy
                            .total_internal_energy_j_per_bulk_m3 *
                        bulk_volume_m3;
                }
            },
            *current[local]);
    }
    return totals;
}

std::array<double, 4>
global_fixed_bhp_well_production_rate(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&
            system,
    Vec state,
    const wdp::
        FixedBhpMultiConnectionWellSourceEvaluatorContext3D&
            well_context,
    double bottom_hole_pressure_pa,
    std::uint64_t*
        authoritative_connection_count) {
    if (state == nullptr ||
        authoritative_connection_count ==
            nullptr ||
        !std::isfinite(
            bottom_hole_pressure_pa) ||
        !(bottom_hole_pressure_pa > 0.0)) {
        throw std::invalid_argument(
            "invalid fixed-BHP well-rate probe input");
    }

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        current;
    std::vector<double>
        porosities;
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    if (system
            .evaluate_local_cells_for_phase_transition(
                state,
                &current,
                &porosities,
                &status) !=
            PETSC_SUCCESS ||
        status !=
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success ||
        current.size() !=
            porosities.size()) {
        throw std::runtime_error(
            "failed to evaluate fixed-BHP well-rate probe state");
    }

    std::array<double, 4>
        local_rate{};
    std::uint64_t
        local_count = 0U;
    for (std::size_t local = 0U;
         local < current.size();
         ++local) {
        const auto& record =
            system.numbering().cell(
                mesh::LocalIndex{
                    static_cast<
                        mesh::LocalIndex::value_type>(
                            local)});
        if (record.owner_rank !=
            system.numbering().local_rank()) {
            continue;
        }
        const auto* connection =
            well_context.find_connection(
                record.cell_global);
        if (connection == nullptr) {
            continue;
        }
        if (!current[local].has_value()) {
            throw std::runtime_error(
                "authoritative fixed-BHP connection has no current cell state");
        }

        const auto source =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    connection
                        ->with_bottom_hole_pressure(
                            bottom_hole_pressure_pa),
                    *current[local]);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            local_rate[component] -=
                source
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component];
        }
        local_rate[3] -=
            source
                .cell_source
                .energy_rate_w;
        ++local_count;
    }

    std::array<double, 4>
        global_rate{};
    std::uint64_t
        global_count = 0U;
    if (MPI_Allreduce(
            local_rate.data(),
            global_rate.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) !=
            MPI_SUCCESS ||
        MPI_Allreduce(
            &local_count,
            &global_count,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            PETSC_COMM_WORLD) !=
            MPI_SUCCESS) {
        throw std::runtime_error(
            "failed to reduce fixed-BHP well-rate probe");
    }

    *authoritative_connection_count =
        global_count;
    return global_rate;
}


void run_rebound_fixed_bhp_physical_timestep(
    int rank,
    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>*
                system,
    Vec* stable_state,
    ControllerFixture* fixture,
    std::size_t expected_phase_count,
    std::size_t expected_scalar_count) {
    if (system == nullptr ||
        *system == nullptr ||
        stable_state == nullptr ||
        *stable_state == nullptr ||
        fixture == nullptr ||
        fixture->well_context == nullptr ||
        expected_phase_count == 0U ||
        expected_phase_count > 3U) {
        throw std::invalid_argument(
            "invalid rebound fixed-BHP physical timestep fixture");
    }

    const auto& target_record =
        (*system)
            ->numbering()
            .cell(
                mesh::LocalIndex{2U});
    require_collective(
        target_record.cell_global ==
                mesh::GlobalEntityId{
                    UINT64_C(30)} &&
            target_record.phase_count ==
                expected_phase_count &&
            target_record.scalar_count ==
                expected_scalar_count,
        "rebound fixed-BHP target chart does not match expected cardinality");

    const bool local_target_owner =
        target_record.owner_rank ==
        (*system)
            ->numbering()
            .local_rank();

    require_collective(
        (*system)
                ->rebase_accepted_timestep(
                    *stable_state,
                    1.0) ==
            PETSC_SUCCESS,
        "failed to rebase stable rebound fixed-BHP state");
    require_collective(
        VecDestroy(
            stable_state) ==
            PETSC_SUCCESS,
        "failed to destroy stable rebound fixed-BHP controller state");

    bool rebased_history_matches = false;
    require_collective(
        (*system)
                ->accepted_history_matches_state(
                    (*system)
                        ->initial_state(),
                    &rebased_history_matches) ==
            PETSC_SUCCESS &&
            rebased_history_matches,
        "rebound fixed-BHP accepted history did not match stable state");

    const auto local_previous_total =
        owned_conserved_totals(
            **system,
            (*system)
                ->initial_state());
    std::array<double, 4>
        global_previous_total{};
    require_collective(
        MPI_Allreduce(
            local_previous_total.data(),
            global_previous_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce rebound fixed-BHP previous conserved totals");

    Vec previous_state = nullptr;
    require_collective(
        VecDuplicate(
            (*system)
                ->initial_state(),
            &previous_state) ==
                PETSC_SUCCESS &&
            VecCopy(
                (*system)
                    ->initial_state(),
                previous_state) ==
                PETSC_SUCCESS,
        "failed to preserve rebound fixed-BHP previous state");

    fixture->well_source_audit.evaluator_calls =
        0U;
    fixture->well_source_audit.target_calls =
        0U;

    fdp::PhysicalTimestepDriverOptions3D
        timestep_options;
    timestep_options.adaptive
        .minimum_timestep_seconds =
        0.125;
    timestep_options.adaptive
        .maximum_timestep_seconds =
        1.0;
    timestep_options.adaptive
        .cutback_factor =
        0.5;
    timestep_options.adaptive
        .growth_factor =
        2.0;
    timestep_options.adaptive
        .maximum_retries =
        4U;
    timestep_options.adaptive
        .growth_nonlinear_iteration_limit =
        20;
    timestep_options.adaptive
        .growth_line_search_direction_change_limit =
        4;
    timestep_options.adaptive
        .growth_transition_restart_limit =
        0U;
    timestep_options.phase_transition
        .max_transition_restarts =
        0U;

    fdp::AcceptedPhysicalTimeClock3D
        timestep_clock{
            0.0,
            1.0};
    FrozenWellTimestepControlAudit
        frozen_scan{
            mesh::LocalIndex{2U},
            mesh::GlobalEntityId{
                UINT64_C(30)},
            expected_phase_count,
            expected_scalar_count,
            0U,
            0U};
    std::optional<
        fdp::PhysicalTimestepDriverReport3D>
        timestep_report;
    const PetscErrorCode error =
        fdp::advance_one_physical_timestep_3d(
            PETSC_COMM_WORLD,
            system,
            0U,
            {
                &frozen_well_timestep_scan,
                &frozen_scan,
                &unexpected_frozen_well_timestep_rebuild,
                &frozen_scan},
            timestep_options,
            &timestep_clock,
            &timestep_report);

    require_collective(
        error == PETSC_SUCCESS &&
            timestep_report.has_value() &&
            timestep_report->accepted() &&
            timestep_report
                    ->adaptive
                    .accepted_timestep_seconds
                    .has_value() &&
            timestep_report
                    ->accepted_record
                    .has_value() &&
            timestep_report
                    ->accepted_record
                    ->phase_transition_restarts ==
                0U &&
            timestep_clock
                    .accepted_step_count() ==
                1U &&
            frozen_scan.scans > 0U &&
            frozen_scan.rebuild_calls == 0U &&
            (*system)
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                expected_phase_count &&
            (local_target_owner
                 ? fixture
                           ->well_source_audit
                           .target_calls >
                       0U
                 : fixture
                           ->well_source_audit
                           .target_calls ==
                       0U),
        "rebound fixed-BHP completion did not complete an owner-only accepted physical timestep");

    const double accepted_dt =
        *timestep_report
             ->adaptive
             .accepted_timestep_seconds;

    Vec state_delta = nullptr;
    PetscReal state_delta_norm = 0.0;
    require_collective(
        VecDuplicate(
            (*system)
                ->initial_state(),
            &state_delta) ==
                PETSC_SUCCESS &&
            VecCopy(
                (*system)
                    ->initial_state(),
                state_delta) ==
                PETSC_SUCCESS &&
            VecAXPY(
                state_delta,
                PetscScalar{-1.0},
                previous_state) ==
                PETSC_SUCCESS &&
            VecNorm(
                state_delta,
                NORM_2,
                &state_delta_norm) ==
                PETSC_SUCCESS &&
            static_cast<double>(
                state_delta_norm) >
                1.0e-10,
        "rebound fixed-BHP physical timestep did not change reservoir state");

    const auto local_final_total =
        owned_conserved_totals(
            **system,
            (*system)
                ->initial_state());
    std::array<double, 4>
        global_final_total{};
    require_collective(
        MPI_Allreduce(
            local_final_total.data(),
            global_final_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce rebound fixed-BHP final conserved totals");

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        final_current;
    std::vector<double>
        final_porosity;
    fdp::NaturalVariableSnesEvaluationStatus3D
        final_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    require_collective(
        (*system)
                ->evaluate_local_cells_for_phase_transition(
                    (*system)
                        ->initial_state(),
                    &final_current,
                    &final_porosity,
                    &final_status) ==
            PETSC_SUCCESS &&
            final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            final_current.size() >
                2U &&
            final_current[2U]
                .has_value(),
        "failed to evaluate accepted rebound fixed-BHP well state");

    std::array<double, 4>
        local_well_production_rate{};
    bool local_final_source_shape_ok =
        true;
    if (local_target_owner) {
        const auto final_well =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    *fixture->well_context,
                    *final_current[2U]);
        local_final_source_shape_ok =
            final_well
                    .cell_source
                    .input_count ==
                expected_scalar_count;
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            local_well_production_rate[
                component] =
                -final_well
                     .cell_source
                     .component_molar_rate_mol_per_s[
                         component];
        }
        local_well_production_rate[3] =
            -final_well
                 .cell_source
                 .energy_rate_w;
    }
    require_collective(
        local_final_source_shape_ok,
        "accepted rebound fixed-BHP source lost expected Jacobian cardinality");

    std::array<double, 4>
        global_well_production_rate{};
    require_collective(
        MPI_Allreduce(
            local_well_production_rate.data(),
            global_well_production_rate.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce accepted rebound fixed-BHP well rates");

    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        require_collective(
            std::isfinite(
                global_well_production_rate[
                    component]) &&
                std::abs(
                    global_well_production_rate[
                        component]) >
                    0.0,
            "accepted rebound fixed-BHP component rate vanished or became non-finite");
        near_collective(
            global_final_total[
                component],
            global_previous_total[
                component] -
                accepted_dt *
                    global_well_production_rate[
                        component],
            2.0e-7,
            2.0e-8);
    }
    require_collective(
        std::isfinite(
            global_well_production_rate[3]) &&
            std::abs(
                global_well_production_rate[3]) >
                0.0,
        "accepted rebound fixed-BHP energy rate vanished or became non-finite");
    near_collective(
        global_final_total[3],
        global_previous_total[3] -
            accepted_dt *
                global_well_production_rate[3],
        2.0e-7,
        2.0e-7);

    bool accepted_history_matches = false;
    require_collective(
        (*system)
                ->accepted_history_matches_state(
                    (*system)
                        ->initial_state(),
                    &accepted_history_matches) ==
            PETSC_SUCCESS &&
            accepted_history_matches,
        "rebound fixed-BHP accepted timestep did not rebase history");

    require_collective(
        VecDestroy(
            &state_delta) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &previous_state) ==
                PETSC_SUCCESS,
        "rebound fixed-BHP physical timestep cleanup failed");

    (void)rank;
}

void run_rebound_multi_connection_fixed_bhp_physical_timestep(
    int rank,
    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>*
                system,
    Vec* stable_state,
    ControllerFixture* fixture,
    std::size_t expected_cell30_phase_count,
    std::size_t expected_cell30_scalar_count,
    PetscInt expected_global_scalar_count) {
    if (system == nullptr ||
        *system == nullptr ||
        stable_state == nullptr ||
        *stable_state == nullptr ||
        fixture == nullptr ||
        fixture->multi_well_context == nullptr) {
        throw std::invalid_argument(
            "invalid rebound multi-connection fixed-BHP physical timestep fixture");
    }

    const auto& record30 =
        (*system)->numbering().cell(
            mesh::LocalIndex{2U});
    const auto& record60 =
        (*system)->numbering().cell(
            mesh::LocalIndex{5U});
    require_collective(
        (*system)->numbering()
                .petsc_global_scalar_count() ==
            expected_global_scalar_count &&
            record30.cell_global ==
                mesh::GlobalEntityId{
                    UINT64_C(30)} &&
            record30.phase_count ==
                expected_cell30_phase_count &&
            record30.scalar_count ==
                expected_cell30_scalar_count &&
            record30.owner_rank.value() == 0U &&
            record60.cell_global ==
                mesh::GlobalEntityId{
                    UINT64_C(60)} &&
            record60.phase_count == 3U &&
            record60.scalar_count == 10U &&
            record60.owner_rank.value() == 1U,
        "rebound multi-connection fixed-BHP target charts or ownership changed unexpectedly");

    require_collective(
        (*system)
                ->rebase_accepted_timestep(
                    *stable_state,
                    1.0) ==
            PETSC_SUCCESS &&
            VecDestroy(
                stable_state) ==
            PETSC_SUCCESS,
        "failed to rebase stable multi-connection fixed-BHP state");

    bool rebased_history_matches = false;
    require_collective(
        (*system)
                ->accepted_history_matches_state(
                    (*system)
                        ->initial_state(),
                    &rebased_history_matches) ==
            PETSC_SUCCESS &&
            rebased_history_matches,
        "rebound multi-connection fixed-BHP history did not match stable state");

    const auto local_previous_total =
        owned_conserved_totals(
            **system,
            (*system)
                ->initial_state());
    std::array<double, 4>
        global_previous_total{};
    require_collective(
        MPI_Allreduce(
            local_previous_total.data(),
            global_previous_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce rebound multi-connection previous conserved totals");

    Vec previous_state = nullptr;
    require_collective(
        VecDuplicate(
            (*system)
                ->initial_state(),
            &previous_state) ==
                PETSC_SUCCESS &&
            VecCopy(
                (*system)
                    ->initial_state(),
                previous_state) ==
                PETSC_SUCCESS,
        "failed to preserve rebound multi-connection previous state");

    fixture
        ->multi_well_source_audit
        .evaluator_calls = 0U;
    fixture
        ->multi_well_source_audit
        .cell30_calls = 0U;
    fixture
        ->multi_well_source_audit
        .cell60_calls = 0U;

    fdp::PhysicalTimestepDriverOptions3D
        timestep_options;
    timestep_options.adaptive
        .minimum_timestep_seconds =
        0.125;
    timestep_options.adaptive
        .maximum_timestep_seconds =
        1.0;
    timestep_options.adaptive
        .cutback_factor =
        0.5;
    timestep_options.adaptive
        .growth_factor =
        2.0;
    timestep_options.adaptive
        .maximum_retries =
        4U;
    timestep_options.adaptive
        .growth_nonlinear_iteration_limit =
        20;
    timestep_options.adaptive
        .growth_line_search_direction_change_limit =
        4;
    timestep_options.adaptive
        .growth_transition_restart_limit =
        0U;
    timestep_options.phase_transition
        .max_transition_restarts =
        0U;

    FrozenMultiWellTimestepControlAudit
        control_audit;
    control_audit
        .expected_cell30_phase_count =
        expected_cell30_phase_count;
    control_audit
        .expected_cell30_scalar_count =
        expected_cell30_scalar_count;
    control_audit
        .expected_global_scalar_count =
        expected_global_scalar_count;

    fdp::AcceptedPhysicalTimeClock3D
        timestep_clock{
            0.0,
            1.0};
    std::optional<
        fdp::PhysicalTimestepDriverReport3D>
        timestep_report;
    const PetscErrorCode error =
        fdp::advance_one_physical_timestep_3d(
            PETSC_COMM_WORLD,
            system,
            0U,
            {
                &frozen_multi_well_timestep_scan,
                &control_audit,
                &unexpected_frozen_multi_well_timestep_rebuild,
                &control_audit},
            timestep_options,
            &timestep_clock,
            &timestep_report);

    require_collective(
        error == PETSC_SUCCESS &&
            timestep_report.has_value() &&
            timestep_report->accepted() &&
            timestep_report
                    ->adaptive
                    .accepted_timestep_seconds
                    .has_value() &&
            timestep_clock
                    .accepted_step_count() ==
                1U &&
            control_audit.scans > 0U &&
            control_audit.rebuild_calls == 0U &&
            fixture
                    ->multi_well_source_audit
                    .evaluator_calls >
                0U &&
            (rank == 0
                 ? fixture
                           ->multi_well_source_audit
                           .cell30_calls >
                       0U &&
                       fixture
                               ->multi_well_source_audit
                               .cell60_calls ==
                           0U
                 : fixture
                           ->multi_well_source_audit
                           .cell30_calls ==
                       0U &&
                       fixture
                               ->multi_well_source_audit
                               .cell60_calls >
                           0U),
        "rebound multi-connection sources were not assembled exactly on authoritative owners");

    const double accepted_dt =
        *timestep_report
             ->adaptive
             .accepted_timestep_seconds;

    Vec state_delta = nullptr;
    PetscReal state_delta_norm = 0.0;
    require_collective(
        VecDuplicate(
            (*system)
                ->initial_state(),
            &state_delta) ==
                PETSC_SUCCESS &&
            VecCopy(
                (*system)
                    ->initial_state(),
                state_delta) ==
                PETSC_SUCCESS &&
            VecAXPY(
                state_delta,
                PetscScalar{-1.0},
                previous_state) ==
                PETSC_SUCCESS &&
            VecNorm(
                state_delta,
                NORM_2,
                &state_delta_norm) ==
                PETSC_SUCCESS &&
            static_cast<double>(
                state_delta_norm) >
                1.0e-10,
        "rebound multi-connection accepted timestep did not change reservoir state");

    const auto local_final_total =
        owned_conserved_totals(
            **system,
            (*system)
                ->initial_state());
    std::array<double, 4>
        global_final_total{};
    require_collective(
        MPI_Allreduce(
            local_final_total.data(),
            global_final_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce rebound multi-connection final conserved totals");

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        final_current;
    std::vector<double>
        final_porosity;
    fdp::NaturalVariableSnesEvaluationStatus3D
        final_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    require_collective(
        (*system)
                ->evaluate_local_cells_for_phase_transition(
                    (*system)
                        ->initial_state(),
                    &final_current,
                    &final_porosity,
                    &final_status) ==
            PETSC_SUCCESS &&
            final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            final_current.size() > 5U &&
            final_current[2U].has_value() &&
            final_current[5U].has_value(),
        "failed to evaluate accepted rebound multi-connection states");

    std::vector<
        wd::
            FixedBhpConnectionCellSourceLinearization3D>
        local_authoritative_sources;
    std::array<double, 4>
        local_connection30_rate{};
    std::array<double, 4>
        local_connection60_rate{};
    bool local_connection_shapes_ok =
        true;

    if (record30.owner_rank ==
        (*system)->numbering().local_rank()) {
        const auto* connection =
            fixture->multi_well_context
                ->find_connection(
                    record30.cell_global);
        if (connection == nullptr) {
            local_connection_shapes_ok =
                false;
        } else {
            const auto source30 =
                wdp::
                    build_fixed_bhp_peaceman_well_source_3d(
                        *connection,
                        *final_current[2U]);
            local_connection_shapes_ok =
                local_connection_shapes_ok &&
                source30.cell_source.input_count ==
                    expected_cell30_scalar_count;
            for (std::size_t component = 0U;
                 component < 3U;
                 ++component) {
                local_connection30_rate[
                    component] =
                    -source30
                         .cell_source
                         .component_molar_rate_mol_per_s[
                             component];
            }
            local_connection30_rate[3] =
                -source30
                     .cell_source
                     .energy_rate_w;
            local_authoritative_sources
                .push_back(source30);
        }
    }

    if (record60.owner_rank ==
        (*system)->numbering().local_rank()) {
        const auto* connection =
            fixture->multi_well_context
                ->find_connection(
                    record60.cell_global);
        if (connection == nullptr) {
            local_connection_shapes_ok =
                false;
        } else {
            const auto source60 =
                wdp::
                    build_fixed_bhp_peaceman_well_source_3d(
                        *connection,
                        *final_current[5U]);
            local_connection_shapes_ok =
                local_connection_shapes_ok &&
                source60.cell_source.input_count ==
                    10U;
            for (std::size_t component = 0U;
                 component < 3U;
                 ++component) {
                local_connection60_rate[
                    component] =
                    -source60
                         .cell_source
                         .component_molar_rate_mol_per_s[
                             component];
            }
            local_connection60_rate[3] =
                -source60
                     .cell_source
                     .energy_rate_w;
            local_authoritative_sources
                .push_back(source60);
        }
    }

    require_collective(
        local_connection_shapes_ok &&
            local_authoritative_sources.size() ==
                1U,
        "rebound multi-connection owner did not expose exactly one authoritative connection");

    const auto local_well_total =
        wd::
            aggregate_fixed_bhp_connection_production_rates_3d(
                local_authoritative_sources);
    std::array<double, 4>
        local_well_rate{};
    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        local_well_rate[component] =
            local_well_total
                .component_molar_rate_mol_per_s[
                    component];
    }
    local_well_rate[3] =
        local_well_total.energy_rate_w;

    std::array<double, 4>
        global_well_rate{};
    std::array<double, 4>
        global_connection30_rate{};
    std::array<double, 4>
        global_connection60_rate{};
    require_collective(
        MPI_Allreduce(
            local_well_rate.data(),
            global_well_rate.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            MPI_Allreduce(
                local_connection30_rate.data(),
                global_connection30_rate.data(),
                4,
                MPI_DOUBLE,
                MPI_SUM,
                PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            MPI_Allreduce(
                local_connection60_rate.data(),
                global_connection60_rate.data(),
                4,
                MPI_DOUBLE,
                MPI_SUM,
                PETSC_COMM_WORLD) ==
                MPI_SUCCESS,
        "failed to reduce rebound multi-connection authoritative rates");

    std::uint64_t local_authoritative_count =
        static_cast<std::uint64_t>(
            local_well_total
                .authoritative_connection_count);
    std::uint64_t global_authoritative_count =
        0U;
    require_collective(
        MPI_Allreduce(
            &local_authoritative_count,
            &global_authoritative_count,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            global_authoritative_count ==
                2U,
        "rebound multi-connection well did not retain exactly two authoritative connections");

    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        require_collective(
            std::isfinite(
                global_connection30_rate[
                    quantity]) &&
                std::isfinite(
                    global_connection60_rate[
                        quantity]) &&
                std::abs(
                    global_connection30_rate[
                        quantity]) >
                    0.0 &&
                std::abs(
                    global_connection60_rate[
                        quantity]) >
                    0.0,
            "rebound multi-connection accepted rate vanished or became non-finite");
        near_collective(
            global_well_rate[quantity],
            global_connection30_rate[
                quantity] +
                global_connection60_rate[
                    quantity],
            2.0e-12,
            2.0e-12);
        near_collective(
            global_final_total[quantity],
            global_previous_total[quantity] -
                accepted_dt *
                    global_well_rate[quantity],
            2.0e-7,
            quantity < 3U
                ? 2.0e-8
                : 2.0e-7);
    }

    bool history_matches = false;
    require_collective(
        (*system)
                ->accepted_history_matches_state(
                    (*system)
                        ->initial_state(),
                    &history_matches) ==
            PETSC_SUCCESS &&
            history_matches,
        "rebound multi-connection accepted timestep did not rebase history");

    require_collective(
        VecDestroy(
            &state_delta) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &previous_state) ==
                PETSC_SUCCESS,
        "rebound multi-connection physical timestep cleanup failed");
}

void run_frozen_fixed_bhp_timestep_case(
    int rank,
    std::uint64_t target_stable,
    std::size_t target_local,
    std::size_t expected_phase_count,
    double bottom_hole_pressure_pa,
    const dp::
        ParallelOwnedConnectionSchedule3D&
            schedule,
    const mesh::PartitionSnapshot&
        partition,
    const dp::
        PetscMpiAijSymbolicPreallocation3D&
            bridge,
    const dp::
        OwnedCellStructuralColumnPatternSnapshot3D&
            pattern,
    DispatchAudit* audit) {
    if (audit == nullptr ||
        expected_phase_count == 0U ||
        expected_phase_count > 3U ||
        target_local >= 6U) {
        throw std::invalid_argument(
            "invalid fixed-BHP cardinality timestep case");
    }

    audit->well_timestep_compressibility =
        true;

    std::vector<double>
        injection_enthalpy;
    injection_enthalpy.reserve(
        expected_phase_count);
    for (std::size_t phase = 0U;
         phase < expected_phase_count;
         ++phase) {
        injection_enthalpy.push_back(
            1000.0 *
            static_cast<double>(
                phase + 1U));
    }

    auto timestep_well_context =
        wdp::
            FixedBhpPeacemanWellSourceEvaluatorContext3D::
                create(
                    mesh::GlobalEntityId{
                        target_stable},
                    well::make_peaceman_well_index_3d(
                        {10.0, 10.0, 5.0},
                        {
                            1.0e-8,
                            1.0e-8,
                            1.0e-8},
                        well::
                            AxisAlignedWellDirection3D::z,
                        0.10,
                        0.0),
                    bottom_hole_pressure_pa,
                    wd::FixedBhpInjectionEnthalpy3D{
                        "fixture/fixed-bhp-variable-cardinality-injection-enthalpy/v1",
                        std::move(
                            injection_enthalpy)},
                    "fixture/fixed-bhp-variable-cardinality-source/v1");

    const auto initial_target_current =
        evaluate_target(
            target_stable,
            audit);
    const auto initial_well =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                timestep_well_context,
                initial_target_current);
    require_collective(
        initial_well
                .phase_volumetric_rate_m3_per_s
                .size() ==
            expected_phase_count,
        "fixed-BHP initial production source changed active phase cardinality");
    for (double phase_rate :
         initial_well
             .phase_volumetric_rate_m3_per_s) {
        require_collective(
            phase_rate > 0.0,
            "fixed-BHP regression was not configured production-positive at the frozen initial state");
    }

    const auto timestep_pr_parameters =
        thermodynamic_adapter_pr_parameters();
    const auto timestep_pr_model =
        th::Pr76Phase<double>::
            from_parameters(
                timestep_pr_parameters);
    auto timestep_provider =
        make_outer_rebuild_pr_provider(
            timestep_pr_model);

    FixedBhpWellSourceAudit
        timestep_source_audit{
            &timestep_well_context};
    FixedBhpWellSourceAudit
        verification_source_audit{
            &timestep_well_context};

    auto timestep_system =
        make_frozen_fixed_bhp_timestep_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &timestep_source_audit,
            &timestep_provider);
    auto verification_system =
        make_frozen_fixed_bhp_timestep_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &verification_source_audit,
            &timestep_provider);

    const auto target_cell =
        mesh::LocalIndex{
            static_cast<
                mesh::LocalIndex::value_type>(
                    target_local)};
    const auto& target_record =
        timestep_system
            ->numbering()
            .cell(target_cell);
    const std::size_t expected_scalar_count =
        expected_phase_count * 3U +
        1U;
    require_collective(
        timestep_system
                ->numbering()
                .petsc_global_scalar_count() ==
            42 &&
            target_record.cell_global ==
                mesh::GlobalEntityId{
                    target_stable} &&
            target_record.phase_count ==
                expected_phase_count &&
            target_record.scalar_count ==
                expected_scalar_count,
        "fixed-BHP timestep fixture changed variable-cardinality reservoir numbering");

    const bool local_target_owner =
        target_record.owner_rank ==
        timestep_system
            ->numbering()
            .local_rank();

    const auto local_previous_total =
        owned_conserved_totals(
            *timestep_system,
            timestep_system
                ->initial_state());
    std::array<double, 4>
        global_previous_total{};
    require_collective(
        MPI_Allreduce(
            local_previous_total.data(),
            global_previous_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce variable-cardinality fixed-BHP previous conserved totals");

    Vec previous_state = nullptr;
    require_collective(
        VecDuplicate(
            timestep_system
                ->initial_state(),
            &previous_state) ==
                PETSC_SUCCESS &&
            VecCopy(
                timestep_system
                    ->initial_state(),
                previous_state) ==
                PETSC_SUCCESS,
        "failed to preserve variable-cardinality fixed-BHP initial state");

    fdp::PhysicalTimestepDriverOptions3D
        timestep_options;
    timestep_options.adaptive
        .minimum_timestep_seconds =
        0.125;
    timestep_options.adaptive
        .maximum_timestep_seconds =
        1.0;
    timestep_options.adaptive
        .cutback_factor =
        0.5;
    timestep_options.adaptive
        .growth_factor =
        2.0;
    timestep_options.adaptive
        .maximum_retries =
        4U;
    timestep_options.adaptive
        .growth_nonlinear_iteration_limit =
        20;
    timestep_options.adaptive
        .growth_line_search_direction_change_limit =
        4;
    timestep_options.adaptive
        .growth_transition_restart_limit =
        0U;
    timestep_options.phase_transition
        .max_transition_restarts =
        0U;

    fdp::AcceptedPhysicalTimeClock3D
        timestep_clock{
            0.0,
            1.0};
    FrozenWellTimestepControlAudit
        timestep_control_audit{
            target_cell,
            mesh::GlobalEntityId{
                target_stable},
            expected_phase_count,
            expected_scalar_count,
            0U,
            0U};
    std::optional<
        fdp::PhysicalTimestepDriverReport3D>
        timestep_report;
    PetscErrorCode error =
        fdp::advance_one_physical_timestep_3d(
            PETSC_COMM_WORLD,
            &timestep_system,
            0U,
            {
                &frozen_well_timestep_scan,
                &timestep_control_audit,
                &unexpected_frozen_well_timestep_rebuild,
                &timestep_control_audit},
            timestep_options,
            &timestep_clock,
            &timestep_report);

    require_collective(
        error == PETSC_SUCCESS &&
            timestep_report.has_value() &&
            timestep_report->accepted() &&
            timestep_report
                ->adaptive
                .accepted_timestep_seconds
                .has_value() &&
            timestep_report
                ->accepted_record
                .has_value() &&
            timestep_report
                ->accepted_record
                ->phase_transition_restarts ==
            0U &&
            timestep_clock
                .accepted_step_count() ==
            1U &&
            timestep_control_audit.scans >
            0U &&
            timestep_control_audit
                .rebuild_calls ==
            0U &&
            timestep_system
                ->numbering()
                .cell(target_cell)
                .phase_count ==
            expected_phase_count &&
            (local_target_owner
                 ? timestep_source_audit
                       .target_calls >
                   0U
                 : timestep_source_audit
                       .target_calls ==
                   0U),
        "PhysicalTimestepDriver did not accept the frozen variable-cardinality fixed-BHP well step owner-only");

    const double accepted_dt =
        *timestep_report
             ->adaptive
             .accepted_timestep_seconds;

    bool history_matches = false;
    require_collective(
        timestep_system
                ->accepted_history_matches_state(
                    timestep_system
                        ->initial_state(),
                    &history_matches) ==
            PETSC_SUCCESS &&
            history_matches,
        "variable-cardinality fixed-BHP accepted timestep did not rebase component/energy history");

    Vec state_delta = nullptr;
    PetscReal state_delta_norm = 0.0;
    require_collective(
        VecDuplicate(
            timestep_system
                ->initial_state(),
            &state_delta) ==
                PETSC_SUCCESS &&
            VecCopy(
                timestep_system
                    ->initial_state(),
                state_delta) ==
                PETSC_SUCCESS &&
            VecAXPY(
                state_delta,
                PetscScalar{-1.0},
                previous_state) ==
                PETSC_SUCCESS &&
            VecNorm(
                state_delta,
                NORM_2,
                &state_delta_norm) ==
                PETSC_SUCCESS &&
            static_cast<double>(
                state_delta_norm) >
                1.0e-8,
        "variable-cardinality fixed-BHP well did not change the accepted reservoir state");

    require_collective(
        verification_system
                ->set_trial_timestep_seconds(
                    accepted_dt) ==
            PETSC_SUCCESS,
        "failed to align independent variable-cardinality fixed-BHP verification timestep");

    Vec verification_residual = nullptr;
    require_collective(
        VecDuplicate(
            timestep_system
                ->initial_state(),
            &verification_residual) ==
                PETSC_SUCCESS &&
            VecSet(
                verification_residual,
                PetscScalar{0.0}) ==
                PETSC_SUCCESS,
        "failed to allocate independent variable-cardinality fixed-BHP residual");

    auto verification_evaluator =
        verification_system
            ->snes_evaluator();
    fdp::NaturalVariableSnesEvaluationStatus3D
        verification_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    error =
        verification_evaluator.function(
            timestep_system
                ->initial_state(),
            verification_residual,
            verification_evaluator
                .user_context,
            &verification_status);
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyBegin(
                verification_residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyEnd(
                verification_residual);
    }
    PetscReal verification_residual_norm =
        0.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                verification_residual,
                NORM_2,
                &verification_residual_norm);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            verification_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            static_cast<double>(
                verification_residual_norm) <=
                1.0e-6,
        "accepted variable-cardinality fixed-BHP timestep failed independent residual reassembly");

    const auto local_final_total =
        owned_conserved_totals(
            *verification_system,
            timestep_system
                ->initial_state());
    std::array<double, 4>
        global_final_total{};
    require_collective(
        MPI_Allreduce(
            local_final_total.data(),
            global_final_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce variable-cardinality fixed-BHP final conserved totals");

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        final_current;
    std::vector<double>
        final_porosities;
    fdp::NaturalVariableSnesEvaluationStatus3D
        final_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    const PetscErrorCode final_error =
        verification_system
            ->evaluate_local_cells_for_phase_transition(
                timestep_system
                    ->initial_state(),
                &final_current,
                &final_porosities,
                &final_status);
    require_collective(
        final_error == PETSC_SUCCESS &&
            final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            final_current.size() >
                target_local &&
            final_porosities.size() ==
                final_current.size() &&
            final_current[target_local]
                .has_value() &&
            std::visit(
                [](const auto& typed) {
                    return typed
                        .transport
                        .state_identity
                        .layout
                        .phase_count();
                },
                *final_current[
                    target_local]) ==
                expected_phase_count,
        "failed to collectively recover accepted variable-cardinality well-cell state");

    std::array<double, 4>
        local_well_production_rate{};
    if (local_target_owner) {
        const auto final_well =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    timestep_well_context,
                    *final_current[
                        target_local]);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            local_well_production_rate[
                component] =
                -final_well
                     .cell_source
                     .component_molar_rate_mol_per_s[
                         component];
        }
        local_well_production_rate[3] =
            -final_well
                 .cell_source
                 .energy_rate_w;
    }

    std::array<double, 4>
        global_well_production_rate{};
    require_collective(
        MPI_Allreduce(
            local_well_production_rate.data(),
            global_well_production_rate.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce variable-cardinality fixed-BHP production rates");

    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        require_collective(
            std::isfinite(
                global_well_production_rate[
                    component]) &&
                std::abs(
                    global_well_production_rate[
                        component]) >
                    0.0,
            "accepted variable-cardinality fixed-BHP component rate vanished or became non-finite");
        near_collective(
            global_final_total[
                component],
            global_previous_total[
                component] -
                accepted_dt *
                    global_well_production_rate[
                        component],
            2.0e-7,
            2.0e-8);
    }
    require_collective(
        std::isfinite(
            global_well_production_rate[3]) &&
            std::abs(
                global_well_production_rate[3]) >
                0.0,
        "accepted variable-cardinality fixed-BHP energy rate vanished or became non-finite");
    near_collective(
        global_final_total[3],
        global_previous_total[3] -
            accepted_dt *
                global_well_production_rate[
                    3],
        2.0e-7,
        2.0e-7);

    require_collective(
        VecDestroy(
            &verification_residual) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &state_delta) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &previous_state) ==
                PETSC_SUCCESS,
        "variable-cardinality fixed-BHP timestep regression cleanup failed");

    audit->well_timestep_compressibility =
        false;
}

void run_multi_connection_fixed_bhp_timestep_case(
    int rank,
    const dp::
        ParallelOwnedConnectionSchedule3D&
            schedule,
    const mesh::PartitionSnapshot&
        partition,
    const dp::
        PetscMpiAijSymbolicPreallocation3D&
            bridge,
    const dp::
        OwnedCellStructuralColumnPatternSnapshot3D&
            pattern,
    DispatchAudit* audit) {
    if (audit == nullptr) {
        throw std::invalid_argument(
            "invalid multi-connection fixed-BHP regression context");
    }

    audit->well_timestep_compressibility =
        true;

    constexpr double shared_bhp_pa = 5.0;
    auto connection30 =
        wdp::
            FixedBhpPeacemanWellSourceEvaluatorContext3D::
                create(
                    mesh::GlobalEntityId{
                        UINT64_C(30)},
                    well::make_peaceman_well_index_3d(
                        {10.0, 10.0, 5.0},
                        {
                            1.0e-8,
                            1.0e-8,
                            1.0e-8},
                        well::
                            AxisAlignedWellDirection3D::z,
                        0.10,
                        0.0),
                    shared_bhp_pa,
                    wd::FixedBhpInjectionEnthalpy3D{
                        "fixture/multi-well-cell30-injection-enthalpy/v1",
                        {1000.0, 2000.0}},
                    "fixture/multi-well-cell30-source/v1");
    auto connection60 =
        wdp::
            FixedBhpPeacemanWellSourceEvaluatorContext3D::
                create(
                    mesh::GlobalEntityId{
                        UINT64_C(60)},
                    well::make_peaceman_well_index_3d(
                        {12.0, 9.0, 4.0},
                        {
                            2.0e-8,
                            1.5e-8,
                            2.5e-8},
                        well::
                            AxisAlignedWellDirection3D::z,
                        0.12,
                        0.1),
                    shared_bhp_pa,
                    wd::FixedBhpInjectionEnthalpy3D{
                        "fixture/multi-well-cell60-injection-enthalpy/v1",
                        {1100.0, 2200.0, 3300.0}},
                    "fixture/multi-well-cell60-source/v1");

    bool duplicate_rejected = false;
    try {
        (void)wdp::
            FixedBhpMultiConnectionWellSourceEvaluatorContext3D::
                create(
                    "fixture/duplicate-completion",
                    {connection30, connection30});
    } catch (const std::invalid_argument&) {
        duplicate_rejected = true;
    }

    bool bhp_mismatch_rejected = false;
    try {
        auto mismatched =
            wdp::
                FixedBhpPeacemanWellSourceEvaluatorContext3D::
                    create(
                        mesh::GlobalEntityId{
                            UINT64_C(60)},
                        connection60.connection(),
                        shared_bhp_pa + 1.0,
                        connection60.injection_enthalpy(),
                        "fixture/multi-well-bhp-mismatch/v1");
        (void)wdp::
            FixedBhpMultiConnectionWellSourceEvaluatorContext3D::
                create(
                    "fixture/mismatched-bhp",
                    {connection30, mismatched});
    } catch (const std::invalid_argument&) {
        bhp_mismatch_rejected = true;
    }

    require_collective(
        duplicate_rejected &&
            bhp_mismatch_rejected,
        "multi-connection fixed-BHP context accepted duplicate stable cell or inconsistent BHP");

    auto multi_context =
        wdp::
            FixedBhpMultiConnectionWellSourceEvaluatorContext3D::
                create(
                    "fixture/single-well-two-connections/v1",
                    {connection60, connection30});

    require_collective(
        multi_context.connection_count() ==
                2U &&
            multi_context
                    .bottom_hole_pressure_pa() ==
                shared_bhp_pa &&
            multi_context
                    .connections()[0U]
                    .target_cell_global() ==
                mesh::GlobalEntityId{
                    UINT64_C(30)} &&
            multi_context
                    .connections()[1U]
                    .target_cell_global() ==
                mesh::GlobalEntityId{
                    UINT64_C(60)} &&
            multi_context.find_connection(
                mesh::GlobalEntityId{
                    UINT64_C(40)}) ==
                nullptr,
        "multi-connection fixed-BHP context did not retain deterministic unique stable-cell lookup");

    auto initial30 =
        evaluate_target(
            UINT64_C(30),
            audit);
    auto initial60 =
        evaluate_target(
            UINT64_C(60),
            audit);
    const auto initial_source30 =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                *multi_context.find_connection(
                    mesh::GlobalEntityId{
                        UINT64_C(30)}),
                initial30);
    const auto initial_source60 =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                *multi_context.find_connection(
                    mesh::GlobalEntityId{
                        UINT64_C(60)}),
                initial60);
    const std::vector<
        wd::
            FixedBhpConnectionCellSourceLinearization3D>
        direct_connection_sources{
            initial_source30,
            initial_source60};
    const auto direct_total =
        wd::
            aggregate_fixed_bhp_connection_production_rates_3d(
                direct_connection_sources);

    bool initial_production_positive =
        initial_source30
                .cell_source
                .input_count ==
            7U &&
        initial_source60
                .cell_source
                .input_count ==
            10U &&
        direct_total
                .authoritative_connection_count ==
            2U &&
        direct_total.bottom_hole_pressure_pa ==
            shared_bhp_pa;
    for (double rate :
         initial_source30
             .phase_volumetric_rate_m3_per_s) {
        initial_production_positive =
            initial_production_positive &&
            rate > 0.0;
    }
    for (double rate :
         initial_source60
             .phase_volumetric_rate_m3_per_s) {
        initial_production_positive =
            initial_production_positive &&
            rate > 0.0;
    }
    require_collective(
        initial_production_positive,
        "multi-connection fixed-BHP fixture was not production-positive with q=7/q=10 connection Jacobians");

    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        near_collective(
            direct_total
                .component_molar_rate_mol_per_s[
                    component],
            -initial_source30
                 .cell_source
                 .component_molar_rate_mol_per_s[
                     component] -
                initial_source60
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component],
            1.0e-12,
            1.0e-12);
    }
    near_collective(
        direct_total.energy_rate_w,
        -initial_source30
             .cell_source
             .energy_rate_w -
            initial_source60
                .cell_source
                .energy_rate_w,
        1.0e-12,
        1.0e-12);

    const auto pr_parameters =
        thermodynamic_adapter_pr_parameters();
    const auto pr_model =
        th::Pr76Phase<double>::
            from_parameters(
                pr_parameters);
    auto provider =
        make_outer_rebuild_pr_provider(
            pr_model);

    FixedBhpMultiConnectionWellSourceAudit
        source_audit{
            &multi_context};
    auto system =
        make_frozen_multi_connection_fixed_bhp_timestep_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &source_audit,
            &provider);

    const auto& record30 =
        system->numbering().cell(
            mesh::LocalIndex{2U});
    const auto& record60 =
        system->numbering().cell(
            mesh::LocalIndex{5U});
    require_collective(
        system->numbering()
                .petsc_global_scalar_count() ==
            42 &&
            record30.phase_count == 2U &&
            record30.scalar_count == 7U &&
            record30.owner_rank.value() == 0U &&
            record60.phase_count == 3U &&
            record60.scalar_count == 10U &&
            record60.owner_rank.value() == 1U,
        "multi-connection fixed-BHP well changed reservoir global numbering or target ownership");

    const auto local_previous_total =
        owned_conserved_totals(
            *system,
            system->initial_state());
    std::array<double, 4>
        global_previous_total{};
    require_collective(
        MPI_Allreduce(
            local_previous_total.data(),
            global_previous_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce multi-connection fixed-BHP previous conserved totals");

    Vec previous_state = nullptr;
    require_collective(
        VecDuplicate(
            system->initial_state(),
            &previous_state) ==
                PETSC_SUCCESS &&
            VecCopy(
                system->initial_state(),
                previous_state) ==
                PETSC_SUCCESS,
        "failed to preserve multi-connection fixed-BHP previous state");

    fdp::PhysicalTimestepDriverOptions3D
        timestep_options;
    timestep_options.adaptive
        .minimum_timestep_seconds =
        0.125;
    timestep_options.adaptive
        .maximum_timestep_seconds =
        1.0;
    timestep_options.adaptive
        .cutback_factor =
        0.5;
    timestep_options.adaptive
        .growth_factor =
        2.0;
    timestep_options.adaptive
        .maximum_retries =
        4U;
    timestep_options.adaptive
        .growth_nonlinear_iteration_limit =
        20;
    timestep_options.adaptive
        .growth_line_search_direction_change_limit =
        4;
    timestep_options.adaptive
        .growth_transition_restart_limit =
        0U;
    timestep_options.phase_transition
        .max_transition_restarts =
        0U;

    FrozenMultiWellTimestepControlAudit
        control_audit;
    fdp::AcceptedPhysicalTimeClock3D
        timestep_clock{
            0.0,
            1.0};
    std::optional<
        fdp::PhysicalTimestepDriverReport3D>
        timestep_report;
    PetscErrorCode error =
        fdp::advance_one_physical_timestep_3d(
            PETSC_COMM_WORLD,
            &system,
            0U,
            {
                &frozen_multi_well_timestep_scan,
                &control_audit,
                &unexpected_frozen_multi_well_timestep_rebuild,
                &control_audit},
            timestep_options,
            &timestep_clock,
            &timestep_report);

    require_collective(
        error == PETSC_SUCCESS &&
            timestep_report.has_value() &&
            timestep_report->accepted() &&
            timestep_report
                    ->adaptive
                    .accepted_timestep_seconds
                    .has_value() &&
            timestep_clock
                    .accepted_step_count() ==
                1U &&
            control_audit.scans > 0U &&
            control_audit.rebuild_calls == 0U &&
            source_audit.evaluator_calls > 0U &&
            (rank == 0
                 ? source_audit.cell30_calls > 0U &&
                       source_audit.cell60_calls == 0U
                 : source_audit.cell30_calls == 0U &&
                       source_audit.cell60_calls > 0U),
        "multi-connection fixed-BHP sources were not assembled exactly on their authoritative owners");

    const double accepted_dt =
        *timestep_report
             ->adaptive
             .accepted_timestep_seconds;

    Vec state_delta = nullptr;
    PetscReal state_delta_norm = 0.0;
    require_collective(
        VecDuplicate(
            system->initial_state(),
            &state_delta) ==
                PETSC_SUCCESS &&
            VecCopy(
                system->initial_state(),
                state_delta) ==
                PETSC_SUCCESS &&
            VecAXPY(
                state_delta,
                PetscScalar{-1.0},
                previous_state) ==
                PETSC_SUCCESS &&
            VecNorm(
                state_delta,
                NORM_2,
                &state_delta_norm) ==
                PETSC_SUCCESS &&
            static_cast<double>(
                state_delta_norm) >
                1.0e-10,
        "multi-connection fixed-BHP accepted timestep did not change reservoir state");

    const auto local_final_total =
        owned_conserved_totals(
            *system,
            system->initial_state());
    std::array<double, 4>
        global_final_total{};
    require_collective(
        MPI_Allreduce(
            local_final_total.data(),
            global_final_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce multi-connection fixed-BHP final conserved totals");

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        final_current;
    std::vector<double>
        final_porosity;
    fdp::NaturalVariableSnesEvaluationStatus3D
        final_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    require_collective(
        system->evaluate_local_cells_for_phase_transition(
                system->initial_state(),
                &final_current,
                &final_porosity,
                &final_status) ==
            PETSC_SUCCESS &&
            final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            final_current.size() > 5U &&
            final_current[2U].has_value() &&
            final_current[5U].has_value(),
        "failed to evaluate accepted multi-connection fixed-BHP states");

    std::vector<
        wd::
            FixedBhpConnectionCellSourceLinearization3D>
        local_authoritative_sources;
    std::array<double, 4>
        local_connection30_rate{};
    std::array<double, 4>
        local_connection60_rate{};
    bool local_connection_shapes_ok =
        true;

    if (record30.owner_rank ==
        system->numbering().local_rank()) {
        const auto source30 =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    *multi_context.find_connection(
                        record30.cell_global),
                    *final_current[2U]);
        local_connection_shapes_ok =
            local_connection_shapes_ok &&
            source30.cell_source.input_count ==
                7U;
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            local_connection30_rate[
                component] =
                -source30
                     .cell_source
                     .component_molar_rate_mol_per_s[
                         component];
        }
        local_connection30_rate[3] =
            -source30.cell_source.energy_rate_w;
        local_authoritative_sources.push_back(
            source30);
    }

    if (record60.owner_rank ==
        system->numbering().local_rank()) {
        const auto source60 =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    *multi_context.find_connection(
                        record60.cell_global),
                    *final_current[5U]);
        local_connection_shapes_ok =
            local_connection_shapes_ok &&
            source60.cell_source.input_count ==
                10U;
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            local_connection60_rate[
                component] =
                -source60
                     .cell_source
                     .component_molar_rate_mol_per_s[
                         component];
        }
        local_connection60_rate[3] =
            -source60.cell_source.energy_rate_w;
        local_authoritative_sources.push_back(
            source60);
    }

    require_collective(
        local_connection_shapes_ok &&
            local_authoritative_sources.size() ==
                1U,
        "multi-connection fixed-BHP authoritative owner did not expose exactly one local connection");

    const auto local_well_total =
        wd::
            aggregate_fixed_bhp_connection_production_rates_3d(
                local_authoritative_sources);
    std::array<double, 4>
        local_well_rate{};
    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        local_well_rate[component] =
            local_well_total
                .component_molar_rate_mol_per_s[
                    component];
    }
    local_well_rate[3] =
        local_well_total.energy_rate_w;

    std::array<double, 4>
        global_well_rate{};
    std::array<double, 4>
        global_connection30_rate{};
    std::array<double, 4>
        global_connection60_rate{};
    require_collective(
        MPI_Allreduce(
            local_well_rate.data(),
            global_well_rate.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            MPI_Allreduce(
                local_connection30_rate.data(),
                global_connection30_rate.data(),
                4,
                MPI_DOUBLE,
                MPI_SUM,
                PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            MPI_Allreduce(
                local_connection60_rate.data(),
                global_connection60_rate.data(),
                4,
                MPI_DOUBLE,
                MPI_SUM,
                PETSC_COMM_WORLD) ==
                MPI_SUCCESS,
        "failed to reduce multi-connection fixed-BHP authoritative rates");

    std::uint64_t local_authoritative_count =
        static_cast<std::uint64_t>(
            local_well_total
                .authoritative_connection_count);
    std::uint64_t global_authoritative_count =
        0U;
    require_collective(
        MPI_Allreduce(
            &local_authoritative_count,
            &global_authoritative_count,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            global_authoritative_count ==
                2U,
        "multi-connection fixed-BHP well did not aggregate exactly two authoritative connections");

    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        require_collective(
            std::isfinite(
                global_connection30_rate[
                    quantity]) &&
                std::isfinite(
                    global_connection60_rate[
                        quantity]) &&
                std::abs(
                    global_connection30_rate[
                        quantity]) >
                    0.0 &&
                std::abs(
                    global_connection60_rate[
                        quantity]) >
                    0.0,
            "multi-connection fixed-BHP accepted connection rate vanished or became non-finite");
        near_collective(
            global_well_rate[quantity],
            global_connection30_rate[
                quantity] +
                global_connection60_rate[
                    quantity],
            2.0e-12,
            2.0e-12);
        near_collective(
            global_final_total[quantity],
            global_previous_total[quantity] -
                accepted_dt *
                    global_well_rate[quantity],
            2.0e-7,
            quantity < 3U
                ? 2.0e-8
                : 2.0e-7);
    }

    bool history_matches = false;
    require_collective(
        system->accepted_history_matches_state(
                system->initial_state(),
                &history_matches) ==
            PETSC_SUCCESS &&
            history_matches,
        "multi-connection fixed-BHP accepted timestep did not rebase history");

    require_collective(
        VecDestroy(
            &state_delta) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &previous_state) ==
                PETSC_SUCCESS,
        "multi-connection fixed-BHP regression cleanup failed");

    audit->well_timestep_compressibility =
        false;
}

void run_multi_connection_local_transition_rebind_case(
    int rank,
    const dp::
        ParallelOwnedConnectionSchedule3D&
            schedule,
    const mesh::PartitionSnapshot&
        partition,
    const dp::
        PetscMpiAijSymbolicPreallocation3D&
            bridge,
    const dp::
        OwnedCellStructuralColumnPatternSnapshot3D&
            pattern,
    DispatchAudit* audit,
    flow::
        Pr76AbsentPhasePotentialExtensionProvider<
            double>* provider) {
    if (audit == nullptr ||
        provider == nullptr) {
        throw std::invalid_argument(
            "invalid local-transition multi-connection fixed-BHP fixture");
    }

    audit->well_timestep_compressibility =
        true;

    const auto aqueous =
        mixed_physical_phase_identity(
            "aqueous");
    const auto hydrocarbon0 =
        mixed_physical_phase_identity(
            "hydrocarbon-0");
    const auto hydrocarbon1 =
        mixed_physical_phase_identity(
            "hydrocarbon-1");

    wdp::FixedBhpPhaseIdentityInjectionEnthalpy3D
        cell30_phase_registry{
            "fixture/multi-transition-cell30-phase-enthalpy/v1",
            {
                {aqueous, 1000.0},
                {hydrocarbon0, 2000.0},
                {hydrocarbon1, 3000.0}}};

    constexpr double shared_bhp_pa =
        5.0;
    auto connection30 =
        wdp::
            FixedBhpPeacemanWellSourceEvaluatorContext3D::
                create_phase_identity_bound(
                    mesh::GlobalEntityId{
                        UINT64_C(30)},
                    well::make_peaceman_well_index_3d(
                        {10.0, 10.0, 5.0},
                        {
                            1.0e-8,
                            1.0e-8,
                            1.0e-8},
                        well::
                            AxisAlignedWellDirection3D::z,
                        0.10,
                        0.0),
                    shared_bhp_pa,
                    cell30_phase_registry,
                    phase_identity_map_for_stable(
                        UINT64_C(30)),
                    "fixture/multi-transition-cell30-source/v1");

    auto connection60 =
        wdp::
            FixedBhpPeacemanWellSourceEvaluatorContext3D::
                create(
                    mesh::GlobalEntityId{
                        UINT64_C(60)},
                    well::make_peaceman_well_index_3d(
                        {12.0, 9.0, 4.0},
                        {
                            2.0e-8,
                            1.5e-8,
                            2.5e-8},
                        well::
                            AxisAlignedWellDirection3D::z,
                        0.12,
                        0.1),
                    shared_bhp_pa,
                    wd::FixedBhpInjectionEnthalpy3D{
                        "fixture/multi-transition-cell60-enthalpy/v1",
                        {1100.0, 2200.0, 3300.0}},
                    "fixture/multi-transition-cell60-source/v1");

    const double frozen60_wi =
        connection60.connection()
            .well_index_m3;
    const auto frozen60_enthalpy =
        connection60
            .injection_enthalpy()
            .specific_enthalpy_j_per_kg;
    const std::string frozen60_source =
        connection60.source_provenance();

    auto multi_context =
        wdp::
            FixedBhpMultiConnectionWellSourceEvaluatorContext3D::
                create(
                    "fixture/multi-transition-single-well/v1",
                    {
                        connection60,
                        connection30});

    ControllerFixture fixture{
        rank,
        &schedule,
        &partition,
        &bridge,
        &pattern,
        audit,
        nullptr,
        provider,
        false,
        0U,
        false,
        false};
    fixture.multi_well_context =
        &multi_context;

    auto initial_system =
        make_controller_initial_system(
            &fixture);
    require_collective(
        initial_system != nullptr &&
            initial_system
                    ->numbering()
                    .petsc_global_scalar_count() ==
                42 &&
            initial_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                2U &&
            initial_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .scalar_count ==
                7U &&
            initial_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{5U})
                    .phase_count ==
                3U &&
            initial_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{5U})
                    .scalar_count ==
                10U,
        "local-transition multi-well fixture did not start on 2P(cell30)+3P(cell60)");

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        appeared_system;
    Vec appeared_state = nullptr;
    std::optional<
        fdp::
            PostSnesPhaseTransitionControllerReport3D>
        appearance_report;
    PetscErrorCode error =
        fdp::
            solve_nonlinear_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                std::move(initial_system),
                {
                    &controller_scan,
                    &fixture,
                    &controller_rebuild,
                    &fixture},
                {4U},
                &appeared_system,
                &appeared_state,
                &appearance_report);

    require_collective(
        error == PETSC_SUCCESS &&
            appeared_system != nullptr &&
            appeared_state != nullptr &&
            appearance_report.has_value() &&
            appearance_report->outcome ==
                fdp::
                    PostSnesPhaseTransitionOutcome3D::
                        stable_phase_set &&
            appearance_report
                    ->transition_restarts ==
                1U &&
            fixture.rebuild_calls ==
                1U &&
            fixture.multi_well_rebound &&
            fixture.multi_well_context !=
                nullptr &&
            appeared_system
                    ->numbering()
                    .petsc_global_scalar_count() ==
                45 &&
            appeared_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                3U &&
            appeared_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .scalar_count ==
                10U &&
            appeared_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{5U})
                    .phase_count ==
                3U &&
            appeared_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{5U})
                    .scalar_count ==
                10U,
        "multi-connection well did not survive local cell30 2P-to-3P restart");

    const auto* appeared30 =
        fixture.multi_well_context
            ->find_connection(
                mesh::GlobalEntityId{
                    UINT64_C(30)});
    const auto* frozen60 =
        fixture.multi_well_context
            ->find_connection(
                mesh::GlobalEntityId{
                    UINT64_C(60)});
    require_collective(
        appeared30 != nullptr &&
            frozen60 != nullptr &&
            appeared30
                    ->active_phase_identities()
                    .has_value() &&
            appeared30
                    ->active_phase_identities()
                    ->phase_count() ==
                3U &&
            appeared30
                    ->injection_enthalpy()
                    .specific_enthalpy_j_per_kg ==
                std::vector<double>{
                    1000.0,
                    2000.0,
                    3000.0} &&
            frozen60
                    ->target_cell_global() ==
                mesh::GlobalEntityId{
                    UINT64_C(60)} &&
            std::abs(
                frozen60
                        ->connection()
                        .well_index_m3 -
                    frozen60_wi) <=
                1.0e-14 *
                    std::max(
                        1.0,
                        std::abs(
                            frozen60_wi)) &&
            frozen60
                    ->bottom_hole_pressure_pa() ==
                shared_bhp_pa &&
            frozen60
                    ->injection_enthalpy()
                    .specific_enthalpy_j_per_kg ==
                frozen60_enthalpy &&
            frozen60
                    ->source_provenance() ==
                frozen60_source &&
            !frozen60
                 ->phase_identity_rebindable(),
        "cell30 local appearance incorrectly mutated frozen cell60 connection");

    run_rebound_multi_connection_fixed_bhp_physical_timestep(
        rank,
        &appeared_system,
        &appeared_state,
        &fixture,
        3U,
        10U,
        45);

    // A second independent controller session removes only cell30's third
    // phase. Cell60 remains a frozen 3P connection throughout the round trip.
    fixture.phase_disappearance_only =
        true;
    fixture.oscillate_after_restart =
        false;
    fixture.rebuild_calls =
        0U;
    fixture.multi_well_rebound =
        false;
    fixture
        .multi_well_source_audit
        .evaluator_calls = 0U;
    fixture
        .multi_well_source_audit
        .cell30_calls = 0U;
    fixture
        .multi_well_source_audit
        .cell60_calls = 0U;

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        disappeared_system;
    Vec disappeared_state = nullptr;
    std::optional<
        fdp::
            PostSnesPhaseTransitionControllerReport3D>
        disappearance_report;
    error =
        fdp::
            solve_nonlinear_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                std::move(appeared_system),
                {
                    &controller_scan,
                    &fixture,
                    &controller_rebuild,
                    &fixture},
                {4U},
                &disappeared_system,
                &disappeared_state,
                &disappearance_report);

    require_collective(
        error == PETSC_SUCCESS &&
            disappeared_system != nullptr &&
            disappeared_state != nullptr &&
            disappearance_report.has_value() &&
            disappearance_report->outcome ==
                fdp::
                    PostSnesPhaseTransitionOutcome3D::
                        stable_phase_set &&
            disappearance_report
                    ->transition_restarts ==
                1U &&
            fixture.rebuild_calls ==
                1U &&
            fixture.multi_well_rebound &&
            fixture.multi_well_context !=
                nullptr &&
            disappeared_system
                    ->numbering()
                    .petsc_global_scalar_count() ==
                42 &&
            disappeared_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                2U &&
            disappeared_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .scalar_count ==
                7U &&
            disappeared_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{5U})
                    .phase_count ==
                3U &&
            disappeared_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{5U})
                    .scalar_count ==
                10U,
        "multi-connection well did not survive local cell30 3P-to-2P restart");

    const auto* disappeared30 =
        fixture.multi_well_context
            ->find_connection(
                mesh::GlobalEntityId{
                    UINT64_C(30)});
    frozen60 =
        fixture.multi_well_context
            ->find_connection(
                mesh::GlobalEntityId{
                    UINT64_C(60)});
    const auto initial30_map =
        phase_identity_map_for_stable(
            UINT64_C(30));
    bool cell30_roundtrip_identity_ok =
        disappeared30 != nullptr &&
        disappeared30
                ->active_phase_identities()
                .has_value() &&
        disappeared30
                ->active_phase_identities()
                ->phase_count() ==
            initial30_map.phase_count();
    if (cell30_roundtrip_identity_ok) {
        for (std::size_t phase = 0U;
             phase < initial30_map.phase_count();
             ++phase) {
            cell30_roundtrip_identity_ok =
                cell30_roundtrip_identity_ok &&
                disappeared30
                        ->active_phase_identities()
                        ->identity(phase) ==
                    initial30_map.identity(
                        phase);
        }
    }

    const auto& retained_registry =
        disappeared30 != nullptr
            ? disappeared30
                  ->phase_identity_injection_enthalpy()
            : std::optional<
                  wdp::
                      FixedBhpPhaseIdentityInjectionEnthalpy3D>{};
    const bool disappeared_phase_retained =
        retained_registry.has_value() &&
        retained_registry->phases.size() ==
            3U &&
        std::any_of(
            retained_registry
                ->phases.begin(),
            retained_registry
                ->phases.end(),
            [&](const auto& entry) {
                return entry.identity ==
                           hydrocarbon1 &&
                    entry
                            .specific_enthalpy_j_per_kg ==
                        3000.0;
            });

    require_collective(
        cell30_roundtrip_identity_ok &&
            disappeared_phase_retained &&
            disappeared30
                    ->injection_enthalpy()
                    .specific_enthalpy_j_per_kg ==
                std::vector<double>{
                    1000.0,
                    2000.0} &&
            frozen60 != nullptr &&
            frozen60
                    ->target_cell_global() ==
                mesh::GlobalEntityId{
                    UINT64_C(60)} &&
            std::abs(
                frozen60
                        ->connection()
                        .well_index_m3 -
                    frozen60_wi) <=
                1.0e-14 *
                    std::max(
                        1.0,
                        std::abs(
                            frozen60_wi)) &&
            frozen60
                    ->bottom_hole_pressure_pa() ==
                shared_bhp_pa &&
            frozen60
                    ->injection_enthalpy()
                    .specific_enthalpy_j_per_kg ==
                frozen60_enthalpy &&
            frozen60
                    ->source_provenance() ==
                frozen60_source &&
            !frozen60
                 ->phase_identity_rebindable(),
        "cell30 local disappearance drifted identity or mutated frozen cell60 connection");

    run_rebound_multi_connection_fixed_bhp_physical_timestep(
        rank,
        &disappeared_system,
        &disappeared_state,
        &fixture,
        2U,
        7U,
        42);

    audit->well_timestep_compressibility =
        false;
}

struct ControlledWellPhaseTransitionFixture {
    int rank{};
    const dp::
        ParallelOwnedConnectionSchedule3D*
            schedule{};
    const mesh::PartitionSnapshot*
        partition{};
    const dp::
        PetscMpiAijSymbolicPreallocation3D*
            bridge{};
    const dp::
        OwnedCellStructuralColumnPatternSnapshot3D*
            pattern{};
    DispatchAudit* audit{};
    flow::
        Pr76AbsentPhasePotentialExtensionProvider<
            double>*
        provider{};
    wdp::
        FixedTotalMolarRateWellSourceEvaluatorContext3D*
        source_context{};
    wdp::
        AcceptedFixedTotalMolarRateWellControlState3D*
        accepted_control{};
    fdp::AcceptedPhysicalTimeClock3D*
        clock{};

    wdp::
        FixedTotalMolarRatePhysicalTimestepControlMode3D
        expected_entry_control{
            wdp::
                FixedTotalMolarRatePhysicalTimestepControlMode3D::
                    fixed_total_molar_rate};
    double expected_entry_bhp_pa{};

    std::size_t scans{};
    std::size_t rebuild_calls{};
    bool accepted_anchor_unchanged_at_rebuild{
        true};
    bool rebuilt_previous_totals_captured{};
    std::array<double, 4>
        rebuilt_previous_totals{};
    bool accepted_anchor_unchanged_at_scan{
        true};
    std::vector<double>
        scanned_source_bottom_hole_pressure_pa;
};

PetscErrorCode
controlled_well_transition_scan(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&
            system,
    Vec,
    const fdp::
        VariableCardinalityNaturalVariableSnesSolveReport3D&
            solve_report,
    void* raw_context,
    fdp::PostSnesPhaseTransitionScanStatus3D*
        scan_status,
    std::vector<
        fdp::PostSnesPhaseTransitionProposal3D>*
            output) {
    if (raw_context == nullptr ||
        scan_status == nullptr ||
        output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* fixture =
        static_cast<
            ControlledWellPhaseTransitionFixture*>(
                raw_context);
    *scan_status =
        fdp::
            PostSnesPhaseTransitionScanStatus3D::
                complete;
    output->clear();

    if (static_cast<int>(
            solve_report.converged_reason) <=
            0 ||
        solve_report
                .locally_owned_solution
                .size() !=
            static_cast<std::size_t>(
                system
                    .numbering()
                    .petsc_local_owned_scalar_count())) {
        return PETSC_ERR_ARG_INCOMP;
    }

    if (fixture->source_context == nullptr ||
        fixture->accepted_control == nullptr ||
        fixture->clock == nullptr) {
        return PETSC_ERR_ARG_INCOMP;
    }
    fixture
        ->accepted_anchor_unchanged_at_scan =
        fixture
            ->accepted_anchor_unchanged_at_scan &&
        fixture
                ->accepted_control
                ->control ==
            fixture
                ->expected_entry_control &&
        fixture
                ->accepted_control
                ->bottom_hole_pressure_pa ==
            fixture
                ->expected_entry_bhp_pa &&
        fixture
                ->clock
                ->accepted_time_seconds() ==
            0.0 &&
        fixture
                ->clock
                ->accepted_step_count() ==
            0U;
    fixture
        ->scanned_source_bottom_hole_pressure_pa
        .push_back(
            fixture
                ->source_context
                ->current_bottom_hole_pressure_pa());

    ++fixture->scans;
    const auto& cell30 =
        system.numbering().cell(
            mesh::LocalIndex{2U});
    if (cell30.owner_rank !=
        system.numbering().local_rank()) {
        return PETSC_SUCCESS;
    }
    if (cell30.cell_global !=
            mesh::GlobalEntityId{
                UINT64_C(30)}) {
        return PETSC_ERR_ARG_INCOMP;
    }

    if (cell30.phase_count == 2U) {
        output->push_back(
            controller_two_to_three_proposal());
    } else if (cell30.phase_count != 3U) {
        return PETSC_ERR_ARG_INCOMP;
    }
    return PETSC_SUCCESS;
}

PetscErrorCode
controlled_well_transition_rebuild(
    const fdp::
        PhaseTransitionRebuiltNaturalVariableSystem3D&
            current_system,
    Vec,
    const fdp::
        VariableCardinalityNaturalVariableSnesSolveReport3D&,
    std::span<
        const fdp::
            PostSnesPhaseTransitionProposal3D>
        local_owned_proposals,
    std::span<
        const fdp::
            AcceptedPhaseTransitionSummary3D>
        accepted_global_batch,
    void* raw_context,
    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>*
        output) {
    if (raw_context == nullptr ||
        output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();
    auto* fixture =
        static_cast<
            ControlledWellPhaseTransitionFixture*>(
                raw_context);
    if (fixture->schedule == nullptr ||
        fixture->partition == nullptr ||
        fixture->bridge == nullptr ||
        fixture->pattern == nullptr ||
        fixture->audit == nullptr ||
        fixture->provider == nullptr ||
        fixture->source_context == nullptr ||
        fixture->accepted_control == nullptr ||
        fixture->clock == nullptr ||
        accepted_global_batch.size() !=
            1U ||
        accepted_global_batch.front()
                .cell_global !=
            mesh::GlobalEntityId{
                UINT64_C(30)} ||
        accepted_global_batch.front()
                .source_phase_count !=
            2U ||
        accepted_global_batch.front()
                .target_phase_count !=
            3U) {
        return PETSC_ERR_ARG_INCOMP;
    }

    if (fixture->rank == 0) {
        if (local_owned_proposals.size() !=
                1U ||
            local_owned_proposals.front()
                    .cell_global !=
                mesh::GlobalEntityId{
                    UINT64_C(30)}) {
            return PETSC_ERR_ARG_INCOMP;
        }
    } else if (
        !local_owned_proposals.empty()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    fixture
        ->accepted_anchor_unchanged_at_rebuild =
        fixture
            ->accepted_anchor_unchanged_at_rebuild &&
        fixture
                ->accepted_control
                ->control ==
            fixture
                ->expected_entry_control &&
        fixture
                ->accepted_control
                ->bottom_hole_pressure_pa ==
            fixture
                ->expected_entry_bhp_pa &&
        fixture
                ->clock
                ->accepted_time_seconds() ==
            0.0 &&
        fixture
                ->clock
                ->accepted_step_count() ==
            0U;

    ++fixture->rebuild_calls;
    auto cells =
        make_controller_cells(
            fixture->rank,
            true,
            fixture->audit);
    auto faces =
        make_face_inputs(
            fixture->rank,
            false);

    PetscErrorCode error =
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_WORLD,
                *fixture->schedule,
                *fixture->partition,
                *fixture->bridge,
                *fixture->pattern,
                current_system
                    .time_step_seconds(),
                std::move(cells),
                std::move(faces),
                {
                    {&evaluate_1p,
                     fixture->audit},
                    {&evaluate_2p,
                     fixture->audit},
                    {&evaluate_3p,
                     fixture->audit}},
                wdp::
                    fixed_total_molar_rate_well_source_binding_3d(
                        fixture
                            ->source_context),
                &fdp::
                    evaluate_absent_phase_thermodynamic_provider_3d<
                        flow::
                            Pr76AbsentPhasePotentialExtensionProvider<
                                double>>,
                fixture->provider,
                output);
    if (error != PETSC_SUCCESS ||
        *output == nullptr) {
        return error != PETSC_SUCCESS
            ? error
            : PETSC_ERR_PLIB;
    }

    const auto local_previous =
        owned_conserved_totals(
            **output,
            (*output)->initial_state());
    if (MPI_Allreduce(
            local_previous.data(),
            fixture
                ->rebuilt_previous_totals
                .data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) !=
        MPI_SUCCESS) {
        output->reset();
        return PETSC_ERR_MPI;
    }
    fixture
        ->rebuilt_previous_totals_captured =
        true;
    return PETSC_SUCCESS;
}

void run_fixed_total_molar_rate_control_case(
    int rank,
    const dp::
        ParallelOwnedConnectionSchedule3D&
            schedule,
    const mesh::PartitionSnapshot&
        partition,
    const dp::
        PetscMpiAijSymbolicPreallocation3D&
            bridge,
    const dp::
        OwnedCellStructuralColumnPatternSnapshot3D&
            pattern,
    DispatchAudit* audit) {
    if (audit == nullptr) {
        throw std::invalid_argument(
            "invalid fixed-total-molar-rate control fixture");
    }

    audit->well_timestep_compressibility =
        true;

    constexpr double initial_bhp_pa = 5.0;
    constexpr double target_reference_bhp_pa = 4.0;

    auto connection30 =
        wdp::
            FixedBhpPeacemanWellSourceEvaluatorContext3D::
                create(
                    mesh::GlobalEntityId{
                        UINT64_C(30)},
                    well::make_peaceman_well_index_3d(
                        {10.0, 10.0, 5.0},
                        {
                            1.0e-8,
                            1.0e-8,
                            1.0e-8},
                        well::
                            AxisAlignedWellDirection3D::z,
                        0.10,
                        0.0),
                    initial_bhp_pa,
                    wd::FixedBhpInjectionEnthalpy3D{
                        "fixture/rate-control-cell30-enthalpy/v1",
                        {1000.0, 2000.0}},
                    "fixture/rate-control-cell30-source/v1");
    auto connection60 =
        wdp::
            FixedBhpPeacemanWellSourceEvaluatorContext3D::
                create(
                    mesh::GlobalEntityId{
                        UINT64_C(60)},
                    well::make_peaceman_well_index_3d(
                        {12.0, 9.0, 4.0},
                        {
                            2.0e-8,
                            1.5e-8,
                            2.5e-8},
                        well::
                            AxisAlignedWellDirection3D::z,
                        0.12,
                        0.1),
                    initial_bhp_pa,
                    wd::FixedBhpInjectionEnthalpy3D{
                        "fixture/rate-control-cell60-enthalpy/v1",
                        {1100.0, 2200.0, 3300.0}},
                    "fixture/rate-control-cell60-source/v1");

    auto multi_context =
        wdp::
            FixedBhpMultiConnectionWellSourceEvaluatorContext3D::
                create(
                    "fixture/rate-controlled-two-connection-well/v1",
                    {
                        connection60,
                        connection30});

    auto initial30 =
        evaluate_target(
            UINT64_C(30),
            audit);
    auto initial60 =
        evaluate_target(
            UINT64_C(60),
            audit);

    const auto target_source30 =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                multi_context
                    .find_connection(
                        mesh::GlobalEntityId{
                            UINT64_C(30)})
                    ->with_bottom_hole_pressure(
                        target_reference_bhp_pa),
                initial30);
    const auto target_source60 =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                multi_context
                    .find_connection(
                        mesh::GlobalEntityId{
                            UINT64_C(60)})
                    ->with_bottom_hole_pressure(
                        target_reference_bhp_pa),
                initial60);

    const auto initial_source30 =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                *multi_context.find_connection(
                    mesh::GlobalEntityId{
                        UINT64_C(30)}),
                initial30);
    const auto initial_source60 =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                *multi_context.find_connection(
                    mesh::GlobalEntityId{
                        UINT64_C(60)}),
                initial60);

    auto total_molar_production =
        [](const auto& first,
           const auto& second) {
            double total = 0.0;
            for (double value :
                 first
                     .cell_source
                     .component_molar_rate_mol_per_s) {
                total -= value;
            }
            for (double value :
                 second
                     .cell_source
                     .component_molar_rate_mol_per_s) {
                total -= value;
            }
            return total;
        };

    const double target_rate =
        total_molar_production(
            target_source30,
            target_source60);
    const double initial_rate =
        total_molar_production(
            initial_source30,
            initial_source60);
    require_collective(
        std::isfinite(target_rate) &&
            target_rate > 0.0 &&
            std::isfinite(initial_rate) &&
            initial_rate > 0.0 &&
            std::abs(
                target_rate -
                initial_rate) >
                1.0e-8 *
                    std::max(
                        target_rate,
                        initial_rate),
        "fixed-total-molar-rate target does not require a nontrivial BHP correction");

    auto source_context =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    multi_context,
                    initial_bhp_pa);

    const auto pr_parameters =
        thermodynamic_adapter_pr_parameters();
    const auto pr_model =
        th::Pr76Phase<double>::
            from_parameters(
                pr_parameters);
    auto provider =
        make_outer_rebuild_pr_provider(
            pr_model);

    auto reservoir_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &source_context,
            &provider);

    require_collective(
        reservoir_system != nullptr &&
            reservoir_system
                    ->numbering()
                    .petsc_global_scalar_count() ==
                42 &&
            reservoir_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                2U &&
            reservoir_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{5U})
                    .phase_count ==
                3U,
        "rate-control reservoir system changed frozen phase cardinality");

    const auto local_previous_total =
        owned_conserved_totals(
            *reservoir_system,
            reservoir_system
                ->initial_state());
    std::array<double, 4>
        global_previous_total{};
    require_collective(
        MPI_Allreduce(
            local_previous_total.data(),
            global_previous_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce rate-control previous conserved totals");

    std::unique_ptr<
        wdp::
            FixedTotalMolarRateWellControlSystem3D>
        control_system;
    PetscErrorCode error =
        wdp::
            FixedTotalMolarRateWellControlSystem3D::
                create(
                    PETSC_COMM_WORLD,
                    reservoir_system.get(),
                    &source_context,
                    initial_bhp_pa,
                    target_rate,
                    &control_system);
    require_collective(
        error == PETSC_SUCCESS &&
            control_system != nullptr &&
            control_system
                    ->reservoir_global_scalar_count() ==
                42 &&
            control_system
                    ->global_scalar_count() ==
                43 &&
            control_system
                    ->well_global_scalar() ==
                42 &&
            control_system
                    ->well_owner_rank() ==
                1,
        "rate-control augmented numbering did not add exactly one authoritative BHP scalar");

    Vec solution = nullptr;
    std::optional<
        wdp::
            FixedTotalMolarRateWellControlSolveReport3D>
        report;
    error =
        control_system->solve(
            &solution,
            &report);

    require_collective(
        error == PETSC_SUCCESS &&
            solution != nullptr &&
            report.has_value() &&
            report->converged() &&
            report->global_scalar_count ==
                43 &&
            report->well_global_scalar ==
                42 &&
            report->well_owner_rank ==
                1 &&
            report->snes_type ==
                SNESNEWTONLS &&
            report->ksp_type ==
                KSPGMRES &&
            report->pc_type ==
                PCASM &&
            std::isfinite(
                report
                    ->bottom_hole_pressure_pa) &&
            report->bottom_hole_pressure_pa >
                0.0 &&
            std::abs(
                report
                    ->bottom_hole_pressure_pa -
                initial_bhp_pa) >
                1.0e-7 &&
            std::isfinite(
                report
                    ->achieved_total_molar_rate_mol_per_s) &&
            std::abs(
                report
                    ->total_molar_rate_residual_mol_per_s()) <=
                1.0e-8 *
                    std::max(
                        1.0,
                        std::abs(
                            target_rate)),
        "fixed-total-molar-rate SNES did not converge to a nontrivial BHP/control root");

    Mat controlled_jacobian = nullptr;
    error =
        MatDuplicate(
            control_system
                ->jacobian_structure(),
            MAT_DO_NOT_COPY_VALUES,
            &controlled_jacobian);
    fdp::NaturalVariableSnesEvaluationStatus3D
        controlled_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    auto controlled_evaluator =
        control_system
            ->snes_evaluator();
    if (error == PETSC_SUCCESS) {
        error =
            controlled_evaluator.jacobian(
                solution,
                controlled_jacobian,
                controlled_evaluator
                    .user_context,
                &controlled_status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyBegin(
                controlled_jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                controlled_jacobian,
                MAT_FINAL_ASSEMBLY);
    }

    bool local_four_blocks_ok =
        error == PETSC_SUCCESS &&
        controlled_status ==
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    const PetscInt well_row =
        control_system
            ->well_global_scalar();

    if (local_four_blocks_ok &&
        rank == 0) {
        const auto& cell30 =
            reservoir_system
                ->numbering()
                .cell(
                    mesh::LocalIndex{2U});
        const PetscInt row =
            cell30
                .petsc_global_scalar_start;
        const PetscInt column =
            row;
        PetscScalar jrr = 0.0;
        PetscScalar jrw = 0.0;
        if (MatGetValues(
                controlled_jacobian,
                1,
                &row,
                1,
                &column,
                &jrr) !=
                PETSC_SUCCESS ||
            MatGetValues(
                controlled_jacobian,
                1,
                &row,
                1,
                &well_row,
                &jrw) !=
                PETSC_SUCCESS ||
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        jrr))) <=
                0.0 ||
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        jrw))) <=
                0.0) {
            local_four_blocks_ok =
                false;
        }
    }

    if (local_four_blocks_ok &&
        rank == 1) {
        const auto& cell30 =
            reservoir_system
                ->numbering()
                .cell(
                    mesh::LocalIndex{2U});
        const auto& cell60 =
            reservoir_system
                ->numbering()
                .cell(
                    mesh::LocalIndex{5U});
        const PetscInt cell30_column =
            cell30
                .petsc_global_scalar_start;
        const PetscInt cell60_column =
            cell60
                .petsc_global_scalar_start;
        const PetscInt cell60_row =
            cell60
                .petsc_global_scalar_start;
        PetscScalar jrr = 0.0;
        PetscScalar jrw = 0.0;
        PetscScalar jwr30 = 0.0;
        PetscScalar jwr60 = 0.0;
        PetscScalar jww = 0.0;
        if (MatGetValues(
                controlled_jacobian,
                1,
                &cell60_row,
                1,
                &cell60_column,
                &jrr) !=
                PETSC_SUCCESS ||
            MatGetValues(
                controlled_jacobian,
                1,
                &cell60_row,
                1,
                &well_row,
                &jrw) !=
                PETSC_SUCCESS ||
            MatGetValues(
                controlled_jacobian,
                1,
                &well_row,
                1,
                &cell30_column,
                &jwr30) !=
                PETSC_SUCCESS ||
            MatGetValues(
                controlled_jacobian,
                1,
                &well_row,
                1,
                &cell60_column,
                &jwr60) !=
                PETSC_SUCCESS ||
            MatGetValues(
                controlled_jacobian,
                1,
                &well_row,
                1,
                &well_row,
                &jww) !=
                PETSC_SUCCESS ||
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        jrr))) <=
                0.0 ||
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        jrw))) <=
                0.0 ||
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        jwr30))) <=
                0.0 ||
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        jwr60))) <=
                0.0 ||
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        jww))) <=
                0.0) {
            local_four_blocks_ok =
                false;
        }
    }
    require_collective(
        local_four_blocks_ok,
        "rate-control Jacobian did not materialize analytic Jrr/Jrw/Jwr/Jww blocks");

    Vec reservoir_solution = nullptr;
    error =
        control_system
            ->copy_reservoir_state(
                solution,
                &reservoir_solution);
    require_collective(
        error == PETSC_SUCCESS &&
            reservoir_solution != nullptr,
        "failed to extract reservoir state from rate-controlled solution");

    const auto local_final_total =
        owned_conserved_totals(
            *reservoir_system,
            reservoir_solution);
    std::array<double, 4>
        global_final_total{};
    require_collective(
        MPI_Allreduce(
            local_final_total.data(),
            global_final_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce rate-control final conserved totals");

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        final_current;
    std::vector<double>
        final_porosity;
    fdp::NaturalVariableSnesEvaluationStatus3D
        final_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    require_collective(
        reservoir_system
                ->evaluate_local_cells_for_phase_transition(
                    reservoir_solution,
                    &final_current,
                    &final_porosity,
                    &final_status) ==
            PETSC_SUCCESS &&
            final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            final_current.size() > 5U &&
            final_current[2U]
                .has_value() &&
            final_current[5U]
                .has_value(),
        "failed to evaluate rate-controlled reservoir solution");

    std::array<double, 4>
        local_well_rate{};
    std::uint64_t
        local_authoritative_count = 0U;

    const auto& record30 =
        reservoir_system
            ->numbering()
            .cell(
                mesh::LocalIndex{2U});
    if (record30.owner_rank ==
        reservoir_system
            ->numbering()
            .local_rank()) {
        const auto source =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    multi_context
                        .find_connection(
                            record30
                                .cell_global)
                        ->with_bottom_hole_pressure(
                            report
                                ->bottom_hole_pressure_pa),
                    *final_current[2U]);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            local_well_rate[
                component] -=
                source
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component];
        }
        local_well_rate[3] -=
            source
                .cell_source
                .energy_rate_w;
        ++local_authoritative_count;
    }

    const auto& record60 =
        reservoir_system
            ->numbering()
            .cell(
                mesh::LocalIndex{5U});
    if (record60.owner_rank ==
        reservoir_system
            ->numbering()
            .local_rank()) {
        const auto source =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    multi_context
                        .find_connection(
                            record60
                                .cell_global)
                        ->with_bottom_hole_pressure(
                            report
                                ->bottom_hole_pressure_pa),
                    *final_current[5U]);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            local_well_rate[
                component] -=
                source
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component];
        }
        local_well_rate[3] -=
            source
                .cell_source
                .energy_rate_w;
        ++local_authoritative_count;
    }

    std::array<double, 4>
        global_well_rate{};
    std::uint64_t
        global_authoritative_count = 0U;
    require_collective(
        MPI_Allreduce(
            local_well_rate.data(),
            global_well_rate.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            MPI_Allreduce(
                &local_authoritative_count,
                &global_authoritative_count,
                1,
                MPI_UINT64_T,
                MPI_SUM,
                PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            global_authoritative_count ==
                2U,
        "rate-controlled well did not retain two unique authoritative connections");

    const double summed_molar_rate =
        global_well_rate[0] +
        global_well_rate[1] +
        global_well_rate[2];
    near_collective(
        summed_molar_rate,
        target_rate,
        1.0e-8,
        1.0e-10);

    const double dt =
        reservoir_system
            ->time_step_seconds();
    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        near_collective(
            global_final_total[
                quantity],
            global_previous_total[
                quantity] -
                dt *
                    global_well_rate[
                        quantity],
            3.0e-7,
            quantity < 3U
                ? 3.0e-8
                : 3.0e-7);
    }

    // Re-run the same physical step through the production rate-control
    // timestep bridge. The earlier direct solve remains an independent
    // algebra/Jacobian/conservation oracle; only this bridge is allowed to
    // advance accepted reservoir history/time and the accepted BHP initial
    // guess.
    source_context.begin_evaluation(
        initial_bhp_pa);
    double accepted_bhp_pa =
        initial_bhp_pa;
    fdp::AcceptedPhysicalTimeClock3D
        rate_control_clock{
            0.0,
            reservoir_system
                ->time_step_seconds()};
    wdp::
        FixedTotalMolarRatePhysicalTimestepDriverOptions3D
        timestep_options;
    timestep_options.adaptive
        .minimum_timestep_seconds =
        0.25;
    timestep_options.adaptive
        .maximum_timestep_seconds =
        4.0;
    timestep_options.adaptive
        .maximum_retries =
        2U;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        timestep_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                reservoir_system.get(),
                &source_context,
                target_rate,
                timestep_options,
                &rate_control_clock,
                &accepted_bhp_pa,
                &timestep_report);

    require_collective(
        error == PETSC_SUCCESS &&
            timestep_report.has_value() &&
            timestep_report->accepted() &&
            timestep_report
                    ->accepted_record
                    ->accepted_step_index ==
                0U &&
            timestep_report
                    ->accepted_record
                    ->time_n_seconds ==
                0.0 &&
            timestep_report
                    ->accepted_record
                    ->time_np1_seconds ==
                1.0 &&
            timestep_report
                    ->accepted_record
                    ->accepted_timestep_seconds ==
                1.0 &&
            timestep_report
                    ->accepted_solve
                    ->global_scalar_count ==
                43 &&
            timestep_report
                    ->accepted_solve
                    ->well_global_scalar ==
                42 &&
            timestep_report
                    ->accepted_solve
                    ->well_owner_rank ==
                1 &&
            std::isfinite(
                accepted_bhp_pa) &&
            accepted_bhp_pa > 0.0 &&
            std::abs(
                accepted_bhp_pa -
                initial_bhp_pa) >
                1.0e-7 &&
            std::abs(
                timestep_report
                    ->accepted_solve
                    ->total_molar_rate_residual_mol_per_s()) <=
                1.0e-8 *
                    std::max(
                        1.0,
                        std::abs(
                            target_rate)) &&
            rate_control_clock
                    .accepted_time_seconds() ==
                1.0 &&
            rate_control_clock
                    .accepted_step_count() ==
                1U &&
            source_context
                    .current_bottom_hole_pressure_pa() ==
                accepted_bhp_pa,
        "fixed-total-molar-rate physical-timestep driver did not commit reservoir time/history and BHP exactly once");

    near_collective(
        accepted_bhp_pa,
        report->bottom_hole_pressure_pa,
        2.0e-8,
        2.0e-10);

    Vec committed_difference = nullptr;
    error =
        VecDuplicate(
            reservoir_system
                ->initial_state(),
            &committed_difference);
    if (error == PETSC_SUCCESS) {
        error =
            VecCopy(
                reservoir_system
                    ->initial_state(),
                committed_difference);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAXPY(
                committed_difference,
                PetscScalar{-1.0},
                reservoir_solution);
    }
    PetscReal committed_difference_norm =
        -1.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                committed_difference,
                NORM_2,
                &committed_difference_norm);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            committed_difference_norm <=
                2.0e-8,
        "rate-control timestep driver committed a different reservoir root from the independently verified augmented solve");

    bool history_matches = false;
    require_collective(
        reservoir_system
                ->accepted_history_matches_state(
                    reservoir_system
                        ->initial_state(),
                    &history_matches) ==
            PETSC_SUCCESS &&
            history_matches,
        "rate-controlled accepted reservoir history did not rebase through the physical-timestep driver");

    Vec accepted_snapshot = nullptr;
    error =
        VecDuplicate(
            reservoir_system
                ->initial_state(),
            &accepted_snapshot);
    if (error == PETSC_SUCCESS) {
        error =
            VecCopy(
                reservoir_system
                    ->initial_state(),
                accepted_snapshot);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            accepted_snapshot != nullptr,
        "failed to snapshot accepted rate-control state for rollback regression");

    const double accepted_time_before_rejection =
        rate_control_clock
            .accepted_time_seconds();
    const std::size_t
        accepted_steps_before_rejection =
            rate_control_clock
                .accepted_step_count();
    const double
        next_dt_before_rejection =
            rate_control_clock
                .next_timestep_seconds();
    const double
        bhp_before_rejection =
            accepted_bhp_pa;

    auto invalid_timestep_options =
        timestep_options;
    invalid_timestep_options.adaptive
        .minimum_timestep_seconds =
        next_dt_before_rejection *
        2.0;
    invalid_timestep_options.adaptive
        .maximum_timestep_seconds =
        next_dt_before_rejection *
        4.0;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        rejected_timestep_report;
    const PetscErrorCode rejection_error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                reservoir_system.get(),
                &source_context,
                target_rate,
                invalid_timestep_options,
                &rate_control_clock,
                &accepted_bhp_pa,
                &rejected_timestep_report);

    Vec rollback_difference = nullptr;
    error =
        VecDuplicate(
            reservoir_system
                ->initial_state(),
            &rollback_difference);
    if (error == PETSC_SUCCESS) {
        error =
            VecCopy(
                reservoir_system
                    ->initial_state(),
                rollback_difference);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAXPY(
                rollback_difference,
                PetscScalar{-1.0},
                accepted_snapshot);
    }
    PetscReal rollback_difference_norm =
        -1.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                rollback_difference,
                NORM_2,
                &rollback_difference_norm);
    }
    require_collective(
        rejection_error ==
                PETSC_ERR_ARG_OUTOFRANGE &&
            !rejected_timestep_report
                 .has_value() &&
            error == PETSC_SUCCESS &&
            rollback_difference_norm <=
                1.0e-14 &&
            rate_control_clock
                    .accepted_time_seconds() ==
                accepted_time_before_rejection &&
            rate_control_clock
                    .accepted_step_count() ==
                accepted_steps_before_rejection &&
            rate_control_clock
                    .next_timestep_seconds() ==
                next_dt_before_rejection &&
            reservoir_system
                    ->time_step_seconds() ==
                next_dt_before_rejection &&
            accepted_bhp_pa ==
                bhp_before_rejection &&
            source_context
                    .current_bottom_hole_pressure_pa() ==
                bhp_before_rejection,
        "rejected rate-control timestep entry mutated accepted state/time/dt/BHP ownership");

    require_collective(
        VecDestroy(
            &rollback_difference) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &accepted_snapshot) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &committed_difference) ==
                PETSC_SUCCESS,
        "rate-control timestep bridge regression cleanup failed");

    // Exercise a real nonlinear retry rather than a test hook.  With a fixed
    // positive whole-well production target, an enormous physical timestep
    // would require negative global molar inventory.  The frozen natural-
    // variable domains cannot represent that state, so the augmented SNES
    // must reject the first trial.  A single cutback returns exactly to the
    // independently validated 1 s problem above.
    auto cutback_source_context =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    multi_context,
                    initial_bhp_pa);
    auto cutback_reservoir_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &cutback_source_context,
            &provider);

    constexpr double impossible_first_dt_seconds =
        1.0e9;
    require_collective(
        cutback_reservoir_system != nullptr &&
            cutback_reservoir_system
                    ->set_trial_timestep_seconds(
                        impossible_first_dt_seconds) ==
                PETSC_SUCCESS,
        "failed to prepare rate-control cutback baseline");

    fdp::AcceptedPhysicalTimeClock3D
        cutback_clock{
            0.0,
            impossible_first_dt_seconds};
    wdp::
        FixedTotalMolarRatePhysicalTimestepDriverOptions3D
        cutback_options;
    cutback_options.adaptive
        .minimum_timestep_seconds =
        1.0;
    cutback_options.adaptive
        .maximum_timestep_seconds =
        impossible_first_dt_seconds;
    cutback_options.adaptive
        .cutback_factor =
        1.0e-9;
    cutback_options.adaptive
        .growth_factor =
        2.0;
    cutback_options.adaptive
        .maximum_retries =
        1U;
    cutback_options.adaptive
        .growth_nonlinear_iteration_limit =
        30;
    cutback_options.adaptive
        .growth_line_search_direction_change_limit =
        30;

    double cutback_accepted_bhp_pa =
        initial_bhp_pa;
    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        cutback_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                cutback_reservoir_system.get(),
                &cutback_source_context,
                target_rate,
                cutback_options,
                &cutback_clock,
                &cutback_accepted_bhp_pa,
                &cutback_report);

    require_collective(
        error == PETSC_SUCCESS &&
            cutback_report.has_value() &&
            cutback_report->accepted() &&
            cutback_report
                    ->adaptive
                    .attempts
                    .size() ==
                2U &&
            cutback_report
                    ->adaptive
                    .retries ==
                1U &&
            cutback_report
                    ->adaptive
                    .accepted_timestep_seconds
                    .has_value() &&
            *cutback_report
                  ->adaptive
                  .accepted_timestep_seconds ==
                1.0 &&
            cutback_report
                    ->accepted_record
                    ->accepted_timestep_seconds ==
                1.0 &&
            cutback_clock
                    .accepted_time_seconds() ==
                1.0 &&
            cutback_clock
                    .accepted_step_count() ==
                1U,
        "rate-control physical-timestep driver did not cut back one real failed augmented solve to the accepted 1 s retry");

    const auto& failed_attempt =
        cutback_report
            ->adaptive
            .attempts[0U];
    const auto& accepted_retry =
        cutback_report
            ->adaptive
            .attempts[1U];
    require_collective(
        failed_attempt
                .request
                .attempt_index ==
            0U &&
            failed_attempt
                .request
                .retry_index ==
            0U &&
            failed_attempt
                .request
                .timestep_seconds ==
            impossible_first_dt_seconds &&
            failed_attempt
                .result
                .outcome ==
            fdp::
                AdaptiveTimestepAttemptOutcome3D::
                    nonlinear_solve_diverged &&
            failed_attempt
                .decision ==
            fdp::
                AdaptiveTimestepDecision3D::
                    reject_and_cutback &&
            failed_attempt
                .next_attempt_timestep_seconds
                .has_value() &&
            *failed_attempt
                  .next_attempt_timestep_seconds ==
                1.0 &&
            accepted_retry
                .request
                .attempt_index ==
            1U &&
            accepted_retry
                .request
                .retry_index ==
            1U &&
            accepted_retry
                .request
                .timestep_seconds ==
            1.0 &&
            accepted_retry
                .result
                .outcome ==
            fdp::
                AdaptiveTimestepAttemptOutcome3D::
                    stable_phase_set,
        "rate-control adaptive report did not preserve failed-first-attempt/cutback/retry provenance");

    near_collective(
        cutback_accepted_bhp_pa,
        report->bottom_hole_pressure_pa,
        2.0e-8,
        2.0e-10);

    Vec cutback_root_difference = nullptr;
    error =
        VecDuplicate(
            cutback_reservoir_system
                ->initial_state(),
            &cutback_root_difference);
    if (error == PETSC_SUCCESS) {
        error =
            VecCopy(
                cutback_reservoir_system
                    ->initial_state(),
                cutback_root_difference);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAXPY(
                cutback_root_difference,
                PetscScalar{-1.0},
                reservoir_solution);
    }
    PetscReal cutback_root_difference_norm =
        -1.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                cutback_root_difference,
                NORM_2,
                &cutback_root_difference_norm);
    }
    bool cutback_history_matches = false;
    if (error == PETSC_SUCCESS) {
        error =
            cutback_reservoir_system
                ->accepted_history_matches_state(
                    cutback_reservoir_system
                        ->initial_state(),
                    &cutback_history_matches);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            cutback_root_difference_norm <=
                2.0e-8 &&
            cutback_history_matches &&
            cutback_source_context
                    .current_bottom_hole_pressure_pa() ==
                cutback_accepted_bhp_pa &&
            std::abs(
                cutback_report
                    ->accepted_solve
                    ->total_molar_rate_residual_mol_per_s()) <=
                1.0e-8 *
                    std::max(
                        1.0,
                        std::abs(
                            target_rate)),
        "rate-control retry reused failed-trial state/BHP instead of rebuilding the validated 1 s accepted root");

    require_collective(
        VecDestroy(
            &cutback_root_difference) ==
            PETSC_SUCCESS,
        "rate-control nonlinear cutback regression cleanup failed");

    // Bind a true producer minimum-BHP constraint between the unconstrained
    // rate-control root and the entry BHP. The fixed-rate candidate must be
    // discarded, then the same 1 s physical timestep must be re-solved from
    // the original accepted reservoir history at exactly p_min.
    require_collective(
        report->bottom_hole_pressure_pa <
            initial_bhp_pa,
        "rate-control switching fixture did not produce a lower-BHP unconstrained rate root");
    const double minimum_bhp_pa =
        0.5 *
        (report->bottom_hole_pressure_pa +
         initial_bhp_pa);
    require_collective(
        minimum_bhp_pa >
                report
                    ->bottom_hole_pressure_pa &&
            minimum_bhp_pa <
                initial_bhp_pa,
        "minimum-BHP switching fixture did not place p_min strictly above the rate-control candidate");

    auto switch_source_context =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    multi_context,
                    initial_bhp_pa);
    auto switch_reservoir_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &switch_source_context,
            &provider);

    const auto switch_local_previous_total =
        owned_conserved_totals(
            *switch_reservoir_system,
            switch_reservoir_system
                ->initial_state());
    std::array<double, 4>
        switch_global_previous_total{};
    require_collective(
        MPI_Allreduce(
            switch_local_previous_total.data(),
            switch_global_previous_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce minimum-BHP switching previous conserved totals");

    fdp::AcceptedPhysicalTimeClock3D
        switch_clock{
            0.0,
            1.0};
    wdp::
        FixedTotalMolarRatePhysicalTimestepDriverOptions3D
        switch_options;
    switch_options.adaptive
        .minimum_timestep_seconds =
        0.25;
    switch_options.adaptive
        .maximum_timestep_seconds =
        4.0;
    switch_options.adaptive
        .maximum_retries =
        2U;
    switch_options.minimum_bottom_hole_pressure_pa =
        minimum_bhp_pa;

    auto switch_control_state =
        wdp::
            AcceptedFixedTotalMolarRateWellControlState3D::
                fixed_total_molar_rate(
                    initial_bhp_pa);
    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        switch_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                switch_reservoir_system.get(),
                &switch_source_context,
                target_rate,
                switch_options,
                &switch_clock,
                &switch_control_state,
                &switch_report);

    require_collective(
        error == PETSC_SUCCESS &&
            switch_report.has_value() &&
            switch_report->accepted() &&
            switch_report
                    ->rate_to_bhp_switch_triggered &&
            switch_report
                    ->accepted_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            !switch_report
                 ->accepted_solve
                 .has_value() &&
            switch_report
                    ->accepted_fixed_bhp_solve
                    .has_value() &&
            switch_report
                    ->discarded_rate_control_candidate
                    .has_value() &&
            switch_report
                    ->discarded_rate_control_candidate
                    ->bottom_hole_pressure_pa <
                minimum_bhp_pa &&
            switch_report
                    ->accepted_record
                    ->accepted_timestep_seconds ==
                1.0 &&
            switch_clock
                    .accepted_time_seconds() ==
                1.0 &&
            switch_clock
                    .accepted_step_count() ==
                1U &&
            switch_control_state
                    .control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            switch_control_state
                    .bottom_hole_pressure_pa ==
                minimum_bhp_pa &&
            switch_source_context
                    .current_bottom_hole_pressure_pa() ==
                minimum_bhp_pa,
        "minimum-BHP constraint did not discard the violating rate candidate and accept a same-timestep fixed-BHP re-solve");

    near_collective(
        switch_report
            ->discarded_rate_control_candidate
            ->bottom_hole_pressure_pa,
        report->bottom_hole_pressure_pa,
        2.0e-8,
        2.0e-10);

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        switch_final_current;
    std::vector<double>
        switch_final_porosity;
    fdp::NaturalVariableSnesEvaluationStatus3D
        switch_final_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    require_collective(
        switch_reservoir_system
                ->evaluate_local_cells_for_phase_transition(
                    switch_reservoir_system
                        ->initial_state(),
                    &switch_final_current,
                    &switch_final_porosity,
                    &switch_final_status) ==
            PETSC_SUCCESS &&
            switch_final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            switch_final_current.size() > 5U &&
            switch_final_current[2U]
                .has_value() &&
            switch_final_current[5U]
                .has_value(),
        "failed to evaluate accepted minimum-BHP reservoir state");

    std::array<double, 4>
        switch_local_well_rate{};
    std::uint64_t
        switch_local_authoritative_count = 0U;
    const auto& switch_record30 =
        switch_reservoir_system
            ->numbering()
            .cell(
                mesh::LocalIndex{2U});
    if (switch_record30.owner_rank ==
        switch_reservoir_system
            ->numbering()
            .local_rank()) {
        const auto source =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    multi_context
                        .find_connection(
                            switch_record30
                                .cell_global)
                        ->with_bottom_hole_pressure(
                            minimum_bhp_pa),
                    *switch_final_current[2U]);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            switch_local_well_rate[
                component] -=
                source
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component];
        }
        switch_local_well_rate[3] -=
            source
                .cell_source
                .energy_rate_w;
        ++switch_local_authoritative_count;
    }

    const auto& switch_record60 =
        switch_reservoir_system
            ->numbering()
            .cell(
                mesh::LocalIndex{5U});
    if (switch_record60.owner_rank ==
        switch_reservoir_system
            ->numbering()
            .local_rank()) {
        const auto source =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    multi_context
                        .find_connection(
                            switch_record60
                                .cell_global)
                        ->with_bottom_hole_pressure(
                            minimum_bhp_pa),
                    *switch_final_current[5U]);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            switch_local_well_rate[
                component] -=
                source
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component];
        }
        switch_local_well_rate[3] -=
            source
                .cell_source
                .energy_rate_w;
        ++switch_local_authoritative_count;
    }

    std::array<double, 4>
        switch_global_well_rate{};
    std::uint64_t
        switch_global_authoritative_count = 0U;
    require_collective(
        MPI_Allreduce(
            switch_local_well_rate.data(),
            switch_global_well_rate.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            MPI_Allreduce(
                &switch_local_authoritative_count,
                &switch_global_authoritative_count,
                1,
                MPI_UINT64_T,
                MPI_SUM,
                PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            switch_global_authoritative_count ==
                2U,
        "minimum-BHP switched well lost owner-only multi-connection aggregation");

    const double switch_total_molar_rate =
        switch_global_well_rate[0] +
        switch_global_well_rate[1] +
        switch_global_well_rate[2];
    require_collective(
        std::isfinite(
            switch_total_molar_rate) &&
            switch_total_molar_rate <
                target_rate -
                    1.0e-10 *
                        std::max(
                            1.0,
                            std::abs(
                                target_rate)),
        "minimum-BHP constrained solution still forced the infeasible fixed-rate target");

    const auto switch_local_final_total =
        owned_conserved_totals(
            *switch_reservoir_system,
            switch_reservoir_system
                ->initial_state());
    std::array<double, 4>
        switch_global_final_total{};
    require_collective(
        MPI_Allreduce(
            switch_local_final_total.data(),
            switch_global_final_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce minimum-BHP switching final conserved totals");

    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        near_collective(
            switch_global_final_total[
                quantity],
            switch_global_previous_total[
                quantity] -
                switch_global_well_rate[
                    quantity],
            3.0e-7,
            quantity < 3U
                ? 3.0e-8
                : 3.0e-7);
    }

    bool switch_history_matches = false;
    require_collective(
        switch_reservoir_system
                ->accepted_history_matches_state(
                    switch_reservoir_system
                        ->initial_state(),
                    &switch_history_matches) ==
                PETSC_SUCCESS &&
            switch_history_matches,
        "minimum-BHP switch did not commit exactly the final fixed-BHP reservoir state/history");

    // Advance a second accepted physical timestep with the same accepted
    // well-control state. It must enter directly in minimum-BHP mode: no
    // augmented rate solve, no discarded rate candidate and no new rate->BHP
    // transition are allowed.
    const double second_step_dt =
        switch_clock
            .next_timestep_seconds();
    const double second_step_time_n =
        switch_clock
            .accepted_time_seconds();
    const auto second_step_local_previous_total =
        owned_conserved_totals(
            *switch_reservoir_system,
            switch_reservoir_system
                ->initial_state());
    std::array<double, 4>
        second_step_global_previous_total{};
    require_collective(
        MPI_Allreduce(
            second_step_local_previous_total.data(),
            second_step_global_previous_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce persisted minimum-BHP second-step previous totals");

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        persisted_bhp_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                switch_reservoir_system.get(),
                &switch_source_context,
                target_rate,
                switch_options,
                &switch_clock,
                &switch_control_state,
                &persisted_bhp_report);

    require_collective(
        error == PETSC_SUCCESS &&
            persisted_bhp_report.has_value() &&
            persisted_bhp_report->accepted() &&
            persisted_bhp_report
                    ->entry_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            persisted_bhp_report
                    ->accepted_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            !persisted_bhp_report
                 ->rate_to_bhp_switch_triggered &&
            !persisted_bhp_report
                 ->discarded_rate_control_candidate
                 .has_value() &&
            !persisted_bhp_report
                 ->accepted_solve
                 .has_value() &&
            persisted_bhp_report
                    ->accepted_fixed_bhp_solve
                    .has_value() &&
            persisted_bhp_report
                    ->accepted_record
                    ->accepted_timestep_seconds ==
                second_step_dt &&
            switch_clock
                    .accepted_time_seconds() ==
                second_step_time_n +
                    second_step_dt &&
            switch_clock
                    .accepted_step_count() ==
                2U &&
            switch_control_state
                    .control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            switch_control_state
                    .bottom_hole_pressure_pa ==
                minimum_bhp_pa &&
            switch_source_context
                    .current_bottom_hole_pressure_pa() ==
                minimum_bhp_pa,
        "accepted minimum-BHP control state did not persist directly into the next physical timestep");

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        second_step_final_current;
    std::vector<double>
        second_step_final_porosity;
    fdp::NaturalVariableSnesEvaluationStatus3D
        second_step_final_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    require_collective(
        switch_reservoir_system
                ->evaluate_local_cells_for_phase_transition(
                    switch_reservoir_system
                        ->initial_state(),
                    &second_step_final_current,
                    &second_step_final_porosity,
                    &second_step_final_status) ==
            PETSC_SUCCESS &&
            second_step_final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            second_step_final_current.size() > 5U &&
            second_step_final_current[2U]
                .has_value() &&
            second_step_final_current[5U]
                .has_value(),
        "failed to evaluate persisted minimum-BHP second-step state");

    std::array<double, 4>
        second_step_local_well_rate{};
    std::uint64_t
        second_step_local_authoritative_count = 0U;
    if (switch_record30.owner_rank ==
        switch_reservoir_system
            ->numbering()
            .local_rank()) {
        const auto source =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    multi_context
                        .find_connection(
                            switch_record30
                                .cell_global)
                        ->with_bottom_hole_pressure(
                            minimum_bhp_pa),
                    *second_step_final_current[2U]);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            second_step_local_well_rate[
                component] -=
                source
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component];
        }
        second_step_local_well_rate[3] -=
            source
                .cell_source
                .energy_rate_w;
        ++second_step_local_authoritative_count;
    }
    if (switch_record60.owner_rank ==
        switch_reservoir_system
            ->numbering()
            .local_rank()) {
        const auto source =
            wdp::
                build_fixed_bhp_peaceman_well_source_3d(
                    multi_context
                        .find_connection(
                            switch_record60
                                .cell_global)
                        ->with_bottom_hole_pressure(
                            minimum_bhp_pa),
                    *second_step_final_current[5U]);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            second_step_local_well_rate[
                component] -=
                source
                    .cell_source
                    .component_molar_rate_mol_per_s[
                        component];
        }
        second_step_local_well_rate[3] -=
            source
                .cell_source
                .energy_rate_w;
        ++second_step_local_authoritative_count;
    }

    std::array<double, 4>
        second_step_global_well_rate{};
    std::uint64_t
        second_step_global_authoritative_count = 0U;
    require_collective(
        MPI_Allreduce(
            second_step_local_well_rate.data(),
            second_step_global_well_rate.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            MPI_Allreduce(
                &second_step_local_authoritative_count,
                &second_step_global_authoritative_count,
                1,
                MPI_UINT64_T,
                MPI_SUM,
                PETSC_COMM_WORLD) ==
                MPI_SUCCESS &&
            second_step_global_authoritative_count ==
                2U,
        "persisted minimum-BHP second step lost owner-only well aggregation");

    const auto second_step_local_final_total =
        owned_conserved_totals(
            *switch_reservoir_system,
            switch_reservoir_system
                ->initial_state());
    std::array<double, 4>
        second_step_global_final_total{};
    require_collective(
        MPI_Allreduce(
            second_step_local_final_total.data(),
            second_step_global_final_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce persisted minimum-BHP second-step final totals");

    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        near_collective(
            second_step_global_final_total[
                quantity],
            second_step_global_previous_total[
                quantity] -
                second_step_dt *
                    second_step_global_well_rate[
                        quantity],
            3.0e-7,
            quantity < 3U
                ? 3.0e-8
                : 3.0e-7);
    }

    // A rejected third entry must not mutate the already accepted BHP control
    // state. This checks control-mode rollback together with the existing
    // reservoir history/time ownership.
    const auto control_before_rejection =
        switch_control_state;
    const double time_before_control_rejection =
        switch_clock
            .accepted_time_seconds();
    const std::size_t steps_before_control_rejection =
        switch_clock
            .accepted_step_count();
    const double dt_before_control_rejection =
        switch_clock
            .next_timestep_seconds();

    auto invalid_persisted_options =
        switch_options;
    invalid_persisted_options.adaptive
        .minimum_timestep_seconds =
        dt_before_control_rejection *
        2.0;
    invalid_persisted_options.adaptive
        .maximum_timestep_seconds =
        dt_before_control_rejection *
        4.0;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        rejected_persisted_report;
    const PetscErrorCode persisted_rejection_error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                switch_reservoir_system.get(),
                &switch_source_context,
                target_rate,
                invalid_persisted_options,
                &switch_clock,
                &switch_control_state,
                &rejected_persisted_report);

    require_collective(
        persisted_rejection_error ==
                PETSC_ERR_ARG_OUTOFRANGE &&
            !rejected_persisted_report
                 .has_value() &&
            switch_control_state.control ==
                control_before_rejection
                    .control &&
            switch_control_state
                    .bottom_hole_pressure_pa ==
                control_before_rejection
                    .bottom_hole_pressure_pa &&
            switch_clock
                    .accepted_time_seconds() ==
                time_before_control_rejection &&
            switch_clock
                    .accepted_step_count() ==
                steps_before_control_rejection &&
            switch_clock
                    .next_timestep_seconds() ==
                dt_before_control_rejection &&
            switch_reservoir_system
                    ->time_step_seconds() ==
                dt_before_control_rejection &&
            switch_source_context
                    .current_bottom_hole_pressure_pa() ==
                minimum_bhp_pa,
        "rejected persisted minimum-BHP timestep mutated accepted control state/time/history");

    // Build an independent accepted minimum-BHP baseline and measure its
    // p_min capacity without committing it. The guarded reactivation target is
    // then chosen strictly below that capacity, so both release gates can be
    // exercised without introducing schedule logic into the production API.
    auto reactivation_source_context =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    multi_context,
                    minimum_bhp_pa);
    auto reactivation_reservoir_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &reactivation_source_context,
            &provider);

    constexpr double reactivation_dt_seconds =
        0.25;
    require_collective(
        reactivation_reservoir_system !=
                nullptr &&
            reactivation_reservoir_system
                    ->set_trial_timestep_seconds(
                        reactivation_dt_seconds) ==
                PETSC_SUCCESS,
        "failed to prepare minimum-BHP reactivation baseline");

    reactivation_source_context
        .begin_fixed_bhp_reservoir_evaluation(
            minimum_bhp_pa);
    Vec reactivation_probe_state = nullptr;
    std::optional<
        fdp::
            VariableCardinalityNaturalVariableSnesSolveReport3D>
        reactivation_probe_solve;
    error =
        reactivation_reservoir_system
            ->solve(
                &reactivation_probe_state,
                &reactivation_probe_solve);
    require_collective(
        error == PETSC_SUCCESS &&
            reactivation_probe_state !=
                nullptr &&
            reactivation_probe_solve
                .has_value(),
        "failed to solve independent minimum-BHP capacity oracle");

    std::uint64_t
        reactivation_probe_authoritative_count = 0U;
    const auto reactivation_probe_rate =
        global_fixed_bhp_well_production_rate(
            *reactivation_reservoir_system,
            reactivation_probe_state,
            multi_context,
            minimum_bhp_pa,
            &reactivation_probe_authoritative_count);
    const double
        reactivation_probe_total_molar_rate =
            reactivation_probe_rate[0] +
            reactivation_probe_rate[1] +
            reactivation_probe_rate[2];
    require_collective(
        reactivation_probe_authoritative_count ==
                2U &&
            std::isfinite(
                reactivation_probe_total_molar_rate) &&
            reactivation_probe_total_molar_rate >
                0.0,
        "minimum-BHP capacity oracle is not a positive two-connection production rate");

    const double reactivation_target_rate =
        0.50 *
        reactivation_probe_total_molar_rate;
    const double reactivation_rate_margin =
        0.10 *
        reactivation_probe_total_molar_rate;
    constexpr double
        reactivation_pressure_margin_pa =
            1.0e-3;

    const auto reactivation_local_previous_total =
        owned_conserved_totals(
            *reactivation_reservoir_system,
            reactivation_reservoir_system
                ->initial_state());
    std::array<double, 4>
        reactivation_global_previous_total{};
    require_collective(
        MPI_Allreduce(
            reactivation_local_previous_total.data(),
            reactivation_global_previous_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce reactivation previous conserved totals");

    auto reactivation_control_state =
        wdp::
            AcceptedFixedTotalMolarRateWellControlState3D::
                minimum_bottom_hole_pressure(
                    minimum_bhp_pa);
    fdp::AcceptedPhysicalTimeClock3D
        reactivation_clock{
            0.0,
            reactivation_dt_seconds};
    wdp::
        FixedTotalMolarRatePhysicalTimestepDriverOptions3D
        reactivation_options;
    reactivation_options.adaptive
        .minimum_timestep_seconds =
        reactivation_dt_seconds;
    reactivation_options.adaptive
        .maximum_timestep_seconds =
        reactivation_dt_seconds;
    reactivation_options.adaptive
        .maximum_retries =
        1U;
    reactivation_options
        .minimum_bottom_hole_pressure_pa =
        minimum_bhp_pa;
    reactivation_options
        .minimum_bhp_release_rate_margin_mol_per_s =
        reactivation_rate_margin;
    reactivation_options
        .minimum_bhp_release_pressure_margin_pa =
        reactivation_pressure_margin_pa;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        reactivation_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                reactivation_reservoir_system.get(),
                &reactivation_source_context,
                reactivation_target_rate,
                reactivation_options,
                &reactivation_clock,
                &reactivation_control_state,
                &reactivation_report);

    require_collective(
        error == PETSC_SUCCESS &&
            reactivation_report.has_value() &&
            reactivation_report->accepted() &&
            reactivation_report
                    ->entry_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            reactivation_report
                    ->accepted_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        fixed_total_molar_rate &&
            reactivation_report
                    ->minimum_bhp_probe_total_molar_rate_mol_per_s
                    .has_value() &&
            reactivation_report
                    ->minimum_bhp_to_rate_reactivation_attempted &&
            reactivation_report
                    ->minimum_bhp_to_rate_reactivation_accepted &&
            reactivation_report
                    ->discarded_minimum_bhp_candidate
                    .has_value() &&
            !reactivation_report
                 ->discarded_rate_reactivation_candidate
                 .has_value() &&
            reactivation_report
                    ->accepted_solve
                    .has_value() &&
            !reactivation_report
                 ->accepted_fixed_bhp_solve
                 .has_value() &&
            reactivation_report
                    ->accepted_solve
                    ->bottom_hole_pressure_pa >=
                minimum_bhp_pa +
                    reactivation_pressure_margin_pa &&
            std::abs(
                reactivation_report
                    ->accepted_solve
                    ->total_molar_rate_residual_mol_per_s()) <=
                1.0e-8 *
                    std::max(
                        1.0,
                        std::abs(
                            reactivation_target_rate)) &&
            reactivation_control_state
                    .control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        fixed_total_molar_rate &&
            reactivation_control_state
                    .bottom_hole_pressure_pa ==
                reactivation_report
                    ->accepted_solve
                    ->bottom_hole_pressure_pa &&
            reactivation_clock
                    .accepted_time_seconds() ==
                reactivation_dt_seconds &&
            reactivation_clock
                    .accepted_step_count() ==
                1U,
        "minimum-BHP feasibility probe did not reactivate fixed-total-molar-rate control through both hysteresis guards");

    near_collective(
        *reactivation_report
             ->minimum_bhp_probe_total_molar_rate_mol_per_s,
        reactivation_probe_total_molar_rate,
        2.0e-8,
        2.0e-10);

    std::uint64_t
        reactivation_final_authoritative_count = 0U;
    const auto reactivation_final_rate =
        global_fixed_bhp_well_production_rate(
            *reactivation_reservoir_system,
            reactivation_reservoir_system
                ->initial_state(),
            multi_context,
            reactivation_control_state
                .bottom_hole_pressure_pa,
            &reactivation_final_authoritative_count);
    const double
        reactivation_final_total_molar_rate =
            reactivation_final_rate[0] +
            reactivation_final_rate[1] +
            reactivation_final_rate[2];
    require_collective(
        reactivation_final_authoritative_count ==
                2U,
        "reactivated rate-control state lost owner-only two-connection aggregation");
    near_collective(
        reactivation_final_total_molar_rate,
        reactivation_target_rate,
        1.0e-8,
        1.0e-10);

    const auto reactivation_local_final_total =
        owned_conserved_totals(
            *reactivation_reservoir_system,
            reactivation_reservoir_system
                ->initial_state());
    std::array<double, 4>
        reactivation_global_final_total{};
    require_collective(
        MPI_Allreduce(
            reactivation_local_final_total.data(),
            reactivation_global_final_total.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce reactivated rate-control final conserved totals");
    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        near_collective(
            reactivation_global_final_total[
                quantity],
            reactivation_global_previous_total[
                quantity] -
                reactivation_dt_seconds *
                    reactivation_final_rate[
                        quantity],
            3.0e-7,
            quantity < 3U
                ? 3.0e-8
                : 3.0e-7);
    }

    Vec reactivation_probe_difference = nullptr;
    error =
        VecDuplicate(
            reactivation_reservoir_system
                ->initial_state(),
            &reactivation_probe_difference);
    if (error == PETSC_SUCCESS) {
        error =
            VecCopy(
                reactivation_reservoir_system
                    ->initial_state(),
                reactivation_probe_difference);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAXPY(
                reactivation_probe_difference,
                PetscScalar{-1.0},
                reactivation_probe_state);
    }
    PetscReal reactivation_probe_difference_norm =
        -1.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                reactivation_probe_difference,
                NORM_2,
                &reactivation_probe_difference_norm);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            reactivation_probe_difference_norm >
                1.0e-9,
        "reactivation committed the discarded minimum-BHP feasibility probe");

    // The accepted reactivated state must persist: the next step enters rate
    // control directly and does not run another fixed-BHP feasibility probe.
    const double reactivated_second_time_n =
        reactivation_clock
            .accepted_time_seconds();
    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        reactivated_second_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                reactivation_reservoir_system.get(),
                &reactivation_source_context,
                reactivation_target_rate,
                reactivation_options,
                &reactivation_clock,
                &reactivation_control_state,
                &reactivated_second_report);
    require_collective(
        error == PETSC_SUCCESS &&
            reactivated_second_report
                .has_value() &&
            reactivated_second_report
                    ->accepted() &&
            reactivated_second_report
                    ->entry_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        fixed_total_molar_rate &&
            reactivated_second_report
                    ->accepted_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        fixed_total_molar_rate &&
            !reactivated_second_report
                 ->minimum_bhp_probe_total_molar_rate_mol_per_s
                 .has_value() &&
            !reactivated_second_report
                 ->minimum_bhp_to_rate_reactivation_attempted &&
            !reactivated_second_report
                 ->minimum_bhp_to_rate_reactivation_accepted &&
            reactivation_clock
                    .accepted_time_seconds() ==
                reactivated_second_time_n +
                    reactivation_dt_seconds &&
            reactivation_clock
                    .accepted_step_count() ==
                2U,
        "reactivated accepted rate control did not persist directly into the next timestep");

    // Rate-capacity deadband: p_min capacity does not clear target + Δq, so
    // no augmented rate candidate may be created.
    auto capacity_deadband_source_context =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    multi_context,
                    minimum_bhp_pa);
    auto capacity_deadband_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &capacity_deadband_source_context,
            &provider);
    require_collective(
        capacity_deadband_system
                    ->set_trial_timestep_seconds(
                        reactivation_dt_seconds) ==
                PETSC_SUCCESS,
        "failed to prepare rate-capacity deadband system");

    auto capacity_deadband_control =
        wdp::
            AcceptedFixedTotalMolarRateWellControlState3D::
                minimum_bottom_hole_pressure(
                    minimum_bhp_pa);
    fdp::AcceptedPhysicalTimeClock3D
        capacity_deadband_clock{
            0.0,
            reactivation_dt_seconds};
    auto capacity_deadband_options =
        reactivation_options;
    const double capacity_deadband_target =
        0.95 *
        reactivation_probe_total_molar_rate;
    capacity_deadband_options
        .minimum_bhp_release_rate_margin_mol_per_s =
        0.10 *
        reactivation_probe_total_molar_rate;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        capacity_deadband_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                capacity_deadband_system.get(),
                &capacity_deadband_source_context,
                capacity_deadband_target,
                capacity_deadband_options,
                &capacity_deadband_clock,
                &capacity_deadband_control,
                &capacity_deadband_report);
    require_collective(
        error == PETSC_SUCCESS &&
            capacity_deadband_report
                .has_value() &&
            capacity_deadband_report
                    ->accepted() &&
            capacity_deadband_report
                    ->accepted_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            capacity_deadband_report
                    ->minimum_bhp_probe_total_molar_rate_mol_per_s
                    .has_value() &&
            !capacity_deadband_report
                 ->minimum_bhp_to_rate_reactivation_attempted &&
            !capacity_deadband_report
                 ->minimum_bhp_to_rate_reactivation_accepted &&
            !capacity_deadband_report
                 ->discarded_minimum_bhp_candidate
                 .has_value() &&
            !capacity_deadband_report
                 ->discarded_rate_reactivation_candidate
                 .has_value() &&
            !capacity_deadband_report
                 ->accepted_solve
                 .has_value() &&
            capacity_deadband_report
                    ->accepted_fixed_bhp_solve
                    .has_value() &&
            capacity_deadband_control
                    .control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure,
        "minimum-BHP rate-capacity hysteresis gate reactivated rate control inside the deadband");

    // Pressure deadband: capacity clears Δq and the augmented rate root
    // converges, but an intentionally wide Δp guard rejects that candidate and
    // keeps the already converged fixed-BHP probe as the accepted state.
    auto pressure_deadband_source_context =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    multi_context,
                    minimum_bhp_pa);
    auto pressure_deadband_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &pressure_deadband_source_context,
            &provider);
    require_collective(
        pressure_deadband_system
                    ->set_trial_timestep_seconds(
                        reactivation_dt_seconds) ==
                PETSC_SUCCESS,
        "failed to prepare pressure-deadband reactivation system");

    auto pressure_deadband_control =
        wdp::
            AcceptedFixedTotalMolarRateWellControlState3D::
                minimum_bottom_hole_pressure(
                    minimum_bhp_pa);
    fdp::AcceptedPhysicalTimeClock3D
        pressure_deadband_clock{
            0.0,
            reactivation_dt_seconds};
    auto pressure_deadband_options =
        reactivation_options;
    pressure_deadband_options
        .minimum_bhp_release_pressure_margin_pa =
        100.0;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        pressure_deadband_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_3d(
                PETSC_COMM_WORLD,
                pressure_deadband_system.get(),
                &pressure_deadband_source_context,
                reactivation_target_rate,
                pressure_deadband_options,
                &pressure_deadband_clock,
                &pressure_deadband_control,
                &pressure_deadband_report);
    require_collective(
        error == PETSC_SUCCESS &&
            pressure_deadband_report
                .has_value() &&
            pressure_deadband_report
                    ->accepted() &&
            pressure_deadband_report
                    ->accepted_control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            pressure_deadband_report
                    ->minimum_bhp_to_rate_reactivation_attempted &&
            !pressure_deadband_report
                 ->minimum_bhp_to_rate_reactivation_accepted &&
            pressure_deadband_report
                    ->discarded_rate_reactivation_candidate
                    .has_value() &&
            !pressure_deadband_report
                 ->discarded_minimum_bhp_candidate
                 .has_value() &&
            !pressure_deadband_report
                 ->accepted_solve
                 .has_value() &&
            pressure_deadband_report
                    ->accepted_fixed_bhp_solve
                    .has_value() &&
            pressure_deadband_report
                    ->discarded_rate_reactivation_candidate
                    ->bottom_hole_pressure_pa <
                minimum_bhp_pa +
                    *pressure_deadband_options
                         .minimum_bhp_release_pressure_margin_pa &&
            pressure_deadband_control
                    .control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            pressure_deadband_control
                    .bottom_hole_pressure_pa ==
                minimum_bhp_pa,
        "minimum-BHP pressure hysteresis gate accepted a rate candidate inside the release deadband");

    require_collective(
        VecDestroy(
            &reactivation_probe_difference) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &reactivation_probe_state) ==
                PETSC_SUCCESS,
        "minimum-BHP reactivation regression cleanup failed");

    Vec discarded_rate_difference = nullptr;
    error =
        VecDuplicate(
            switch_reservoir_system
                ->initial_state(),
            &discarded_rate_difference);
    if (error == PETSC_SUCCESS) {
        error =
            VecCopy(
                switch_reservoir_system
                    ->initial_state(),
                discarded_rate_difference);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAXPY(
                discarded_rate_difference,
                PetscScalar{-1.0},
                reservoir_solution);
    }
    PetscReal discarded_rate_difference_norm =
        -1.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                discarded_rate_difference,
                NORM_2,
                &discarded_rate_difference_norm);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            discarded_rate_difference_norm >
                1.0e-9,
        "minimum-BHP switch appears to have committed the discarded rate-control reservoir candidate");

    require_collective(
        VecDestroy(
            &discarded_rate_difference) ==
            PETSC_SUCCESS,
        "minimum-BHP switching regression cleanup failed");

    // Controlled post-SNES phase-transition restart: cell30 starts 2P,
    // cell60 stays frozen 3P. The first selected rate candidate is scanned,
    // discarded, rebuilt to 3P by stable identity, and the same 1 s physical
    // timestep is solved again from the entry accepted rate-control state.
    const auto aqueous_identity =
        mixed_physical_phase_identity(
            "aqueous");
    const auto hc0_identity =
        mixed_physical_phase_identity(
            "hydrocarbon-0");
    const auto hc1_identity =
        mixed_physical_phase_identity(
            "hydrocarbon-1");
    const wdp::
        FixedBhpPhaseIdentityInjectionEnthalpy3D
        controlled_phase_enthalpies{
            "fixture/controlled-transition-enthalpy/v1",
            {
                {aqueous_identity, 1000.0},
                {hc0_identity, 2000.0},
                {hc1_identity, 3000.0}}};

    const auto make_transition_well =
        [&](double bhp_pa) {
            auto reboundable30 =
                wdp::
                    FixedBhpPeacemanWellSourceEvaluatorContext3D::
                        create_phase_identity_bound(
                            mesh::GlobalEntityId{
                                UINT64_C(30)},
                            well::
                                make_peaceman_well_index_3d(
                                    {10.0, 10.0, 5.0},
                                    {
                                        1.0e-8,
                                        1.0e-8,
                                        1.0e-8},
                                    well::
                                        AxisAlignedWellDirection3D::z,
                                    0.10,
                                    0.0),
                            bhp_pa,
                            controlled_phase_enthalpies,
                            phase_identity_map_for_stable(
                                UINT64_C(30)),
                            "fixture/controlled-transition-cell30/v1");
            auto frozen60 =
                wdp::
                    FixedBhpPeacemanWellSourceEvaluatorContext3D::
                        create(
                            mesh::GlobalEntityId{
                                UINT64_C(60)},
                            well::
                                make_peaceman_well_index_3d(
                                    {12.0, 9.0, 4.0},
                                    {
                                        2.0e-8,
                                        1.5e-8,
                                        2.5e-8},
                                    well::
                                        AxisAlignedWellDirection3D::z,
                                    0.12,
                                    0.1),
                            bhp_pa,
                            wd::
                                FixedBhpInjectionEnthalpy3D{
                                    "fixture/controlled-transition-cell60-enthalpy/v1",
                                    {
                                        1100.0,
                                        2200.0,
                                        3300.0}},
                            "fixture/controlled-transition-cell60/v1");
            return wdp::
                FixedBhpMultiConnectionWellSourceEvaluatorContext3D::
                    create(
                        "fixture/controlled-transition-well/v1",
                        {
                            frozen60,
                            reboundable30});
        };

    const flow::FrozenActivePhaseIdentityMap
        controlled_transition_three_phase_map{
            {
                aqueous_identity,
                hc0_identity,
                hc1_identity}};

    // Construct a transition-specific control target from the 2P chart, then
    // prove directly with independent augmented solves that the same target
    // has a valid root on both frozen phase charts. Transactional restart does
    // not require those roots to occur in the same BHP interval.
    auto transition_fixed_bhp_rate =
        [&](double bhp_pa,
            bool cell30_is_three_phase)
            -> std::optional<double> {
            auto reference_source =
                wdp::
                    FixedTotalMolarRateWellSourceEvaluatorContext3D::
                        create(
                            make_transition_well(
                                bhp_pa),
                            bhp_pa);
            if (cell30_is_three_phase) {
                reference_source
                    .rebind_connection(
                        mesh::GlobalEntityId{
                            UINT64_C(30)},
                        controlled_transition_three_phase_map);
            }
            reference_source
                .begin_fixed_bhp_reservoir_evaluation(
                    bhp_pa);
            auto reference_system =
                make_fixed_total_molar_rate_reservoir_system(
                    rank,
                    schedule,
                    partition,
                    bridge,
                    pattern,
                    audit,
                    &reference_source,
                    &provider,
                    cell30_is_three_phase);

            Vec reference_state = nullptr;
            std::optional<
                fdp::
                    VariableCardinalityNaturalVariableSnesSolveReport3D>
                reference_solve;
            const PetscErrorCode reference_error =
                reference_system
                    ->solve(
                        &reference_state,
                        &reference_solve);
            const bool solved =
                reference_error ==
                    PETSC_SUCCESS &&
                reference_state !=
                    nullptr &&
                reference_solve
                    .has_value() &&
                static_cast<int>(
                    reference_solve
                        ->converged_reason) >
                    0;
            int local_solved =
                solved ? 1 : 0;
            int global_solved = 0;
            if (MPI_Allreduce(
                    &local_solved,
                    &global_solved,
                    1,
                    MPI_INT,
                    MPI_MIN,
                    PETSC_COMM_WORLD) !=
                MPI_SUCCESS) {
                if (reference_state !=
                    nullptr) {
                    (void)VecDestroy(
                        &reference_state);
                }
                throw std::runtime_error(
                    "failed to reduce transition fixed-BHP solve status");
            }
            if (global_solved == 0) {
                if (reference_state !=
                    nullptr) {
                    (void)VecDestroy(
                        &reference_state);
                }
                return std::nullopt;
            }

            std::uint64_t
                authoritative_count = 0U;
            const auto rate =
                global_fixed_bhp_well_production_rate(
                    *reference_system,
                    reference_state,
                    reference_source.well(),
                    bhp_pa,
                    &authoritative_count);
            const double total_molar_rate =
                rate[0] +
                rate[1] +
                rate[2];
            const PetscErrorCode destroy_error =
                VecDestroy(
                    &reference_state);
            if (destroy_error !=
                PETSC_SUCCESS ||
                authoritative_count !=
                    2U ||
                !std::isfinite(
                    total_molar_rate) ||
                !(total_molar_rate > 0.0)) {
                return std::nullopt;
            }
            return total_molar_rate;
        };

    auto transition_rate_root =
        [&](double target,
            bool cell30_is_three_phase)
            -> std::optional<double> {
            auto oracle_source =
                wdp::
                    FixedTotalMolarRateWellSourceEvaluatorContext3D::
                        create(
                            make_transition_well(
                                initial_bhp_pa),
                            initial_bhp_pa);
            if (cell30_is_three_phase) {
                oracle_source
                    .rebind_connection(
                        mesh::GlobalEntityId{
                            UINT64_C(30)},
                        controlled_transition_three_phase_map);
            }
            auto oracle_system =
                make_fixed_total_molar_rate_reservoir_system(
                    rank,
                    schedule,
                    partition,
                    bridge,
                    pattern,
                    audit,
                    &oracle_source,
                    &provider,
                    cell30_is_three_phase);
            std::unique_ptr<
                wdp::
                    FixedTotalMolarRateWellControlSystem3D>
                oracle_control;
            PetscErrorCode oracle_error =
                wdp::
                    FixedTotalMolarRateWellControlSystem3D::
                        create(
                            PETSC_COMM_WORLD,
                            oracle_system.get(),
                            &oracle_source,
                            initial_bhp_pa,
                            target,
                            &oracle_control);

            Vec oracle_state = nullptr;
            std::optional<
                wdp::
                    FixedTotalMolarRateWellControlSolveReport3D>
                oracle_report;
            if (oracle_error ==
                PETSC_SUCCESS) {
                oracle_error =
                    oracle_control
                        ->solve(
                            &oracle_state,
                            &oracle_report);
            }

            bool root_valid =
                oracle_error ==
                    PETSC_SUCCESS &&
                oracle_state !=
                    nullptr &&
                oracle_report
                    .has_value() &&
                oracle_report
                    ->converged() &&
                std::isfinite(
                    oracle_report
                        ->bottom_hole_pressure_pa) &&
                oracle_report
                        ->bottom_hole_pressure_pa >
                    0.0 &&
                std::abs(
                    oracle_report
                        ->total_molar_rate_residual_mol_per_s()) <=
                    1.0e-8 *
                        std::max(
                            1.0,
                            std::abs(
                                target));
            int local_root_valid =
                root_valid ? 1 : 0;
            int global_root_valid = 0;
            if (MPI_Allreduce(
                    &local_root_valid,
                    &global_root_valid,
                    1,
                    MPI_INT,
                    MPI_MIN,
                    PETSC_COMM_WORLD) !=
                MPI_SUCCESS) {
                if (oracle_state !=
                    nullptr) {
                    (void)VecDestroy(
                        &oracle_state);
                }
                throw std::runtime_error(
                    "failed to reduce transition rate-root status");
            }

            double root_bhp = 0.0;
            if (global_root_valid != 0) {
                root_bhp =
                    oracle_report
                        ->bottom_hole_pressure_pa;
            }
            if (oracle_state !=
                nullptr) {
                if (VecDestroy(
                        &oracle_state) !=
                    PETSC_SUCCESS) {
                    throw std::runtime_error(
                        "transition rate-root cleanup failed");
                }
            }
            if (global_root_valid == 0) {
                return std::nullopt;
            }
            return root_bhp;
        };

    const auto
        transition_two_phase_rate_at_entry =
            transition_fixed_bhp_rate(
                initial_bhp_pa,
                false);
    const auto
        transition_two_phase_rate_at_reference =
            transition_fixed_bhp_rate(
                target_reference_bhp_pa,
                false);
    require_collective(
        transition_two_phase_rate_at_entry
                .has_value() &&
            transition_two_phase_rate_at_reference
                .has_value() &&
            target_reference_bhp_pa <
                initial_bhp_pa &&
            *transition_two_phase_rate_at_reference >
                *transition_two_phase_rate_at_entry,
        "2P transition chart did not provide a lower-BHP production interval");

    constexpr std::array<double, 6>
        transition_target_fractions{
            0.01,
            0.025,
            0.05,
            0.10,
            0.25,
            0.50};

    std::optional<double>
        transition_target_rate_candidate;
    std::optional<double>
        transition_two_phase_rate_root_bhp;
    std::optional<double>
        transition_three_phase_rate_root_bhp;

    const double transition_rate_span =
        *transition_two_phase_rate_at_reference -
        *transition_two_phase_rate_at_entry;
    for (double fraction :
         transition_target_fractions) {
        const double candidate =
            *transition_two_phase_rate_at_entry +
            fraction *
                transition_rate_span;
        const auto two_phase_root =
            transition_rate_root(
                candidate,
                false);
        if (!two_phase_root.has_value() ||
            !(*two_phase_root <
              initial_bhp_pa)) {
            continue;
        }
        const auto three_phase_root =
            transition_rate_root(
                candidate,
                true);
        if (!three_phase_root.has_value()) {
            continue;
        }

        transition_target_rate_candidate =
            candidate;
        transition_two_phase_rate_root_bhp =
            two_phase_root;
        transition_three_phase_rate_root_bhp =
            three_phase_root;
        break;
    }

    require_collective(
        transition_target_rate_candidate
                .has_value() &&
            transition_two_phase_rate_root_bhp
                .has_value() &&
            transition_three_phase_rate_root_bhp
                .has_value(),
        "no transition-specific target has valid frozen 2P and 3P augmented rate roots");

    const double transition_target_rate =
        *transition_target_rate_candidate;
    const double
        transition_two_phase_rate_root_bhp_pa =
            *transition_two_phase_rate_root_bhp;

    auto transition_rate_source =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    make_transition_well(
                        initial_bhp_pa),
                    initial_bhp_pa);
    auto transition_rate_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &transition_rate_source,
            &provider);

    auto transition_rate_control =
        wdp::
            AcceptedFixedTotalMolarRateWellControlState3D::
                fixed_total_molar_rate(
                    initial_bhp_pa);
    fdp::AcceptedPhysicalTimeClock3D
        transition_rate_clock{
            0.0,
            1.0};
    ControlledWellPhaseTransitionFixture
        transition_rate_fixture{
            rank,
            &schedule,
            &partition,
            &bridge,
            &pattern,
            audit,
            &provider,
            &transition_rate_source,
            &transition_rate_control,
            &transition_rate_clock,
            wdp::
                FixedTotalMolarRatePhysicalTimestepControlMode3D::
                    fixed_total_molar_rate,
            initial_bhp_pa,
            0U,
            0U,
            true,
            false,
            {},
            true,
            {}};

    wdp::
        FixedTotalMolarRatePhysicalTimestepDriverOptions3D
        transition_rate_options;
    transition_rate_options.adaptive
        .minimum_timestep_seconds =
        1.0;
    transition_rate_options.adaptive
        .maximum_timestep_seconds =
        1.0;
    transition_rate_options.adaptive
        .maximum_retries =
        1U;
    transition_rate_options.phase_transition
        .max_transition_restarts =
        2U;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        transition_rate_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                &transition_rate_system,
                &transition_rate_source,
                transition_target_rate,
                {
                    &controlled_well_transition_scan,
                    &transition_rate_fixture,
                    &controlled_well_transition_rebuild,
                    &transition_rate_fixture},
                transition_rate_options,
                &transition_rate_clock,
                &transition_rate_control,
                &transition_rate_report);

    const auto* rebound_rate_cell30 =
        transition_rate_source
            .well()
            .find_connection(
                mesh::GlobalEntityId{
                    UINT64_C(30)});
    const bool transition_rate_report_present =
        transition_rate_report.has_value();
    const bool transition_rate_accepted =
        transition_rate_report_present &&
        transition_rate_report->accepted();
    const bool transition_rate_record_present =
        transition_rate_report_present &&
        transition_rate_report
            ->accepted_record
            .has_value();
    const bool transition_rate_solve_present =
        transition_rate_report_present &&
        transition_rate_report
            ->accepted_solve
            .has_value();
    const bool transition_rate_reservoir_report_present =
        transition_rate_solve_present &&
        transition_rate_report
            ->accepted_solve
            ->reservoir_nonlinear_solve
            .has_value();
    const std::string transition_rate_restart_diagnostic =
        std::string{"error="} +
        std::to_string(
            static_cast<int>(error)) +
        ", report=" +
        std::to_string(
            transition_rate_report_present ? 1 : 0) +
        ", accepted=" +
        std::to_string(
            transition_rate_accepted ? 1 : 0) +
        ", adaptive_outcome=" +
        std::to_string(
            transition_rate_report_present
                ? static_cast<int>(
                      transition_rate_report
                          ->adaptive
                          .outcome)
                : -999) +
        ", adaptive_retries=" +
        std::to_string(
            transition_rate_report_present
                ? transition_rate_report
                      ->adaptive
                      .retries
                : 999U) +
        ", adaptive_attempts=" +
        std::to_string(
            transition_rate_report_present
                ? transition_rate_report
                      ->adaptive
                      .attempts
                      .size()
                : 999U) +
        ", attempt_outcome=" +
        std::to_string(
            transition_rate_report_present &&
                    !transition_rate_report
                         ->adaptive
                         .attempts
                         .empty()
                ? static_cast<int>(
                      transition_rate_report
                          ->adaptive
                          .attempts
                          .back()
                          .result
                          .outcome)
                : -999) +
        ", attempt_decision=" +
        std::to_string(
            transition_rate_report_present &&
                    !transition_rate_report
                         ->adaptive
                         .attempts
                         .empty()
                ? static_cast<int>(
                      transition_rate_report
                          ->adaptive
                          .attempts
                          .back()
                          .decision)
                : -999) +
        ", attempt_nonlinear_iterations=" +
        std::to_string(
            transition_rate_report_present &&
                    !transition_rate_report
                         ->adaptive
                         .attempts
                         .empty()
                ? transition_rate_report
                      ->adaptive
                      .attempts
                      .back()
                      .result
                      .nonlinear_iterations
                : -999) +
        ", attempt_function_domain_errors=" +
        std::to_string(
            transition_rate_report_present &&
                    !transition_rate_report
                         ->adaptive
                         .attempts
                         .empty()
                ? transition_rate_report
                      ->adaptive
                      .attempts
                      .back()
                      .result
                      .function_domain_errors
                : -999) +
        ", attempt_jacobian_domain_errors=" +
        std::to_string(
            transition_rate_report_present &&
                    !transition_rate_report
                         ->adaptive
                         .attempts
                         .empty()
                ? transition_rate_report
                      ->adaptive
                      .attempts
                      .back()
                      .result
                      .jacobian_domain_errors
                : -999) +
        ", attempt_line_search_changes=" +
        std::to_string(
            transition_rate_report_present &&
                    !transition_rate_report
                         ->adaptive
                         .attempts
                         .empty()
                ? transition_rate_report
                      ->adaptive
                      .attempts
                      .back()
                      .result
                      .line_search_direction_changes
                : -999) +
        ", attempt_transition_restarts=" +
        std::to_string(
            transition_rate_report_present &&
                    !transition_rate_report
                         ->adaptive
                         .attempts
                         .empty()
                ? transition_rate_report
                      ->adaptive
                      .attempts
                      .back()
                      .result
                      .phase_transition_restarts
                : 999U) +
        ", report_restarts=" +
        std::to_string(
            transition_rate_report_present
                ? transition_rate_report
                      ->phase_transition_restarts
                : 999U) +
        ", record=" +
        std::to_string(
            transition_rate_record_present ? 1 : 0) +
        ", record_restarts=" +
        std::to_string(
            transition_rate_record_present
                ? transition_rate_report
                      ->accepted_record
                      ->phase_transition_restarts
                : 999U) +
        ", scans=" +
        std::to_string(
            transition_rate_fixture.scans) +
        ", rebuilds=" +
        std::to_string(
            transition_rate_fixture.rebuild_calls) +
        ", anchor_rebuild=" +
        std::to_string(
            transition_rate_fixture
                    .accepted_anchor_unchanged_at_rebuild
                ? 1
                : 0) +
        ", anchor_scan=" +
        std::to_string(
            transition_rate_fixture
                    .accepted_anchor_unchanged_at_scan
                ? 1
                : 0) +
        ", previous_totals=" +
        std::to_string(
            transition_rate_fixture
                    .rebuilt_previous_totals_captured
                ? 1
                : 0) +
        ", cell30_phases=" +
        std::to_string(
            transition_rate_system
                ->numbering()
                .cell(
                    mesh::LocalIndex{2U})
                .phase_count) +
        ", cell60_phases=" +
        std::to_string(
            transition_rate_system
                ->numbering()
                .cell(
                    mesh::LocalIndex{5U})
                .phase_count) +
        ", control=" +
        std::to_string(
            static_cast<int>(
                transition_rate_control
                    .control)) +
        ", solve=" +
        std::to_string(
            transition_rate_solve_present ? 1 : 0) +
        ", reservoir_report=" +
        std::to_string(
            transition_rate_reservoir_report_present
                ? 1
                : 0) +
        ", clock=" +
        std::to_string(
            transition_rate_clock
                .accepted_time_seconds()) +
        ", steps=" +
        std::to_string(
            transition_rate_clock
                .accepted_step_count());

    require_collective(
        error == PETSC_SUCCESS &&
            transition_rate_report_present &&
            transition_rate_accepted &&
            transition_rate_report
                    ->phase_transition_restarts ==
                1U &&
            transition_rate_record_present &&
            transition_rate_report
                    ->accepted_record
                    ->phase_transition_restarts ==
                1U &&
            transition_rate_fixture.scans ==
                2U &&
            transition_rate_fixture.rebuild_calls ==
                1U &&
            transition_rate_fixture
                .accepted_anchor_unchanged_at_rebuild &&
            transition_rate_fixture
                .accepted_anchor_unchanged_at_scan &&
            transition_rate_fixture
                .rebuilt_previous_totals_captured &&
            transition_rate_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                3U &&
            transition_rate_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{5U})
                    .phase_count ==
                3U &&
            transition_rate_control
                    .control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        fixed_total_molar_rate &&
            transition_rate_solve_present &&
            transition_rate_reservoir_report_present &&
            transition_rate_report
                    ->accepted_solve
                    ->global_scalar_count ==
                transition_rate_system
                    ->numbering()
                    .petsc_global_scalar_count() +
                    1 &&
            transition_rate_clock
                    .accepted_time_seconds() ==
                1.0 &&
            transition_rate_clock
                    .accepted_step_count() ==
                1U &&
            rebound_rate_cell30 !=
                nullptr &&
            rebound_rate_cell30
                    ->active_phase_identities()
                    .has_value() &&
            rebound_rate_cell30
                    ->active_phase_identities()
                    ->phase_count() ==
                3U,
        std::string{
            "rate-controlled post-SNES transition invariant failed: "} +
            transition_rate_restart_diagnostic);

    std::uint64_t
        transition_rate_authoritative_count = 0U;
    const auto transition_rate_well_rate =
        global_fixed_bhp_well_production_rate(
            *transition_rate_system,
            transition_rate_system
                ->initial_state(),
            transition_rate_source.well(),
            transition_rate_control
                .bottom_hole_pressure_pa,
            &transition_rate_authoritative_count);
    const double
        transition_rate_total_molar =
            transition_rate_well_rate[0] +
            transition_rate_well_rate[1] +
            transition_rate_well_rate[2];
    require_collective(
        transition_rate_authoritative_count ==
            2U,
        "rate-controlled transition lost owner-only two-connection aggregation");
    const double
        transition_rate_residual_tolerance =
            1.0e-8 *
            std::max(
                1.0,
                std::abs(
                    transition_target_rate));
    const double
        transition_report_rate_residual =
            transition_rate_report
                ->accepted_solve
                ->total_molar_rate_residual_mol_per_s();
    require_collective(
        std::abs(
            transition_report_rate_residual) <=
            transition_rate_residual_tolerance,
        std::string{
            "transitioned augmented rate solve violated control residual: residual="} +
            std::to_string(
                transition_report_rate_residual) +
            ", target=" +
            std::to_string(
                transition_target_rate));
    require_collective(
        std::abs(
            transition_rate_total_molar -
            transition_target_rate) <=
            transition_rate_residual_tolerance,
        std::string{
            "transitioned committed whole-well rate drifted from control root: actual="} +
            std::to_string(
                transition_rate_total_molar) +
            ", target=" +
            std::to_string(
                transition_target_rate) +
            ", accepted_bhp=" +
            std::to_string(
                transition_rate_control
                    .bottom_hole_pressure_pa));

    const auto transition_rate_local_final =
        owned_conserved_totals(
            *transition_rate_system,
            transition_rate_system
                ->initial_state());
    std::array<double, 4>
        transition_rate_global_final{};
    require_collective(
        MPI_Allreduce(
            transition_rate_local_final.data(),
            transition_rate_global_final.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce transitioned rate-control final totals");
    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        near_collective(
            transition_rate_global_final[
                quantity],
            transition_rate_fixture
                    .rebuilt_previous_totals[
                        quantity] -
                transition_rate_well_rate[
                    quantity],
            4.0e-7,
            quantity < 3U
                ? 4.0e-8
                : 4.0e-7);
    }

    bool transition_rate_history_matches =
        false;
    require_collective(
        transition_rate_system
                ->accepted_history_matches_state(
                    transition_rate_system
                        ->initial_state(),
                    &transition_rate_history_matches) ==
                PETSC_SUCCESS &&
            transition_rate_history_matches,
        "transitioned rate-control accepted history did not commit exactly once");

    // Cross the two state machines in one physical timestep. Entry control is
    // accepted RATE, but the first 2P rate candidate is known to violate the
    // already validated p_min and therefore selects a tentative fixed-BHP
    // candidate. That selected BHP candidate triggers 2P -> 3P. The phase
    // restart must discard the complete trial control generation and restart
    // the rebuilt system from the accepted RATE/BHP anchor, not from p_min.
    auto transition_switch_source =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    make_transition_well(
                        initial_bhp_pa),
                    initial_bhp_pa);
    auto transition_switch_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &transition_switch_source,
            &provider);
    auto transition_switch_control =
        wdp::
            AcceptedFixedTotalMolarRateWellControlState3D::
                fixed_total_molar_rate(
                    initial_bhp_pa);
    fdp::AcceptedPhysicalTimeClock3D
        transition_switch_clock{
            0.0,
            1.0};
    ControlledWellPhaseTransitionFixture
        transition_switch_fixture{
            rank,
            &schedule,
            &partition,
            &bridge,
            &pattern,
            audit,
            &provider,
            &transition_switch_source,
            &transition_switch_control,
            &transition_switch_clock,
            wdp::
                FixedTotalMolarRatePhysicalTimestepControlMode3D::
                    fixed_total_molar_rate,
            initial_bhp_pa,
            0U,
            0U,
            true,
            false,
            {},
            true,
            {}};

    const double
        transition_switch_minimum_bhp_pa =
            0.5 *
            (transition_two_phase_rate_root_bhp_pa +
             initial_bhp_pa);
    require_collective(
        transition_switch_minimum_bhp_pa >
                transition_two_phase_rate_root_bhp_pa &&
            transition_switch_minimum_bhp_pa <
                initial_bhp_pa,
        "transactional switch p_min is not strictly above the feasible 2P rate root");

    auto transition_switch_options =
        transition_rate_options;
    transition_switch_options
        .minimum_bottom_hole_pressure_pa =
        transition_switch_minimum_bhp_pa;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        transition_switch_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                &transition_switch_system,
                &transition_switch_source,
                transition_target_rate,
                {
                    &controlled_well_transition_scan,
                    &transition_switch_fixture,
                    &controlled_well_transition_rebuild,
                    &transition_switch_fixture},
                transition_switch_options,
                &transition_switch_clock,
                &transition_switch_control,
                &transition_switch_report);

    const wdp::
        FixedTotalMolarRateWellControlSolveReport3D*
        restarted_rate_candidate =
            nullptr;
    if (transition_switch_report
            .has_value()) {
        if (transition_switch_report
                ->accepted_solve
                .has_value()) {
            restarted_rate_candidate =
                &*transition_switch_report
                      ->accepted_solve;
        } else if (
            transition_switch_report
                ->discarded_rate_control_candidate
                .has_value()) {
            restarted_rate_candidate =
                &*transition_switch_report
                      ->discarded_rate_control_candidate;
        }
    }

    const auto* rebound_switch_cell30 =
        transition_switch_source
            .well()
            .find_connection(
                mesh::GlobalEntityId{
                    UINT64_C(30)});
    require_collective(
        error == PETSC_SUCCESS &&
            transition_switch_report
                .has_value() &&
            transition_switch_report
                    ->accepted() &&
            transition_switch_report
                    ->phase_transition_restarts ==
                1U &&
            transition_switch_report
                    ->accepted_record
                    ->phase_transition_restarts ==
                1U &&
            transition_switch_fixture.scans ==
                2U &&
            transition_switch_fixture
                    .scanned_source_bottom_hole_pressure_pa
                    .size() ==
                2U &&
            transition_switch_fixture
                    .scanned_source_bottom_hole_pressure_pa[
                        0U] ==
                transition_switch_minimum_bhp_pa &&
            transition_switch_fixture
                .accepted_anchor_unchanged_at_rebuild &&
            transition_switch_fixture
                .accepted_anchor_unchanged_at_scan &&
            restarted_rate_candidate !=
                nullptr &&
            restarted_rate_candidate
                    ->initial_bottom_hole_pressure_pa ==
                initial_bhp_pa &&
            transition_switch_clock
                    .accepted_time_seconds() ==
                1.0 &&
            transition_switch_clock
                    .accepted_step_count() ==
                1U &&
            rebound_switch_cell30 !=
                nullptr &&
            rebound_switch_cell30
                    ->active_phase_identities()
                    .has_value() &&
            rebound_switch_cell30
                    ->active_phase_identities()
                    ->phase_count() ==
                3U,
        "rate-to-BHP trial control generation leaked across the phase-transition restart");

    std::uint64_t
        transition_switch_authoritative_count =
            0U;
    const auto transition_switch_well_rate =
        global_fixed_bhp_well_production_rate(
            *transition_switch_system,
            transition_switch_system
                ->initial_state(),
            transition_switch_source.well(),
            transition_switch_control
                .bottom_hole_pressure_pa,
            &transition_switch_authoritative_count);
    require_collective(
        transition_switch_authoritative_count ==
            2U,
        "control/phase transactional restart lost owner-only two-connection aggregation");

    const double
        transition_switch_total_molar_rate =
            transition_switch_well_rate[0] +
            transition_switch_well_rate[1] +
            transition_switch_well_rate[2];
    if (transition_switch_control.control ==
        wdp::
            FixedTotalMolarRatePhysicalTimestepControlMode3D::
                fixed_total_molar_rate) {
        require_collective(
            transition_switch_report
                    ->accepted_solve
                    .has_value() &&
                std::abs(
                    transition_switch_total_molar_rate -
                    transition_target_rate) <=
                    transition_rate_residual_tolerance &&
                std::abs(
                    transition_switch_report
                        ->accepted_solve
                        ->total_molar_rate_residual_mol_per_s()) <=
                    transition_rate_residual_tolerance,
            "restarted RATE generation did not satisfy the rate residual contract");
    } else {
        require_collective(
            transition_switch_control.control ==
                    wdp::
                        FixedTotalMolarRatePhysicalTimestepControlMode3D::
                            minimum_bottom_hole_pressure &&
                transition_switch_control
                        .bottom_hole_pressure_pa ==
                    transition_switch_minimum_bhp_pa &&
                transition_switch_report
                        ->accepted_fixed_bhp_solve
                        .has_value(),
            "restarted arbitration produced an invalid final minimum-BHP control state");
    }

    const auto transition_switch_local_final =
        owned_conserved_totals(
            *transition_switch_system,
            transition_switch_system
                ->initial_state());
    std::array<double, 4>
        transition_switch_global_final{};
    require_collective(
        MPI_Allreduce(
            transition_switch_local_final.data(),
            transition_switch_global_final.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce control/phase transactional final totals");
    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        near_collective(
            transition_switch_global_final[
                quantity],
            transition_switch_fixture
                    .rebuilt_previous_totals[
                        quantity] -
                transition_switch_well_rate[
                    quantity],
            4.0e-7,
            quantity < 3U
                ? 4.0e-8
                : 4.0e-7);
    }

    bool transition_switch_history_matches =
        false;
    require_collective(
        transition_switch_system
                ->accepted_history_matches_state(
                    transition_switch_system
                        ->initial_state(),
                    &transition_switch_history_matches) ==
                PETSC_SUCCESS &&
            transition_switch_history_matches,
        "control/phase transactional restart did not commit exactly one final history");

    // Repeat from an already accepted minimum-BHP mode. The 2P fixed-BHP
    // candidate is discarded by the phase transition, the completion is
    // rebound to 3P, and the same timestep restarts from the accepted
    // minimum-BHP mode rather than carrying any trial control decision.
    constexpr double
        transition_minimum_bhp_pa =
            5.0;
    auto transition_bhp_source =
        wdp::
            FixedTotalMolarRateWellSourceEvaluatorContext3D::
                create(
                    make_transition_well(
                        transition_minimum_bhp_pa),
                    transition_minimum_bhp_pa);
    auto transition_bhp_system =
        make_fixed_total_molar_rate_reservoir_system(
            rank,
            schedule,
            partition,
            bridge,
            pattern,
            audit,
            &transition_bhp_source,
            &provider);
    auto transition_bhp_control =
        wdp::
            AcceptedFixedTotalMolarRateWellControlState3D::
                minimum_bottom_hole_pressure(
                    transition_minimum_bhp_pa);
    fdp::AcceptedPhysicalTimeClock3D
        transition_bhp_clock{
            0.0,
            1.0};
    ControlledWellPhaseTransitionFixture
        transition_bhp_fixture{
            rank,
            &schedule,
            &partition,
            &bridge,
            &pattern,
            audit,
            &provider,
            &transition_bhp_source,
            &transition_bhp_control,
            &transition_bhp_clock,
            wdp::
                FixedTotalMolarRatePhysicalTimestepControlMode3D::
                    minimum_bottom_hole_pressure,
            transition_minimum_bhp_pa,
            0U,
            0U,
            true,
            false,
            {},
            true,
            {}};

    auto transition_bhp_options =
        transition_rate_options;
    transition_bhp_options
        .minimum_bottom_hole_pressure_pa =
        transition_minimum_bhp_pa;

    std::optional<
        wdp::
            FixedTotalMolarRatePhysicalTimestepDriverReport3D>
        transition_bhp_report;
    error =
        wdp::
            advance_fixed_total_molar_rate_controlled_physical_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                &transition_bhp_system,
                &transition_bhp_source,
                transition_target_rate,
                {
                    &controlled_well_transition_scan,
                    &transition_bhp_fixture,
                    &controlled_well_transition_rebuild,
                    &transition_bhp_fixture},
                transition_bhp_options,
                &transition_bhp_clock,
                &transition_bhp_control,
                &transition_bhp_report);

    const auto* rebound_bhp_cell30 =
        transition_bhp_source
            .well()
            .find_connection(
                mesh::GlobalEntityId{
                    UINT64_C(30)});
    require_collective(
        error == PETSC_SUCCESS &&
            transition_bhp_report
                .has_value() &&
            transition_bhp_report
                    ->accepted() &&
            transition_bhp_report
                    ->phase_transition_restarts ==
                1U &&
            transition_bhp_report
                    ->accepted_record
                    ->phase_transition_restarts ==
                1U &&
            transition_bhp_fixture.scans ==
                2U &&
            transition_bhp_fixture.rebuild_calls ==
                1U &&
            transition_bhp_fixture
                .accepted_anchor_unchanged_at_rebuild &&
            transition_bhp_fixture
                .accepted_anchor_unchanged_at_scan &&
            transition_bhp_control
                    .control ==
                wdp::
                    FixedTotalMolarRatePhysicalTimestepControlMode3D::
                        minimum_bottom_hole_pressure &&
            transition_bhp_control
                    .bottom_hole_pressure_pa ==
                transition_minimum_bhp_pa &&
            !transition_bhp_report
                 ->accepted_solve
                 .has_value() &&
            transition_bhp_report
                    ->accepted_fixed_bhp_solve
                    .has_value() &&
            transition_bhp_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                3U &&
            transition_bhp_clock
                    .accepted_time_seconds() ==
                1.0 &&
            transition_bhp_clock
                    .accepted_step_count() ==
                1U &&
            rebound_bhp_cell30 !=
                nullptr &&
            rebound_bhp_cell30
                    ->active_phase_identities()
                    .has_value() &&
            rebound_bhp_cell30
                    ->active_phase_identities()
                    ->phase_count() ==
                3U,
        "minimum-BHP post-SNES transition did not restart from accepted BHP control after rebind");

    std::uint64_t
        transition_bhp_authoritative_count = 0U;
    const auto transition_bhp_well_rate =
        global_fixed_bhp_well_production_rate(
            *transition_bhp_system,
            transition_bhp_system
                ->initial_state(),
            transition_bhp_source.well(),
            transition_minimum_bhp_pa,
            &transition_bhp_authoritative_count);
    require_collective(
        transition_bhp_authoritative_count ==
            2U,
        "minimum-BHP transition lost owner-only two-connection aggregation");

    const auto transition_bhp_local_final =
        owned_conserved_totals(
            *transition_bhp_system,
            transition_bhp_system
                ->initial_state());
    std::array<double, 4>
        transition_bhp_global_final{};
    require_collective(
        MPI_Allreduce(
            transition_bhp_local_final.data(),
            transition_bhp_global_final.data(),
            4,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce transitioned minimum-BHP final totals");
    for (std::size_t quantity = 0U;
         quantity < 4U;
         ++quantity) {
        near_collective(
            transition_bhp_global_final[
                quantity],
            transition_bhp_fixture
                    .rebuilt_previous_totals[
                        quantity] -
                transition_bhp_well_rate[
                    quantity],
            4.0e-7,
            quantity < 3U
                ? 4.0e-8
                : 4.0e-7);
    }

    bool transition_bhp_history_matches =
        false;
    require_collective(
        transition_bhp_system
                ->accepted_history_matches_state(
                    transition_bhp_system
                        ->initial_state(),
                    &transition_bhp_history_matches) ==
                PETSC_SUCCESS &&
            transition_bhp_history_matches,
        "transitioned minimum-BHP accepted history did not commit exactly once");

    require_collective(
        MatDestroy(
            &controlled_jacobian) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &reservoir_solution) ==
                PETSC_SUCCESS &&
            VecDestroy(
                &solution) ==
                PETSC_SUCCESS,
        "rate-control regression cleanup failed");

    audit->well_timestep_compressibility =
        false;
}

void run_phase_transition_rebound_fixed_bhp_case(
    int rank,
    const dp::
        ParallelOwnedConnectionSchedule3D&
            schedule,
    const mesh::PartitionSnapshot&
        partition,
    const dp::
        PetscMpiAijSymbolicPreallocation3D&
            bridge,
    const dp::
        OwnedCellStructuralColumnPatternSnapshot3D&
            pattern,
    DispatchAudit* audit,
    flow::
        Pr76AbsentPhasePotentialExtensionProvider<
            double>* provider) {
    if (audit == nullptr ||
        provider == nullptr) {
        throw std::invalid_argument(
            "invalid phase-transition fixed-BHP well fixture");
    }

    audit->well_timestep_compressibility =
        true;

    const auto aqueous =
        mixed_physical_phase_identity(
            "aqueous");
    const auto hydrocarbon0 =
        mixed_physical_phase_identity(
            "hydrocarbon-0");
    const auto hydrocarbon1 =
        mixed_physical_phase_identity(
            "hydrocarbon-1");

    wdp::FixedBhpPhaseIdentityInjectionEnthalpy3D
        phase_enthalpy_registry{
            "fixture/fixed-bhp-phase-identity-enthalpy/v1",
            {
                {aqueous, 1100.0},
                {hydrocarbon0, 2200.0},
                {hydrocarbon1, 3300.0}}};

    // Prove the registry is identity-driven rather than slot-driven before it
    // enters the nonlinear system: a deliberately reordered active map must
    // reorder the resolved enthalpy payload accordingly.
    const auto reordered_enthalpy =
        wdp::
            resolve_fixed_bhp_injection_enthalpy_by_phase_identity_3d(
                phase_enthalpy_registry,
                flow::FrozenActivePhaseIdentityMap{
                    {hydrocarbon0, aqueous}});
    require_collective(
        reordered_enthalpy
                .specific_enthalpy_j_per_kg
                .size() ==
            2U &&
            reordered_enthalpy
                    .specific_enthalpy_j_per_kg[
                        0U] ==
                2200.0 &&
            reordered_enthalpy
                    .specific_enthalpy_j_per_kg[
                        1U] ==
                1100.0,
        "fixed-BHP injection enthalpy followed slot order instead of stable phase identity");

    auto well_context =
        wdp::
            FixedBhpPeacemanWellSourceEvaluatorContext3D::
                create_phase_identity_bound(
                    mesh::GlobalEntityId{
                        UINT64_C(30)},
                    well::make_peaceman_well_index_3d(
                        {10.0, 10.0, 5.0},
                        {
                            1.0e-8,
                            1.0e-8,
                            1.0e-8},
                        well::
                            AxisAlignedWellDirection3D::z,
                        0.10,
                        0.0),
                    15.0,
                    phase_enthalpy_registry,
                    phase_identity_map_for_stable(
                        UINT64_C(30)),
                    "fixture/fixed-bhp-transition-cell30/v1");

    require_collective(
        well_context.phase_identity_rebindable() &&
            well_context
                    .active_phase_identities()
                    .has_value() &&
            well_context
                    .active_phase_identities()
                    ->phase_count() ==
                2U &&
            well_context
                    .injection_enthalpy()
                    .specific_enthalpy_j_per_kg ==
                std::vector<double>{
                    1100.0,
                    2200.0},
        "fixed-BHP initial 2P phase-identity binding is malformed");

    const double initial_well_index =
        well_context.connection()
            .well_index_m3;
    const double initial_bhp =
        well_context
            .bottom_hole_pressure_pa();

    ControllerFixture fixture{
        rank,
        &schedule,
        &partition,
        &bridge,
        &pattern,
        audit,
        nullptr,
        provider,
        false,
        0U,
        false};
    fixture.well_context =
        &well_context;

    auto initial_system =
        make_controller_initial_system(
            &fixture);
    require_collective(
        initial_system != nullptr &&
            initial_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .cell_global ==
                mesh::GlobalEntityId{
                    UINT64_C(30)} &&
            initial_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                2U &&
            initial_system
                    ->coordinate_registry()
                    .cell(
                        mesh::GlobalEntityId{
                            UINT64_C(30)})
                    .active_phases
                    .phase_count() ==
                2U,
        "fixed-BHP transition fixture did not start on the stable 2P completion");

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        stable_system;
    Vec stable_state = nullptr;
    std::optional<
        fdp::
            PostSnesPhaseTransitionControllerReport3D>
        transition_report;
    PetscErrorCode error =
        fdp::
            solve_nonlinear_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                std::move(initial_system),
                {
                    &controller_scan,
                    &fixture,
                    &controller_rebuild,
                    &fixture},
                {4U},
                &stable_system,
                &stable_state,
                &transition_report);

    require_collective(
        error == PETSC_SUCCESS &&
            stable_system != nullptr &&
            stable_state != nullptr &&
            transition_report.has_value() &&
            transition_report->outcome ==
                fdp::
                    PostSnesPhaseTransitionOutcome3D::
                        stable_phase_set &&
            transition_report
                    ->transition_restarts ==
                1U &&
            transition_report
                    ->generations.size() ==
                2U &&
            fixture.rebuild_calls ==
                1U &&
            fixture.well_rebound &&
            fixture.well_context !=
                nullptr &&
            stable_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                3U &&
            stable_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .scalar_count ==
                10U,
        "fixed-BHP completion did not survive the accepted 2P-to-3P rebuild and SNES restart");

    const auto& rebound_context =
        *fixture.well_context;
    const auto& rebuilt_phase_map =
        stable_system
            ->coordinate_registry()
            .cell(
                mesh::GlobalEntityId{
                    UINT64_C(30)})
            .active_phases;
    bool rebound_identity_matches =
        rebound_context
                .active_phase_identities()
                .has_value() &&
        rebound_context
                .active_phase_identities()
                ->phase_count() ==
            rebuilt_phase_map.phase_count();
    if (rebound_identity_matches) {
        for (std::size_t phase = 0U;
             phase < rebuilt_phase_map.phase_count();
             ++phase) {
            rebound_identity_matches =
                rebound_identity_matches &&
                rebound_context
                        .active_phase_identities()
                        ->identity(phase) ==
                    rebuilt_phase_map.identity(
                        phase);
        }
    }

    require_collective(
        rebound_identity_matches &&
            rebound_context
                    .target_cell_global() ==
                mesh::GlobalEntityId{
                    UINT64_C(30)} &&
            std::abs(
                rebound_context
                        .connection()
                        .well_index_m3 -
                    initial_well_index) <=
                1.0e-14 *
                    std::max(
                        1.0,
                        std::abs(
                            initial_well_index)) &&
            rebound_context
                    .bottom_hole_pressure_pa() ==
                initial_bhp &&
            rebound_context
                    .injection_enthalpy()
                    .specific_enthalpy_j_per_kg ==
                std::vector<double>{
                    1100.0,
                    2200.0,
                    3300.0},
        "fixed-BHP completion rebinding changed stable cell/WI/BHP or failed stable phase-identity enthalpy mapping");

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        rebound_current;
    std::vector<double>
        rebound_porosity;
    fdp::NaturalVariableSnesEvaluationStatus3D
        rebound_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    require_collective(
        stable_system
                ->evaluate_local_cells_for_phase_transition(
                    stable_state,
                    &rebound_current,
                    &rebound_porosity,
                    &rebound_status) ==
            PETSC_SUCCESS &&
            rebound_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            rebound_current.size() >
                2U &&
            rebound_current[2U]
                .has_value(),
        "failed to evaluate rebuilt fixed-BHP completion state");

    const auto rebound_well =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                rebound_context,
                *rebound_current[2U]);
    require_collective(
        rebound_well
                .cell_source
                .input_count ==
            10U,
        "fixed-BHP rebound source did not rebuild its Jacobian on q=3Nc+1");

    // First prove the post-appearance 3P completion can run a real accepted
    // physical timestep with owner-only well source and conservation.
    run_rebound_fixed_bhp_physical_timestep(
        rank,
        &stable_system,
        &stable_state,
        &fixture,
        3U,
        10U);

    // Start a second, independent controller session from the already accepted
    // 3P system. This preserves the controller's cycle-detection contract: the
    // reverse 3P->2P transition is not a rollback inside the first session.
    fixture.phase_disappearance_only =
        true;
    fixture.oscillate_after_restart =
        false;
    fixture.rebuild_calls =
        0U;
    fixture.well_rebound =
        false;
    fixture.well_source_audit.evaluator_calls =
        0U;
    fixture.well_source_audit.target_calls =
        0U;

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        roundtrip_system;
    Vec roundtrip_state = nullptr;
    std::optional<
        fdp::
            PostSnesPhaseTransitionControllerReport3D>
        disappearance_report;
    error =
        fdp::
            solve_nonlinear_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                std::move(stable_system),
                {
                    &controller_scan,
                    &fixture,
                    &controller_rebuild,
                    &fixture},
                {4U},
                &roundtrip_system,
                &roundtrip_state,
                &disappearance_report);

    require_collective(
        error == PETSC_SUCCESS &&
            roundtrip_system != nullptr &&
            roundtrip_state != nullptr &&
            disappearance_report.has_value() &&
            disappearance_report->outcome ==
                fdp::
                    PostSnesPhaseTransitionOutcome3D::
                        stable_phase_set &&
            disappearance_report
                    ->transition_restarts ==
                1U &&
            disappearance_report
                    ->generations.size() ==
                2U &&
            fixture.rebuild_calls ==
                1U &&
            fixture.well_rebound &&
            fixture.well_context !=
                nullptr &&
            roundtrip_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .cell_global ==
                mesh::GlobalEntityId{
                    UINT64_C(30)} &&
            roundtrip_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                2U &&
            roundtrip_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .scalar_count ==
                7U &&
            roundtrip_system
                    ->numbering()
                    .petsc_global_scalar_count() ==
                42,
        "fixed-BHP completion did not survive the independent 3P-to-2P phase-disappearance restart");

    const auto initial_two_phase_map =
        phase_identity_map_for_stable(
            UINT64_C(30));
    const auto& roundtrip_phase_map =
        roundtrip_system
            ->coordinate_registry()
            .cell(
                mesh::GlobalEntityId{
                    UINT64_C(30)})
            .active_phases;
    bool roundtrip_identity_matches =
        roundtrip_phase_map.phase_count() ==
            initial_two_phase_map.phase_count() &&
        fixture.well_context
                ->active_phase_identities()
                .has_value() &&
        fixture.well_context
                ->active_phase_identities()
                ->phase_count() ==
            initial_two_phase_map.phase_count();
    if (roundtrip_identity_matches) {
        for (std::size_t phase = 0U;
             phase <
             initial_two_phase_map.phase_count();
             ++phase) {
            roundtrip_identity_matches =
                roundtrip_identity_matches &&
                roundtrip_phase_map.identity(
                    phase) ==
                    initial_two_phase_map.identity(
                        phase) &&
                fixture.well_context
                        ->active_phase_identities()
                        ->identity(phase) ==
                    initial_two_phase_map.identity(
                        phase);
        }
    }

    const auto& retained_registry =
        fixture.well_context
            ->phase_identity_injection_enthalpy();
    const bool retained_disappeared_phase =
        retained_registry.has_value() &&
        retained_registry->phases.size() ==
            3U &&
        std::any_of(
            retained_registry
                ->phases.begin(),
            retained_registry
                ->phases.end(),
            [&](const auto& entry) {
                return entry.identity ==
                           hydrocarbon1 &&
                    entry
                            .specific_enthalpy_j_per_kg ==
                        3300.0;
            });

    require_collective(
        roundtrip_identity_matches &&
            retained_disappeared_phase &&
            fixture.well_context
                    ->target_cell_global() ==
                mesh::GlobalEntityId{
                    UINT64_C(30)} &&
            std::abs(
                fixture.well_context
                        ->connection()
                        .well_index_m3 -
                    initial_well_index) <=
                1.0e-14 *
                    std::max(
                        1.0,
                        std::abs(
                            initial_well_index)) &&
            fixture.well_context
                    ->bottom_hole_pressure_pa() ==
                initial_bhp &&
            fixture.well_context
                    ->injection_enthalpy()
                    .specific_enthalpy_j_per_kg ==
                std::vector<double>{
                    1100.0,
                    2200.0},
        "fixed-BHP 2P->3P->2P round-trip drifted stable phase identity or discarded inactive enthalpy registry data");

    std::vector<std::optional<
        fdp::
            MixedCardinalityPhysicalCurrentCellLinearization3D>>
        roundtrip_current;
    std::vector<double>
        roundtrip_porosity;
    fdp::NaturalVariableSnesEvaluationStatus3D
        roundtrip_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    require_collective(
        roundtrip_system
                ->evaluate_local_cells_for_phase_transition(
                    roundtrip_state,
                    &roundtrip_current,
                    &roundtrip_porosity,
                    &roundtrip_status) ==
            PETSC_SUCCESS &&
            roundtrip_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            roundtrip_current.size() >
                2U &&
            roundtrip_current[2U]
                .has_value(),
        "failed to evaluate round-trip 2P fixed-BHP completion state");

    const auto roundtrip_well =
        wdp::
            build_fixed_bhp_peaceman_well_source_3d(
                *fixture.well_context,
                *roundtrip_current[2U]);
    require_collective(
        roundtrip_well
                .cell_source
                .input_count ==
            7U,
        "phase-disappearance rebinding did not shrink fixed-BHP Jacobian back to q=2Nc+1");

    run_rebound_fixed_bhp_physical_timestep(
        rank,
        &roundtrip_system,
        &roundtrip_state,
        &fixture,
        2U,
        7U);

    audit->well_timestep_compressibility =
        false;
}

} // namespace

void mixed_cardinality_physical_snes_assembly_test() {
    int rank = -1;
    int size = -1;
    if (MPI_Comm_rank(
            PETSC_COMM_WORLD,
            &rank) != MPI_SUCCESS ||
        MPI_Comm_size(
            PETSC_COMM_WORLD,
            &size) != MPI_SUCCESS) {
        throw std::runtime_error(
            "failed to query MPI rank/size for mixed physical dispatcher test");
    }
    require_collective(
        size == 2,
        "mixed physical dispatcher regression requires two ranks");

    const auto partition =
        make_partition(rank);
    const std::vector<std::size_t>
        phase_counts{
            1U, 1U, 2U, 2U, 3U, 3U};
    std::optional<
        fdp::
            VariableCardinalityNaturalVariableNumbering3D>
        numbering;
    PetscErrorCode error =
        fdp::
            make_variable_cardinality_natural_variable_numbering_3d(
                PETSC_COMM_WORLD,
                partition,
                3U,
                phase_counts,
                &numbering);
    require_collective(
        error == PETSC_SUCCESS &&
            numbering.has_value() &&
            numbering->petsc_global_scalar_count() ==
                42,
        "failed to build 1P/2P/3P mixed physical numbering");

    const auto bridge =
        make_cell_bridge(rank);
    const auto pattern =
        make_cell_pattern(rank);
    const auto schedule =
        make_schedule(
            rank,
            false);

    DispatchAudit audit;

    auto fixed_bhp_well_context =
        wdp::
            FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D::
                create(
                    mesh::GlobalEntityId{
                        UINT64_C(60)},
                    well::make_peaceman_well_index_3d(
                        {10.0, 10.0, 5.0},
                        {
                            1.0e-12,
                            1.0e-12,
                            1.0e-12},
                        well::
                            AxisAlignedWellDirection3D::z,
                        0.10,
                        0.0),
                    25.0,
                    wd::InjectionPhaseSpecificEnthalpy3P{
                        "fixture/fixed-bhp-injection-enthalpy/v1",
                        {1000.0, 2000.0, 3000.0}},
                    "fixture/fixed-bhp-cell60-source/v1");

    // The first bridge intentionally supports only a frozen 3P target. A 1P
    // target is rejected explicitly rather than padded with fictitious phases.
    {
        auto unsupported_context =
            wdp::
                FixedBhpThreePhasePeacemanWellSourceEvaluatorContext3D::
                    create(
                        mesh::GlobalEntityId{
                            UINT64_C(20)},
                        fixed_bhp_well_context
                            .connection(),
                        5.0,
                        wd::InjectionPhaseSpecificEnthalpy3P{
                            "fixture/fixed-bhp-unsupported/v1",
                            {1000.0, 2000.0, 3000.0}},
                        "fixture/fixed-bhp-unsupported-source/v1");
        auto current_1p =
            evaluate_target(
                UINT64_C(20),
                &audit);
        std::optional<
            fd::CellSourceLinearization3D>
            unsupported_source;
        fdp::NaturalVariableSnesEvaluationStatus3D
            unsupported_status =
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success;
        const auto q1 =
            target_1p();
        const PetscErrorCode unsupported_error =
            wdp::
                evaluate_fixed_bhp_three_phase_peaceman_well_source_3d(
                    mesh::LocalIndex{1U},
                    mesh::GlobalEntityId{
                        UINT64_C(20)},
                    q1,
                    current_1p,
                    &unsupported_context,
                    &unsupported_source,
                    &unsupported_status);
        require_collective(
            unsupported_error ==
                    PETSC_ERR_SUP &&
                !unsupported_source.has_value(),
            "fixed-BHP well bridge fabricated inactive phases for a non-3P target");
    }

    auto cell_inputs =
        make_cell_inputs(
            rank,
            &audit);
    auto face_inputs =
        make_face_inputs(
            rank,
            false);

    std::optional<
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D>
        context;
    error =
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D::
                create(
                    PETSC_COMM_WORLD,
                    schedule,
                    partition,
                    *numbering,
                    bridge,
                    pattern,
                    1.0,
                    std::move(cell_inputs),
                    make_phase_identity_maps(),
                    std::move(face_inputs),
                    {
                        {&evaluate_1p, &audit},
                        {&evaluate_2p, &audit},
                        {&evaluate_3p, &audit}},
                    {},
                    {
                        &evaluate_explicit_cell_source,
                        &audit},
                    &context);
    require_collective(
        error == PETSC_SUCCESS &&
            context.has_value(),
        "failed to create mixed physical production dispatcher");

    Mat jacobian = nullptr;
    error =
        context->create_jacobian_structure(
            &jacobian);
    require_collective(
        error == PETSC_SUCCESS &&
            jacobian != nullptr,
        "failed to create mixed physical ragged MPIAIJ structure");

    Vec state = nullptr;
    error =
        fdp::
            create_variable_cardinality_natural_variable_vec_3d(
                PETSC_COMM_WORLD,
                *numbering,
                &state);
    require_collective(
        error == PETSC_SUCCESS &&
            state != nullptr,
        "failed to create mixed physical state vector");
    insert_target(
        state,
        *numbering);

    Vec residual = nullptr;
    error =
        VecDuplicate(
            state,
            &residual);
    if (error == PETSC_SUCCESS) {
        error =
            VecSet(
                residual,
                PetscScalar{0.0});
    }
    auto evaluator =
        context->snes_evaluator();
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    if (error == PETSC_SUCCESS) {
        error =
            evaluator.function(
                state,
                residual,
                evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyBegin(
                residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyEnd(
                residual);
    }
    PetscReal residual_norm = -1.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                residual,
                NORM_2,
                &residual_norm);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            residual_norm < 1.0e-12,
        "target mixed physical residual is not zero");

    error =
        MatZeroEntries(
            jacobian);
    status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    if (error == PETSC_SUCCESS) {
        error =
            evaluator.jacobian(
                state,
                jacobian,
                evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyBegin(
                jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success,
        "failed to assemble mixed physical Jacobian");

    bool local_spatial_blocks_ok = true;
    if (rank == 0) {
        const std::array<
            std::pair<PetscInt, PetscInt>,
            3>
            probes{{
                {0, 21},
                {4, 25},
                {11, 32}}};
        for (const auto& [row, column] :
             probes) {
            PetscScalar value = 0.0;
            const auto get_error =
                MatGetValues(
                    jacobian,
                    1,
                    &row,
                    1,
                    &column,
                    &value);
            local_spatial_blocks_ok =
                local_spatial_blocks_ok &&
                get_error ==
                    PETSC_SUCCESS &&
                std::isfinite(
                    static_cast<double>(
                        PetscRealPart(
                            value))) &&
                std::abs(
                    static_cast<double>(
                        PetscRealPart(
                            value))) >
                    0.0;
        }
    }
    require_collective(
        local_spatial_blocks_ok,
        "1P/2P/3P real TPFA off-diagonal blocks were not assembled");

    require_collective(
        audit.single_calls > 0U &&
            audit.two_calls > 0U &&
            audit.three_calls > 0U,
        "mixed physical dispatcher did not exercise all three cell closures");

    PetscScalar baseline_component_d = 0.0;
    PetscScalar baseline_energy_d = 0.0;
    if (rank == 1) {
        const auto& source_record =
            numbering->cell(
                mesh::LocalIndex{1U});
        const PetscInt row0 =
            source_record.petsc_global_scalar_start;
        const PetscInt energy_row =
            row0 + 3;
        const PetscInt column = row0;
        require_collective(
            MatGetValues(
                jacobian,
                1,
                &row0,
                1,
                &column,
                &baseline_component_d) ==
                PETSC_SUCCESS &&
            MatGetValues(
                jacobian,
                1,
                &energy_row,
                1,
                &column,
                &baseline_energy_d) ==
                PETSC_SUCCESS,
            "failed to capture source-disabled diagonal Jacobian");
    } else {
        require_collective(
            true,
            "source-disabled diagonal Jacobian");
    }

    PetscScalar baseline_well_component_d = 0.0;
    PetscScalar baseline_well_energy_d = 0.0;
    if (rank == 1) {
        const auto& well_record =
            numbering->cell(
                mesh::LocalIndex{5U});
        const PetscInt row0 =
            well_record.petsc_global_scalar_start;
        const PetscInt energy_row =
            row0 + 3;
        const PetscInt column =
            row0;
        require_collective(
            MatGetValues(
                jacobian,
                1,
                &row0,
                1,
                &column,
                &baseline_well_component_d) ==
                PETSC_SUCCESS &&
            MatGetValues(
                jacobian,
                1,
                &energy_row,
                1,
                &column,
                &baseline_well_energy_d) ==
                PETSC_SUCCESS,
            "failed to capture fixed-BHP well target baseline Jacobian");
    } else {
        require_collective(
            true,
            "fixed-BHP well target baseline Jacobian");
    }

    audit.source_enabled = true;
    error = VecSet(
        residual,
        PetscScalar{0.0});
    status =
        fdp::NaturalVariableSnesEvaluationStatus3D::
            success;
    if (error == PETSC_SUCCESS) {
        error =
            evaluator.function(
                state,
                residual,
                evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error = VecAssemblyBegin(residual);
    }
    if (error == PETSC_SUCCESS) {
        error = VecAssemblyEnd(residual);
    }
    if (error == PETSC_SUCCESS) {
        error = MatZeroEntries(jacobian);
    }
    if (error == PETSC_SUCCESS) {
        error =
            evaluator.jacobian(
                state,
                jacobian,
                evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyBegin(
                jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                jacobian,
                MAT_FINAL_ASSEMBLY);
    }

    bool local_source_ok = true;
    if (rank == 1) {
        const auto& source_record =
            numbering->cell(
                mesh::LocalIndex{1U});
        const PetscInt row0 =
            source_record.petsc_global_scalar_start;
        const PetscInt energy_row =
            row0 + 3;
        const PetscInt column = row0;
        std::array<PetscInt, 4> rows{
            row0,
            row0 + 1,
            row0 + 2,
            energy_row};
        std::array<PetscScalar, 4> values{};
        local_source_ok =
            VecGetValues(
                residual,
                static_cast<PetscInt>(rows.size()),
                rows.data(),
                values.data()) ==
            PETSC_SUCCESS;
        const std::array<double, 4> expected{
            -2.0 / 3.0,
            1.0 / 3.0,
            -0.5 / 3.0,
            -20.0 / 3.0};
        for (std::size_t i = 0U;
             i < expected.size();
             ++i) {
            local_source_ok =
                local_source_ok &&
                std::abs(
                    static_cast<double>(
                        PetscRealPart(values[i])) -
                    expected[i]) <
                    1.0e-11;
        }

        PetscScalar sourced_component_d = 0.0;
        PetscScalar sourced_energy_d = 0.0;
        local_source_ok =
            local_source_ok &&
            MatGetValues(
                jacobian,
                1,
                &row0,
                1,
                &column,
                &sourced_component_d) ==
                PETSC_SUCCESS &&
            MatGetValues(
                jacobian,
                1,
                &energy_row,
                1,
                &column,
                &sourced_energy_d) ==
                PETSC_SUCCESS &&
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        sourced_component_d -
                        baseline_component_d)) +
                0.1 / 3.0) <
                1.0e-11 &&
            std::abs(
                static_cast<double>(
                    PetscRealPart(
                        sourced_energy_d -
                        baseline_energy_d)) +
                0.2 / 3.0) <
                1.0e-11;
    }
    require_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    success &&
            local_source_ok &&
            (rank == 0
                 ? audit.source_cell20_calls == 0U
                 : audit.source_cell20_calls > 0U),
        "explicit cell source did not map to owner component/energy residual and diagonal Jacobian");

    // The actual fixed-BHP well bridge now enters the same owner-only source
    // callback and therefore the fully implicit mixed-cardinality residual and
    // diagonal Jacobian. The connection is attached to stable cell60, owned by
    // rank1 and ghosted on rank0.
    FixedBhpWellSourceAudit
        fixed_bhp_source_audit{
            &fixed_bhp_well_context};
    auto well_cells =
        make_cell_inputs(
            rank,
            &audit);
    auto well_faces =
        make_face_inputs(
            rank,
            false);
    std::optional<
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D>
        well_assembly_context;
    error =
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D::
                create(
                    PETSC_COMM_WORLD,
                    schedule,
                    partition,
                    *numbering,
                    bridge,
                    pattern,
                    1.0,
                    std::move(well_cells),
                    make_phase_identity_maps(),
                    std::move(well_faces),
                    {
                        {&evaluate_1p, &audit},
                        {&evaluate_2p, &audit},
                        {&evaluate_3p, &audit}},
                    {},
                    {
                        &evaluate_audited_fixed_bhp_well_source,
                        &fixed_bhp_source_audit},
                    &well_assembly_context);
    require_collective(
        error == PETSC_SUCCESS &&
            well_assembly_context.has_value(),
        "failed to create fixed-BHP well mixed physical assembly context");

    auto direct_current =
        evaluate_target(
            UINT64_C(60),
            &audit);
    const auto& direct_three_phase =
        std::get<
            fdp::
                FixedThreePhaseCurrentCellLinearization3D>(
                    direct_current);
    const auto expected_well_source =
        wdp::
            build_fixed_bhp_three_phase_peaceman_well_source_3d(
                fixed_bhp_well_context,
                direct_three_phase);
    const auto expected_well_normalized =
        fd::normalize_cell_source_by_bulk_volume(
            expected_well_source.cell_source,
            7.0);

    Mat well_jacobian = nullptr;
    Vec well_residual = nullptr;
    error =
        well_assembly_context
            ->create_jacobian_structure(
                &well_jacobian);
    if (error == PETSC_SUCCESS) {
        error =
            VecDuplicate(
                state,
                &well_residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecSet(
                well_residual,
                PetscScalar{0.0});
    }
    auto well_evaluator =
        well_assembly_context
            ->snes_evaluator();
    status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    if (error == PETSC_SUCCESS) {
        error =
            well_evaluator.function(
                state,
                well_residual,
                well_evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyBegin(
                well_residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyEnd(
                well_residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatZeroEntries(
                well_jacobian);
    }
    if (error == PETSC_SUCCESS) {
        error =
            well_evaluator.jacobian(
                state,
                well_jacobian,
                well_evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyBegin(
                well_jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                well_jacobian,
                MAT_FINAL_ASSEMBLY);
    }

    bool local_well_source_ok =
        error == PETSC_SUCCESS &&
        status ==
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    if (rank == 1 &&
        local_well_source_ok) {
        const auto& well_record =
            numbering->cell(
                mesh::LocalIndex{5U});
        const PetscInt row0 =
            well_record.petsc_global_scalar_start;
        const PetscInt energy_row =
            row0 + 3;
        const PetscInt column =
            row0;
        std::array<PetscInt, 4>
            rows{
                row0,
                row0 + 1,
                row0 + 2,
                energy_row};
        std::array<PetscScalar, 4>
            values{};
        local_well_source_ok =
            VecGetValues(
                well_residual,
                static_cast<PetscInt>(
                    rows.size()),
                rows.data(),
                values.data()) ==
            PETSC_SUCCESS;

        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            const double actual =
                static_cast<double>(
                    PetscRealPart(
                        values[component]));
            const double expected =
                expected_well_normalized
                    .component_residual_mol_per_bulk_m3_s[
                        component];
            local_well_source_ok =
                local_well_source_ok &&
                std::abs(
                    actual -
                    expected) <=
                    1.0e-10 *
                        std::max(
                            1.0e-30,
                            std::abs(expected)) +
                    1.0e-20;
        }
        {
            const double actual =
                static_cast<double>(
                    PetscRealPart(
                        values[3]));
            const double expected =
                expected_well_normalized
                    .energy_residual_w_per_bulk_m3;
            local_well_source_ok =
                local_well_source_ok &&
                std::abs(
                    actual -
                    expected) <=
                    1.0e-10 *
                        std::max(
                            1.0e-30,
                            std::abs(expected)) +
                    1.0e-20;
        }

        PetscScalar well_component_d = 0.0;
        PetscScalar well_energy_d = 0.0;
        local_well_source_ok =
            local_well_source_ok &&
            MatGetValues(
                well_jacobian,
                1,
                &row0,
                1,
                &column,
                &well_component_d) ==
                PETSC_SUCCESS &&
            MatGetValues(
                well_jacobian,
                1,
                &energy_row,
                1,
                &column,
                &well_energy_d) ==
                PETSC_SUCCESS;

        const double component_delta =
            static_cast<double>(
                PetscRealPart(
                    well_component_d -
                    baseline_well_component_d));
        const double energy_delta =
            static_cast<double>(
                PetscRealPart(
                    well_energy_d -
                    baseline_well_energy_d));
        const double expected_component_delta =
            expected_well_normalized
                .d_component_residual(
                    0U,
                    0U);
        const double expected_energy_delta =
            expected_well_normalized
                .d_energy_residual(
                    0U);
        local_well_source_ok =
            local_well_source_ok &&
            std::abs(
                component_delta -
                expected_component_delta) <=
                1.0e-10 *
                    std::max(
                        1.0e-30,
                        std::abs(
                            expected_component_delta)) +
                1.0e-20 &&
            std::abs(
                energy_delta -
                expected_energy_delta) <=
                1.0e-10 *
                    std::max(
                        1.0e-30,
                        std::abs(
                            expected_energy_delta)) +
                1.0e-20;
    }

    require_collective(
        local_well_source_ok &&
            (rank == 0
                 ? fixed_bhp_source_audit
                           .target_calls ==
                       0U
                 : fixed_bhp_source_audit
                           .target_calls ==
                       2U),
        "fixed-BHP production well source did not enter owner-only residual/Jacobian correctly");

    require_collective(
        (well_residual == nullptr ||
         VecDestroy(
             &well_residual) ==
             PETSC_SUCCESS) &&
            (well_jacobian == nullptr ||
             MatDestroy(
                 &well_jacobian) ==
                 PETSC_SUCCESS),
        "fixed-BHP well source assembly cleanup failed");

    audit.source_enabled = false;


    // Exercise the same production fixed-BHP source/driver/conservation path
    // on one frozen target cell of each supported phase cardinality. The mixed
    // reservoir around the completion is unchanged; only the well target chart
    // and its exact active-phase injection-enthalpy payload differ.
    run_frozen_fixed_bhp_timestep_case(
        rank,
        UINT64_C(20),
        1U,
        1U,
        5.0,
        schedule,
        partition,
        bridge,
        pattern,
        &audit);
    run_frozen_fixed_bhp_timestep_case(
        rank,
        UINT64_C(40),
        3U,
        2U,
        15.0,
        schedule,
        partition,
        bridge,
        pattern,
        &audit);
    run_frozen_fixed_bhp_timestep_case(
        rank,
        UINT64_C(60),
        5U,
        3U,
        25.0,
        schedule,
        partition,
        bridge,
        pattern,
        &audit);

    run_multi_connection_fixed_bhp_timestep_case(
        rank,
        schedule,
        partition,
        bridge,
        pattern,
        &audit);

    run_fixed_total_molar_rate_control_case(
        rank,
        schedule,
        partition,
        bridge,
        pattern,
        &audit);

    Vec solution = nullptr;
    std::optional<
        fdp::
            VariableCardinalityNaturalVariableSnesSolveReport3D>
        report;
    error =
        fdp::
            solve_variable_cardinality_natural_variable_snes_3d(
                PETSC_COMM_WORLD,
                *numbering,
                state,
                jacobian,
                context->snes_evaluator(),
                &solution,
                &report);
    require_collective(
        error == PETSC_SUCCESS &&
            solution != nullptr &&
            report.has_value() &&
            static_cast<int>(
                report->converged_reason) >
                0 &&
            report->snes_type ==
                SNESNEWTONLS &&
            report->ksp_type ==
                KSPGMRES &&
            report->pc_type ==
                PCASM &&
            local_solution_matches_target(
                *report),
        "mixed physical system did not pass PETSc SNES production solve");

    // A cross-cardinality face must never be guessed from the fixed-cardinality
    // TPFA kernels. Without an explicit bridge evaluator, construction is
    // collectively rejected before nonlinear callbacks are installed.
    const auto cross_schedule =
        make_schedule(
            rank,
            true);
    auto cross_cells =
        make_cell_inputs(
            rank,
            &audit);
    auto cross_faces =
        make_face_inputs(
            rank,
            true);
    std::optional<
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D>
        rejected;
    error =
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D::
                create(
                    PETSC_COMM_WORLD,
                    cross_schedule,
                    partition,
                    *numbering,
                    bridge,
                    pattern,
                    1.0,
                    std::move(cross_cells),
                    make_phase_identity_maps(),
                    std::move(cross_faces),
                    {
                        {&evaluate_1p, &audit},
                        {&evaluate_2p, &audit},
                        {&evaluate_3p, &audit}},
                    {},
                    &rejected);
    require_collective(
        error == PETSC_ERR_SUP &&
            !rejected.has_value(),
        "cross-cardinality TPFA face was not rejected without an explicit bridge");

    // Explicit potential extensions enable the standard physical bridge for
    // all three mixed-cardinality pair classes. The bridge computes actual
    // TPFA component/energy flux from active upwind payloads and never requests
    // absent-side composition or enthalpy.
    check_direct_cross_cardinality_bridge(
        UINT64_C(20),
        UINT64_C(30),
        3U,
        &audit);
    check_direct_cross_cardinality_bridge(
        UINT64_C(20),
        UINT64_C(50),
        4U,
        &audit);
    check_direct_cross_cardinality_bridge(
        UINT64_C(40),
        UINT64_C(50),
        5U,
        &audit);

    // The 1P<->2P face is now admitted through the standard bridge into the
    // real distributed dispatcher. Rank 0 owns the authoritative face while
    // cell20 is owned by rank 1, so this exercises off-process residual and
    // rectangular Jacobian insertion, not just a direct local bridge call.
    auto bridged_cells =
        make_cell_inputs(
            rank,
            &audit);
    auto bridged_faces =
        make_face_inputs(
            rank,
            true);
    const auto pr_parameters =
        thermodynamic_adapter_pr_parameters();
    const auto pr_model =
        th::Pr76Phase<double>::
            from_parameters(
                pr_parameters);
    auto pr_provider =
        make_thermodynamic_adapter_pr_provider(
            pr_model);
    ThermodynamicCoordinateResolverAudit
        coordinate_audit;
    fdp::
        ThermodynamicAbsentPhaseExtensionAdapterBinding3D
        thermodynamic_adapter{
            &resolve_absent_phase_thermodynamic_coordinates,
            &coordinate_audit,
            &fdp::
                evaluate_absent_phase_thermodynamic_provider_3d<
                    flow::
                        Pr76AbsentPhasePotentialExtensionProvider<
                            double>>,
            &pr_provider};
    fdp::CrossCardinalityTpfaBridgeBinding3D
        standard_bridge{
            &fdp::
                evaluate_thermodynamic_absent_phase_extension_3d,
            &thermodynamic_adapter};
    std::optional<
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D>
        bridged_context;
    error =
        fdp::
            MixedCardinalityPhysicalSnesAssemblyContext3D::
                create(
                    PETSC_COMM_WORLD,
                    cross_schedule,
                    partition,
                    *numbering,
                    bridge,
                    pattern,
                    1.0,
                    std::move(bridged_cells),
                    make_phase_identity_maps(),
                    std::move(bridged_faces),
                    {
                        {&evaluate_1p, &audit},
                        {&evaluate_2p, &audit},
                        {&evaluate_3p, &audit}},
                    {
                        &fdp::
                            evaluate_standard_cross_cardinality_tpfa_face_3d,
                        &standard_bridge},
                    &bridged_context);
    require_collective(
        error == PETSC_SUCCESS &&
            bridged_context.has_value(),
        "standard cross-cardinality bridge was not admitted by distributed dispatcher");

    Mat bridged_jacobian = nullptr;
    error =
        bridged_context
            ->create_jacobian_structure(
                &bridged_jacobian);
    require_collective(
        error == PETSC_SUCCESS &&
            bridged_jacobian != nullptr,
        "failed to create bridged cross-cardinality MPIAIJ structure");

    Vec bridged_state = nullptr;
    error =
        fdp::
            create_variable_cardinality_natural_variable_vec_3d(
                PETSC_COMM_WORLD,
                *numbering,
                &bridged_state);
    require_collective(
        error == PETSC_SUCCESS &&
            bridged_state != nullptr,
        "failed to create bridged cross-cardinality state");
    insert_target(
        bridged_state,
        *numbering);

    Vec bridged_residual = nullptr;
    error =
        VecDuplicate(
            bridged_state,
            &bridged_residual);
    if (error == PETSC_SUCCESS) {
        error =
            VecSet(
                bridged_residual,
                PetscScalar{0.0});
    }
    auto bridged_evaluator =
        bridged_context->snes_evaluator();
    status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    if (error == PETSC_SUCCESS) {
        error =
            bridged_evaluator.function(
                bridged_state,
                bridged_residual,
                bridged_evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyBegin(
                bridged_residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyEnd(
                bridged_residual);
    }
    PetscReal bridged_norm = 0.0;
    if (error == PETSC_SUCCESS) {
        error =
            VecNorm(
                bridged_residual,
                NORM_2,
                &bridged_norm);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            bridged_norm > 0.0,
        "distributed cross-cardinality face did not contribute a physical residual");

    error =
        MatZeroEntries(
            bridged_jacobian);
    status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    if (error == PETSC_SUCCESS) {
        error =
            bridged_evaluator.jacobian(
                bridged_state,
                bridged_jacobian,
                bridged_evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyBegin(
                bridged_jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                bridged_jacobian,
                MAT_FINAL_ASSEMBLY);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success,
        "distributed cross-cardinality Jacobian assembly failed");

    std::array<double, 4>
        local_cross_rate{};
    if (rank == 0) {
        const auto& cell30 =
            numbering->cell(
                mesh::LocalIndex{2U});
        for (std::size_t row = 0U;
             row < local_cross_rate.size();
             ++row) {
            const PetscInt index =
                cell30.petsc_global_scalar_start +
                static_cast<PetscInt>(row);
            PetscScalar value = 0.0;
            if (VecGetValues(
                    bridged_residual,
                    1,
                    &index,
                    &value) !=
                PETSC_SUCCESS) {
                throw std::runtime_error(
                    "failed to read cell30 bridged residual");
            }
            local_cross_rate[row] =
                static_cast<double>(
                    PetscRealPart(value)) *
                4.0;
        }
    } else {
        const auto& cell20 =
            numbering->cell(
                mesh::LocalIndex{1U});
        for (std::size_t row = 0U;
             row < local_cross_rate.size();
             ++row) {
            const PetscInt index =
                cell20.petsc_global_scalar_start +
                static_cast<PetscInt>(row);
            PetscScalar value = 0.0;
            if (VecGetValues(
                    bridged_residual,
                    1,
                    &index,
                    &value) !=
                PETSC_SUCCESS) {
                throw std::runtime_error(
                    "failed to read cell20 bridged residual");
            }
            local_cross_rate[row] =
                static_cast<double>(
                    PetscRealPart(value)) *
                3.0;
        }
    }

    std::array<double, 4>
        global_cross_rate{};
    require_collective(
        MPI_Allreduce(
            local_cross_rate.data(),
            global_cross_rate.data(),
            static_cast<int>(
                global_cross_rate.size()),
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce bridged cross-face conservation rows");
    bool conserved = true;
    for (const double value :
         global_cross_rate) {
        conserved =
            conserved &&
            std::abs(value) <
                1.0e-10;
    }
    require_collective(
        conserved,
        "distributed cross-cardinality component/energy face rows are not conservative");

    bool cross_block_nonzero = true;
    {
        const auto& row_cell =
            numbering->cell(
                rank == 0
                    ? mesh::LocalIndex{2U}
                    : mesh::LocalIndex{1U});
        const auto& column_cell =
            numbering->cell(
                rank == 0
                    ? mesh::LocalIndex{1U}
                    : mesh::LocalIndex{2U});
        const PetscInt row =
            row_cell
                .petsc_global_scalar_start;
        const PetscInt column =
            column_cell
                .petsc_global_scalar_start;
        PetscScalar value = 0.0;
        cross_block_nonzero =
            MatGetValues(
                bridged_jacobian,
                1,
                &row,
                1,
                &column,
                &value) ==
                PETSC_SUCCESS &&
            std::isfinite(
                static_cast<double>(
                    PetscRealPart(value))) &&
            std::abs(
                static_cast<double>(
                    PetscRealPart(value))) >
                0.0;
    }
    require_collective(
        cross_block_nonzero,
        "distributed 4x7/7x4 cross-cardinality Jacobian block is missing");
    std::uint64_t global_coordinate_calls = 0U;
    require_collective(
        MPI_Allreduce(
            &coordinate_audit.calls,
            &global_coordinate_calls,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            PETSC_COMM_WORLD) ==
            MPI_SUCCESS,
        "failed to reduce thermodynamic coordinate resolver call count");
    require_collective(
        global_coordinate_calls > 0U,
        "authoritative distributed bridge did not invoke the thermodynamic coordinate resolver");

    require_collective(
        VecDestroy(&bridged_residual) ==
                PETSC_SUCCESS &&
            VecDestroy(&bridged_state) ==
                PETSC_SUCCESS &&
            MatDestroy(&bridged_jacobian) ==
                PETSC_SUCCESS,
        "bridged cross-cardinality fixture cleanup failed");

    // Accepted phase transitions are applied only outside SNES.  Rebuild the
    // complete q-ragged production system after cell30 changes 2P -> 3P.
    // This must recreate global numbering/state, the PetscSF-backed physical
    // context and MPIAIJ structure while preserving BE histories.
    auto rebuild_cells =
        make_outer_rebuild_cells(
            rank,
            &audit);
    auto rebuild_faces =
        make_face_inputs(
            rank,
            true);
    auto outer_provider =
        make_outer_rebuild_pr_provider(
            pr_model);

    // A real fixed-BHP completion now follows stable cell30 through the
    // existing 2P->3P outer restart, resolves its well-side enthalpy by stable
    // phase identity, and then advances one post-rebuild physical timestep.
    run_phase_transition_rebound_fixed_bhp_case(
        rank,
        schedule,
        partition,
        bridge,
        pattern,
        &audit,
        &outer_provider);

    run_multi_connection_local_transition_rebind_case(
        rank,
        schedule,
        partition,
        bridge,
        pattern,
        &audit,
        &outer_provider);

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        rebuilt_system;
    error =
        fdp::
            rebuild_phase_transition_natural_variable_system_3d(
                PETSC_COMM_WORLD,
                cross_schedule,
                partition,
                bridge,
                pattern,
                1.0,
                std::move(rebuild_cells),
                std::move(rebuild_faces),
                {
                    {&evaluate_1p, &audit},
                    {&evaluate_2p, &audit},
                    {&evaluate_3p, &audit}},
                &fdp::
                    evaluate_absent_phase_thermodynamic_provider_3d<
                        flow::
                            Pr76AbsentPhasePotentialExtensionProvider<
                                double>>,
                &outer_provider,
                &rebuilt_system);
    require_collective(
        error == PETSC_SUCCESS &&
            rebuilt_system != nullptr &&
            rebuilt_system
                    ->numbering()
                    .petsc_global_scalar_count() ==
                45 &&
            rebuilt_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                3U &&
            rebuilt_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .scalar_count ==
                10U,
        "accepted 2P->3P transition did not rebuild q-ragged numbering");

    Vec rebuild_residual = nullptr;
    error =
        VecDuplicate(
            rebuilt_system
                ->initial_state(),
            &rebuild_residual);
    if (error == PETSC_SUCCESS) {
        error =
            VecSet(
                rebuild_residual,
                PetscScalar{0.0});
    }
    auto rebuild_evaluator =
        rebuilt_system
            ->snes_evaluator();
    status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    if (error == PETSC_SUCCESS) {
        error =
            rebuild_evaluator.function(
                rebuilt_system
                    ->initial_state(),
                rebuild_residual,
                rebuild_evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyBegin(
                rebuild_residual);
    }
    if (error == PETSC_SUCCESS) {
        error =
            VecAssemblyEnd(
                rebuild_residual);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success,
        "rebuilt phase-transition residual callback failed");

    error =
        MatZeroEntries(
            rebuilt_system
                ->jacobian_structure());
    status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    if (error == PETSC_SUCCESS) {
        error =
            rebuild_evaluator.jacobian(
                rebuilt_system
                    ->initial_state(),
                rebuilt_system
                    ->jacobian_structure(),
                rebuild_evaluator.user_context,
                &status);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyBegin(
                rebuilt_system
                    ->jacobian_structure(),
                MAT_FINAL_ASSEMBLY);
    }
    if (error == PETSC_SUCCESS) {
        error =
            MatAssemblyEnd(
                rebuilt_system
                    ->jacobian_structure(),
                MAT_FINAL_ASSEMBLY);
    }
    require_collective(
        error == PETSC_SUCCESS &&
            status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success,
        "rebuilt phase-transition Jacobian callback failed");

    require_collective(
        rebuilt_system
                ->coordinate_registry()
                .reference_entry(
                    mesh::GlobalEntityId{
                        UINT64_C(20)},
                    mixed_physical_phase_identity(
                        "hydrocarbon-1"))
                .selected_branch_provenance ==
            "fixture/hydrocarbon1-root0",
        "outer rebuild did not retain frozen absent-phase branch provenance");

    require_collective(
        VecDestroy(&rebuild_residual) ==
            PETSC_SUCCESS,
        "phase-transition rebuild residual cleanup failed");
    rebuilt_system.reset();

    // Full post-SNES controller: stable after one accepted restart.
    ControllerFixture stable_controller{
        rank,
        &schedule,
        &partition,
        &bridge,
        &pattern,
        &audit,
        &pr_model,
        &outer_provider,
        false,
        0U};
    run_controller_case(
        &stable_controller,
        4U,
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                stable_phase_set,
        1U);

    // Preserved proposal handoff: the first outer-controller scan is not
    // recomputed. It replays the exact local-owned batch captured by the
    // preceding fixed-cardinality solve, then delegates subsequent generations
    // to the normal controller scanner.
    stable_controller.rebuild_calls = 0U;
    auto handoff_initial =
        make_controller_initial_system(
            &stable_controller);
    fdp::PostSnesPhaseTransitionHandoffScannerContext3D
        handoff_scanner;
    if (rank == 0) {
        handoff_scanner
            .initial_local_owned_proposals
            .push_back(
                controller_two_to_three_proposal());
    }
    handoff_scanner.subsequent_scanner =
        &controller_scan;
    handoff_scanner.subsequent_scanner_context =
        &stable_controller;

    std::unique_ptr<
        fdp::
            PhaseTransitionRebuiltNaturalVariableSystem3D>
        handoff_final_system;
    Vec handoff_final_state = nullptr;
    std::optional<
        fdp::
            PostSnesPhaseTransitionControllerReport3D>
        handoff_controller_report;
    error =
        fdp::
            solve_nonlinear_timestep_with_phase_transitions_3d(
                PETSC_COMM_WORLD,
                std::move(
                    handoff_initial),
                {
                    &fdp::
                        scan_post_snes_phase_transition_handoff_3d,
                    &handoff_scanner,
                    &controller_rebuild,
                    &stable_controller},
                {4U},
                &handoff_final_system,
                &handoff_final_state,
                &handoff_controller_report);
    require_collective(
        error == PETSC_SUCCESS &&
            handoff_scanner.initial_batch_consumed &&
            handoff_controller_report.has_value() &&
            handoff_controller_report->outcome ==
                fdp::
                    PostSnesPhaseTransitionOutcome3D::
                        stable_phase_set &&
            handoff_controller_report
                    ->transition_restarts ==
                1U &&
            handoff_controller_report
                    ->generations.size() ==
                2U &&
            handoff_controller_report
                    ->generations.front()
                    .accepted_transition_batch
                    .size() ==
                1U &&
            stable_controller.rebuild_calls ==
                1U &&
            handoff_final_system != nullptr &&
            handoff_final_system
                    ->numbering()
                    .cell(
                        mesh::LocalIndex{2U})
                    .phase_count ==
                3U,
        "preserved proposal did not drive one outer-controller topology restart");
    require_collective(
        VecDestroy(
            &handoff_final_state) ==
            PETSC_SUCCESS,
        "handoff outer-controller final state cleanup failed");

    // Reverse 3P->2P proposal returns to the previously visited global
    // signature and must be rejected as a phase-set cycle before rebuilding.
    ControllerFixture cycle_controller{
        rank,
        &schedule,
        &partition,
        &bridge,
        &pattern,
        &audit,
        &pr_model,
        &outer_provider,
        true,
        0U};
    run_controller_case(
        &cycle_controller,
        4U,
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                phase_set_cycle_detected,
        1U);

    // A zero restart budget permits the converged scan but never mutates the
    // system; the timestep therefore remains unaccepted.
    ControllerFixture budget_controller{
        rank,
        &schedule,
        &partition,
        &bridge,
        &pattern,
        &audit,
        &pr_model,
        &outer_provider,
        false,
        0U};
    run_controller_case(
        &budget_controller,
        0U,
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                transition_restart_budget_exhausted,
        0U);

    ControllerFixture indeterminate_controller{
        rank,
        &schedule,
        &partition,
        &bridge,
        &pattern,
        &audit,
        &pr_model,
        &outer_provider,
        false,
        0U};
    indeterminate_controller.scan_indeterminate =
        true;
    run_controller_case(
        &indeterminate_controller,
        4U,
        fdp::
            PostSnesPhaseTransitionOutcome3D::
                phase_set_scan_indeterminate,
        0U);

    require_collective(
        VecDestroy(&solution) ==
                PETSC_SUCCESS &&
            VecDestroy(&residual) ==
                PETSC_SUCCESS &&
            VecDestroy(&state) ==
                PETSC_SUCCESS &&
            MatDestroy(&jacobian) ==
                PETSC_SUCCESS,
        "mixed physical dispatcher fixture cleanup failed");
}
