#include <mpmc/flow_discretization_petsc/physical_timestep_driver.hpp>
#include <mpmc/flow_discretization_petsc/post_snes_phase_transition_controller.hpp>
#include <mpmc/flow_discretization_petsc/post_snes_phase_transition_handoff_scanner.hpp>
#include <mpmc/well_discretization_petsc/fixed_bhp_well_source_evaluator.hpp>
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
        input.phase_properties =
            flow::PhasePropertyPrerequisiteInput{
                5.0,
                2.0,
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
        auto molar =
            flow::
                make_single_phase_molar_density_linearization(
                    state,
                    zero);
        const flow::TransportPropertyProvenance
            provenance{
                "mixed-physical-dispatch",
                "controlled-regression",
                "v1"};
        auto transport =
            flow::
                make_single_phase_transport_linearization(
                    state,
                    zero,
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

        input.phase_properties[0] =
            flow::PhasePropertyPrerequisiteInput{
                6.0,
                2.0,
                1.0,
                temperature,
                0.5 * temperature};
        input.phase_properties[1] =
            flow::PhasePropertyPrerequisiteInput{
                3.0,
                1.5,
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
        auto molar =
            flow::
                make_two_phase_molar_density_linearization(
                    state,
                    zero);
        const flow::TransportPropertyProvenance
            provenance{
                "mixed-physical-dispatch",
                "controlled-regression",
                "v1"};
        auto transport =
            flow::
                make_two_phase_transport_linearization(
                    state,
                    zero,
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

        const std::array<double, 3>
            molar_density{
                2.0, 4.0, 7.0};
        const std::array<double, 3>
            mass_density{
                1.0, 2.0, 3.0};
        const std::array<double, 3>
            viscosity{
                1.0, 1.2, 1.4};
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
            dh;
        std::array<std::vector<double>, 3>
            du;
        for (std::size_t phase = 0U;
             phase < 3U;
             ++phase) {
            zero[phase].assign(
                q,
                0.0);
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
                zero};
        const flow::TransportPropertyProvenance
            provenance{
                "mixed-physical-dispatch",
                "controlled-regression",
                "v1"};
        auto transport =
            flow::
                make_phase_transport_property_linearization(
                    state,
                    zero,
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
        output->push_back(
            controller_two_to_three_proposal());
        return PETSC_SUCCESS;
    }

    if (record.phase_count == 3U &&
        fixture
            ->oscillate_after_restart) {
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
            global_well_production_rate[
                component] >
                0.0,
            "variable-cardinality fixed-BHP regression did not remain production-positive");
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
        global_well_production_rate[3] >
            0.0,
        "variable-cardinality fixed-BHP energy rate did not remain production-positive");
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
        1.0e-12,
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
        1.0e-12,
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
        1.0e-12,
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
