#include <mpmc/flow_discretization_petsc/fixed_three_phase_snes_assembly.hpp>
#include <mpmc/flow_discretization_petsc/complete_natural_variable_petsc_materialization.hpp>

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace mesh = mpmc::mesh;
namespace flow = mpmc::flow;
namespace disc = mpmc::discretization;
namespace dp = mpmc::discretization_petsc;
namespace fdp = mpmc::flow_discretization_petsc;

void require_collective(
    bool local_condition,
    std::string_view message) {
    int local =
        local_condition ? 1 : 0;
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
            "MPI_Allreduce failed in fixed-three-phase SNES regression");
    }
    if (global == 0) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_collective(
    double actual,
    double expected,
    double relative = 2.0e-8,
    double absolute = 2.0e-9) {
    const bool local =
        std::isfinite(actual) &&
        std::isfinite(expected) &&
        std::abs(actual - expected) <=
            absolute +
                relative *
                    std::max(
                        std::abs(actual),
                        std::abs(expected));
    require_collective(
        local,
        "numeric mismatch in fixed-three-phase SNES regression");
}

flow::NaturalVariableLayout3P
target_layout(std::uint64_t stable_cell) {
    if (stable_cell == UINT64_C(10)) {
        return flow::NaturalVariableLayout3P{
            flow::NaturalVariableCompositionPivot3P::
                from_dependent_components(
                    3U,
                    {1U, 0U, 2U})};
    }
    if (stable_cell == UINT64_C(20)) {
        return flow::NaturalVariableLayout3P{
            flow::NaturalVariableCompositionPivot3P::
                from_dependent_components(
                    3U,
                    {0U, 2U, 1U})};
    }
    throw std::invalid_argument(
        "unknown controlled stable cell");
}

std::array<std::vector<double>, 3>
target_full_composition(
    std::uint64_t stable_cell) {
    if (stable_cell == UINT64_C(10)) {
        return {
            std::vector<double>{0.10, 0.70, 0.20},
            std::vector<double>{0.60, 0.20, 0.20},
            std::vector<double>{0.20, 0.30, 0.50}};
    }
    if (stable_cell == UINT64_C(20)) {
        return {
            std::vector<double>{0.55, 0.25, 0.20},
            std::vector<double>{0.25, 0.25, 0.50},
            std::vector<double>{0.20, 0.60, 0.20}};
    }
    throw std::invalid_argument(
        "unknown controlled stable cell");
}

std::vector<double>
target_natural_variables(
    std::uint64_t stable_cell) {
    const auto layout =
        target_layout(
            stable_cell);
    const auto composition =
        target_full_composition(
            stable_cell);
    std::vector<double> q(
        layout.unknown_count(),
        0.0);

    q[layout.pressure_unknown_index()] =
        stable_cell == UINT64_C(10)
            ? 10.0
            : 12.0;
    q[layout.temperature_unknown_index()] =
        stable_cell == UINT64_C(10)
            ? 8.0
            : 9.0;

    const auto s0 =
        layout.independent_saturation_unknown_index(
            flow::PhaseSlot3::phase0);
    const auto s1 =
        layout.independent_saturation_unknown_index(
            flow::PhaseSlot3::phase1);
    if (!s0 || !s1) {
        throw std::logic_error(
            "controlled saturation layout malformed");
    }
    q[*s0] =
        stable_cell == UINT64_C(10)
            ? 0.20
            : 0.25;
    q[*s1] =
        stable_cell == UINT64_C(10)
            ? 0.30
            : 0.35;

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
                    composition[phase][component];
            }
        }
    }
    return q;
}

std::vector<double>
initial_natural_variables(
    std::uint64_t stable_cell) {
    auto values =
        target_natural_variables(
            stable_cell);
    values[0U] *= 1.02;
    values[1U] *= 1.015;
    values[2U] += 0.008;
    values[3U] -= 0.004;
    for (std::size_t column = 4U;
         column < values.size();
         ++column) {
        values[column] +=
            (column % 2U == 0U)
                ? 0.002
                : -0.002;
    }
    return values;
}

struct ZeroRelativePermeability3P {
    template <typename Number>
    flow::RelativePermeabilityEvaluation3P<Number>
    operator()(
        const flow::ThreePhaseSaturationState3P<Number>&)
        const {
        return {
            std::array<Number, 3>{
                Number{},
                Number{},
                Number{}}};
    }
};

struct ControlledClosureAudit {
    std::uint64_t calls{};
    std::array<std::vector<double>, 2>
        last_local_state;
};

std::size_t controlled_cell_slot(
    mesh::GlobalEntityId global) {
    if (global.value() == UINT64_C(10)) {
        return 0U;
    }
    if (global.value() == UINT64_C(20)) {
        return 1U;
    }
    throw std::invalid_argument(
        "unexpected controlled stable cell");
}

PetscErrorCode evaluate_controlled_cell(
    mesh::LocalIndex,
    mesh::GlobalEntityId cell_global,
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
        static_cast<ControlledClosureAudit*>(
            raw_context);
    ++audit->calls;

    try {
        if (natural_variables.size() !=
                frozen_layout.unknown_count() ||
            component_ids.size() !=
                frozen_layout.component_count()) {
            return PETSC_ERR_ARG_INCOMP;
        }

        const std::size_t q =
            frozen_layout.unknown_count();
        const double pressure =
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
            !std::isfinite(pressure) ||
            !(pressure > 0.0) ||
            !std::isfinite(temperature) ||
            !(temperature > 0.0)) {
            *status =
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        domain_error;
            return PETSC_SUCCESS;
        }

        flow::NaturalVariableCellStateInput3P
            state_input;
        state_input.component_ids.assign(
            component_ids.begin(),
            component_ids.end());
        state_input.reference_pressure_pa =
            pressure;
        state_input.temperature_k =
            temperature;
        state_input.independent_saturations = {
            natural_variables[*s0],
            natural_variables[*s1]};
        state_input.composition_pivot =
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
                const std::size_t component =
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
                state_input
                    .independent_phase_compositions[
                        phase]
                    .push_back(
                        natural_variables[*column]);
            }
        }

        const std::uint64_t stable =
            cell_global.value();
        const auto target =
            target_natural_variables(
                stable);
        const double target_pressure =
            target[0U];
        const double target_temperature =
            target[1U];
        constexpr std::array<double, 3>
            base_molar_density{
                2.0, 4.0, 7.0};
        constexpr std::array<double, 3>
            pressure_density_slope{
                0.20, 0.30, 0.40};
        constexpr std::array<double, 3>
            temperature_density_slope{
                0.05, 0.08, 0.11};
        constexpr std::array<double, 3>
            mass_density{
                1.0, 2.0, 3.0};
        constexpr std::array<double, 3>
            viscosity{
                1.0, 1.2, 1.4};

        std::array<double, 3>
            molar_density{};
        for (std::size_t phase = 0U;
             phase < 3U;
             ++phase) {
            molar_density[phase] =
                base_molar_density[phase] +
                pressure_density_slope[phase] *
                    (pressure -
                     target_pressure) +
                temperature_density_slope[phase] *
                    (temperature -
                     target_temperature);
            if (!std::isfinite(
                    molar_density[phase]) ||
                !(molar_density[phase] >
                  0.0)) {
                *status =
                    fdp::
                        NaturalVariableSnesEvaluationStatus3D::
                            domain_error;
                return PETSC_SUCCESS;
            }
            const double enthalpy =
                static_cast<double>(
                    phase + 1U) *
                temperature;
            const double internal_energy =
                0.5 *
                static_cast<double>(
                    phase + 1U) *
                temperature;
            state_input.phase_properties[phase] =
                flow::PhasePropertyPrerequisiteInput{
                    molar_density[phase],
                    mass_density[phase],
                    viscosity[phase],
                    enthalpy,
                    internal_energy};
        }

        auto state =
            flow::NaturalVariableCellState3P::
                create(
                    std::move(state_input));

        std::array<std::vector<double>, 3>
            molar_gradient;
        std::array<std::vector<double>, 3>
            zero_gradient;
        std::array<std::vector<double>, 3>
            enthalpy_gradient;
        std::array<std::vector<double>, 3>
            internal_energy_gradient;

        for (std::size_t phase = 0U;
             phase < 3U;
             ++phase) {
            molar_gradient[phase].assign(
                q,
                0.0);
            zero_gradient[phase].assign(
                q,
                0.0);
            enthalpy_gradient[phase].assign(
                q,
                0.0);
            internal_energy_gradient[phase]
                .assign(
                    q,
                    0.0);

            molar_gradient[phase][
                frozen_layout
                    .pressure_unknown_index()] =
                pressure_density_slope[phase];
            molar_gradient[phase][
                frozen_layout
                    .temperature_unknown_index()] =
                temperature_density_slope[phase];
            enthalpy_gradient[phase][
                frozen_layout
                    .temperature_unknown_index()] =
                static_cast<double>(
                    phase + 1U);
            internal_energy_gradient[phase][
                frozen_layout
                    .temperature_unknown_index()] =
                0.5 *
                static_cast<double>(
                    phase + 1U);
        }

        flow::PhaseMolarDensityNaturalVariableLinearization3P
            molar_linearization{
                frozen_layout,
                molar_density,
                molar_gradient};

        const flow::TransportPropertyProvenance
            controlled_provenance{
                "controlled-snes-assembly",
                "software-regression",
                "v1"};

        auto transport =
            flow::
                make_phase_transport_property_linearization(
                    state,
                    zero_gradient,
                    zero_gradient,
                    controlled_provenance,
                    controlled_provenance);
        auto caloric =
            flow::
                make_phase_caloric_property_linearization(
                    state,
                    enthalpy_gradient,
                    internal_energy_gradient,
                    controlled_provenance,
                    controlled_provenance);

        std::vector<double>
            rock_gradient(
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
                    2.0 *
                        temperature,
                    std::move(
                        rock_gradient),
                    controlled_provenance);

        auto saturation_primal =
            flow::
                evaluate_three_phase_saturation_constitutive(
                    pressure,
                    natural_variables[*s0],
                    natural_variables[*s1],
                    ZeroRelativePermeability3P{},
                    flow::NoCapillaryPressure3P{});
        flow::ThreePhaseSaturationCoordinateDerivatives3P
            saturation_derivatives{};
        auto saturation_linearization =
            flow::
                make_saturation_constitutive_natural_variable_linearization(
                    state,
                    saturation_primal,
                    saturation_derivatives);

        std::vector<double>
            fugacity_values(
                6U,
                0.0);
        std::vector<double>
            fugacity_jacobian(
                6U * q,
                0.0);
        for (std::size_t row = 0U;
             row < 6U;
             ++row) {
            const std::size_t column =
                4U + row;
            fugacity_values[row] =
                natural_variables[column] -
                target[column];
            fugacity_jacobian[
                row * q +
                column] =
                1.0;
        }
        flow::FugacityEquilibriumResidualLinearization3P
            fugacity{
                frozen_layout,
                std::vector<std::string>{
                    component_ids.begin(),
                    component_ids.end()},
                flow::FugacityEquilibriumResidual3P<double>{
                    3U,
                    std::move(
                        fugacity_values)},
                q,
                std::move(
                    fugacity_jacobian)};

        audit->last_local_state[
            controlled_cell_slot(
                cell_global)] =
            std::vector<double>{
                natural_variables.begin(),
                natural_variables.end()};

        output->emplace(
            fdp::FixedThreePhaseCurrentCellLinearization3D{
                std::move(state),
                std::move(
                    molar_linearization),
                std::move(transport),
                std::move(caloric),
                std::move(rock),
                std::move(
                    saturation_linearization),
                std::move(fugacity)});
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

mesh::Topology
make_topology(int rank) {
    mesh::Topology::EntityIds ids;
    ids.faces = {
        mesh::GlobalEntityId{
            UINT64_C(100)}};
    ids.cells =
        rank == 0
            ? std::vector<mesh::GlobalEntityId>{
                  mesh::GlobalEntityId{
                      UINT64_C(10)},
                  mesh::GlobalEntityId{
                      UINT64_C(20)}}
            : std::vector<mesh::GlobalEntityId>{
                  mesh::GlobalEntityId{
                      UINT64_C(20)},
                  mesh::GlobalEntityId{
                      UINT64_C(10)}};
    return {
        std::move(ids),
        {}};
}

mesh::PartitionSnapshot
make_partition(int rank) {
    const auto topology =
        make_topology(rank);
    mesh::EntityOwnerRanks owners;
    owners.faces = {
        mesh::PartitionRank{0U}};
    owners.cells =
        rank == 0
            ? std::vector<mesh::PartitionRank>{
                  mesh::PartitionRank{0U},
                  mesh::PartitionRank{1U}}
            : std::vector<mesh::PartitionRank>{
                  mesh::PartitionRank{1U},
                  mesh::PartitionRank{0U}};
    return mesh::PartitionSnapshot::create(
        topology,
        mesh::PartitionRank{
            static_cast<
                mesh::PartitionRank::value_type>(
                    rank)},
        2U,
        std::move(owners));
}

dp::ParallelOwnedConnectionSchedule3D
make_schedule(int rank) {
    const dp::AssemblyReadyInternalConnectionRow3D
        row{
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{
                UINT64_C(100)},
            rank == 0
                ? mesh::LocalIndex{0U}
                : mesh::LocalIndex{1U},
            mesh::GlobalEntityId{
                UINT64_C(10)},
            rank == 0
                ? mesh::LocalIndex{1U}
                : mesh::LocalIndex{0U},
            mesh::GlobalEntityId{
                UINT64_C(20)},
            2.0e-12};

    if (rank == 0) {
        return {
            mesh::PartitionRank{0U},
            2U,
            {row},
            {}};
    }
    return {
        mesh::PartitionRank{1U},
        2U,
        {},
        {row}};
}

mesh::DofLayout
make_dof_layout(int rank) {
    return mesh::DofLayout::create(
        make_topology(rank),
        {
            mesh::DofVariable{
                "natural_state",
                mesh::EntityKind::cell,
                10U}
        });
}

mesh::DofNumberingSnapshot
make_numbering(
    int rank,
    const mesh::DofLayout& layout,
    const mesh::PartitionSnapshot& partition) {
    mesh::GlobalEntityNumberingInput input;
    input.global_cell_count = 2U;
    input.global_face_count = 1U;
    input.faces = {
        {
            mesh::GlobalEntityId{
                UINT64_C(100)},
            mesh::GlobalEntityOrdinal{
                UINT64_C(0)}}
    };
    input.cells =
        rank == 0
            ? std::vector<
                  mesh::GlobalEntityOrdinalRecord>{
                  {
                      mesh::GlobalEntityId{
                          UINT64_C(10)},
                      mesh::GlobalEntityOrdinal{
                          UINT64_C(1)}},
                  {
                      mesh::GlobalEntityId{
                          UINT64_C(20)},
                      mesh::GlobalEntityOrdinal{
                          UINT64_C(0)}}
              }
            : std::vector<
                  mesh::GlobalEntityOrdinalRecord>{
                  {
                      mesh::GlobalEntityId{
                          UINT64_C(20)},
                      mesh::GlobalEntityOrdinal{
                          UINT64_C(0)}},
                  {
                      mesh::GlobalEntityId{
                          UINT64_C(10)},
                      mesh::GlobalEntityOrdinal{
                          UINT64_C(1)}}
              };
    return mesh::DofNumberingSnapshot::
        create_local(
            layout,
            partition,
            std::move(input));
}

dp::PetscMpiAijSymbolicPreallocation3D
make_cell_bridge(int rank) {
    const PetscInt row =
        static_cast<PetscInt>(
            rank);
    return {
        mesh::PartitionRank{
            static_cast<
                mesh::PartitionRank::value_type>(
                    rank)},
        2U,
        row,
        row + 1,
        2,
        {mesh::LocalIndex{0U}},
        {
            mesh::GlobalEntityId{
                rank == 0
                    ? UINT64_C(10)
                    : UINT64_C(20)}
        },
        {row},
        {1},
        {1},
        rank == 0
            ? std::vector<PetscInt>{
                  0, 1}
            : std::vector<PetscInt>{
                  1, 0}};
}

dp::OwnedCellStructuralColumnPatternSnapshot3D
make_cell_pattern(int rank) {
    const PetscInt row =
        static_cast<PetscInt>(
            rank);
    return {
        mesh::PartitionRank{
            static_cast<
                mesh::PartitionRank::value_type>(
                    rank)},
        2U,
        2U,
        row,
        row + 1,
        2,
        {
            dp::
                OwnedCellStructuralColumnPatternRow3D{
                    mesh::LocalIndex{0U},
                    mesh::GlobalEntityId{
                        rank == 0
                            ? UINT64_C(10)
                            : UINT64_C(20)},
                    row,
                    0U,
                    1U,
                    0U,
                    1U}
        },
        {row},
        {
            static_cast<PetscInt>(
                rank == 0
                    ? 1
                    : 0)}
    };
}

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

disc::TpfaInternalFaceTransmissibilityEntry3D
materialized_entry() {
    return {
        mesh::LocalIndex{0U},
        disc::
            TpfaInternalFaceTransmissibilityDisposition3D::
                materialized,
        direct_admissibility(),
        disc::TpfaStaticFaceTransmissibility3D{
            disc::
                TpfaStaticFaceTransmissibilityDisposition3D::
                    positive_harmonic_combination,
            1.5,
            2.0e-12}};
}

std::optional<
    fdp::FixedThreePhaseCurrentCellLinearization3D>
evaluate_direct(
    std::uint64_t stable_cell,
    const std::vector<double>& q,
    ControlledClosureAudit* audit) {
    std::optional<
        fdp::FixedThreePhaseCurrentCellLinearization3D>
        output;
    auto layout =
        target_layout(
            stable_cell);
    const std::vector<std::string>
        ids{"A", "B", "C"};
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    const PetscErrorCode error =
        evaluate_controlled_cell(
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{
                stable_cell},
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
                    success) {
        throw std::runtime_error(
            "controlled direct cell evaluation failed");
    }
    return output;
}

std::vector<
    fdp::FixedThreePhaseSnesCellInput3D>
make_cell_inputs(
    int rank,
    ControlledClosureAudit* audit) {
    std::vector<
        fdp::FixedThreePhaseSnesCellInput3D>
        result;
    result.reserve(2U);

    for (std::size_t local = 0U;
         local < 2U;
         ++local) {
        const std::uint64_t stable =
            rank == 0
                ? (local == 0U
                       ? UINT64_C(10)
                       : UINT64_C(20))
                : (local == 0U
                       ? UINT64_C(20)
                       : UINT64_C(10));
        const bool owned =
            (rank == 0 &&
             stable == UINT64_C(10)) ||
            (rank == 1 &&
             stable == UINT64_C(20));
        auto target =
            target_natural_variables(
                stable);
        auto evaluation =
            evaluate_direct(
                stable,
                target,
                audit);

        fdp::FixedThreePhaseSnesCellInput3D
            input;
        input.cell =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)};
        input.cell_global =
            mesh::GlobalEntityId{
                stable};
        input.bulk_volume_m3 =
            stable == UINT64_C(10)
                ? 2.0
                : 5.0;
        input.porosity =
            stable == UINT64_C(10)
                ? 0.25
                : 0.30;
        input.frozen_layout =
            target_layout(stable);
        input.component_ids =
            {"A", "B", "C"};

        if (owned) {
            input.previous_component_accumulation =
                flow::
                    build_pore_volume_component_accumulation(
                        evaluation->state,
                        input.porosity);
            input.previous_energy_accumulation =
                flow::
                    build_pore_volume_energy_accumulation_snapshot(
                        evaluation->state,
                        input.porosity,
                        evaluation->transport,
                        evaluation->caloric,
                        evaluation->rock);
        }
        result.push_back(
            std::move(input));
    }
    return result;
}

std::vector<
    fdp::FixedThreePhaseSnesAuthoritativeFaceInput3D>
make_face_inputs(
    int rank) {
    if (rank != 0) {
        return {};
    }
    return {
        {
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{
                UINT64_C(100)},
            materialized_entry(),
            flow::GravityVector3D{
                0.0, 0.0, -9.81},
            flow::OwnerToNeighbourDisplacement3D{
                0.0, 0.0, 1.0},
            mpmc::flow_discretization::
                StaticThermalFaceConductance3D{
                    0.0}}
    };
}

void set_owned_state(
    Vec state,
    int rank,
    bool target) {
    const std::uint64_t stable =
        rank == 0
            ? UINT64_C(10)
            : UINT64_C(20);
    auto values =
        target
            ? target_natural_variables(
                  stable)
            : initial_natural_variables(
                  stable);

    const PetscInt start =
        static_cast<PetscInt>(
            rank * 10);
    for (std::size_t slot = 0U;
         slot < values.size();
         ++slot) {
        const PetscInt index =
            start +
            static_cast<PetscInt>(
                slot);
        const PetscScalar value =
            static_cast<PetscScalar>(
                values[slot]);
        require_collective(
            VecSetValues(
                state,
                1,
                &index,
                &value,
                INSERT_VALUES) ==
                PETSC_SUCCESS,
            "failed to insert controlled owned state");
    }
    const PetscErrorCode begin =
        VecAssemblyBegin(state);
    const PetscErrorCode end =
        VecAssemblyEnd(state);
    require_collective(
        begin == PETSC_SUCCESS &&
            end == PETSC_SUCCESS,
        "failed to assemble controlled owned state");
}

} // namespace

void fixed_three_phase_snes_physical_assembly_test() {
    int rank = -1;
    int size = -1;
    if (MPI_Comm_rank(
            PETSC_COMM_WORLD,
            &rank) != MPI_SUCCESS ||
        MPI_Comm_size(
            PETSC_COMM_WORLD,
            &size) != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI rank/size query failed");
    }
    require_collective(
        size == 2,
        "fixed-three-phase SNES assembly test requires two ranks");

    auto partition =
        make_partition(rank);
    auto schedule =
        make_schedule(rank);
    auto dof_layout =
        make_dof_layout(rank);
    auto dof_numbering =
        make_numbering(
            rank,
            dof_layout,
            partition);
    auto cell_bridge =
        make_cell_bridge(rank);
    auto cell_pattern =
        make_cell_pattern(rank);

    ControlledClosureAudit audit;
    auto cells =
        make_cell_inputs(
            rank,
            &audit);
    auto faces =
        make_face_inputs(rank);

    std::optional<
        fdp::FixedThreePhaseSnesAssemblyContext3D>
        context;
    PetscErrorCode error =
        fdp::FixedThreePhaseSnesAssemblyContext3D::
            create(
                PETSC_COMM_WORLD,
                schedule,
                partition,
                dof_layout,
                dof_numbering,
                cell_bridge,
                cell_pattern,
                "natural_state",
                1.0,
                std::move(cells),
                std::move(faces),
                {
                    evaluate_controlled_cell,
                    &audit},
                &context);
    require_collective(
        error == PETSC_SUCCESS &&
            context.has_value(),
        "failed to create fixed-three-phase SNES assembly context");

    Vec initial = nullptr;
    require_collective(
        VecCreateMPI(
            PETSC_COMM_WORLD,
            10,
            20,
            &initial) ==
                PETSC_SUCCESS,
        "failed to create fixed-three-phase initial Vec");
    set_owned_state(
        initial,
        rank,
        false);

    std::optional<
        fdp::CompleteNaturalVariableAssemblySnapshot3D>
        initial_assembly;
    fdp::NaturalVariableSnesEvaluationStatus3D
        initial_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    error =
        context
            ->evaluate_complete_assembly(
                initial,
                &initial_assembly,
                &initial_status);
    require_collective(
        error == PETSC_SUCCESS &&
            initial_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            initial_assembly.has_value(),
        "initial fixed-three-phase production assembly failed");

    const std::uint64_t remote_stable =
        rank == 0
            ? UINT64_C(20)
            : UINT64_C(10);
    const auto remote_initial =
        initial_natural_variables(
            remote_stable);
    const auto& observed_remote =
        audit.last_local_state[
            controlled_cell_slot(
                mesh::GlobalEntityId{
                    remote_stable})];
    require_collective(
        observed_remote.size() ==
            remote_initial.size(),
        "PetscSF ghost state was not evaluated");
    for (std::size_t slot = 0U;
         slot < remote_initial.size();
         ++slot) {
        near_collective(
            observed_remote[slot],
            remote_initial[slot],
            0.0,
            0.0);
    }

    Vec unused_residual = nullptr;
    Mat jacobian_template = nullptr;
    error =
        fdp::
            materialize_complete_natural_variable_petsc_system_3d(
                PETSC_COMM_WORLD,
                *initial_assembly,
                cell_bridge,
                &unused_residual,
                &jacobian_template);
    require_collective(
        error == PETSC_SUCCESS &&
            unused_residual != nullptr &&
            jacobian_template != nullptr,
        "failed to materialize fixed-three-phase SNES Jacobian template");

    const std::uint64_t calls_before_snes =
        audit.calls;

    Vec solution = nullptr;
    std::optional<
        fdp::NaturalVariableSnesSolveReport3D>
        report;
    error =
        fdp::solve_natural_variable_snes_3d(
            PETSC_COMM_WORLD,
            *initial_assembly,
            initial,
            jacobian_template,
            context->snes_evaluator(),
            &solution,
            &report);
    require_collective(
        error == PETSC_SUCCESS &&
            solution != nullptr &&
            report.has_value(),
        "PETSc SNES failed on fixed-three-phase production assembly");

    require_collective(
        audit.calls >
                calls_before_snes &&
            report->function_evaluations() >
                0 &&
            report->jacobian_evaluations() >
                0 &&
            report->function_domain_errors() ==
                0 &&
            report->jacobian_domain_errors() ==
                0 &&
            static_cast<int>(
                report->converged_reason()) >
                0,
        "SNES did not repeatedly evaluate production flow assembly");

    const std::uint64_t stable =
        rank == 0
            ? UINT64_C(10)
            : UINT64_C(20);
    const auto target =
        target_natural_variables(
            stable);
    const PetscInt start =
        static_cast<PetscInt>(
            rank * 10);

    for (std::size_t slot = 0U;
         slot < target.size();
         ++slot) {
        const PetscInt index =
            start +
            static_cast<PetscInt>(
                slot);
        PetscScalar value{};
        require_collective(
            VecGetValues(
                solution,
                1,
                &index,
                &value) ==
                PETSC_SUCCESS,
            "failed to read converged fixed-three-phase state");
        near_collective(
            static_cast<double>(
                PetscRealPart(value)),
            target[slot]);
    }

    std::optional<
        fdp::CompleteNaturalVariableAssemblySnapshot3D>
        final_assembly;
    fdp::NaturalVariableSnesEvaluationStatus3D
        final_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    error =
        context
            ->evaluate_complete_assembly(
                solution,
                &final_assembly,
                &final_status);
    require_collective(
        error == PETSC_SUCCESS &&
            final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success &&
            final_assembly.has_value(),
        "independent final production assembly failed");

    double local_squared = 0.0;
    for (const auto& entry :
         final_assembly->residual_entries()) {
        local_squared +=
            entry.native_value *
            entry.native_value;
    }
    double global_squared = 0.0;
    if (MPI_Allreduce(
            &local_squared,
            &global_squared,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            PETSC_COMM_WORLD) !=
        MPI_SUCCESS) {
        throw std::runtime_error(
            "failed to reduce final assembly residual norm");
    }
    require_collective(
        std::sqrt(global_squared) <=
            2.0e-8 &&
            report->final_function_l2_norm() <=
                2.0e-8,
        "fixed-three-phase production assembly did not converge to target");

    for (const auto& entry :
         final_assembly->jacobian_entries()) {
        require_collective(
            std::isfinite(
                entry.value),
            "final production Jacobian contains non-finite value");
    }

    const PetscErrorCode solution_destroy =
        VecDestroy(&solution);
    const PetscErrorCode unused_destroy =
        VecDestroy(&unused_residual);
    const PetscErrorCode matrix_destroy =
        MatDestroy(&jacobian_template);
    const PetscErrorCode initial_destroy =
        VecDestroy(&initial);
    require_collective(
        solution_destroy ==
                PETSC_SUCCESS &&
            unused_destroy ==
                PETSC_SUCCESS &&
            matrix_destroy ==
                PETSC_SUCCESS &&
            initial_destroy ==
                PETSC_SUCCESS,
        "fixed-three-phase SNES fixture cleanup failed");
}
