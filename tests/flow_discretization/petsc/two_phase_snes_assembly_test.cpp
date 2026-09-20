#include <mpmc/flow_discretization_petsc/two_phase_snes_assembly.hpp>
#include <mpmc/flow_discretization_petsc/complete_natural_variable_petsc_materialization.hpp>

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
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

void require_collective_2p(
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
            PETSC_COMM_WORLD) != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI_Allreduce failed in two-phase SNES regression");
    }
    if (global == 0) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_collective_2p(
    double actual,
    double expected,
    double relative = 5.0e-8,
    double absolute = 5.0e-9) {
    const bool local =
        std::isfinite(actual) &&
        std::isfinite(expected) &&
        std::abs(actual - expected) <=
            absolute +
            relative *
                std::max(
                    std::abs(actual),
                    std::abs(expected));
    require_collective_2p(
        local,
        "two-phase distributed numeric mismatch");
}

flow::NaturalVariableLayout2P
target_layout_2p(std::uint64_t stable_cell) {
    if (stable_cell != UINT64_C(10) &&
        stable_cell != UINT64_C(20)) {
        throw std::invalid_argument(
            "unknown two-phase controlled cell");
    }
    return flow::NaturalVariableLayout2P{
        flow::NaturalVariableCompositionPivot2P::
            from_dependent_components(
                3U,
                {1U, 2U})};
}

std::vector<double>
target_state_2p(std::uint64_t stable_cell) {
    const auto layout =
        target_layout_2p(stable_cell);
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
    q[layout.independent_saturation_unknown_index()] =
        stable_cell == UINT64_C(10)
            ? 0.40
            : 0.45;

    q[*layout.independent_composition_unknown_index(
        0U,
        0U)] =
        stable_cell == UINT64_C(10)
            ? 0.20
            : 0.25;
    q[*layout.independent_composition_unknown_index(
        0U,
        2U)] =
        stable_cell == UINT64_C(10)
            ? 0.30
            : 0.25;

    q[*layout.independent_composition_unknown_index(
        1U,
        0U)] =
        stable_cell == UINT64_C(10)
            ? 0.35
            : 0.30;
    q[*layout.independent_composition_unknown_index(
        1U,
        1U)] =
        stable_cell == UINT64_C(10)
            ? 0.25
            : 0.30;
    return q;
}

std::vector<double>
initial_state_2p(std::uint64_t stable_cell) {
    auto q =
        target_state_2p(stable_cell);
    q[0] *= 1.03;
    q[1] *= 1.02;
    q[2] += 0.01;
    q[3] += 0.005;
    q[4] -= 0.004;
    q[5] -= 0.006;
    q[6] += 0.005;
    return q;
}

std::array<std::array<double, 3>, 2>
target_full_composition_2p(
    std::uint64_t stable_cell) {
    if (stable_cell == UINT64_C(10)) {
        return {{
            {{0.20, 0.50, 0.30}},
            {{0.35, 0.25, 0.40}}
        }};
    }
    if (stable_cell == UINT64_C(20)) {
        return {{
            {{0.25, 0.50, 0.25}},
            {{0.30, 0.30, 0.40}}
        }};
    }
    throw std::invalid_argument(
        "unknown two-phase target composition cell");
}

struct TwoPhaseControlledAudit {
    std::uint64_t calls{};
    std::vector<double> cell10;
    std::vector<double> cell20;
};

std::vector<double>& audit_cell_2p(
    TwoPhaseControlledAudit& audit,
    std::uint64_t stable) {
    if (stable == UINT64_C(10)) {
        return audit.cell10;
    }
    if (stable == UINT64_C(20)) {
        return audit.cell20;
    }
    throw std::invalid_argument(
        "unexpected two-phase stable cell");
}

PetscErrorCode evaluate_cell_2p(
    mesh::LocalIndex,
    mesh::GlobalEntityId cell_global,
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
        static_cast<TwoPhaseControlledAudit*>(
            raw_context);
    ++audit->calls;

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
            state_input;
        state_input.component_ids.assign(
            component_ids.begin(),
            component_ids.end());
        state_input.reference_pressure_pa = p;
        state_input.temperature_k =
            temperature;
        state_input.independent_saturation =
            s0;
        state_input.composition_pivot =
            frozen_layout.composition_pivot();

        for (std::size_t phase = 0U;
             phase < 2U;
             ++phase) {
            for (std::size_t rank = 0U;
                 rank <
                     frozen_layout.component_count() -
                         1U;
                 ++rank) {
                const std::size_t component =
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
                state_input
                    .independent_phase_compositions[
                        phase]
                    .push_back(
                        natural_variables[*column]);
            }
        }

        const auto target =
            target_state_2p(
                cell_global.value());
        const double target_p =
            target[0];
        const bool first =
            cell_global.value() == UINT64_C(10);
        const std::array<double, 2>
            base_density{
                first ? 6.0 : 7.0,
                first ? 3.0 : 4.0};
        const std::array<double, 2>
            density_slope{
                first ? 0.40 : 0.35,
                first ? 0.20 : 0.15};
        std::array<double, 2>
            molar_density{};
        for (std::size_t phase = 0U;
             phase < 2U;
             ++phase) {
            molar_density[phase] =
                base_density[phase] +
                density_slope[phase] *
                    (p - target_p);
            if (!std::isfinite(
                    molar_density[phase]) ||
                !(molar_density[phase] > 0.0)) {
                *status =
                    fdp::NaturalVariableSnesEvaluationStatus3D::
                        domain_error;
                return PETSC_SUCCESS;
            }
        }

        state_input.phase_properties[0] =
            flow::PhasePropertyPrerequisiteInput{
                molar_density[0],
                2.0,
                1.0,
                temperature,
                0.5 * temperature};
        state_input.phase_properties[1] =
            flow::PhasePropertyPrerequisiteInput{
                molar_density[1],
                1.5,
                1.0,
                1.2 * temperature,
                0.8 * temperature};

        auto state =
            flow::NaturalVariableCellState2P::
                create(std::move(state_input));

        const std::size_t q =
            frozen_layout.unknown_count();
        std::array<std::vector<double>, 2>
            dc{
                std::vector<double>(q, 0.0),
                std::vector<double>(q, 0.0)};
        dc[0][frozen_layout
                  .pressure_unknown_index()] =
            density_slope[0];
        dc[1][frozen_layout
                  .pressure_unknown_index()] =
            density_slope[1];
        auto molar =
            flow::make_two_phase_molar_density_linearization(
                state,
                dc);

        std::array<std::vector<double>, 2>
            zero{
                std::vector<double>(q, 0.0),
                std::vector<double>(q, 0.0)};
        const flow::TransportPropertyProvenance
            provenance{
                "controlled-two-phase-snes",
                "software-regression",
                "v1"};
        auto transport =
            flow::make_two_phase_transport_linearization(
                state,
                zero,
                zero,
                {0.0, 0.0},
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
            flow::make_two_phase_caloric_linearization(
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
            flow::make_two_phase_rock_thermal_storage_linearization(
                state,
                2.0 * temperature,
                rock_gradient,
                provenance);

        const auto target_x =
            target_full_composition_2p(
                cell_global.value());
        const std::array<double, 3>
            pressure_coupling{
                0.01, -0.015, 0.02};
        const std::array<double, 3>
            temperature_coupling{
                0.02, 0.01, -0.01};
        std::vector<double> fugacity_residual(
            3U,
            0.0);
        std::vector<double> fugacity_jacobian(
            3U * q,
            0.0);
        const auto x0 =
            state.phase_composition(0U);
        const auto x1 =
            state.phase_composition(1U);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            fugacity_residual[component] =
                std::log(
                    x0[component] /
                    x1[component]) -
                std::log(
                    target_x[0][component] /
                    target_x[1][component]) +
                pressure_coupling[component] *
                    (p - target[0]) +
                temperature_coupling[component] *
                    (temperature - target[1]);

            for (std::size_t column = 0U;
                 column < q;
                 ++column) {
                double derivative = 0.0;
                if (column ==
                    frozen_layout
                        .pressure_unknown_index()) {
                    derivative +=
                        pressure_coupling[
                            component];
                }
                if (column ==
                    frozen_layout
                        .temperature_unknown_index()) {
                    derivative +=
                        temperature_coupling[
                            component];
                }
                derivative +=
                    flow::two_phase_detail::
                        d_composition(
                            frozen_layout,
                            0U,
                            component,
                            column) /
                    x0[component];
                derivative -=
                    flow::two_phase_detail::
                        d_composition(
                            frozen_layout,
                            1U,
                            component,
                            column) /
                    x1[component];
                fugacity_jacobian[
                    component * q +
                    column] =
                    derivative;
            }
        }
        auto fugacity =
            flow::make_two_phase_fugacity_equilibrium_linearization(
                state,
                std::move(fugacity_residual),
                std::move(fugacity_jacobian));

        audit_cell_2p(
            *audit,
            cell_global.value()) =
            std::vector<double>{
                natural_variables.begin(),
                natural_variables.end()};

        output->emplace(
            fdp::TwoPhaseCurrentCellLinearization3D{
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

mesh::Topology make_topology_2p(int rank) {
    mesh::Topology::EntityIds ids;
    ids.faces = {
        mesh::GlobalEntityId{
            UINT64_C(100)}};
    ids.cells =
        rank == 0
            ? std::vector<mesh::GlobalEntityId>{
                  mesh::GlobalEntityId{UINT64_C(10)},
                  mesh::GlobalEntityId{UINT64_C(20)}}
            : std::vector<mesh::GlobalEntityId>{
                  mesh::GlobalEntityId{UINT64_C(20)},
                  mesh::GlobalEntityId{UINT64_C(10)}};
    return {
        std::move(ids),
        {}};
}

mesh::PartitionSnapshot
make_partition_2p(int rank) {
    auto topology =
        make_topology_2p(rank);
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
make_schedule_2p(int rank) {
    const dp::AssemblyReadyInternalConnectionRow3D
        row{
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{UINT64_C(100)},
            rank == 0
                ? mesh::LocalIndex{0U}
                : mesh::LocalIndex{1U},
            mesh::GlobalEntityId{UINT64_C(10)},
            rank == 0
                ? mesh::LocalIndex{1U}
                : mesh::LocalIndex{0U},
            mesh::GlobalEntityId{UINT64_C(20)},
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
make_dof_layout_2p(int rank) {
    return mesh::DofLayout::create(
        make_topology_2p(rank),
        {
            mesh::DofVariable{
                "natural_state_2p",
                mesh::EntityKind::cell,
                7U}
        });
}

mesh::DofNumberingSnapshot
make_numbering_2p(
    int rank,
    const mesh::DofLayout& layout,
    const mesh::PartitionSnapshot& partition) {
    mesh::GlobalEntityNumberingInput input;
    input.global_cell_count = 2U;
    input.global_face_count = 1U;
    input.faces = {
        {
            mesh::GlobalEntityId{UINT64_C(100)},
            mesh::GlobalEntityOrdinal{UINT64_C(0)}}
    };
    input.cells =
        rank == 0
            ? std::vector<
                  mesh::GlobalEntityOrdinalRecord>{
                  {
                      mesh::GlobalEntityId{UINT64_C(10)},
                      mesh::GlobalEntityOrdinal{UINT64_C(1)}},
                  {
                      mesh::GlobalEntityId{UINT64_C(20)},
                      mesh::GlobalEntityOrdinal{UINT64_C(0)}}}
            : std::vector<
                  mesh::GlobalEntityOrdinalRecord>{
                  {
                      mesh::GlobalEntityId{UINT64_C(20)},
                      mesh::GlobalEntityOrdinal{UINT64_C(0)}},
                  {
                      mesh::GlobalEntityId{UINT64_C(10)},
                      mesh::GlobalEntityOrdinal{UINT64_C(1)}}};
    return mesh::DofNumberingSnapshot::
        create_local(
            layout,
            partition,
            std::move(input));
}

dp::PetscMpiAijSymbolicPreallocation3D
make_cell_bridge_2p(int rank) {
    const PetscInt row =
        static_cast<PetscInt>(rank);
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
            ? std::vector<PetscInt>{0, 1}
            : std::vector<PetscInt>{1, 0}};
}

dp::OwnedCellStructuralColumnPatternSnapshot3D
make_cell_pattern_2p(int rank) {
    const PetscInt row =
        static_cast<PetscInt>(rank);
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
            dp::OwnedCellStructuralColumnPatternRow3D{
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
                rank == 0 ? 1 : 0)}
    };
}

disc::CombinedTransmissibilityAdmissibility3D
direct_admissibility_2p() {
    return {
        disc::
            CombinedTransmissibilityAdmissibilityDisposition3D::
                direct_normal_projection_k_orthogonal_candidate,
        {
            disc::TransmissibilityGeometryDisposition3D::
                direct_normal_projection_allowed,
            0.0,
            0.05},
        {
            disc::KOrthogonalityDisposition3D::
                k_orthogonal_within_policy,
            std::nullopt,
            std::nullopt,
            {
                std::nullopt,
                std::nullopt},
            0.05}};
}

disc::TpfaInternalFaceTransmissibilityEntry3D
materialized_entry_2p() {
    return {
        mesh::LocalIndex{0U},
        disc::TpfaInternalFaceTransmissibilityDisposition3D::
            materialized,
        direct_admissibility_2p(),
        disc::TpfaStaticFaceTransmissibility3D{
            disc::TpfaStaticFaceTransmissibilityDisposition3D::
                positive_harmonic_combination,
            1.5,
            2.0e-12}};
}

std::optional<
    fdp::TwoPhaseCurrentCellLinearization3D>
evaluate_direct_2p(
    std::uint64_t stable_cell,
    const std::vector<double>& q,
    TwoPhaseControlledAudit* audit) {
    std::optional<
        fdp::TwoPhaseCurrentCellLinearization3D>
        output;
    auto layout =
        target_layout_2p(stable_cell);
    const std::vector<std::string>
        ids{"A", "B", "C"};
    fdp::NaturalVariableSnesEvaluationStatus3D
        status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    const PetscErrorCode error =
        evaluate_cell_2p(
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{stable_cell},
            q,
            layout,
            ids,
            audit,
            &output,
            &status);
    if (error != PETSC_SUCCESS ||
        status !=
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success ||
        !output.has_value()) {
        throw std::runtime_error(
            "two-phase controlled direct evaluation failed");
    }
    return output;
}

std::vector<fdp::TwoPhaseSnesCellInput3D>
make_cells_2p(
    int rank,
    TwoPhaseControlledAudit* audit) {
    std::vector<fdp::TwoPhaseSnesCellInput3D>
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
            target_state_2p(stable);
        auto evaluation =
            evaluate_direct_2p(
                stable,
                target,
                audit);

        fdp::TwoPhaseSnesCellInput3D input;
        input.cell =
            mesh::LocalIndex{
                static_cast<
                    mesh::LocalIndex::value_type>(
                        local)};
        input.cell_global =
            mesh::GlobalEntityId{stable};
        input.bulk_volume_m3 =
            stable == UINT64_C(10)
                ? 2.0
                : 5.0;
        input.porosity =
            stable == UINT64_C(10)
                ? 0.25
                : 0.30;
        input.frozen_layout =
            target_layout_2p(stable);
        input.component_ids =
            {"A", "B", "C"};

        if (owned) {
            input.previous_component_accumulation =
                flow::build_two_phase_component_accumulation(
                    evaluation->state,
                    input.porosity);
            input.previous_energy_accumulation =
                flow::build_two_phase_energy_accumulation_snapshot(
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
    fdp::TwoPhaseSnesAuthoritativeFaceInput3D>
make_faces_2p(int rank) {
    if (rank != 0) {
        return {};
    }
    return {
        {
            mesh::LocalIndex{0U},
            mesh::GlobalEntityId{UINT64_C(100)},
            materialized_entry_2p(),
            flow::GravityVector3D{
                0.0, 0.0, -9.81},
            flow::OwnerToNeighbourDisplacement3D{
                0.0, 0.0, 1.0},
            mpmc::flow_discretization::
                StaticThermalFaceConductance3D{
                    0.0}}
    };
}

void set_owned_state_2p(
    Vec state,
    int rank,
    bool use_target) {
    const std::uint64_t stable =
        rank == 0
            ? UINT64_C(10)
            : UINT64_C(20);
    const auto values =
        use_target
            ? target_state_2p(stable)
            : initial_state_2p(stable);
    const PetscInt start =
        static_cast<PetscInt>(
            rank * 7);
    for (std::size_t slot = 0U;
         slot < values.size();
         ++slot) {
        const PetscInt index =
            start +
            static_cast<PetscInt>(slot);
        const PetscScalar value =
            static_cast<PetscScalar>(
                values[slot]);
        require_collective_2p(
            VecSetValues(
                state,
                1,
                &index,
                &value,
                INSERT_VALUES) ==
                PETSC_SUCCESS,
            "failed to insert two-phase owned state");
    }
    require_collective_2p(
        VecAssemblyBegin(state) ==
                PETSC_SUCCESS &&
            VecAssemblyEnd(state) ==
                PETSC_SUCCESS,
        "failed to assemble two-phase owned state");
}

} // namespace

void two_phase_snes_physical_assembly_test() {
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
    require_collective_2p(
        size == 2,
        "two-phase SNES test requires two ranks");

    auto partition =
        make_partition_2p(rank);
    auto schedule =
        make_schedule_2p(rank);
    auto dof_layout =
        make_dof_layout_2p(rank);
    auto dof_numbering =
        make_numbering_2p(
            rank,
            dof_layout,
            partition);
    auto cell_bridge =
        make_cell_bridge_2p(rank);
    auto cell_pattern =
        make_cell_pattern_2p(rank);

    TwoPhaseControlledAudit audit;
    auto cells =
        make_cells_2p(
            rank,
            &audit);
    auto faces =
        make_faces_2p(rank);

    std::optional<
        fdp::TwoPhaseSnesAssemblyContext3D>
        context;
    PetscErrorCode error =
        fdp::TwoPhaseSnesAssemblyContext3D::
            create(
                PETSC_COMM_WORLD,
                schedule,
                partition,
                dof_layout,
                dof_numbering,
                cell_bridge,
                cell_pattern,
                "natural_state_2p",
                1.0,
                std::move(cells),
                std::move(faces),
                {evaluate_cell_2p, &audit},
                &context);
    require_collective_2p(
        error == PETSC_SUCCESS &&
            context.has_value(),
        "failed to create two-phase SNES assembly context");

    Vec initial = nullptr;
    require_collective_2p(
        VecCreateMPI(
            PETSC_COMM_WORLD,
            7,
            14,
            &initial) ==
            PETSC_SUCCESS,
        "failed to create two-phase initial Vec");
    set_owned_state_2p(
        initial,
        rank,
        false);

    std::optional<
        fdp::CompleteNaturalVariableAssemblySnapshot3D>
        initial_assembly;
    fdp::NaturalVariableSnesEvaluationStatus3D
        initial_status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    error =
        context->evaluate_complete_assembly(
            initial,
            &initial_assembly,
            &initial_status);
    require_collective_2p(
        error == PETSC_SUCCESS &&
            initial_status ==
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    success &&
            initial_assembly.has_value(),
        "two-phase initial production assembly failed");
    require_collective_2p(
        initial_assembly->component_count() == 3U &&
            initial_assembly
                    ->natural_variable_count() ==
                7U &&
            initial_assembly
                    ->residual_entries()
                    .size() ==
                7U,
        "two-phase complete assembly did not reduce to 2*Nc+1");

    std::size_t fugacity_rows = 0U;
    for (const auto& entry :
         initial_assembly->residual_entries()) {
        if (entry.equation_kind ==
            fdp::NaturalVariableEquationKind3D::
                fugacity_equilibrium) {
            ++fugacity_rows;
            require_collective_2p(
                entry.equation_slot >= 4U &&
                    entry.equation_slot <= 6U,
                "two-phase fugacity row escaped active phase1 block");
        }
    }
    require_collective_2p(
        fugacity_rows == 3U,
        "two-phase assembly did not expose exactly Nc fugacity rows");

    const auto remote =
        initial_state_2p(
            rank == 0
                ? UINT64_C(20)
                : UINT64_C(10));
    const auto& observed_remote =
        rank == 0
            ? audit.cell20
            : audit.cell10;
    require_collective_2p(
        observed_remote.size() ==
            remote.size(),
        "two-phase PetscSF ghost state was not evaluated");
    for (std::size_t slot = 0U;
         slot < remote.size();
         ++slot) {
        near_collective_2p(
            observed_remote[slot],
            remote[slot],
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
    require_collective_2p(
        error == PETSC_SUCCESS &&
            unused_residual != nullptr &&
            jacobian_template != nullptr,
        "failed to materialize two-phase PETSc system");

    const std::uint64_t calls_before =
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
    require_collective_2p(
        error == PETSC_SUCCESS &&
            solution != nullptr &&
            report.has_value() &&
            report->natural_variable_count() ==
                7U &&
            static_cast<int>(
                report->converged_reason()) >
                0 &&
            report->function_domain_errors() ==
                0 &&
            report->jacobian_domain_errors() ==
                0 &&
            audit.calls > calls_before,
        "PETSc SNES failed on reduced two-phase production assembly");

    require_collective_2p(
        report->snes_type() ==
                std::string_view{SNESNEWTONLS} &&
            report->line_search_type() ==
                std::string_view{SNESLINESEARCHBT} &&
            report->ksp_type() ==
                std::string_view{KSPGMRES} &&
            report->pc_type() ==
                std::string_view{PCASM} &&
            report->default_asm_overlap() == 1,
        "two-phase solve did not use NewtonLS/BT + GMRES + ASM(1)");

    const auto target =
        target_state_2p(
            rank == 0
                ? UINT64_C(10)
                : UINT64_C(20));
    const PetscInt start =
        static_cast<PetscInt>(
            rank * 7);
    for (std::size_t slot = 0U;
         slot < target.size();
         ++slot) {
        const PetscInt index =
            start +
            static_cast<PetscInt>(slot);
        PetscScalar value{};
        require_collective_2p(
            VecGetValues(
                solution,
                1,
                &index,
                &value) ==
                PETSC_SUCCESS,
            "failed to read converged two-phase state");
        near_collective_2p(
            static_cast<double>(
                PetscRealPart(value)),
            target[slot]);
    }

    std::optional<
        fdp::CompleteNaturalVariableAssemblySnapshot3D>
        final_assembly;
    fdp::NaturalVariableSnesEvaluationStatus3D
        final_status =
            fdp::NaturalVariableSnesEvaluationStatus3D::
                success;
    error =
        context->evaluate_complete_assembly(
            solution,
            &final_assembly,
            &final_status);
    require_collective_2p(
        error == PETSC_SUCCESS &&
            final_status ==
                fdp::NaturalVariableSnesEvaluationStatus3D::
                    success &&
            final_assembly.has_value() &&
            final_assembly
                    ->natural_variable_count() ==
                7U,
        "independent final two-phase assembly failed");

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
            "failed to reduce two-phase final residual");
    }
    require_collective_2p(
        std::sqrt(global_squared) <=
            5.0e-8 &&
            report->final_function_l2_norm() <=
                5.0e-8,
        "two-phase production assembly did not converge to target");

    for (const auto& entry :
         final_assembly->jacobian_entries()) {
        require_collective_2p(
            std::isfinite(entry.value) &&
                entry.natural_variable_column <
                    7U,
            "two-phase final Jacobian is invalid");
    }

    const PetscErrorCode solution_destroy =
        VecDestroy(&solution);
    const PetscErrorCode residual_destroy =
        VecDestroy(&unused_residual);
    const PetscErrorCode matrix_destroy =
        MatDestroy(&jacobian_template);
    const PetscErrorCode initial_destroy =
        VecDestroy(&initial);
    require_collective_2p(
        solution_destroy == PETSC_SUCCESS &&
            residual_destroy == PETSC_SUCCESS &&
            matrix_destroy == PETSC_SUCCESS &&
            initial_destroy == PETSC_SUCCESS,
        "two-phase SNES fixture cleanup failed");
}
