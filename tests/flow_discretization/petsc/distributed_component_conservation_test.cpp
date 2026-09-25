#include <mpmc/flow_discretization_petsc/distributed_component_conservation.hpp>
#include <mpmc/flow_discretization_petsc/distributed_energy_conservation.hpp>
#include <mpmc/flow_discretization_petsc/energy_global_assembly_mapping.hpp>
#include <mpmc/flow_discretization_petsc/complete_natural_variable_assembly_snapshot.hpp>
#include <mpmc/flow_discretization_petsc/complete_natural_variable_petsc_materialization.hpp>
#include <mpmc/flow_discretization_petsc/natural_variable_newton_linear_system.hpp>
#include <mpmc/flow_discretization_petsc/natural_variable_snes_solver.hpp>
#include <mpmc/flow_discretization_petsc/global_component_assembly_mapping.hpp>
#include <mpmc/flow_discretization_petsc/fugacity_equilibrium_global_assembly_mapping.hpp>

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

bool distributed_component_conservation_header();
bool distributed_energy_conservation_header();
bool energy_global_assembly_mapping_header();
bool complete_natural_variable_assembly_snapshot_header();
bool complete_natural_variable_petsc_materialization_header();
bool natural_variable_newton_linear_system_header();
bool natural_variable_snes_solver_header();
bool fixed_three_phase_snes_assembly_header();
void fixed_three_phase_snes_physical_assembly_test();
bool single_phase_snes_assembly_header();
void single_phase_snes_physical_assembly_test();
bool two_phase_snes_assembly_header();
void two_phase_snes_physical_assembly_test();
bool variable_cardinality_natural_variable_numbering_header();
void variable_cardinality_natural_variable_numbering_test();
bool variable_cardinality_snes_assembly_header();
void variable_cardinality_snes_assembly_test();
bool mixed_cardinality_physical_snes_assembly_header();
void mixed_cardinality_physical_snes_assembly_test();
bool fixed_bhp_well_source_evaluator_header();
bool pr76_production_cell_evaluator_header();
bool pr76_single_phase_transition_target_rebuild_header();
bool pr76_transition_rebuild_materialization_header();
void pr76_production_fully_implicit_transient_test();
void pr76_li_firoozabadi_sour_gas_short_step_test();
bool adaptive_timestep_controller_header();
bool physical_timestep_driver_header();
bool physical_time_loop_header();
void adaptive_timestep_controller_test();
bool single_phase_adaptive_timestep_attempt_header();
bool single_phase_handoff_initial_system_header();
bool cross_cardinality_tpfa_bridge_header();
bool thermodynamic_cross_cardinality_tpfa_adapter_header();
bool phase_transition_outer_rebuild_header();
bool post_snes_phase_transition_controller_header();
bool post_snes_phase_transition_handoff_scanner_header();
bool post_snes_pt_flash_phase_transition_scanner_header();
void post_snes_pt_flash_phase_transition_scanner_test();
bool global_component_assembly_mapping_header();
bool fugacity_equilibrium_global_assembly_mapping_header();

namespace {

namespace dp =
    mpmc::discretization_petsc;
namespace fd =
    mpmc::flow_discretization;
namespace fdp =
    mpmc::flow_discretization_petsc;
namespace flow = mpmc::flow;
namespace mesh = mpmc::mesh;

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
            PETSC_COMM_WORLD) != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI_Allreduce failed in test assertion");
    }
    if (global == 0) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_collective(
    double actual,
    double expected,
    double relative = 5.0e-12,
    double absolute = 5.0e-12) {
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
        "distributed numeric mismatch");
}

flow::NaturalVariableStateIdentity3P
identity(
    std::uint64_t stable_cell) {
    if (stable_cell == 10U) {
        return {
            flow::NaturalVariableLayout3P{
                flow::NaturalVariableCompositionPivot3P::
                    from_dependent_components(
                        3U,
                        {1U, 0U, 2U})},
            {"A", "B", "C"},
            1.0e6,
            350.0,
            {0.20, 0.30, 0.50},
            {
                std::vector<double>{0.10, 0.70, 0.20},
                std::vector<double>{0.60, 0.20, 0.20},
                std::vector<double>{0.20, 0.30, 0.50}
            }};
    }
    if (stable_cell == 20U) {
        return {
            flow::NaturalVariableLayout3P{
                flow::NaturalVariableCompositionPivot3P::
                    from_dependent_components(
                        3U,
                        {0U, 2U, 1U})},
            {"A", "B", "C"},
            0.99e6,
            351.0,
            {0.25, 0.35, 0.40},
            {
                std::vector<double>{0.55, 0.25, 0.20},
                std::vector<double>{0.25, 0.25, 0.50},
                std::vector<double>{0.20, 0.60, 0.20}
            }};
    }
    throw std::invalid_argument(
        "unknown stable cell fixture");
}

std::vector<double> component_jacobian(
    std::size_t component_count,
    std::size_t column_count,
    double scale) {
    std::vector<double> result(
        component_count * column_count,
        0.0);
    for (std::size_t component = 0U;
         component < component_count;
         ++component) {
        for (std::size_t column = 0U;
             column < column_count;
             ++column) {
            result[
                component * column_count +
                column] =
                scale *
                static_cast<double>(
                    1U +
                    component * 10U +
                    column);
        }
    }
    return result;
}

std::vector<double> total_gradient(
    const std::vector<double>& values,
    std::size_t component_count,
    std::size_t column_count) {
    std::vector<double> result(
        column_count,
        0.0);
    for (std::size_t column = 0U;
         column < column_count;
         ++column) {
        for (std::size_t component = 0U;
             component < component_count;
             ++component) {
            result[column] +=
                values[
                    component * column_count +
                    column];
        }
    }
    return result;
}

flow::BackwardEulerComponentAccumulationResidual3P
accumulation(
    std::uint64_t stable_cell) {
    auto state =
        identity(stable_cell);
    const std::size_t q =
        state.layout.unknown_count();
    const double scale =
        stable_cell == 10U
            ? 1.0
            : 1.5;

    std::vector<double> residual{
        0.10 * scale,
        -0.04 * scale,
        0.06 * scale};
    double total = 0.0;
    for (double value : residual) {
        total += value;
    }
    const auto jacobian =
        component_jacobian(
            3U,
            q,
            0.002 * scale);

    return {
        state.layout,
        stable_cell == 10U
            ? 0.20
            : 0.25,
        120.0,
        {"A", "B", "C"},
        std::move(residual),
        total,
        q,
        jacobian,
        total_gradient(
            jacobian,
            3U,
            q)};
}

fd::ConservativeComponentFaceRateScatterLinearization3D
face_scatter() {
    auto owner =
        identity(10U);
    auto neighbour =
        identity(20U);
    const std::size_t owner_q =
        owner.layout.unknown_count();
    const std::size_t neighbour_q =
        neighbour.layout.unknown_count();

    const std::vector<double>
        owner_rate{2.0, -1.0, 3.0};
    std::vector<double>
        neighbour_rate{-2.0, 1.0, -3.0};

    const auto owner_owner =
        component_jacobian(
            3U,
            owner_q,
            0.01);
    const auto owner_neighbour =
        component_jacobian(
            3U,
            neighbour_q,
            -0.02);

    std::vector<double> neighbour_owner(
        owner_owner.size());
    std::vector<double> neighbour_neighbour(
        owner_neighbour.size());
    for (std::size_t index = 0U;
         index < owner_owner.size();
         ++index) {
        neighbour_owner[index] =
            -owner_owner[index];
    }
    for (std::size_t index = 0U;
         index < owner_neighbour.size();
         ++index) {
        neighbour_neighbour[index] =
            -owner_neighbour[index];
    }

    const auto owner_total_owner =
        total_gradient(
            owner_owner,
            3U,
            owner_q);
    const auto owner_total_neighbour =
        total_gradient(
            owner_neighbour,
            3U,
            neighbour_q);

    std::vector<double> neighbour_total_owner(
        owner_q);
    std::vector<double> neighbour_total_neighbour(
        neighbour_q);
    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        neighbour_total_owner[column] =
            -owner_total_owner[column];
    }
    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        neighbour_total_neighbour[column] =
            -owner_total_neighbour[column];
    }

    return {
        mesh::LocalIndex{0U},
        {"A", "B", "C"},
        std::move(owner),
        std::move(neighbour),
        owner_rate,
        std::move(neighbour_rate),
        owner_owner,
        owner_neighbour,
        std::move(neighbour_owner),
        std::move(neighbour_neighbour),
        4.0,
        -4.0,
        owner_total_owner,
        owner_total_neighbour,
        std::move(neighbour_total_owner),
        std::move(neighbour_total_neighbour)};
}

fd::NormalizedComponentFaceContributionLinearization3D
normalized_face() {
    return fd::
        normalize_component_face_rate_by_bulk_volume(
            face_scatter(),
            {2.0, 5.0});
}

mesh::Topology
make_topology(int rank) {
    mesh::Topology::EntityIds ids;
    ids.faces = {
        mesh::GlobalEntityId{100U}};
    ids.cells =
        rank == 0
            ? std::vector<mesh::GlobalEntityId>{
                  mesh::GlobalEntityId{10U},
                  mesh::GlobalEntityId{20U}}
            : std::vector<mesh::GlobalEntityId>{
                  mesh::GlobalEntityId{20U},
                  mesh::GlobalEntityId{10U}};
    return mesh::Topology{
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
        row =
            rank == 0
                ? dp::AssemblyReadyInternalConnectionRow3D{
                      mesh::LocalIndex{0U},
                      mesh::GlobalEntityId{100U},
                      mesh::LocalIndex{0U},
                      mesh::GlobalEntityId{10U},
                      mesh::LocalIndex{1U},
                      mesh::GlobalEntityId{20U},
                      1.0}
                : dp::AssemblyReadyInternalConnectionRow3D{
                      mesh::LocalIndex{0U},
                      mesh::GlobalEntityId{100U},
                      mesh::LocalIndex{1U},
                      mesh::GlobalEntityId{10U},
                      mesh::LocalIndex{0U},
                      mesh::GlobalEntityId{20U},
                      1.0};

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

struct LocalInputs {
    flow::BackwardEulerComponentAccumulationResidual3P
        accumulation10;
    flow::BackwardEulerComponentAccumulationResidual3P
        accumulation20;
    std::vector<fdp::DistributedCellStateBinding3D>
        cells;
    fd::NormalizedComponentFaceContributionLinearization3D
        face;
    std::vector<
        fd::AuthoritativeNormalizedFaceContributionBinding3D>
        faces;

    explicit LocalInputs(int rank)
        : accumulation10(
              accumulation(10U)),
          accumulation20(
              accumulation(20U)),
          face(normalized_face()) {
        if (rank == 0) {
            cells = {
                {
                    mesh::LocalIndex{0U},
                    mesh::GlobalEntityId{10U},
                    2.0,
                    identity(10U),
                    &accumulation10},
                {
                    mesh::LocalIndex{1U},
                    mesh::GlobalEntityId{20U},
                    5.0,
                    identity(20U),
                    nullptr}
            };
            faces = {
                {
                    mesh::LocalIndex{0U},
                    mesh::GlobalEntityId{100U},
                    &face}
            };
        } else {
            cells = {
                {
                    mesh::LocalIndex{0U},
                    mesh::GlobalEntityId{20U},
                    5.0,
                    identity(20U),
                    &accumulation20},
                {
                    mesh::LocalIndex{1U},
                    mesh::GlobalEntityId{10U},
                    2.0,
                    identity(10U),
                    nullptr}
            };
        }
    }
};

std::optional<
    fdp::DistributedOwnedMultiCellComponentConservationSnapshot3D>
build(
    int rank,
    mesh::PartitionSnapshot& partition,
    dp::ParallelOwnedConnectionSchedule3D& schedule,
    LocalInputs& inputs,
    PetscErrorCode* error_out) {
    std::optional<
        fdp::DistributedOwnedMultiCellComponentConservationSnapshot3D>
        output;
    const PetscErrorCode error =
        fdp::
            make_distributed_owned_component_conservation_snapshot_3d(
                PETSC_COMM_WORLD,
                schedule,
                partition,
                inputs.cells,
                inputs.faces,
                &output);
    if (error_out != nullptr) {
        *error_out = error;
    }
    (void)rank;
    return output;
}

const fd::OwnedCellPairComponentJacobianBlock3D*
find_block(
    const fd::OwnedCellComponentConservationRow3D& row,
    std::uint64_t column_global) {
    const auto found =
        std::find_if(
            row.off_diagonal_cell_pair_blocks.begin(),
            row.off_diagonal_cell_pair_blocks.end(),
            [column_global](const auto& block) {
                return block.column_cell_global ==
                    mesh::GlobalEntityId{
                        column_global};
            });
    return found ==
               row.off_diagonal_cell_pair_blocks.end()
        ? nullptr
        : &*found;
}

void audit_cross_rank_jacobian(
    int rank,
    const fdp::
        DistributedOwnedMultiCellComponentConservationSnapshot3D&
            snapshot,
    const LocalInputs& inputs) {
    const auto& row =
        snapshot.owned_rows().front();
    const auto& owned_binding =
        inputs.cells.front();
    const auto* owned_accumulation =
        owned_binding.owned_accumulation;
    require_collective(
        owned_accumulation != nullptr,
        "owned accumulation missing in Jacobian audit");

    for (const std::uint64_t source_global :
         {UINT64_C(10), UINT64_C(20)}) {
        const std::size_t q =
            identity(source_global)
                .layout.unknown_count();

        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            for (std::size_t column = 0U;
                 column < q;
                 ++column) {
                const bool source_row =
                    row.cell_global ==
                    mesh::GlobalEntityId{
                        source_global};
                const auto* block =
                    source_row
                        ? nullptr
                        : find_block(
                              row,
                              source_global);
                require_collective(
                    source_row ||
                        block != nullptr,
                    "missing cross-rank off-diagonal block");

                double local = 0.0;
                if (source_row) {
                    local =
                        owned_binding.bulk_volume_m3 *
                        (row.local_residual
                             .d_local(
                                 component,
                                 column) -
                         owned_accumulation
                             ->d_residual(
                                 component,
                                 column));
                } else {
                    local =
                        owned_binding.bulk_volume_m3 *
                        block->d_component(
                            component,
                            column);
                }

                double global = 0.0;
                if (MPI_Allreduce(
                        &local,
                        &global,
                        1,
                        MPI_DOUBLE,
                        MPI_SUM,
                        PETSC_COMM_WORLD) !=
                    MPI_SUCCESS) {
                    throw std::runtime_error(
                        "MPI_Allreduce failed in Jacobian audit");
                }
                near_collective(
                    global,
                    0.0,
                    0.0,
                    2.0e-13);
            }
        }

        for (std::size_t column = 0U;
             column < q;
             ++column) {
            const bool source_row =
                row.cell_global ==
                mesh::GlobalEntityId{
                    source_global};
            const auto* block =
                source_row
                    ? nullptr
                    : find_block(
                          row,
                          source_global);
            require_collective(
                source_row ||
                    block != nullptr,
                "missing total cross-rank block");

            double local = 0.0;
            if (source_row) {
                local =
                    owned_binding.bulk_volume_m3 *
                    (row.local_residual
                         .d_total_local(column) -
                     owned_accumulation
                         ->d_total(column));
            } else {
                local =
                    owned_binding.bulk_volume_m3 *
                    block->d_total(column);
            }

            double global = 0.0;
            if (MPI_Allreduce(
                    &local,
                    &global,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    PETSC_COMM_WORLD) !=
                MPI_SUCCESS) {
                throw std::runtime_error(
                    "MPI_Allreduce failed in total Jacobian audit");
            }
            near_collective(
                global,
                0.0,
                0.0,
                2.0e-13);
        }
    }

    (void)rank;
}

void distributed_exchange() {
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
        "distributed conservation test requires exactly two ranks");

    auto partition =
        make_partition(rank);
    auto schedule =
        make_schedule(rank);
    LocalInputs inputs{rank};

    PetscErrorCode error =
        PETSC_SUCCESS;
    auto output =
        build(
            rank,
            partition,
            schedule,
            inputs,
            &error);

    require_collective(
        error == PETSC_SUCCESS &&
            output.has_value(),
        "distributed owner-targeted exchange failed");

    const auto& snapshot =
        *output;
    require_collective(
        snapshot.owned_rows().size() == 1U &&
            snapshot.received_face_side_count() == 1U &&
            snapshot.local_authoritative_face_count() ==
                (rank == 0 ? 1U : 0U),
        "distributed snapshot cardinality mismatch");

    const auto& row =
        snapshot.owned_rows().front();
    const std::uint64_t expected_cell =
        rank == 0
            ? UINT64_C(10)
            : UINT64_C(20);
    const std::uint64_t expected_column =
        rank == 0
            ? UINT64_C(20)
            : UINT64_C(10);

    require_collective(
        row.cell ==
                mesh::LocalIndex{0U} &&
            row.cell_global ==
                mesh::GlobalEntityId{
                    expected_cell} &&
            row.off_diagonal_cell_pair_blocks.size() ==
                1U,
        "owned row identity/local reordering mismatch");

    const auto& block =
        row.off_diagonal_cell_pair_blocks.front();
    require_collective(
        block.column_cell ==
                mesh::LocalIndex{1U} &&
            block.column_cell_global ==
                mesh::GlobalEntityId{
                    expected_column} &&
            block.contributing_face_global_ids ==
                std::vector<mesh::GlobalEntityId>{
                    mesh::GlobalEntityId{100U}},
        "stable-ID off-diagonal block mapping mismatch");

    const auto face =
        normalized_face();
    const auto& accumulation_ref =
        rank == 0
            ? inputs.accumulation10
            : inputs.accumulation20;

    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        const double spatial =
            rank == 0
                ? face.owner_contribution(component)
                : face.neighbour_contribution(component);
        near_collective(
            row.local_residual.residual(component),
            accumulation_ref.residual(component) +
                spatial);
    }

    for (double value :
         snapshot
             .global_volume_weighted_spatial_component_balance_mol_per_s()) {
        near_collective(
            value,
            0.0,
            0.0,
            2.0e-13);
    }
    near_collective(
        snapshot
            .global_volume_weighted_spatial_total_balance_mol_per_s(),
        0.0,
        0.0,
        2.0e-13);

    audit_cross_rank_jacobian(
        rank,
        snapshot,
        inputs);
}

mesh::DofLayout
make_natural_variable_dof_layout(
    int rank,
    std::size_t q) {
    const auto topology =
        make_topology(rank);
    return mesh::DofLayout::create(
        topology,
        {
            mesh::DofVariable{
                "natural_state",
                mesh::EntityKind::cell,
                q}
        });
}

mesh::DofNumberingSnapshot
make_reversed_mesh_dof_numbering(
    int rank,
    const mesh::DofLayout& layout,
    const mesh::PartitionSnapshot& partition) {
    mesh::GlobalEntityNumberingInput numbering;
    numbering.global_cell_count = 2U;
    numbering.global_face_count = 1U;
    numbering.faces = {
        {
            mesh::GlobalEntityId{100U},
            mesh::GlobalEntityOrdinal{0U}}
    };
    numbering.cells =
        rank == 0
            ? std::vector<
                  mesh::GlobalEntityOrdinalRecord>{
                  {
                      mesh::GlobalEntityId{10U},
                      mesh::GlobalEntityOrdinal{1U}},
                  {
                      mesh::GlobalEntityId{20U},
                      mesh::GlobalEntityOrdinal{0U}}
              }
            : std::vector<
                  mesh::GlobalEntityOrdinalRecord>{
                  {
                      mesh::GlobalEntityId{20U},
                      mesh::GlobalEntityOrdinal{0U}},
                  {
                      mesh::GlobalEntityId{10U},
                      mesh::GlobalEntityOrdinal{1U}}
              };

    return mesh::DofNumberingSnapshot::
        create_local(
            layout,
            partition,
            std::move(numbering));
}

dp::PetscMpiAijSymbolicPreallocation3D
make_cell_row_bridge(int rank) {
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
make_cell_column_pattern(int rank) {
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

dp::OwnedCellStructuralColumnPatternSnapshot3D
make_incomplete_cell_column_pattern(
    int rank) {
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
                0U}
        },
        {row},
        {}
    };
}

void global_component_assembly_mapping() {
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
        "global mapping test requires exactly two ranks");

    auto partition =
        make_partition(rank);
    auto schedule =
        make_schedule(rank);
    LocalInputs inputs{rank};

    // Force one exact-zero diagonal Jacobian scalar on rank 0. The mapping
    // bridge must retain it because numerical zero is not structural zero.
    if (rank == 0) {
        const double face_derivative =
            inputs.face
                .d_owner_contribution_wrt_owner(
                    0U,
                    0U);
        const double old =
            inputs.accumulation10
                .component_jacobian[0U];
        inputs.accumulation10
            .component_jacobian[0U] =
            -face_derivative;
        inputs.accumulation10
            .total_residual_gradient[0U] +=
            -face_derivative -
            old;
    }

    PetscErrorCode error =
        PETSC_SUCCESS;
    auto conservation =
        build(
            rank,
            partition,
            schedule,
            inputs,
            &error);
    require_collective(
        error == PETSC_SUCCESS &&
            conservation.has_value(),
        "distributed snapshot for global mapping failed");

    const std::size_t q =
        identity(10U)
            .layout.unknown_count();
    auto layout =
        make_natural_variable_dof_layout(
            rank,
            q);
    auto numbering =
        make_reversed_mesh_dof_numbering(
            rank,
            layout,
            partition);
    auto cell_bridge =
        make_cell_row_bridge(
            rank);
    auto cell_pattern =
        make_cell_column_pattern(
            rank);

    std::optional<
        fdp::ComponentConservationGlobalAssemblyEntries3D>
        mapping;
    error =
        fdp::
            make_component_conservation_global_assembly_entries_3d(
                PETSC_COMM_WORLD,
                *conservation,
                partition,
                layout,
                numbering,
                cell_bridge,
                cell_pattern,
                "natural_state",
                &mapping);

    require_collective(
        error == PETSC_SUCCESS &&
            mapping.has_value(),
        "global component assembly mapping failed");

    const auto& result =
        *mapping;
    const PetscInt expected_start =
        static_cast<PetscInt>(
            rank) *
        static_cast<PetscInt>(q);
    const PetscInt expected_end =
        expected_start +
        static_cast<PetscInt>(q);
    const PetscInt expected_count =
        2 *
        static_cast<PetscInt>(q);

    require_collective(
        result.petsc_scalar_row_start() ==
                expected_start &&
            result.petsc_scalar_row_end() ==
                expected_end &&
            result.petsc_scalar_row_count() ==
                expected_count &&
            result.component_count() == 3U &&
            result.natural_variable_count() ==
                q &&
            result.residual_entries().size() ==
                3U &&
            result.jacobian_entries().size() ==
                3U * q * 2U,
        "global mapping scalar range/cardinality mismatch");

    const auto& owned_row =
        conservation->owned_rows().front();
    const auto* off_block =
        find_block(
            owned_row,
            rank == 0
                ? UINT64_C(20)
                : UINT64_C(10));
    require_collective(
        off_block != nullptr,
        "global mapping fixture missing off-diagonal block");

    const std::uint64_t row_cell_global =
        rank == 0
            ? UINT64_C(10)
            : UINT64_C(20);
    const std::uint64_t other_cell_global =
        rank == 0
            ? UINT64_C(20)
            : UINT64_C(10);

    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        const auto& residual =
            result.residual_entries()[
                component];
        const PetscInt expected_row =
            expected_start +
            static_cast<PetscInt>(
                component);
        const std::uint64_t
            expected_mesh_row =
                (rank == 0
                     ? static_cast<std::uint64_t>(q)
                     : UINT64_C(0)) +
                static_cast<std::uint64_t>(
                    component);

        require_collective(
            residual.petsc_global_row ==
                    expected_row &&
                residual.mesh_global_row_dof ==
                    mesh::GlobalDofIndex{
                        expected_mesh_row} &&
                residual.row_cell_global ==
                    mesh::GlobalEntityId{
                        row_cell_global} &&
                residual.component_row ==
                    component,
            "residual scalar mapping/provenance mismatch");
        near_collective(
            residual.value_mol_per_bulk_m3_s,
            owned_row.local_residual
                .residual(component));

        // Deliberately reversed mesh-global entity ordinals prove that the
        // mesh-global DoF provenance is not being used as PETSc row ownership.
        require_collective(
            static_cast<std::uint64_t>(
                residual.petsc_global_row) !=
                residual
                    .mesh_global_row_dof
                    .value(),
            "mesh-global DoF was incorrectly conflated with PETSc scalar row");

        for (PetscInt global_column = 0;
             global_column <
             expected_count;
             ++global_column) {
            const std::size_t entry_index =
                component *
                    static_cast<std::size_t>(
                        expected_count) +
                static_cast<std::size_t>(
                    global_column);
            const auto& entry =
                result.jacobian_entries()[
                    entry_index];

            const bool diagonal =
                (rank == 0 &&
                 global_column <
                     static_cast<PetscInt>(
                         q)) ||
                (rank == 1 &&
                 global_column >=
                     static_cast<PetscInt>(
                         q));
            const std::size_t column =
                static_cast<std::size_t>(
                    global_column %
                    static_cast<PetscInt>(
                        q));
            const std::uint64_t
                column_cell_global =
                    diagonal
                        ? row_cell_global
                        : other_cell_global;
            const double expected_value =
                diagonal
                    ? owned_row.local_residual
                          .d_local(
                              component,
                              column)
                    : off_block->d_component(
                          component,
                          column);

            const std::uint64_t
                expected_mesh_column =
                    column_cell_global ==
                            UINT64_C(10)
                        ? static_cast<std::uint64_t>(
                              q) +
                              static_cast<std::uint64_t>(
                                  column)
                        : static_cast<std::uint64_t>(
                              column);

            require_collective(
                entry.petsc_global_row ==
                        expected_row &&
                    entry.petsc_global_column ==
                        global_column &&
                    entry.mesh_global_row_dof ==
                        residual
                            .mesh_global_row_dof &&
                    entry.mesh_global_column_dof ==
                        mesh::GlobalDofIndex{
                            expected_mesh_column} &&
                    entry.row_cell_global ==
                        mesh::GlobalEntityId{
                            row_cell_global} &&
                    entry.column_cell_global ==
                        mesh::GlobalEntityId{
                            column_cell_global} &&
                    entry.component_row ==
                        component &&
                    entry.natural_variable_column ==
                        column &&
                    entry.block_kind ==
                        (diagonal
                             ? fdp::
                                   ComponentJacobianCellBlockKind3D::
                                       diagonal_cell
                             : fdp::
                                   ComponentJacobianCellBlockKind3D::
                                       off_diagonal_cell),
                "Jacobian scalar row/column mapping mismatch");
            near_collective(
                entry.value,
                expected_value);
        }
    }

    std::uint64_t local_zero_count = 0U;
    for (const auto& entry :
         result.jacobian_entries()) {
        if (entry.value == 0.0) {
            ++local_zero_count;
        }
    }
    std::uint64_t global_zero_count = 0U;
    if (MPI_Allreduce(
            &local_zero_count,
            &global_zero_count,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            PETSC_COMM_WORLD) !=
        MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI_Allreduce failed for zero-entry audit");
    }
    require_collective(
        global_zero_count > 0U,
        "zero Jacobian scalar was dropped from assembly-ready entries");

    {
        const auto bad_layout =
            make_natural_variable_dof_layout(
                rank,
                q - 1U);
        const auto bad_numbering =
            make_reversed_mesh_dof_numbering(
                rank,
                bad_layout,
                partition);
        std::optional<
            fdp::ComponentConservationGlobalAssemblyEntries3D>
            bad_output;
        const PetscErrorCode bad_error =
            fdp::
                make_component_conservation_global_assembly_entries_3d(
                    PETSC_COMM_WORLD,
                    *conservation,
                    partition,
                    bad_layout,
                    bad_numbering,
                    cell_bridge,
                    cell_pattern,
                    "natural_state",
                    &bad_output);
        require_collective(
            bad_error ==
                    PETSC_ERR_ARG_INCOMP &&
                !bad_output.has_value(),
            "wrong natural-variable DoF width was not rejected collectively");
    }

    {
        auto bad_pattern =
            make_incomplete_cell_column_pattern(
                rank);
        std::optional<
            fdp::ComponentConservationGlobalAssemblyEntries3D>
            bad_output;
        const PetscErrorCode bad_error =
            fdp::
                make_component_conservation_global_assembly_entries_3d(
                    PETSC_COMM_WORLD,
                    *conservation,
                    partition,
                    layout,
                    numbering,
                    cell_bridge,
                    bad_pattern,
                    "natural_state",
                    &bad_output);
        require_collective(
            bad_error ==
                    PETSC_ERR_ARG_INCOMP &&
                !bad_output.has_value(),
            "incomplete cell structural pattern was not rejected collectively");
    }
}

flow::FugacityEquilibriumResidualLinearization3P
make_fugacity_linearization(
    std::uint64_t stable_cell) {
    const auto state =
        identity(stable_cell);
    const auto chart =
        state.layout;
    const std::size_t n =
        chart.component_count();
    const std::size_t q =
        chart.unknown_count();
    const std::size_t rows =
        2U * n;
    const double factor =
        stable_cell == UINT64_C(10)
            ? 1.0
            : 1.5;

    std::vector<double> values(
        rows,
        0.0);
    std::vector<double> derivatives(
        rows * q,
        0.0);

    for (std::size_t row = 0U;
         row < rows;
         ++row) {
        values[row] =
            factor *
            0.01 *
            static_cast<double>(
                row + 1U);
        for (std::size_t column = 0U;
             column < q;
             ++column) {
            // pc=none local fugacity rows currently do not depend on the two
            // independent saturation coordinates. Keep those exact zeros in
            // the dense assembly-ready block.
            if (column == 2U ||
                column == 3U) {
                derivatives[
                    row * q +
                    column] =
                    0.0;
            } else {
                derivatives[
                    row * q +
                    column] =
                    factor *
                    0.001 *
                    static_cast<double>(
                        1U +
                        row * 100U +
                        column);
            }
        }
    }

    return {
        chart,
        state.component_ids,
        flow::FugacityEquilibriumResidual3P<double>{
            n,
            std::move(values)},
        q,
        std::move(derivatives)};
}

void fugacity_global_assembly_mapping() {
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
        "fugacity global mapping test requires exactly two ranks");

    auto partition =
        make_partition(rank);
    auto schedule =
        make_schedule(rank);
    LocalInputs inputs{rank};

    PetscErrorCode error =
        PETSC_SUCCESS;
    auto conservation =
        build(
            rank,
            partition,
            schedule,
            inputs,
            &error);
    require_collective(
        error == PETSC_SUCCESS &&
            conservation.has_value(),
        "distributed snapshot for fugacity mapping failed");

    const std::size_t q =
        identity(10U)
            .layout.unknown_count();
    auto layout =
        make_natural_variable_dof_layout(
            rank,
            q);
    auto numbering =
        make_reversed_mesh_dof_numbering(
            rank,
            layout,
            partition);
    auto cell_bridge =
        make_cell_row_bridge(
            rank);
    auto cell_pattern =
        make_cell_column_pattern(
            rank);

    const std::uint64_t stable_cell =
        rank == 0
            ? UINT64_C(10)
            : UINT64_C(20);
    auto fugacity =
        make_fugacity_linearization(
            stable_cell);

    const std::array<
        fdp::OwnedCellFugacityEquilibriumLinearizationBinding3D,
        1>
        bindings{{
            {
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{
                    stable_cell},
                &fugacity}
        }};

    std::optional<
        fdp::FugacityEquilibriumGlobalAssemblyEntries3D>
        mapping;
    error =
        fdp::
            make_fugacity_equilibrium_global_assembly_entries_3d(
                PETSC_COMM_WORLD,
                *conservation,
                partition,
                layout,
                numbering,
                cell_bridge,
                cell_pattern,
                "natural_state",
                bindings,
                &mapping);

    require_collective(
        error == PETSC_SUCCESS &&
            mapping.has_value(),
        "fugacity global assembly mapping failed");

    const auto& result =
        *mapping;
    const PetscInt expected_start =
        static_cast<PetscInt>(
            rank) *
        static_cast<PetscInt>(
            q);
    const PetscInt expected_end =
        expected_start +
        static_cast<PetscInt>(
            q);
    const PetscInt expected_count =
        2 *
        static_cast<PetscInt>(
            q);

    require_collective(
        result.petsc_scalar_row_start() ==
                expected_start &&
            result.petsc_scalar_row_end() ==
                expected_end &&
            result.petsc_scalar_row_count() ==
                expected_count &&
            result.component_count() == 3U &&
            result.natural_variable_count() ==
                q &&
            result.residual_entries().size() ==
                6U &&
            result.jacobian_entries().size() ==
                6U * q,
        "fugacity global mapping range/cardinality mismatch");

    const PetscInt energy_row =
        expected_start + 3;
    require_collective(
        std::none_of(
            result.residual_entries().begin(),
            result.residual_entries().end(),
            [energy_row](const auto& entry) {
                return entry.petsc_global_row ==
                    energy_row;
            }),
        "fugacity mapping occupied the reserved energy row");

    for (std::size_t local_row = 0U;
         local_row < 6U;
         ++local_row) {
        const std::size_t phase =
            1U +
            local_row / 3U;
        const std::size_t component =
            local_row % 3U;
        const auto slot =
            static_cast<flow::PhaseSlot3>(
                phase);
        const std::size_t equation_slot =
            fugacity.equation_index(
                slot,
                component);
        require_collective(
            equation_slot ==
                4U + local_row,
            "fugacity equation slot ordering changed");

        const auto& residual =
            result.residual_entries()[
                local_row];
        const PetscInt expected_row =
            expected_start +
            static_cast<PetscInt>(
                equation_slot);
        const std::uint64_t
            expected_mesh_row =
                (stable_cell ==
                         UINT64_C(10)
                     ? static_cast<
                           std::uint64_t>(
                               q)
                     : UINT64_C(0)) +
                static_cast<std::uint64_t>(
                    equation_slot);

        require_collective(
            residual.petsc_global_row ==
                    expected_row &&
                residual.mesh_global_row_dof ==
                    mesh::GlobalDofIndex{
                        expected_mesh_row} &&
                residual.row_cell_global ==
                    mesh::GlobalEntityId{
                        stable_cell} &&
                residual.non_reference_phase ==
                    slot &&
                residual.component ==
                    component,
            "fugacity residual row/provenance mapping mismatch");
        near_collective(
            residual.value,
            fugacity.residual(
                slot,
                component));

        require_collective(
            static_cast<std::uint64_t>(
                residual.petsc_global_row) !=
                residual
                    .mesh_global_row_dof
                    .value(),
            "fugacity mesh-global DoF was conflated with PETSc row");

        for (std::size_t column = 0U;
             column < q;
             ++column) {
            const std::size_t entry_index =
                local_row * q +
                column;
            const auto& entry =
                result.jacobian_entries()[
                    entry_index];
            const PetscInt expected_column =
                expected_start +
                static_cast<PetscInt>(
                    column);
            const std::uint64_t
                expected_mesh_column =
                    (stable_cell ==
                             UINT64_C(10)
                         ? static_cast<
                               std::uint64_t>(
                                   q)
                         : UINT64_C(0)) +
                    static_cast<std::uint64_t>(
                        column);

            require_collective(
                entry.petsc_global_row ==
                        expected_row &&
                    entry.petsc_global_column ==
                        expected_column &&
                    entry.mesh_global_row_dof ==
                        residual
                            .mesh_global_row_dof &&
                    entry.mesh_global_column_dof ==
                        mesh::GlobalDofIndex{
                            expected_mesh_column} &&
                    entry.row_cell_global ==
                        mesh::GlobalEntityId{
                            stable_cell} &&
                    entry.non_reference_phase ==
                        slot &&
                    entry.component ==
                        component &&
                    entry.natural_variable_column ==
                        column,
                "fugacity Jacobian scalar mapping mismatch");
            near_collective(
                entry.value,
                fugacity.d_residual(
                    slot,
                    component,
                    column));
        }
    }

    std::uint64_t local_zero_count = 0U;
    for (const auto& entry :
         result.jacobian_entries()) {
        if (entry.value == 0.0) {
            ++local_zero_count;
        }
    }
    std::uint64_t global_zero_count = 0U;
    if (MPI_Allreduce(
            &local_zero_count,
            &global_zero_count,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            PETSC_COMM_WORLD) !=
        MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI_Allreduce failed for fugacity zero-entry audit");
    }
    require_collective(
        global_zero_count == 24U,
        "fugacity zero saturation-column entries were not retained");

    {
        auto wrong_chart =
            make_fugacity_linearization(
                rank == 0
                    ? UINT64_C(20)
                    : UINT64_C(10));
        const std::array<
            fdp::OwnedCellFugacityEquilibriumLinearizationBinding3D,
            1>
            bad_bindings{{
                {
                    mesh::LocalIndex{0U},
                    mesh::GlobalEntityId{
                        stable_cell},
                    &wrong_chart}
            }};
        std::optional<
            fdp::FugacityEquilibriumGlobalAssemblyEntries3D>
            bad_output;
        const PetscErrorCode bad_error =
            fdp::
                make_fugacity_equilibrium_global_assembly_entries_3d(
                    PETSC_COMM_WORLD,
                    *conservation,
                    partition,
                    layout,
                    numbering,
                    cell_bridge,
                    cell_pattern,
                    "natural_state",
                    bad_bindings,
                    &bad_output);
        require_collective(
            bad_error ==
                    PETSC_ERR_ARG_INCOMP &&
                !bad_output.has_value(),
            "fugacity pivot/chart mismatch was not rejected collectively");
    }

    {
        const std::span<
            const fdp::
                OwnedCellFugacityEquilibriumLinearizationBinding3D>
            missing{};
        std::optional<
            fdp::FugacityEquilibriumGlobalAssemblyEntries3D>
            bad_output;
        const PetscErrorCode bad_error =
            fdp::
                make_fugacity_equilibrium_global_assembly_entries_3d(
                    PETSC_COMM_WORLD,
                    *conservation,
                    partition,
                    layout,
                    numbering,
                    cell_bridge,
                    cell_pattern,
                    "natural_state",
                    missing,
                    &bad_output);
        require_collective(
            bad_error ==
                    PETSC_ERR_ARG_INCOMP &&
                !bad_output.has_value(),
            "missing fugacity binding was not rejected collectively");
    }
}

flow::BackwardEulerEnergyAccumulationResidual3P
energy_accumulation(
    std::uint64_t stable_cell) {
    const auto state =
        identity(stable_cell);
    const std::size_t q =
        state.layout.unknown_count();
    std::vector<double> gradient(
        q,
        0.0);
    const double scale =
        stable_cell == UINT64_C(10)
            ? 10.0
            : 20.0;
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        gradient[column] =
            scale *
            static_cast<double>(
                column + 1U);
    }
    return {
        state,
        0.25,
        10.0,
        stable_cell == UINT64_C(10)
            ? 100.0
            : 200.0,
        q,
        std::move(gradient)};
}

fd::NormalizedEnergyFaceContributionLinearization3D
normalized_energy_face() {
    const auto owner =
        identity(UINT64_C(10));
    const auto neighbour =
        identity(UINT64_C(20));
    const std::size_t owner_q =
        owner.layout.unknown_count();
    const std::size_t neighbour_q =
        neighbour.layout.unknown_count();

    std::vector<double> owner_source_gradient(
        owner_q,
        0.0);
    std::vector<double> neighbour_source_gradient(
        neighbour_q,
        0.0);
    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        owner_source_gradient[column] =
            3.0 *
            static_cast<double>(
                column + 1U);
    }
    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        neighbour_source_gradient[column] =
            -4.0 *
            static_cast<double>(
                column + 1U);
    }

    std::vector<double> owner_owner(
        owner_q,
        0.0);
    std::vector<double> neighbour_owner(
        owner_q,
        0.0);
    std::vector<double> owner_neighbour(
        neighbour_q,
        0.0);
    std::vector<double> neighbour_neighbour(
        neighbour_q,
        0.0);

    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        owner_owner[column] =
            owner_source_gradient[column] /
            2.0;
        neighbour_owner[column] =
            -owner_source_gradient[column] /
            5.0;
    }
    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        owner_neighbour[column] =
            neighbour_source_gradient[column] /
            2.0;
        neighbour_neighbour[column] =
            -neighbour_source_gradient[column] /
            5.0;
    }

    return {
        mesh::LocalIndex{0U},
        {2.0, 5.0},
        owner,
        neighbour,
        300.0,
        -120.0,
        std::move(owner_owner),
        std::move(owner_neighbour),
        std::move(neighbour_owner),
        std::move(neighbour_neighbour)};
}

struct EnergyLocalInputs {
    flow::BackwardEulerEnergyAccumulationResidual3P
        accumulation10;
    flow::BackwardEulerEnergyAccumulationResidual3P
        accumulation20;
    std::vector<
        fdp::DistributedEnergyCellStateBinding3D>
        cells;
    fd::NormalizedEnergyFaceContributionLinearization3D
        face;
    std::vector<
        fdp::AuthoritativeNormalizedEnergyFaceBinding3D>
        faces;

    explicit EnergyLocalInputs(int rank)
        : accumulation10(
              energy_accumulation(
                  UINT64_C(10))),
          accumulation20(
              energy_accumulation(
                  UINT64_C(20))),
          face(
              normalized_energy_face()) {
        if (rank == 0) {
            cells = {
                {
                    mesh::LocalIndex{0U},
                    mesh::GlobalEntityId{10U},
                    2.0,
                    identity(UINT64_C(10)),
                    &accumulation10},
                {
                    mesh::LocalIndex{1U},
                    mesh::GlobalEntityId{20U},
                    5.0,
                    identity(UINT64_C(20)),
                    nullptr}
            };
            faces = {
                {
                    mesh::LocalIndex{0U},
                    mesh::GlobalEntityId{100U},
                    &face}
            };
        } else {
            cells = {
                {
                    mesh::LocalIndex{0U},
                    mesh::GlobalEntityId{20U},
                    5.0,
                    identity(UINT64_C(20)),
                    &accumulation20},
                {
                    mesh::LocalIndex{1U},
                    mesh::GlobalEntityId{10U},
                    2.0,
                    identity(UINT64_C(10)),
                    nullptr}
            };
        }
    }
};

std::optional<
    fdp::DistributedOwnedMultiCellEnergyConservationSnapshot3D>
build_energy(
    mesh::PartitionSnapshot& partition,
    dp::ParallelOwnedConnectionSchedule3D& schedule,
    EnergyLocalInputs& inputs,
    PetscErrorCode* error_out) {
    std::optional<
        fdp::DistributedOwnedMultiCellEnergyConservationSnapshot3D>
        output;
    const PetscErrorCode error =
        fdp::
            make_distributed_owned_energy_conservation_snapshot_3d(
                PETSC_COMM_WORLD,
                schedule,
                partition,
                inputs.cells,
                inputs.faces,
                &output);
    if (error_out != nullptr) {
        *error_out = error;
    }
    return output;
}

const fdp::OwnedCellPairEnergyJacobianBlock3D*
find_energy_block(
    const fdp::OwnedCellEnergyConservationRow3D& row,
    std::uint64_t column_global) {
    const auto found =
        std::find_if(
            row.off_diagonal_cell_pair_blocks.begin(),
            row.off_diagonal_cell_pair_blocks.end(),
            [column_global](const auto& block) {
                return block.column_cell_global ==
                    mesh::GlobalEntityId{
                        column_global};
            });
    return found ==
               row.off_diagonal_cell_pair_blocks.end()
        ? nullptr
        : &*found;
}

void audit_cross_rank_energy_jacobian(
    const fdp::
        DistributedOwnedMultiCellEnergyConservationSnapshot3D&
            snapshot,
    const EnergyLocalInputs& inputs) {
    const auto& row =
        snapshot.owned_rows().front();
    const auto& owned_binding =
        inputs.cells.front();
    const auto* owned_accumulation =
        owned_binding.owned_accumulation;
    require_collective(
        owned_accumulation != nullptr,
        "owned energy accumulation missing");

    for (const std::uint64_t source_global :
         {UINT64_C(10), UINT64_C(20)}) {
        const std::size_t q =
            identity(source_global)
                .layout.unknown_count();

        for (std::size_t column = 0U;
             column < q;
             ++column) {
            const bool source_row =
                row.cell_global ==
                mesh::GlobalEntityId{
                    source_global};
            const auto* block =
                source_row
                    ? nullptr
                    : find_energy_block(
                          row,
                          source_global);

            require_collective(
                source_row ||
                    block != nullptr,
                "missing cross-rank energy off-diagonal block");

            double local = 0.0;
            if (source_row) {
                local =
                    owned_binding.bulk_volume_m3 *
                    (row.local_residual
                         .d_local(column) -
                     owned_accumulation
                         ->d_residual(column));
            } else {
                local =
                    owned_binding.bulk_volume_m3 *
                    block->d_residual(
                        column);
            }

            double global = 0.0;
            if (MPI_Allreduce(
                    &local,
                    &global,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    PETSC_COMM_WORLD) !=
                MPI_SUCCESS) {
                throw std::runtime_error(
                    "MPI_Allreduce failed in energy Jacobian audit");
            }
            near_collective(
                global,
                0.0,
                0.0,
                2.0e-13);
        }
    }
}

void distributed_energy_exchange() {
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
        "distributed energy exchange test requires exactly two ranks");

    auto partition =
        make_partition(rank);
    auto schedule =
        make_schedule(rank);
    EnergyLocalInputs inputs{rank};

    PetscErrorCode error =
        PETSC_SUCCESS;
    auto output =
        build_energy(
            partition,
            schedule,
            inputs,
            &error);

    require_collective(
        error == PETSC_SUCCESS &&
            output.has_value(),
        "distributed owner-targeted energy exchange failed");

    const auto& snapshot =
        *output;
    require_collective(
        snapshot.owned_rows().size() == 1U &&
            snapshot.received_face_side_count() == 1U &&
            snapshot.local_authoritative_face_count() ==
                (rank == 0 ? 1U : 0U),
        "distributed energy snapshot cardinality mismatch");

    const auto& row =
        snapshot.owned_rows().front();
    const std::uint64_t expected_cell =
        rank == 0
            ? UINT64_C(10)
            : UINT64_C(20);
    const std::uint64_t expected_column =
        rank == 0
            ? UINT64_C(20)
            : UINT64_C(10);

    require_collective(
        row.cell ==
                mesh::LocalIndex{0U} &&
            row.cell_global ==
                mesh::GlobalEntityId{
                    expected_cell} &&
            row.local_residual
                    .neighbour_blocks.size() ==
                1U &&
            row.off_diagonal_cell_pair_blocks.size() ==
                1U,
        "owned energy row identity/local reordering mismatch");

    const auto& block =
        row.off_diagonal_cell_pair_blocks.front();
    require_collective(
        block.column_cell ==
                mesh::LocalIndex{1U} &&
            block.column_cell_global ==
                mesh::GlobalEntityId{
                    expected_column} &&
            block.contributing_face_global_ids ==
                std::vector<mesh::GlobalEntityId>{
                    mesh::GlobalEntityId{100U}},
        "stable-ID energy off-diagonal block mapping mismatch");

    const auto& accumulation_ref =
        rank == 0
            ? inputs.accumulation10
            : inputs.accumulation20;
    const double spatial =
        rank == 0
            ? inputs.face
                  .owner_contribution_w_per_bulk_m3
            : inputs.face
                  .neighbour_contribution_w_per_bulk_m3;

    near_collective(
        row.local_residual
            .residual_w_per_bulk_m3,
        accumulation_ref
                .residual_w_per_bulk_m3 +
            spatial);
    near_collective(
        snapshot
            .global_volume_weighted_spatial_balance_w(),
        0.0,
        0.0,
        2.0e-13);

    audit_cross_rank_energy_jacobian(
        snapshot,
        inputs);
}

void energy_global_assembly_mapping() {
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
        "energy global mapping test requires exactly two ranks");

    auto partition =
        make_partition(rank);
    auto schedule =
        make_schedule(rank);
    EnergyLocalInputs energy_inputs{rank};

    // Force one exact-zero final diagonal scalar on rank 0. The mapping must
    // retain it because matrix structure is defined by the cell pattern.
    if (rank == 0) {
        energy_inputs.accumulation10
            .gradient[0U] =
            -energy_inputs.face
                 .owner_row_owner_column_gradient[
                     0U];
    }

    PetscErrorCode error =
        PETSC_SUCCESS;
    auto energy_snapshot =
        build_energy(
            partition,
            schedule,
            energy_inputs,
            &error);
    require_collective(
        error == PETSC_SUCCESS &&
            energy_snapshot.has_value(),
        "distributed energy snapshot for global mapping failed");

    const std::size_t q =
        identity(UINT64_C(10))
            .layout.unknown_count();
    auto layout =
        make_natural_variable_dof_layout(
            rank,
            q);
    auto numbering =
        make_reversed_mesh_dof_numbering(
            rank,
            layout,
            partition);
    auto cell_bridge =
        make_cell_row_bridge(
            rank);
    auto cell_pattern =
        make_cell_column_pattern(
            rank);

    std::optional<
        fdp::EnergyConservationGlobalAssemblyEntries3D>
        mapping;
    error =
        fdp::
            make_energy_conservation_global_assembly_entries_3d(
                PETSC_COMM_WORLD,
                *energy_snapshot,
                partition,
                layout,
                numbering,
                cell_bridge,
                cell_pattern,
                "natural_state",
                &mapping);

    require_collective(
        error == PETSC_SUCCESS &&
            mapping.has_value(),
        "energy global assembly mapping failed");

    const auto& result =
        *mapping;
    const PetscInt expected_start =
        static_cast<PetscInt>(rank) *
        static_cast<PetscInt>(q);
    const PetscInt expected_end =
        expected_start +
        static_cast<PetscInt>(q);
    const PetscInt expected_count =
        2 *
        static_cast<PetscInt>(q);
    constexpr std::size_t component_count = 3U;
    const std::size_t energy_slot =
        component_count;

    require_collective(
        result.petsc_scalar_row_start() ==
                expected_start &&
            result.petsc_scalar_row_end() ==
                expected_end &&
            result.petsc_scalar_row_count() ==
                expected_count &&
            result.component_count() ==
                component_count &&
            result.natural_variable_count() ==
                q &&
            result.residual_entries().size() ==
                1U &&
            result.jacobian_entries().size() ==
                2U * q,
        "energy global mapping range/cardinality mismatch");

    const auto& owned_row =
        energy_snapshot->owned_rows().front();
    const auto* off_block =
        find_energy_block(
            owned_row,
            rank == 0
                ? UINT64_C(20)
                : UINT64_C(10));
    require_collective(
        off_block != nullptr,
        "energy global mapping fixture missing off-diagonal block");

    const std::uint64_t row_cell_global =
        rank == 0
            ? UINT64_C(10)
            : UINT64_C(20);
    const auto& residual =
        result.residual_entries().front();
    const PetscInt expected_row =
        expected_start +
        static_cast<PetscInt>(
            energy_slot);
    const std::uint64_t expected_mesh_row =
        (row_cell_global == UINT64_C(10)
             ? static_cast<std::uint64_t>(q)
             : UINT64_C(0)) +
        static_cast<std::uint64_t>(
            energy_slot);

    require_collective(
        residual.petsc_global_row ==
                expected_row &&
            residual.mesh_global_row_dof ==
                mesh::GlobalDofIndex{
                    expected_mesh_row} &&
            residual.row_cell_global ==
                mesh::GlobalEntityId{
                    row_cell_global},
        "energy residual row/provenance mapping mismatch");
    near_collective(
        residual.value_w_per_bulk_m3,
        owned_row.local_residual
            .residual_w_per_bulk_m3);

    require_collective(
        static_cast<std::uint64_t>(
            residual.petsc_global_row) !=
            residual.mesh_global_row_dof.value(),
        "energy mesh-global DoF was conflated with PETSc row");

    for (PetscInt global_column = 0;
         global_column < expected_count;
         ++global_column) {
        const auto& entry =
            result.jacobian_entries()[
                static_cast<std::size_t>(
                    global_column)];
        const std::size_t column =
            static_cast<std::size_t>(
                global_column %
                static_cast<PetscInt>(q));
        const std::uint64_t column_cell_global =
            global_column <
                    static_cast<PetscInt>(q)
                ? UINT64_C(10)
                : UINT64_C(20);
        const bool diagonal =
            column_cell_global ==
            row_cell_global;
        const double expected_value =
            diagonal
                ? owned_row.local_residual
                      .d_local(column)
                : off_block->d_residual(
                      column);
        const std::uint64_t expected_mesh_column =
            (column_cell_global ==
                     UINT64_C(10)
                 ? static_cast<std::uint64_t>(
                       q)
                 : UINT64_C(0)) +
            static_cast<std::uint64_t>(
                column);

        require_collective(
            entry.petsc_global_row ==
                    expected_row &&
                entry.petsc_global_column ==
                    global_column &&
                entry.mesh_global_row_dof ==
                    residual.mesh_global_row_dof &&
                entry.mesh_global_column_dof ==
                    mesh::GlobalDofIndex{
                        expected_mesh_column} &&
                entry.row_cell_global ==
                    mesh::GlobalEntityId{
                        row_cell_global} &&
                entry.column_cell_global ==
                    mesh::GlobalEntityId{
                        column_cell_global} &&
                entry.natural_variable_column ==
                    column &&
                entry.block_kind ==
                    (diagonal
                         ? fdp::
                               EnergyJacobianCellBlockKind3D::
                                   diagonal_cell
                         : fdp::
                               EnergyJacobianCellBlockKind3D::
                                   off_diagonal_cell),
            "energy Jacobian scalar mapping mismatch");
        near_collective(
            entry.value,
            expected_value);
    }

    std::uint64_t local_zero_count = 0U;
    for (const auto& entry :
         result.jacobian_entries()) {
        if (entry.value == 0.0) {
            ++local_zero_count;
        }
    }
    std::uint64_t global_zero_count = 0U;
    if (MPI_Allreduce(
            &local_zero_count,
            &global_zero_count,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            PETSC_COMM_WORLD) !=
        MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI_Allreduce failed for energy zero-entry audit");
    }
    require_collective(
        global_zero_count > 0U,
        "zero energy Jacobian scalar was dropped");

    {
        const auto bad_layout =
            make_natural_variable_dof_layout(
                rank,
                q - 1U);
        const auto bad_numbering =
            make_reversed_mesh_dof_numbering(
                rank,
                bad_layout,
                partition);
        std::optional<
            fdp::EnergyConservationGlobalAssemblyEntries3D>
            bad_output;
        const PetscErrorCode bad_error =
            fdp::
                make_energy_conservation_global_assembly_entries_3d(
                    PETSC_COMM_WORLD,
                    *energy_snapshot,
                    partition,
                    bad_layout,
                    bad_numbering,
                    cell_bridge,
                    cell_pattern,
                    "natural_state",
                    &bad_output);
        require_collective(
            bad_error ==
                    PETSC_ERR_ARG_INCOMP &&
                !bad_output.has_value(),
            "wrong energy natural-variable DoF width was not rejected collectively");
    }

    {
        auto bad_pattern =
            make_incomplete_cell_column_pattern(
                rank);
        std::optional<
            fdp::EnergyConservationGlobalAssemblyEntries3D>
            bad_output;
        const PetscErrorCode bad_error =
            fdp::
                make_energy_conservation_global_assembly_entries_3d(
                    PETSC_COMM_WORLD,
                    *energy_snapshot,
                    partition,
                    layout,
                    numbering,
                    cell_bridge,
                    bad_pattern,
                    "natural_state",
                    &bad_output);
        require_collective(
            bad_error ==
                    PETSC_ERR_ARG_INCOMP &&
                !bad_output.has_value(),
            "incomplete energy structural pattern was not rejected collectively");
    }

    // Final row-coverage gate: component + energy + fugacity residual rows must
    // occupy every scalar equation slot exactly once for the owned cell.
    LocalInputs component_inputs{rank};
    auto component_snapshot =
        build(
            rank,
            partition,
            schedule,
            component_inputs,
            &error);
    require_collective(
        error == PETSC_SUCCESS &&
            component_snapshot.has_value(),
        "component snapshot for full-row coverage failed");

    std::optional<
        fdp::ComponentConservationGlobalAssemblyEntries3D>
        component_mapping;
    error =
        fdp::
            make_component_conservation_global_assembly_entries_3d(
                PETSC_COMM_WORLD,
                *component_snapshot,
                partition,
                layout,
                numbering,
                cell_bridge,
                cell_pattern,
                "natural_state",
                &component_mapping);
    require_collective(
        error == PETSC_SUCCESS &&
            component_mapping.has_value(),
        "component mapping for full-row coverage failed");

    auto fugacity =
        make_fugacity_linearization(
            row_cell_global);
    const std::array<
        fdp::OwnedCellFugacityEquilibriumLinearizationBinding3D,
        1>
        fugacity_bindings{{
            {
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{
                    row_cell_global},
                &fugacity}
        }};
    std::optional<
        fdp::FugacityEquilibriumGlobalAssemblyEntries3D>
        fugacity_mapping;
    error =
        fdp::
            make_fugacity_equilibrium_global_assembly_entries_3d(
                PETSC_COMM_WORLD,
                *component_snapshot,
                partition,
                layout,
                numbering,
                cell_bridge,
                cell_pattern,
                "natural_state",
                fugacity_bindings,
                &fugacity_mapping);
    require_collective(
        error == PETSC_SUCCESS &&
            fugacity_mapping.has_value(),
        "fugacity mapping for full-row coverage failed");

    std::vector<PetscInt> row_slots;
    row_slots.reserve(q);
    for (const auto& entry :
         component_mapping->residual_entries()) {
        row_slots.push_back(
            entry.petsc_global_row -
            expected_start);
    }
    for (const auto& entry :
         result.residual_entries()) {
        row_slots.push_back(
            entry.petsc_global_row -
            expected_start);
    }
    for (const auto& entry :
         fugacity_mapping->residual_entries()) {
        row_slots.push_back(
            entry.petsc_global_row -
            expected_start);
    }
    std::sort(
        row_slots.begin(),
        row_slots.end());

    require_collective(
        row_slots.size() == q,
        "complete equation block row count mismatch");
    for (std::size_t slot = 0U;
         slot < q;
         ++slot) {
        require_collective(
            row_slots[slot] ==
                static_cast<PetscInt>(slot),
            "component+energy+fugacity rows do not cover the full natural-variable block exactly once");
    }
}

void complete_natural_variable_assembly_snapshot() {
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
        "complete assembly snapshot test requires exactly two ranks");

    auto partition =
        make_partition(rank);
    auto schedule =
        make_schedule(rank);

    LocalInputs component_inputs{rank};
    PetscErrorCode error =
        PETSC_SUCCESS;
    auto component_snapshot =
        build(
            rank,
            partition,
            schedule,
            component_inputs,
            &error);
    require_collective(
        error == PETSC_SUCCESS &&
            component_snapshot.has_value(),
        "component snapshot for complete assembly failed");

    EnergyLocalInputs energy_inputs{rank};
    auto energy_snapshot =
        build_energy(
            partition,
            schedule,
            energy_inputs,
            &error);
    require_collective(
        error == PETSC_SUCCESS &&
            energy_snapshot.has_value(),
        "energy snapshot for complete assembly failed");

    const std::size_t q =
        identity(UINT64_C(10))
            .layout.unknown_count();
    constexpr std::size_t n = 3U;
    auto layout =
        make_natural_variable_dof_layout(
            rank,
            q);
    auto numbering =
        make_reversed_mesh_dof_numbering(
            rank,
            layout,
            partition);
    auto cell_bridge =
        make_cell_row_bridge(
            rank);
    auto cell_pattern =
        make_cell_column_pattern(
            rank);

    std::optional<
        fdp::ComponentConservationGlobalAssemblyEntries3D>
        component_mapping;
    error =
        fdp::
            make_component_conservation_global_assembly_entries_3d(
                PETSC_COMM_WORLD,
                *component_snapshot,
                partition,
                layout,
                numbering,
                cell_bridge,
                cell_pattern,
                "natural_state",
                &component_mapping);
    require_collective(
        error == PETSC_SUCCESS &&
            component_mapping.has_value(),
        "component mapping for complete assembly failed");

    std::optional<
        fdp::EnergyConservationGlobalAssemblyEntries3D>
        energy_mapping;
    error =
        fdp::
            make_energy_conservation_global_assembly_entries_3d(
                PETSC_COMM_WORLD,
                *energy_snapshot,
                partition,
                layout,
                numbering,
                cell_bridge,
                cell_pattern,
                "natural_state",
                &energy_mapping);
    require_collective(
        error == PETSC_SUCCESS &&
            energy_mapping.has_value(),
        "energy mapping for complete assembly failed");

    const std::uint64_t stable_cell =
        rank == 0
            ? UINT64_C(10)
            : UINT64_C(20);
    auto fugacity =
        make_fugacity_linearization(
            stable_cell);
    const std::array<
        fdp::OwnedCellFugacityEquilibriumLinearizationBinding3D,
        1>
        fugacity_bindings{{
            {
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{
                    stable_cell},
                &fugacity}
        }};
    std::optional<
        fdp::FugacityEquilibriumGlobalAssemblyEntries3D>
        fugacity_mapping;
    error =
        fdp::
            make_fugacity_equilibrium_global_assembly_entries_3d(
                PETSC_COMM_WORLD,
                *component_snapshot,
                partition,
                layout,
                numbering,
                cell_bridge,
                cell_pattern,
                "natural_state",
                fugacity_bindings,
                &fugacity_mapping);
    require_collective(
        error == PETSC_SUCCESS &&
            fugacity_mapping.has_value(),
        "fugacity mapping for complete assembly failed");

    std::optional<
        fdp::CompleteNaturalVariableAssemblySnapshot3D>
        complete;
    error =
        fdp::
            make_complete_natural_variable_assembly_snapshot_3d(
                PETSC_COMM_WORLD,
                *component_mapping,
                *energy_mapping,
                *fugacity_mapping,
                partition,
                cell_bridge,
                cell_pattern,
                &complete);
    require_collective(
        error == PETSC_SUCCESS &&
            complete.has_value(),
        "complete natural-variable assembly snapshot failed");

    const PetscInt expected_start =
        static_cast<PetscInt>(rank) *
        static_cast<PetscInt>(q);
    const PetscInt expected_end =
        expected_start +
        static_cast<PetscInt>(q);
    const PetscInt global_count =
        2 *
        static_cast<PetscInt>(q);

    require_collective(
        complete->petsc_scalar_row_start() ==
                expected_start &&
            complete->petsc_scalar_row_end() ==
                expected_end &&
            complete->petsc_scalar_row_count() ==
                global_count &&
            complete->component_count() == n &&
            complete->natural_variable_count() ==
                q &&
            complete->residual_entries().size() ==
                q &&
            complete->jacobian_entries().size() ==
                14U * q,
        "complete assembly cardinality/metadata mismatch");

    std::size_t component_rows = 0U;
    std::size_t energy_rows = 0U;
    std::size_t fugacity_rows = 0U;
    for (std::size_t slot = 0U;
         slot < q;
         ++slot) {
        const auto& residual =
            complete->residual_entries()[slot];
        require_collective(
            residual.petsc_global_row ==
                    expected_start +
                        static_cast<PetscInt>(
                            slot) &&
                residual.equation_slot ==
                    slot &&
                residual.row_cell_global ==
                    mesh::GlobalEntityId{
                        stable_cell},
            "complete residual row coverage/stable identity mismatch");

        if (slot < n) {
            require_collective(
                residual.equation_kind ==
                    fdp::NaturalVariableEquationKind3D::
                        component_conservation,
                "component equation kind mismatch in complete snapshot");
            ++component_rows;
        } else if (slot == n) {
            require_collective(
                residual.equation_kind ==
                    fdp::NaturalVariableEquationKind3D::
                        energy_conservation,
                "energy equation kind mismatch in complete snapshot");
            ++energy_rows;
        } else {
            require_collective(
                residual.equation_kind ==
                    fdp::NaturalVariableEquationKind3D::
                        fugacity_equilibrium,
                "fugacity equation kind mismatch in complete snapshot");
            ++fugacity_rows;
        }
    }
    require_collective(
        component_rows == n &&
            energy_rows == 1U &&
            fugacity_rows == 2U * n,
        "complete residual equation-family counts changed");

    std::size_t jacobian_index = 0U;
    std::uint64_t local_zero_count = 0U;
    for (const auto& residual :
         complete->residual_entries()) {
        const std::size_t expected_columns =
            residual.equation_kind ==
                    fdp::NaturalVariableEquationKind3D::
                        fugacity_equilibrium
                ? q
                : 2U * q;
        std::size_t actual_columns = 0U;
        PetscInt previous_column = -1;
        while (jacobian_index <
                   complete->jacobian_entries().size() &&
               complete->jacobian_entries()[
                       jacobian_index]
                       .petsc_global_row ==
                   residual.petsc_global_row) {
            const auto& entry =
                complete->jacobian_entries()[
                    jacobian_index];
            require_collective(
                entry.petsc_global_column >=
                        0 &&
                    entry.petsc_global_column <
                        global_count &&
                    entry.petsc_global_column >
                        previous_column &&
                    entry.equation_slot ==
                        residual.equation_slot &&
                    entry.equation_kind ==
                        residual.equation_kind &&
                    entry.row_cell_global ==
                        residual.row_cell_global,
                "complete Jacobian ordering/provenance mismatch");
            if (entry.value == 0.0) {
                ++local_zero_count;
            }
            previous_column =
                entry.petsc_global_column;
            ++actual_columns;
            ++jacobian_index;
        }
        require_collective(
            actual_columns ==
                expected_columns,
            "complete Jacobian per-row structural cardinality mismatch");
    }
    require_collective(
        jacobian_index ==
            complete->jacobian_entries().size(),
        "complete Jacobian contained orphan entries");

    std::uint64_t global_zero_count = 0U;
    if (MPI_Allreduce(
            &local_zero_count,
            &global_zero_count,
            1,
            MPI_UINT64_T,
            MPI_SUM,
            PETSC_COMM_WORLD) !=
        MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI_Allreduce failed for complete assembly zero-entry audit");
    }
    require_collective(
        global_zero_count >= 24U,
        "complete assembly dropped exact-zero structural Jacobian entries");

    // Mapping metadata mismatch must reject collectively.
    {
        std::vector<
            fdp::AssemblyReadyEnergyResidualEntry3D>
            residuals{
                energy_mapping
                    ->residual_entries()
                    .begin(),
                energy_mapping
                    ->residual_entries()
                    .end()};
        std::vector<
            fdp::AssemblyReadyEnergyJacobianEntry3D>
            jacobians{
                energy_mapping
                    ->jacobian_entries()
                    .begin(),
                energy_mapping
                    ->jacobian_entries()
                    .end()};
        fdp::EnergyConservationGlobalAssemblyEntries3D
            bad_energy{
                energy_mapping->local_rank(),
                energy_mapping->rank_count(),
                "different_state",
                energy_mapping->component_count(),
                energy_mapping
                    ->natural_variable_count(),
                energy_mapping
                    ->petsc_scalar_row_start(),
                energy_mapping
                    ->petsc_scalar_row_end(),
                energy_mapping
                    ->petsc_scalar_row_count(),
                std::move(residuals),
                std::move(jacobians)};

        std::optional<
            fdp::CompleteNaturalVariableAssemblySnapshot3D>
            bad_output;
        const PetscErrorCode bad_error =
            fdp::
                make_complete_natural_variable_assembly_snapshot_3d(
                    PETSC_COMM_WORLD,
                    *component_mapping,
                    bad_energy,
                    *fugacity_mapping,
                    partition,
                    cell_bridge,
                    cell_pattern,
                    &bad_output);
        require_collective(
            bad_error ==
                    PETSC_ERR_ARG_INCOMP &&
                !bad_output.has_value(),
            "complete assembly metadata mismatch was not rejected collectively");
    }

    // A syntactically valid energy mapping that collides with component row 0
    // must be rejected by the complete-row coverage gate.
    {
        std::vector<
            fdp::AssemblyReadyEnergyResidualEntry3D>
            residuals{
                energy_mapping
                    ->residual_entries()
                    .begin(),
                energy_mapping
                    ->residual_entries()
                    .end()};
        std::vector<
            fdp::AssemblyReadyEnergyJacobianEntry3D>
            jacobians{
                energy_mapping
                    ->jacobian_entries()
                    .begin(),
                energy_mapping
                    ->jacobian_entries()
                    .end()};
        residuals.front().petsc_global_row =
            expected_start;
        for (auto& entry : jacobians) {
            entry.petsc_global_row =
                expected_start;
        }
        fdp::EnergyConservationGlobalAssemblyEntries3D
            colliding_energy{
                energy_mapping->local_rank(),
                energy_mapping->rank_count(),
                std::string{
                    energy_mapping
                        ->natural_variable_id()},
                energy_mapping->component_count(),
                energy_mapping
                    ->natural_variable_count(),
                energy_mapping
                    ->petsc_scalar_row_start(),
                energy_mapping
                    ->petsc_scalar_row_end(),
                energy_mapping
                    ->petsc_scalar_row_count(),
                std::move(residuals),
                std::move(jacobians)};

        std::optional<
            fdp::CompleteNaturalVariableAssemblySnapshot3D>
            bad_output;
        const PetscErrorCode bad_error =
            fdp::
                make_complete_natural_variable_assembly_snapshot_3d(
                    PETSC_COMM_WORLD,
                    *component_mapping,
                    colliding_energy,
                    *fugacity_mapping,
                    partition,
                    cell_bridge,
                    cell_pattern,
                    &bad_output);
        require_collective(
            bad_error ==
                    PETSC_ERR_ARG_INCOMP &&
                !bad_output.has_value(),
            "complete assembly row collision was not rejected collectively");
    }

    // Exact cell pattern is an input to the complete snapshot too; dropping
    // the remote neighbour block must be rejected.
    {
        auto bad_pattern =
            make_incomplete_cell_column_pattern(
                rank);
        std::optional<
            fdp::CompleteNaturalVariableAssemblySnapshot3D>
            bad_output;
        const PetscErrorCode bad_error =
            fdp::
                make_complete_natural_variable_assembly_snapshot_3d(
                    PETSC_COMM_WORLD,
                    *component_mapping,
                    *energy_mapping,
                    *fugacity_mapping,
                    partition,
                    cell_bridge,
                    bad_pattern,
                    &bad_output);
        require_collective(
            bad_error ==
                    PETSC_ERR_ARG_INCOMP &&
                !bad_output.has_value(),
            "complete assembly symbolic mismatch was not rejected collectively");
    }
}

void complete_natural_variable_petsc_materialization() {
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
        "PETSc materialization test requires exactly two ranks");

    auto partition =
        make_partition(rank);
    auto schedule =
        make_schedule(rank);

    LocalInputs component_inputs{rank};
    PetscErrorCode error =
        PETSC_SUCCESS;
    auto component_snapshot =
        build(
            rank,
            partition,
            schedule,
            component_inputs,
            &error);
    require_collective(
        error == PETSC_SUCCESS &&
            component_snapshot.has_value(),
        "component snapshot for PETSc materialization failed");

    EnergyLocalInputs energy_inputs{rank};
    auto energy_snapshot =
        build_energy(
            partition,
            schedule,
            energy_inputs,
            &error);
    require_collective(
        error == PETSC_SUCCESS &&
            energy_snapshot.has_value(),
        "energy snapshot for PETSc materialization failed");

    const std::size_t q =
        identity(UINT64_C(10))
            .layout.unknown_count();
    constexpr std::size_t n = 3U;
    auto layout =
        make_natural_variable_dof_layout(
            rank,
            q);
    auto numbering =
        make_reversed_mesh_dof_numbering(
            rank,
            layout,
            partition);
    auto cell_bridge =
        make_cell_row_bridge(
            rank);
    auto cell_pattern =
        make_cell_column_pattern(
            rank);

    std::optional<
        fdp::ComponentConservationGlobalAssemblyEntries3D>
        component_mapping;
    error =
        fdp::
            make_component_conservation_global_assembly_entries_3d(
                PETSC_COMM_WORLD,
                *component_snapshot,
                partition,
                layout,
                numbering,
                cell_bridge,
                cell_pattern,
                "natural_state",
                &component_mapping);
    require_collective(
        error == PETSC_SUCCESS &&
            component_mapping.has_value(),
        "component mapping for PETSc materialization failed");

    std::optional<
        fdp::EnergyConservationGlobalAssemblyEntries3D>
        energy_mapping;
    error =
        fdp::
            make_energy_conservation_global_assembly_entries_3d(
                PETSC_COMM_WORLD,
                *energy_snapshot,
                partition,
                layout,
                numbering,
                cell_bridge,
                cell_pattern,
                "natural_state",
                &energy_mapping);
    require_collective(
        error == PETSC_SUCCESS &&
            energy_mapping.has_value(),
        "energy mapping for PETSc materialization failed");

    const std::uint64_t stable_cell =
        rank == 0
            ? UINT64_C(10)
            : UINT64_C(20);
    auto fugacity =
        make_fugacity_linearization(
            stable_cell);
    const std::array<
        fdp::OwnedCellFugacityEquilibriumLinearizationBinding3D,
        1>
        fugacity_bindings{{
            {
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{
                    stable_cell},
                &fugacity}
        }};
    std::optional<
        fdp::FugacityEquilibriumGlobalAssemblyEntries3D>
        fugacity_mapping;
    error =
        fdp::
            make_fugacity_equilibrium_global_assembly_entries_3d(
                PETSC_COMM_WORLD,
                *component_snapshot,
                partition,
                layout,
                numbering,
                cell_bridge,
                cell_pattern,
                "natural_state",
                fugacity_bindings,
                &fugacity_mapping);
    require_collective(
        error == PETSC_SUCCESS &&
            fugacity_mapping.has_value(),
        "fugacity mapping for PETSc materialization failed");

    std::optional<
        fdp::CompleteNaturalVariableAssemblySnapshot3D>
        complete;
    error =
        fdp::
            make_complete_natural_variable_assembly_snapshot_3d(
                PETSC_COMM_WORLD,
                *component_mapping,
                *energy_mapping,
                *fugacity_mapping,
                partition,
                cell_bridge,
                cell_pattern,
                &complete);
    require_collective(
        error == PETSC_SUCCESS &&
            complete.has_value(),
        "complete snapshot for PETSc materialization failed");

    Vec residual = nullptr;
    Mat jacobian = nullptr;
    error =
        fdp::
            materialize_complete_natural_variable_petsc_system_3d(
                PETSC_COMM_WORLD,
                *complete,
                cell_bridge,
                &residual,
                &jacobian);
    require_collective(
        error == PETSC_SUCCESS &&
            residual != nullptr &&
            jacobian != nullptr,
        "PETSc complete-system materialization failed");

    const PetscInt expected_start =
        static_cast<PetscInt>(rank) *
        static_cast<PetscInt>(q);
    const PetscInt expected_end =
        expected_start +
        static_cast<PetscInt>(q);
    const PetscInt global_count =
        2 *
        static_cast<PetscInt>(q);

    PetscInt vec_local = -1;
    PetscInt vec_global = -1;
    PetscInt vec_start = -1;
    PetscInt vec_end = -1;
    PetscInt mat_local_rows = -1;
    PetscInt mat_local_columns = -1;
    PetscInt mat_global_rows = -1;
    PetscInt mat_global_columns = -1;
    PetscInt mat_start = -1;
    PetscInt mat_end = -1;
    PetscInt mat_column_start = -1;
    PetscInt mat_column_end = -1;
    PetscBool is_mpiaij =
        PETSC_FALSE;

    require_collective(
        VecGetLocalSize(
            residual,
            &vec_local) ==
                PETSC_SUCCESS &&
            VecGetSize(
                residual,
                &vec_global) ==
                PETSC_SUCCESS &&
            VecGetOwnershipRange(
                residual,
                &vec_start,
                &vec_end) ==
                PETSC_SUCCESS &&
            MatGetLocalSize(
                jacobian,
                &mat_local_rows,
                &mat_local_columns) ==
                PETSC_SUCCESS &&
            MatGetSize(
                jacobian,
                &mat_global_rows,
                &mat_global_columns) ==
                PETSC_SUCCESS &&
            MatGetOwnershipRange(
                jacobian,
                &mat_start,
                &mat_end) ==
                PETSC_SUCCESS &&
            MatGetOwnershipRangeColumn(
                jacobian,
                &mat_column_start,
                &mat_column_end) ==
                PETSC_SUCCESS &&
            PetscObjectTypeCompare(
                reinterpret_cast<PetscObject>(
                    jacobian),
                MATMPIAIJ,
                &is_mpiaij) ==
                PETSC_SUCCESS,
        "PETSc materialized system metadata query failed");

    require_collective(
        vec_local ==
                static_cast<PetscInt>(q) &&
            vec_global ==
                global_count &&
            vec_start ==
                expected_start &&
            vec_end ==
                expected_end &&
            mat_local_rows ==
                static_cast<PetscInt>(q) &&
            mat_local_columns ==
                static_cast<PetscInt>(q) &&
            mat_global_rows ==
                global_count &&
            mat_global_columns ==
                global_count &&
            mat_start ==
                expected_start &&
            mat_end ==
                expected_end &&
            mat_column_start ==
                expected_start &&
            mat_column_end ==
                expected_end &&
            is_mpiaij ==
                PETSC_TRUE,
        "PETSc materialized system size/ownership/type mismatch");

    for (const auto& entry :
         complete->residual_entries()) {
        const PetscInt row =
            entry.petsc_global_row;
        PetscScalar actual{};
        require_collective(
            VecGetValues(
                residual,
                1,
                &row,
                &actual) ==
                PETSC_SUCCESS,
            "VecGetValues failed for materialized residual");
        near_collective(
            static_cast<double>(
                PetscRealPart(actual)),
            entry.native_value);
    }

    for (const auto& entry :
         complete->jacobian_entries()) {
        const PetscInt row =
            entry.petsc_global_row;
        const PetscInt column =
            entry.petsc_global_column;
        PetscScalar actual{};
        require_collective(
            MatGetValues(
                jacobian,
                1,
                &row,
                1,
                &column,
                &actual) ==
                PETSC_SUCCESS,
            "MatGetValues failed for materialized Jacobian");
        near_collective(
            static_cast<double>(
                PetscRealPart(actual)),
            entry.value);
    }

    MatInfo info{};
    require_collective(
        MatGetInfo(
            jacobian,
            MAT_LOCAL,
            &info) ==
                PETSC_SUCCESS &&
            info.mallocs == 0.0,
        "PETSc materialized Jacobian required dynamic nonzero allocation");

    const auto zero_entry =
        std::find_if(
            complete->jacobian_entries().begin(),
            complete->jacobian_entries().end(),
            [](const auto& entry) {
                return entry.value == 0.0;
            });
    require_collective(
        zero_entry !=
            complete->jacobian_entries().end(),
        "PETSc materialization fixture lacks exact-zero structural entry");

    const PetscInt zero_row =
        zero_entry->petsc_global_row;
    const PetscInt zero_column =
        zero_entry->petsc_global_column;
    const PetscScalar probe_value =
        static_cast<PetscScalar>(
            1.234567);
    error =
        MatSetValues(
            jacobian,
            1,
            &zero_row,
            1,
            &zero_column,
            &probe_value,
            INSERT_VALUES);
    require_collective(
        error == PETSC_SUCCESS,
        "existing exact-zero Jacobian location was not materialized");
    const PetscErrorCode probe_begin_error =
        MatAssemblyBegin(
            jacobian,
            MAT_FINAL_ASSEMBLY);
    const PetscErrorCode probe_end_error =
        MatAssemblyEnd(
            jacobian,
            MAT_FINAL_ASSEMBLY);
    require_collective(
        probe_begin_error ==
                PETSC_SUCCESS &&
            probe_end_error ==
                PETSC_SUCCESS,
        "assembly after exact-zero location probe failed");

    PetscScalar probed{};
    require_collective(
        MatGetValues(
            jacobian,
            1,
            &zero_row,
            1,
            &zero_column,
            &probed) ==
                PETSC_SUCCESS,
        "readback after exact-zero location probe failed");
    near_collective(
        static_cast<double>(
            PetscRealPart(probed)),
        1.234567);

    const PetscScalar restore_zero =
        PetscScalar{0.0};
    const PetscErrorCode restore_set_error =
        MatSetValues(
            jacobian,
            1,
            &zero_row,
            1,
            &zero_column,
            &restore_zero,
            INSERT_VALUES);
    require_collective(
        restore_set_error ==
            PETSC_SUCCESS,
        "restoring exact-zero Jacobian value insertion failed");
    const PetscErrorCode restore_begin_error =
        MatAssemblyBegin(
            jacobian,
            MAT_FINAL_ASSEMBLY);
    const PetscErrorCode restore_end_error =
        MatAssemblyEnd(
            jacobian,
            MAT_FINAL_ASSEMBLY);
    require_collective(
        restore_begin_error ==
                PETSC_SUCCESS &&
            restore_end_error ==
                PETSC_SUCCESS,
        "restoring exact-zero Jacobian assembly failed");

    // Cell-level preallocation includes the neighbour block for every scalar
    // row, but fugacity physics does not materialize that location. The
    // location-freeze option must reject such an insertion.
    const PetscInt fugacity_row =
        expected_start +
        static_cast<PetscInt>(
            n + 1U);
    const PetscInt remote_cell_row =
        rank == 0
            ? 1
            : 0;
    const PetscInt forbidden_column =
        remote_cell_row *
        static_cast<PetscInt>(q);
    const PetscScalar forbidden_value =
        PetscScalar{7.0};

    require_collective(
        PetscPushErrorHandler(
            PetscReturnErrorHandler,
            nullptr) ==
            PETSC_SUCCESS,
        "failed to install PETSc return-error handler");
    const PetscErrorCode forbidden_error =
        MatSetValues(
            jacobian,
            1,
            &fugacity_row,
            1,
            &forbidden_column,
            &forbidden_value,
            INSERT_VALUES);
    require_collective(
        PetscPopErrorHandler() ==
            PETSC_SUCCESS,
        "failed to restore PETSc error handler");
    require_collective(
        forbidden_error !=
            PETSC_SUCCESS,
        "PETSc accepted a non-materialized fugacity neighbour location");

    // Materializer metadata mismatch must be collectively rejected before
    // creating any output objects.
    {
        const PetscInt cell_start =
            cell_bridge.global_row_start();
        const PetscInt cell_end =
            cell_bridge.global_row_end();
        auto bad_bridge =
            dp::PetscMpiAijSymbolicPreallocation3D{
                cell_bridge.local_rank(),
                cell_bridge.rank_count(),
                cell_start,
                cell_end,
                3,
                std::vector<mesh::LocalIndex>(
                    cell_bridge
                        .owned_cells_in_petsc_row_order()
                        .begin(),
                    cell_bridge
                        .owned_cells_in_petsc_row_order()
                        .end()),
                std::vector<mesh::GlobalEntityId>(
                    cell_bridge
                        .owned_cell_global_ids()
                        .begin(),
                    cell_bridge
                        .owned_cell_global_ids()
                        .end()),
                std::vector<PetscInt>(
                    cell_bridge
                        .owned_global_rows()
                        .begin(),
                    cell_bridge
                        .owned_global_rows()
                        .end()),
                std::vector<PetscInt>(
                    cell_bridge
                        .diagonal_nnz()
                        .begin(),
                    cell_bridge
                        .diagonal_nnz()
                        .end()),
                std::vector<PetscInt>(
                    cell_bridge
                        .off_diagonal_nnz()
                        .begin(),
                    cell_bridge
                        .off_diagonal_nnz()
                        .end()),
                rank == 0
                    ? std::vector<PetscInt>{0, 1}
                    : std::vector<PetscInt>{1, 0}};

        Vec bad_residual = nullptr;
        Mat bad_jacobian = nullptr;
        const PetscErrorCode bad_error =
            fdp::
                materialize_complete_natural_variable_petsc_system_3d(
                    PETSC_COMM_WORLD,
                    *complete,
                    bad_bridge,
                    &bad_residual,
                    &bad_jacobian);
        require_collective(
            bad_error ==
                    PETSC_ERR_ARG_INCOMP &&
                bad_residual == nullptr &&
                bad_jacobian == nullptr,
            "PETSc materializer metadata mismatch was not rejected collectively");
    }

    const PetscErrorCode vec_destroy_error =
        VecDestroy(
            &residual);
    const PetscErrorCode mat_destroy_error =
        MatDestroy(
            &jacobian);
    require_collective(
        vec_destroy_error ==
                PETSC_SUCCESS &&
            mat_destroy_error ==
                PETSC_SUCCESS,
        "PETSc materialized system destroy failed");
}


double manufactured_reference_delta(
    PetscInt global_scalar) {
    if (global_scalar < 0) {
        throw std::invalid_argument(
            "negative manufactured scalar index");
    }
    const double magnitude =
        0.01 *
        static_cast<double>(
            global_scalar + 1);
    return global_scalar % 2 == 0
        ? magnitude
        : -magnitude;
}

fdp::CompleteNaturalVariableAssemblySnapshot3D
make_manufactured_solvable_snapshot(
    const fdp::
        CompleteNaturalVariableAssemblySnapshot3D&
            source) {
    std::vector<
        fdp::CompleteNaturalVariableResidualEntry3D>
        residuals{
            source.residual_entries().begin(),
            source.residual_entries().end()};
    std::vector<
        fdp::CompleteNaturalVariableJacobianEntry3D>
        jacobian{
            source.jacobian_entries().begin(),
            source.jacobian_entries().end()};

    for (auto& entry : jacobian) {
        if (entry.petsc_global_row ==
            entry.petsc_global_column) {
            entry.value =
                10.0 +
                0.1 *
                    static_cast<double>(
                        entry.petsc_global_row +
                        1);
        } else {
            const PetscInt pattern =
                (entry.petsc_global_row +
                 entry.petsc_global_column) %
                7;
            entry.value =
                0.001 *
                static_cast<double>(
                    pattern + 1);
        }
    }

    for (auto& residual : residuals) {
        double action = 0.0;
        for (const auto& entry : jacobian) {
            if (entry.petsc_global_row !=
                residual.petsc_global_row) {
                continue;
            }
            action +=
                entry.value *
                manufactured_reference_delta(
                    entry.petsc_global_column);
        }
        residual.native_value =
            -action;
    }

    return {
        source.local_rank(),
        source.rank_count(),
        std::string{
            source.natural_variable_id()},
        source.component_count(),
        source.natural_variable_count(),
        source.petsc_scalar_row_start(),
        source.petsc_scalar_row_end(),
        source.petsc_scalar_row_count(),
        std::move(residuals),
        std::move(jacobian)};
}

void natural_variable_newton_linear_system() {
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
        "frozen Newton linear-system test requires exactly two ranks");

    auto partition =
        make_partition(rank);
    auto schedule =
        make_schedule(rank);

    LocalInputs component_inputs{rank};
    PetscErrorCode error =
        PETSC_SUCCESS;
    auto component_snapshot =
        build(
            rank,
            partition,
            schedule,
            component_inputs,
            &error);
    require_collective(
        error == PETSC_SUCCESS &&
            component_snapshot.has_value(),
        "component snapshot for frozen Newton test failed");

    EnergyLocalInputs energy_inputs{rank};
    auto energy_snapshot =
        build_energy(
            partition,
            schedule,
            energy_inputs,
            &error);
    require_collective(
        error == PETSC_SUCCESS &&
            energy_snapshot.has_value(),
        "energy snapshot for frozen Newton test failed");

    const std::size_t q =
        identity(UINT64_C(10))
            .layout.unknown_count();
    auto layout =
        make_natural_variable_dof_layout(
            rank,
            q);
    auto numbering =
        make_reversed_mesh_dof_numbering(
            rank,
            layout,
            partition);
    auto cell_bridge =
        make_cell_row_bridge(
            rank);
    auto cell_pattern =
        make_cell_column_pattern(
            rank);

    std::optional<
        fdp::ComponentConservationGlobalAssemblyEntries3D>
        component_mapping;
    error =
        fdp::
            make_component_conservation_global_assembly_entries_3d(
                PETSC_COMM_WORLD,
                *component_snapshot,
                partition,
                layout,
                numbering,
                cell_bridge,
                cell_pattern,
                "natural_state",
                &component_mapping);
    require_collective(
        error == PETSC_SUCCESS &&
            component_mapping.has_value(),
        "component mapping for frozen Newton test failed");

    std::optional<
        fdp::EnergyConservationGlobalAssemblyEntries3D>
        energy_mapping;
    error =
        fdp::
            make_energy_conservation_global_assembly_entries_3d(
                PETSC_COMM_WORLD,
                *energy_snapshot,
                partition,
                layout,
                numbering,
                cell_bridge,
                cell_pattern,
                "natural_state",
                &energy_mapping);
    require_collective(
        error == PETSC_SUCCESS &&
            energy_mapping.has_value(),
        "energy mapping for frozen Newton test failed");

    const std::uint64_t stable_cell =
        rank == 0
            ? UINT64_C(10)
            : UINT64_C(20);
    auto fugacity =
        make_fugacity_linearization(
            stable_cell);
    const std::array<
        fdp::OwnedCellFugacityEquilibriumLinearizationBinding3D,
        1>
        fugacity_bindings{{
            {
                mesh::LocalIndex{0U},
                mesh::GlobalEntityId{
                    stable_cell},
                &fugacity}
        }};
    std::optional<
        fdp::FugacityEquilibriumGlobalAssemblyEntries3D>
        fugacity_mapping;
    error =
        fdp::
            make_fugacity_equilibrium_global_assembly_entries_3d(
                PETSC_COMM_WORLD,
                *component_snapshot,
                partition,
                layout,
                numbering,
                cell_bridge,
                cell_pattern,
                "natural_state",
                fugacity_bindings,
                &fugacity_mapping);
    require_collective(
        error == PETSC_SUCCESS &&
            fugacity_mapping.has_value(),
        "fugacity mapping for frozen Newton test failed");

    std::optional<
        fdp::CompleteNaturalVariableAssemblySnapshot3D>
        complete;
    error =
        fdp::
            make_complete_natural_variable_assembly_snapshot_3d(
                PETSC_COMM_WORLD,
                *component_mapping,
                *energy_mapping,
                *fugacity_mapping,
                partition,
                cell_bridge,
                cell_pattern,
                &complete);
    require_collective(
        error == PETSC_SUCCESS &&
            complete.has_value(),
        "complete snapshot for frozen Newton test failed");

    auto manufactured =
        make_manufactured_solvable_snapshot(
            *complete);

    Vec residual = nullptr;
    Mat jacobian = nullptr;
    error =
        fdp::
            materialize_complete_natural_variable_petsc_system_3d(
                PETSC_COMM_WORLD,
                manufactured,
                cell_bridge,
                &residual,
                &jacobian);
    require_collective(
        error == PETSC_SUCCESS &&
            residual != nullptr &&
            jacobian != nullptr,
        "manufactured PETSc system materialization failed");

    Vec reference_delta = nullptr;
    Vec manufactured_balance = nullptr;
    const PetscErrorCode reference_duplicate_error =
        VecDuplicate(
            residual,
            &reference_delta);
    const PetscErrorCode balance_duplicate_error =
        VecDuplicate(
            residual,
            &manufactured_balance);
    require_collective(
        reference_duplicate_error ==
                PETSC_SUCCESS &&
            balance_duplicate_error ==
                PETSC_SUCCESS,
        "manufactured Newton reference vector allocation failed");

    for (const auto& entry :
         manufactured.residual_entries()) {
        const PetscInt index =
            entry.petsc_global_row;
        const PetscScalar value =
            static_cast<PetscScalar>(
                manufactured_reference_delta(
                    index));
        require_collective(
            VecSetValues(
                reference_delta,
                1,
                &index,
                &value,
                INSERT_VALUES) ==
                PETSC_SUCCESS,
            "manufactured Newton reference correction insertion failed");
    }
    const PetscErrorCode reference_begin =
        VecAssemblyBegin(
            reference_delta);
    const PetscErrorCode reference_end =
        VecAssemblyEnd(
            reference_delta);
    require_collective(
        reference_begin ==
                PETSC_SUCCESS &&
            reference_end ==
                PETSC_SUCCESS,
        "manufactured Newton reference correction assembly failed");

    const PetscErrorCode manufactured_matmult_error =
        MatMult(
            jacobian,
            reference_delta,
            manufactured_balance);
    const PetscErrorCode manufactured_axpy_error =
        VecAXPY(
            manufactured_balance,
            PetscScalar{1.0},
            residual);
    require_collective(
        manufactured_matmult_error ==
                PETSC_SUCCESS &&
            manufactured_axpy_error ==
                PETSC_SUCCESS,
        "manufactured Newton exact-balance evaluation failed");

    PetscReal manufactured_balance_norm = 0.0;
    require_collective(
        VecNorm(
            manufactured_balance,
            NORM_2,
            &manufactured_balance_norm) ==
                PETSC_SUCCESS,
        "manufactured Newton exact-balance norm failed");
    require_collective(
        static_cast<double>(
            manufactured_balance_norm) <=
            5.0e-12,
        "manufactured Newton fixture does not satisfy J*delta_ref + R = 0");

    Vec residual_before = nullptr;
    const PetscErrorCode residual_before_duplicate_error =
        VecDuplicate(
            residual,
            &residual_before);
    const PetscErrorCode residual_before_copy_error =
        VecCopy(
            residual,
            residual_before);
    require_collective(
        residual_before_duplicate_error ==
                PETSC_SUCCESS &&
            residual_before_copy_error ==
                PETSC_SUCCESS,
        "frozen Newton residual preservation copy failed");

    Vec correction = nullptr;
    std::optional<
        fdp::FrozenNaturalVariableNewtonLinearCorrection3D>
        report;
    error =
        fdp::
            solve_frozen_natural_variable_newton_linear_system_3d(
                PETSC_COMM_WORLD,
                manufactured,
                jacobian,
                residual,
                &correction,
                &report);
    require_collective(
        error == PETSC_SUCCESS &&
            correction != nullptr &&
            report.has_value(),
        "frozen Newton linear-system contract failed");

    require_collective(
        report->ksp_type() ==
                std::string_view{
                    KSPGMRES} &&
            report->pc_type() ==
                std::string_view{
                    PCNONE} &&
            static_cast<int>(
                report->converged_reason()) >
                0 &&
            report->iteration_count() > 0 &&
            report->iteration_count() <=
                manufactured
                    .petsc_scalar_row_count(),
        "frozen Newton baseline KSP configuration/convergence mismatch");

    near_collective(
        report->residual_l2_norm(),
        report->rhs_l2_norm(),
        2.0e-13,
        2.0e-13);
    require_collective(
        report->linear_residual_l2_norm() <=
            5.0e-10 *
                std::max(
                    1.0,
                    report->rhs_l2_norm()),
        "frozen Newton reported linear residual is too large");

    const PetscInt expected_start =
        manufactured
            .petsc_scalar_row_start();
    const auto& cell_layout =
        identity(stable_cell).layout;
    require_collective(
        report->locally_owned_correction().size() ==
                q &&
            cell_layout.unknown_count() ==
                q,
        "frozen Newton local correction/layout cardinality mismatch");

    for (std::size_t slot = 0U;
         slot < q;
         ++slot) {
        const auto& entry =
            report
                ->locally_owned_correction()[
                    slot];
        const PetscInt global_scalar =
            expected_start +
            static_cast<PetscInt>(
                slot);
        const std::uint64_t expected_mesh_global =
            stable_cell == UINT64_C(10)
                ? static_cast<std::uint64_t>(
                      q + slot)
                : static_cast<std::uint64_t>(
                      slot);

        require_collective(
            entry.petsc_global_scalar ==
                    global_scalar &&
                entry.mesh_global_dof ==
                    mesh::GlobalDofIndex{
                        expected_mesh_global} &&
                entry.cell_global ==
                    mesh::GlobalEntityId{
                        stable_cell} &&
                entry.natural_variable_slot ==
                    slot,
            "frozen Newton correction global/local numbering mismatch");

        const bool composition_slot =
            slot >= 4U;
        require_collective(
            cell_layout
                    .composition_unknown_identity(
                        slot)
                    .has_value() ==
                composition_slot,
            "frozen Newton correction slot does not match natural-variable layout semantics");

        const double expected_value =
            manufactured_reference_delta(
                global_scalar);
        near_collective(
            entry.value,
            expected_value,
            2.0e-10,
            2.0e-12);

        PetscScalar vector_value{};
        require_collective(
            VecGetValues(
                correction,
                1,
                &global_scalar,
                &vector_value) ==
                PETSC_SUCCESS,
            "frozen Newton correction Vec readback failed");
        near_collective(
            static_cast<double>(
                PetscRealPart(
                    vector_value)),
            expected_value,
            2.0e-10,
            2.0e-12);
    }

    Vec independent_linear_residual = nullptr;
    const PetscErrorCode independent_duplicate_error =
        VecDuplicate(
            residual,
            &independent_linear_residual);
    const PetscErrorCode independent_matmult_error =
        MatMult(
            jacobian,
            correction,
            independent_linear_residual);
    const PetscErrorCode independent_axpy_error =
        VecAXPY(
            independent_linear_residual,
            PetscScalar{1.0},
            residual);
    require_collective(
        independent_duplicate_error ==
                PETSC_SUCCESS &&
            independent_matmult_error ==
                PETSC_SUCCESS &&
            independent_axpy_error ==
                PETSC_SUCCESS,
        "independent frozen Newton linear residual construction failed");

    PetscReal independent_linear_norm = 0.0;
    require_collective(
        VecNorm(
            independent_linear_residual,
            NORM_2,
            &independent_linear_norm) ==
                PETSC_SUCCESS,
        "independent frozen Newton linear residual norm failed");
    near_collective(
        static_cast<double>(
            independent_linear_norm),
        report->linear_residual_l2_norm(),
        5.0e-10,
        5.0e-12);

    Vec residual_difference = nullptr;
    const PetscErrorCode residual_difference_duplicate_error =
        VecDuplicate(
            residual,
            &residual_difference);
    const PetscErrorCode residual_difference_copy_error =
        VecCopy(
            residual_before,
            residual_difference);
    const PetscErrorCode residual_difference_axpy_error =
        VecAXPY(
            residual_difference,
            PetscScalar{-1.0},
            residual);
    require_collective(
        residual_difference_duplicate_error ==
                PETSC_SUCCESS &&
            residual_difference_copy_error ==
                PETSC_SUCCESS &&
            residual_difference_axpy_error ==
                PETSC_SUCCESS,
        "frozen Newton residual preservation audit failed");

    PetscReal residual_difference_norm = 0.0;
    require_collective(
        VecNorm(
            residual_difference,
            NORM_INFINITY,
            &residual_difference_norm) ==
                PETSC_SUCCESS &&
            residual_difference_norm ==
                PetscReal{0.0},
        "frozen Newton solve modified the input residual Vec");

    // Wrong residual ownership/size must be rejected collectively before a
    // correction Vec is returned.
    Vec wrong_residual = nullptr;
    require_collective(
        VecCreateMPI(
            PETSC_COMM_WORLD,
            static_cast<PetscInt>(
                q + 1U),
            PETSC_DECIDE,
            &wrong_residual) ==
                PETSC_SUCCESS,
        "wrong-size residual fixture allocation failed");
    Vec wrong_correction = nullptr;
    std::optional<
        fdp::FrozenNaturalVariableNewtonLinearCorrection3D>
        wrong_report;
    const PetscErrorCode wrong_error =
        fdp::
            solve_frozen_natural_variable_newton_linear_system_3d(
                PETSC_COMM_WORLD,
                manufactured,
                jacobian,
                wrong_residual,
                &wrong_correction,
                &wrong_report);
    require_collective(
        wrong_error ==
                PETSC_ERR_ARG_INCOMP &&
            wrong_correction ==
                nullptr &&
            !wrong_report.has_value(),
        "wrong residual size/ownership was not rejected collectively");

    const PetscErrorCode wrong_destroy =
        VecDestroy(
            &wrong_residual);
    const PetscErrorCode residual_difference_destroy =
        VecDestroy(
            &residual_difference);
    const PetscErrorCode independent_destroy =
        VecDestroy(
            &independent_linear_residual);
    const PetscErrorCode correction_destroy =
        VecDestroy(
            &correction);
    const PetscErrorCode residual_before_destroy =
        VecDestroy(
            &residual_before);
    const PetscErrorCode manufactured_balance_destroy =
        VecDestroy(
            &manufactured_balance);
    const PetscErrorCode reference_destroy =
        VecDestroy(
            &reference_delta);
    const PetscErrorCode residual_destroy =
        VecDestroy(
            &residual);
    const PetscErrorCode jacobian_destroy =
        MatDestroy(
            &jacobian);

    require_collective(
        wrong_destroy ==
                PETSC_SUCCESS &&
            residual_difference_destroy ==
                PETSC_SUCCESS &&
            independent_destroy ==
                PETSC_SUCCESS &&
            correction_destroy ==
                PETSC_SUCCESS &&
            residual_before_destroy ==
                PETSC_SUCCESS &&
            manufactured_balance_destroy ==
                PETSC_SUCCESS &&
            reference_destroy ==
                PETSC_SUCCESS &&
            residual_destroy ==
                PETSC_SUCCESS &&
            jacobian_destroy ==
                PETSC_SUCCESS,
        "frozen Newton PETSc fixture cleanup failed");
}


fdp::CompleteNaturalVariableAssemblySnapshot3D
make_snes_numbering_snapshot(
    int rank) {
    constexpr std::size_t n = 3U;
    constexpr std::size_t q = 10U;
    const PetscInt start =
        static_cast<PetscInt>(
            rank) *
        static_cast<PetscInt>(
            q);
    const std::uint64_t stable_cell =
        rank == 0
            ? UINT64_C(10)
            : UINT64_C(20);

    std::vector<
        fdp::CompleteNaturalVariableResidualEntry3D>
        residuals;
    residuals.reserve(q);
    for (std::size_t slot = 0U;
         slot < q;
         ++slot) {
        const std::uint64_t mesh_global =
            stable_cell == UINT64_C(10)
                ? static_cast<std::uint64_t>(
                      q + slot)
                : static_cast<std::uint64_t>(
                      slot);
        const auto kind =
            slot < n
                ? fdp::NaturalVariableEquationKind3D::
                      component_conservation
                : (slot == n
                       ? fdp::NaturalVariableEquationKind3D::
                             energy_conservation
                       : fdp::NaturalVariableEquationKind3D::
                             fugacity_equilibrium);
        residuals.push_back(
            fdp::CompleteNaturalVariableResidualEntry3D{
                start +
                    static_cast<PetscInt>(
                        slot),
                mesh::GlobalDofIndex{
                    mesh_global},
                mesh::GlobalEntityId{
                    stable_cell},
                slot,
                kind,
                0.0});
    }

    std::vector<
        fdp::CompleteNaturalVariableJacobianEntry3D>
        jacobian;
    jacobian.reserve(14U * q);
    for (std::size_t row_slot = 0U;
         row_slot < q;
         ++row_slot) {
        const auto kind =
            row_slot < n
                ? fdp::NaturalVariableEquationKind3D::
                      component_conservation
                : (row_slot == n
                       ? fdp::NaturalVariableEquationKind3D::
                             energy_conservation
                       : fdp::NaturalVariableEquationKind3D::
                             fugacity_equilibrium);
        const PetscInt first_column =
            row_slot <= n
                ? PetscInt{0}
                : start;
        const PetscInt last_column =
            row_slot <= n
                ? PetscInt{20}
                : start +
                      static_cast<PetscInt>(
                          q);

        for (PetscInt global_column =
                 first_column;
             global_column < last_column;
             ++global_column) {
            const std::size_t column_slot =
                static_cast<std::size_t>(
                    global_column %
                    static_cast<PetscInt>(
                        q));
            const std::uint64_t column_cell =
                global_column <
                        static_cast<PetscInt>(
                            q)
                    ? UINT64_C(10)
                    : UINT64_C(20);
            const std::uint64_t
                mesh_global_column =
                    column_cell ==
                            UINT64_C(10)
                        ? static_cast<
                              std::uint64_t>(
                                  q +
                                  column_slot)
                        : static_cast<
                              std::uint64_t>(
                                  column_slot);
            jacobian.push_back(
                fdp::CompleteNaturalVariableJacobianEntry3D{
                    start +
                        static_cast<PetscInt>(
                            row_slot),
                    global_column,
                    residuals[row_slot]
                        .mesh_global_row_dof,
                    mesh::GlobalDofIndex{
                        mesh_global_column},
                    mesh::GlobalEntityId{
                        stable_cell},
                    mesh::GlobalEntityId{
                        column_cell},
                    row_slot,
                    column_slot,
                    kind,
                    global_column ==
                            start +
                                static_cast<
                                    PetscInt>(
                                    row_slot)
                        ? 1.0
                        : 0.0});
        }
    }

    return {
        mesh::PartitionRank{
            static_cast<
                mesh::PartitionRank::value_type>(
                    rank)},
        2U,
        "natural_state",
        n,
        q,
        start,
        start +
            static_cast<PetscInt>(
                q),
        PetscInt{20},
        std::move(residuals),
        std::move(jacobian)};
}

double snes_target_value(
    PetscInt global_scalar) {
    if (global_scalar < 0) {
        throw std::invalid_argument(
            "negative SNES manufactured scalar index");
    }
    return 1.0 +
        0.01 *
            static_cast<double>(
                global_scalar + 1);
}

struct SnesManufacturedContext {
    PetscInt function_calls{};
    PetscInt jacobian_calls{};
    PetscInt precheck_calls{};
    PetscInt precheck_changes{};
};

PetscErrorCode snes_manufactured_function(
    Vec state,
    Vec residual,
    void* raw_context,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (state == nullptr ||
        residual == nullptr ||
        raw_context == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* context =
        static_cast<
            SnesManufacturedContext*>(
                raw_context);
    ++context->function_calls;

    PetscInt start = -1;
    PetscInt end = -1;
    PetscErrorCode error =
        VecGetOwnershipRange(
            state,
            &start,
            &end);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    const PetscScalar* x = nullptr;
    error =
        VecGetArrayRead(
            state,
            &x);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    bool domain_error = false;
    for (PetscInt global = start;
         global < end;
         ++global) {
        const std::size_t local =
            static_cast<std::size_t>(
                global - start);
        const double value =
            static_cast<double>(
                PetscRealPart(
                    x[local]));
        if (!std::isfinite(value) ||
            value <= 0.0) {
            domain_error = true;
            break;
        }
    }

    const PetscErrorCode restore_error =
        VecRestoreArrayRead(
            state,
            &x);
    if (restore_error !=
        PETSC_SUCCESS) {
        return restore_error;
    }

    if (domain_error) {
        *status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    domain_error;
        return PETSC_SUCCESS;
    }

    PetscScalar* f = nullptr;
    error =
        VecGetArray(
            residual,
            &f);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    error =
        VecGetArrayRead(
            state,
            &x);
    if (error != PETSC_SUCCESS) {
        (void)VecRestoreArray(
            residual,
            &f);
        return error;
    }

    for (PetscInt global = start;
         global < end;
         ++global) {
        const std::size_t local =
            static_cast<std::size_t>(
                global - start);
        const double value =
            static_cast<double>(
                PetscRealPart(
                    x[local]));
        const double target =
            snes_target_value(
                global);
        f[local] =
            static_cast<PetscScalar>(
                std::log(
                    value /
                    target));
    }

    const PetscErrorCode x_restore_error =
        VecRestoreArrayRead(
            state,
            &x);
    const PetscErrorCode f_restore_error =
        VecRestoreArray(
            residual,
            &f);
    if (x_restore_error !=
        PETSC_SUCCESS) {
        return x_restore_error;
    }
    if (f_restore_error !=
        PETSC_SUCCESS) {
        return f_restore_error;
    }

    *status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    return PETSC_SUCCESS;
}

PetscErrorCode snes_manufactured_precheck(
    Vec state,
    Vec search_direction,
    void* raw_context,
    PetscBool* changed_direction) {
    if (state == nullptr ||
        search_direction == nullptr ||
        raw_context == nullptr ||
        changed_direction == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    auto* context =
        static_cast<
            SnesManufacturedContext*>(
                raw_context);
    ++context->precheck_calls;
    *changed_direction =
        PETSC_FALSE;

    PetscInt start = -1;
    PetscInt end = -1;
    PetscErrorCode error =
        VecGetOwnershipRange(
            state,
            &start,
            &end);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    const PetscScalar* x = nullptr;
    const PetscScalar* y = nullptr;
    error =
        VecGetArrayRead(
            state,
            &x);
    if (error == PETSC_SUCCESS) {
        error =
            VecGetArrayRead(
                search_direction,
                &y);
    }
    if (error != PETSC_SUCCESS) {
        if (x != nullptr) {
            (void)VecRestoreArrayRead(
                state,
                &x);
        }
        return error;
    }

    double local_scale = 1.0;
    for (PetscInt global = start;
         global < end;
         ++global) {
        const std::size_t local =
            static_cast<std::size_t>(
                global - start);
        const double state_value =
            static_cast<double>(
                PetscRealPart(
                    x[local]));
        const double direction_value =
            static_cast<double>(
                PetscRealPart(
                    y[local]));
        if (!std::isfinite(state_value) ||
            !std::isfinite(direction_value) ||
            state_value <= 0.0) {
            (void)VecRestoreArrayRead(
                search_direction,
                &y);
            (void)VecRestoreArrayRead(
                state,
                &x);
            return PETSC_ERR_ARG_OUTOFRANGE;
        }

        if (direction_value > 0.0 &&
            state_value -
                    direction_value <=
                0.0) {
            local_scale =
                std::min(
                    local_scale,
                    0.8 *
                        state_value /
                        direction_value);
        }
    }

    const PetscErrorCode y_restore =
        VecRestoreArrayRead(
            search_direction,
            &y);
    const PetscErrorCode x_restore =
        VecRestoreArrayRead(
            state,
            &x);
    if (y_restore != PETSC_SUCCESS) {
        return y_restore;
    }
    if (x_restore != PETSC_SUCCESS) {
        return x_restore;
    }

    double global_scale = 1.0;
    if (MPI_Allreduce(
            &local_scale,
            &global_scale,
            1,
            MPI_DOUBLE,
            MPI_MIN,
            PetscObjectComm(
                reinterpret_cast<PetscObject>(
                    state))) !=
        MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    if (!std::isfinite(global_scale) ||
        global_scale <= 0.0 ||
        global_scale > 1.0) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    if (global_scale < 1.0) {
        error =
            VecScale(
                search_direction,
                static_cast<PetscScalar>(
                    global_scale));
        if (error != PETSC_SUCCESS) {
            return error;
        }
        ++context->precheck_changes;
        *changed_direction =
            PETSC_TRUE;
    }

    return PETSC_SUCCESS;
}

PetscErrorCode snes_manufactured_jacobian(
    Vec state,
    Mat jacobian,
    void* raw_context,
    fdp::NaturalVariableSnesEvaluationStatus3D*
        status) {
    if (state == nullptr ||
        jacobian == nullptr ||
        raw_context == nullptr ||
        status == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    auto* context =
        static_cast<
            SnesManufacturedContext*>(
                raw_context);
    ++context->jacobian_calls;

    PetscInt start = -1;
    PetscInt end = -1;
    PetscErrorCode error =
        VecGetOwnershipRange(
            state,
            &start,
            &end);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    const PetscScalar* x = nullptr;
    error =
        VecGetArrayRead(
            state,
            &x);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    for (PetscInt global = start;
         global < end;
         ++global) {
        const std::size_t local =
            static_cast<std::size_t>(
                global - start);
        const double value =
            static_cast<double>(
                PetscRealPart(
                    x[local]));
        if (!std::isfinite(value) ||
            value <= 0.0) {
            const PetscErrorCode restore =
                VecRestoreArrayRead(
                    state,
                    &x);
            if (restore != PETSC_SUCCESS) {
                return restore;
            }
            *status =
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        domain_error;
            return PETSC_SUCCESS;
        }

        const PetscScalar derivative =
            static_cast<PetscScalar>(
                1.0 /
                value);
        error =
            MatSetValues(
                jacobian,
                1,
                &global,
                1,
                &global,
                &derivative,
                INSERT_VALUES);
        if (error != PETSC_SUCCESS) {
            (void)VecRestoreArrayRead(
                state,
                &x);
            return error;
        }
    }

    error =
        VecRestoreArrayRead(
            state,
            &x);
    if (error != PETSC_SUCCESS) {
        return error;
    }

    *status =
        fdp::
            NaturalVariableSnesEvaluationStatus3D::
                success;
    return PETSC_SUCCESS;
}

void natural_variable_snes_solver() {
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
        "SNES natural-variable test requires exactly two ranks");

    constexpr std::size_t q = 10U;
    auto numbering =
        make_snes_numbering_snapshot(
            rank);
    auto cell_bridge =
        make_cell_row_bridge(
            rank);

    Vec unused_residual = nullptr;
    Mat jacobian_template = nullptr;
    PetscErrorCode error =
        fdp::
            materialize_complete_natural_variable_petsc_system_3d(
                PETSC_COMM_WORLD,
                numbering,
                cell_bridge,
                &unused_residual,
                &jacobian_template);
    require_collective(
        error == PETSC_SUCCESS &&
            unused_residual != nullptr &&
            jacobian_template != nullptr,
        "SNES Jacobian structure template materialization failed");

    Vec initial_state = nullptr;
    const PetscErrorCode initial_duplicate_error =
        VecDuplicate(
            unused_residual,
            &initial_state);
    require_collective(
        initial_duplicate_error ==
                PETSC_SUCCESS,
        "SNES initial-state allocation failed");

    const PetscInt start =
        numbering.petsc_scalar_row_start();
    for (std::size_t slot = 0U;
         slot < q;
         ++slot) {
        const PetscInt global =
            start +
            static_cast<PetscInt>(
                slot);
        const PetscScalar value =
            static_cast<PetscScalar>(
                10.0 *
                snes_target_value(
                    global));
        require_collective(
            VecSetValues(
                initial_state,
                1,
                &global,
                &value,
                INSERT_VALUES) ==
                PETSC_SUCCESS,
            "SNES initial-state insertion failed");
    }
    const PetscErrorCode initial_begin =
        VecAssemblyBegin(
            initial_state);
    const PetscErrorCode initial_end =
        VecAssemblyEnd(
            initial_state);
    require_collective(
        initial_begin ==
                PETSC_SUCCESS &&
            initial_end ==
                PETSC_SUCCESS,
        "SNES initial-state assembly failed");

    Vec initial_before = nullptr;
    const PetscErrorCode initial_before_duplicate =
        VecDuplicate(
            initial_state,
            &initial_before);
    const PetscErrorCode initial_before_copy =
        VecCopy(
            initial_state,
            initial_before);
    require_collective(
        initial_before_duplicate ==
                PETSC_SUCCESS &&
            initial_before_copy ==
                PETSC_SUCCESS,
        "SNES initial-state preservation copy failed");

    SnesManufacturedContext manufactured_context;
    const fdp::NaturalVariableSnesEvaluator3D
        evaluator{
            snes_manufactured_function,
            snes_manufactured_jacobian,
            snes_manufactured_precheck,
            &manufactured_context};

    Vec solution = nullptr;
    std::optional<
        fdp::NaturalVariableSnesSolveReport3D>
        report;
    error =
        fdp::
            solve_natural_variable_snes_3d(
                PETSC_COMM_WORLD,
                numbering,
                initial_state,
                jacobian_template,
                evaluator,
                &solution,
                &report);
    require_collective(
        error == PETSC_SUCCESS &&
            solution != nullptr &&
            report.has_value(),
        "PETSc SNES natural-variable solve failed");

    require_collective(
        report->snes_type() ==
                std::string_view{
                    SNESNEWTONLS} &&
            report->line_search_type() ==
                std::string_view{
                    SNESLINESEARCHBT} &&
            report->ksp_type() ==
                std::string_view{
                    KSPGMRES} &&
            report->pc_type() ==
                std::string_view{
                    PCASM} &&
            report->default_asm_overlap() ==
                1 &&
            static_cast<int>(
                report->converged_reason()) >
                0 &&
            report->nonlinear_iterations() >
                0 &&
            report->function_evaluations() >=
                report->nonlinear_iterations() &&
            report->jacobian_evaluations() >
                0 &&
            report->function_domain_errors() ==
                0 &&
            report->jacobian_domain_errors() ==
                0 &&
            report->line_search_prechecks() >
                0 &&
            report->line_search_direction_changes() >
                0,
        "SNES/GMRES/ASM/BT configuration or domain handling mismatch");

    require_collective(
        manufactured_context.function_calls ==
                report->function_evaluations() &&
            manufactured_context.jacobian_calls ==
                report->jacobian_evaluations() &&
            manufactured_context.precheck_calls ==
                report->line_search_prechecks() &&
            manufactured_context.precheck_changes ==
                report
                    ->line_search_direction_changes(),
        "SNES adapter callback/precheck accounting mismatch");

    require_collective(
        report->final_function_l2_norm() <=
            2.0e-9,
        "SNES final manufactured function norm is too large");

    const std::uint64_t stable_cell =
        rank == 0
            ? UINT64_C(10)
            : UINT64_C(20);
    const auto& cell_layout =
        identity(stable_cell).layout;
    require_collective(
        report->locally_owned_solution().size() ==
                q &&
            cell_layout.unknown_count() ==
                q,
        "SNES solution/layout cardinality mismatch");

    for (std::size_t slot = 0U;
         slot < q;
         ++slot) {
        const PetscInt global =
            start +
            static_cast<PetscInt>(
                slot);
        const double target =
            snes_target_value(
                global);
        PetscScalar actual{};
        require_collective(
            VecGetValues(
                solution,
                1,
                &global,
                &actual) ==
                PETSC_SUCCESS,
            "SNES solution Vec readback failed");
        near_collective(
            static_cast<double>(
                PetscRealPart(
                    actual)),
            target,
            5.0e-9,
            5.0e-11);

        const auto& entry =
            report
                ->locally_owned_solution()[
                    slot];
        const std::uint64_t
            expected_mesh_global =
                stable_cell ==
                        UINT64_C(10)
                    ? static_cast<
                          std::uint64_t>(
                              q +
                              slot)
                    : static_cast<
                          std::uint64_t>(
                              slot);
        require_collective(
            entry.petsc_global_scalar ==
                    global &&
                entry.mesh_global_dof ==
                    mesh::GlobalDofIndex{
                        expected_mesh_global} &&
                entry.cell_global ==
                    mesh::GlobalEntityId{
                        stable_cell} &&
                entry.natural_variable_slot ==
                    slot,
            "SNES solution global/local natural-variable numbering mismatch");
        near_collective(
            entry.value,
            target,
            5.0e-9,
            5.0e-11);
    }

    Vec independent_residual = nullptr;
    require_collective(
        VecDuplicate(
            solution,
            &independent_residual) ==
                PETSC_SUCCESS,
        "SNES independent residual allocation failed");
    fdp::NaturalVariableSnesEvaluationStatus3D
        final_status =
            fdp::
                NaturalVariableSnesEvaluationStatus3D::
                    success;
    error =
        snes_manufactured_function(
            solution,
            independent_residual,
            &manufactured_context,
            &final_status);
    require_collective(
        error == PETSC_SUCCESS &&
            final_status ==
                fdp::
                    NaturalVariableSnesEvaluationStatus3D::
                        success,
        "SNES independent final residual evaluation failed");
    const PetscErrorCode final_begin =
        VecAssemblyBegin(
            independent_residual);
    const PetscErrorCode final_end =
        VecAssemblyEnd(
            independent_residual);
    require_collective(
        final_begin ==
                PETSC_SUCCESS &&
            final_end ==
                PETSC_SUCCESS,
        "SNES independent final residual assembly failed");

    PetscReal independent_norm = 0.0;
    require_collective(
        VecNorm(
            independent_residual,
            NORM_2,
            &independent_norm) ==
                PETSC_SUCCESS,
        "SNES independent final residual norm failed");
    near_collective(
        static_cast<double>(
            independent_norm),
        report->final_function_l2_norm(),
        5.0e-8,
        5.0e-11);

    Vec initial_difference = nullptr;
    const PetscErrorCode initial_difference_duplicate =
        VecDuplicate(
            initial_state,
            &initial_difference);
    const PetscErrorCode initial_difference_copy =
        VecCopy(
            initial_before,
            initial_difference);
    const PetscErrorCode initial_difference_axpy =
        VecAXPY(
            initial_difference,
            PetscScalar{-1.0},
            initial_state);
    require_collective(
        initial_difference_duplicate ==
                PETSC_SUCCESS &&
            initial_difference_copy ==
                PETSC_SUCCESS &&
            initial_difference_axpy ==
                PETSC_SUCCESS,
        "SNES initial-state preservation audit failed");

    PetscReal initial_difference_norm = 0.0;
    require_collective(
        VecNorm(
            initial_difference,
            NORM_INFINITY,
            &initial_difference_norm) ==
                PETSC_SUCCESS &&
            initial_difference_norm ==
                PetscReal{0.0},
        "SNES adapter modified caller initial-state Vec");

    Vec immediate_state = nullptr;
    require_collective(
        VecDuplicate(
            unused_residual,
            &immediate_state) ==
                PETSC_SUCCESS,
        "SNES immediate-convergence state allocation failed");
    for (std::size_t slot = 0U;
         slot < q;
         ++slot) {
        const PetscInt global =
            start +
            static_cast<PetscInt>(
                slot);
        const PetscScalar value =
            static_cast<PetscScalar>(
                snes_target_value(
                    global));
        require_collective(
            VecSetValues(
                immediate_state,
                1,
                &global,
                &value,
                INSERT_VALUES) ==
                PETSC_SUCCESS,
            "SNES immediate-convergence state insertion failed");
    }
    require_collective(
        VecAssemblyBegin(
            immediate_state) ==
                PETSC_SUCCESS &&
            VecAssemblyEnd(
                immediate_state) ==
                PETSC_SUCCESS,
        "SNES immediate-convergence state assembly failed");

    SnesManufacturedContext
        immediate_context;
    const fdp::NaturalVariableSnesEvaluator3D
        immediate_evaluator{
            snes_manufactured_function,
            snes_manufactured_jacobian,
            snes_manufactured_precheck,
            &immediate_context};
    Vec immediate_solution = nullptr;
    std::optional<
        fdp::NaturalVariableSnesSolveReport3D>
        immediate_report;
    error =
        fdp::
            solve_natural_variable_snes_3d(
                PETSC_COMM_WORLD,
                numbering,
                immediate_state,
                jacobian_template,
                immediate_evaluator,
                &immediate_solution,
                &immediate_report);
    require_collective(
        error == PETSC_SUCCESS &&
            immediate_solution != nullptr &&
            immediate_report.has_value() &&
            static_cast<int>(
                immediate_report
                    ->converged_reason()) >
                0 &&
            immediate_report
                    ->nonlinear_iterations() ==
                0 &&
            immediate_report
                    ->jacobian_evaluations() ==
                0 &&
            immediate_report
                    ->line_search_prechecks() ==
                0 &&
            immediate_context
                    .jacobian_calls ==
                0 &&
            immediate_context
                    .precheck_calls ==
                0 &&
            immediate_report
                    ->function_evaluations() >=
                1 &&
            immediate_report
                    ->final_function_l2_norm() <=
                1.0e-14,
        "SNES immediate convergence must publish a valid zero-iteration/zero-Jacobian report");
    const PetscErrorCode
        immediate_solution_destroy =
            VecDestroy(
                &immediate_solution);
    const PetscErrorCode
        immediate_state_destroy =
            VecDestroy(
                &immediate_state);

    const PetscErrorCode initial_difference_destroy =
        VecDestroy(
            &initial_difference);
    const PetscErrorCode independent_residual_destroy =
        VecDestroy(
            &independent_residual);
    const PetscErrorCode solution_destroy =
        VecDestroy(
            &solution);
    const PetscErrorCode initial_before_destroy =
        VecDestroy(
            &initial_before);
    const PetscErrorCode initial_destroy =
        VecDestroy(
            &initial_state);
    const PetscErrorCode unused_residual_destroy =
        VecDestroy(
            &unused_residual);
    const PetscErrorCode jacobian_template_destroy =
        MatDestroy(
            &jacobian_template);
    require_collective(
        initial_difference_destroy ==
                PETSC_SUCCESS &&
            independent_residual_destroy ==
                PETSC_SUCCESS &&
            solution_destroy ==
                PETSC_SUCCESS &&
            immediate_solution_destroy ==
                PETSC_SUCCESS &&
            immediate_state_destroy ==
                PETSC_SUCCESS &&
            initial_before_destroy ==
                PETSC_SUCCESS &&
            initial_destroy ==
                PETSC_SUCCESS &&
            unused_residual_destroy ==
                PETSC_SUCCESS &&
            jacobian_template_destroy ==
                PETSC_SUCCESS,
        "SNES manufactured PETSc fixture cleanup failed");
}

void invalid_collective_energy_inputs() {
    int rank = -1;
    if (MPI_Comm_rank(
            PETSC_COMM_WORLD,
            &rank) != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI rank query failed");
    }

    {
        auto partition =
            make_partition(rank);
        auto schedule =
            make_schedule(rank);
        EnergyLocalInputs inputs{rank};

        if (rank == 1) {
            inputs.cells[1]
                .state_identity
                .temperature_k +=
                0.25;
        }

        PetscErrorCode error =
            PETSC_SUCCESS;
        auto output =
            build_energy(
                partition,
                schedule,
                inputs,
                &error);
        require_collective(
            error == PETSC_ERR_ARG_INCOMP &&
                !output.has_value(),
            "stale ghost energy state was not rejected collectively");
    }

    {
        auto partition =
            make_partition(rank);
        auto schedule =
            make_schedule(rank);
        EnergyLocalInputs inputs{rank};

        if (rank == 1) {
            inputs.cells[0]
                .bulk_volume_m3 =
                6.0;
        }

        PetscErrorCode error =
            PETSC_SUCCESS;
        auto output =
            build_energy(
                partition,
                schedule,
                inputs,
                &error);
        require_collective(
            error == PETSC_ERR_ARG_INCOMP &&
                !output.has_value(),
            "energy endpoint volume mismatch was not rejected collectively");
    }

    {
        auto partition =
            make_partition(rank);
        auto schedule =
            make_schedule(rank);
        EnergyLocalInputs inputs{rank};

        if (rank == 0) {
            inputs.faces.clear();
        }

        PetscErrorCode error =
            PETSC_SUCCESS;
        auto output =
            build_energy(
                partition,
                schedule,
                inputs,
                &error);
        require_collective(
            error == PETSC_ERR_ARG_INCOMP &&
                !output.has_value(),
            "missing authoritative energy face binding was not rejected collectively");
    }

    {
        auto partition =
            make_partition(rank);
        auto schedule =
            make_schedule(rank);
        EnergyLocalInputs inputs{rank};

        if (rank == 0) {
            inputs.accumulation10.porosity =
                1.25;
        }

        PetscErrorCode error =
            PETSC_SUCCESS;
        auto output =
            build_energy(
                partition,
                schedule,
                inputs,
                &error);
        require_collective(
            error == PETSC_ERR_ARG_INCOMP &&
                !output.has_value(),
            "invalid owned energy accumulation was not rejected collectively");
    }
}

void invalid_collective_inputs() {
    int rank = -1;
    if (MPI_Comm_rank(
            PETSC_COMM_WORLD,
            &rank) != MPI_SUCCESS) {
        throw std::runtime_error(
            "MPI rank query failed");
    }

    {
        auto partition =
            make_partition(rank);
        auto schedule =
            make_schedule(rank);
        LocalInputs inputs{rank};

        if (rank == 1) {
            inputs.cells[1]
                .state_identity
                .temperature_k +=
                0.5;
        }

        PetscErrorCode error =
            PETSC_SUCCESS;
        auto output =
            build(
                rank,
                partition,
                schedule,
                inputs,
                &error);
        require_collective(
            error == PETSC_ERR_ARG_INCOMP &&
                !output.has_value(),
            "state fingerprint mismatch was not rejected collectively");
    }

    {
        auto partition =
            make_partition(rank);
        auto schedule =
            make_schedule(rank);
        LocalInputs inputs{rank};

        if (rank == 1) {
            inputs.cells[0]
                .bulk_volume_m3 =
                6.0;
        }

        PetscErrorCode error =
            PETSC_SUCCESS;
        auto output =
            build(
                rank,
                partition,
                schedule,
                inputs,
                &error);
        require_collective(
            error == PETSC_ERR_ARG_INCOMP &&
                !output.has_value(),
            "endpoint volume mismatch was not rejected collectively");
    }

    {
        auto partition =
            make_partition(rank);
        auto schedule =
            make_schedule(rank);
        LocalInputs inputs{rank};

        if (rank == 0) {
            inputs.faces.clear();
        }

        PetscErrorCode error =
            PETSC_SUCCESS;
        auto output =
            build(
                rank,
                partition,
                schedule,
                inputs,
                &error);
        require_collective(
            error == PETSC_ERR_ARG_INCOMP &&
                !output.has_value(),
            "missing authoritative face binding was not rejected collectively");
    }
}

void headers() {
    require_collective(
        distributed_component_conservation_header(),
        "distributed component conservation header probe failed");
    require_collective(
        distributed_energy_conservation_header(),
        "distributed energy conservation header probe failed");
    require_collective(
        energy_global_assembly_mapping_header(),
        "energy global assembly mapping header probe failed");
    require_collective(
        complete_natural_variable_assembly_snapshot_header(),
        "complete natural-variable assembly snapshot header probe failed");
    require_collective(
        complete_natural_variable_petsc_materialization_header(),
        "complete natural-variable PETSc materialization header probe failed");
    require_collective(
        natural_variable_newton_linear_system_header(),
        "natural-variable Newton linear-system header probe failed");
    require_collective(
        natural_variable_snes_solver_header(),
        "natural-variable SNES solver header probe failed");
    require_collective(
        fixed_three_phase_snes_assembly_header(),
        "fixed-three-phase SNES assembly header probe failed");
    require_collective(
        single_phase_snes_assembly_header(),
        "single-phase SNES assembly header probe failed");
    require_collective(
        two_phase_snes_assembly_header(),
        "two-phase SNES assembly header probe failed");
    require_collective(
        variable_cardinality_natural_variable_numbering_header(),
        "variable-cardinality natural-variable numbering header probe failed");
    require_collective(
        variable_cardinality_snes_assembly_header(),
        "variable-cardinality SNES assembly header probe failed");
    require_collective(
        mixed_cardinality_physical_snes_assembly_header(),
        "mixed-cardinality physical SNES assembly header probe failed");
    require_collective(
        fixed_bhp_well_source_evaluator_header(),
        "fixed-BHP well source evaluator header probe failed");
    require_collective(
        pr76_production_cell_evaluator_header(),
        "PR76 production cell evaluator header probe failed");
    require_collective(
        pr76_single_phase_transition_target_rebuild_header(),
        "PR76 single-phase transition target rebuild header probe failed");
    require_collective(
        pr76_transition_rebuild_materialization_header(),
        "PR76 transition rebuild materialization header probe failed");
    require_collective(
        adaptive_timestep_controller_header(),
        "adaptive timestep controller header probe failed");
    require_collective(
        physical_timestep_driver_header(),
        "physical timestep driver header probe failed");
    require_collective(
        physical_time_loop_header(),
        "physical time loop header probe failed");
    require_collective(
        single_phase_adaptive_timestep_attempt_header(),
        "single-phase adaptive timestep attempt header probe failed");
    require_collective(
        single_phase_handoff_initial_system_header(),
        "single-phase handoff initial-system header probe failed");
    require_collective(
        cross_cardinality_tpfa_bridge_header(),
        "cross-cardinality TPFA bridge header probe failed");
    require_collective(
        thermodynamic_cross_cardinality_tpfa_adapter_header(),
        "thermodynamic cross-cardinality TPFA adapter header probe failed");
    require_collective(
        phase_transition_outer_rebuild_header(),
        "phase-transition outer rebuild header probe failed");
    require_collective(
        post_snes_phase_transition_controller_header(),
        "post-SNES phase-transition controller header probe failed");
    require_collective(
        post_snes_phase_transition_handoff_scanner_header(),
        "post-SNES phase-transition handoff scanner header probe failed");
    require_collective(
        post_snes_pt_flash_phase_transition_scanner_header(),
        "post-SNES PT flash phase-transition scanner header probe failed");
    require_collective(
        global_component_assembly_mapping_header(),
        "global component assembly mapping header probe failed");
    require_collective(
        fugacity_equilibrium_global_assembly_mapping_header(),
        "fugacity equilibrium global assembly mapping header probe failed");
}

} // namespace

int main(int argc, char** argv) {
    PetscErrorCode error =
        PetscInitialize(
            &argc,
            &argv,
            nullptr,
            nullptr);
    if (error != PETSC_SUCCESS) {
        return static_cast<int>(error);
    }

    int result = 0;
    try {
        int size = 0;
        if (MPI_Comm_size(
                PETSC_COMM_WORLD,
                &size) != MPI_SUCCESS ||
            size != 2) {
            throw std::runtime_error(
                "test must run with two MPI ranks");
        }

        distributed_exchange();
        distributed_energy_exchange();
        energy_global_assembly_mapping();
        complete_natural_variable_assembly_snapshot();
        complete_natural_variable_petsc_materialization();
        natural_variable_newton_linear_system();
        natural_variable_snes_solver();
        fixed_three_phase_snes_physical_assembly_test();
        single_phase_snes_physical_assembly_test();
        two_phase_snes_physical_assembly_test();
        variable_cardinality_natural_variable_numbering_test();
        variable_cardinality_snes_assembly_test();
        mixed_cardinality_physical_snes_assembly_test();
        pr76_production_fully_implicit_transient_test();
        pr76_li_firoozabadi_sour_gas_short_step_test();
        adaptive_timestep_controller_test();
        post_snes_pt_flash_phase_transition_scanner_test();
        global_component_assembly_mapping();
        fugacity_global_assembly_mapping();
        invalid_collective_inputs();
        invalid_collective_energy_inputs();
        headers();

        int rank = -1;
        if (MPI_Comm_rank(
                PETSC_COMM_WORLD,
                &rank) == MPI_SUCCESS &&
            rank == 0) {
            std::cout
                << "[PASS] PETSc SNES post-convergence phase-transition controller lifecycle\n";
        }
    } catch (const std::exception& exception) {
        int rank = -1;
        (void)MPI_Comm_rank(
            PETSC_COMM_WORLD,
            &rank);
        std::cerr
            << "[rank "
            << rank
            << "] [FAIL] "
            << exception.what()
            << '\n';
        result = 1;
    }

    error =
        PetscFinalize();
    if (error != PETSC_SUCCESS) {
        return static_cast<int>(error);
    }
    return result;
}
