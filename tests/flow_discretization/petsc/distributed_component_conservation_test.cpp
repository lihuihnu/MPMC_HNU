#include <mpmc/flow_discretization_petsc/distributed_component_conservation.hpp>

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

mesh::PartitionSnapshot
make_partition(int rank) {
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

    mesh::Topology topology{
        std::move(ids),
        {}};

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
        invalid_collective_inputs();
        headers();

        int rank = -1;
        if (MPI_Comm_rank(
                PETSC_COMM_WORLD,
                &rank) == MPI_SUCCESS &&
            rank == 0) {
            std::cout
                << "[PASS] distributed owner-targeted component conservation\n";
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
