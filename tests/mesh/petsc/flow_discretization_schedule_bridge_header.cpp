#include <mpmc/discretization_petsc/adapter.hpp>
#include <mpmc/flow_discretization/owned_multi_cell_component_conservation.hpp>

#include <vector>

bool flow_discretization_schedule_bridge_header() {
    mpmc::discretization_petsc::
        ParallelOwnedConnectionSchedule3D schedule{
            mpmc::mesh::PartitionRank{0U},
            1U,
            {},
            {}};

    const auto rows =
        mpmc::flow_discretization::
            copy_authoritative_internal_connection_rows(
                schedule);
    return rows.empty();
}
