#ifndef MPMC_MESH_PETSC_ADAPTER_HPP
#define MPMC_MESH_PETSC_ADAPTER_HPP

#include <mpmc/mesh/dof_layout.hpp>
#include <mpmc/mesh/dof_numbering.hpp>
#include <mpmc/mesh/partition_snapshot.hpp>
#include <mpmc/mesh/shared_entity_plan.hpp>

#include <petscsection.h>
#include <petscsf.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace mpmc::mesh_petsc {

namespace detail {

inline PetscErrorCode checked_petsc_int_size(std::size_t value, PetscInt* output) {
    if (output == nullptr) return PETSC_ERR_ARG_NULL;
    if (value > static_cast<std::size_t>(std::numeric_limits<PetscInt>::max())) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    *output = static_cast<PetscInt>(value);
    return PETSC_SUCCESS;
}

inline PetscErrorCode checked_petsc_int_u64(std::uint64_t value, PetscInt* output) {
    if (output == nullptr) return PETSC_ERR_ARG_NULL;
    const auto maximum =
        static_cast<std::uint64_t>(std::numeric_limits<PetscInt>::max());
    if (value > maximum) return PETSC_ERR_ARG_OUTOFRANGE;
    *output = static_cast<PetscInt>(value);
    return PETSC_SUCCESS;
}

inline PetscErrorCode checked_mpi_rank(
    mpmc::mesh::PartitionRank rank,
    PetscMPIInt* output) {
    if (output == nullptr) return PETSC_ERR_ARG_NULL;
    if (rank.value() >
        static_cast<std::uint32_t>(
            std::numeric_limits<PetscMPIInt>::max())) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    *output = static_cast<PetscMPIInt>(rank.value());
    return PETSC_SUCCESS;
}

inline PetscErrorCode validate_communicator(
    MPI_Comm comm,
    mpmc::mesh::PartitionRank local_rank,
    std::uint32_t rank_count) {
    int mpi_rank = -1;
    int mpi_size = -1;
    if (MPI_Comm_rank(comm, &mpi_rank) != MPI_SUCCESS) return PETSC_ERR_MPI;
    if (MPI_Comm_size(comm, &mpi_size) != MPI_SUCCESS) return PETSC_ERR_MPI;
    if (mpi_rank < 0 || mpi_size < 0) return PETSC_ERR_MPI;
    if (static_cast<std::uint32_t>(mpi_rank) != local_rank.value() ||
        static_cast<std::uint32_t>(mpi_size) != rank_count) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }
    return PETSC_SUCCESS;
}

struct PointRanges {
    std::size_t cell_begin{0U};
    std::size_t face_begin{0U};
    std::size_t vertex_begin{0U};
    std::size_t end{0U};
};

inline PetscErrorCode point_ranges(
    const mpmc::mesh::DofLayout& layout,
    PointRanges* ranges) {
    if (ranges == nullptr) return PETSC_ERR_ARG_NULL;

    const std::size_t cells = layout.entity_count(mpmc::mesh::EntityKind::cell);
    const std::size_t faces = layout.entity_count(mpmc::mesh::EntityKind::face);
    const std::size_t vertices =
        layout.entity_count(mpmc::mesh::EntityKind::vertex);

    if (faces > std::numeric_limits<std::size_t>::max() - cells) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    const std::size_t face_begin = cells;
    const std::size_t vertex_begin = cells + faces;
    if (vertices >
        std::numeric_limits<std::size_t>::max() - vertex_begin) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    ranges->cell_begin = 0U;
    ranges->face_begin = face_begin;
    ranges->vertex_begin = vertex_begin;
    ranges->end = vertex_begin + vertices;
    return PETSC_SUCCESS;
}

inline std::size_t point_base(
    const PointRanges& ranges,
    mpmc::mesh::EntityKind kind) {
    switch (kind) {
    case mpmc::mesh::EntityKind::cell:
        return ranges.cell_begin;
    case mpmc::mesh::EntityKind::face:
        return ranges.face_begin;
    case mpmc::mesh::EntityKind::vertex:
        return ranges.vertex_begin;
    case mpmc::mesh::EntityKind::edge:
        break;
    }
    return std::numeric_limits<std::size_t>::max();
}

} // namespace detail

inline PetscErrorCode create_section_mapping(
    MPI_Comm comm,
    const mpmc::mesh::DofLayout& layout,
    const mpmc::mesh::DofNumberingSnapshot& numbering,
    PetscSection* section,
    std::vector<PetscInt>* local_to_global) {
    if (section == nullptr || local_to_global == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    *section = nullptr;
    local_to_global->clear();

    if (numbering.local_dof_count() != layout.total_dof_count()) {
        return PETSC_ERR_ARG_SIZ;
    }

    PetscErrorCode error = detail::validate_communicator(
        comm, numbering.local_rank(), numbering.rank_count());
    if (error != PETSC_SUCCESS) return error;

    local_to_global->reserve(numbering.local_dof_count());
    for (std::size_t local = 0U;
         local < numbering.local_dof_count();
         ++local) {
        PetscInt global = 0;
        error = detail::checked_petsc_int_u64(
            numbering.global_index(local).value(), &global);
        if (error != PETSC_SUCCESS) {
            local_to_global->clear();
            return error;
        }
        local_to_global->push_back(global);
    }

    detail::PointRanges ranges;
    error = detail::point_ranges(layout, &ranges);
    if (error != PETSC_SUCCESS) return error;

    PetscInt point_end = 0;
    PetscInt field_count = 0;
    error = detail::checked_petsc_int_size(ranges.end, &point_end);
    if (error != PETSC_SUCCESS) return error;
    error = detail::checked_petsc_int_size(layout.variable_count(), &field_count);
    if (error != PETSC_SUCCESS) return error;

    PetscSection local_section = nullptr;
    error = PetscSectionCreate(comm, &local_section);
    if (error != PETSC_SUCCESS) return error;

    error = PetscSectionSetNumFields(local_section, field_count);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&local_section);
        return error;
    }

    for (std::size_t variable_index = 0U;
         variable_index < layout.variable_count();
         ++variable_index) {
        PetscInt field = 0;
        PetscInt components = 0;
        error = detail::checked_petsc_int_size(variable_index, &field);
        if (error != PETSC_SUCCESS) {
            PetscSectionDestroy(&local_section);
            return error;
        }
        error = detail::checked_petsc_int_size(
            layout.variable(variable_index).component_count, &components);
        if (error != PETSC_SUCCESS) {
            PetscSectionDestroy(&local_section);
            return error;
        }
        error = PetscSectionSetFieldName(
            local_section, field,
            layout.variable(variable_index).id.c_str());
        if (error != PETSC_SUCCESS) {
            PetscSectionDestroy(&local_section);
            return error;
        }
        error = PetscSectionSetFieldComponents(
            local_section, field, components);
        if (error != PETSC_SUCCESS) {
            PetscSectionDestroy(&local_section);
            return error;
        }
    }

    error = PetscSectionSetPointMajor(local_section, PETSC_TRUE);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&local_section);
        return error;
    }

    error = PetscSectionSetChart(local_section, 0, point_end);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&local_section);
        return error;
    }

    for (const auto kind :
         {mpmc::mesh::EntityKind::cell,
          mpmc::mesh::EntityKind::face,
          mpmc::mesh::EntityKind::vertex}) {
        const std::size_t base = detail::point_base(ranges, kind);
        const std::size_t entity_count = layout.entity_count(kind);

        PetscInt total_dof = 0;
        error = detail::checked_petsc_int_size(
            layout.dofs_per_entity(kind), &total_dof);
        if (error != PETSC_SUCCESS) {
            PetscSectionDestroy(&local_section);
            return error;
        }

        for (std::size_t entity = 0U;
             entity < entity_count;
             ++entity) {
            PetscInt point = 0;
            error = detail::checked_petsc_int_size(base + entity, &point);
            if (error != PETSC_SUCCESS) {
                PetscSectionDestroy(&local_section);
                return error;
            }

            error = PetscSectionSetDof(local_section, point, total_dof);
            if (error != PETSC_SUCCESS) {
                PetscSectionDestroy(&local_section);
                return error;
            }

            for (std::size_t variable_index = 0U;
                 variable_index < layout.variable_count();
                 ++variable_index) {
                const auto& variable = layout.variable(variable_index);
                if (variable.location != kind) continue;

                PetscInt field = 0;
                PetscInt field_dof = 0;
                error = detail::checked_petsc_int_size(
                    variable_index, &field);
                if (error != PETSC_SUCCESS) {
                    PetscSectionDestroy(&local_section);
                    return error;
                }
                error = detail::checked_petsc_int_size(
                    variable.component_count, &field_dof);
                if (error != PETSC_SUCCESS) {
                    PetscSectionDestroy(&local_section);
                    return error;
                }
                error = PetscSectionSetFieldDof(
                    local_section, point, field, field_dof);
                if (error != PETSC_SUCCESS) {
                    PetscSectionDestroy(&local_section);
                    return error;
                }
            }
        }
    }

    error = PetscSectionSetUp(local_section);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&local_section);
        return error;
    }

    PetscInt storage_size = 0;
    error = PetscSectionGetStorageSize(local_section, &storage_size);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&local_section);
        return error;
    }
    PetscInt expected_storage = 0;
    error = detail::checked_petsc_int_size(
        layout.total_dof_count(), &expected_storage);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&local_section);
        return error;
    }
    if (storage_size != expected_storage) {
        PetscSectionDestroy(&local_section);
        return PETSC_ERR_PLIB;
    }

    *section = local_section;
    return PETSC_SUCCESS;
}

inline PetscErrorCode create_entity_sf(
    MPI_Comm comm,
    const mpmc::mesh::PartitionSnapshot& partition,
    const mpmc::mesh::SharedEntityPlan& plan,
    mpmc::mesh::EntityKind kind,
    PetscSF* sf) {
    if (sf == nullptr) return PETSC_ERR_ARG_NULL;
    *sf = nullptr;

    if (plan.local_rank() != partition.local_rank() ||
        plan.rank_count() != partition.rank_count()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    PetscErrorCode error = detail::validate_communicator(
        comm, partition.local_rank(), partition.rank_count());
    if (error != PETSC_SUCCESS) return error;

    PetscInt nroots = 0;
    error = detail::checked_petsc_int_size(
        partition.entity_count(kind), &nroots);
    if (error != PETSC_SUCCESS) return error;

    std::vector<PetscInt> local_leaves;
    std::vector<PetscSFNode> remote_roots;
    local_leaves.reserve(plan.receive_count());
    remote_roots.reserve(plan.receive_count());

    for (const auto& neighbor : plan.neighbors()) {
        PetscMPIInt remote_rank = 0;
        error = detail::checked_mpi_rank(neighbor.rank, &remote_rank);
        if (error != PETSC_SUCCESS) return error;

        for (const auto& entity :
             plan.receive_entities_from(neighbor.rank)) {
            if (entity.kind != kind) continue;

            PetscInt local_leaf = 0;
            PetscInt remote_root = 0;
            error = detail::checked_petsc_int_size(
                static_cast<std::size_t>(entity.local.value()),
                &local_leaf);
            if (error != PETSC_SUCCESS) return error;
            error = detail::checked_petsc_int_size(
                static_cast<std::size_t>(entity.remote_local.value()),
                &remote_root);
            if (error != PETSC_SUCCESS) return error;

            if (!partition.is_ghost(kind, entity.local) ||
                partition.owner_rank(kind, entity.local) != neighbor.rank ||
                partition.global_id(kind, entity.local) != entity.global_id) {
                return PETSC_ERR_ARG_INCOMP;
            }

            local_leaves.push_back(local_leaf);
            remote_roots.push_back(PetscSFNode{remote_rank, remote_root});
        }
    }

    PetscInt nleaves = 0;
    error = detail::checked_petsc_int_size(
        local_leaves.size(), &nleaves);
    if (error != PETSC_SUCCESS) return error;

    PetscSF local_sf = nullptr;
    error = PetscSFCreate(comm, &local_sf);
    if (error != PETSC_SUCCESS) return error;

    error = PetscSFSetGraph(
        local_sf,
        nroots,
        nleaves,
        local_leaves.empty() ? nullptr : local_leaves.data(),
        PETSC_COPY_VALUES,
        remote_roots.empty() ? nullptr : remote_roots.data(),
        PETSC_COPY_VALUES);
    if (error != PETSC_SUCCESS) {
        PetscSFDestroy(&local_sf);
        return error;
    }

    error = PetscSFSetUp(local_sf);
    if (error != PETSC_SUCCESS) {
        PetscSFDestroy(&local_sf);
        return error;
    }

    *sf = local_sf;
    return PETSC_SUCCESS;
}

} // namespace mpmc::mesh_petsc

#endif // MPMC_MESH_PETSC_ADAPTER_HPP
