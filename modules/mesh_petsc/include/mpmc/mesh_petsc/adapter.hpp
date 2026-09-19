#ifndef MPMC_MESH_PETSC_ADAPTER_HPP
#define MPMC_MESH_PETSC_ADAPTER_HPP

#include <mpmc/mesh/dense_field.hpp>
#include <mpmc/mesh/dof_layout.hpp>
#include <mpmc/mesh/dof_numbering.hpp>
#include <mpmc/mesh/face_boundary.hpp>
#include <mpmc/mesh/geometry_2d.hpp>
#include <mpmc/mesh/partition_snapshot.hpp>
#include <mpmc/mesh/shared_entity_plan.hpp>
#include <mpmc/mesh/topology.hpp>

#include <petscdmplex.h>
#include <petscsection.h>
#include <petscsf.h>
#include <petscvec.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
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

/// Create one flattened point PetscSF matching create_section_mapping()'s chart.
///
/// The local point space is [cell points][face points][vertex points]. SharedEntityPlan
/// stores remote indices within each EntityKind, so the adapter collectively gathers
/// each rank's local point counts to recover the remote chart base without changing
/// the core halo contract.
inline PetscErrorCode create_point_sf(
    MPI_Comm comm,
    const mpmc::mesh::DofLayout& layout,
    const mpmc::mesh::PartitionSnapshot& partition,
    const mpmc::mesh::SharedEntityPlan& plan,
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

    for (const auto kind :
         {mpmc::mesh::EntityKind::cell,
          mpmc::mesh::EntityKind::face,
          mpmc::mesh::EntityKind::vertex}) {
        if (layout.entity_count(kind) != partition.entity_count(kind)) {
            return PETSC_ERR_ARG_SIZ;
        }
    }

    detail::PointRanges local_ranges;
    error = detail::point_ranges(layout, &local_ranges);
    if (error != PETSC_SUCCESS) return error;

    const std::array<std::uint64_t, 3> local_counts{
        static_cast<std::uint64_t>(
            partition.entity_count(mpmc::mesh::EntityKind::cell)),
        static_cast<std::uint64_t>(
            partition.entity_count(mpmc::mesh::EntityKind::face)),
        static_cast<std::uint64_t>(
            partition.entity_count(mpmc::mesh::EntityKind::vertex))};

    const std::size_t rank_count =
        static_cast<std::size_t>(partition.rank_count());
    if (rank_count >
        std::numeric_limits<std::size_t>::max() / std::size_t{3U}) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    std::vector<std::uint64_t> all_counts(rank_count * std::size_t{3U});

    if (MPI_Allgather(
            local_counts.data(), 3, MPI_UINT64_T,
            all_counts.data(), 3, MPI_UINT64_T,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    PetscInt nroots = 0;
    error = detail::checked_petsc_int_size(local_ranges.end, &nroots);
    if (error != PETSC_SUCCESS) return error;

    std::vector<PetscInt> local_leaves;
    std::vector<PetscSFNode> remote_roots;
    local_leaves.reserve(plan.receive_count());
    remote_roots.reserve(plan.receive_count());

    for (const auto& neighbor : plan.neighbors()) {
        PetscMPIInt remote_rank = 0;
        error = detail::checked_mpi_rank(neighbor.rank, &remote_rank);
        if (error != PETSC_SUCCESS) return error;

        const std::size_t remote_slot =
            static_cast<std::size_t>(neighbor.rank.value()) *
            std::size_t{3U};
        const std::uint64_t remote_cells = all_counts[remote_slot];
        const std::uint64_t remote_faces = all_counts[remote_slot + 1U];
        const std::uint64_t remote_vertices = all_counts[remote_slot + 2U];

        if (remote_faces >
            std::numeric_limits<std::uint64_t>::max() - remote_cells) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }
        const std::uint64_t remote_face_base = remote_cells;
        const std::uint64_t remote_vertex_base =
            remote_cells + remote_faces;

        for (const auto& entity :
             plan.receive_entities_from(neighbor.rank)) {
            std::size_t local_base = 0U;
            std::uint64_t remote_base = 0U;
            std::uint64_t remote_count = 0U;

            switch (entity.kind) {
            case mpmc::mesh::EntityKind::cell:
                local_base = local_ranges.cell_begin;
                remote_base = 0U;
                remote_count = remote_cells;
                break;
            case mpmc::mesh::EntityKind::face:
                local_base = local_ranges.face_begin;
                remote_base = remote_face_base;
                remote_count = remote_faces;
                break;
            case mpmc::mesh::EntityKind::vertex:
                local_base = local_ranges.vertex_begin;
                remote_base = remote_vertex_base;
                remote_count = remote_vertices;
                break;
            case mpmc::mesh::EntityKind::edge:
                continue;
            }

            if (!partition.is_ghost(entity.kind, entity.local) ||
                partition.owner_rank(entity.kind, entity.local) !=
                    neighbor.rank ||
                partition.global_id(entity.kind, entity.local) !=
                    entity.global_id) {
                return PETSC_ERR_ARG_INCOMP;
            }

            const std::uint64_t remote_local =
                static_cast<std::uint64_t>(
                    entity.remote_local.value());
            if (remote_local >= remote_count) {
                return PETSC_ERR_ARG_OUTOFRANGE;
            }
            if (remote_local >
                std::numeric_limits<std::uint64_t>::max() - remote_base) {
                return PETSC_ERR_ARG_OUTOFRANGE;
            }

            const std::size_t local_point =
                local_base +
                static_cast<std::size_t>(entity.local.value());
            const std::uint64_t remote_point =
                remote_base + remote_local;

            PetscInt local_leaf = 0;
            PetscInt remote_root = 0;
            error = detail::checked_petsc_int_size(
                local_point, &local_leaf);
            if (error != PETSC_SUCCESS) return error;
            error = detail::checked_petsc_int_u64(
                remote_point, &remote_root);
            if (error != PETSC_SUCCESS) return error;

            local_leaves.push_back(local_leaf);
            remote_roots.push_back(
                PetscSFNode{remote_rank, remote_root});
        }
    }

    PetscInt nleaves = 0;
    error = detail::checked_petsc_int_size(
        local_leaves.size(), &nleaves);
    if (error != PETSC_SUCCESS) return error;

    PetscSF point_sf = nullptr;
    error = PetscSFCreate(comm, &point_sf);
    if (error != PETSC_SUCCESS) return error;

    error = PetscSFSetGraph(
        point_sf,
        nroots,
        nleaves,
        local_leaves.empty() ? nullptr : local_leaves.data(),
        PETSC_COPY_VALUES,
        remote_roots.empty() ? nullptr : remote_roots.data(),
        PETSC_COPY_VALUES);
    if (error != PETSC_SUCCESS) {
        PetscSFDestroy(&point_sf);
        return error;
    }

    error = PetscSFSetUp(point_sf);
    if (error != PETSC_SUCCESS) {
        PetscSFDestroy(&point_sf);
        return error;
    }

    *sf = point_sf;
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


struct DMPlexPointIdentity {
    PetscInt point;
    mpmc::mesh::EntityKind kind;
    mpmc::mesh::LocalIndex local;
    mpmc::mesh::GlobalEntityId global;
};

inline PetscErrorCode create_serial_dmplex_topology(
    const mpmc::mesh::Topology& topology,
    DM* dm,
    std::vector<DMPlexPointIdentity>* identities) {
    if (dm == nullptr || identities == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    *dm = nullptr;
    identities->clear();

    if (topology.entity_count(mpmc::mesh::EntityKind::edge) != 0U) {
        return PETSC_ERR_SUP;
    }
    if (!topology.has_relation(
            mpmc::mesh::EntityKind::cell,
            mpmc::mesh::EntityKind::face) ||
        !topology.has_relation(
            mpmc::mesh::EntityKind::cell,
            mpmc::mesh::EntityKind::vertex) ||
        !topology.has_relation(
            mpmc::mesh::EntityKind::face,
            mpmc::mesh::EntityKind::vertex) ||
        !topology.has_relation(
            mpmc::mesh::EntityKind::face,
            mpmc::mesh::EntityKind::cell)) {
        return PETSC_ERR_ARG_INCOMP;
    }

    const std::size_t cell_count =
        topology.entity_count(mpmc::mesh::EntityKind::cell);
    const std::size_t face_count =
        topology.entity_count(mpmc::mesh::EntityKind::face);
    const std::size_t vertex_count =
        topology.entity_count(mpmc::mesh::EntityKind::vertex);
    if (cell_count == 0U ||
        face_count == 0U ||
        vertex_count == 0U) {
        return PETSC_ERR_ARG_SIZ;
    }

    const auto& cell_faces = topology.relation(
        mpmc::mesh::EntityKind::cell,
        mpmc::mesh::EntityKind::face);
    const auto& cell_vertices = topology.relation(
        mpmc::mesh::EntityKind::cell,
        mpmc::mesh::EntityKind::vertex);
    const auto& face_vertices = topology.relation(
        mpmc::mesh::EntityKind::face,
        mpmc::mesh::EntityKind::vertex);
    const auto& face_cells = topology.relation(
        mpmc::mesh::EntityKind::face,
        mpmc::mesh::EntityKind::cell);

    for (std::size_t cell = 0U; cell < cell_count; ++cell) {
        const auto local = mpmc::mesh::LocalIndex{
            static_cast<mpmc::mesh::LocalIndex::value_type>(cell)};
        if (cell_faces.adjacent(local).size() != 4U ||
            cell_vertices.adjacent(local).size() != 4U) {
            return PETSC_ERR_SUP;
        }
    }
    for (std::size_t face = 0U; face < face_count; ++face) {
        const auto local = mpmc::mesh::LocalIndex{
            static_cast<mpmc::mesh::LocalIndex::value_type>(face)};
        const auto vertices = face_vertices.adjacent(local);
        const auto cells = face_cells.adjacent(local);
        if (vertices.size() != 2U ||
            (cells.size() != 1U && cells.size() != 2U)) {
            return PETSC_ERR_SUP;
        }
    }

    std::vector<std::vector<mpmc::mesh::LocalIndex>> derived_face_cells(
        face_count);
    for (std::size_t cell = 0U; cell < cell_count; ++cell) {
        const auto cell_local = mpmc::mesh::LocalIndex{
            static_cast<mpmc::mesh::LocalIndex::value_type>(cell)};
        const auto faces = cell_faces.adjacent(cell_local);
        std::array<mpmc::mesh::LocalIndex, 4> union_vertices{
            faces[0], faces[0], faces[0], faces[0]};
        std::size_t unique_vertex_count = 0U;
        for (const auto face_local : faces) {
            const std::size_t face_position =
                static_cast<std::size_t>(face_local.value());
            if (face_position >= face_count) {
                return PETSC_ERR_ARG_OUTOFRANGE;
            }
            derived_face_cells[face_position].push_back(cell_local);
            for (const auto vertex_local :
                 face_vertices.adjacent(face_local)) {
                const bool already_present =
                    std::find(
                        union_vertices.begin(),
                        union_vertices.begin() +
                            static_cast<std::ptrdiff_t>(
                                unique_vertex_count),
                        vertex_local) !=
                    union_vertices.begin() +
                        static_cast<std::ptrdiff_t>(
                            unique_vertex_count);
                if (!already_present) {
                    if (unique_vertex_count >=
                        union_vertices.size()) {
                        return PETSC_ERR_ARG_INCOMP;
                    }
                    union_vertices[unique_vertex_count++] =
                        vertex_local;
                }
            }
        }
        if (unique_vertex_count != 4U) {
            return PETSC_ERR_ARG_INCOMP;
        }

        auto expected_vertices = cell_vertices.adjacent(cell_local);
        std::array<mpmc::mesh::LocalIndex, 4> expected{
            expected_vertices[0],
            expected_vertices[1],
            expected_vertices[2],
            expected_vertices[3]};
        std::sort(expected.begin(), expected.end());
        std::sort(union_vertices.begin(), union_vertices.end());
        if (expected != union_vertices) {
            return PETSC_ERR_ARG_INCOMP;
        }
    }

    for (std::size_t face = 0U; face < face_count; ++face) {
        const auto face_local = mpmc::mesh::LocalIndex{
            static_cast<mpmc::mesh::LocalIndex::value_type>(face)};
        auto expected = derived_face_cells[face];
        auto actual_span = face_cells.adjacent(face_local);
        std::vector<mpmc::mesh::LocalIndex> actual(
            actual_span.begin(), actual_span.end());
        std::sort(expected.begin(), expected.end());
        std::sort(actual.begin(), actual.end());
        if (expected != actual) {
            return PETSC_ERR_ARG_INCOMP;
        }
    }

    if (face_count >
        std::numeric_limits<std::size_t>::max() - cell_count) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    const std::size_t face_base = cell_count;
    const std::size_t vertex_base = cell_count + face_count;
    if (vertex_count >
        std::numeric_limits<std::size_t>::max() - vertex_base) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    const std::size_t point_count = vertex_base + vertex_count;

    PetscInt chart_end = 0;
    PetscErrorCode error =
        detail::checked_petsc_int_size(point_count, &chart_end);
    if (error != PETSC_SUCCESS) return error;

    DM plex = nullptr;
    error = DMPlexCreate(PETSC_COMM_SELF, &plex);
    if (error != PETSC_SUCCESS) return error;

    error = DMSetDimension(plex, 2);
    if (error != PETSC_SUCCESS) {
        DMDestroy(&plex);
        return error;
    }
    error = DMPlexSetChart(plex, 0, chart_end);
    if (error != PETSC_SUCCESS) {
        DMDestroy(&plex);
        return error;
    }

    for (std::size_t cell = 0U; cell < cell_count; ++cell) {
        PetscInt point = 0;
        error = detail::checked_petsc_int_size(cell, &point);
        if (error != PETSC_SUCCESS) {
            DMDestroy(&plex);
            return error;
        }
        error = DMPlexSetConeSize(plex, point, 4);
        if (error != PETSC_SUCCESS) {
            DMDestroy(&plex);
            return error;
        }
    }
    for (std::size_t face = 0U; face < face_count; ++face) {
        PetscInt point = 0;
        error = detail::checked_petsc_int_size(
            face_base + face, &point);
        if (error != PETSC_SUCCESS) {
            DMDestroy(&plex);
            return error;
        }
        error = DMPlexSetConeSize(plex, point, 2);
        if (error != PETSC_SUCCESS) {
            DMDestroy(&plex);
            return error;
        }
    }

    error = DMSetUp(plex);
    if (error != PETSC_SUCCESS) {
        DMDestroy(&plex);
        return error;
    }

    for (std::size_t cell = 0U; cell < cell_count; ++cell) {
        const auto local = mpmc::mesh::LocalIndex{
            static_cast<mpmc::mesh::LocalIndex::value_type>(cell)};
        const auto faces = cell_faces.adjacent(local);
        std::array<PetscInt, 4> cone{};
        for (std::size_t i = 0U; i < cone.size(); ++i) {
            error = detail::checked_petsc_int_size(
                face_base +
                    static_cast<std::size_t>(
                        faces[i].value()),
                &cone[i]);
            if (error != PETSC_SUCCESS) {
                DMDestroy(&plex);
                return error;
            }
        }
        PetscInt point = 0;
        error = detail::checked_petsc_int_size(cell, &point);
        if (error != PETSC_SUCCESS) {
            DMDestroy(&plex);
            return error;
        }
        error = DMPlexSetCone(plex, point, cone.data());
        if (error != PETSC_SUCCESS) {
            DMDestroy(&plex);
            return error;
        }
    }

    for (std::size_t face = 0U; face < face_count; ++face) {
        const auto local = mpmc::mesh::LocalIndex{
            static_cast<mpmc::mesh::LocalIndex::value_type>(face)};
        const auto vertices = face_vertices.adjacent(local);
        std::array<PetscInt, 2> cone{};
        for (std::size_t i = 0U; i < cone.size(); ++i) {
            error = detail::checked_petsc_int_size(
                vertex_base +
                    static_cast<std::size_t>(
                        vertices[i].value()),
                &cone[i]);
            if (error != PETSC_SUCCESS) {
                DMDestroy(&plex);
                return error;
            }
        }
        PetscInt point = 0;
        error = detail::checked_petsc_int_size(
            face_base + face, &point);
        if (error != PETSC_SUCCESS) {
            DMDestroy(&plex);
            return error;
        }
        error = DMPlexSetCone(plex, point, cone.data());
        if (error != PETSC_SUCCESS) {
            DMDestroy(&plex);
            return error;
        }
    }

    error = DMPlexSymmetrize(plex);
    if (error != PETSC_SUCCESS) {
        DMDestroy(&plex);
        return error;
    }
    error = DMPlexStratify(plex);
    if (error != PETSC_SUCCESS) {
        DMDestroy(&plex);
        return error;
    }
    error = DMPlexComputeCellTypes(plex);
    if (error != PETSC_SUCCESS) {
        DMDestroy(&plex);
        return error;
    }

    identities->reserve(point_count);
    for (std::size_t cell = 0U; cell < cell_count; ++cell) {
        PetscInt point = 0;
        error = detail::checked_petsc_int_size(cell, &point);
        if (error != PETSC_SUCCESS) {
            identities->clear();
            DMDestroy(&plex);
            return error;
        }
        const auto local = mpmc::mesh::LocalIndex{
            static_cast<mpmc::mesh::LocalIndex::value_type>(cell)};
        identities->push_back(DMPlexPointIdentity{
            point,
            mpmc::mesh::EntityKind::cell,
            local,
            topology.global_id(
                mpmc::mesh::EntityKind::cell, local)});
    }
    for (std::size_t face = 0U; face < face_count; ++face) {
        PetscInt point = 0;
        error = detail::checked_petsc_int_size(
            face_base + face, &point);
        if (error != PETSC_SUCCESS) {
            identities->clear();
            DMDestroy(&plex);
            return error;
        }
        const auto local = mpmc::mesh::LocalIndex{
            static_cast<mpmc::mesh::LocalIndex::value_type>(face)};
        identities->push_back(DMPlexPointIdentity{
            point,
            mpmc::mesh::EntityKind::face,
            local,
            topology.global_id(
                mpmc::mesh::EntityKind::face, local)});
    }
    for (std::size_t vertex = 0U;
         vertex < vertex_count;
         ++vertex) {
        PetscInt point = 0;
        error = detail::checked_petsc_int_size(
            vertex_base + vertex, &point);
        if (error != PETSC_SUCCESS) {
            identities->clear();
            DMDestroy(&plex);
            return error;
        }
        const auto local = mpmc::mesh::LocalIndex{
            static_cast<mpmc::mesh::LocalIndex::value_type>(vertex)};
        identities->push_back(DMPlexPointIdentity{
            point,
            mpmc::mesh::EntityKind::vertex,
            local,
            topology.global_id(
                mpmc::mesh::EntityKind::vertex, local)});
    }

    *dm = plex;
    return PETSC_SUCCESS;
}

inline PetscErrorCode attach_root_geometry2d_coordinates(
    DM dm,
    PetscMPIInt root_rank,
    const mpmc::mesh::Geometry2D* root_geometry,
    std::span<const DMPlexPointIdentity> identities) {
    if (dm == nullptr) return PETSC_ERR_ARG_NULL;

    MPI_Comm comm =
        PetscObjectComm(reinterpret_cast<PetscObject>(dm));
    int mpi_rank = -1;
    int mpi_size = -1;
    if (MPI_Comm_rank(comm, &mpi_rank) != MPI_SUCCESS ||
        MPI_Comm_size(comm, &mpi_size) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    if (root_rank < 0 || root_rank >= mpi_size) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    PetscInt chart_start = 0;
    PetscInt chart_end = 0;
    PetscErrorCode error =
        DMPlexGetChart(dm, &chart_start, &chart_end);
    if (error != PETSC_SUCCESS) return error;

    PetscErrorCode local_validation = PETSC_SUCCESS;
    if (mpi_rank == root_rank) {
        if (root_geometry == nullptr) {
            local_validation = PETSC_ERR_ARG_NULL;
        } else {
            std::size_t vertex_identity_count = 0U;
            for (const auto& identity : identities) {
                if (identity.point < chart_start ||
                    identity.point >= chart_end) {
                    local_validation = PETSC_ERR_ARG_OUTOFRANGE;
                    break;
                }
                if (identity.kind ==
                    mpmc::mesh::EntityKind::vertex) {
                    ++vertex_identity_count;
                    if (static_cast<std::size_t>(
                            identity.local.value()) >=
                        root_geometry->vertex_count()) {
                        local_validation = PETSC_ERR_ARG_SIZ;
                        break;
                    }
                }
            }
            if (local_validation == PETSC_SUCCESS &&
                vertex_identity_count !=
                    root_geometry->vertex_count()) {
                local_validation = PETSC_ERR_ARG_SIZ;
            }
        }
    } else if (!identities.empty()) {
        local_validation = PETSC_ERR_ARG_INCOMP;
    }

    int validation_code =
        static_cast<int>(local_validation);
    if (MPI_Bcast(
            &validation_code,
            1,
            MPI_INT,
            root_rank,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    if (validation_code !=
        static_cast<int>(PETSC_SUCCESS)) {
        return static_cast<PetscErrorCode>(
            validation_code);
    }

    error = DMSetCoordinateDim(dm, 2);
    if (error != PETSC_SUCCESS) return error;

    PetscSection coordinate_section = nullptr;
    error = PetscSectionCreate(comm, &coordinate_section);
    if (error != PETSC_SUCCESS) return error;

    error = PetscSectionSetChart(
        coordinate_section,
        chart_start,
        chart_end);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&coordinate_section);
        return error;
    }

    for (const auto& identity : identities) {
        if (identity.kind !=
            mpmc::mesh::EntityKind::vertex) {
            continue;
        }
        error = PetscSectionSetDof(
            coordinate_section,
            identity.point,
            2);
        if (error != PETSC_SUCCESS) {
            PetscSectionDestroy(&coordinate_section);
            return error;
        }
    }
    error = PetscSectionSetUp(coordinate_section);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&coordinate_section);
        return error;
    }

    PetscInt coordinate_storage = 0;
    error = PetscSectionGetStorageSize(
        coordinate_section,
        &coordinate_storage);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&coordinate_section);
        return error;
    }
    if (coordinate_storage < 0) {
        PetscSectionDestroy(&coordinate_section);
        return PETSC_ERR_PLIB;
    }

    Vec coordinates = nullptr;
    error = VecCreateSeq(
        PETSC_COMM_SELF,
        coordinate_storage,
        &coordinates);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&coordinate_section);
        return error;
    }
    error = VecSetBlockSize(coordinates, 2);
    if (error != PETSC_SUCCESS) {
        VecDestroy(&coordinates);
        PetscSectionDestroy(&coordinate_section);
        return error;
    }

    if (mpi_rank == root_rank) {
        PetscScalar* values = nullptr;
        error = VecGetArray(coordinates, &values);
        if (error != PETSC_SUCCESS) {
            VecDestroy(&coordinates);
            PetscSectionDestroy(&coordinate_section);
            return error;
        }
        for (const auto& identity : identities) {
            if (identity.kind !=
                mpmc::mesh::EntityKind::vertex) {
                continue;
            }
            PetscInt offset = -1;
            error = PetscSectionGetOffset(
                coordinate_section,
                identity.point,
                &offset);
            if (error != PETSC_SUCCESS ||
                offset < 0 ||
                offset + 1 >= coordinate_storage) {
                VecRestoreArray(coordinates, &values);
                VecDestroy(&coordinates);
                PetscSectionDestroy(&coordinate_section);
                return error != PETSC_SUCCESS
                           ? error
                           : PETSC_ERR_PLIB;
            }
            const auto coordinate =
                root_geometry->vertex_coordinate_m(
                    identity.local);
            values[static_cast<std::size_t>(offset)] =
                static_cast<PetscScalar>(coordinate.x_m);
            values[static_cast<std::size_t>(offset + 1)] =
                static_cast<PetscScalar>(coordinate.y_m);
        }
        error = VecRestoreArray(coordinates, &values);
        if (error != PETSC_SUCCESS) {
            VecDestroy(&coordinates);
            PetscSectionDestroy(&coordinate_section);
            return error;
        }
    }

    error = DMSetCoordinateSection(
        dm, 2, coordinate_section);
    if (error != PETSC_SUCCESS) {
        VecDestroy(&coordinates);
        PetscSectionDestroy(&coordinate_section);
        return error;
    }
    error = DMSetCoordinatesLocal(dm, coordinates);
    if (error != PETSC_SUCCESS) {
        VecDestroy(&coordinates);
        PetscSectionDestroy(&coordinate_section);
        return error;
    }

    const PetscErrorCode vector_destroy_error =
        VecDestroy(&coordinates);
    const PetscErrorCode section_destroy_error =
        PetscSectionDestroy(&coordinate_section);
    if (vector_destroy_error != PETSC_SUCCESS) {
        return vector_destroy_error;
    }
    return section_destroy_error;
}

inline PetscErrorCode create_root_dmplex_topology(
    MPI_Comm comm,
    PetscMPIInt root_rank,
    const mpmc::mesh::Topology* root_topology,
    DM* dm,
    std::vector<DMPlexPointIdentity>* identities) {
    if (dm == nullptr || identities == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    *dm = nullptr;
    identities->clear();

    int mpi_rank = -1;
    int mpi_size = -1;
    if (MPI_Comm_rank(comm, &mpi_rank) != MPI_SUCCESS ||
        MPI_Comm_size(comm, &mpi_size) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    if (root_rank < 0 || root_rank >= mpi_size) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    DM serial_dm = nullptr;
    std::vector<DMPlexPointIdentity> serial_identities;
    PetscErrorCode root_error = PETSC_SUCCESS;
    if (mpi_rank == root_rank) {
        if (root_topology == nullptr) {
            root_error = PETSC_ERR_ARG_NULL;
        } else {
            root_error = create_serial_dmplex_topology(
                *root_topology, &serial_dm, &serial_identities);
        }
    }

    int broadcast_error = static_cast<int>(root_error);
    if (MPI_Bcast(
            &broadcast_error, 1, MPI_INT,
            root_rank, comm) != MPI_SUCCESS) {
        if (serial_dm != nullptr) DMDestroy(&serial_dm);
        return PETSC_ERR_MPI;
    }
    if (broadcast_error != static_cast<int>(PETSC_SUCCESS)) {
        if (serial_dm != nullptr) DMDestroy(&serial_dm);
        return static_cast<PetscErrorCode>(broadcast_error);
    }

    PetscErrorCode error = PETSC_SUCCESS;
    PetscInt chart_start = 0;
    PetscInt chart_end = 0;
    if (mpi_rank == root_rank) {
        error = DMPlexGetChart(
            serial_dm, &chart_start, &chart_end);
        if (error != PETSC_SUCCESS) {
            DMDestroy(&serial_dm);
            return error;
        }
        if (chart_start != 0) {
            DMDestroy(&serial_dm);
            return PETSC_ERR_ARG_INCOMP;
        }
    }

    DM rooted_dm = nullptr;
    error = DMPlexCreate(comm, &rooted_dm);
    if (error != PETSC_SUCCESS) {
        if (serial_dm != nullptr) DMDestroy(&serial_dm);
        return error;
    }

    error = DMSetDimension(rooted_dm, 2);
    if (error != PETSC_SUCCESS) {
        DMDestroy(&rooted_dm);
        if (serial_dm != nullptr) DMDestroy(&serial_dm);
        return error;
    }
    error = DMPlexSetChart(rooted_dm, 0, chart_end);
    if (error != PETSC_SUCCESS) {
        DMDestroy(&rooted_dm);
        if (serial_dm != nullptr) DMDestroy(&serial_dm);
        return error;
    }

    if (mpi_rank == root_rank) {
        for (PetscInt point = 0; point < chart_end; ++point) {
            PetscInt cone_size = 0;
            error = DMPlexGetConeSize(
                serial_dm, point, &cone_size);
            if (error != PETSC_SUCCESS) {
                DMDestroy(&rooted_dm);
                DMDestroy(&serial_dm);
                return error;
            }
            error = DMPlexSetConeSize(
                rooted_dm, point, cone_size);
            if (error != PETSC_SUCCESS) {
                DMDestroy(&rooted_dm);
                DMDestroy(&serial_dm);
                return error;
            }
        }
    }

    error = DMSetUp(rooted_dm);
    if (error != PETSC_SUCCESS) {
        DMDestroy(&rooted_dm);
        if (serial_dm != nullptr) DMDestroy(&serial_dm);
        return error;
    }

    if (mpi_rank == root_rank) {
        for (PetscInt point = 0; point < chart_end; ++point) {
            PetscInt cone_size = 0;
            const PetscInt* cone = nullptr;
            const PetscInt* orientation = nullptr;
            error = DMPlexGetConeSize(
                serial_dm, point, &cone_size);
            if (error != PETSC_SUCCESS) {
                DMDestroy(&rooted_dm);
                DMDestroy(&serial_dm);
                return error;
            }
            if (cone_size == 0) continue;

            error = DMPlexGetCone(
                serial_dm, point, &cone);
            if (error != PETSC_SUCCESS) {
                DMDestroy(&rooted_dm);
                DMDestroy(&serial_dm);
                return error;
            }
            error = DMPlexGetConeOrientation(
                serial_dm, point, &orientation);
            if (error != PETSC_SUCCESS) {
                DMDestroy(&rooted_dm);
                DMDestroy(&serial_dm);
                return error;
            }
            error = DMPlexSetCone(
                rooted_dm, point, cone);
            if (error != PETSC_SUCCESS) {
                DMDestroy(&rooted_dm);
                DMDestroy(&serial_dm);
                return error;
            }
            error = DMPlexSetConeOrientation(
                rooted_dm, point, orientation);
            if (error != PETSC_SUCCESS) {
                DMDestroy(&rooted_dm);
                DMDestroy(&serial_dm);
                return error;
            }
        }
    }

    error = DMPlexSymmetrize(rooted_dm);
    if (error != PETSC_SUCCESS) {
        DMDestroy(&rooted_dm);
        if (serial_dm != nullptr) DMDestroy(&serial_dm);
        return error;
    }
    error = DMPlexStratify(rooted_dm);
    if (error != PETSC_SUCCESS) {
        DMDestroy(&rooted_dm);
        if (serial_dm != nullptr) DMDestroy(&serial_dm);
        return error;
    }
    error = DMPlexComputeCellTypes(rooted_dm);
    if (error != PETSC_SUCCESS) {
        DMDestroy(&rooted_dm);
        if (serial_dm != nullptr) DMDestroy(&serial_dm);
        return error;
    }

    if (mpi_rank == root_rank) {
        *identities = std::move(serial_identities);
    }

    if (serial_dm != nullptr) {
        error = DMDestroy(&serial_dm);
        if (error != PETSC_SUCCESS) {
            identities->clear();
            DMDestroy(&rooted_dm);
            return error;
        }
    }

    *dm = rooted_dm;
    return PETSC_SUCCESS;
}

inline PetscErrorCode migrate_dmplex_identities(
    DM source_dm,
    PetscSF migration_sf,
    std::span<const DMPlexPointIdentity> source_identities,
    DM target_dm,
    std::vector<DMPlexPointIdentity>* target_identities) {
    if (source_dm == nullptr ||
        migration_sf == nullptr ||
        target_dm == nullptr ||
        target_identities == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    target_identities->clear();

    PetscInt source_start = 0;
    PetscInt source_end = 0;
    PetscInt target_start = 0;
    PetscInt target_end = 0;
    PetscErrorCode error = DMPlexGetChart(
        source_dm, &source_start, &source_end);
    if (error != PETSC_SUCCESS) return error;
    error = DMPlexGetChart(
        target_dm, &target_start, &target_end);
    if (error != PETSC_SUCCESS) return error;
    if (source_start < 0 || source_end < source_start ||
        target_start < 0 || target_end < target_start) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    const std::size_t source_count =
        static_cast<std::size_t>(
            source_end - source_start);
    if (source_identities.size() != source_count) {
        return PETSC_ERR_ARG_SIZ;
    }

    PetscSection source_section = nullptr;
    PetscSection target_section = nullptr;
    error = PetscSectionCreate(
        PetscObjectComm(
            reinterpret_cast<PetscObject>(source_dm)),
        &source_section);
    if (error != PETSC_SUCCESS) return error;
    error = PetscSectionCreate(
        PetscObjectComm(
            reinterpret_cast<PetscObject>(target_dm)),
        &target_section);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&source_section);
        return error;
    }

    error = PetscSectionSetChart(
        source_section, source_start, source_end);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return error;
    }
    for (PetscInt point = source_start;
         point < source_end;
         ++point) {
        error = PetscSectionSetDof(
            source_section, point, 2);
        if (error != PETSC_SUCCESS) {
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return error;
        }
    }
    error = PetscSectionSetUp(source_section);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return error;
    }

    PetscInt source_storage = 0;
    error = PetscSectionGetStorageSize(
        source_section, &source_storage);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return error;
    }
    if (source_storage < 0 ||
        source_storage !=
            static_cast<PetscInt>(
                source_count * std::size_t{2U})) {
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return PETSC_ERR_PLIB;
    }

    std::vector<std::uint8_t> seen(
        source_count, std::uint8_t{0U});
    std::vector<std::uint64_t> source_data(
        static_cast<std::size_t>(source_storage), 0U);
    for (const auto& identity : source_identities) {
        if (identity.point < source_start ||
            identity.point >= source_end) {
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return PETSC_ERR_ARG_OUTOFRANGE;
        }
        const std::size_t position =
            static_cast<std::size_t>(
                identity.point - source_start);
        if (seen[position] != 0U) {
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return PETSC_ERR_ARG_INCOMP;
        }
        seen[position] = std::uint8_t{1U};

        PetscInt offset = 0;
        error = PetscSectionGetOffset(
            source_section, identity.point, &offset);
        if (error != PETSC_SUCCESS) {
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return error;
        }
        if (offset < 0 ||
            offset + 1 >= source_storage) {
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return PETSC_ERR_PLIB;
        }
        source_data[
            static_cast<std::size_t>(offset)] =
            static_cast<std::uint64_t>(identity.kind);
        source_data[
            static_cast<std::size_t>(offset + 1)] =
            identity.global.value();
    }
    if (std::find(
            seen.begin(), seen.end(),
            std::uint8_t{0U}) != seen.end()) {
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return PETSC_ERR_ARG_INCOMP;
    }

    void* target_raw = nullptr;
    error = DMPlexDistributeData(
        source_dm,
        migration_sf,
        source_section,
        MPI_UINT64_T,
        source_data.empty()
            ? nullptr
            : source_data.data(),
        target_section,
        &target_raw);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return error;
    }

    PetscInt target_storage = 0;
    error = PetscSectionGetStorageSize(
        target_section, &target_storage);
    if (error != PETSC_SUCCESS) {
        PetscFree(target_raw);
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return error;
    }

    const std::size_t target_count =
        static_cast<std::size_t>(
            target_end - target_start);
    if (target_storage < 0 ||
        target_storage !=
            static_cast<PetscInt>(
                target_count * std::size_t{2U})) {
        PetscFree(target_raw);
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return PETSC_ERR_PLIB;
    }

    PetscInt cell_start = -1;
    PetscInt cell_end = -1;
    PetscInt face_start = -1;
    PetscInt face_end = -1;
    PetscInt vertex_start = -1;
    PetscInt vertex_end = -1;
    error = DMPlexGetHeightStratum(
        target_dm, 0, &cell_start, &cell_end);
    if (error == PETSC_SUCCESS) {
        error = DMPlexGetHeightStratum(
            target_dm, 1, &face_start, &face_end);
    }
    if (error == PETSC_SUCCESS) {
        error = DMPlexGetDepthStratum(
            target_dm, 0, &vertex_start, &vertex_end);
    }
    if (error != PETSC_SUCCESS) {
        PetscFree(target_raw);
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return error;
    }

    auto* target_data =
        static_cast<std::uint64_t*>(target_raw);
    target_identities->reserve(target_count);
    for (PetscInt point = target_start;
         point < target_end;
         ++point) {
        PetscInt dof = 0;
        PetscInt offset = 0;
        error = PetscSectionGetDof(
            target_section, point, &dof);
        if (error == PETSC_SUCCESS) {
            error = PetscSectionGetOffset(
                target_section, point, &offset);
        }
        if (error != PETSC_SUCCESS) {
            target_identities->clear();
            PetscFree(target_raw);
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return error;
        }
        if (dof != 2 || offset < 0 ||
            offset + 1 >= target_storage) {
            target_identities->clear();
            PetscFree(target_raw);
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return PETSC_ERR_PLIB;
        }

        const std::uint64_t raw_kind =
            target_data[
                static_cast<std::size_t>(offset)];
        mpmc::mesh::EntityKind kind;
        PetscInt kind_start = -1;
        PetscInt kind_end = -1;
        switch (raw_kind) {
        case static_cast<std::uint64_t>(
                 mpmc::mesh::EntityKind::cell):
            kind = mpmc::mesh::EntityKind::cell;
            kind_start = cell_start;
            kind_end = cell_end;
            break;
        case static_cast<std::uint64_t>(
                 mpmc::mesh::EntityKind::face):
            kind = mpmc::mesh::EntityKind::face;
            kind_start = face_start;
            kind_end = face_end;
            break;
        case static_cast<std::uint64_t>(
                 mpmc::mesh::EntityKind::vertex):
            kind = mpmc::mesh::EntityKind::vertex;
            kind_start = vertex_start;
            kind_end = vertex_end;
            break;
        default:
            target_identities->clear();
            PetscFree(target_raw);
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return PETSC_ERR_ARG_INCOMP;
        }

        if (point < kind_start || point >= kind_end) {
            target_identities->clear();
            PetscFree(target_raw);
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return PETSC_ERR_ARG_INCOMP;
        }
        const PetscInt local_value =
            point - kind_start;
        if (local_value < 0 ||
            static_cast<std::uint64_t>(local_value) >
                static_cast<std::uint64_t>(
                    std::numeric_limits<
                        mpmc::mesh::LocalIndex::value_type>::max())) {
            target_identities->clear();
            PetscFree(target_raw);
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return PETSC_ERR_ARG_OUTOFRANGE;
        }

        target_identities->push_back(
            DMPlexPointIdentity{
                point,
                kind,
                mpmc::mesh::LocalIndex{
                    static_cast<
                        mpmc::mesh::LocalIndex::value_type>(
                            local_value)},
                mpmc::mesh::GlobalEntityId{
                    target_data[
                        static_cast<std::size_t>(
                            offset + 1)]}});
    }

    const PetscErrorCode free_error =
        PetscFree(target_raw);
    const PetscErrorCode target_destroy_error =
        PetscSectionDestroy(&target_section);
    const PetscErrorCode source_destroy_error =
        PetscSectionDestroy(&source_section);
    if (free_error != PETSC_SUCCESS) return free_error;
    if (target_destroy_error != PETSC_SUCCESS) {
        return target_destroy_error;
    }
    return source_destroy_error;
}

inline PetscErrorCode migrate_face_boundary_snapshot(
    DM source_dm,
    PetscSF migration_sf,
    const mpmc::mesh::FaceBoundarySnapshot& source_boundary,
    std::span<const DMPlexPointIdentity> source_identities,
    DM target_dm,
    std::span<const DMPlexPointIdentity> target_identities,
    std::optional<mpmc::mesh::FaceBoundarySnapshot>* target_boundary) {
    if (source_dm == nullptr ||
        migration_sf == nullptr ||
        target_dm == nullptr ||
        target_boundary == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    target_boundary->reset();

    PetscInt source_start = 0;
    PetscInt source_end = 0;
    PetscErrorCode error =
        DMPlexGetChart(
            source_dm, &source_start, &source_end);
    if (error != PETSC_SUCCESS) return error;

    PetscSection source_section = nullptr;
    PetscSection target_section = nullptr;
    error = PetscSectionCreate(
        PetscObjectComm(
            reinterpret_cast<PetscObject>(source_dm)),
        &source_section);
    if (error != PETSC_SUCCESS) return error;
    error = PetscSectionCreate(
        PetscObjectComm(
            reinterpret_cast<PetscObject>(target_dm)),
        &target_section);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&source_section);
        return error;
    }

    error = PetscSectionSetChart(
        source_section, source_start, source_end);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return error;
    }
    for (const auto& identity : source_identities) {
        if (identity.kind !=
            mpmc::mesh::EntityKind::face) {
            continue;
        }
        if (identity.point < source_start ||
            identity.point >= source_end ||
            static_cast<std::size_t>(
                identity.local.value()) >=
                source_boundary.face_count()) {
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return PETSC_ERR_ARG_SIZ;
        }
        error = PetscSectionSetDof(
            source_section, identity.point, 2);
        if (error != PETSC_SUCCESS) {
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return error;
        }
    }
    error = PetscSectionSetUp(source_section);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return error;
    }

    PetscInt source_storage = 0;
    error = PetscSectionGetStorageSize(
        source_section, &source_storage);
    if (error != PETSC_SUCCESS || source_storage < 0) {
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return error != PETSC_SUCCESS
                   ? error
                   : PETSC_ERR_PLIB;
    }

    std::vector<std::uint32_t> source_data(
        static_cast<std::size_t>(source_storage), 0U);
    for (const auto& identity : source_identities) {
        if (identity.kind !=
            mpmc::mesh::EntityKind::face) {
            continue;
        }
        PetscInt offset = -1;
        error = PetscSectionGetOffset(
            source_section,
            identity.point,
            &offset);
        if (error != PETSC_SUCCESS ||
            offset < 0 ||
            offset + 1 >= source_storage) {
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return error != PETSC_SUCCESS
                       ? error
                       : PETSC_ERR_PLIB;
        }
        const auto classification =
            source_boundary.classification(
                identity.local);
        const auto tag =
            source_boundary.physical_tag(
                identity.local);
        source_data[
            static_cast<std::size_t>(offset)] =
            static_cast<std::uint32_t>(
                classification);
        source_data[
            static_cast<std::size_t>(offset + 1)] =
            tag.value();
    }

    void* target_raw = nullptr;
    error = DMPlexDistributeData(
        source_dm,
        migration_sf,
        source_section,
        MPI_UINT32_T,
        source_data.empty()
            ? nullptr
            : source_data.data(),
        target_section,
        &target_raw);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return error;
    }

    PetscInt face_start = -1;
    PetscInt face_end = -1;
    error = DMPlexGetHeightStratum(
        target_dm, 1, &face_start, &face_end);
    if (error != PETSC_SUCCESS ||
        face_start < 0 ||
        face_end < face_start) {
        PetscFree(target_raw);
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return error != PETSC_SUCCESS
                   ? error
                   : PETSC_ERR_PLIB;
    }
    const std::size_t face_count =
        static_cast<std::size_t>(
            face_end - face_start);

    std::vector<mpmc::mesh::FaceClassification>
        classifications(
            face_count,
            mpmc::mesh::FaceClassification::interior);
    std::vector<mpmc::mesh::PhysicalTag>
        physical_tags(
            face_count,
            mpmc::mesh::PhysicalTag{0U});
    std::vector<std::uint8_t> seen(
        face_count, std::uint8_t{0U});

    auto* target_data =
        static_cast<std::uint32_t*>(target_raw);
    for (const auto& identity : target_identities) {
        if (identity.kind !=
            mpmc::mesh::EntityKind::face) {
            continue;
        }
        const std::size_t local =
            static_cast<std::size_t>(
                identity.local.value());
        if (local >= face_count ||
            seen[local] != 0U) {
            PetscFree(target_raw);
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return PETSC_ERR_ARG_INCOMP;
        }

        PetscInt dof = 0;
        PetscInt offset = -1;
        error = PetscSectionGetDof(
            target_section,
            identity.point,
            &dof);
        if (error == PETSC_SUCCESS) {
            error = PetscSectionGetOffset(
                target_section,
                identity.point,
                &offset);
        }
        if (error != PETSC_SUCCESS ||
            dof != 2 ||
            offset < 0) {
            PetscFree(target_raw);
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return error != PETSC_SUCCESS
                       ? error
                       : PETSC_ERR_PLIB;
        }

        const std::uint32_t raw_classification =
            target_data[
                static_cast<std::size_t>(offset)];
        mpmc::mesh::FaceClassification classification;
        if (raw_classification ==
            static_cast<std::uint32_t>(
                mpmc::mesh::FaceClassification::interior)) {
            classification =
                mpmc::mesh::FaceClassification::interior;
        } else if (
            raw_classification ==
            static_cast<std::uint32_t>(
                mpmc::mesh::FaceClassification::boundary)) {
            classification =
                mpmc::mesh::FaceClassification::boundary;
        } else {
            PetscFree(target_raw);
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return PETSC_ERR_ARG_INCOMP;
        }

        const auto tag =
            mpmc::mesh::PhysicalTag{
                target_data[
                    static_cast<std::size_t>(
                        offset + 1)]};
        if (classification ==
                mpmc::mesh::FaceClassification::interior &&
            tag.is_tagged()) {
            PetscFree(target_raw);
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return PETSC_ERR_ARG_INCOMP;
        }

        classifications[local] = classification;
        physical_tags[local] = tag;
        seen[local] = std::uint8_t{1U};
    }

    if (std::find(
            seen.begin(),
            seen.end(),
            std::uint8_t{0U}) != seen.end()) {
        PetscFree(target_raw);
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return PETSC_ERR_ARG_INCOMP;
    }

    target_boundary->emplace(
        std::move(classifications),
        std::move(physical_tags));

    const PetscErrorCode free_error =
        PetscFree(target_raw);
    const PetscErrorCode target_destroy_error =
        PetscSectionDestroy(&target_section);
    const PetscErrorCode source_destroy_error =
        PetscSectionDestroy(&source_section);
    if (free_error != PETSC_SUCCESS) return free_error;
    if (target_destroy_error != PETSC_SUCCESS) {
        return target_destroy_error;
    }
    return source_destroy_error;
}

inline PetscErrorCode migrate_dense_field_snapshot(
    DM source_dm,
    PetscSF migration_sf,
    const mpmc::mesh::DenseFieldSnapshot& source_field,
    std::span<const DMPlexPointIdentity> source_identities,
    DM target_dm,
    std::span<const DMPlexPointIdentity> target_identities,
    std::optional<mpmc::mesh::DenseFieldSnapshot>* target_field) {
    if (source_dm == nullptr ||
        migration_sf == nullptr ||
        target_dm == nullptr ||
        target_field == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    target_field->reset();

    const auto location = source_field.location();
    if (location == mpmc::mesh::EntityKind::edge) {
        return PETSC_ERR_SUP;
    }
    const std::size_t component_count =
        source_field.component_count();
    PetscInt component_count_petsc = 0;
    PetscErrorCode error =
        detail::checked_petsc_int_size(
            component_count,
            &component_count_petsc);
    if (error != PETSC_SUCCESS) return error;

    PetscInt source_start = 0;
    PetscInt source_end = 0;
    error = DMPlexGetChart(
        source_dm, &source_start, &source_end);
    if (error != PETSC_SUCCESS) return error;

    PetscSection source_section = nullptr;
    PetscSection target_section = nullptr;
    error = PetscSectionCreate(
        PetscObjectComm(
            reinterpret_cast<PetscObject>(source_dm)),
        &source_section);
    if (error != PETSC_SUCCESS) return error;
    error = PetscSectionCreate(
        PetscObjectComm(
            reinterpret_cast<PetscObject>(target_dm)),
        &target_section);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&source_section);
        return error;
    }

    error = PetscSectionSetChart(
        source_section, source_start, source_end);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return error;
    }
    for (const auto& identity : source_identities) {
        if (identity.kind != location) continue;
        if (identity.point < source_start ||
            identity.point >= source_end ||
            static_cast<std::size_t>(
                identity.local.value()) >=
                source_field.entity_count()) {
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return PETSC_ERR_ARG_SIZ;
        }
        error = PetscSectionSetDof(
            source_section,
            identity.point,
            component_count_petsc);
        if (error != PETSC_SUCCESS) {
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return error;
        }
    }
    error = PetscSectionSetUp(source_section);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return error;
    }

    PetscInt source_storage = 0;
    error = PetscSectionGetStorageSize(
        source_section, &source_storage);
    if (error != PETSC_SUCCESS ||
        source_storage < 0) {
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return error != PETSC_SUCCESS
                   ? error
                   : PETSC_ERR_PLIB;
    }

    std::vector<double> source_data(
        static_cast<std::size_t>(source_storage),
        0.0);
    for (const auto& identity : source_identities) {
        if (identity.kind != location) continue;
        PetscInt offset = -1;
        error = PetscSectionGetOffset(
            source_section,
            identity.point,
            &offset);
        if (error != PETSC_SUCCESS ||
            offset < 0 ||
            component_count_petsc >
                source_storage - offset) {
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return error != PETSC_SUCCESS
                       ? error
                       : PETSC_ERR_PLIB;
        }
        const auto values =
            source_field.entity_values(
                identity.local);
        for (std::size_t component = 0U;
             component < component_count;
             ++component) {
            source_data[
                static_cast<std::size_t>(offset) +
                component] =
                values[component];
        }
    }

    void* target_raw = nullptr;
    error = DMPlexDistributeData(
        source_dm,
        migration_sf,
        source_section,
        MPI_DOUBLE,
        source_data.empty()
            ? nullptr
            : source_data.data(),
        target_section,
        &target_raw);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return error;
    }

    PetscInt kind_start = -1;
    PetscInt kind_end = -1;
    switch (location) {
    case mpmc::mesh::EntityKind::cell:
        error = DMPlexGetHeightStratum(
            target_dm, 0, &kind_start, &kind_end);
        break;
    case mpmc::mesh::EntityKind::face:
        error = DMPlexGetHeightStratum(
            target_dm, 1, &kind_start, &kind_end);
        break;
    case mpmc::mesh::EntityKind::vertex:
        error = DMPlexGetDepthStratum(
            target_dm, 0, &kind_start, &kind_end);
        break;
    case mpmc::mesh::EntityKind::edge:
        error = PETSC_ERR_SUP;
        break;
    }
    if (error != PETSC_SUCCESS ||
        kind_start < 0 ||
        kind_end < kind_start) {
        PetscFree(target_raw);
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return error != PETSC_SUCCESS
                   ? error
                   : PETSC_ERR_PLIB;
    }

    const std::size_t entity_count =
        static_cast<std::size_t>(
            kind_end - kind_start);
    if (entity_count != 0U &&
        component_count >
            std::numeric_limits<std::size_t>::max() /
                entity_count) {
        PetscFree(target_raw);
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    std::vector<double> values(
        entity_count * component_count,
        0.0);
    std::vector<std::uint8_t> seen(
        entity_count, std::uint8_t{0U});

    auto* target_data =
        static_cast<double*>(target_raw);
    for (const auto& identity : target_identities) {
        if (identity.kind != location) continue;
        const std::size_t local =
            static_cast<std::size_t>(
                identity.local.value());
        if (local >= entity_count ||
            seen[local] != 0U) {
            PetscFree(target_raw);
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return PETSC_ERR_ARG_INCOMP;
        }

        PetscInt dof = 0;
        PetscInt offset = -1;
        error = PetscSectionGetDof(
            target_section,
            identity.point,
            &dof);
        if (error == PETSC_SUCCESS) {
            error = PetscSectionGetOffset(
                target_section,
                identity.point,
                &offset);
        }
        if (error != PETSC_SUCCESS ||
            dof != component_count_petsc ||
            offset < 0) {
            PetscFree(target_raw);
            PetscSectionDestroy(&target_section);
            PetscSectionDestroy(&source_section);
            return error != PETSC_SUCCESS
                       ? error
                       : PETSC_ERR_PLIB;
        }

        for (std::size_t component = 0U;
             component < component_count;
             ++component) {
            values[
                local * component_count +
                component] =
                target_data[
                    static_cast<std::size_t>(offset) +
                    component];
        }
        seen[local] = std::uint8_t{1U};
    }

    if (std::find(
            seen.begin(),
            seen.end(),
            std::uint8_t{0U}) != seen.end()) {
        PetscFree(target_raw);
        PetscSectionDestroy(&target_section);
        PetscSectionDestroy(&source_section);
        return PETSC_ERR_ARG_INCOMP;
    }

    mpmc::mesh::Topology::EntityIds ids;
    ids.cells.assign(
        location == mpmc::mesh::EntityKind::cell
            ? entity_count
            : 0U,
        mpmc::mesh::GlobalEntityId{0U});
    ids.faces.assign(
        location == mpmc::mesh::EntityKind::face
            ? entity_count
            : 0U,
        mpmc::mesh::GlobalEntityId{0U});
    ids.vertices.assign(
        location == mpmc::mesh::EntityKind::vertex
            ? entity_count
            : 0U,
        mpmc::mesh::GlobalEntityId{0U});

    auto assign_id =
        [&](const DMPlexPointIdentity& identity) {
            const std::size_t local =
                static_cast<std::size_t>(
                    identity.local.value());
            switch (identity.kind) {
            case mpmc::mesh::EntityKind::cell:
                if (location ==
                    mpmc::mesh::EntityKind::cell) {
                    ids.cells[local] =
                        identity.global;
                }
                break;
            case mpmc::mesh::EntityKind::face:
                if (location ==
                    mpmc::mesh::EntityKind::face) {
                    ids.faces[local] =
                        identity.global;
                }
                break;
            case mpmc::mesh::EntityKind::vertex:
                if (location ==
                    mpmc::mesh::EntityKind::vertex) {
                    ids.vertices[local] =
                        identity.global;
                }
                break;
            case mpmc::mesh::EntityKind::edge:
                break;
            }
        };
    for (const auto& identity : target_identities) {
        if (identity.kind == location) {
            assign_id(identity);
        }
    }

    const mpmc::mesh::Topology target_topology{
        std::move(ids), {}};
    target_field->emplace(
        mpmc::mesh::DenseFieldSnapshot::create(
            target_topology,
            location,
            component_count,
            std::move(values),
            source_field.metadata()));

    const PetscErrorCode free_error =
        PetscFree(target_raw);
    const PetscErrorCode target_destroy_error =
        PetscSectionDestroy(&target_section);
    const PetscErrorCode source_destroy_error =
        PetscSectionDestroy(&source_section);
    if (free_error != PETSC_SUCCESS) return free_error;
    if (target_destroy_error != PETSC_SUCCESS) {
        return target_destroy_error;
    }
    return source_destroy_error;
}

inline PetscErrorCode create_section_vecs(
    MPI_Comm comm,
    PetscSection local_section,
    PetscSection global_section,
    Vec* global_vec,
    Vec* local_vec) {
    if (global_vec == nullptr || local_vec == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    *global_vec = nullptr;
    *local_vec = nullptr;
    if (local_section == nullptr || global_section == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    PetscInt local_size = 0;
    PetscInt owned_size = 0;
    PetscErrorCode error =
        PetscSectionGetStorageSize(local_section, &local_size);
    if (error != PETSC_SUCCESS) return error;
    error = PetscSectionGetConstrainedStorageSize(
        global_section, &owned_size);
    if (error != PETSC_SUCCESS) return error;
    if (local_size < 0 || owned_size < 0) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    Vec global = nullptr;
    error = VecCreateMPI(
        comm, owned_size, PETSC_DETERMINE, &global);
    if (error != PETSC_SUCCESS) return error;

    Vec local = nullptr;
    error = VecCreateSeq(PETSC_COMM_SELF, local_size, &local);
    if (error != PETSC_SUCCESS) {
        VecDestroy(&global);
        return error;
    }

    *global_vec = global;
    *local_vec = local;
    return PETSC_SUCCESS;
}

inline PetscErrorCode global_to_local(
    PetscSF section_sf,
    Vec global_vec,
    Vec local_vec) {
    if (section_sf == nullptr ||
        global_vec == nullptr ||
        local_vec == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    PetscInt nroots = 0;
    PetscInt nleaves = 0;
    PetscErrorCode error =
        PetscSFGetGraph(
            section_sf, &nroots, &nleaves, nullptr, nullptr);
    if (error != PETSC_SUCCESS) return error;
    if (nroots < 0 || nleaves < 0) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }

    PetscInt global_local_size = 0;
    PetscInt local_size = 0;
    error = VecGetLocalSize(global_vec, &global_local_size);
    if (error != PETSC_SUCCESS) return error;
    error = VecGetLocalSize(local_vec, &local_size);
    if (error != PETSC_SUCCESS) return error;
    if (global_local_size != nroots || local_size != nleaves) {
        return PETSC_ERR_ARG_SIZ;
    }

    const PetscScalar* roots = nullptr;
    PetscScalar* leaves = nullptr;
    error = VecGetArrayRead(global_vec, &roots);
    if (error != PETSC_SUCCESS) return error;
    error = VecGetArray(local_vec, &leaves);
    if (error != PETSC_SUCCESS) {
        VecRestoreArrayRead(global_vec, &roots);
        return error;
    }

    PetscErrorCode communication_error =
        PetscSFBcastBegin(
            section_sf,
            MPIU_SCALAR,
            roots,
            leaves,
            MPI_REPLACE);
    if (communication_error == PETSC_SUCCESS) {
        communication_error =
            PetscSFBcastEnd(
                section_sf,
                MPIU_SCALAR,
                roots,
                leaves,
                MPI_REPLACE);
    }

    const PetscErrorCode local_restore_error =
        VecRestoreArray(local_vec, &leaves);
    const PetscErrorCode global_restore_error =
        VecRestoreArrayRead(global_vec, &roots);

    if (communication_error != PETSC_SUCCESS) {
        return communication_error;
    }
    if (local_restore_error != PETSC_SUCCESS) {
        return local_restore_error;
    }
    return global_restore_error;
}

inline PetscErrorCode local_to_global_add(
    PetscSF section_sf,
    Vec local_vec,
    Vec global_vec) {
    if (section_sf == nullptr ||
        local_vec == nullptr ||
        global_vec == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }

    PetscInt nroots = 0;
    PetscInt nleaves = 0;
    PetscErrorCode error =
        PetscSFGetGraph(
            section_sf, &nroots, &nleaves, nullptr, nullptr);
    if (error != PETSC_SUCCESS) return error;
    if (nroots < 0 || nleaves < 0) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }

    PetscInt global_local_size = 0;
    PetscInt local_size = 0;
    error = VecGetLocalSize(global_vec, &global_local_size);
    if (error != PETSC_SUCCESS) return error;
    error = VecGetLocalSize(local_vec, &local_size);
    if (error != PETSC_SUCCESS) return error;
    if (global_local_size != nroots || local_size != nleaves) {
        return PETSC_ERR_ARG_SIZ;
    }

    Vec accumulated = nullptr;
    error = VecDuplicate(global_vec, &accumulated);
    if (error != PETSC_SUCCESS) return error;
    error = VecSet(accumulated, PetscScalar{0.0});
    if (error != PETSC_SUCCESS) {
        VecDestroy(&accumulated);
        return error;
    }

    const PetscScalar* leaves = nullptr;
    PetscScalar* roots = nullptr;
    error = VecGetArrayRead(local_vec, &leaves);
    if (error != PETSC_SUCCESS) {
        VecDestroy(&accumulated);
        return error;
    }
    error = VecGetArray(accumulated, &roots);
    if (error != PETSC_SUCCESS) {
        VecRestoreArrayRead(local_vec, &leaves);
        VecDestroy(&accumulated);
        return error;
    }

    PetscErrorCode communication_error =
        PetscSFReduceBegin(
            section_sf,
            MPIU_SCALAR,
            leaves,
            roots,
            MPIU_SUM);
    if (communication_error == PETSC_SUCCESS) {
        communication_error =
            PetscSFReduceEnd(
                section_sf,
                MPIU_SCALAR,
                leaves,
                roots,
                MPIU_SUM);
    }

    const PetscErrorCode root_restore_error =
        VecRestoreArray(accumulated, &roots);
    const PetscErrorCode leaf_restore_error =
        VecRestoreArrayRead(local_vec, &leaves);

    if (communication_error == PETSC_SUCCESS &&
        root_restore_error == PETSC_SUCCESS &&
        leaf_restore_error == PETSC_SUCCESS) {
        communication_error =
            VecAXPY(
                global_vec,
                PetscScalar{1.0},
                accumulated);
    }

    const PetscErrorCode destroy_error =
        VecDestroy(&accumulated);
    if (communication_error != PETSC_SUCCESS) {
        return communication_error;
    }
    if (root_restore_error != PETSC_SUCCESS) {
        return root_restore_error;
    }
    if (leaf_restore_error != PETSC_SUCCESS) {
        return leaf_restore_error;
    }
    return destroy_error;
}

} // namespace mpmc::mesh_petsc

#endif // MPMC_MESH_PETSC_ADAPTER_HPP
