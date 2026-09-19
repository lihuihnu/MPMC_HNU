#ifndef MPMC_MESH_PETSC_ADAPTER_HPP
#define MPMC_MESH_PETSC_ADAPTER_HPP

#include <mpmc/mesh/corner_point_geometry_3d.hpp>
#include <mpmc/mesh/dense_field.hpp>
#include <mpmc/mesh/dof_layout.hpp>
#include <mpmc/mesh/dof_numbering.hpp>
#include <mpmc/mesh/face_boundary.hpp>
#include <mpmc/mesh/face_geometry_3d.hpp>
#include <mpmc/mesh/geometry_2d.hpp>
#include <mpmc/mesh/partition_snapshot.hpp>
#include <mpmc/mesh/shared_entity_plan.hpp>
#include <mpmc/mesh/topology.hpp>
#include <mpmc/mesh/tpfa_internal_face_transmissibility_snapshot_3d.hpp>

#include <petscdmplex.h>
#include <petscsection.h>
#include <petscsf.h>
#include <petscvec.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cmath>
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

    std::size_t face_vertex_degree = 0U;
    for (std::size_t face = 0U; face < face_count; ++face) {
        const auto local = mpmc::mesh::LocalIndex{
            static_cast<mpmc::mesh::LocalIndex::value_type>(face)};
        const auto vertices = face_vertices.adjacent(local);
        const auto cells = face_cells.adjacent(local);
        if ((vertices.size() != 2U && vertices.size() != 4U) ||
            (cells.size() != 1U && cells.size() != 2U)) {
            return PETSC_ERR_SUP;
        }
        if (face_vertex_degree == 0U) {
            face_vertex_degree = vertices.size();
        } else if (vertices.size() != face_vertex_degree) {
            return PETSC_ERR_SUP;
        }
    }

    const PetscInt dm_dimension =
        face_vertex_degree == 2U ? PetscInt{2} : PetscInt{3};
    for (std::size_t cell = 0U; cell < cell_count; ++cell) {
        const auto local = mpmc::mesh::LocalIndex{
            static_cast<mpmc::mesh::LocalIndex::value_type>(cell)};
        const std::size_t face_degree =
            cell_faces.adjacent(local).size();
        const std::size_t vertex_degree =
            cell_vertices.adjacent(local).size();
        if (dm_dimension == PetscInt{2}) {
            if ((face_degree != 3U && face_degree != 4U) ||
                vertex_degree != face_degree) {
                return PETSC_ERR_SUP;
            }
        } else if (face_degree != 6U || vertex_degree != 8U) {
            return PETSC_ERR_SUP;
        }
    }

    std::vector<std::vector<mpmc::mesh::LocalIndex>> derived_face_cells(
        face_count);
    for (std::size_t cell = 0U; cell < cell_count; ++cell) {
        const auto cell_local = mpmc::mesh::LocalIndex{
            static_cast<mpmc::mesh::LocalIndex::value_type>(cell)};
        const auto faces = cell_faces.adjacent(cell_local);
        std::vector<mpmc::mesh::LocalIndex> union_vertices;
        union_vertices.reserve(faces.size());
        for (const auto face_local : faces) {
            const std::size_t face_position =
                static_cast<std::size_t>(face_local.value());
            if (face_position >= face_count) {
                return PETSC_ERR_ARG_OUTOFRANGE;
            }
            derived_face_cells[face_position].push_back(cell_local);
            for (const auto vertex_local :
                 face_vertices.adjacent(face_local)) {
                if (std::find(
                        union_vertices.begin(),
                        union_vertices.end(),
                        vertex_local) ==
                    union_vertices.end()) {
                    union_vertices.push_back(vertex_local);
                }
            }
        }
        const auto expected_span =
            cell_vertices.adjacent(cell_local);
        std::vector<mpmc::mesh::LocalIndex> expected(
            expected_span.begin(), expected_span.end());
        std::sort(expected.begin(), expected.end());
        std::sort(
            union_vertices.begin(), union_vertices.end());
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

    error = DMSetDimension(plex, dm_dimension);
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
        PetscInt cone_size = 0;
        error = detail::checked_petsc_int_size(
            cell_faces.adjacent(
                mpmc::mesh::LocalIndex{
                    static_cast<
                        mpmc::mesh::LocalIndex::value_type>(
                            cell)})
                .size(),
            &cone_size);
        if (error == PETSC_SUCCESS) {
            error = DMPlexSetConeSize(
                plex, point, cone_size);
        }
        if (error != PETSC_SUCCESS) {
            DMDestroy(&plex);
            return error;
        }
    }
    for (std::size_t face = 0U; face < face_count; ++face) {
        const auto local = mpmc::mesh::LocalIndex{
            static_cast<mpmc::mesh::LocalIndex::value_type>(face)};
        PetscInt point = 0;
        PetscInt cone_size = 0;
        error = detail::checked_petsc_int_size(
            face_base + face, &point);
        if (error == PETSC_SUCCESS) {
            error = detail::checked_petsc_int_size(
                face_vertices.adjacent(local).size(),
                &cone_size);
        }
        if (error == PETSC_SUCCESS) {
            error = DMPlexSetConeSize(
                plex, point, cone_size);
        }
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
        std::vector<PetscInt> cone(faces.size());
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
        std::vector<PetscInt> cone(vertices.size());
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

inline PetscErrorCode attach_serial_vertex_coordinates_3d(
    DM dm,
    std::span<const mpmc::mesh::Coordinate3D> vertex_coordinates_m,
    std::span<const DMPlexPointIdentity> identities) {
    if (dm == nullptr) return PETSC_ERR_ARG_NULL;

    MPI_Comm comm =
        PetscObjectComm(reinterpret_cast<PetscObject>(dm));
    int mpi_size = -1;
    if (MPI_Comm_size(comm, &mpi_size) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    if (mpi_size != 1) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }

    PetscInt dimension = -1;
    PetscErrorCode error =
        DMGetDimension(dm, &dimension);
    if (error != PETSC_SUCCESS) return error;
    if (dimension != PetscInt{3}) {
        return PETSC_ERR_ARG_INCOMP;
    }

    PetscInt chart_start = 0;
    PetscInt chart_end = 0;
    error = DMPlexGetChart(
        dm, &chart_start, &chart_end);
    if (error != PETSC_SUCCESS) return error;

    PetscInt vertex_start = 0;
    PetscInt vertex_end = 0;
    error = DMPlexGetDepthStratum(
        dm, 0, &vertex_start, &vertex_end);
    if (error != PETSC_SUCCESS) return error;
    if (vertex_start < chart_start ||
        vertex_end < vertex_start ||
        vertex_end > chart_end) {
        return PETSC_ERR_PLIB;
    }

    const auto vertex_count =
        static_cast<std::size_t>(
            vertex_end - vertex_start);
    if (vertex_coordinates_m.size() !=
        vertex_count) {
        return PETSC_ERR_ARG_SIZ;
    }

    std::vector<std::uint8_t>
        seen_local(vertex_count, std::uint8_t{0U});
    std::vector<std::uint8_t>
        seen_point(vertex_count, std::uint8_t{0U});
    std::size_t vertex_identity_count = 0U;
    for (const auto& identity : identities) {
        if (identity.point < chart_start ||
            identity.point >= chart_end) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }
        const bool point_is_vertex =
            identity.point >= vertex_start &&
            identity.point < vertex_end;
        if (identity.kind !=
            mpmc::mesh::EntityKind::vertex) {
            if (point_is_vertex) {
                return PETSC_ERR_ARG_INCOMP;
            }
            continue;
        }
        if (!point_is_vertex) {
            return PETSC_ERR_ARG_INCOMP;
        }

        const std::size_t local =
            static_cast<std::size_t>(
                identity.local.value());
        if (local >= vertex_count) {
            return PETSC_ERR_ARG_SIZ;
        }
        const std::size_t point_slot =
            static_cast<std::size_t>(
                identity.point - vertex_start);
        if (seen_local[local] != std::uint8_t{0U} ||
            seen_point[point_slot] != std::uint8_t{0U}) {
            return PETSC_ERR_ARG_INCOMP;
        }
        seen_local[local] = std::uint8_t{1U};
        seen_point[point_slot] = std::uint8_t{1U};
        ++vertex_identity_count;

        const auto coordinate =
            vertex_coordinates_m[local];
        if (!std::isfinite(coordinate.x_m) ||
            !std::isfinite(coordinate.y_m) ||
            !std::isfinite(coordinate.z_m)) {
            return PETSC_ERR_FP;
        }
    }
    if (vertex_identity_count != vertex_count) {
        return PETSC_ERR_ARG_SIZ;
    }

    error = DMSetCoordinateDim(dm, 3);
    if (error != PETSC_SUCCESS) return error;

    PetscSection coordinate_section = nullptr;
    error = PetscSectionCreate(
        comm, &coordinate_section);
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
            3);
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
    if (vertex_count >
        static_cast<std::size_t>(
            std::numeric_limits<PetscInt>::max()) /
            std::size_t{3U}) {
        PetscSectionDestroy(&coordinate_section);
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    const PetscInt expected_storage =
        static_cast<PetscInt>(
            vertex_count * std::size_t{3U});
    if (coordinate_storage != expected_storage) {
        PetscSectionDestroy(&coordinate_section);
        return PETSC_ERR_PLIB;
    }

    Vec coordinates = nullptr;
    error = VecCreateSeq(
        comm,
        coordinate_storage,
        &coordinates);
    if (error != PETSC_SUCCESS) {
        PetscSectionDestroy(&coordinate_section);
        return error;
    }
    error = VecSetBlockSize(coordinates, 3);
    if (error != PETSC_SUCCESS) {
        VecDestroy(&coordinates);
        PetscSectionDestroy(&coordinate_section);
        return error;
    }

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
            offset + 2 >= coordinate_storage) {
            VecRestoreArray(coordinates, &values);
            VecDestroy(&coordinates);
            PetscSectionDestroy(&coordinate_section);
            return error != PETSC_SUCCESS
                       ? error
                       : PETSC_ERR_PLIB;
        }

        const auto coordinate =
            vertex_coordinates_m[
                static_cast<std::size_t>(
                    identity.local.value())];
        values[static_cast<std::size_t>(offset)] =
            static_cast<PetscScalar>(
                coordinate.x_m);
        values[static_cast<std::size_t>(offset + 1)] =
            static_cast<PetscScalar>(
                coordinate.y_m);
        values[static_cast<std::size_t>(offset + 2)] =
            static_cast<PetscScalar>(
                coordinate.z_m);
    }
    error = VecRestoreArray(coordinates, &values);
    if (error != PETSC_SUCCESS) {
        VecDestroy(&coordinates);
        PetscSectionDestroy(&coordinate_section);
        return error;
    }

    error = DMSetCoordinateSection(
        dm, 3, coordinate_section);
    if (error != PETSC_SUCCESS) {
        VecDestroy(&coordinates);
        PetscSectionDestroy(&coordinate_section);
        return error;
    }
    error = DMSetCoordinatesLocal(
        dm, coordinates);
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

inline PetscErrorCode attach_root_vertex_coordinates_3d(
    DM dm,
    PetscMPIInt root_rank,
    std::span<const mpmc::mesh::Coordinate3D> root_vertex_coordinates_m,
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

    PetscInt dimension = -1;
    PetscErrorCode error =
        DMGetDimension(dm, &dimension);
    if (error != PETSC_SUCCESS) return error;
    if (dimension != PetscInt{3}) {
        return PETSC_ERR_ARG_INCOMP;
    }

    PetscInt chart_start = 0;
    PetscInt chart_end = 0;
    error = DMPlexGetChart(
        dm, &chart_start, &chart_end);
    if (error != PETSC_SUCCESS) return error;

    PetscErrorCode local_validation = PETSC_SUCCESS;
    if (mpi_rank == root_rank) {
        PetscInt vertex_start = 0;
        PetscInt vertex_end = 0;
        local_validation = DMPlexGetDepthStratum(
            dm, 0, &vertex_start, &vertex_end);
        if (local_validation == PETSC_SUCCESS &&
            (vertex_start < chart_start ||
             vertex_end < vertex_start ||
             vertex_end > chart_end)) {
            local_validation = PETSC_ERR_PLIB;
        }

        std::size_t vertex_identity_count = 0U;
        if (local_validation == PETSC_SUCCESS) {
            const auto vertex_count =
                static_cast<std::size_t>(
                    vertex_end - vertex_start);
            if (root_vertex_coordinates_m.size() !=
                vertex_count) {
                local_validation = PETSC_ERR_ARG_SIZ;
            }

            std::vector<std::uint8_t>
                seen_local(vertex_count, std::uint8_t{0U});
            std::vector<std::uint8_t>
                seen_point(vertex_count, std::uint8_t{0U});
            for (const auto& identity : identities) {
                if (identity.point < chart_start ||
                    identity.point >= chart_end) {
                    local_validation =
                        PETSC_ERR_ARG_OUTOFRANGE;
                    break;
                }

                const bool point_is_vertex =
                    identity.point >= vertex_start &&
                    identity.point < vertex_end;
                if (identity.kind !=
                    mpmc::mesh::EntityKind::vertex) {
                    if (point_is_vertex) {
                        local_validation =
                            PETSC_ERR_ARG_INCOMP;
                        break;
                    }
                    continue;
                }
                if (!point_is_vertex) {
                    local_validation =
                        PETSC_ERR_ARG_INCOMP;
                    break;
                }

                const std::size_t local =
                    static_cast<std::size_t>(
                        identity.local.value());
                const std::size_t point_slot =
                    static_cast<std::size_t>(
                        identity.point - vertex_start);
                if (local >= vertex_count ||
                    seen_local[local] !=
                        std::uint8_t{0U} ||
                    seen_point[point_slot] !=
                        std::uint8_t{0U}) {
                    local_validation =
                        PETSC_ERR_ARG_INCOMP;
                    break;
                }

                const auto coordinate =
                    root_vertex_coordinates_m[local];
                if (!std::isfinite(coordinate.x_m) ||
                    !std::isfinite(coordinate.y_m) ||
                    !std::isfinite(coordinate.z_m)) {
                    local_validation = PETSC_ERR_FP;
                    break;
                }

                seen_local[local] = std::uint8_t{1U};
                seen_point[point_slot] =
                    std::uint8_t{1U};
                ++vertex_identity_count;
            }
            if (local_validation == PETSC_SUCCESS &&
                vertex_identity_count != vertex_count) {
                local_validation = PETSC_ERR_ARG_SIZ;
            }
        }
    } else if (!identities.empty() ||
               !root_vertex_coordinates_m.empty()) {
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

    error = DMSetCoordinateDim(dm, 3);
    if (error != PETSC_SUCCESS) return error;

    PetscSection coordinate_section = nullptr;
    error = PetscSectionCreate(
        comm, &coordinate_section);
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
            3);
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
    if (error != PETSC_SUCCESS ||
        coordinate_storage < 0) {
        PetscSectionDestroy(&coordinate_section);
        return error != PETSC_SUCCESS
                   ? error
                   : PETSC_ERR_PLIB;
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
    error = VecSetBlockSize(coordinates, 3);
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
                offset + 2 >= coordinate_storage) {
                VecRestoreArray(coordinates, &values);
                VecDestroy(&coordinates);
                PetscSectionDestroy(&coordinate_section);
                return error != PETSC_SUCCESS
                           ? error
                           : PETSC_ERR_PLIB;
            }

            const auto coordinate =
                root_vertex_coordinates_m[
                    static_cast<std::size_t>(
                        identity.local.value())];
            values[static_cast<std::size_t>(offset)] =
                static_cast<PetscScalar>(
                    coordinate.x_m);
            values[static_cast<std::size_t>(
                offset + 1)] =
                static_cast<PetscScalar>(
                    coordinate.y_m);
            values[static_cast<std::size_t>(
                offset + 2)] =
                static_cast<PetscScalar>(
                    coordinate.z_m);
        }

        error = VecRestoreArray(
            coordinates, &values);
        if (error != PETSC_SUCCESS) {
            VecDestroy(&coordinates);
            PetscSectionDestroy(&coordinate_section);
            return error;
        }
    }

    error = DMSetCoordinateSection(
        dm, 3, coordinate_section);
    if (error != PETSC_SUCCESS) {
        VecDestroy(&coordinates);
        PetscSectionDestroy(&coordinate_section);
        return error;
    }
    error = DMSetCoordinatesLocal(
        dm, coordinates);
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
    PetscInt root_dimension = 0;
    if (mpi_rank == root_rank) {
        error = DMPlexGetChart(
            serial_dm, &chart_start, &chart_end);
        if (error == PETSC_SUCCESS) {
            error = DMGetDimension(
                serial_dm, &root_dimension);
        }
        if (error != PETSC_SUCCESS) {
            DMDestroy(&serial_dm);
            return error;
        }
        if (chart_start != 0 ||
            root_dimension <= 0) {
            DMDestroy(&serial_dm);
            return PETSC_ERR_ARG_INCOMP;
        }
    }

    if (MPI_Bcast(
            &root_dimension,
            1,
            MPIU_INT,
            root_rank,
            comm) != MPI_SUCCESS) {
        if (serial_dm != nullptr) DMDestroy(&serial_dm);
        return PETSC_ERR_MPI;
    }

    DM rooted_dm = nullptr;
    error = DMPlexCreate(comm, &rooted_dm);
    if (error != PETSC_SUCCESS) {
        if (serial_dm != nullptr) DMDestroy(&serial_dm);
        return error;
    }

    error = DMSetDimension(rooted_dm, root_dimension);
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

struct StableOwnerFaceGeometry3D {
    std::vector<mpmc::mesh::Coordinate3D>
        face_centroids_m;
    std::vector<double>
        face_areas_m2;
    std::vector<mpmc::mesh::GlobalEntityId>
        face_owner_global_ids;
    std::vector<mpmc::mesh::UnitVector3D>
        face_owner_unit_normals;

    [[nodiscard]] std::size_t
    face_count() const noexcept {
        return face_centroids_m.size();
    }
};

inline PetscErrorCode make_stable_owner_face_geometry_3d(
    const mpmc::mesh::FaceGeometry3D& source_geometry,
    std::span<const DMPlexPointIdentity> source_identities,
    StableOwnerFaceGeometry3D* output) {
    if (output == nullptr) return PETSC_ERR_ARG_NULL;
    *output = StableOwnerFaceGeometry3D{};

    const std::size_t face_count =
        source_geometry.face_count();
    if (source_geometry.face_centroids_m().size() != face_count ||
        source_geometry.face_areas_m2().size() != face_count ||
        source_geometry.face_owners().size() != face_count ||
        source_geometry.face_owner_unit_normals().size() != face_count) {
        return PETSC_ERR_ARG_SIZ;
    }

    std::vector<mpmc::mesh::GlobalEntityId>
        cell_global_ids(
            source_geometry.cell_count(),
            mpmc::mesh::GlobalEntityId{0U});
    std::vector<std::uint8_t> cell_seen(
        source_geometry.cell_count(),
        std::uint8_t{0U});
    std::vector<std::uint8_t> face_seen(
        face_count,
        std::uint8_t{0U});

    for (const auto& identity : source_identities) {
        if (identity.kind ==
            mpmc::mesh::EntityKind::cell) {
            const std::size_t local =
                static_cast<std::size_t>(
                    identity.local.value());
            if (local >= cell_global_ids.size() ||
                cell_seen[local] != std::uint8_t{0U}) {
                return PETSC_ERR_ARG_INCOMP;
            }
            cell_seen[local] = std::uint8_t{1U};
            cell_global_ids[local] = identity.global;
        } else if (identity.kind ==
                   mpmc::mesh::EntityKind::face) {
            const std::size_t local =
                static_cast<std::size_t>(
                    identity.local.value());
            if (local >= face_count ||
                face_seen[local] != std::uint8_t{0U}) {
                return PETSC_ERR_ARG_INCOMP;
            }
            face_seen[local] = std::uint8_t{1U};
        }
    }

    if (std::find(
            cell_seen.begin(),
            cell_seen.end(),
            std::uint8_t{0U}) != cell_seen.end() ||
        std::find(
            face_seen.begin(),
            face_seen.end(),
            std::uint8_t{0U}) != face_seen.end()) {
        return PETSC_ERR_ARG_SIZ;
    }

    output->face_centroids_m.assign(
        source_geometry.face_centroids_m().begin(),
        source_geometry.face_centroids_m().end());
    output->face_areas_m2.assign(
        source_geometry.face_areas_m2().begin(),
        source_geometry.face_areas_m2().end());
    output->face_owner_unit_normals.assign(
        source_geometry.face_owner_unit_normals().begin(),
        source_geometry.face_owner_unit_normals().end());
    output->face_owner_global_ids.reserve(face_count);

    for (std::size_t face = 0U;
         face < face_count;
         ++face) {
        const auto owner =
            source_geometry.face_owner(
                mpmc::mesh::LocalIndex{
                    static_cast<
                        mpmc::mesh::LocalIndex::value_type>(
                            face)});
        const std::size_t owner_local =
            static_cast<std::size_t>(
                owner.value());
        if (owner_local >= cell_global_ids.size()) {
            *output = StableOwnerFaceGeometry3D{};
            return PETSC_ERR_ARG_OUTOFRANGE;
        }
        output->face_owner_global_ids.push_back(
            cell_global_ids[owner_local]);
    }

    return PETSC_SUCCESS;
}

inline PetscErrorCode migrate_stable_owner_face_geometry_3d(
    DM source_dm,
    PetscSF migration_sf,
    const StableOwnerFaceGeometry3D& source_geometry,
    std::span<const DMPlexPointIdentity> source_identities,
    DM target_dm,
    std::span<const DMPlexPointIdentity> target_identities,
    StableOwnerFaceGeometry3D* target_geometry) {
    if (source_dm == nullptr ||
        migration_sf == nullptr ||
        target_dm == nullptr ||
        target_geometry == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    *target_geometry = StableOwnerFaceGeometry3D{};

    const std::size_t source_face_count =
        source_geometry.face_count();
    if (source_geometry.face_areas_m2.size() != source_face_count ||
        source_geometry.face_owner_global_ids.size() != source_face_count ||
        source_geometry.face_owner_unit_normals.size() != source_face_count) {
        return PETSC_ERR_ARG_SIZ;
    }

    PetscInt source_start = 0;
    PetscInt source_end = 0;
    PetscErrorCode error =
        DMPlexGetChart(
            source_dm, &source_start, &source_end);
    if (error != PETSC_SUCCESS) return error;

    PetscSection metric_source_section = nullptr;
    PetscSection metric_target_section = nullptr;
    PetscSection owner_source_section = nullptr;
    PetscSection owner_target_section = nullptr;

    auto destroy_sections = [&]() {
        PetscSectionDestroy(&owner_target_section);
        PetscSectionDestroy(&owner_source_section);
        PetscSectionDestroy(&metric_target_section);
        PetscSectionDestroy(&metric_source_section);
    };

    error = PetscSectionCreate(
        PetscObjectComm(
            reinterpret_cast<PetscObject>(source_dm)),
        &metric_source_section);
    if (error != PETSC_SUCCESS) return error;
    error = PetscSectionCreate(
        PetscObjectComm(
            reinterpret_cast<PetscObject>(target_dm)),
        &metric_target_section);
    if (error != PETSC_SUCCESS) {
        destroy_sections();
        return error;
    }
    error = PetscSectionCreate(
        PetscObjectComm(
            reinterpret_cast<PetscObject>(source_dm)),
        &owner_source_section);
    if (error != PETSC_SUCCESS) {
        destroy_sections();
        return error;
    }
    error = PetscSectionCreate(
        PetscObjectComm(
            reinterpret_cast<PetscObject>(target_dm)),
        &owner_target_section);
    if (error != PETSC_SUCCESS) {
        destroy_sections();
        return error;
    }

    for (PetscSection section :
         {metric_source_section,
          owner_source_section}) {
        error = PetscSectionSetChart(
            section, source_start, source_end);
        if (error != PETSC_SUCCESS) {
            destroy_sections();
            return error;
        }
    }

    std::vector<std::uint8_t>
        source_face_seen(
            source_face_count,
            std::uint8_t{0U});
    for (const auto& identity : source_identities) {
        if (identity.kind !=
            mpmc::mesh::EntityKind::face) {
            continue;
        }
        const std::size_t local =
            static_cast<std::size_t>(
                identity.local.value());
        if (identity.point < source_start ||
            identity.point >= source_end ||
            local >= source_face_count ||
            source_face_seen[local] != std::uint8_t{0U}) {
            destroy_sections();
            return PETSC_ERR_ARG_INCOMP;
        }
        source_face_seen[local] = std::uint8_t{1U};

        error = PetscSectionSetDof(
            metric_source_section,
            identity.point,
            7);
        if (error == PETSC_SUCCESS) {
            error = PetscSectionSetDof(
                owner_source_section,
                identity.point,
                1);
        }
        if (error != PETSC_SUCCESS) {
            destroy_sections();
            return error;
        }
    }
    if (std::find(
            source_face_seen.begin(),
            source_face_seen.end(),
            std::uint8_t{0U}) !=
        source_face_seen.end()) {
        destroy_sections();
        return PETSC_ERR_ARG_SIZ;
    }

    error = PetscSectionSetUp(metric_source_section);
    if (error == PETSC_SUCCESS) {
        error = PetscSectionSetUp(owner_source_section);
    }
    if (error != PETSC_SUCCESS) {
        destroy_sections();
        return error;
    }

    PetscInt metric_source_storage = 0;
    PetscInt owner_source_storage = 0;
    error = PetscSectionGetStorageSize(
        metric_source_section,
        &metric_source_storage);
    if (error == PETSC_SUCCESS) {
        error = PetscSectionGetStorageSize(
            owner_source_section,
            &owner_source_storage);
    }
    if (error != PETSC_SUCCESS ||
        metric_source_storage < 0 ||
        owner_source_storage < 0) {
        destroy_sections();
        return error != PETSC_SUCCESS
                   ? error
                   : PETSC_ERR_PLIB;
    }

    std::vector<double> metric_source_data(
        static_cast<std::size_t>(
            metric_source_storage),
        0.0);
    std::vector<std::uint64_t> owner_source_data(
        static_cast<std::size_t>(
            owner_source_storage),
        0U);

    for (const auto& identity : source_identities) {
        if (identity.kind !=
            mpmc::mesh::EntityKind::face) {
            continue;
        }
        const std::size_t local =
            static_cast<std::size_t>(
                identity.local.value());
        PetscInt metric_offset = -1;
        PetscInt owner_offset = -1;
        error = PetscSectionGetOffset(
            metric_source_section,
            identity.point,
            &metric_offset);
        if (error == PETSC_SUCCESS) {
            error = PetscSectionGetOffset(
                owner_source_section,
                identity.point,
                &owner_offset);
        }
        if (error != PETSC_SUCCESS ||
            metric_offset < 0 ||
            metric_offset + 6 >= metric_source_storage ||
            owner_offset < 0 ||
            owner_offset >= owner_source_storage) {
            destroy_sections();
            return error != PETSC_SUCCESS
                       ? error
                       : PETSC_ERR_PLIB;
        }

        const auto centroid =
            source_geometry.face_centroids_m[local];
        const double area =
            source_geometry.face_areas_m2[local];
        const auto normal =
            source_geometry
                .face_owner_unit_normals[local];

        metric_source_data[
            static_cast<std::size_t>(
                metric_offset)] = centroid.x_m;
        metric_source_data[
            static_cast<std::size_t>(
                metric_offset + 1)] = centroid.y_m;
        metric_source_data[
            static_cast<std::size_t>(
                metric_offset + 2)] = centroid.z_m;
        metric_source_data[
            static_cast<std::size_t>(
                metric_offset + 3)] = area;
        metric_source_data[
            static_cast<std::size_t>(
                metric_offset + 4)] = normal.x;
        metric_source_data[
            static_cast<std::size_t>(
                metric_offset + 5)] = normal.y;
        metric_source_data[
            static_cast<std::size_t>(
                metric_offset + 6)] = normal.z;
        owner_source_data[
            static_cast<std::size_t>(
                owner_offset)] =
            source_geometry
                .face_owner_global_ids[local]
                .value();
    }

    void* metric_target_raw = nullptr;
    void* owner_target_raw = nullptr;
    error = DMPlexDistributeData(
        source_dm,
        migration_sf,
        metric_source_section,
        MPI_DOUBLE,
        metric_source_data.empty()
            ? nullptr
            : metric_source_data.data(),
        metric_target_section,
        &metric_target_raw);
    if (error == PETSC_SUCCESS) {
        error = DMPlexDistributeData(
            source_dm,
            migration_sf,
            owner_source_section,
            MPI_UINT64_T,
            owner_source_data.empty()
                ? nullptr
                : owner_source_data.data(),
            owner_target_section,
            &owner_target_raw);
    }
    if (error != PETSC_SUCCESS) {
        PetscFree(metric_target_raw);
        PetscFree(owner_target_raw);
        destroy_sections();
        return error;
    }

    PetscInt face_start = -1;
    PetscInt face_end = -1;
    error = DMPlexGetHeightStratum(
        target_dm, 1, &face_start, &face_end);
    if (error != PETSC_SUCCESS ||
        face_start < 0 ||
        face_end < face_start) {
        PetscFree(metric_target_raw);
        PetscFree(owner_target_raw);
        destroy_sections();
        return error != PETSC_SUCCESS
                   ? error
                   : PETSC_ERR_PLIB;
    }
    const std::size_t target_face_count =
        static_cast<std::size_t>(
            face_end - face_start);

    target_geometry->face_centroids_m.resize(
        target_face_count);
    target_geometry->face_areas_m2.resize(
        target_face_count);
    target_geometry->face_owner_global_ids.assign(
        target_face_count,
        mpmc::mesh::GlobalEntityId{0U});
    target_geometry->face_owner_unit_normals.resize(
        target_face_count);
    std::vector<std::uint8_t> target_seen(
        target_face_count,
        std::uint8_t{0U});

    auto* metric_target_data =
        static_cast<double*>(
            metric_target_raw);
    auto* owner_target_data =
        static_cast<std::uint64_t*>(
            owner_target_raw);

    for (const auto& identity : target_identities) {
        if (identity.kind !=
            mpmc::mesh::EntityKind::face) {
            continue;
        }
        const std::size_t local =
            static_cast<std::size_t>(
                identity.local.value());
        if (local >= target_face_count ||
            target_seen[local] != std::uint8_t{0U}) {
            *target_geometry =
                StableOwnerFaceGeometry3D{};
            PetscFree(metric_target_raw);
            PetscFree(owner_target_raw);
            destroy_sections();
            return PETSC_ERR_ARG_INCOMP;
        }

        PetscInt metric_dof = 0;
        PetscInt owner_dof = 0;
        PetscInt metric_offset = -1;
        PetscInt owner_offset = -1;
        error = PetscSectionGetDof(
            metric_target_section,
            identity.point,
            &metric_dof);
        if (error == PETSC_SUCCESS) {
            error = PetscSectionGetDof(
                owner_target_section,
                identity.point,
                &owner_dof);
        }
        if (error == PETSC_SUCCESS) {
            error = PetscSectionGetOffset(
                metric_target_section,
                identity.point,
                &metric_offset);
        }
        if (error == PETSC_SUCCESS) {
            error = PetscSectionGetOffset(
                owner_target_section,
                identity.point,
                &owner_offset);
        }
        if (error != PETSC_SUCCESS ||
            metric_dof != 7 ||
            owner_dof != 1 ||
            metric_offset < 0 ||
            owner_offset < 0) {
            *target_geometry =
                StableOwnerFaceGeometry3D{};
            PetscFree(metric_target_raw);
            PetscFree(owner_target_raw);
            destroy_sections();
            return error != PETSC_SUCCESS
                       ? error
                       : PETSC_ERR_PLIB;
        }

        const auto metric_base =
            static_cast<std::size_t>(
                metric_offset);
        target_geometry
            ->face_centroids_m[local] =
            mpmc::mesh::Coordinate3D{
                metric_target_data[
                    metric_base],
                metric_target_data[
                    metric_base + 1U],
                metric_target_data[
                    metric_base + 2U]};
        target_geometry
            ->face_areas_m2[local] =
            metric_target_data[
                metric_base + 3U];
        target_geometry
            ->face_owner_unit_normals[local] =
            mpmc::mesh::UnitVector3D{
                metric_target_data[
                    metric_base + 4U],
                metric_target_data[
                    metric_base + 5U],
                metric_target_data[
                    metric_base + 6U]};
        target_geometry
            ->face_owner_global_ids[local] =
            mpmc::mesh::GlobalEntityId{
                owner_target_data[
                    static_cast<std::size_t>(
                        owner_offset)]};
        target_seen[local] = std::uint8_t{1U};
    }

    if (std::find(
            target_seen.begin(),
            target_seen.end(),
            std::uint8_t{0U}) !=
        target_seen.end()) {
        *target_geometry =
            StableOwnerFaceGeometry3D{};
        PetscFree(metric_target_raw);
        PetscFree(owner_target_raw);
        destroy_sections();
        return PETSC_ERR_ARG_INCOMP;
    }

    for (std::size_t face = 0U;
         face < target_face_count;
         ++face) {
        const auto centroid =
            target_geometry
                ->face_centroids_m[face];
        const double area =
            target_geometry
                ->face_areas_m2[face];
        const auto normal =
            target_geometry
                ->face_owner_unit_normals[face];
        const double normal_magnitude =
            std::sqrt(
                normal.x * normal.x +
                normal.y * normal.y +
                normal.z * normal.z);
        if (!std::isfinite(centroid.x_m) ||
            !std::isfinite(centroid.y_m) ||
            !std::isfinite(centroid.z_m) ||
            !std::isfinite(area) ||
            area <= 0.0 ||
            !std::isfinite(normal_magnitude) ||
            std::abs(normal_magnitude - 1.0) >
                128.0 *
                std::numeric_limits<double>::epsilon()) {
            *target_geometry =
                StableOwnerFaceGeometry3D{};
            PetscFree(metric_target_raw);
            PetscFree(owner_target_raw);
            destroy_sections();
            return PETSC_ERR_FP;
        }
    }

    const PetscErrorCode metric_free_error =
        PetscFree(metric_target_raw);
    const PetscErrorCode owner_free_error =
        PetscFree(owner_target_raw);
    destroy_sections();
    if (metric_free_error != PETSC_SUCCESS) {
        return metric_free_error;
    }
    return owner_free_error;
}

inline PetscErrorCode materialize_face_geometry_3d(
    const StableOwnerFaceGeometry3D& source_geometry,
    std::span<const DMPlexPointIdentity> target_identities,
    std::optional<mpmc::mesh::FaceGeometry3D>* target_geometry) {
    if (target_geometry == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    target_geometry->reset();

    const std::size_t face_count =
        source_geometry.face_count();
    if (source_geometry.face_areas_m2.size() != face_count ||
        source_geometry.face_owner_global_ids.size() != face_count ||
        source_geometry.face_owner_unit_normals.size() != face_count) {
        return PETSC_ERR_ARG_SIZ;
    }

    std::size_t cell_count = 0U;
    for (const auto& identity : target_identities) {
        if (identity.kind ==
            mpmc::mesh::EntityKind::cell) {
            cell_count = std::max(
                cell_count,
                static_cast<std::size_t>(
                    identity.local.value()) +
                    std::size_t{1U});
        }
    }

    std::vector<mpmc::mesh::LocalIndex>
        owners;
    owners.reserve(face_count);
    for (const auto owner_global :
         source_geometry.face_owner_global_ids) {
        const auto found =
            std::find_if(
                target_identities.begin(),
                target_identities.end(),
                [owner_global](
                    const DMPlexPointIdentity& identity) {
                    return identity.kind ==
                               mpmc::mesh::EntityKind::cell &&
                           identity.global ==
                               owner_global;
                });
        if (found == target_identities.end()) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }
        owners.push_back(found->local);
    }

    try {
        target_geometry->emplace(
            cell_count,
            source_geometry.face_centroids_m,
            source_geometry.face_areas_m2,
            std::move(owners),
            source_geometry.face_owner_unit_normals);
    } catch (...) {
        target_geometry->reset();
        return PETSC_ERR_ARG_INCOMP;
    }
    return PETSC_SUCCESS;
}


struct StableFaceGatedTpfaSnapshot3D {
    mpmc::mesh::TransmissibilityGeometryAdmissibilityPolicy3D
        geometry_policy{0.0};
    mpmc::mesh::KOrthogonalityAdmissibilityPolicy3D
        k_policy{0.0};
    std::vector<mpmc::mesh::GlobalEntityId>
        face_global_ids;
    std::vector<
        mpmc::mesh::TpfaInternalFaceTransmissibilityDisposition3D>
        dispositions;
    std::vector<std::optional<double>>
        materialized_face_transmissibilities_m3;

    [[nodiscard]] std::size_t
    entry_count() const noexcept {
        return face_global_ids.size();
    }
};

inline PetscErrorCode make_root_stable_face_gated_tpfa_snapshot_3d(
    MPI_Comm comm,
    PetscMPIInt root_rank,
    const mpmc::mesh::TpfaInternalFaceTransmissibilitySnapshot3D*
        root_snapshot,
    std::span<const DMPlexPointIdentity> root_identities,
    StableFaceGatedTpfaSnapshot3D* output) {
    if (output == nullptr) return PETSC_ERR_ARG_NULL;
    *output = StableFaceGatedTpfaSnapshot3D{};

    int mpi_rank = -1;
    int mpi_size = -1;
    if (MPI_Comm_rank(comm, &mpi_rank) != MPI_SUCCESS ||
        MPI_Comm_size(comm, &mpi_size) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    if (root_rank < 0 || root_rank >= mpi_size) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    std::array<double, 2> policies{0.0, 0.0};
    PetscErrorCode local_error = PETSC_SUCCESS;
    if (mpi_rank == root_rank) {
        if (root_snapshot == nullptr) {
            local_error = PETSC_ERR_ARG_NULL;
        } else {
            policies[0] =
                root_snapshot
                    ->geometry_policy()
                    .max_direct_normal_projection_angle_rad;
            policies[1] =
                root_snapshot
                    ->k_policy()
                    .max_half_face_co_normal_angle_rad;
            if (!std::isfinite(policies[0]) ||
                !std::isfinite(policies[1])) {
                local_error = PETSC_ERR_FP;
            }
        }
    } else if (root_snapshot != nullptr ||
               !root_identities.empty()) {
        local_error = PETSC_ERR_ARG_INCOMP;
    }

    int local_code =
        static_cast<int>(local_error);
    if (MPI_Bcast(
            &local_code,
            1,
            MPI_INT,
            root_rank,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    if (local_code !=
        static_cast<int>(PETSC_SUCCESS)) {
        return static_cast<PetscErrorCode>(
            local_code);
    }
    if (MPI_Bcast(
            policies.data(),
            static_cast<int>(policies.size()),
            MPI_DOUBLE,
            root_rank,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    output->geometry_policy =
        mpmc::mesh::
            TransmissibilityGeometryAdmissibilityPolicy3D{
                policies[0]};
    output->k_policy =
        mpmc::mesh::
            KOrthogonalityAdmissibilityPolicy3D{
                policies[1]};

    if (mpi_rank != root_rank) {
        return PETSC_SUCCESS;
    }

    output->face_global_ids.reserve(
        root_snapshot->internal_face_count());
    output->dispositions.reserve(
        root_snapshot->internal_face_count());
    output->materialized_face_transmissibilities_m3.reserve(
        root_snapshot->internal_face_count());

    for (const auto& entry :
         root_snapshot->entries()) {
        const auto identity =
            std::find_if(
                root_identities.begin(),
                root_identities.end(),
                [&entry](
                    const DMPlexPointIdentity& candidate) {
                    return candidate.kind ==
                               mpmc::mesh::EntityKind::face &&
                           candidate.local ==
                               entry.face;
                });
        if (identity == root_identities.end()) {
            *output = StableFaceGatedTpfaSnapshot3D{};
            return PETSC_ERR_ARG_INCOMP;
        }

        output->face_global_ids.push_back(
            identity->global);
        output->dispositions.push_back(
            entry.disposition);
        if (entry.static_transmissibility.has_value()) {
            const double value =
                entry.static_transmissibility
                    ->face_transmissibility_m3;
            if (!std::isfinite(value) ||
                value <= 0.0) {
                *output =
                    StableFaceGatedTpfaSnapshot3D{};
                return PETSC_ERR_FP;
            }
            output
                ->materialized_face_transmissibilities_m3
                .push_back(value);
        } else {
            output
                ->materialized_face_transmissibilities_m3
                .push_back(std::nullopt);
        }
    }

    return PETSC_SUCCESS;
}

inline PetscErrorCode validate_stable_face_gated_tpfa_snapshot_3d(
    const StableFaceGatedTpfaSnapshot3D& snapshot) {
    const std::size_t count =
        snapshot.entry_count();
    if (snapshot.dispositions.size() != count ||
        snapshot
                .materialized_face_transmissibilities_m3
                .size() !=
            count) {
        return PETSC_ERR_ARG_SIZ;
    }
    const double half_pi =
        0.5 * std::acos(-1.0);
    if (!std::isfinite(
            snapshot.geometry_policy
                .max_direct_normal_projection_angle_rad) ||
        snapshot.geometry_policy
                .max_direct_normal_projection_angle_rad <
            0.0 ||
        snapshot.geometry_policy
                .max_direct_normal_projection_angle_rad >=
            half_pi ||
        !std::isfinite(
            snapshot.k_policy
                .max_half_face_co_normal_angle_rad) ||
        snapshot.k_policy
                .max_half_face_co_normal_angle_rad <
            0.0 ||
        snapshot.k_policy
                .max_half_face_co_normal_angle_rad >=
            half_pi) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    for (std::size_t i = 0U; i < count; ++i) {
        for (std::size_t j = 0U; j < i; ++j) {
            if (snapshot.face_global_ids[i] ==
                snapshot.face_global_ids[j]) {
                return PETSC_ERR_ARG_INCOMP;
            }
        }

        const bool has_value =
            snapshot
                .materialized_face_transmissibilities_m3[i]
                .has_value();
        switch (snapshot.dispositions[i]) {
        case mpmc::mesh::
            TpfaInternalFaceTransmissibilityDisposition3D::
                materialized:
            if (!has_value ||
                !std::isfinite(
                    *snapshot
                         .materialized_face_transmissibilities_m3[i]) ||
                *snapshot
                     .materialized_face_transmissibilities_m3[i] <=
                    0.0) {
                return PETSC_ERR_ARG_INCOMP;
            }
            break;
        case mpmc::mesh::
            TpfaInternalFaceTransmissibilityDisposition3D::
                blocked_geometry_non_orthogonal:
        case mpmc::mesh::
            TpfaInternalFaceTransmissibilityDisposition3D::
                blocked_k_non_orthogonal:
        case mpmc::mesh::
            TpfaInternalFaceTransmissibilityDisposition3D::
                blocked_geometry_and_k_non_orthogonal:
        case mpmc::mesh::
            TpfaInternalFaceTransmissibilityDisposition3D::
                blocked_degenerate_permeability_direction:
            if (has_value) {
                return PETSC_ERR_ARG_INCOMP;
            }
            break;
        default:
            return PETSC_ERR_ARG_INCOMP;
        }
    }
    return PETSC_SUCCESS;
}

inline PetscErrorCode migrate_stable_face_gated_tpfa_snapshot_3d(
    DM source_dm,
    PetscSF migration_sf,
    const StableFaceGatedTpfaSnapshot3D& source_snapshot,
    std::span<const DMPlexPointIdentity> source_identities,
    DM target_dm,
    std::span<const DMPlexPointIdentity> target_identities,
    StableFaceGatedTpfaSnapshot3D* target_snapshot) {
    if (source_dm == nullptr ||
        migration_sf == nullptr ||
        target_dm == nullptr ||
        target_snapshot == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    *target_snapshot = StableFaceGatedTpfaSnapshot3D{};

    PetscErrorCode error =
        validate_stable_face_gated_tpfa_snapshot_3d(
            source_snapshot);
    if (error != PETSC_SUCCESS) return error;

    MPI_Comm comm =
        PetscObjectComm(
            reinterpret_cast<PetscObject>(source_dm));
    std::array<double, 2> local_policy{
        source_snapshot.geometry_policy
            .max_direct_normal_projection_angle_rad,
        source_snapshot.k_policy
            .max_half_face_co_normal_angle_rad};
    std::array<double, 2> minimum_policy{};
    std::array<double, 2> maximum_policy{};
    if (MPI_Allreduce(
            local_policy.data(),
            minimum_policy.data(),
            static_cast<int>(local_policy.size()),
            MPI_DOUBLE,
            MPI_MIN,
            comm) != MPI_SUCCESS ||
        MPI_Allreduce(
            local_policy.data(),
            maximum_policy.data(),
            static_cast<int>(local_policy.size()),
            MPI_DOUBLE,
            MPI_MAX,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    for (std::size_t i = 0U;
         i < local_policy.size();
         ++i) {
        if (minimum_policy[i] !=
            maximum_policy[i]) {
            return PETSC_ERR_ARG_INCOMP;
        }
    }

    PetscInt source_start = 0;
    PetscInt source_end = 0;
    error = DMPlexGetChart(
        source_dm,
        &source_start,
        &source_end);
    if (error != PETSC_SUCCESS) return error;

    PetscSection disposition_source_section = nullptr;
    PetscSection disposition_target_section = nullptr;
    PetscSection value_source_section = nullptr;
    PetscSection value_target_section = nullptr;
    auto destroy_sections = [&]() {
        PetscSectionDestroy(&value_target_section);
        PetscSectionDestroy(&value_source_section);
        PetscSectionDestroy(&disposition_target_section);
        PetscSectionDestroy(&disposition_source_section);
    };

    error = PetscSectionCreate(
        comm,
        &disposition_source_section);
    if (error != PETSC_SUCCESS) return error;
    error = PetscSectionCreate(
        PetscObjectComm(
            reinterpret_cast<PetscObject>(target_dm)),
        &disposition_target_section);
    if (error != PETSC_SUCCESS) {
        destroy_sections();
        return error;
    }
    error = PetscSectionCreate(
        comm,
        &value_source_section);
    if (error != PETSC_SUCCESS) {
        destroy_sections();
        return error;
    }
    error = PetscSectionCreate(
        PetscObjectComm(
            reinterpret_cast<PetscObject>(target_dm)),
        &value_target_section);
    if (error != PETSC_SUCCESS) {
        destroy_sections();
        return error;
    }

    for (PetscSection section :
         {disposition_source_section,
          value_source_section}) {
        error = PetscSectionSetChart(
            section,
            source_start,
            source_end);
        if (error != PETSC_SUCCESS) {
            destroy_sections();
            return error;
        }
    }

    std::vector<PetscInt> source_points(
        source_snapshot.entry_count(),
        PetscInt{-1});
    for (std::size_t entry = 0U;
         entry < source_snapshot.entry_count();
         ++entry) {
        const auto global =
            source_snapshot.face_global_ids[entry];
        const auto identity =
            std::find_if(
                source_identities.begin(),
                source_identities.end(),
                [global](
                    const DMPlexPointIdentity& candidate) {
                    return candidate.kind ==
                               mpmc::mesh::EntityKind::face &&
                           candidate.global == global;
                });
        if (identity == source_identities.end() ||
            identity->point < source_start ||
            identity->point >= source_end) {
            destroy_sections();
            return PETSC_ERR_ARG_INCOMP;
        }
        source_points[entry] =
            identity->point;

        error = PetscSectionSetDof(
            disposition_source_section,
            identity->point,
            1);
        if (error == PETSC_SUCCESS &&
            source_snapshot
                .materialized_face_transmissibilities_m3[
                    entry]
                .has_value()) {
            error = PetscSectionSetDof(
                value_source_section,
                identity->point,
                1);
        }
        if (error != PETSC_SUCCESS) {
            destroy_sections();
            return error;
        }
    }

    error = PetscSectionSetUp(
        disposition_source_section);
    if (error == PETSC_SUCCESS) {
        error = PetscSectionSetUp(
            value_source_section);
    }
    if (error != PETSC_SUCCESS) {
        destroy_sections();
        return error;
    }

    PetscInt disposition_storage = 0;
    PetscInt value_storage = 0;
    error = PetscSectionGetStorageSize(
        disposition_source_section,
        &disposition_storage);
    if (error == PETSC_SUCCESS) {
        error = PetscSectionGetStorageSize(
            value_source_section,
            &value_storage);
    }
    if (error != PETSC_SUCCESS ||
        disposition_storage < 0 ||
        value_storage < 0) {
        destroy_sections();
        return error != PETSC_SUCCESS
                   ? error
                   : PETSC_ERR_PLIB;
    }

    std::vector<std::uint32_t> disposition_data(
        static_cast<std::size_t>(
            disposition_storage),
        0U);
    std::vector<double> value_data(
        static_cast<std::size_t>(
            value_storage),
        0.0);

    for (std::size_t entry = 0U;
         entry < source_snapshot.entry_count();
         ++entry) {
        PetscInt disposition_offset = -1;
        error = PetscSectionGetOffset(
            disposition_source_section,
            source_points[entry],
            &disposition_offset);
        if (error != PETSC_SUCCESS ||
            disposition_offset < 0 ||
            disposition_offset >=
                disposition_storage) {
            destroy_sections();
            return error != PETSC_SUCCESS
                       ? error
                       : PETSC_ERR_PLIB;
        }
        disposition_data[
            static_cast<std::size_t>(
                disposition_offset)] =
            static_cast<std::uint32_t>(
                source_snapshot.dispositions[entry]);

        if (source_snapshot
                .materialized_face_transmissibilities_m3[
                    entry]
                .has_value()) {
            PetscInt value_offset = -1;
            error = PetscSectionGetOffset(
                value_source_section,
                source_points[entry],
                &value_offset);
            if (error != PETSC_SUCCESS ||
                value_offset < 0 ||
                value_offset >= value_storage) {
                destroy_sections();
                return error != PETSC_SUCCESS
                           ? error
                           : PETSC_ERR_PLIB;
            }
            value_data[
                static_cast<std::size_t>(
                    value_offset)] =
                *source_snapshot
                     .materialized_face_transmissibilities_m3[
                         entry];
        }
    }

    void* disposition_target_raw = nullptr;
    void* value_target_raw = nullptr;
    error = DMPlexDistributeData(
        source_dm,
        migration_sf,
        disposition_source_section,
        MPI_UINT32_T,
        disposition_data.empty()
            ? nullptr
            : disposition_data.data(),
        disposition_target_section,
        &disposition_target_raw);
    if (error == PETSC_SUCCESS) {
        error = DMPlexDistributeData(
            source_dm,
            migration_sf,
            value_source_section,
            MPI_DOUBLE,
            value_data.empty()
                ? nullptr
                : value_data.data(),
            value_target_section,
            &value_target_raw);
    }
    if (error != PETSC_SUCCESS) {
        PetscFree(disposition_target_raw);
        PetscFree(value_target_raw);
        destroy_sections();
        return error;
    }

    target_snapshot->geometry_policy =
        source_snapshot.geometry_policy;
    target_snapshot->k_policy =
        source_snapshot.k_policy;

    auto* disposition_target_data =
        static_cast<std::uint32_t*>(
            disposition_target_raw);
    auto* value_target_data =
        static_cast<double*>(
            value_target_raw);

    for (const auto& identity :
         target_identities) {
        if (identity.kind !=
            mpmc::mesh::EntityKind::face) {
            continue;
        }

        PetscInt disposition_dof = 0;
        PetscInt disposition_offset = -1;
        error = PetscSectionGetDof(
            disposition_target_section,
            identity.point,
            &disposition_dof);
        if (error != PETSC_SUCCESS) {
            *target_snapshot =
                StableFaceGatedTpfaSnapshot3D{};
            PetscFree(disposition_target_raw);
            PetscFree(value_target_raw);
            destroy_sections();
            return error;
        }
        if (disposition_dof == 0) {
            continue;
        }
        if (disposition_dof != 1) {
            *target_snapshot =
                StableFaceGatedTpfaSnapshot3D{};
            PetscFree(disposition_target_raw);
            PetscFree(value_target_raw);
            destroy_sections();
            return PETSC_ERR_PLIB;
        }
        error = PetscSectionGetOffset(
            disposition_target_section,
            identity.point,
            &disposition_offset);
        if (error != PETSC_SUCCESS ||
            disposition_offset < 0) {
            *target_snapshot =
                StableFaceGatedTpfaSnapshot3D{};
            PetscFree(disposition_target_raw);
            PetscFree(value_target_raw);
            destroy_sections();
            return error != PETSC_SUCCESS
                       ? error
                       : PETSC_ERR_PLIB;
        }

        const auto raw_disposition =
            disposition_target_data[
                static_cast<std::size_t>(
                    disposition_offset)];
        mpmc::mesh::
            TpfaInternalFaceTransmissibilityDisposition3D
                disposition;
        switch (raw_disposition) {
        case static_cast<std::uint32_t>(
                 mpmc::mesh::
                     TpfaInternalFaceTransmissibilityDisposition3D::
                         materialized):
            disposition =
                mpmc::mesh::
                    TpfaInternalFaceTransmissibilityDisposition3D::
                        materialized;
            break;
        case static_cast<std::uint32_t>(
                 mpmc::mesh::
                     TpfaInternalFaceTransmissibilityDisposition3D::
                         blocked_geometry_non_orthogonal):
            disposition =
                mpmc::mesh::
                    TpfaInternalFaceTransmissibilityDisposition3D::
                        blocked_geometry_non_orthogonal;
            break;
        case static_cast<std::uint32_t>(
                 mpmc::mesh::
                     TpfaInternalFaceTransmissibilityDisposition3D::
                         blocked_k_non_orthogonal):
            disposition =
                mpmc::mesh::
                    TpfaInternalFaceTransmissibilityDisposition3D::
                        blocked_k_non_orthogonal;
            break;
        case static_cast<std::uint32_t>(
                 mpmc::mesh::
                     TpfaInternalFaceTransmissibilityDisposition3D::
                         blocked_geometry_and_k_non_orthogonal):
            disposition =
                mpmc::mesh::
                    TpfaInternalFaceTransmissibilityDisposition3D::
                        blocked_geometry_and_k_non_orthogonal;
            break;
        case static_cast<std::uint32_t>(
                 mpmc::mesh::
                     TpfaInternalFaceTransmissibilityDisposition3D::
                         blocked_degenerate_permeability_direction):
            disposition =
                mpmc::mesh::
                    TpfaInternalFaceTransmissibilityDisposition3D::
                        blocked_degenerate_permeability_direction;
            break;
        default:
            *target_snapshot =
                StableFaceGatedTpfaSnapshot3D{};
            PetscFree(disposition_target_raw);
            PetscFree(value_target_raw);
            destroy_sections();
            return PETSC_ERR_ARG_INCOMP;
        }

        PetscInt value_dof = 0;
        error = PetscSectionGetDof(
            value_target_section,
            identity.point,
            &value_dof);
        if (error != PETSC_SUCCESS) {
            *target_snapshot =
                StableFaceGatedTpfaSnapshot3D{};
            PetscFree(disposition_target_raw);
            PetscFree(value_target_raw);
            destroy_sections();
            return error;
        }

        std::optional<double> value;
        if (disposition ==
            mpmc::mesh::
                TpfaInternalFaceTransmissibilityDisposition3D::
                    materialized) {
            if (value_dof != 1) {
                *target_snapshot =
                    StableFaceGatedTpfaSnapshot3D{};
                PetscFree(disposition_target_raw);
                PetscFree(value_target_raw);
                destroy_sections();
                return PETSC_ERR_ARG_INCOMP;
            }
            PetscInt value_offset = -1;
            error = PetscSectionGetOffset(
                value_target_section,
                identity.point,
                &value_offset);
            if (error != PETSC_SUCCESS ||
                value_offset < 0) {
                *target_snapshot =
                    StableFaceGatedTpfaSnapshot3D{};
                PetscFree(disposition_target_raw);
                PetscFree(value_target_raw);
                destroy_sections();
                return error != PETSC_SUCCESS
                           ? error
                           : PETSC_ERR_PLIB;
            }
            const double transported =
                value_target_data[
                    static_cast<std::size_t>(
                        value_offset)];
            if (!std::isfinite(transported) ||
                transported <= 0.0) {
                *target_snapshot =
                    StableFaceGatedTpfaSnapshot3D{};
                PetscFree(disposition_target_raw);
                PetscFree(value_target_raw);
                destroy_sections();
                return PETSC_ERR_FP;
            }
            value = transported;
        } else if (value_dof != 0) {
            *target_snapshot =
                StableFaceGatedTpfaSnapshot3D{};
            PetscFree(disposition_target_raw);
            PetscFree(value_target_raw);
            destroy_sections();
            return PETSC_ERR_ARG_INCOMP;
        }

        target_snapshot->face_global_ids.push_back(
            identity.global);
        target_snapshot->dispositions.push_back(
            disposition);
        target_snapshot
            ->materialized_face_transmissibilities_m3
            .push_back(value);
    }

    error =
        validate_stable_face_gated_tpfa_snapshot_3d(
            *target_snapshot);
    const PetscErrorCode disposition_free_error =
        PetscFree(disposition_target_raw);
    const PetscErrorCode value_free_error =
        PetscFree(value_target_raw);
    destroy_sections();

    if (error != PETSC_SUCCESS) return error;
    if (disposition_free_error != PETSC_SUCCESS) {
        return disposition_free_error;
    }
    return value_free_error;
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
