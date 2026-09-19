#ifndef MPMC_MESH_PETSC_ADAPTER_HPP
#define MPMC_MESH_PETSC_ADAPTER_HPP

#include <mpmc/mesh/dof_layout.hpp>
#include <mpmc/mesh/dof_numbering.hpp>
#include <mpmc/mesh/partition_snapshot.hpp>
#include <mpmc/mesh/shared_entity_plan.hpp>
#include <mpmc/mesh/topology.hpp>

#include <petscdmplex.h>
#include <petscsection.h>
#include <petscsf.h>
#include <petscvec.h>

#include <array>
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
