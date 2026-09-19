#ifndef MPMC_DISCRETIZATION_PETSC_ADAPTER_HPP
#define MPMC_DISCRETIZATION_PETSC_ADAPTER_HPP

#include <mpmc/mesh/tpfa_internal_face_transmissibility_snapshot_3d.hpp>
#include <mpmc/mesh_petsc/adapter.hpp>

#include <petscdmplex.h>
#include <petscmat.h>
#include <petscsection.h>
#include <petscsf.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::discretization_petsc {

using mpmc::mesh_petsc::DMPlexPointIdentity;
using mpmc::mesh_petsc::StableOwnerFaceGeometry3D;

namespace detail {

inline PetscErrorCode checked_petsc_int_size(
    std::size_t value,
    PetscInt* output) {
    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    if (value >
        static_cast<std::size_t>(
            std::numeric_limits<PetscInt>::max())) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    *output =
        static_cast<PetscInt>(value);
    return PETSC_SUCCESS;
}

inline PetscErrorCode validate_communicator(
    MPI_Comm comm,
    mpmc::mesh::PartitionRank local_rank,
    std::uint32_t rank_count) {
    int mpi_rank = -1;
    int mpi_size = -1;
    if (MPI_Comm_rank(
            comm,
            &mpi_rank) != MPI_SUCCESS ||
        MPI_Comm_size(
            comm,
            &mpi_size) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    if (mpi_rank < 0 ||
        mpi_size < 0) {
        return PETSC_ERR_MPI;
    }
    if (static_cast<std::uint32_t>(
            mpi_rank) !=
            local_rank.value() ||
        static_cast<std::uint32_t>(
            mpi_size) !=
            rank_count) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }
    return PETSC_SUCCESS;
}

} // namespace detail

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


struct TargetLocalGatedTpfaTransmissibilityEntry3D {
    mpmc::mesh::LocalIndex face;
    mpmc::mesh::GlobalEntityId global;
    mpmc::mesh::TpfaInternalFaceTransmissibilityDisposition3D
        disposition;
    std::optional<double> transmissibility_m3;
};

/// Read-only target-local view for assembly-facing static TPFA data.
///
/// Membership comes from the stable transported internal-face snapshot, not
/// target-local DMPlex support cardinality. This is required at overlap=0,
/// where a globally internal shared face can have only one local supporting
/// cell because the remote cell is absent from the rank-local DM.
class TargetLocalGatedTpfaTransmissibilityView3D {
public:
    TargetLocalGatedTpfaTransmissibilityView3D(
        std::size_t target_face_count,
        mpmc::mesh::
            TransmissibilityGeometryAdmissibilityPolicy3D
                geometry_policy,
        mpmc::mesh::KOrthogonalityAdmissibilityPolicy3D
            k_policy,
        std::vector<
            TargetLocalGatedTpfaTransmissibilityEntry3D>
                entries)
        : target_face_count_(target_face_count),
          geometry_policy_(geometry_policy),
          k_policy_(k_policy),
          entries_(std::move(entries)),
          face_to_entry_(target_face_count) {
        validate_and_index();
    }

    TargetLocalGatedTpfaTransmissibilityView3D(
        const TargetLocalGatedTpfaTransmissibilityView3D&) =
        default;
    TargetLocalGatedTpfaTransmissibilityView3D(
        TargetLocalGatedTpfaTransmissibilityView3D&&) noexcept =
        default;
    TargetLocalGatedTpfaTransmissibilityView3D& operator=(
        const TargetLocalGatedTpfaTransmissibilityView3D&) =
        delete;
    TargetLocalGatedTpfaTransmissibilityView3D& operator=(
        TargetLocalGatedTpfaTransmissibilityView3D&&) =
        delete;
    ~TargetLocalGatedTpfaTransmissibilityView3D() = default;

    [[nodiscard]] std::size_t
    target_face_count() const noexcept {
        return target_face_count_;
    }

    [[nodiscard]] std::size_t
    internal_face_count() const noexcept {
        return entries_.size();
    }

    [[nodiscard]] std::size_t
    materialized_face_count() const noexcept {
        return materialized_face_count_;
    }

    [[nodiscard]] std::size_t
    blocked_face_count() const noexcept {
        return entries_.size() -
               materialized_face_count_;
    }

    [[nodiscard]] mpmc::mesh::
        TransmissibilityGeometryAdmissibilityPolicy3D
    geometry_policy() const noexcept {
        return geometry_policy_;
    }

    [[nodiscard]]
    mpmc::mesh::KOrthogonalityAdmissibilityPolicy3D
    k_policy() const noexcept {
        return k_policy_;
    }

    [[nodiscard]] std::span<
        const TargetLocalGatedTpfaTransmissibilityEntry3D>
    entries() const noexcept {
        return entries_;
    }

    [[nodiscard]] bool contains_internal_face(
        mpmc::mesh::LocalIndex face) const {
        const std::size_t local =
            static_cast<std::size_t>(
                face.value());
        if (local >= target_face_count_) {
            throw std::out_of_range(
                "mpmc::discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D: face index out of range");
        }
        return face_to_entry_[local].has_value();
    }

    [[nodiscard]]
    const TargetLocalGatedTpfaTransmissibilityEntry3D&
    entry(mpmc::mesh::LocalIndex face) const {
        const std::size_t local =
            static_cast<std::size_t>(
                face.value());
        if (local >= target_face_count_) {
            throw std::out_of_range(
                "mpmc::discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D: face index out of range");
        }
        const auto mapped =
            face_to_entry_[local];
        if (!mapped.has_value()) {
            throw std::invalid_argument(
                "mpmc::discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D: face is not a transported internal face");
        }
        return entries_.at(*mapped);
    }

    [[nodiscard]] mpmc::mesh::GlobalEntityId
    global_id(mpmc::mesh::LocalIndex face) const {
        return entry(face).global;
    }

    [[nodiscard]] mpmc::mesh::
        TpfaInternalFaceTransmissibilityDisposition3D
    disposition(mpmc::mesh::LocalIndex face) const {
        return entry(face).disposition;
    }

    [[nodiscard]] std::optional<double>
    optional_transmissibility_m3(
        mpmc::mesh::LocalIndex face) const {
        return entry(face).transmissibility_m3;
    }

    /// Return the assembly-usable static T_f [m3].
    ///
    /// Blocked faces deliberately reject this accessor rather than returning
    /// zero or another sentinel that could accidentally enter an assembly.
    [[nodiscard]] double transmissibility_m3(
        mpmc::mesh::LocalIndex face) const {
        const auto& local_entry =
            entry(face);
        if (!local_entry
                 .transmissibility_m3
                 .has_value()) {
            throw std::invalid_argument(
                "mpmc::discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D: blocked face has no assembly-usable transmissibility");
        }
        return *local_entry.transmissibility_m3;
    }

private:
    void validate_and_index() {
        const double half_pi =
            0.5 * std::acos(-1.0);
        if (!std::isfinite(
                geometry_policy_
                    .max_direct_normal_projection_angle_rad) ||
            geometry_policy_
                    .max_direct_normal_projection_angle_rad <
                0.0 ||
            geometry_policy_
                    .max_direct_normal_projection_angle_rad >=
                half_pi ||
            !std::isfinite(
                k_policy_
                    .max_half_face_co_normal_angle_rad) ||
            k_policy_
                    .max_half_face_co_normal_angle_rad <
                0.0 ||
            k_policy_
                    .max_half_face_co_normal_angle_rad >=
                half_pi) {
            throw std::invalid_argument(
                "mpmc::discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D: policy angles must be finite in [0, pi/2)");
        }

        materialized_face_count_ = 0U;
        for (std::size_t index = 0U;
             index < entries_.size();
             ++index) {
            const auto& local_entry =
                entries_[index];
            const std::size_t local =
                static_cast<std::size_t>(
                    local_entry.face.value());
            if (local >= target_face_count_) {
                throw std::out_of_range(
                    "mpmc::discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D: entry face index out of range");
            }
            if (face_to_entry_[local].has_value()) {
                throw std::invalid_argument(
                    "mpmc::discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D: duplicate local face entry");
            }
            for (std::size_t prior = 0U;
                 prior < index;
                 ++prior) {
                if (entries_[prior].global ==
                    local_entry.global) {
                    throw std::invalid_argument(
                        "mpmc::discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D: duplicate stable face GlobalEntityId");
                }
            }

            const bool has_value =
                local_entry
                    .transmissibility_m3
                    .has_value();
            switch (local_entry.disposition) {
            case mpmc::mesh::
                TpfaInternalFaceTransmissibilityDisposition3D::
                    materialized:
                if (!has_value ||
                    !std::isfinite(
                        *local_entry
                             .transmissibility_m3) ||
                    *local_entry
                         .transmissibility_m3 <=
                        0.0) {
                    throw std::invalid_argument(
                        "mpmc::discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D: materialized face requires finite positive T_f");
                }
                ++materialized_face_count_;
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
                    throw std::invalid_argument(
                        "mpmc::discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D: blocked face must not expose T_f");
                }
                break;
            default:
                throw std::invalid_argument(
                    "mpmc::discretization_petsc::TargetLocalGatedTpfaTransmissibilityView3D: invalid disposition");
            }

            face_to_entry_[local] = index;
        }
    }

    std::size_t target_face_count_;
    mpmc::mesh::
        TransmissibilityGeometryAdmissibilityPolicy3D
            geometry_policy_;
    mpmc::mesh::KOrthogonalityAdmissibilityPolicy3D
        k_policy_;
    std::vector<
        TargetLocalGatedTpfaTransmissibilityEntry3D>
            entries_;
    std::vector<std::optional<std::size_t>>
        face_to_entry_;
    std::size_t materialized_face_count_{0U};
};

inline PetscErrorCode
make_target_local_gated_tpfa_transmissibility_view_3d(
    DM target_dm,
    const StableFaceGatedTpfaSnapshot3D& transport,
    std::span<const DMPlexPointIdentity> target_identities,
    std::optional<
        TargetLocalGatedTpfaTransmissibilityView3D>* output) {
    if (target_dm == nullptr ||
        output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();

    PetscErrorCode error =
        validate_stable_face_gated_tpfa_snapshot_3d(
            transport);
    if (error != PETSC_SUCCESS) return error;

    PetscInt face_start = -1;
    PetscInt face_end = -1;
    error = DMPlexGetHeightStratum(
        target_dm,
        1,
        &face_start,
        &face_end);
    if (error != PETSC_SUCCESS ||
        face_start < 0 ||
        face_end < face_start) {
        return error != PETSC_SUCCESS
                   ? error
                   : PETSC_ERR_PLIB;
    }
    const std::size_t face_count =
        static_cast<std::size_t>(
            face_end - face_start);

    std::vector<std::uint8_t> identity_face_seen(
        face_count,
        std::uint8_t{0U});
    for (const auto& identity :
         target_identities) {
        if (identity.kind !=
            mpmc::mesh::EntityKind::face) {
            continue;
        }
        if (identity.point < face_start ||
            identity.point >= face_end) {
            return PETSC_ERR_ARG_INCOMP;
        }
        const std::size_t local =
            static_cast<std::size_t>(
                identity.local.value());
        const std::size_t expected_local =
            static_cast<std::size_t>(
                identity.point -
                face_start);
        if (local >= face_count ||
            local != expected_local ||
            identity_face_seen[local] !=
                std::uint8_t{0U}) {
            return PETSC_ERR_ARG_INCOMP;
        }
        identity_face_seen[local] =
            std::uint8_t{1U};
    }
    if (std::find(
            identity_face_seen.begin(),
            identity_face_seen.end(),
            std::uint8_t{0U}) !=
        identity_face_seen.end()) {
        return PETSC_ERR_ARG_SIZ;
    }

    std::vector<
        TargetLocalGatedTpfaTransmissibilityEntry3D>
            entries;
    entries.reserve(
        transport.entry_count());

    for (std::size_t entry = 0U;
         entry < transport.entry_count();
         ++entry) {
        const auto global =
            transport.face_global_ids[entry];
        const auto found =
            std::find_if(
                target_identities.begin(),
                target_identities.end(),
                [global](
                    const DMPlexPointIdentity& identity) {
                    return identity.kind ==
                               mpmc::mesh::EntityKind::face &&
                           identity.global ==
                               global;
                });
        if (found == target_identities.end()) {
            return PETSC_ERR_ARG_INCOMP;
        }

        entries.push_back(
            TargetLocalGatedTpfaTransmissibilityEntry3D{
                found->local,
                global,
                transport.dispositions[entry],
                transport
                    .materialized_face_transmissibilities_m3[
                        entry]});
    }

    try {
        output->emplace(
            face_count,
            transport.geometry_policy,
            transport.k_policy,
            std::move(entries));
    } catch (...) {
        output->reset();
        return PETSC_ERR_ARG_INCOMP;
    }

    return PETSC_SUCCESS;
}


struct AssemblyReadyInternalConnectionRow3D {
    mpmc::mesh::LocalIndex face;
    mpmc::mesh::GlobalEntityId face_global;
    mpmc::mesh::LocalIndex owner_cell;
    mpmc::mesh::GlobalEntityId owner_cell_global;
    mpmc::mesh::LocalIndex neighbour_cell;
    mpmc::mesh::GlobalEntityId neighbour_cell_global;
    double transmissibility_m3;
};

/// Immutable target-local active internal-connection table.
///
/// Rows contain only materialized gated TPFA faces. Blocked faces remain
/// available exclusively through TargetLocalGatedTpfaTransmissibilityView3D
/// and never receive an active connection row.
///
/// "owner" is the canonical geometric owner cell carried by
/// StableOwnerFaceGeometry3D; it is unrelated to MPI point ownership.
class AssemblyReadyInternalConnectionTable3D {
public:
    AssemblyReadyInternalConnectionTable3D(
        std::size_t target_face_count,
        std::size_t target_cell_count,
        std::vector<AssemblyReadyInternalConnectionRow3D>
            rows)
        : target_face_count_(target_face_count),
          target_cell_count_(target_cell_count),
          rows_(std::move(rows)),
          face_to_row_(target_face_count) {
        validate_and_index();
    }

    AssemblyReadyInternalConnectionTable3D(
        const AssemblyReadyInternalConnectionTable3D&) =
        default;
    AssemblyReadyInternalConnectionTable3D(
        AssemblyReadyInternalConnectionTable3D&&) noexcept =
        default;
    AssemblyReadyInternalConnectionTable3D& operator=(
        const AssemblyReadyInternalConnectionTable3D&) =
        delete;
    AssemblyReadyInternalConnectionTable3D& operator=(
        AssemblyReadyInternalConnectionTable3D&&) =
        delete;
    ~AssemblyReadyInternalConnectionTable3D() = default;

    [[nodiscard]] std::size_t
    target_face_count() const noexcept {
        return target_face_count_;
    }

    [[nodiscard]] std::size_t
    target_cell_count() const noexcept {
        return target_cell_count_;
    }

    [[nodiscard]] std::size_t
    row_count() const noexcept {
        return rows_.size();
    }

    [[nodiscard]] std::span<
        const AssemblyReadyInternalConnectionRow3D>
    rows() const noexcept {
        return rows_;
    }

    [[nodiscard]] bool contains_face(
        mpmc::mesh::LocalIndex face) const {
        const std::size_t local =
            static_cast<std::size_t>(
                face.value());
        if (local >= target_face_count_) {
            throw std::out_of_range(
                "mpmc::discretization_petsc::AssemblyReadyInternalConnectionTable3D: face index out of range");
        }
        return face_to_row_[local].has_value();
    }

    [[nodiscard]]
    const AssemblyReadyInternalConnectionRow3D&
    row(mpmc::mesh::LocalIndex face) const {
        const std::size_t local =
            static_cast<std::size_t>(
                face.value());
        if (local >= target_face_count_) {
            throw std::out_of_range(
                "mpmc::discretization_petsc::AssemblyReadyInternalConnectionTable3D: face index out of range");
        }
        const auto mapped =
            face_to_row_[local];
        if (!mapped.has_value()) {
            throw std::invalid_argument(
                "mpmc::discretization_petsc::AssemblyReadyInternalConnectionTable3D: face has no active materialized connection row");
        }
        return rows_.at(*mapped);
    }

    [[nodiscard]] double transmissibility_m3(
        mpmc::mesh::LocalIndex face) const {
        return row(face).transmissibility_m3;
    }

private:
    void validate_and_index() {
        for (std::size_t index = 0U;
             index < rows_.size();
             ++index) {
            const auto& connection =
                rows_[index];
            const std::size_t face =
                static_cast<std::size_t>(
                    connection.face.value());
            const std::size_t owner =
                static_cast<std::size_t>(
                    connection.owner_cell.value());
            const std::size_t neighbour =
                static_cast<std::size_t>(
                    connection.neighbour_cell.value());

            if (face >= target_face_count_) {
                throw std::out_of_range(
                    "mpmc::discretization_petsc::AssemblyReadyInternalConnectionTable3D: row face index out of range");
            }
            if (owner >= target_cell_count_ ||
                neighbour >= target_cell_count_ ||
                owner == neighbour) {
                throw std::invalid_argument(
                    "mpmc::discretization_petsc::AssemblyReadyInternalConnectionTable3D: row must reference two distinct target-local cells");
            }
            if (!std::isfinite(
                    connection.transmissibility_m3) ||
                connection.transmissibility_m3 <=
                    0.0) {
                throw std::invalid_argument(
                    "mpmc::discretization_petsc::AssemblyReadyInternalConnectionTable3D: row transmissibility must be finite and positive");
            }
            if (face_to_row_[face].has_value()) {
                throw std::invalid_argument(
                    "mpmc::discretization_petsc::AssemblyReadyInternalConnectionTable3D: duplicate target-local face row");
            }
            for (std::size_t prior = 0U;
                 prior < index;
                 ++prior) {
                if (rows_[prior].face_global ==
                    connection.face_global) {
                    throw std::invalid_argument(
                        "mpmc::discretization_petsc::AssemblyReadyInternalConnectionTable3D: duplicate stable face GlobalEntityId");
                }
            }

            face_to_row_[face] = index;
        }
    }

    std::size_t target_face_count_;
    std::size_t target_cell_count_;
    std::vector<AssemblyReadyInternalConnectionRow3D>
        rows_;
    std::vector<std::optional<std::size_t>>
        face_to_row_;
};

/// Project materialized target-local gated faces into active connection rows.
///
/// Every materialized face must have both adjacent cells present in target_dm.
/// At overlap=0 a cross-rank internal face normally has only one local support,
/// so this routine returns PETSC_ERR_ARG_WRONGSTATE rather than fabricating a
/// remote neighbour LocalIndex or silently omitting the connection.
///
/// Blocked faces are skipped unconditionally and therefore never enter the
/// active table.
inline PetscErrorCode
make_assembly_ready_internal_connection_table_3d(
    DM target_dm,
    const TargetLocalGatedTpfaTransmissibilityView3D& view,
    const StableOwnerFaceGeometry3D& stable_face_geometry,
    std::span<const DMPlexPointIdentity> target_identities,
    std::optional<AssemblyReadyInternalConnectionTable3D>*
        output) {
    if (target_dm == nullptr ||
        output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();

    const std::size_t face_count =
        view.target_face_count();
    if (stable_face_geometry.face_count() !=
            face_count ||
        stable_face_geometry.face_areas_m2.size() !=
            face_count ||
        stable_face_geometry.face_owner_global_ids.size() !=
            face_count ||
        stable_face_geometry.face_owner_unit_normals.size() !=
            face_count) {
        return PETSC_ERR_ARG_SIZ;
    }

    PetscInt cell_start = -1;
    PetscInt cell_end = -1;
    PetscInt face_start = -1;
    PetscInt face_end = -1;
    PetscErrorCode error =
        DMPlexGetHeightStratum(
            target_dm,
            0,
            &cell_start,
            &cell_end);
    if (error == PETSC_SUCCESS) {
        error = DMPlexGetHeightStratum(
            target_dm,
            1,
            &face_start,
            &face_end);
    }
    if (error != PETSC_SUCCESS ||
        cell_start < 0 ||
        cell_end < cell_start ||
        face_start < 0 ||
        face_end < face_start) {
        return error != PETSC_SUCCESS
                   ? error
                   : PETSC_ERR_PLIB;
    }

    const std::size_t cell_count =
        static_cast<std::size_t>(
            cell_end - cell_start);
    if (static_cast<std::size_t>(
            face_end - face_start) !=
        face_count) {
        return PETSC_ERR_ARG_SIZ;
    }

    auto identity_for_point =
        [&](PetscInt point)
            -> const DMPlexPointIdentity* {
        const auto found =
            std::find_if(
                target_identities.begin(),
                target_identities.end(),
                [point](
                    const DMPlexPointIdentity& identity) {
                    return identity.point ==
                           point;
                });
        return found == target_identities.end()
                   ? nullptr
                   : &*found;
    };

    std::vector<AssemblyReadyInternalConnectionRow3D>
        rows;
    rows.reserve(
        view.materialized_face_count());

    for (const auto& gated :
         view.entries()) {
        if (gated.disposition !=
            mpmc::mesh::
                TpfaInternalFaceTransmissibilityDisposition3D::
                    materialized) {
            continue;
        }
        if (!gated.transmissibility_m3.has_value() ||
            !std::isfinite(
                *gated.transmissibility_m3) ||
            *gated.transmissibility_m3 <= 0.0) {
            return PETSC_ERR_ARG_INCOMP;
        }

        const std::size_t face_local =
            static_cast<std::size_t>(
                gated.face.value());
        if (face_local >= face_count) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }
        const PetscInt face_point =
            face_start +
            static_cast<PetscInt>(
                face_local);

        PetscInt support_size = 0;
        const PetscInt* support = nullptr;
        error = DMPlexGetSupportSize(
            target_dm,
            face_point,
            &support_size);
        if (error == PETSC_SUCCESS) {
            error = DMPlexGetSupport(
                target_dm,
                face_point,
                &support);
        }
        if (error != PETSC_SUCCESS) {
            return error;
        }
        if (support_size != 2 ||
            support == nullptr) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }

        const auto* first_identity =
            identity_for_point(
                support[0]);
        const auto* second_identity =
            identity_for_point(
                support[1]);
        if (first_identity == nullptr ||
            second_identity == nullptr ||
            first_identity->kind !=
                mpmc::mesh::EntityKind::cell ||
            second_identity->kind !=
                mpmc::mesh::EntityKind::cell ||
            first_identity->point < cell_start ||
            first_identity->point >= cell_end ||
            second_identity->point < cell_start ||
            second_identity->point >= cell_end ||
            first_identity->local ==
                second_identity->local) {
            return PETSC_ERR_ARG_INCOMP;
        }

        const auto canonical_owner_global =
            stable_face_geometry
                .face_owner_global_ids[
                    face_local];

        const DMPlexPointIdentity* owner = nullptr;
        const DMPlexPointIdentity* neighbour = nullptr;
        if (first_identity->global ==
            canonical_owner_global) {
            owner = first_identity;
            neighbour = second_identity;
        } else if (
            second_identity->global ==
            canonical_owner_global) {
            owner = second_identity;
            neighbour = first_identity;
        } else {
            return PETSC_ERR_ARG_INCOMP;
        }

        rows.push_back(
            AssemblyReadyInternalConnectionRow3D{
                gated.face,
                gated.global,
                owner->local,
                owner->global,
                neighbour->local,
                neighbour->global,
                *gated.transmissibility_m3});
    }

    if (rows.size() !=
        view.materialized_face_count()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    try {
        output->emplace(
            face_count,
            cell_count,
            std::move(rows));
    } catch (...) {
        output->reset();
        return PETSC_ERR_ARG_INCOMP;
    }
    return PETSC_SUCCESS;
}


class ParallelOwnedConnectionSchedule3D {
public:
    ParallelOwnedConnectionSchedule3D(
        mpmc::mesh::PartitionRank local_rank,
        std::uint32_t rank_count,
        std::vector<AssemblyReadyInternalConnectionRow3D>
            authoritative_rows,
        std::vector<AssemblyReadyInternalConnectionRow3D>
            ghost_rows)
        : local_rank_(local_rank),
          rank_count_(rank_count),
          authoritative_rows_(
              std::move(authoritative_rows)),
          ghost_rows_(
              std::move(ghost_rows)) {
        validate();
    }

    ParallelOwnedConnectionSchedule3D(
        const ParallelOwnedConnectionSchedule3D&) =
        default;
    ParallelOwnedConnectionSchedule3D(
        ParallelOwnedConnectionSchedule3D&&) noexcept =
        default;
    ParallelOwnedConnectionSchedule3D& operator=(
        const ParallelOwnedConnectionSchedule3D&) =
        delete;
    ParallelOwnedConnectionSchedule3D& operator=(
        ParallelOwnedConnectionSchedule3D&&) =
        delete;
    ~ParallelOwnedConnectionSchedule3D() = default;

    [[nodiscard]] mpmc::mesh::PartitionRank
    local_rank() const noexcept {
        return local_rank_;
    }

    [[nodiscard]] std::uint32_t
    rank_count() const noexcept {
        return rank_count_;
    }

    /// The only rows permitted to participate in a future assembly loop.
    [[nodiscard]] std::span<
        const AssemblyReadyInternalConnectionRow3D>
    assembly_rows() const noexcept {
        return authoritative_rows_;
    }

    /// Duplicate non-owned rows retained only for stable-ID/value consistency
    /// diagnostics. They are intentionally excluded from assembly_rows().
    [[nodiscard]] std::span<
        const AssemblyReadyInternalConnectionRow3D>
    ghost_rows() const noexcept {
        return ghost_rows_;
    }

    [[nodiscard]] std::size_t
    authoritative_row_count() const noexcept {
        return authoritative_rows_.size();
    }

    [[nodiscard]] std::size_t
    ghost_row_count() const noexcept {
        return ghost_rows_.size();
    }

    [[nodiscard]] std::size_t
    local_copy_count() const noexcept {
        return authoritative_rows_.size() +
               ghost_rows_.size();
    }

private:
    static void validate_row(
        const AssemblyReadyInternalConnectionRow3D&
            row) {
        if (row.owner_cell ==
                row.neighbour_cell ||
            row.owner_cell_global ==
                row.neighbour_cell_global ||
            !std::isfinite(
                row.transmissibility_m3) ||
            row.transmissibility_m3 <= 0.0) {
            throw std::invalid_argument(
                "mpmc::discretization_petsc::ParallelOwnedConnectionSchedule3D: invalid connection row");
        }
    }

    void validate() {
        if (rank_count_ == 0U ||
            local_rank_.value() >= rank_count_) {
            throw std::invalid_argument(
                "mpmc::discretization_petsc::ParallelOwnedConnectionSchedule3D: invalid rank metadata");
        }

        for (const auto& row :
             authoritative_rows_) {
            validate_row(row);
        }
        for (const auto& row :
             ghost_rows_) {
            validate_row(row);
        }

        for (std::size_t i = 0U;
             i < authoritative_rows_.size();
             ++i) {
            for (std::size_t j = 0U;
                 j < i;
                 ++j) {
                if (authoritative_rows_[i].face_global ==
                    authoritative_rows_[j].face_global) {
                    throw std::invalid_argument(
                        "mpmc::discretization_petsc::ParallelOwnedConnectionSchedule3D: duplicate authoritative stable face ID");
                }
            }
            for (const auto& ghost :
                 ghost_rows_) {
                if (authoritative_rows_[i].face_global ==
                    ghost.face_global) {
                    throw std::invalid_argument(
                        "mpmc::discretization_petsc::ParallelOwnedConnectionSchedule3D: one local stable face copy cannot be both authoritative and ghost");
                }
            }
        }
        for (std::size_t i = 0U;
             i < ghost_rows_.size();
             ++i) {
            for (std::size_t j = 0U;
                 j < i;
                 ++j) {
                if (ghost_rows_[i].face_global ==
                    ghost_rows_[j].face_global) {
                    throw std::invalid_argument(
                        "mpmc::discretization_petsc::ParallelOwnedConnectionSchedule3D: duplicate ghost stable face ID");
                }
            }
        }
    }

    mpmc::mesh::PartitionRank local_rank_;
    std::uint32_t rank_count_;
    std::vector<AssemblyReadyInternalConnectionRow3D>
        authoritative_rows_;
    std::vector<AssemblyReadyInternalConnectionRow3D>
        ghost_rows_;
};

/// Split a target-local active connection table by PETSc/core face ownership.
///
/// The authoritative criterion is exclusively PartitionSnapshot ownership of
/// the FACE represented by the row. Canonical geometric owner_cell is not used
/// to choose the MPI assembly rank.
///
/// Future assembly code must iterate schedule.assembly_rows(); ghost_rows()
/// exist only to validate duplicate copies after overlap construction.
inline PetscErrorCode
make_parallel_owned_connection_schedule_3d(
    const AssemblyReadyInternalConnectionTable3D& table,
    const mpmc::mesh::PartitionSnapshot& partition,
    std::optional<ParallelOwnedConnectionSchedule3D>*
        output) {
    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();

    if (partition.entity_count(
            mpmc::mesh::EntityKind::face) !=
            table.target_face_count() ||
        partition.entity_count(
            mpmc::mesh::EntityKind::cell) !=
            table.target_cell_count()) {
        return PETSC_ERR_ARG_SIZ;
    }

    std::vector<AssemblyReadyInternalConnectionRow3D>
        authoritative_rows;
    std::vector<AssemblyReadyInternalConnectionRow3D>
        ghost_rows;
    authoritative_rows.reserve(
        table.row_count());
    ghost_rows.reserve(
        table.row_count());

    for (const auto& row :
         table.rows()) {
        try {
            if (partition.global_id(
                    mpmc::mesh::EntityKind::face,
                    row.face) !=
                    row.face_global ||
                partition.global_id(
                    mpmc::mesh::EntityKind::cell,
                    row.owner_cell) !=
                    row.owner_cell_global ||
                partition.global_id(
                    mpmc::mesh::EntityKind::cell,
                    row.neighbour_cell) !=
                    row.neighbour_cell_global) {
                return PETSC_ERR_ARG_INCOMP;
            }

            if (partition.is_owned(
                    mpmc::mesh::EntityKind::face,
                    row.face)) {
                authoritative_rows.push_back(
                    row);
            } else if (
                partition.is_ghost(
                    mpmc::mesh::EntityKind::face,
                    row.face)) {
                ghost_rows.push_back(
                    row);
            } else {
                return PETSC_ERR_ARG_INCOMP;
            }
        } catch (...) {
            return PETSC_ERR_ARG_INCOMP;
        }
    }

    if (authoritative_rows.size() +
            ghost_rows.size() !=
        table.row_count()) {
        return PETSC_ERR_PLIB;
    }

    try {
        output->emplace(
            partition.local_rank(),
            partition.rank_count(),
            std::move(authoritative_rows),
            std::move(ghost_rows));
    } catch (...) {
        output->reset();
        return PETSC_ERR_ARG_INCOMP;
    }

    return PETSC_SUCCESS;
}


struct CellPairCoupling3D {
    mpmc::mesh::LocalIndex first_cell;
    mpmc::mesh::GlobalEntityId first_cell_global;
    mpmc::mesh::LocalIndex second_cell;
    mpmc::mesh::GlobalEntityId second_cell_global;
};

struct OwnedCellStructuralCounts3D {
    mpmc::mesh::LocalIndex cell;
    mpmc::mesh::GlobalEntityId cell_global;

    /// Structural columns in the local diagonal block:
    /// one self column plus unique neighbouring cells owned by this rank.
    std::size_t diagonal_block_nnz;

    /// Structural columns in the off-diagonal block:
    /// unique neighbouring cells owned by remote ranks and present as ghosts.
    std::size_t off_diagonal_block_nnz;
};

/// Immutable local cell-pair sparsity/stencil snapshot.
///
/// Only schedule.assembly_rows() contribute structural pairs. Each pair is
/// owner-targeted to its endpoint cell owner rank(s), then locally deduplicated;
/// schedule.ghost_rows() are deliberately not consumed. A rank therefore stores
/// only couplings needed by its locally owned cell rows, never a replicated
/// global coupling graph.
///
/// No transmissibility value, pressure variable, matrix coefficient, flux, or
/// residual is stored in this snapshot.
class CellPairSparsityStencilSnapshot3D {
public:
    CellPairSparsityStencilSnapshot3D(
        mpmc::mesh::PartitionRank local_rank,
        std::uint32_t rank_count,
        std::size_t local_cell_count,
        std::vector<CellPairCoupling3D> couplings,
        std::vector<OwnedCellStructuralCounts3D>
            owned_cell_counts)
        : local_rank_(local_rank),
          rank_count_(rank_count),
          local_cell_count_(local_cell_count),
          couplings_(std::move(couplings)),
          owned_cell_counts_(
              std::move(owned_cell_counts)),
          owned_cell_to_counts_(local_cell_count) {
        validate_and_index();
    }

    CellPairSparsityStencilSnapshot3D(
        const CellPairSparsityStencilSnapshot3D&) =
        default;
    CellPairSparsityStencilSnapshot3D(
        CellPairSparsityStencilSnapshot3D&&) noexcept =
        default;
    CellPairSparsityStencilSnapshot3D& operator=(
        const CellPairSparsityStencilSnapshot3D&) =
        delete;
    CellPairSparsityStencilSnapshot3D& operator=(
        CellPairSparsityStencilSnapshot3D&&) =
        delete;
    ~CellPairSparsityStencilSnapshot3D() = default;

    [[nodiscard]] mpmc::mesh::PartitionRank
    local_rank() const noexcept {
        return local_rank_;
    }

    [[nodiscard]] std::uint32_t
    rank_count() const noexcept {
        return rank_count_;
    }

    [[nodiscard]] std::size_t
    local_cell_count() const noexcept {
        return local_cell_count_;
    }

    [[nodiscard]] std::size_t
    coupling_count() const noexcept {
        return couplings_.size();
    }

    [[nodiscard]] std::size_t
    owned_cell_count() const noexcept {
        return owned_cell_counts_.size();
    }

    [[nodiscard]] std::span<const CellPairCoupling3D>
    couplings() const noexcept {
        return couplings_;
    }

    [[nodiscard]] std::span<
        const OwnedCellStructuralCounts3D>
    owned_cell_counts() const noexcept {
        return owned_cell_counts_;
    }

    [[nodiscard]] bool contains_owned_cell(
        mpmc::mesh::LocalIndex cell) const {
        const std::size_t local =
            static_cast<std::size_t>(
                cell.value());
        if (local >= local_cell_count_) {
            throw std::out_of_range(
                "mpmc::discretization_petsc::CellPairSparsityStencilSnapshot3D: cell index out of range");
        }
        return owned_cell_to_counts_[local].has_value();
    }

    [[nodiscard]]
    const OwnedCellStructuralCounts3D&
    structural_counts(
        mpmc::mesh::LocalIndex cell) const {
        const std::size_t local =
            static_cast<std::size_t>(
                cell.value());
        if (local >= local_cell_count_) {
            throw std::out_of_range(
                "mpmc::discretization_petsc::CellPairSparsityStencilSnapshot3D: cell index out of range");
        }
        const auto mapped =
            owned_cell_to_counts_[local];
        if (!mapped.has_value()) {
            throw std::invalid_argument(
                "mpmc::discretization_petsc::CellPairSparsityStencilSnapshot3D: structural counts are defined only for locally owned cells");
        }
        return owned_cell_counts_.at(*mapped);
    }

private:
    void validate_and_index() {
        if (rank_count_ == 0U ||
            local_rank_.value() >= rank_count_) {
            throw std::invalid_argument(
                "mpmc::discretization_petsc::CellPairSparsityStencilSnapshot3D: invalid rank metadata");
        }

        for (std::size_t index = 0U;
             index < couplings_.size();
             ++index) {
            const auto& pair =
                couplings_[index];
            const std::size_t first =
                static_cast<std::size_t>(
                    pair.first_cell.value());
            const std::size_t second =
                static_cast<std::size_t>(
                    pair.second_cell.value());
            if (first >= local_cell_count_ ||
                second >= local_cell_count_ ||
                first == second ||
                pair.first_cell_global ==
                    pair.second_cell_global ||
                !(pair.first_cell_global <
                  pair.second_cell_global)) {
                throw std::invalid_argument(
                    "mpmc::discretization_petsc::CellPairSparsityStencilSnapshot3D: coupling must use two distinct local cells ordered by stable GlobalEntityId");
            }
            for (std::size_t prior = 0U;
                 prior < index;
                 ++prior) {
                if (couplings_[prior]
                            .first_cell_global ==
                        pair.first_cell_global &&
                    couplings_[prior]
                            .second_cell_global ==
                        pair.second_cell_global) {
                    throw std::invalid_argument(
                        "mpmc::discretization_petsc::CellPairSparsityStencilSnapshot3D: duplicate stable cell pair");
                }
            }
        }

        for (std::size_t index = 0U;
             index < owned_cell_counts_.size();
             ++index) {
            const auto& counts =
                owned_cell_counts_[index];
            const std::size_t local =
                static_cast<std::size_t>(
                    counts.cell.value());
            if (local >= local_cell_count_ ||
                counts.diagonal_block_nnz == 0U ||
                owned_cell_to_counts_[local]
                    .has_value()) {
                throw std::invalid_argument(
                    "mpmc::discretization_petsc::CellPairSparsityStencilSnapshot3D: invalid or duplicate owned-cell structural counts");
            }
            owned_cell_to_counts_[local] =
                index;
        }
    }

    mpmc::mesh::PartitionRank local_rank_;
    std::uint32_t rank_count_;
    std::size_t local_cell_count_;
    std::vector<CellPairCoupling3D> couplings_;
    std::vector<OwnedCellStructuralCounts3D>
        owned_cell_counts_;
    std::vector<std::optional<std::size_t>>
        owned_cell_to_counts_;
};

namespace cell_pair_sparsity_detail {

struct StableCellPair {
    std::uint64_t first;
    std::uint64_t second;
};

[[nodiscard]] inline StableCellPair
canonical_pair(
    mpmc::mesh::GlobalEntityId left,
    mpmc::mesh::GlobalEntityId right) {
    const std::uint64_t a = left.value();
    const std::uint64_t b = right.value();
    if (a == b) {
        throw std::invalid_argument(
            "mpmc::discretization_petsc::cell_pair_sparsity: self coupling is invalid");
    }
    return a < b
               ? StableCellPair{a, b}
               : StableCellPair{b, a};
}

[[nodiscard]] inline bool pair_less(
    const StableCellPair& left,
    const StableCellPair& right) noexcept {
    return left.first < right.first ||
           (left.first == right.first &&
            left.second < right.second);
}

[[nodiscard]] inline bool pair_equal(
    const StableCellPair& left,
    const StableCellPair& right) noexcept {
    return left.first == right.first &&
           left.second == right.second;
}

struct TargetedStableCellPair {
    std::uint32_t target_rank;
    StableCellPair pair;
};

[[nodiscard]] inline bool targeted_pair_less(
    const TargetedStableCellPair& left,
    const TargetedStableCellPair& right) noexcept {
    return left.target_rank < right.target_rank ||
           (left.target_rank == right.target_rank &&
            pair_less(left.pair, right.pair));
}

[[nodiscard]] inline bool targeted_pair_equal(
    const TargetedStableCellPair& left,
    const TargetedStableCellPair& right) noexcept {
    return left.target_rank == right.target_rank &&
           pair_equal(left.pair, right.pair);
}

} // namespace cell_pair_sparsity_detail

/// Collectively build the local owned-row cell-pair sparsity snapshot.
///
/// Only schedule.assembly_rows() are consumed. Every authoritative stable cell
/// pair is routed only to the owner rank of each endpoint cell; when both cells
/// share one owner rank the pair is sent once. schedule.ghost_rows() remain
/// diagnostics-only and never participate in the structural graph.
///
/// The exchange is owner-targeted: MPI_Alltoall communicates only per-rank word
/// counts, and MPI_Alltoallv transfers pair payloads solely to endpoint owner
/// ranks. No rank reconstructs or temporarily stores the global cell-pair graph.
///
/// PETSc/MPI diagonal-block counts follow standard distributed sparse-matrix
/// preallocation semantics:
///   diagonal_block_nnz = 1 self + unique neighbours owned by this rank
///   off_diagonal_block_nnz = unique neighbours owned by remote ranks
inline PetscErrorCode
make_cell_pair_sparsity_stencil_snapshot_3d(
    MPI_Comm comm,
    const ParallelOwnedConnectionSchedule3D& schedule,
    const mpmc::mesh::PartitionSnapshot& partition,
    std::optional<CellPairSparsityStencilSnapshot3D>*
        output) {
    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();

    PetscErrorCode error =
        detail::validate_communicator(
            comm,
            schedule.local_rank(),
            schedule.rank_count());
    if (error != PETSC_SUCCESS) return error;
    if (partition.local_rank() !=
            schedule.local_rank() ||
        partition.rank_count() !=
            schedule.rank_count()) {
        return PETSC_ERR_ARG_INCOMP;
    }

    int mpi_size = 0;
    if (MPI_Comm_size(
            comm,
            &mpi_size) != MPI_SUCCESS ||
        mpi_size <= 0 ||
        static_cast<std::uint32_t>(mpi_size) !=
            partition.rank_count()) {
        return PETSC_ERR_MPI;
    }
    const std::size_t rank_count =
        static_cast<std::size_t>(mpi_size);

    const std::size_t local_pair_count =
        schedule.assembly_rows().size();
    if (local_pair_count >
        std::numeric_limits<std::size_t>::max() /
            2U) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    std::vector<
        cell_pair_sparsity_detail::
            TargetedStableCellPair>
        outbound_pairs;
    outbound_pairs.reserve(
        local_pair_count * 2U);

    try {
        for (const auto& row :
             schedule.assembly_rows()) {
            if (partition.global_id(
                    mpmc::mesh::EntityKind::cell,
                    row.owner_cell) !=
                    row.owner_cell_global ||
                partition.global_id(
                    mpmc::mesh::EntityKind::cell,
                    row.neighbour_cell) !=
                    row.neighbour_cell_global) {
                return PETSC_ERR_ARG_INCOMP;
            }

            const auto pair =
                cell_pair_sparsity_detail::
                    canonical_pair(
                        row.owner_cell_global,
                        row.neighbour_cell_global);
            const auto owner_rank =
                partition.owner_rank(
                    mpmc::mesh::EntityKind::cell,
                    row.owner_cell);
            const auto neighbour_rank =
                partition.owner_rank(
                    mpmc::mesh::EntityKind::cell,
                    row.neighbour_cell);

            outbound_pairs.push_back(
                cell_pair_sparsity_detail::
                    TargetedStableCellPair{
                        owner_rank.value(),
                        pair});
            if (neighbour_rank !=
                owner_rank) {
                outbound_pairs.push_back(
                    cell_pair_sparsity_detail::
                        TargetedStableCellPair{
                            neighbour_rank.value(),
                            pair});
            }
        }
    } catch (...) {
        return PETSC_ERR_ARG_INCOMP;
    }

    std::sort(
        outbound_pairs.begin(),
        outbound_pairs.end(),
        cell_pair_sparsity_detail::
            targeted_pair_less);
    outbound_pairs.erase(
        std::unique(
            outbound_pairs.begin(),
            outbound_pairs.end(),
            cell_pair_sparsity_detail::
                targeted_pair_equal),
        outbound_pairs.end());

    if (outbound_pairs.size() >
        static_cast<std::size_t>(
            std::numeric_limits<int>::max() /
            2)) {
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    std::vector<int> send_word_counts(
        rank_count,
        0);
    for (const auto& targeted :
         outbound_pairs) {
        const std::size_t target =
            static_cast<std::size_t>(
                targeted.target_rank);
        if (target >= rank_count ||
            send_word_counts[target] >
                std::numeric_limits<int>::max() -
                    2) {
            return PETSC_ERR_ARG_INCOMP;
        }
        send_word_counts[target] += 2;
    }

    std::vector<int> send_displacements(
        rank_count,
        0);
    int total_send_words = 0;
    for (std::size_t rank = 0U;
         rank < rank_count;
         ++rank) {
        const int count =
            send_word_counts[rank];
        if (count < 0 ||
            count % 2 != 0 ||
            count >
                std::numeric_limits<int>::max() -
                    total_send_words) {
            return PETSC_ERR_ARG_INCOMP;
        }
        send_displacements[rank] =
            total_send_words;
        total_send_words += count;
    }

    std::vector<std::uint64_t> send_words;
    send_words.reserve(
        static_cast<std::size_t>(
            total_send_words));
    for (const auto& targeted :
         outbound_pairs) {
        send_words.push_back(
            targeted.pair.first);
        send_words.push_back(
            targeted.pair.second);
    }
    if (send_words.size() !=
        static_cast<std::size_t>(
            total_send_words)) {
        return PETSC_ERR_PLIB;
    }

    std::vector<int> receive_word_counts(
        rank_count,
        0);
    if (MPI_Alltoall(
            send_word_counts.data(),
            1,
            MPI_INT,
            receive_word_counts.data(),
            1,
            MPI_INT,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    std::vector<int> receive_displacements(
        rank_count,
        0);
    int total_receive_words = 0;
    for (std::size_t rank = 0U;
         rank < rank_count;
         ++rank) {
        const int count =
            receive_word_counts[rank];
        if (count < 0 ||
            count % 2 != 0 ||
            count >
                std::numeric_limits<int>::max() -
                    total_receive_words) {
            return PETSC_ERR_ARG_INCOMP;
        }
        receive_displacements[rank] =
            total_receive_words;
        total_receive_words += count;
    }

    std::vector<std::uint64_t> received_words(
        static_cast<std::size_t>(
            total_receive_words));
    if (MPI_Alltoallv(
            send_words.empty()
                ? nullptr
                : send_words.data(),
            send_word_counts.data(),
            send_displacements.data(),
            MPI_UINT64_T,
            received_words.empty()
                ? nullptr
                : received_words.data(),
            receive_word_counts.data(),
            receive_displacements.data(),
            MPI_UINT64_T,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }

    std::vector<
        cell_pair_sparsity_detail::StableCellPair>
            unique_pairs;
    unique_pairs.reserve(
        received_words.size() / 2U);
    for (std::size_t offset = 0U;
         offset < received_words.size();
         offset += 2U) {
        const auto first =
            received_words[offset];
        const auto second =
            received_words[offset + 1U];
        if (first >= second) {
            return PETSC_ERR_ARG_INCOMP;
        }
        unique_pairs.push_back(
            cell_pair_sparsity_detail::
                StableCellPair{
                    first, second});
    }
    std::sort(
        unique_pairs.begin(),
        unique_pairs.end(),
        cell_pair_sparsity_detail::
            pair_less);
    unique_pairs.erase(
        std::unique(
            unique_pairs.begin(),
            unique_pairs.end(),
            cell_pair_sparsity_detail::
                pair_equal),
        unique_pairs.end());

    const std::size_t local_cell_count =
        partition.entity_count(
            mpmc::mesh::EntityKind::cell);
    std::vector<CellPairCoupling3D>
        local_couplings;
    local_couplings.reserve(
        unique_pairs.size());

    for (const auto pair :
         unique_pairs) {
        const auto first_global =
            mpmc::mesh::GlobalEntityId{
                pair.first};
        const auto second_global =
            mpmc::mesh::GlobalEntityId{
                pair.second};

        if (!partition.contains_global(
                mpmc::mesh::EntityKind::cell,
                first_global) ||
            !partition.contains_global(
                mpmc::mesh::EntityKind::cell,
                second_global)) {
            return PETSC_ERR_ARG_WRONGSTATE;
        }

        mpmc::mesh::LocalIndex first_local_index{
            0U};
        mpmc::mesh::LocalIndex second_local_index{
            0U};
        bool first_owned = false;
        bool second_owned = false;
        try {
            first_local_index =
                partition.local_index(
                    mpmc::mesh::EntityKind::cell,
                    first_global);
            second_local_index =
                partition.local_index(
                    mpmc::mesh::EntityKind::cell,
                    second_global);
            first_owned =
                partition.is_owned(
                    mpmc::mesh::EntityKind::cell,
                    first_local_index);
            second_owned =
                partition.is_owned(
                    mpmc::mesh::EntityKind::cell,
                    second_local_index);
        } catch (...) {
            return PETSC_ERR_ARG_INCOMP;
        }

        if (!first_owned &&
            !second_owned) {
            return PETSC_ERR_ARG_INCOMP;
        }

        local_couplings.push_back(
            CellPairCoupling3D{
                first_local_index,
                first_global,
                second_local_index,
                second_global});
    }

    std::vector<OwnedCellStructuralCounts3D>
        owned_counts;
    owned_counts.reserve(
        partition.owned_count(
            mpmc::mesh::EntityKind::cell));

    for (std::size_t local = 0U;
         local < local_cell_count;
         ++local) {
        const auto cell =
            mpmc::mesh::LocalIndex{
                static_cast<
                    mpmc::mesh::LocalIndex::value_type>(
                        local)};
        if (!partition.is_owned(
                mpmc::mesh::EntityKind::cell,
                cell)) {
            continue;
        }

        std::size_t diagonal_block_nnz = 1U;
        std::size_t off_diagonal_block_nnz = 0U;

        for (const auto& pair :
             local_couplings) {
            mpmc::mesh::LocalIndex neighbour{
                0U};
            bool incident = false;
            if (pair.first_cell == cell) {
                neighbour = pair.second_cell;
                incident = true;
            } else if (
                pair.second_cell == cell) {
                neighbour = pair.first_cell;
                incident = true;
            }
            if (!incident) continue;

            if (partition.is_owned(
                    mpmc::mesh::EntityKind::cell,
                    neighbour)) {
                ++diagonal_block_nnz;
            } else if (
                partition.is_ghost(
                    mpmc::mesh::EntityKind::cell,
                    neighbour)) {
                ++off_diagonal_block_nnz;
            } else {
                return PETSC_ERR_ARG_INCOMP;
            }
        }

        owned_counts.push_back(
            OwnedCellStructuralCounts3D{
                cell,
                partition.global_id(
                    mpmc::mesh::EntityKind::cell,
                    cell),
                diagonal_block_nnz,
                off_diagonal_block_nnz});
    }

    if (owned_counts.size() !=
        partition.owned_count(
            mpmc::mesh::EntityKind::cell)) {
        return PETSC_ERR_PLIB;
    }

    try {
        output->emplace(
            partition.local_rank(),
            partition.rank_count(),
            local_cell_count,
            std::move(local_couplings),
            std::move(owned_counts));
    } catch (...) {
        output->reset();
        return PETSC_ERR_ARG_INCOMP;
    }

    return PETSC_SUCCESS;
}

class PetscMpiAijSymbolicPreallocation3D {
public:
    PetscMpiAijSymbolicPreallocation3D(
        mpmc::mesh::PartitionRank local_rank,
        std::uint32_t rank_count,
        PetscInt global_row_start,
        PetscInt global_row_end,
        PetscInt global_row_count,
        std::vector<mpmc::mesh::LocalIndex>
            owned_cells_in_petsc_row_order,
        std::vector<mpmc::mesh::GlobalEntityId>
            owned_cell_global_ids,
        std::vector<PetscInt> owned_global_rows,
        std::vector<PetscInt> diagonal_nnz,
        std::vector<PetscInt> off_diagonal_nnz,
        std::vector<PetscInt> local_cell_global_rows)
        : local_rank_(local_rank),
          rank_count_(rank_count),
          global_row_start_(global_row_start),
          global_row_end_(global_row_end),
          global_row_count_(global_row_count),
          owned_cells_in_petsc_row_order_(
              std::move(
                  owned_cells_in_petsc_row_order)),
          owned_cell_global_ids_(
              std::move(owned_cell_global_ids)),
          owned_global_rows_(
              std::move(owned_global_rows)),
          diagonal_nnz_(
              std::move(diagonal_nnz)),
          off_diagonal_nnz_(
              std::move(off_diagonal_nnz)),
          local_cell_global_rows_(
              std::move(local_cell_global_rows)) {
        validate();
    }

    PetscMpiAijSymbolicPreallocation3D(
        const PetscMpiAijSymbolicPreallocation3D&) =
        default;
    PetscMpiAijSymbolicPreallocation3D(
        PetscMpiAijSymbolicPreallocation3D&&) noexcept =
        default;
    PetscMpiAijSymbolicPreallocation3D& operator=(
        const PetscMpiAijSymbolicPreallocation3D&) =
        delete;
    PetscMpiAijSymbolicPreallocation3D& operator=(
        PetscMpiAijSymbolicPreallocation3D&&) =
        delete;
    ~PetscMpiAijSymbolicPreallocation3D() = default;

    [[nodiscard]] mpmc::mesh::PartitionRank
    local_rank() const noexcept {
        return local_rank_;
    }

    [[nodiscard]] std::uint32_t
    rank_count() const noexcept {
        return rank_count_;
    }

    [[nodiscard]] PetscInt
    global_row_start() const noexcept {
        return global_row_start_;
    }

    [[nodiscard]] PetscInt
    global_row_end() const noexcept {
        return global_row_end_;
    }

    [[nodiscard]] PetscInt
    global_row_count() const noexcept {
        return global_row_count_;
    }

    [[nodiscard]] PetscInt
    local_owned_row_count() const noexcept {
        return global_row_end_ -
               global_row_start_;
    }

    [[nodiscard]] std::span<
        const mpmc::mesh::LocalIndex>
    owned_cells_in_petsc_row_order() const noexcept {
        return owned_cells_in_petsc_row_order_;
    }

    [[nodiscard]] std::span<
        const mpmc::mesh::GlobalEntityId>
    owned_cell_global_ids() const noexcept {
        return owned_cell_global_ids_;
    }

    [[nodiscard]] std::span<const PetscInt>
    owned_global_rows() const noexcept {
        return owned_global_rows_;
    }

    /// Direct inputs for a future MatMPIAIJSetPreallocation() call.
    [[nodiscard]] std::span<const PetscInt>
    diagonal_nnz() const noexcept {
        return diagonal_nnz_;
    }

    [[nodiscard]] std::span<const PetscInt>
    off_diagonal_nnz() const noexcept {
        return off_diagonal_nnz_;
    }

    /// PETSc global row for any local cell, owned or ghost, resolved through
    /// the collectively verified stable-cell numbering.
    [[nodiscard]] PetscInt global_row(
        mpmc::mesh::LocalIndex local_cell) const {
        const std::size_t local =
            static_cast<std::size_t>(
                local_cell.value());
        if (local >=
            local_cell_global_rows_.size()) {
            throw std::out_of_range(
                "mpmc::discretization_petsc::PetscMpiAijSymbolicPreallocation3D: local cell index out of range");
        }
        return local_cell_global_rows_[local];
    }

private:
    void validate() {
        if (rank_count_ == 0U ||
            local_rank_.value() >= rank_count_ ||
            global_row_start_ < 0 ||
            global_row_end_ < global_row_start_ ||
            global_row_count_ < global_row_end_) {
            throw std::invalid_argument(
                "mpmc::discretization_petsc::PetscMpiAijSymbolicPreallocation3D: invalid PETSc row ownership metadata");
        }

        const std::size_t owned_count =
            owned_cells_in_petsc_row_order_.size();
        if (owned_cell_global_ids_.size() !=
                owned_count ||
            owned_global_rows_.size() !=
                owned_count ||
            diagonal_nnz_.size() !=
                owned_count ||
            off_diagonal_nnz_.size() !=
                owned_count ||
            static_cast<std::size_t>(
                global_row_end_ -
                global_row_start_) !=
                owned_count) {
            throw std::invalid_argument(
                "mpmc::discretization_petsc::PetscMpiAijSymbolicPreallocation3D: preallocation arrays do not match owned PETSc row count");
        }

        for (std::size_t i = 0U;
             i < owned_count;
             ++i) {
            if (owned_global_rows_[i] !=
                    global_row_start_ +
                        static_cast<PetscInt>(i) ||
                diagonal_nnz_[i] <= 0 ||
                off_diagonal_nnz_[i] < 0) {
                throw std::invalid_argument(
                    "mpmc::discretization_petsc::PetscMpiAijSymbolicPreallocation3D: invalid row order or structural nnz count");
            }
            if (i > 0U &&
                !(owned_cell_global_ids_[i - 1U] <
                  owned_cell_global_ids_[i])) {
                throw std::invalid_argument(
                    "mpmc::discretization_petsc::PetscMpiAijSymbolicPreallocation3D: owned stable cell IDs must be strictly ordered within PETSc ownership range");
            }
        }

        for (const PetscInt row :
             local_cell_global_rows_) {
            if (row < 0 ||
                row >= global_row_count_) {
                throw std::invalid_argument(
                    "mpmc::discretization_petsc::PetscMpiAijSymbolicPreallocation3D: local stable cell did not resolve to a valid PETSc global row");
            }
        }
    }

    mpmc::mesh::PartitionRank local_rank_;
    std::uint32_t rank_count_;
    PetscInt global_row_start_;
    PetscInt global_row_end_;
    PetscInt global_row_count_;
    std::vector<mpmc::mesh::LocalIndex>
        owned_cells_in_petsc_row_order_;
    std::vector<mpmc::mesh::GlobalEntityId>
        owned_cell_global_ids_;
    std::vector<PetscInt>
        owned_global_rows_;
    std::vector<PetscInt>
        diagonal_nnz_;
    std::vector<PetscInt>
        off_diagonal_nnz_;
    std::vector<PetscInt>
        local_cell_global_rows_;
};

/// Convert owned-cell structural counts into PETSc MPIAIJ symbolic
/// preallocation arrays and verify a one-to-one stable-cell <-> PETSc-row map.
///
/// This routine creates only a transient PetscLayout to obtain PETSc's actual
/// contiguous rank ownership ranges. It does NOT create a Mat, insert matrix
/// values, define a pressure variable, or consume any transmissibility.
inline PetscErrorCode
make_petsc_mpiaij_symbolic_preallocation_3d(
    DM target_dm,
    const CellPairSparsityStencilSnapshot3D& sparsity,
    const mpmc::mesh::PartitionSnapshot& partition,
    std::span<const DMPlexPointIdentity> target_identities,
    std::optional<PetscMpiAijSymbolicPreallocation3D>*
        output) {
    if (target_dm == nullptr ||
        output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();

    MPI_Comm comm =
        PetscObjectComm(
            reinterpret_cast<PetscObject>(
                target_dm));
    PetscErrorCode error =
        detail::validate_communicator(
            comm,
            partition.local_rank(),
            partition.rank_count());
    if (error != PETSC_SUCCESS) return error;

    if (sparsity.local_rank() !=
            partition.local_rank() ||
        sparsity.rank_count() !=
            partition.rank_count() ||
        sparsity.local_cell_count() !=
            partition.entity_count(
                mpmc::mesh::EntityKind::cell) ||
        sparsity.owned_cell_count() !=
            partition.owned_count(
                mpmc::mesh::EntityKind::cell)) {
        return PETSC_ERR_ARG_INCOMP;
    }

    PetscInt cell_start = -1;
    PetscInt cell_end = -1;
    error = DMPlexGetHeightStratum(
        target_dm,
        0,
        &cell_start,
        &cell_end);
    if (error != PETSC_SUCCESS ||
        cell_start < 0 ||
        cell_end < cell_start ||
        static_cast<std::size_t>(
            cell_end - cell_start) !=
            partition.entity_count(
                mpmc::mesh::EntityKind::cell)) {
        return error != PETSC_SUCCESS
                   ? error
                   : PETSC_ERR_ARG_SIZ;
    }

    std::vector<PetscInt>
        cell_points(
            partition.entity_count(
                mpmc::mesh::EntityKind::cell),
            PetscInt{-1});
    std::vector<std::uint8_t>
        cell_identity_seen(
            cell_points.size(),
            std::uint8_t{0U});

    for (const auto& identity :
         target_identities) {
        if (identity.kind !=
            mpmc::mesh::EntityKind::cell) {
            continue;
        }
        const std::size_t local =
            static_cast<std::size_t>(
                identity.local.value());
        if (local >= cell_points.size() ||
            cell_identity_seen[local] !=
                std::uint8_t{0U} ||
            identity.point < cell_start ||
            identity.point >= cell_end ||
            identity.point - cell_start !=
                static_cast<PetscInt>(local) ||
            partition.global_id(
                mpmc::mesh::EntityKind::cell,
                identity.local) !=
                identity.global) {
            return PETSC_ERR_ARG_INCOMP;
        }
        cell_points[local] =
            identity.point;
        cell_identity_seen[local] =
            std::uint8_t{1U};
    }
    if (std::find(
            cell_identity_seen.begin(),
            cell_identity_seen.end(),
            std::uint8_t{0U}) !=
        cell_identity_seen.end()) {
        return PETSC_ERR_ARG_SIZ;
    }

    PetscInt local_owned_rows = 0;
    error = detail::checked_petsc_int_size(
        partition.owned_count(
            mpmc::mesh::EntityKind::cell),
        &local_owned_rows);
    if (error != PETSC_SUCCESS) return error;

    PetscInt global_rows = 0;
    if (MPI_Allreduce(
            &local_owned_rows,
            &global_rows,
            1,
            MPIU_INT,
            MPI_SUM,
            comm) != MPI_SUCCESS) {
        return PETSC_ERR_MPI;
    }
    if (global_rows < local_owned_rows) {
        return PETSC_ERR_PLIB;
    }

    PetscLayout layout = nullptr;
    error = PetscLayoutCreate(
        comm,
        &layout);
    if (error != PETSC_SUCCESS) return error;

    error = PetscLayoutSetLocalSize(
        layout,
        local_owned_rows);
    if (error == PETSC_SUCCESS) {
        error = PetscLayoutSetSize(
            layout,
            global_rows);
    }
    if (error == PETSC_SUCCESS) {
        error = PetscLayoutSetUp(
            layout);
    }

    PetscInt row_start = -1;
    PetscInt row_end = -1;
    if (error == PETSC_SUCCESS) {
        error = PetscLayoutGetRange(
            layout,
            &row_start,
            &row_end);
    }
    if (error != PETSC_SUCCESS ||
        row_start < 0 ||
        row_end < row_start ||
        row_end - row_start !=
            local_owned_rows) {
        PetscLayoutDestroy(&layout);
        return error != PETSC_SUCCESS
                   ? error
                   : PETSC_ERR_PLIB;
    }

    std::vector<
        OwnedCellStructuralCounts3D>
            ordered_counts(
                sparsity
                    .owned_cell_counts()
                    .begin(),
                sparsity
                    .owned_cell_counts()
                    .end());
    std::sort(
        ordered_counts.begin(),
        ordered_counts.end(),
        [](const auto& left,
           const auto& right) {
            return left.cell_global <
                   right.cell_global;
        });

    std::vector<mpmc::mesh::LocalIndex>
        owned_cells;
    std::vector<mpmc::mesh::GlobalEntityId>
        owned_global_ids;
    std::vector<PetscInt>
        owned_global_rows;
    std::vector<PetscInt>
        diagonal_nnz;
    std::vector<PetscInt>
        off_diagonal_nnz;

    const std::size_t owned_count =
        ordered_counts.size();
    owned_cells.reserve(owned_count);
    owned_global_ids.reserve(owned_count);
    owned_global_rows.reserve(owned_count);
    diagonal_nnz.reserve(owned_count);
    off_diagonal_nnz.reserve(owned_count);

    std::vector<PetscInt>
        local_cell_global_rows(
            partition.entity_count(
                mpmc::mesh::EntityKind::cell),
            PetscInt{-1});

    PetscInt point_sf_roots = -1;
    PetscInt point_sf_leaves = -1;
    const PetscInt* point_sf_ilocal = nullptr;
    const PetscSFNode* point_sf_remote = nullptr;
    PetscSF point_sf = nullptr;
    error = DMGetPointSF(
        target_dm,
        &point_sf);
    if (error == PETSC_SUCCESS) {
        error = PetscSFGetGraph(
            point_sf,
            &point_sf_roots,
            &point_sf_leaves,
            &point_sf_ilocal,
            &point_sf_remote);
    }
    if (error != PETSC_SUCCESS ||
        point_sf == nullptr ||
        point_sf_roots < 0 ||
        point_sf_leaves < 0) {
        PetscLayoutDestroy(&layout);
        return error != PETSC_SUCCESS
                   ? error
                   : PETSC_ERR_ARG_WRONGSTATE;
    }

    std::vector<PetscInt>
        root_rows(
            static_cast<std::size_t>(
                point_sf_roots),
            PetscInt{-1});
    std::vector<PetscInt>
        leaf_rows(
            static_cast<std::size_t>(
                point_sf_roots),
            PetscInt{-1});

    for (std::size_t i = 0U;
         i < owned_count;
         ++i) {
        const auto& counts =
            ordered_counts[i];
        try {
            if (!partition.is_owned(
                    mpmc::mesh::EntityKind::cell,
                    counts.cell) ||
                partition.global_id(
                    mpmc::mesh::EntityKind::cell,
                    counts.cell) !=
                    counts.cell_global) {
                PetscLayoutDestroy(&layout);
                return PETSC_ERR_ARG_INCOMP;
            }
        } catch (...) {
            PetscLayoutDestroy(&layout);
            return PETSC_ERR_ARG_INCOMP;
        }

        const std::size_t local =
            static_cast<std::size_t>(
                counts.cell.value());
        const PetscInt point =
            cell_points[local];
        if (point < 0 ||
            point >= point_sf_roots) {
            PetscLayoutDestroy(&layout);
            return PETSC_ERR_ARG_INCOMP;
        }

        PetscInt d_nnz = 0;
        PetscInt o_nnz = 0;
        error = detail::checked_petsc_int_size(
            counts.diagonal_block_nnz,
            &d_nnz);
        if (error == PETSC_SUCCESS) {
            error = detail::checked_petsc_int_size(
                counts.off_diagonal_block_nnz,
                &o_nnz);
        }
        if (error != PETSC_SUCCESS ||
            d_nnz <= 0 ||
            d_nnz > local_owned_rows ||
            o_nnz < 0 ||
            o_nnz >
                global_rows -
                    local_owned_rows) {
            PetscLayoutDestroy(&layout);
            return error != PETSC_SUCCESS
                       ? error
                       : PETSC_ERR_ARG_INCOMP;
        }

        const PetscInt global_row =
            row_start +
            static_cast<PetscInt>(i);
        root_rows[
            static_cast<std::size_t>(
                point)] =
            global_row;
        local_cell_global_rows[local] =
            global_row;

        owned_cells.push_back(
            counts.cell);
        owned_global_ids.push_back(
            counts.cell_global);
        owned_global_rows.push_back(
            global_row);
        diagonal_nnz.push_back(
            d_nnz);
        off_diagonal_nnz.push_back(
            o_nnz);
    }

    if (owned_count !=
        static_cast<std::size_t>(
            local_owned_rows)) {
        PetscLayoutDestroy(&layout);
        return PETSC_ERR_PLIB;
    }

    error = PetscSFBcastBegin(
        point_sf,
        MPIU_INT,
        root_rows.data(),
        leaf_rows.data(),
        MPI_REPLACE);
    if (error == PETSC_SUCCESS) {
        error = PetscSFBcastEnd(
            point_sf,
            MPIU_INT,
            root_rows.data(),
            leaf_rows.data(),
            MPI_REPLACE);
    }
    if (error != PETSC_SUCCESS) {
        PetscLayoutDestroy(&layout);
        return error;
    }

    std::vector<std::uint8_t>
        ghost_leaf_seen(
            cell_points.size(),
            std::uint8_t{0U});
    for (PetscInt leaf = 0;
         leaf < point_sf_leaves;
         ++leaf) {
        const PetscInt point =
            point_sf_ilocal != nullptr
                ? point_sf_ilocal[leaf]
                : leaf;
        if (point < cell_start ||
            point >= cell_end) {
            continue;
        }
        const std::size_t local =
            static_cast<std::size_t>(
                point - cell_start);
        if (local >=
                local_cell_global_rows.size() ||
            !partition.is_ghost(
                mpmc::mesh::EntityKind::cell,
                mpmc::mesh::LocalIndex{
                    static_cast<
                        mpmc::mesh::LocalIndex::value_type>(
                            local)}) ||
            point_sf_remote == nullptr ||
            point_sf_remote[leaf].rank !=
                static_cast<PetscMPIInt>(
                    partition.owner_rank(
                        mpmc::mesh::EntityKind::cell,
                        mpmc::mesh::LocalIndex{
                            static_cast<
                                mpmc::mesh::LocalIndex::value_type>(
                                    local)})
                        .value()) ||
            leaf_rows[
                static_cast<std::size_t>(
                    point)] < 0) {
            PetscLayoutDestroy(&layout);
            return PETSC_ERR_ARG_INCOMP;
        }

        local_cell_global_rows[local] =
            leaf_rows[
                static_cast<std::size_t>(
                    point)];
        ghost_leaf_seen[local] =
            std::uint8_t{1U};
    }

    for (std::size_t local = 0U;
         local <
             local_cell_global_rows.size();
         ++local) {
        const auto cell =
            mpmc::mesh::LocalIndex{
                static_cast<
                    mpmc::mesh::LocalIndex::value_type>(
                        local)};
        const PetscInt global_row =
            local_cell_global_rows[local];
        if (global_row < 0 ||
            global_row >= global_rows) {
            PetscLayoutDestroy(&layout);
            return PETSC_ERR_ARG_INCOMP;
        }

        const bool row_is_local =
            global_row >= row_start &&
            global_row < row_end;
        if (partition.is_owned(
                mpmc::mesh::EntityKind::cell,
                cell)) {
            if (!row_is_local) {
                PetscLayoutDestroy(&layout);
                return PETSC_ERR_ARG_INCOMP;
            }
        } else {
            if (!partition.is_ghost(
                    mpmc::mesh::EntityKind::cell,
                    cell) ||
                row_is_local ||
                ghost_leaf_seen[local] ==
                    std::uint8_t{0U}) {
                PetscLayoutDestroy(&layout);
                return PETSC_ERR_ARG_INCOMP;
            }
        }
    }

    // Stable numbering is deterministic within each PETSc ownership range:
    // owned cells are strictly ordered by stable GlobalEntityId.
    for (std::size_t i = 1U;
         i < owned_global_ids.size();
         ++i) {
        if (!(owned_global_ids[i - 1U] <
              owned_global_ids[i])) {
            PetscLayoutDestroy(&layout);
            return PETSC_ERR_ARG_INCOMP;
        }
    }

    const PetscErrorCode destroy_error =
        PetscLayoutDestroy(&layout);
    if (destroy_error != PETSC_SUCCESS) {
        return destroy_error;
    }

    try {
        output->emplace(
            partition.local_rank(),
            partition.rank_count(),
            row_start,
            row_end,
            global_rows,
            std::move(owned_cells),
            std::move(owned_global_ids),
            std::move(owned_global_rows),
            std::move(diagonal_nnz),
            std::move(off_diagonal_nnz),
            std::move(local_cell_global_rows));
    } catch (...) {
        output->reset();
        return PETSC_ERR_ARG_INCOMP;
    }

    return PETSC_SUCCESS;
}

struct OwnedCellStructuralColumnPatternRow3D {
    mpmc::mesh::LocalIndex cell;
    mpmc::mesh::GlobalEntityId cell_global;
    PetscInt global_row;
    std::size_t diagonal_column_offset;
    std::size_t diagonal_column_count;
    std::size_t off_diagonal_column_offset;
    std::size_t off_diagonal_column_count;
};

/// Immutable exact PETSc global-column pattern for locally owned cell rows.
///
/// Rows are stored in the same PETSc/stable-cell order as
/// PetscMpiAijSymbolicPreallocation3D. Finalized column storage is compact:
/// row metadata carries offsets/counts into one flat diagonal-column array and
/// one flat off-diagonal-column array. There is no per-row heap allocation.
///
/// Diagonal columns contain self plus locally owned neighbours. Off-diagonal
/// columns contain remote-owned neighbours visible as local ghosts. Each row's
/// two spans are strictly sorted and deduplicated.
///
/// This snapshot is structural only: it does not own or modify a Mat and does
/// not store any numerical matrix coefficient or physical quantity.
class OwnedCellStructuralColumnPatternSnapshot3D {
public:
    OwnedCellStructuralColumnPatternSnapshot3D(
        mpmc::mesh::PartitionRank local_rank,
        std::uint32_t rank_count,
        std::size_t local_cell_count,
        PetscInt global_row_start,
        PetscInt global_row_end,
        PetscInt global_row_count,
        std::vector<
            OwnedCellStructuralColumnPatternRow3D>
            rows,
        std::vector<PetscInt>
            diagonal_global_columns,
        std::vector<PetscInt>
            off_diagonal_global_columns)
        : local_rank_(local_rank),
          rank_count_(rank_count),
          local_cell_count_(local_cell_count),
          global_row_start_(global_row_start),
          global_row_end_(global_row_end),
          global_row_count_(global_row_count),
          rows_(std::move(rows)),
          diagonal_global_columns_(
              std::move(
                  diagonal_global_columns)),
          off_diagonal_global_columns_(
              std::move(
                  off_diagonal_global_columns)),
          local_cell_to_row_(local_cell_count) {
        validate_and_index();
    }

    OwnedCellStructuralColumnPatternSnapshot3D(
        const OwnedCellStructuralColumnPatternSnapshot3D&) =
        default;
    OwnedCellStructuralColumnPatternSnapshot3D(
        OwnedCellStructuralColumnPatternSnapshot3D&&) noexcept =
        default;
    OwnedCellStructuralColumnPatternSnapshot3D& operator=(
        const OwnedCellStructuralColumnPatternSnapshot3D&) =
        delete;
    OwnedCellStructuralColumnPatternSnapshot3D& operator=(
        OwnedCellStructuralColumnPatternSnapshot3D&&) =
        delete;
    ~OwnedCellStructuralColumnPatternSnapshot3D() =
        default;

    [[nodiscard]] mpmc::mesh::PartitionRank
    local_rank() const noexcept {
        return local_rank_;
    }

    [[nodiscard]] std::uint32_t
    rank_count() const noexcept {
        return rank_count_;
    }

    [[nodiscard]] std::size_t
    local_cell_count() const noexcept {
        return local_cell_count_;
    }

    [[nodiscard]] PetscInt
    global_row_start() const noexcept {
        return global_row_start_;
    }

    [[nodiscard]] PetscInt
    global_row_end() const noexcept {
        return global_row_end_;
    }

    [[nodiscard]] PetscInt
    global_row_count() const noexcept {
        return global_row_count_;
    }

    [[nodiscard]] std::size_t
    row_count() const noexcept {
        return rows_.size();
    }

    [[nodiscard]] std::span<
        const OwnedCellStructuralColumnPatternRow3D>
    rows() const noexcept {
        return rows_;
    }

    [[nodiscard]] std::span<const PetscInt>
    diagonal_column_storage() const noexcept {
        return diagonal_global_columns_;
    }

    [[nodiscard]] std::span<const PetscInt>
    off_diagonal_column_storage() const noexcept {
        return off_diagonal_global_columns_;
    }

    [[nodiscard]] bool contains_owned_cell(
        mpmc::mesh::LocalIndex cell) const {
        const std::size_t local =
            static_cast<std::size_t>(
                cell.value());
        if (local >= local_cell_count_) {
            throw std::out_of_range(
                "mpmc::discretization_petsc::OwnedCellStructuralColumnPatternSnapshot3D: cell index out of range");
        }
        return local_cell_to_row_[local]
            .has_value();
    }

    [[nodiscard]]
    const OwnedCellStructuralColumnPatternRow3D&
    row(mpmc::mesh::LocalIndex cell) const {
        return rows_.at(
            row_index(cell));
    }

    [[nodiscard]] std::span<const PetscInt>
    diagonal_global_columns(
        mpmc::mesh::LocalIndex cell) const {
        const auto& pattern =
            rows_.at(
                row_index(cell));
        return std::span<const PetscInt>{
            diagonal_global_columns_}
            .subspan(
                pattern.diagonal_column_offset,
                pattern.diagonal_column_count);
    }

    [[nodiscard]] std::span<const PetscInt>
    off_diagonal_global_columns(
        mpmc::mesh::LocalIndex cell) const {
        const auto& pattern =
            rows_.at(
                row_index(cell));
        return std::span<const PetscInt>{
            off_diagonal_global_columns_}
            .subspan(
                pattern.off_diagonal_column_offset,
                pattern.off_diagonal_column_count);
    }

private:
    [[nodiscard]] std::size_t row_index(
        mpmc::mesh::LocalIndex cell) const {
        const std::size_t local =
            static_cast<std::size_t>(
                cell.value());
        if (local >= local_cell_count_) {
            throw std::out_of_range(
                "mpmc::discretization_petsc::OwnedCellStructuralColumnPatternSnapshot3D: cell index out of range");
        }
        const auto mapped =
            local_cell_to_row_[local];
        if (!mapped.has_value()) {
            throw std::invalid_argument(
                "mpmc::discretization_petsc::OwnedCellStructuralColumnPatternSnapshot3D: column pattern is defined only for locally owned cells");
        }
        return *mapped;
    }

    void validate_and_index() {
        if (rank_count_ == 0U ||
            local_rank_.value() >= rank_count_ ||
            global_row_start_ < 0 ||
            global_row_end_ < global_row_start_ ||
            global_row_count_ < global_row_end_ ||
            static_cast<std::size_t>(
                global_row_end_ -
                global_row_start_) !=
                rows_.size()) {
            throw std::invalid_argument(
                "mpmc::discretization_petsc::OwnedCellStructuralColumnPatternSnapshot3D: invalid PETSc row ownership metadata");
        }

        std::size_t expected_diagonal_offset = 0U;
        std::size_t expected_off_diagonal_offset = 0U;

        for (std::size_t index = 0U;
             index < rows_.size();
             ++index) {
            const auto& pattern =
                rows_[index];
            const std::size_t local =
                static_cast<std::size_t>(
                    pattern.cell.value());
            const PetscInt expected_row =
                global_row_start_ +
                static_cast<PetscInt>(
                    index);

            if (local >= local_cell_count_ ||
                local_cell_to_row_[local]
                    .has_value() ||
                pattern.global_row !=
                    expected_row ||
                pattern.global_row < 0 ||
                pattern.global_row >=
                    global_row_count_ ||
                pattern.diagonal_column_offset !=
                    expected_diagonal_offset ||
                pattern.off_diagonal_column_offset !=
                    expected_off_diagonal_offset ||
                pattern.diagonal_column_count == 0U ||
                pattern.diagonal_column_count >
                    diagonal_global_columns_.size() -
                        pattern.diagonal_column_offset ||
                pattern.off_diagonal_column_count >
                    off_diagonal_global_columns_.size() -
                        pattern.off_diagonal_column_offset) {
                throw std::invalid_argument(
                    "mpmc::discretization_petsc::OwnedCellStructuralColumnPatternSnapshot3D: invalid row identity or compact column range");
            }

            const auto diagonal_columns =
                std::span<const PetscInt>{
                    diagonal_global_columns_}
                    .subspan(
                        pattern.diagonal_column_offset,
                        pattern.diagonal_column_count);
            const auto off_diagonal_columns =
                std::span<const PetscInt>{
                    off_diagonal_global_columns_}
                    .subspan(
                        pattern.off_diagonal_column_offset,
                        pattern.off_diagonal_column_count);

            if (!std::is_sorted(
                    diagonal_columns.begin(),
                    diagonal_columns.end()) ||
                std::adjacent_find(
                    diagonal_columns.begin(),
                    diagonal_columns.end()) !=
                    diagonal_columns.end() ||
                !std::binary_search(
                    diagonal_columns.begin(),
                    diagonal_columns.end(),
                    pattern.global_row)) {
                throw std::invalid_argument(
                    "mpmc::discretization_petsc::OwnedCellStructuralColumnPatternSnapshot3D: diagonal columns must be sorted, unique, and include the self row");
            }

            for (const PetscInt column :
                 diagonal_columns) {
                if (column <
                        global_row_start_ ||
                    column >=
                        global_row_end_) {
                    throw std::invalid_argument(
                        "mpmc::discretization_petsc::OwnedCellStructuralColumnPatternSnapshot3D: diagonal column is outside local PETSc ownership range");
                }
            }

            if (!std::is_sorted(
                    off_diagonal_columns.begin(),
                    off_diagonal_columns.end()) ||
                std::adjacent_find(
                    off_diagonal_columns.begin(),
                    off_diagonal_columns.end()) !=
                    off_diagonal_columns.end()) {
                throw std::invalid_argument(
                    "mpmc::discretization_petsc::OwnedCellStructuralColumnPatternSnapshot3D: off-diagonal columns must be sorted and unique");
            }

            for (const PetscInt column :
                 off_diagonal_columns) {
                if (column < 0 ||
                    column >=
                        global_row_count_ ||
                    (column >=
                         global_row_start_ &&
                     column <
                         global_row_end_)) {
                    throw std::invalid_argument(
                        "mpmc::discretization_petsc::OwnedCellStructuralColumnPatternSnapshot3D: off-diagonal column must be remote-owned");
                }
            }

            if (index > 0U &&
                !(rows_[index - 1U]
                      .cell_global <
                  pattern.cell_global)) {
                throw std::invalid_argument(
                    "mpmc::discretization_petsc::OwnedCellStructuralColumnPatternSnapshot3D: owned stable cell IDs must follow PETSc row order");
            }

            expected_diagonal_offset +=
                pattern.diagonal_column_count;
            expected_off_diagonal_offset +=
                pattern.off_diagonal_column_count;
            local_cell_to_row_[local] =
                index;
        }

        if (expected_diagonal_offset !=
                diagonal_global_columns_.size() ||
            expected_off_diagonal_offset !=
                off_diagonal_global_columns_.size()) {
            throw std::invalid_argument(
                "mpmc::discretization_petsc::OwnedCellStructuralColumnPatternSnapshot3D: compact column storage has unreferenced tail data");
        }
    }

    mpmc::mesh::PartitionRank local_rank_;
    std::uint32_t rank_count_;
    std::size_t local_cell_count_;
    PetscInt global_row_start_;
    PetscInt global_row_end_;
    PetscInt global_row_count_;
    std::vector<
        OwnedCellStructuralColumnPatternRow3D>
        rows_;
    std::vector<PetscInt>
        diagonal_global_columns_;
    std::vector<PetscInt>
        off_diagonal_global_columns_;
    std::vector<std::optional<std::size_t>>
        local_cell_to_row_;
};

/// Materialize exact structural PETSc global columns for each locally owned
/// cell row from the existing sparsity snapshot and symbolic row map.
///
/// self is always a diagonal-block column. Every unique sparsity neighbour is
/// mapped through bridge.global_row(); locally owned neighbours enter the
/// diagonal block and ghost/remote-owned neighbours enter the off-diagonal
/// block. Each row's column spans are sorted and deduplicated before being
/// appended to compact flat storage.
///
/// The final row cardinalities must exactly equal the frozen d_nnz/o_nnz
/// counts. No Mat is created or modified and no transmissibility/value/physics
/// quantity is consumed.
inline PetscErrorCode
make_owned_cell_structural_column_pattern_snapshot_3d(
    const CellPairSparsityStencilSnapshot3D& sparsity,
    const PetscMpiAijSymbolicPreallocation3D& bridge,
    const mpmc::mesh::PartitionSnapshot& partition,
    std::optional<
        OwnedCellStructuralColumnPatternSnapshot3D>*
        output) {
    if (output == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    output->reset();

    if (sparsity.local_rank() !=
            partition.local_rank() ||
        bridge.local_rank() !=
            partition.local_rank() ||
        sparsity.rank_count() !=
            partition.rank_count() ||
        bridge.rank_count() !=
            partition.rank_count() ||
        sparsity.local_cell_count() !=
            partition.entity_count(
                mpmc::mesh::EntityKind::cell) ||
        bridge.local_owned_row_count() < 0 ||
        static_cast<std::size_t>(
            bridge.local_owned_row_count()) !=
            sparsity.owned_cell_count() ||
        sparsity.owned_cell_count() !=
            partition.owned_count(
                mpmc::mesh::EntityKind::cell)) {
        return PETSC_ERR_ARG_INCOMP;
    }

    const auto owned_cells =
        bridge
            .owned_cells_in_petsc_row_order();
    const auto owned_ids =
        bridge.owned_cell_global_ids();
    const auto owned_rows =
        bridge.owned_global_rows();
    const auto diagonal_nnz =
        bridge.diagonal_nnz();
    const auto off_diagonal_nnz =
        bridge.off_diagonal_nnz();

    const std::size_t owned_count =
        owned_cells.size();
    if (owned_ids.size() != owned_count ||
        owned_rows.size() != owned_count ||
        diagonal_nnz.size() != owned_count ||
        off_diagonal_nnz.size() !=
            owned_count) {
        return PETSC_ERR_ARG_INCOMP;
    }

    std::size_t total_diagonal_columns = 0U;
    std::size_t total_off_diagonal_columns = 0U;
    for (std::size_t index = 0U;
         index < owned_count;
         ++index) {
        if (diagonal_nnz[index] <= 0 ||
            off_diagonal_nnz[index] < 0) {
            return PETSC_ERR_ARG_INCOMP;
        }
        const auto diagonal_count =
            static_cast<std::size_t>(
                diagonal_nnz[index]);
        const auto off_diagonal_count =
            static_cast<std::size_t>(
                off_diagonal_nnz[index]);
        if (diagonal_count >
                std::numeric_limits<std::size_t>::max() -
                    total_diagonal_columns ||
            off_diagonal_count >
                std::numeric_limits<std::size_t>::max() -
                    total_off_diagonal_columns) {
            return PETSC_ERR_ARG_OUTOFRANGE;
        }
        total_diagonal_columns +=
            diagonal_count;
        total_off_diagonal_columns +=
            off_diagonal_count;
    }

    std::vector<
        OwnedCellStructuralColumnPatternRow3D>
        patterns;
    std::vector<PetscInt>
        flat_diagonal_columns;
    std::vector<PetscInt>
        flat_off_diagonal_columns;
    patterns.reserve(owned_count);
    flat_diagonal_columns.reserve(
        total_diagonal_columns);
    flat_off_diagonal_columns.reserve(
        total_off_diagonal_columns);

    for (std::size_t index = 0U;
         index < owned_count;
         ++index) {
        const auto cell =
            owned_cells[index];
        const auto cell_global =
            owned_ids[index];
        const PetscInt global_row =
            owned_rows[index];

        try {
            if (!partition.is_owned(
                    mpmc::mesh::EntityKind::cell,
                    cell) ||
                partition.global_id(
                    mpmc::mesh::EntityKind::cell,
                    cell) !=
                    cell_global ||
                bridge.global_row(cell) !=
                    global_row) {
                return PETSC_ERR_ARG_INCOMP;
            }
        } catch (...) {
            return PETSC_ERR_ARG_INCOMP;
        }

        std::vector<PetscInt>
            diagonal_columns{
                global_row};
        std::vector<PetscInt>
            off_diagonal_columns;

        for (const auto& coupling :
             sparsity.couplings()) {
            try {
                if (partition.global_id(
                        mpmc::mesh::EntityKind::cell,
                        coupling.first_cell) !=
                        coupling.first_cell_global ||
                    partition.global_id(
                        mpmc::mesh::EntityKind::cell,
                        coupling.second_cell) !=
                        coupling.second_cell_global) {
                    return PETSC_ERR_ARG_INCOMP;
                }
            } catch (...) {
                return PETSC_ERR_ARG_INCOMP;
            }

            mpmc::mesh::LocalIndex neighbour{
                0U};
            bool incident = false;
            if (coupling.first_cell ==
                cell) {
                neighbour =
                    coupling.second_cell;
                incident = true;
            } else if (
                coupling.second_cell ==
                cell) {
                neighbour =
                    coupling.first_cell;
                incident = true;
            }
            if (!incident) {
                continue;
            }

            PetscInt neighbour_row = -1;
            try {
                neighbour_row =
                    bridge.global_row(
                        neighbour);
                if (partition.is_owned(
                        mpmc::mesh::EntityKind::cell,
                        neighbour)) {
                    diagonal_columns.push_back(
                        neighbour_row);
                } else if (
                    partition.is_ghost(
                        mpmc::mesh::EntityKind::cell,
                        neighbour)) {
                    off_diagonal_columns.push_back(
                        neighbour_row);
                } else {
                    return PETSC_ERR_ARG_INCOMP;
                }
            } catch (...) {
                return PETSC_ERR_ARG_INCOMP;
            }
        }

        std::sort(
            diagonal_columns.begin(),
            diagonal_columns.end());
        diagonal_columns.erase(
            std::unique(
                diagonal_columns.begin(),
                diagonal_columns.end()),
            diagonal_columns.end());

        std::sort(
            off_diagonal_columns.begin(),
            off_diagonal_columns.end());
        off_diagonal_columns.erase(
            std::unique(
                off_diagonal_columns.begin(),
                off_diagonal_columns.end()),
            off_diagonal_columns.end());

        const OwnedCellStructuralCounts3D*
            structural = nullptr;
        try {
            structural =
                &sparsity.structural_counts(
                    cell);
        } catch (...) {
            return PETSC_ERR_ARG_INCOMP;
        }
        if (structural->cell !=
                cell ||
            structural->cell_global !=
                cell_global ||
            diagonal_columns.size() !=
                structural
                    ->diagonal_block_nnz ||
            off_diagonal_columns.size() !=
                structural
                    ->off_diagonal_block_nnz ||
            diagonal_columns.size() !=
                static_cast<std::size_t>(
                    diagonal_nnz[index]) ||
            off_diagonal_columns.size() !=
                static_cast<std::size_t>(
                    off_diagonal_nnz[index]) ||
            diagonal_columns.size() +
                    off_diagonal_columns.size() !=
                static_cast<std::size_t>(
                    diagonal_nnz[index] +
                    off_diagonal_nnz[index])) {
            return PETSC_ERR_ARG_INCOMP;
        }

        const std::size_t diagonal_offset =
            flat_diagonal_columns.size();
        const std::size_t off_diagonal_offset =
            flat_off_diagonal_columns.size();

        flat_diagonal_columns.insert(
            flat_diagonal_columns.end(),
            diagonal_columns.begin(),
            diagonal_columns.end());
        flat_off_diagonal_columns.insert(
            flat_off_diagonal_columns.end(),
            off_diagonal_columns.begin(),
            off_diagonal_columns.end());

        patterns.push_back(
            OwnedCellStructuralColumnPatternRow3D{
                cell,
                cell_global,
                global_row,
                diagonal_offset,
                diagonal_columns.size(),
                off_diagonal_offset,
                off_diagonal_columns.size()});
    }

    if (flat_diagonal_columns.size() !=
            total_diagonal_columns ||
        flat_off_diagonal_columns.size() !=
            total_off_diagonal_columns) {
        return PETSC_ERR_PLIB;
    }

    try {
        output->emplace(
            partition.local_rank(),
            partition.rank_count(),
            partition.entity_count(
                mpmc::mesh::EntityKind::cell),
            bridge.global_row_start(),
            bridge.global_row_end(),
            bridge.global_row_count(),
            std::move(patterns),
            std::move(flat_diagonal_columns),
            std::move(flat_off_diagonal_columns));
    } catch (...) {
        output->reset();
        return PETSC_ERR_ARG_INCOMP;
    }

    return PETSC_SUCCESS;
}


/// Create an empty square MPIAIJ matrix using only the frozen symbolic
/// preallocation bridge.
///
/// This gate materializes PETSc row/column ownership and calls
/// MatMPIAIJSetPreallocation() with the bridge d_nnz/o_nnz arrays. It does not
/// insert any matrix values, define any physical unknown, or consume a
/// transmissibility/flux/residual quantity. The returned Mat remains value-free
/// and is owned by the caller.
inline PetscErrorCode
create_empty_petsc_mpiaij_symbolic_matrix_3d(
    MPI_Comm comm,
    const PetscMpiAijSymbolicPreallocation3D& bridge,
    Mat* matrix) {
    if (matrix == nullptr) {
        return PETSC_ERR_ARG_NULL;
    }
    if (*matrix != nullptr) {
        return PETSC_ERR_ARG_WRONGSTATE;
    }

    PetscErrorCode error =
        detail::validate_communicator(
            comm,
            bridge.local_rank(),
            bridge.rank_count());
    if (error != PETSC_SUCCESS) return error;

    const PetscInt local_rows =
        bridge.local_owned_row_count();
    const PetscInt global_rows =
        bridge.global_row_count();
    if (local_rows < 0 ||
        global_rows < local_rows) {
        return PETSC_ERR_ARG_INCOMP;
    }

    Mat created = nullptr;
    error = MatCreate(
        comm,
        &created);
    if (error == PETSC_SUCCESS) {
        error = MatSetSizes(
            created,
            local_rows,
            local_rows,
            global_rows,
            global_rows);
    }
    if (error == PETSC_SUCCESS) {
        error = MatSetType(
            created,
            MATMPIAIJ);
    }
    if (error == PETSC_SUCCESS) {
        const auto diagonal_nnz =
            bridge.diagonal_nnz();
        const auto off_diagonal_nnz =
            bridge.off_diagonal_nnz();
        error = MatMPIAIJSetPreallocation(
            created,
            0,
            diagonal_nnz.empty()
                ? nullptr
                : diagonal_nnz.data(),
            0,
            off_diagonal_nnz.empty()
                ? nullptr
                : off_diagonal_nnz.data());
    }
    if (error != PETSC_SUCCESS) {
        MatDestroy(&created);
        return error;
    }

    PetscBool is_mpiaij = PETSC_FALSE;
    error = PetscObjectTypeCompare(
        reinterpret_cast<PetscObject>(
            created),
        MATMPIAIJ,
        &is_mpiaij);

    PetscInt actual_local_rows = -1;
    PetscInt actual_local_columns = -1;
    PetscInt actual_global_rows = -1;
    PetscInt actual_global_columns = -1;
    PetscInt row_start = -1;
    PetscInt row_end = -1;
    PetscInt column_start = -1;
    PetscInt column_end = -1;

    if (error == PETSC_SUCCESS) {
        error = MatGetLocalSize(
            created,
            &actual_local_rows,
            &actual_local_columns);
    }
    if (error == PETSC_SUCCESS) {
        error = MatGetSize(
            created,
            &actual_global_rows,
            &actual_global_columns);
    }
    if (error == PETSC_SUCCESS) {
        error = MatGetOwnershipRange(
            created,
            &row_start,
            &row_end);
    }
    if (error == PETSC_SUCCESS) {
        error = MatGetOwnershipRangeColumn(
            created,
            &column_start,
            &column_end);
    }
    if (error != PETSC_SUCCESS) {
        MatDestroy(&created);
        return error;
    }

    if (is_mpiaij != PETSC_TRUE ||
        actual_local_rows != local_rows ||
        actual_local_columns != local_rows ||
        actual_global_rows != global_rows ||
        actual_global_columns != global_rows ||
        row_start != bridge.global_row_start() ||
        row_end != bridge.global_row_end() ||
        column_start != row_start ||
        column_end != row_end) {
        MatDestroy(&created);
        return PETSC_ERR_ARG_INCOMP;
    }

    const auto owned_cells =
        bridge.owned_cells_in_petsc_row_order();
    const auto owned_global_rows =
        bridge.owned_global_rows();
    const auto owned_global_ids =
        bridge.owned_cell_global_ids();
    if (owned_cells.size() !=
            owned_global_rows.size() ||
        owned_cells.size() !=
            owned_global_ids.size() ||
        owned_cells.size() !=
            static_cast<std::size_t>(
                local_rows)) {
        MatDestroy(&created);
        return PETSC_ERR_ARG_INCOMP;
    }

    for (std::size_t index = 0U;
         index < owned_cells.size();
         ++index) {
        const PetscInt expected_row =
            row_start +
            static_cast<PetscInt>(
                index);
        if (owned_global_rows[index] !=
                expected_row ||
            bridge.global_row(
                owned_cells[index]) !=
                expected_row ||
            (index > 0U &&
             !(owned_global_ids[index - 1U] <
               owned_global_ids[index]))) {
            MatDestroy(&created);
            return PETSC_ERR_ARG_INCOMP;
        }
    }

    *matrix = created;
    return PETSC_SUCCESS;
}


} // namespace mpmc::discretization_petsc

#endif // MPMC_DISCRETIZATION_PETSC_ADAPTER_HPP
