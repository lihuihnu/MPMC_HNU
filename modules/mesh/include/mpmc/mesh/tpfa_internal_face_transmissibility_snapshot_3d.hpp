#ifndef MPMC_MESH_TPFA_INTERNAL_FACE_TRANSMISSIBILITY_SNAPSHOT_3D_HPP
#define MPMC_MESH_TPFA_INTERNAL_FACE_TRANSMISSIBILITY_SNAPSHOT_3D_HPP

#include <mpmc/mesh/permeability_tensor_3d.hpp>
#include <mpmc/mesh/tpfa_static_face_transmissibility_3d.hpp>

#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

/// Gated state for one internal face.
///
/// Only materialized entries carry a static face transmissibility. Every blocked
/// state preserves why the upstream geometry/K admissibility gate refused direct
/// TPFA use.
enum class TpfaInternalFaceTransmissibilityDisposition3D {
    materialized,
    blocked_geometry_non_orthogonal,
    blocked_k_non_orthogonal,
    blocked_geometry_and_k_non_orthogonal,
    blocked_degenerate_permeability_direction
};

struct TpfaInternalFaceTransmissibilityEntry3D {
    LocalIndex face;
    TpfaInternalFaceTransmissibilityDisposition3D disposition;
    CombinedTransmissibilityAdmissibility3D admissibility;
    std::optional<TpfaStaticFaceTransmissibility3D>
        static_transmissibility;
};

/// Immutable snapshot of admissibility-gated static TPFA transmissibilities.
///
/// Boundary faces are deliberately absent from entries(). Every internal face
/// appears exactly once. A direct geometry/K candidate is materialized to one
/// static T_f [m3]; all other internal faces retain a blocked disposition and
/// expose no usable transmissibility.
///
/// The explicit geometry/K angular policies are retained with the snapshot so
/// the materialization decision remains auditable.
class TpfaInternalFaceTransmissibilitySnapshot3D {
public:
    TpfaInternalFaceTransmissibilitySnapshot3D(
        std::size_t total_face_count,
        TransmissibilityGeometryAdmissibilityPolicy3D
            geometry_policy,
        KOrthogonalityAdmissibilityPolicy3D
            k_policy,
        std::vector<TpfaInternalFaceTransmissibilityEntry3D>
            entries)
        : total_face_count_(total_face_count),
          geometry_policy_(geometry_policy),
          k_policy_(k_policy),
          entries_(std::move(entries)),
          face_to_entry_(total_face_count_) {
        validate_and_index();
    }

    TpfaInternalFaceTransmissibilitySnapshot3D(
        const TpfaInternalFaceTransmissibilitySnapshot3D&) = default;
    TpfaInternalFaceTransmissibilitySnapshot3D(
        TpfaInternalFaceTransmissibilitySnapshot3D&&) noexcept = default;
    TpfaInternalFaceTransmissibilitySnapshot3D& operator=(
        const TpfaInternalFaceTransmissibilitySnapshot3D&) = delete;
    TpfaInternalFaceTransmissibilitySnapshot3D& operator=(
        TpfaInternalFaceTransmissibilitySnapshot3D&&) = delete;
    ~TpfaInternalFaceTransmissibilitySnapshot3D() = default;

    [[nodiscard]] std::size_t total_face_count() const noexcept {
        return total_face_count_;
    }

    [[nodiscard]] std::size_t internal_face_count() const noexcept {
        return entries_.size();
    }

    [[nodiscard]] std::size_t materialized_face_count() const noexcept {
        return materialized_face_count_;
    }

    [[nodiscard]] std::size_t blocked_face_count() const noexcept {
        return entries_.size() - materialized_face_count_;
    }

    [[nodiscard]] TransmissibilityGeometryAdmissibilityPolicy3D
    geometry_policy() const noexcept {
        return geometry_policy_;
    }

    [[nodiscard]] KOrthogonalityAdmissibilityPolicy3D
    k_policy() const noexcept {
        return k_policy_;
    }

    [[nodiscard]] std::span<
        const TpfaInternalFaceTransmissibilityEntry3D>
    entries() const noexcept {
        return entries_;
    }

    [[nodiscard]] bool contains_internal_face(
        LocalIndex face) const {
        const std::size_t local =
            static_cast<std::size_t>(face.value());
        if (local >= total_face_count_) {
            throw std::out_of_range(
                "mpmc::mesh::TpfaInternalFaceTransmissibilitySnapshot3D: face index out of range");
        }
        return face_to_entry_[local].has_value();
    }

    [[nodiscard]] const TpfaInternalFaceTransmissibilityEntry3D&
    entry(LocalIndex face) const {
        const std::size_t local =
            static_cast<std::size_t>(face.value());
        if (local >= total_face_count_) {
            throw std::out_of_range(
                "mpmc::mesh::TpfaInternalFaceTransmissibilitySnapshot3D: face index out of range");
        }
        const auto mapped =
            face_to_entry_[local];
        if (!mapped.has_value()) {
            throw std::invalid_argument(
                "mpmc::mesh::TpfaInternalFaceTransmissibilitySnapshot3D: requested face is not an internal face in this snapshot");
        }
        return entries_.at(*mapped);
    }

private:
    void validate_and_index() {
        materialized_face_count_ = 0U;

        for (std::size_t index = 0U;
             index < entries_.size();
             ++index) {
            const auto& entry = entries_[index];
            const std::size_t face =
                static_cast<std::size_t>(
                    entry.face.value());
            if (face >= total_face_count_) {
                throw std::out_of_range(
                    "mpmc::mesh::TpfaInternalFaceTransmissibilitySnapshot3D: internal face index out of range");
            }
            if (face_to_entry_[face].has_value()) {
                throw std::invalid_argument(
                    "mpmc::mesh::TpfaInternalFaceTransmissibilitySnapshot3D: duplicate internal face entry");
            }

            switch (entry.disposition) {
            case TpfaInternalFaceTransmissibilityDisposition3D::
                materialized:
                if (entry.admissibility.disposition !=
                        CombinedTransmissibilityAdmissibilityDisposition3D::
                            direct_normal_projection_k_orthogonal_candidate ||
                    !entry.static_transmissibility.has_value() ||
                    entry.static_transmissibility->disposition !=
                        TpfaStaticFaceTransmissibilityDisposition3D::
                            positive_harmonic_combination ||
                    entry.static_transmissibility
                            ->face_transmissibility_m3 <=
                        0.0) {
                    throw std::invalid_argument(
                        "mpmc::mesh::TpfaInternalFaceTransmissibilitySnapshot3D: materialized entry must be a direct candidate with positive static transmissibility");
                }
                ++materialized_face_count_;
                break;
            case TpfaInternalFaceTransmissibilityDisposition3D::
                blocked_geometry_non_orthogonal:
                require_blocked(
                    entry,
                    CombinedTransmissibilityAdmissibilityDisposition3D::
                        requires_geometry_non_orthogonal_treatment);
                break;
            case TpfaInternalFaceTransmissibilityDisposition3D::
                blocked_k_non_orthogonal:
                require_blocked(
                    entry,
                    CombinedTransmissibilityAdmissibilityDisposition3D::
                        requires_k_non_orthogonal_treatment);
                break;
            case TpfaInternalFaceTransmissibilityDisposition3D::
                blocked_geometry_and_k_non_orthogonal:
                require_blocked(
                    entry,
                    CombinedTransmissibilityAdmissibilityDisposition3D::
                        requires_geometry_and_k_non_orthogonal_treatment);
                break;
            case TpfaInternalFaceTransmissibilityDisposition3D::
                blocked_degenerate_permeability_direction:
                require_blocked(
                    entry,
                    CombinedTransmissibilityAdmissibilityDisposition3D::
                        degenerate_permeability_direction);
                break;
            default:
                throw std::invalid_argument(
                    "mpmc::mesh::TpfaInternalFaceTransmissibilitySnapshot3D: invalid internal face disposition");
            }

            face_to_entry_[face] = index;
        }
    }

    static void require_blocked(
        const TpfaInternalFaceTransmissibilityEntry3D& entry,
        CombinedTransmissibilityAdmissibilityDisposition3D
            expected) {
        if (entry.admissibility.disposition != expected ||
            entry.static_transmissibility.has_value()) {
            throw std::invalid_argument(
                "mpmc::mesh::TpfaInternalFaceTransmissibilitySnapshot3D: blocked entry must preserve matching admissibility reason and expose no transmissibility");
        }
    }

    std::size_t total_face_count_;
    TransmissibilityGeometryAdmissibilityPolicy3D
        geometry_policy_;
    KOrthogonalityAdmissibilityPolicy3D
        k_policy_;
    std::vector<TpfaInternalFaceTransmissibilityEntry3D>
        entries_;
    std::vector<std::optional<std::size_t>>
        face_to_entry_;
    std::size_t materialized_face_count_{0U};
};

[[nodiscard]] inline TpfaInternalFaceTransmissibilityDisposition3D
gated_disposition(
    CombinedTransmissibilityAdmissibilityDisposition3D
        disposition) {
    switch (disposition) {
    case CombinedTransmissibilityAdmissibilityDisposition3D::
        direct_normal_projection_k_orthogonal_candidate:
        return TpfaInternalFaceTransmissibilityDisposition3D::
            materialized;
    case CombinedTransmissibilityAdmissibilityDisposition3D::
        requires_geometry_non_orthogonal_treatment:
        return TpfaInternalFaceTransmissibilityDisposition3D::
            blocked_geometry_non_orthogonal;
    case CombinedTransmissibilityAdmissibilityDisposition3D::
        requires_k_non_orthogonal_treatment:
        return TpfaInternalFaceTransmissibilityDisposition3D::
            blocked_k_non_orthogonal;
    case CombinedTransmissibilityAdmissibilityDisposition3D::
        requires_geometry_and_k_non_orthogonal_treatment:
        return TpfaInternalFaceTransmissibilityDisposition3D::
            blocked_geometry_and_k_non_orthogonal;
    case CombinedTransmissibilityAdmissibilityDisposition3D::
        degenerate_permeability_direction:
        return TpfaInternalFaceTransmissibilityDisposition3D::
            blocked_degenerate_permeability_direction;
    }
    throw std::invalid_argument(
        "mpmc::mesh::gated_disposition: invalid combined admissibility disposition");
}

/// Build an internal-face snapshot using the supplied explicit admissibility
/// policies.
///
/// Crucially, static transmissibility construction is invoked only for direct
/// geometry/K candidates. Blocked faces therefore never manufacture or expose
/// a T_f that a downstream flow layer could accidentally consume.
[[nodiscard]] inline TpfaInternalFaceTransmissibilitySnapshot3D
make_admissibility_gated_internal_face_transmissibility_snapshot_3d(
    const CellFaceGeometricOperator3D& geometry,
    const CellCartesianDiagonalPermeability3D& permeability,
    TransmissibilityGeometryAdmissibilityPolicy3D
        geometry_policy,
    KOrthogonalityAdmissibilityPolicy3D
        k_policy) {
    if (geometry.cell_count() !=
        permeability.cell_count()) {
        throw std::invalid_argument(
            "mpmc::mesh::make_admissibility_gated_internal_face_transmissibility_snapshot_3d: geometry/permeability cell counts do not match");
    }

    std::vector<TpfaInternalFaceTransmissibilityEntry3D>
        entries;
    entries.reserve(geometry.face_count());

    for (std::size_t face = 0U;
         face < geometry.face_count();
         ++face) {
        const LocalIndex local{
            static_cast<LocalIndex::value_type>(
                face)};
        if (!geometry.face_neighbour(local).has_value()) {
            continue;
        }

        const auto admissibility =
            classify_internal_face_transmissibility_admissibility(
                geometry,
                permeability,
                local,
                geometry_policy,
                k_policy);
        const auto disposition =
            gated_disposition(
                admissibility.disposition);

        std::optional<TpfaStaticFaceTransmissibility3D>
            static_transmissibility;
        if (disposition ==
            TpfaInternalFaceTransmissibilityDisposition3D::
                materialized) {
            static_transmissibility =
                make_internal_face_static_tpfa_transmissibility_3d(
                    geometry,
                    permeability,
                    local);
        }

        entries.push_back(
            TpfaInternalFaceTransmissibilityEntry3D{
                local,
                disposition,
                admissibility,
                static_transmissibility});
    }

    return TpfaInternalFaceTransmissibilitySnapshot3D{
        geometry.face_count(),
        geometry_policy,
        k_policy,
        std::move(entries)};
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_TPFA_INTERNAL_FACE_TRANSMISSIBILITY_SNAPSHOT_3D_HPP
