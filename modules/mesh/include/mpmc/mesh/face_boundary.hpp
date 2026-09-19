#ifndef MPMC_MESH_FACE_BOUNDARY_HPP
#define MPMC_MESH_FACE_BOUNDARY_HPP

#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/topology.hpp>

#include <compare>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

enum class FaceClassification : std::uint8_t {
    interior = 0,
    boundary = 1,
};

/// Stable integer identity for a boundary/physical group.
///
/// Value zero is reserved for "untagged". Positive values are opaque identifiers:
/// this layer assigns no pressure, flux, wall, well, material or other physics.
class PhysicalTag {
public:
    using value_type = std::uint32_t;
    static constexpr value_type untagged_value = 0U;

    explicit constexpr PhysicalTag(value_type value) noexcept : value_(value) {}

    [[nodiscard]] constexpr value_type value() const noexcept { return value_; }
    [[nodiscard]] constexpr bool is_tagged() const noexcept {
        return value_ != untagged_value;
    }

    [[nodiscard]] friend constexpr auto operator<=>(const PhysicalTag&,
                                                     const PhysicalTag&) noexcept = default;

private:
    value_type value_;
};

/// Immutable per-face classification and physical-tag snapshot.
///
/// The arrays align exactly with Topology face local ordering. Classification is
/// derived from face->cell cardinality: one adjacent cell is boundary, two is
/// interior. A boundary face may be untagged; an interior face must be untagged.
class FaceBoundarySnapshot {
public:
    FaceBoundarySnapshot(std::vector<FaceClassification> classifications,
                         std::vector<PhysicalTag> physical_tags)
        : classifications_(std::move(classifications)),
          physical_tags_(std::move(physical_tags)) {
        validate();
    }

    FaceBoundarySnapshot(const FaceBoundarySnapshot&) = default;
    FaceBoundarySnapshot(FaceBoundarySnapshot&&) noexcept = default;
    FaceBoundarySnapshot& operator=(const FaceBoundarySnapshot&) = delete;
    FaceBoundarySnapshot& operator=(FaceBoundarySnapshot&&) = delete;
    ~FaceBoundarySnapshot() = default;

    [[nodiscard]] std::size_t face_count() const noexcept {
        return classifications_.size();
    }
    [[nodiscard]] std::size_t boundary_face_count() const noexcept {
        return boundary_face_count_;
    }
    [[nodiscard]] std::size_t interior_face_count() const noexcept {
        return classifications_.size() - boundary_face_count_;
    }

    [[nodiscard]] std::span<const FaceClassification> classifications() const noexcept {
        return classifications_;
    }
    [[nodiscard]] std::span<const PhysicalTag> physical_tags() const noexcept {
        return physical_tags_;
    }

    [[nodiscard]] FaceClassification classification(LocalIndex face) const {
        return classifications_.at(static_cast<std::size_t>(face.value()));
    }
    [[nodiscard]] bool is_boundary(LocalIndex face) const {
        return classification(face) == FaceClassification::boundary;
    }
    [[nodiscard]] PhysicalTag physical_tag(LocalIndex face) const {
        return physical_tags_.at(static_cast<std::size_t>(face.value()));
    }
    [[nodiscard]] bool has_physical_tag(LocalIndex face) const {
        return physical_tag(face).is_tagged();
    }

private:
    void validate() {
        if (classifications_.size() != physical_tags_.size()) {
            throw std::invalid_argument(
                "mpmc::mesh::FaceBoundarySnapshot: classification/tag array sizes do not match");
        }

        boundary_face_count_ = 0U;
        for (std::size_t face = 0; face < classifications_.size(); ++face) {
            switch (classifications_[face]) {
            case FaceClassification::interior:
                if (physical_tags_[face].is_tagged()) {
                    throw std::invalid_argument(
                        "mpmc::mesh::FaceBoundarySnapshot: interior face cannot carry a physical tag");
                }
                break;
            case FaceClassification::boundary:
                ++boundary_face_count_;
                break;
            default:
                throw std::invalid_argument(
                    "mpmc::mesh::FaceBoundarySnapshot: invalid face classification");
            }
        }
    }

    std::vector<FaceClassification> classifications_;
    std::vector<PhysicalTag> physical_tags_;
    std::size_t boundary_face_count_{0U};
};

/// Classify all faces from Topology and attach one optional integer tag per face.
///
/// physical_tags must have exactly Topology::face count entries. Zero means
/// untagged. Nonzero tags are accepted only on boundary faces.
[[nodiscard]] inline FaceBoundarySnapshot make_face_boundary_snapshot(
    const Topology& topology,
    std::span<const PhysicalTag> physical_tags) {
    const std::size_t face_count = topology.entity_count(EntityKind::face);
    if (physical_tags.size() != face_count) {
        throw std::invalid_argument(
            "mpmc::mesh::make_face_boundary_snapshot: physical tag count does not match faces");
    }
    if (!topology.has_relation(EntityKind::face, EntityKind::cell)) {
        throw std::invalid_argument(
            "mpmc::mesh::make_face_boundary_snapshot: face-to-cell relation is required");
    }

    const auto& face_cell = topology.relation(EntityKind::face, EntityKind::cell);
    std::vector<FaceClassification> classifications;
    classifications.reserve(face_count);
    std::vector<PhysicalTag> tags;
    tags.reserve(face_count);

    for (std::size_t face = 0; face < face_count; ++face) {
        const auto adjacent =
            face_cell.adjacent(LocalIndex{static_cast<LocalIndex::value_type>(face)});
        if (adjacent.size() == 1U) {
            classifications.push_back(FaceClassification::boundary);
        } else if (adjacent.size() == 2U) {
            classifications.push_back(FaceClassification::interior);
            if (physical_tags[face].is_tagged()) {
                throw std::invalid_argument(
                    "mpmc::mesh::make_face_boundary_snapshot: interior face cannot carry a physical tag");
            }
        } else {
            throw std::invalid_argument(
                "mpmc::mesh::make_face_boundary_snapshot: each face must have one or two adjacent cells");
        }
        tags.push_back(physical_tags[face]);
    }

    return FaceBoundarySnapshot{std::move(classifications), std::move(tags)};
}

/// Convenience classification with every boundary face initially untagged.
[[nodiscard]] inline FaceBoundarySnapshot make_face_boundary_snapshot(
    const Topology& topology) {
    const std::size_t face_count = topology.entity_count(EntityKind::face);
    const std::vector<PhysicalTag> untagged(face_count, PhysicalTag{0U});
    return make_face_boundary_snapshot(topology, untagged);
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_FACE_BOUNDARY_HPP
