#ifndef MPMC_MESH_CSR_ADJACENCY_HPP
#define MPMC_MESH_CSR_ADJACENCY_HPP

#include <mpmc/mesh/entity.hpp>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

/// Compact immutable-after-construction adjacency relation.
///
/// Rows are source entities, values are target LocalIndex values. Storage consists
/// of exactly two contiguous arrays: CSR offsets and target indices. The relation
/// owns no per-row heap objects and performs no allocation during indexed reads.
class CsrAdjacency {
public:
    using Offset = std::uint32_t;

    CsrAdjacency(EntityKind source_kind,
                 EntityKind target_kind,
                 std::size_t target_count,
                 std::vector<Offset> offsets,
                 std::vector<LocalIndex> indices)
        : source_kind_(source_kind),
          target_kind_(target_kind),
          target_count_(target_count),
          offsets_(std::move(offsets)),
          indices_(std::move(indices)) {
        validate();
    }

    [[nodiscard]] EntityKind source_kind() const noexcept { return source_kind_; }
    [[nodiscard]] EntityKind target_kind() const noexcept { return target_kind_; }
    [[nodiscard]] std::size_t source_count() const noexcept { return offsets_.size() - 1U; }
    [[nodiscard]] std::size_t target_count() const noexcept { return target_count_; }
    [[nodiscard]] std::size_t entry_count() const noexcept { return indices_.size(); }

    [[nodiscard]] std::span<const Offset> offsets() const noexcept { return offsets_; }
    [[nodiscard]] std::span<const LocalIndex> indices() const noexcept { return indices_; }

    /// Bounds-checked source-row lookup. Returned storage aliases this object.
    [[nodiscard]] std::span<const LocalIndex> adjacent(LocalIndex source) const {
        const std::size_t row = static_cast<std::size_t>(source.value());
        if (row >= source_count()) {
            throw std::out_of_range("mpmc::mesh::CsrAdjacency: source index out of range");
        }
        const std::size_t begin = static_cast<std::size_t>(offsets_[row]);
        const std::size_t end = static_cast<std::size_t>(offsets_[row + 1U]);
        return std::span<const LocalIndex>{indices_}.subspan(begin, end - begin);
    }

private:
    void validate() const {
        if (offsets_.empty()) {
            throw std::invalid_argument("mpmc::mesh::CsrAdjacency: offsets must contain row zero");
        }
        if (offsets_.front() != Offset{0}) {
            throw std::invalid_argument("mpmc::mesh::CsrAdjacency: first offset must be zero");
        }

        constexpr auto local_max = std::numeric_limits<LocalIndex::value_type>::max();
        const std::uint64_t local_capacity = static_cast<std::uint64_t>(local_max) + 1ULL;
        if (static_cast<std::uint64_t>(source_count()) > local_capacity ||
            static_cast<std::uint64_t>(target_count_) > local_capacity) {
            throw std::length_error(
                "mpmc::mesh::CsrAdjacency: local entity count exceeds LocalIndex capacity");
        }
        if (indices_.size() > static_cast<std::size_t>(std::numeric_limits<Offset>::max())) {
            throw std::length_error(
                "mpmc::mesh::CsrAdjacency: adjacency entries exceed CSR offset capacity");
        }

        Offset previous = Offset{0};
        for (const Offset offset : offsets_) {
            if (offset < previous || static_cast<std::size_t>(offset) > indices_.size()) {
                throw std::invalid_argument("mpmc::mesh::CsrAdjacency: invalid CSR offsets");
            }
            previous = offset;
        }
        if (static_cast<std::size_t>(offsets_.back()) != indices_.size()) {
            throw std::invalid_argument(
                "mpmc::mesh::CsrAdjacency: final offset does not match entry count");
        }

        for (const LocalIndex target : indices_) {
            if (static_cast<std::size_t>(target.value()) >= target_count_) {
                throw std::out_of_range("mpmc::mesh::CsrAdjacency: target index out of range");
            }
        }
    }

    EntityKind source_kind_;
    EntityKind target_kind_;
    std::size_t target_count_;
    std::vector<Offset> offsets_;
    std::vector<LocalIndex> indices_;
};

} // namespace mpmc::mesh

#endif // MPMC_MESH_CSR_ADJACENCY_HPP
