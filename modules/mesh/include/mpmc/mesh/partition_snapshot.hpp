#ifndef MPMC_MESH_PARTITION_SNAPSHOT_HPP
#define MPMC_MESH_PARTITION_SNAPSHOT_HPP

#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <array>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

class PartitionRank {
public:
    using value_type = std::uint32_t;

    explicit constexpr PartitionRank(value_type value) noexcept : value_(value) {}

    [[nodiscard]] constexpr value_type value() const noexcept { return value_; }

    [[nodiscard]] friend constexpr auto operator<=>(const PartitionRank&,
                                                     const PartitionRank&) noexcept = default;

private:
    value_type value_;
};

enum class EntityOwnership : std::uint8_t {
    owned = 0,
    ghost = 1,
};

/// Per-kind owner-rank input aligned to Topology local ordering.
struct EntityOwnerRanks {
    std::vector<PartitionRank> vertices;
    std::vector<PartitionRank> edges;
    std::vector<PartitionRank> faces;
    std::vector<PartitionRank> cells;
};

/// Immutable local partition/ownership snapshot.
///
/// Each local entity carries the stable GlobalEntityId copied from Topology and one
/// owner rank. Ownership is derived, never duplicated:
/// owner_rank == local_rank -> owned, otherwise ghost.
///
/// The snapshot owns all mapping data and retains no reference to Topology. Local
/// ordering is exactly the source Topology ordering. Global-to-local lookup uses a
/// compact LocalIndex permutation sorted by GlobalEntityId, not per-entity nodes.
class PartitionSnapshot {
public:
    PartitionSnapshot(const PartitionSnapshot&) = default;
    PartitionSnapshot(PartitionSnapshot&&) noexcept = default;
    PartitionSnapshot& operator=(const PartitionSnapshot&) = delete;
    PartitionSnapshot& operator=(PartitionSnapshot&&) = delete;
    ~PartitionSnapshot() = default;

    [[nodiscard]] static PartitionSnapshot create(
        const Topology& topology,
        PartitionRank local_rank,
        std::uint32_t rank_count,
        EntityOwnerRanks owner_ranks) {
        if (rank_count == 0U) {
            throw std::invalid_argument(
                "mpmc::mesh::PartitionSnapshot: rank_count must be positive");
        }
        if (local_rank.value() >= rank_count) {
            throw std::out_of_range(
                "mpmc::mesh::PartitionSnapshot: local_rank is outside rank_count");
        }

        std::array<KindPartition, 4> kinds{};
        initialize_kind(topology, EntityKind::vertex, local_rank, rank_count,
                        std::move(owner_ranks.vertices), kinds[0]);
        initialize_kind(topology, EntityKind::edge, local_rank, rank_count,
                        std::move(owner_ranks.edges), kinds[1]);
        initialize_kind(topology, EntityKind::face, local_rank, rank_count,
                        std::move(owner_ranks.faces), kinds[2]);
        initialize_kind(topology, EntityKind::cell, local_rank, rank_count,
                        std::move(owner_ranks.cells), kinds[3]);

        return PartitionSnapshot{local_rank, rank_count, std::move(kinds)};
    }

    [[nodiscard]] PartitionRank local_rank() const noexcept { return local_rank_; }
    [[nodiscard]] std::uint32_t rank_count() const noexcept { return rank_count_; }
    [[nodiscard]] bool is_serial() const noexcept { return rank_count_ == 1U; }

    [[nodiscard]] std::size_t entity_count(EntityKind kind) const {
        return kind_partition(kind).global_ids.size();
    }
    [[nodiscard]] std::size_t owned_count(EntityKind kind) const {
        return kind_partition(kind).owned_count;
    }
    [[nodiscard]] std::size_t ghost_count(EntityKind kind) const {
        const auto& partition = kind_partition(kind);
        return partition.global_ids.size() - partition.owned_count;
    }

    [[nodiscard]] std::span<const GlobalEntityId> global_ids(EntityKind kind) const {
        return kind_partition(kind).global_ids;
    }
    [[nodiscard]] std::span<const PartitionRank> owner_ranks(EntityKind kind) const {
        return kind_partition(kind).owner_ranks;
    }

    [[nodiscard]] GlobalEntityId global_id(EntityKind kind, LocalIndex local) const {
        return kind_partition(kind).global_ids.at(local_position(kind, local));
    }
    [[nodiscard]] PartitionRank owner_rank(EntityKind kind, LocalIndex local) const {
        return kind_partition(kind).owner_ranks.at(local_position(kind, local));
    }
    [[nodiscard]] EntityOwnership ownership(EntityKind kind, LocalIndex local) const {
        return owner_rank(kind, local) == local_rank_
                   ? EntityOwnership::owned
                   : EntityOwnership::ghost;
    }
    [[nodiscard]] bool is_owned(EntityKind kind, LocalIndex local) const {
        return ownership(kind, local) == EntityOwnership::owned;
    }
    [[nodiscard]] bool is_ghost(EntityKind kind, LocalIndex local) const {
        return ownership(kind, local) == EntityOwnership::ghost;
    }

    [[nodiscard]] bool contains_global(EntityKind kind, GlobalEntityId global) const {
        const auto& partition = kind_partition(kind);
        const auto found = lower_bound_global(partition, global);
        return found != partition.global_order.end() &&
               partition.global_ids[static_cast<std::size_t>(found->value())] == global;
    }

    [[nodiscard]] LocalIndex local_index(EntityKind kind, GlobalEntityId global) const {
        const auto& partition = kind_partition(kind);
        const auto found = lower_bound_global(partition, global);
        if (found == partition.global_order.end() ||
            partition.global_ids[static_cast<std::size_t>(found->value())] != global) {
            throw std::out_of_range(
                "mpmc::mesh::PartitionSnapshot: global entity ID is not local");
        }
        return *found;
    }

private:
    struct KindPartition {
        std::vector<GlobalEntityId> global_ids;
        std::vector<PartitionRank> owner_ranks;
        std::vector<LocalIndex> global_order;
        std::size_t owned_count{0U};
    };

    PartitionSnapshot(PartitionRank local_rank,
                      std::uint32_t rank_count,
                      std::array<KindPartition, 4> kinds)
        : local_rank_(local_rank),
          rank_count_(rank_count),
          kinds_(std::move(kinds)) {}

    [[nodiscard]] static std::size_t kind_index(EntityKind kind) {
        switch (kind) {
        case EntityKind::vertex: return 0U;
        case EntityKind::edge: return 1U;
        case EntityKind::face: return 2U;
        case EntityKind::cell: return 3U;
        }
        throw std::invalid_argument(
            "mpmc::mesh::PartitionSnapshot: invalid entity kind");
    }

    [[nodiscard]] const KindPartition& kind_partition(EntityKind kind) const {
        return kinds_[kind_index(kind)];
    }

    [[nodiscard]] std::size_t local_position(EntityKind kind, LocalIndex local) const {
        const auto& partition = kind_partition(kind);
        const std::size_t index = static_cast<std::size_t>(local.value());
        if (index >= partition.global_ids.size()) {
            throw std::out_of_range(
                "mpmc::mesh::PartitionSnapshot: local entity index out of range");
        }
        return index;
    }

    static void initialize_kind(const Topology& topology,
                                EntityKind kind,
                                PartitionRank local_rank,
                                std::uint32_t rank_count,
                                std::vector<PartitionRank> owner_ranks,
                                KindPartition& output) {
        const auto topology_ids = topology.global_ids(kind);
        if (owner_ranks.size() != topology_ids.size()) {
            throw std::invalid_argument(
                "mpmc::mesh::PartitionSnapshot: owner-rank count does not match topology");
        }

        output.global_ids.assign(topology_ids.begin(), topology_ids.end());
        output.owner_ranks = std::move(owner_ranks);
        output.global_order.reserve(output.global_ids.size());

        for (std::size_t local = 0; local < output.global_ids.size(); ++local) {
            const PartitionRank owner = output.owner_ranks[local];
            if (owner.value() >= rank_count) {
                throw std::out_of_range(
                    "mpmc::mesh::PartitionSnapshot: owner rank is outside rank_count");
            }
            if (owner == local_rank) {
                ++output.owned_count;
            }
            output.global_order.emplace_back(
                static_cast<LocalIndex::value_type>(local));
        }

        std::sort(
            output.global_order.begin(), output.global_order.end(),
            [&](LocalIndex left, LocalIndex right) {
                const auto left_id =
                    output.global_ids[static_cast<std::size_t>(left.value())];
                const auto right_id =
                    output.global_ids[static_cast<std::size_t>(right.value())];
                return left_id < right_id;
            });
    }

    [[nodiscard]] static std::vector<LocalIndex>::const_iterator lower_bound_global(
        const KindPartition& partition,
        GlobalEntityId global) {
        return std::lower_bound(
            partition.global_order.begin(),
            partition.global_order.end(),
            global,
            [&](LocalIndex local, GlobalEntityId target) {
                return partition.global_ids[
                           static_cast<std::size_t>(local.value())] < target;
            });
    }

    PartitionRank local_rank_;
    std::uint32_t rank_count_;
    std::array<KindPartition, 4> kinds_;
};

/// Strict single-rank baseline: rank 0 owns every local entity and there are no ghosts.
[[nodiscard]] inline PartitionSnapshot make_serial_partition_snapshot(
    const Topology& topology) {
    EntityOwnerRanks owners;
    owners.vertices.assign(
        topology.entity_count(EntityKind::vertex), PartitionRank{0U});
    owners.edges.assign(
        topology.entity_count(EntityKind::edge), PartitionRank{0U});
    owners.faces.assign(
        topology.entity_count(EntityKind::face), PartitionRank{0U});
    owners.cells.assign(
        topology.entity_count(EntityKind::cell), PartitionRank{0U});
    return PartitionSnapshot::create(
        topology, PartitionRank{0U}, 1U, std::move(owners));
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_PARTITION_SNAPSHOT_HPP
