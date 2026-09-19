#ifndef MPMC_MESH_SHARED_ENTITY_PLAN_HPP
#define MPMC_MESH_SHARED_ENTITY_PLAN_HPP

#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/partition_snapshot.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

/// Canonical description of one owned entity mirrored as one ghost on another rank.
///
/// owner_local and ghost_local are rank-local indices in their respective rank
/// topologies. GlobalEntityId is the stable cross-rank identity used to reconcile
/// the two endpoints.
struct SharedEntityLink {
    EntityKind kind;
    GlobalEntityId global_id;
    PartitionRank owner_rank;
    LocalIndex owner_local;
    PartitionRank ghost_rank;
    LocalIndex ghost_local;
};

/// One local endpoint in a compact halo send/receive list.
///
/// local is this rank's entity index. remote_local is the corresponding entity
/// index on the neighbor rank. For receive entries, remote_local is owner-local.
struct HaloEntityRef {
    EntityKind kind;
    GlobalEntityId global_id;
    LocalIndex local;
    LocalIndex remote_local;
};

/// Offsets into flat send/receive arrays for one neighbor rank.
struct NeighborExchangeRange {
    PartitionRank rank;
    std::size_t send_begin;
    std::size_t send_count;
    std::size_t receive_begin;
    std::size_t receive_count;
};

/// Immutable non-MPI halo exchange contract derived from canonical shared links.
///
/// Every local ghost entity must appear in exactly one SharedEntityLink whose
/// owner_rank matches PartitionSnapshot::owner_rank(). Locally owned entities may
/// appear in zero or more links because multiple neighbor ranks may mirror them.
///
/// Neighbor ranges are sorted by rank. Send and receive arrays are flat and grouped
/// by those ranges; entries inside a neighbor are deterministically sorted by
/// entity kind, GlobalEntityId and local index.
class SharedEntityPlan {
public:
    SharedEntityPlan(const SharedEntityPlan&) = default;
    SharedEntityPlan(SharedEntityPlan&&) noexcept = default;
    SharedEntityPlan& operator=(const SharedEntityPlan&) = delete;
    SharedEntityPlan& operator=(SharedEntityPlan&&) = delete;
    ~SharedEntityPlan() = default;

    [[nodiscard]] static SharedEntityPlan create(
        const PartitionSnapshot& partition,
        std::span<const SharedEntityLink> links) {
        std::array<std::vector<std::uint8_t>, 4> ghost_seen;
        for (const EntityKind kind :
             {EntityKind::vertex, EntityKind::edge,
              EntityKind::face, EntityKind::cell}) {
            ghost_seen[kind_index(kind)].assign(
                partition.entity_count(kind), std::uint8_t{0U});
        }

        std::vector<Pending> pending_sends;
        std::vector<Pending> pending_receives;
        pending_sends.reserve(links.size());
        pending_receives.reserve(links.size());

        for (const SharedEntityLink& link : links) {
            validate_rank(link.owner_rank, partition.rank_count());
            validate_rank(link.ghost_rank, partition.rank_count());
            if (link.owner_rank == link.ghost_rank) {
                throw std::invalid_argument(
                    "mpmc::mesh::SharedEntityPlan: owner and ghost ranks must differ");
            }

            const bool local_is_owner =
                link.owner_rank == partition.local_rank();
            const bool local_is_ghost =
                link.ghost_rank == partition.local_rank();
            if (!local_is_owner && !local_is_ghost) {
                throw std::invalid_argument(
                    "mpmc::mesh::SharedEntityPlan: shared link does not touch local rank");
            }

            if (local_is_owner) {
                validate_local_endpoint(
                    partition, link.kind, link.owner_local, link.global_id,
                    EntityOwnership::owned);
                pending_sends.push_back(Pending{
                    link.ghost_rank,
                    HaloEntityRef{
                        link.kind,
                        link.global_id,
                        link.owner_local,
                        link.ghost_local}});
            }

            if (local_is_ghost) {
                validate_local_endpoint(
                    partition, link.kind, link.ghost_local, link.global_id,
                    EntityOwnership::ghost);
                if (partition.owner_rank(link.kind, link.ghost_local) !=
                    link.owner_rank) {
                    throw std::invalid_argument(
                        "mpmc::mesh::SharedEntityPlan: ghost owner rank disagrees with partition");
                }

                auto& seen = ghost_seen[kind_index(link.kind)];
                const std::size_t ghost_position =
                    static_cast<std::size_t>(link.ghost_local.value());
                if (seen[ghost_position] != 0U) {
                    throw std::invalid_argument(
                        "mpmc::mesh::SharedEntityPlan: local ghost appears in more than one shared link");
                }
                seen[ghost_position] = std::uint8_t{1U};

                pending_receives.push_back(Pending{
                    link.owner_rank,
                    HaloEntityRef{
                        link.kind,
                        link.global_id,
                        link.ghost_local,
                        link.owner_local}});
            }
        }

        for (const EntityKind kind :
             {EntityKind::vertex, EntityKind::edge,
              EntityKind::face, EntityKind::cell}) {
            const auto& seen = ghost_seen[kind_index(kind)];
            for (std::size_t local = 0; local < seen.size(); ++local) {
                const LocalIndex index{
                    static_cast<LocalIndex::value_type>(local)};
                if (partition.is_ghost(kind, index) && seen[local] == 0U) {
                    throw std::invalid_argument(
                        "mpmc::mesh::SharedEntityPlan: local ghost is missing an owner link");
                }
            }
        }

        sort_pending(pending_sends);
        sort_pending(pending_receives);

        std::vector<NeighborExchangeRange> neighbors;
        std::vector<HaloEntityRef> sends;
        std::vector<HaloEntityRef> receives;
        sends.reserve(pending_sends.size());
        receives.reserve(pending_receives.size());

        std::size_t send_position = 0U;
        std::size_t receive_position = 0U;
        while (send_position < pending_sends.size() ||
               receive_position < pending_receives.size()) {
            const PartitionRank neighbor = next_neighbor(
                pending_sends, send_position,
                pending_receives, receive_position);

            const std::size_t send_begin = sends.size();
            while (send_position < pending_sends.size() &&
                   pending_sends[send_position].neighbor == neighbor) {
                sends.push_back(pending_sends[send_position].entity);
                ++send_position;
            }

            const std::size_t receive_begin = receives.size();
            while (receive_position < pending_receives.size() &&
                   pending_receives[receive_position].neighbor == neighbor) {
                receives.push_back(pending_receives[receive_position].entity);
                ++receive_position;
            }

            neighbors.push_back(NeighborExchangeRange{
                neighbor,
                send_begin,
                sends.size() - send_begin,
                receive_begin,
                receives.size() - receive_begin});
        }

        return SharedEntityPlan{
            partition.local_rank(),
            partition.rank_count(),
            std::move(neighbors),
            std::move(sends),
            std::move(receives)};
    }

    [[nodiscard]] PartitionRank local_rank() const noexcept { return local_rank_; }
    [[nodiscard]] std::uint32_t rank_count() const noexcept { return rank_count_; }
    [[nodiscard]] bool is_serial() const noexcept { return rank_count_ == 1U; }

    [[nodiscard]] std::size_t neighbor_count() const noexcept {
        return neighbors_.size();
    }
    [[nodiscard]] std::size_t send_count() const noexcept {
        return sends_.size();
    }
    [[nodiscard]] std::size_t receive_count() const noexcept {
        return receives_.size();
    }

    [[nodiscard]] std::span<const NeighborExchangeRange> neighbors() const noexcept {
        return neighbors_;
    }
    [[nodiscard]] std::span<const HaloEntityRef> send_entities() const noexcept {
        return sends_;
    }
    [[nodiscard]] std::span<const HaloEntityRef> receive_entities() const noexcept {
        return receives_;
    }

    [[nodiscard]] bool contains_neighbor(PartitionRank rank) const {
        const auto found = find_neighbor(rank);
        return found != neighbors_.end() && found->rank == rank;
    }

    [[nodiscard]] std::span<const HaloEntityRef> send_entities_to(
        PartitionRank rank) const {
        const auto found = require_neighbor(rank);
        return std::span<const HaloEntityRef>{sends_}.subspan(
            found->send_begin, found->send_count);
    }

    [[nodiscard]] std::span<const HaloEntityRef> receive_entities_from(
        PartitionRank rank) const {
        const auto found = require_neighbor(rank);
        return std::span<const HaloEntityRef>{receives_}.subspan(
            found->receive_begin, found->receive_count);
    }

private:
    struct Pending {
        PartitionRank neighbor;
        HaloEntityRef entity;
    };

    SharedEntityPlan(PartitionRank local_rank,
                     std::uint32_t rank_count,
                     std::vector<NeighborExchangeRange> neighbors,
                     std::vector<HaloEntityRef> sends,
                     std::vector<HaloEntityRef> receives)
        : local_rank_(local_rank),
          rank_count_(rank_count),
          neighbors_(std::move(neighbors)),
          sends_(std::move(sends)),
          receives_(std::move(receives)) {}

    [[nodiscard]] static std::size_t kind_index(EntityKind kind) {
        switch (kind) {
        case EntityKind::vertex: return 0U;
        case EntityKind::edge: return 1U;
        case EntityKind::face: return 2U;
        case EntityKind::cell: return 3U;
        }
        throw std::invalid_argument(
            "mpmc::mesh::SharedEntityPlan: invalid entity kind");
    }

    static void validate_rank(PartitionRank rank, std::uint32_t rank_count) {
        if (rank.value() >= rank_count) {
            throw std::out_of_range(
                "mpmc::mesh::SharedEntityPlan: rank is outside partition rank_count");
        }
    }

    static void validate_local_endpoint(
        const PartitionSnapshot& partition,
        EntityKind kind,
        LocalIndex local,
        GlobalEntityId global,
        EntityOwnership expected) {
        const GlobalEntityId actual = partition.global_id(kind, local);
        if (actual != global) {
            throw std::invalid_argument(
                "mpmc::mesh::SharedEntityPlan: local endpoint GlobalEntityId mismatch");
        }
        if (partition.ownership(kind, local) != expected) {
            throw std::invalid_argument(
                "mpmc::mesh::SharedEntityPlan: local endpoint ownership mismatch");
        }
    }

    [[nodiscard]] static std::uint8_t kind_order(EntityKind kind) {
        return static_cast<std::uint8_t>(kind_index(kind));
    }

    [[nodiscard]] static bool entity_less(
        const HaloEntityRef& left,
        const HaloEntityRef& right) {
        const std::uint8_t left_kind = kind_order(left.kind);
        const std::uint8_t right_kind = kind_order(right.kind);
        if (left_kind != right_kind) {
            return left_kind < right_kind;
        }
        if (left.global_id != right.global_id) {
            return left.global_id < right.global_id;
        }
        if (left.local != right.local) {
            return left.local < right.local;
        }
        return left.remote_local < right.remote_local;
    }

    static void sort_pending(std::vector<Pending>& pending) {
        std::sort(
            pending.begin(), pending.end(),
            [](const Pending& left, const Pending& right) {
                if (left.neighbor != right.neighbor) {
                    return left.neighbor < right.neighbor;
                }
                return entity_less(left.entity, right.entity);
            });

        for (std::size_t i = 1U; i < pending.size(); ++i) {
            const Pending& previous = pending[i - 1U];
            const Pending& current = pending[i];
            if (previous.neighbor == current.neighbor &&
                previous.entity.kind == current.entity.kind &&
                previous.entity.global_id == current.entity.global_id) {
                throw std::invalid_argument(
                    "mpmc::mesh::SharedEntityPlan: duplicate shared entity for one neighbor");
            }
        }
    }

    [[nodiscard]] static PartitionRank next_neighbor(
        const std::vector<Pending>& sends,
        std::size_t send_position,
        const std::vector<Pending>& receives,
        std::size_t receive_position) {
        if (send_position == sends.size()) {
            return receives[receive_position].neighbor;
        }
        if (receive_position == receives.size()) {
            return sends[send_position].neighbor;
        }
        return sends[send_position].neighbor < receives[receive_position].neighbor
                   ? sends[send_position].neighbor
                   : receives[receive_position].neighbor;
    }

    [[nodiscard]] std::vector<NeighborExchangeRange>::const_iterator find_neighbor(
        PartitionRank rank) const {
        return std::lower_bound(
            neighbors_.begin(), neighbors_.end(), rank,
            [](const NeighborExchangeRange& neighbor, PartitionRank target) {
                return neighbor.rank < target;
            });
    }

    [[nodiscard]] std::vector<NeighborExchangeRange>::const_iterator require_neighbor(
        PartitionRank rank) const {
        const auto found = find_neighbor(rank);
        if (found == neighbors_.end() || found->rank != rank) {
            throw std::out_of_range(
                "mpmc::mesh::SharedEntityPlan: rank is not a halo neighbor");
        }
        return found;
    }

    PartitionRank local_rank_;
    std::uint32_t rank_count_;
    std::vector<NeighborExchangeRange> neighbors_;
    std::vector<HaloEntityRef> sends_;
    std::vector<HaloEntityRef> receives_;
};

} // namespace mpmc::mesh

#endif // MPMC_MESH_SHARED_ENTITY_PLAN_HPP
