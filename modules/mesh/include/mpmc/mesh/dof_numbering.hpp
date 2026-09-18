#ifndef MPMC_MESH_DOF_NUMBERING_HPP
#define MPMC_MESH_DOF_NUMBERING_HPP

#include <mpmc/mesh/dof_layout.hpp>
#include <mpmc/mesh/partition_snapshot.hpp>

#include <algorithm>
#include <array>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

class GlobalDofIndex {
public:
    using value_type = std::uint64_t;

    explicit constexpr GlobalDofIndex(value_type value) noexcept : value_(value) {}

    [[nodiscard]] constexpr value_type value() const noexcept { return value_; }

    [[nodiscard]] friend constexpr auto operator<=>(const GlobalDofIndex&,
                                                     const GlobalDofIndex&) noexcept = default;

private:
    value_type value_;
};

class GlobalEntityOrdinal {
public:
    using value_type = std::uint64_t;

    explicit constexpr GlobalEntityOrdinal(value_type value) noexcept : value_(value) {}

    [[nodiscard]] constexpr value_type value() const noexcept { return value_; }

    [[nodiscard]] friend constexpr auto operator<=>(const GlobalEntityOrdinal&,
                                                     const GlobalEntityOrdinal&) noexcept = default;

private:
    value_type value_;
};

struct GlobalEntityOrdinalRecord {
    GlobalEntityId global_id;
    GlobalEntityOrdinal ordinal;
};

struct GlobalEntityNumberingInput {
    std::uint64_t global_vertex_count{0U};
    std::uint64_t global_face_count{0U};
    std::uint64_t global_cell_count{0U};

    std::vector<GlobalEntityOrdinalRecord> vertices;
    std::vector<GlobalEntityOrdinalRecord> faces;
    std::vector<GlobalEntityOrdinalRecord> cells;
};

/// Immutable partition-aware mapping between local and contiguous global scalar DoFs.
///
/// Global ordering mirrors DofLayout:
///   [cell block][face block][vertex block]
///
/// Within a global location block:
///   global entity ordinal -> variable declaration order -> component.
///
/// A serial snapshot can derive ordinals directly from local ordering. A generic
/// local snapshot cannot infer contiguous cross-rank ordinals from owner ranks alone,
/// so create_local() requires explicit global entity ordinals supplied by a future
/// distributor/MPI layer and validates them against stable GlobalEntityId values.
class DofNumberingSnapshot {
private:
    struct LocationNumbering {
        EntityKind kind{EntityKind::cell};
        std::size_t local_entity_count{0U};
        std::size_t dofs_per_entity{0U};
        std::size_t local_scalar_base{0U};
        std::size_t local_scalar_count{0U};

        std::uint64_t global_entity_count{0U};
        std::uint64_t global_scalar_base{0U};
        std::uint64_t global_scalar_count{0U};

        std::vector<GlobalEntityOrdinal> local_ordinals;
        std::vector<PartitionRank> owner_ranks;
        std::vector<LocalIndex> ordinal_order;
    };

public:
    DofNumberingSnapshot(const DofNumberingSnapshot&) = default;
    DofNumberingSnapshot(DofNumberingSnapshot&&) noexcept = default;
    DofNumberingSnapshot& operator=(const DofNumberingSnapshot&) = delete;
    DofNumberingSnapshot& operator=(DofNumberingSnapshot&&) = delete;
    ~DofNumberingSnapshot() = default;

    [[nodiscard]] static DofNumberingSnapshot create_serial(
        const DofLayout& layout,
        const PartitionSnapshot& partition) {
        if (!partition.is_serial() ||
            partition.local_rank().value() != 0U ||
            partition.rank_count() != 1U) {
            throw std::invalid_argument(
                "mpmc::mesh::DofNumberingSnapshot: serial numbering requires strict rank-0 serial partition");
        }

        GlobalEntityNumberingInput numbering;
        numbering.global_cell_count =
            static_cast<std::uint64_t>(partition.entity_count(EntityKind::cell));
        numbering.global_face_count =
            static_cast<std::uint64_t>(partition.entity_count(EntityKind::face));
        numbering.global_vertex_count =
            static_cast<std::uint64_t>(partition.entity_count(EntityKind::vertex));

        append_serial_records(partition, EntityKind::cell, numbering.cells);
        append_serial_records(partition, EntityKind::face, numbering.faces);
        append_serial_records(partition, EntityKind::vertex, numbering.vertices);

        return create_local(layout, partition, std::move(numbering));
    }

    [[nodiscard]] static DofNumberingSnapshot create_local(
        const DofLayout& layout,
        const PartitionSnapshot& partition,
        GlobalEntityNumberingInput numbering) {
        std::array<LocationNumbering, 3> locations{};
        std::uint64_t global_scalar_base = 0U;
        std::size_t owned_dof_count = 0U;
        std::size_t ghost_dof_count = 0U;

        initialize_location(
            layout, partition, EntityKind::cell,
            numbering.global_cell_count, std::move(numbering.cells),
            global_scalar_base, owned_dof_count, ghost_dof_count, locations[0]);
        initialize_location(
            layout, partition, EntityKind::face,
            numbering.global_face_count, std::move(numbering.faces),
            global_scalar_base, owned_dof_count, ghost_dof_count, locations[1]);
        initialize_location(
            layout, partition, EntityKind::vertex,
            numbering.global_vertex_count, std::move(numbering.vertices),
            global_scalar_base, owned_dof_count, ghost_dof_count, locations[2]);

        if (owned_dof_count >
            std::numeric_limits<std::size_t>::max() - ghost_dof_count) {
            throw std::length_error(
                "mpmc::mesh::DofNumberingSnapshot: local DoF ownership count overflow");
        }
        if (owned_dof_count + ghost_dof_count != layout.total_dof_count()) {
            throw std::invalid_argument(
                "mpmc::mesh::DofNumberingSnapshot: layout and partition local DoF counts disagree");
        }

        return DofNumberingSnapshot{
            partition.local_rank(),
            partition.rank_count(),
            std::move(locations),
            layout.total_dof_count(),
            global_scalar_base,
            owned_dof_count,
            ghost_dof_count};
    }

    [[nodiscard]] PartitionRank local_rank() const noexcept { return local_rank_; }
    [[nodiscard]] std::uint32_t rank_count() const noexcept { return rank_count_; }
    [[nodiscard]] bool is_serial() const noexcept { return rank_count_ == 1U; }

    [[nodiscard]] std::size_t local_dof_count() const noexcept {
        return local_dof_count_;
    }
    [[nodiscard]] std::uint64_t global_dof_count() const noexcept {
        return global_dof_count_;
    }
    [[nodiscard]] std::size_t owned_dof_count() const noexcept {
        return owned_dof_count_;
    }
    [[nodiscard]] std::size_t ghost_dof_count() const noexcept {
        return ghost_dof_count_;
    }

    [[nodiscard]] EntityOwnership ownership(std::size_t local_scalar) const {
        const auto decoded = decode_local(local_scalar);
        const auto& location = locations_[decoded.location_slot];
        return location.owner_ranks[decoded.local_entity] == local_rank_
                   ? EntityOwnership::owned
                   : EntityOwnership::ghost;
    }
    [[nodiscard]] bool is_owned(std::size_t local_scalar) const {
        return ownership(local_scalar) == EntityOwnership::owned;
    }
    [[nodiscard]] bool is_ghost(std::size_t local_scalar) const {
        return ownership(local_scalar) == EntityOwnership::ghost;
    }

    [[nodiscard]] GlobalDofIndex global_index(std::size_t local_scalar) const {
        const auto decoded = decode_local(local_scalar);
        const auto& location = locations_[decoded.location_slot];
        const std::uint64_t ordinal =
            location.local_ordinals[decoded.local_entity].value();
        const std::uint64_t dofs_per_entity =
            static_cast<std::uint64_t>(location.dofs_per_entity);
        const std::uint64_t intra =
            static_cast<std::uint64_t>(decoded.intra_entity_dof);
        return GlobalDofIndex{
            location.global_scalar_base + ordinal * dofs_per_entity + intra};
    }

    [[nodiscard]] bool contains_global(GlobalDofIndex global) const {
        return find_local_scalar(global).found;
    }

    [[nodiscard]] std::size_t local_scalar(GlobalDofIndex global) const {
        const auto found = find_local_scalar(global);
        if (!found.found) {
            throw std::out_of_range(
                "mpmc::mesh::DofNumberingSnapshot: global DoF is not local");
        }
        return found.local_scalar;
    }

    [[nodiscard]] std::uint64_t global_entity_count(EntityKind location) const {
        return location_numbering(location).global_entity_count;
    }
    [[nodiscard]] std::uint64_t global_location_offset(EntityKind location) const {
        return location_numbering(location).global_scalar_base;
    }
    [[nodiscard]] std::uint64_t global_location_dof_count(EntityKind location) const {
        return location_numbering(location).global_scalar_count;
    }

private:
    struct DecodedLocal {
        std::size_t location_slot;
        std::size_t local_entity;
        std::size_t intra_entity_dof;
    };

    struct LocalSearchResult {
        bool found;
        std::size_t local_scalar;
    };

    DofNumberingSnapshot(
        PartitionRank local_rank,
        std::uint32_t rank_count,
        std::array<LocationNumbering, 3> locations,
        std::size_t local_dof_count,
        std::uint64_t global_dof_count,
        std::size_t owned_dof_count,
        std::size_t ghost_dof_count)
        : local_rank_(local_rank),
          rank_count_(rank_count),
          locations_(std::move(locations)),
          local_dof_count_(local_dof_count),
          global_dof_count_(global_dof_count),
          owned_dof_count_(owned_dof_count),
          ghost_dof_count_(ghost_dof_count) {}

    [[nodiscard]] static std::size_t location_slot(EntityKind location) {
        switch (location) {
        case EntityKind::cell: return 0U;
        case EntityKind::face: return 1U;
        case EntityKind::vertex: return 2U;
        case EntityKind::edge:
            throw std::invalid_argument(
                "mpmc::mesh::DofNumberingSnapshot: edge DoFs are not supported by DofLayout");
        }
        throw std::invalid_argument(
            "mpmc::mesh::DofNumberingSnapshot: invalid DoF location");
    }

    [[nodiscard]] const LocationNumbering& location_numbering(
        EntityKind location) const {
        return locations_[location_slot(location)];
    }

    static void append_serial_records(
        const PartitionSnapshot& partition,
        EntityKind kind,
        std::vector<GlobalEntityOrdinalRecord>& output) {
        const auto ids = partition.global_ids(kind);
        output.reserve(ids.size());
        for (std::size_t local = 0; local < ids.size(); ++local) {
            output.push_back(GlobalEntityOrdinalRecord{
                ids[local],
                GlobalEntityOrdinal{static_cast<std::uint64_t>(local)}});
        }
    }

    static void initialize_location(
        const DofLayout& layout,
        const PartitionSnapshot& partition,
        EntityKind kind,
        std::uint64_t global_entity_count,
        std::vector<GlobalEntityOrdinalRecord> records,
        std::uint64_t& global_scalar_base,
        std::size_t& owned_dof_count,
        std::size_t& ghost_dof_count,
        LocationNumbering& output) {
        const std::size_t local_entity_count = layout.entity_count(kind);
        if (partition.entity_count(kind) != local_entity_count) {
            throw std::invalid_argument(
                "mpmc::mesh::DofNumberingSnapshot: DofLayout and PartitionSnapshot entity counts disagree");
        }
        if (records.size() != local_entity_count) {
            throw std::invalid_argument(
                "mpmc::mesh::DofNumberingSnapshot: global ordinal record count does not match local entities");
        }

        output.kind = kind;
        output.local_entity_count = local_entity_count;
        output.dofs_per_entity = layout.dofs_per_entity(kind);
        output.local_scalar_base = layout.location_offset(kind);
        output.local_scalar_count = layout.location_scalar_count(kind);
        output.global_entity_count = global_entity_count;
        output.global_scalar_base = global_scalar_base;
        output.local_ordinals.reserve(local_entity_count);
        output.owner_ranks.assign(
            partition.owner_ranks(kind).begin(), partition.owner_ranks(kind).end());
        output.ordinal_order.reserve(local_entity_count);

        const std::uint64_t dofs_per_entity =
            checked_u64(output.dofs_per_entity,
                        "mpmc::mesh::DofNumberingSnapshot: DoFs-per-entity exceed global index capacity");
        if (global_entity_count != 0U &&
            dofs_per_entity >
                std::numeric_limits<std::uint64_t>::max() / global_entity_count) {
            throw std::length_error(
                "mpmc::mesh::DofNumberingSnapshot: global location DoF count overflow");
        }
        output.global_scalar_count = global_entity_count * dofs_per_entity;
        if (output.global_scalar_count >
            std::numeric_limits<std::uint64_t>::max() - global_scalar_base) {
            throw std::length_error(
                "mpmc::mesh::DofNumberingSnapshot: global DoF count overflow");
        }
        global_scalar_base += output.global_scalar_count;

        for (std::size_t local = 0; local < local_entity_count; ++local) {
            const auto expected_id =
                partition.global_id(
                    kind, LocalIndex{static_cast<LocalIndex::value_type>(local)});
            const auto& record = records[local];
            if (record.global_id != expected_id) {
                throw std::invalid_argument(
                    "mpmc::mesh::DofNumberingSnapshot: ordinal record GlobalEntityId does not match partition local ordering");
            }
            if (record.ordinal.value() >= global_entity_count) {
                throw std::out_of_range(
                    "mpmc::mesh::DofNumberingSnapshot: global entity ordinal is outside global entity count");
            }
            output.local_ordinals.push_back(record.ordinal);
            output.ordinal_order.emplace_back(
                static_cast<LocalIndex::value_type>(local));
        }

        std::sort(
            output.ordinal_order.begin(), output.ordinal_order.end(),
            [&](LocalIndex left, LocalIndex right) {
                return output.local_ordinals[
                           static_cast<std::size_t>(left.value())] <
                       output.local_ordinals[
                           static_cast<std::size_t>(right.value())];
            });
        for (std::size_t i = 1U; i < output.ordinal_order.size(); ++i) {
            const auto previous =
                output.local_ordinals[static_cast<std::size_t>(
                    output.ordinal_order[i - 1U].value())];
            const auto current =
                output.local_ordinals[static_cast<std::size_t>(
                    output.ordinal_order[i].value())];
            if (previous == current) {
                throw std::invalid_argument(
                    "mpmc::mesh::DofNumberingSnapshot: duplicate local global entity ordinal");
            }
        }

        if (output.dofs_per_entity != 0U) {
            const std::size_t owned_entities = partition.owned_count(kind);
            const std::size_t ghost_entities = partition.ghost_count(kind);
            if (owned_entities >
                std::numeric_limits<std::size_t>::max() / output.dofs_per_entity) {
                throw std::length_error(
                    "mpmc::mesh::DofNumberingSnapshot: owned local DoF count overflow");
            }
            if (ghost_entities >
                std::numeric_limits<std::size_t>::max() / output.dofs_per_entity) {
                throw std::length_error(
                    "mpmc::mesh::DofNumberingSnapshot: ghost local DoF count overflow");
            }
            const std::size_t owned_here =
                owned_entities * output.dofs_per_entity;
            const std::size_t ghost_here =
                ghost_entities * output.dofs_per_entity;
            if (owned_here >
                std::numeric_limits<std::size_t>::max() - owned_dof_count ||
                ghost_here >
                std::numeric_limits<std::size_t>::max() - ghost_dof_count) {
                throw std::length_error(
                    "mpmc::mesh::DofNumberingSnapshot: accumulated local DoF count overflow");
            }
            owned_dof_count += owned_here;
            ghost_dof_count += ghost_here;
        }
    }

    [[nodiscard]] static std::uint64_t checked_u64(
        std::size_t value,
        const char* message) {
        if constexpr (sizeof(std::size_t) > sizeof(std::uint64_t)) {
            if (value > static_cast<std::size_t>(
                            std::numeric_limits<std::uint64_t>::max())) {
                throw std::length_error(message);
            }
        }
        return static_cast<std::uint64_t>(value);
    }

    [[nodiscard]] DecodedLocal decode_local(std::size_t local_scalar) const {
        if (local_scalar >= local_dof_count_) {
            throw std::out_of_range(
                "mpmc::mesh::DofNumberingSnapshot: local scalar index out of range");
        }

        for (std::size_t slot = 0U; slot < locations_.size(); ++slot) {
            const auto& location = locations_[slot];
            if (local_scalar >= location.local_scalar_base &&
                local_scalar - location.local_scalar_base <
                    location.local_scalar_count) {
                const std::size_t relative =
                    local_scalar - location.local_scalar_base;
                if (location.dofs_per_entity == 0U) {
                    throw std::logic_error(
                        "mpmc::mesh::DofNumberingSnapshot: nonempty zero-width local DoF block");
                }
                return DecodedLocal{
                    slot,
                    relative / location.dofs_per_entity,
                    relative % location.dofs_per_entity};
            }
        }

        throw std::logic_error(
            "mpmc::mesh::DofNumberingSnapshot: local scalar does not belong to a layout block");
    }

    [[nodiscard]] LocalSearchResult find_local_scalar(GlobalDofIndex global) const {
        if (global.value() >= global_dof_count_) {
            return LocalSearchResult{false, 0U};
        }

        for (const auto& location : locations_) {
            if (global.value() >= location.global_scalar_base &&
                global.value() - location.global_scalar_base <
                    location.global_scalar_count) {
                if (location.dofs_per_entity == 0U) {
                    return LocalSearchResult{false, 0U};
                }
                const std::uint64_t relative =
                    global.value() - location.global_scalar_base;
                const std::uint64_t width =
                    static_cast<std::uint64_t>(location.dofs_per_entity);
                const GlobalEntityOrdinal target{relative / width};
                const std::size_t intra =
                    static_cast<std::size_t>(relative % width);

                const auto found = std::lower_bound(
                    location.ordinal_order.begin(),
                    location.ordinal_order.end(),
                    target,
                    [&](LocalIndex local, GlobalEntityOrdinal ordinal) {
                        return location.local_ordinals[
                                   static_cast<std::size_t>(local.value())] <
                               ordinal;
                    });
                if (found == location.ordinal_order.end()) {
                    return LocalSearchResult{false, 0U};
                }
                const std::size_t local_entity =
                    static_cast<std::size_t>(found->value());
                if (location.local_ordinals[local_entity] != target) {
                    return LocalSearchResult{false, 0U};
                }
                return LocalSearchResult{
                    true,
                    location.local_scalar_base +
                        local_entity * location.dofs_per_entity + intra};
            }
        }
        return LocalSearchResult{false, 0U};
    }

    PartitionRank local_rank_;
    std::uint32_t rank_count_;
    std::array<LocationNumbering, 3> locations_;
    std::size_t local_dof_count_;
    std::uint64_t global_dof_count_;
    std::size_t owned_dof_count_;
    std::size_t ghost_dof_count_;
};

} // namespace mpmc::mesh

#endif // MPMC_MESH_DOF_NUMBERING_HPP
