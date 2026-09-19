#ifndef MPMC_MESH_TOPOLOGY_HPP
#define MPMC_MESH_TOPOLOGY_HPP

#include <mpmc/mesh/csr_adjacency.hpp>
#include <mpmc/mesh/entity.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

/// Immutable owning snapshot of mesh entity identities and supplied topology relations.
///
/// Global IDs are unique within each EntityKind. Together, (kind, GlobalEntityId)
/// forms the stable entity identity. Relations occupy fixed kind-pair slots, so
/// relation lookup never requires a string or hash lookup.
class Topology {
public:
    static constexpr std::size_t entity_kind_count = 4U;
    static constexpr std::size_t relation_slot_count = entity_kind_count * entity_kind_count;

    using GlobalIdArray = std::vector<GlobalEntityId>;

    struct EntityIds {
        GlobalIdArray vertices;
        GlobalIdArray edges;
        GlobalIdArray faces;
        GlobalIdArray cells;
    };

    Topology(EntityIds entity_ids, std::vector<CsrAdjacency> relations)
        : entity_ids_(std::move(entity_ids)) {
        validate_global_ids();
        for (auto& relation : relations) {
            insert_relation(std::move(relation));
        }
    }

    Topology(const Topology&) = default;
    Topology(Topology&&) noexcept = default;
    Topology& operator=(const Topology&) = delete;
    Topology& operator=(Topology&&) = delete;
    ~Topology() = default;

    [[nodiscard]] std::size_t entity_count(EntityKind kind) const {
        return ids_for(kind).size();
    }

    [[nodiscard]] std::span<const GlobalEntityId> global_ids(EntityKind kind) const {
        return ids_for(kind);
    }

    [[nodiscard]] GlobalEntityId global_id(EntityKind kind, LocalIndex local) const {
        const auto ids = global_ids(kind);
        const std::size_t index = static_cast<std::size_t>(local.value());
        if (index >= ids.size()) {
            throw std::out_of_range("mpmc::mesh::Topology: local entity index out of range");
        }
        return ids[index];
    }

    [[nodiscard]] std::size_t relation_count() const noexcept { return relation_count_; }

    [[nodiscard]] bool has_relation(EntityKind source, EntityKind target) const {
        return relations_[relation_slot(source, target)].has_value();
    }

    [[nodiscard]] const CsrAdjacency& relation(EntityKind source, EntityKind target) const {
        const auto& stored = relations_[relation_slot(source, target)];
        if (!stored) {
            throw std::out_of_range("mpmc::mesh::Topology: requested relation is absent");
        }
        return *stored;
    }

private:
    [[nodiscard]] static std::size_t kind_index(EntityKind kind) {
        switch (kind) {
        case EntityKind::vertex: return 0U;
        case EntityKind::edge: return 1U;
        case EntityKind::face: return 2U;
        case EntityKind::cell: return 3U;
        }
        throw std::invalid_argument("mpmc::mesh::Topology: invalid entity kind");
    }

    [[nodiscard]] static std::size_t relation_slot(EntityKind source, EntityKind target) {
        return kind_index(source) * entity_kind_count + kind_index(target);
    }

    [[nodiscard]] const GlobalIdArray& ids_for(EntityKind kind) const {
        switch (kind) {
        case EntityKind::vertex: return entity_ids_.vertices;
        case EntityKind::edge: return entity_ids_.edges;
        case EntityKind::face: return entity_ids_.faces;
        case EntityKind::cell: return entity_ids_.cells;
        }
        throw std::invalid_argument("mpmc::mesh::Topology: invalid entity kind");
    }

    void validate_global_ids() const {
        constexpr auto local_max = std::numeric_limits<LocalIndex::value_type>::max();
        const std::uint64_t local_capacity = static_cast<std::uint64_t>(local_max) + 1ULL;
        const std::array<const GlobalIdArray*, entity_kind_count> arrays{
            &entity_ids_.vertices, &entity_ids_.edges, &entity_ids_.faces, &entity_ids_.cells};

        for (const GlobalIdArray* ids : arrays) {
            if (static_cast<std::uint64_t>(ids->size()) > local_capacity) {
                throw std::length_error(
                    "mpmc::mesh::Topology: entity count exceeds LocalIndex capacity");
            }

            // Preserve caller ordering: duplicate detection uses a temporary scalar copy.
            std::vector<GlobalEntityId::value_type> sorted;
            sorted.reserve(ids->size());
            for (const GlobalEntityId id : *ids) {
                sorted.push_back(id.value());
            }
            std::sort(sorted.begin(), sorted.end());
            if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) {
                throw std::invalid_argument(
                    "mpmc::mesh::Topology: duplicate global ID within entity kind");
            }
        }
    }

    void insert_relation(CsrAdjacency relation) {
        const std::size_t source_count = ids_for(relation.source_kind()).size();
        const std::size_t target_count = ids_for(relation.target_kind()).size();
        if (relation.source_count() != source_count) {
            throw std::invalid_argument(
                "mpmc::mesh::Topology: relation source count does not match entity IDs");
        }
        if (relation.target_count() != target_count) {
            throw std::invalid_argument(
                "mpmc::mesh::Topology: relation target count does not match entity IDs");
        }

        const std::size_t slot = relation_slot(relation.source_kind(), relation.target_kind());
        if (relations_[slot]) {
            throw std::invalid_argument(
                "mpmc::mesh::Topology: duplicate relation for entity-kind pair");
        }
        relations_[slot].emplace(std::move(relation));
        ++relation_count_;
    }

    EntityIds entity_ids_;
    std::array<std::optional<CsrAdjacency>, relation_slot_count> relations_{};
    std::size_t relation_count_{0U};
};

} // namespace mpmc::mesh

#endif // MPMC_MESH_TOPOLOGY_HPP
