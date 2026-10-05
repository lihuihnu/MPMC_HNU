#ifndef MPMC_MESH_MESH_EXCHANGE_GROUP_HPP
#define MPMC_MESH_MESH_EXCHANGE_GROUP_HPP

#include <mpmc/mesh/topology.hpp>
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::mesh {

struct MeshExchangeGroup {
    EntityKind location;
    std::uint32_t tag;
    std::string name;
    std::vector<GlobalEntityId> members;
};

namespace mesh_exchange_detail {

// Stable IDs may be sparse and unrelated to local order. The caller owns the
// ID span for this lookup's lifetime. Build only on first query; empty groups
// allocate no index. Keep missing-ID policy and diagnostics in the caller.
class GroupEntityLookup {
public:
    explicit GroupEntityLookup(std::span<const GlobalEntityId> ids)
        : ids_(ids) {}

    [[nodiscard]] std::optional<std::size_t> find(GlobalEntityId id) {
        if (sorted_.empty()) {
            sorted_.reserve(ids_.size());
            for (std::size_t local = 0U; local < ids_.size(); ++local) {
                sorted_.emplace_back(ids_[local].value(), local);
            }
            std::sort(sorted_.begin(), sorted_.end());
        }
        const auto found = std::lower_bound(
            sorted_.begin(), sorted_.end(), id.value(),
            [](const auto& entry, GlobalEntityId::value_type value) {
                return entry.first < value;
            });
        if (found == sorted_.end() || found->first != id.value()) {
            return std::nullopt;
        }
        return found->second;
    }

private:
    std::span<const GlobalEntityId> ids_;
    std::vector<std::pair<GlobalEntityId::value_type, std::size_t>> sorted_;
};

} // namespace mesh_exchange_detail

inline void validate_mesh_exchange_groups(
        const Topology& topology,
        std::span<const MeshExchangeGroup> groups) {
        // Reuse one index per entity kind across all groups in this validation.
        std::array<std::optional<mesh_exchange_detail::GroupEntityLookup>, 4> lookups;
        for (std::size_t index = 0U;
             index < groups.size();
             ++index) {
            const auto& group = groups[index];
            if (group.tag == 0U) {
                throw std::invalid_argument(
                    "mpmc::mesh::MeshExchangeDocument: group tag zero is reserved");
            }
            if (group.name.find('\0') !=
                std::string::npos) {
                throw std::invalid_argument(
                    "mpmc::mesh::MeshExchangeDocument: group name cannot contain NUL");
            }
            for (std::size_t previous = 0U;
                 previous < index;
                 ++previous) {
                if (groups[previous].location ==
                        group.location &&
                    groups[previous].tag ==
                        group.tag) {
                    throw std::invalid_argument(
                        "mpmc::mesh::MeshExchangeDocument: duplicate group key");
                }
            }

            const auto ids =
                topology.global_ids(
                    group.location);
            // global_ids above validates the kind before it is used as an index.
            auto& lookup = lookups[static_cast<std::size_t>(group.location)];
            if (!lookup.has_value()) {
                lookup.emplace(ids);
            }
            for (const auto member : group.members) {
                if (!lookup->find(member).has_value()) {
                    throw std::invalid_argument(
                        "mpmc::mesh::MeshExchangeDocument: group member is absent from topology");
                }
            }
        }
    }


} // namespace mpmc::mesh

#endif // MPMC_MESH_MESH_EXCHANGE_GROUP_HPP
