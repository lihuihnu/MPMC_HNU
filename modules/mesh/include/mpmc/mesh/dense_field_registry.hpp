#ifndef MPMC_MESH_DENSE_FIELD_REGISTRY_HPP
#define MPMC_MESH_DENSE_FIELD_REGISTRY_HPP

#include <mpmc/mesh/dense_field.hpp>
#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/topology.hpp>

#include <array>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::mesh {

/// Immutable owning collection of dense fields keyed by (location, field id).
///
/// Registry lookup is intended for setup/control paths. Numerical hot loops
/// should retain a DenseFieldSnapshot reference and use LocalIndex/component
/// access directly. The same field id may exist at different entity locations;
/// duplicate ids at the same location are rejected.
///
/// DenseFieldSnapshot does not retain a Topology object, so create() rechecks
/// location entity counts against the supplied topology but cannot prove that a
/// same-sized field originated from that exact topology instance.
class DenseFieldRegistry {
public:
    DenseFieldRegistry(
        const DenseFieldRegistry&) = default;
    DenseFieldRegistry(
        DenseFieldRegistry&&) noexcept = default;
    DenseFieldRegistry& operator=(
        const DenseFieldRegistry&) = delete;
    DenseFieldRegistry& operator=(
        DenseFieldRegistry&&) = delete;
    ~DenseFieldRegistry() = default;

    [[nodiscard]] static DenseFieldRegistry create(
        const Topology& topology,
        std::vector<DenseFieldSnapshot> fields) {
        std::array<std::size_t, 4>
            counts{0U, 0U, 0U, 0U};

        for (std::size_t index = 0U;
             index < fields.size();
             ++index) {
            const auto location =
                fields[index].location();
            const std::size_t slot =
                location_index(location);
            if (fields[index].entity_count() !=
                topology.entity_count(location)) {
                throw std::invalid_argument(
                    "mpmc::mesh::DenseFieldRegistry: field entity count does not match registry topology");
            }

            for (std::size_t previous = 0U;
                 previous < index;
                 ++previous) {
                if (fields[previous].location() ==
                        location &&
                    fields[previous].metadata().id ==
                        fields[index].metadata().id) {
                    throw std::invalid_argument(
                        "mpmc::mesh::DenseFieldRegistry: duplicate field id at one entity location");
                }
            }
            ++counts[slot];
        }

        return DenseFieldRegistry{
            std::move(fields),
            counts};
    }

    [[nodiscard]] std::size_t size()
        const noexcept {
        return fields_.size();
    }

    [[nodiscard]] bool empty()
        const noexcept {
        return fields_.empty();
    }

    [[nodiscard]] std::span<
        const DenseFieldSnapshot>
    fields() const noexcept {
        return fields_;
    }

    [[nodiscard]] std::size_t field_count(
        EntityKind location) const {
        return counts_[
            location_index(location)];
    }

    [[nodiscard]] const DenseFieldSnapshot*
    find(
        EntityKind location,
        std::string_view id) const {
        (void)location_index(location);
        for (const auto& field : fields_) {
            if (field.location() == location &&
                field.metadata().id == id) {
                return &field;
            }
        }
        return nullptr;
    }

    [[nodiscard]] bool contains(
        EntityKind location,
        std::string_view id) const {
        return find(location, id) != nullptr;
    }

    [[nodiscard]] const DenseFieldSnapshot&
    at(
        EntityKind location,
        std::string_view id) const {
        const auto* field =
            find(location, id);
        if (field == nullptr) {
            throw std::out_of_range(
                "mpmc::mesh::DenseFieldRegistry: field key not found");
        }
        return *field;
    }

private:
    DenseFieldRegistry(
        std::vector<DenseFieldSnapshot> fields,
        std::array<std::size_t, 4> counts)
        : fields_(std::move(fields)),
          counts_(counts) {}

    [[nodiscard]] static std::size_t
    location_index(EntityKind location) {
        switch (location) {
        case EntityKind::vertex:
            return 0U;
        case EntityKind::edge:
            return 1U;
        case EntityKind::face:
            return 2U;
        case EntityKind::cell:
            return 3U;
        }
        throw std::invalid_argument(
            "mpmc::mesh::DenseFieldRegistry: invalid field location");
    }

    std::vector<DenseFieldSnapshot> fields_;
    std::array<std::size_t, 4> counts_;
};

} // namespace mpmc::mesh

#endif // MPMC_MESH_DENSE_FIELD_REGISTRY_HPP
