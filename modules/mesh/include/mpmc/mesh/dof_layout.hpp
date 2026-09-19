#ifndef MPMC_MESH_DOF_LAYOUT_HPP
#define MPMC_MESH_DOF_LAYOUT_HPP

#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <map>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::mesh {

struct DofVariable {
    std::string id;
    EntityKind location;
    std::size_t component_count;
};

/// Immutable local scalar layout for logical variables attached to mesh entities.
///
/// Scalar ordering is fixed as:
///   [cell block][face block][edge block][vertex block]
///
/// Within one location the order is:
///   entity-major -> variable declaration order -> component.
///
/// The layout contains no physics and no numeric solution values. String lookup is
/// a control-path convenience; hot assembly should resolve a variable index once
/// and use scalar_offset(variable_index, entity, component), which is O(1).
class DofLayout {
public:
    DofLayout(const DofLayout&) = default;
    DofLayout(DofLayout&&) noexcept = default;
    DofLayout& operator=(const DofLayout&) = delete;
    DofLayout& operator=(DofLayout&&) = delete;
    ~DofLayout() = default;

    [[nodiscard]] static DofLayout create(
        const Topology& topology,
        std::vector<DofVariable> variables) {
        if (variables.empty()) {
            throw std::invalid_argument(
                "mpmc::mesh::DofLayout: at least one variable is required");
        }

        std::map<std::string, std::size_t, std::less<>> indices;
        std::vector<std::size_t> component_offsets;
        component_offsets.reserve(variables.size());
        std::array<LocationLayout, 4> locations{};

        for (std::size_t variable_index = 0; variable_index < variables.size();
             ++variable_index) {
            const auto& variable = variables[variable_index];
            require_text(
                variable.id,
                "mpmc::mesh::DofLayout: variable id must be nonblank and contain no NUL");
            const std::size_t location_slot = location_index(variable.location);
            if (variable.component_count == 0U) {
                throw std::invalid_argument(
                    "mpmc::mesh::DofLayout: component_count must be positive");
            }
            if (!indices.emplace(variable.id, variable_index).second) {
                throw std::invalid_argument(
                    "mpmc::mesh::DofLayout: duplicate variable id");
            }

            auto& location = locations[location_slot];
            component_offsets.push_back(location.dofs_per_entity);
            if (variable.component_count >
                std::numeric_limits<std::size_t>::max() - location.dofs_per_entity) {
                throw std::length_error(
                    "mpmc::mesh::DofLayout: per-entity DoF count overflow");
            }
            location.dofs_per_entity += variable.component_count;
        }

        std::size_t total_dof_count = 0U;
        for (const EntityKind location_kind :
             {EntityKind::cell, EntityKind::face,
              EntityKind::edge, EntityKind::vertex}) {
            auto& location = locations[location_index(location_kind)];
            location.entity_count = topology.entity_count(location_kind);
            location.base_offset = total_dof_count;

            if (location.entity_count != 0U &&
                location.dofs_per_entity >
                    std::numeric_limits<std::size_t>::max() / location.entity_count) {
                throw std::length_error(
                    "mpmc::mesh::DofLayout: location scalar count overflow");
            }
            location.scalar_count =
                location.entity_count * location.dofs_per_entity;
            if (location.scalar_count >
                std::numeric_limits<std::size_t>::max() - total_dof_count) {
                throw std::length_error(
                    "mpmc::mesh::DofLayout: total scalar count overflow");
            }
            total_dof_count += location.scalar_count;
        }

        return DofLayout{
            std::move(variables),
            std::move(component_offsets),
            std::move(indices),
            locations,
            total_dof_count};
    }

    [[nodiscard]] std::size_t variable_count() const noexcept {
        return variables_.size();
    }
    [[nodiscard]] std::size_t total_dof_count() const noexcept {
        return total_dof_count_;
    }

    [[nodiscard]] std::span<const DofVariable> variables() const noexcept {
        return variables_;
    }
    [[nodiscard]] const DofVariable& variable(std::size_t variable_index) const {
        return variables_.at(variable_index);
    }

    [[nodiscard]] bool contains(std::string_view id) const {
        return variable_indices_.contains(id);
    }
    [[nodiscard]] std::size_t variable_index(std::string_view id) const {
        const auto found = variable_indices_.find(id);
        if (found == variable_indices_.end()) {
            throw std::out_of_range(
                "mpmc::mesh::DofLayout: unknown variable id");
        }
        return found->second;
    }

    [[nodiscard]] std::size_t entity_count(EntityKind location) const {
        return layout_for(location).entity_count;
    }
    [[nodiscard]] std::size_t dofs_per_entity(EntityKind location) const {
        return layout_for(location).dofs_per_entity;
    }
    [[nodiscard]] std::size_t location_offset(EntityKind location) const {
        return layout_for(location).base_offset;
    }
    [[nodiscard]] std::size_t location_scalar_count(EntityKind location) const {
        return layout_for(location).scalar_count;
    }

    [[nodiscard]] std::size_t entity_offset(
        EntityKind location,
        LocalIndex entity) const {
        const auto& layout = layout_for(location);
        const std::size_t entity_index =
            static_cast<std::size_t>(entity.value());
        if (entity_index >= layout.entity_count) {
            throw std::out_of_range(
                "mpmc::mesh::DofLayout: entity index out of range");
        }
        return layout.base_offset + entity_index * layout.dofs_per_entity;
    }

    [[nodiscard]] std::size_t scalar_offset(
        std::size_t variable_index,
        LocalIndex entity,
        std::size_t component) const {
        const auto& variable_spec = variable(variable_index);
        if (component >= variable_spec.component_count) {
            throw std::out_of_range(
                "mpmc::mesh::DofLayout: component index out of range");
        }

        return entity_offset(variable_spec.location, entity) +
               component_offsets_.at(variable_index) + component;
    }

    [[nodiscard]] std::size_t scalar_offset(
        std::string_view variable_id,
        LocalIndex entity,
        std::size_t component) const {
        return scalar_offset(variable_index(variable_id), entity, component);
    }

private:
    struct LocationLayout {
        std::size_t entity_count{0U};
        std::size_t dofs_per_entity{0U};
        std::size_t base_offset{0U};
        std::size_t scalar_count{0U};
    };

    DofLayout(std::vector<DofVariable> variables,
              std::vector<std::size_t> component_offsets,
              std::map<std::string, std::size_t, std::less<>> variable_indices,
              std::array<LocationLayout, 4> locations,
              std::size_t total_dof_count)
        : variables_(std::move(variables)),
          component_offsets_(std::move(component_offsets)),
          variable_indices_(std::move(variable_indices)),
          locations_(locations),
          total_dof_count_(total_dof_count) {}

    static void require_text(const std::string& text, const char* message) {
        const bool nonblank =
            std::any_of(text.begin(), text.end(),
                        [](unsigned char character) { return character > 32U; });
        if (!nonblank || text.find('\0') != std::string::npos) {
            throw std::invalid_argument(message);
        }
    }

    [[nodiscard]] static std::size_t location_index(EntityKind location) {
        switch (location) {
        case EntityKind::cell: return 0U;
        case EntityKind::face: return 1U;
        case EntityKind::edge: return 2U;
        case EntityKind::vertex: return 3U;
        }
        throw std::invalid_argument(
            "mpmc::mesh::DofLayout: invalid variable location");
    }

    [[nodiscard]] const LocationLayout& layout_for(EntityKind location) const {
        return locations_[location_index(location)];
    }

    std::vector<DofVariable> variables_;
    std::vector<std::size_t> component_offsets_;
    std::map<std::string, std::size_t, std::less<>> variable_indices_;
    std::array<LocationLayout, 4> locations_;
    std::size_t total_dof_count_;
};

} // namespace mpmc::mesh

#endif // MPMC_MESH_DOF_LAYOUT_HPP
