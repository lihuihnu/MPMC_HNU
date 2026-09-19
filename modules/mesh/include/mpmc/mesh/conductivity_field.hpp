#ifndef MPMC_MESH_CONDUCTIVITY_FIELD_HPP
#define MPMC_MESH_CONDUCTIVITY_FIELD_HPP

#include <mpmc/mesh/dense_field.hpp>
#include <mpmc/mesh/dense_field_registry.hpp>
#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>

namespace mpmc::mesh {

/// Contract for a scalar conductivity field.
///
/// Conductivity is intentionally not classified here as thermal, electrical or
/// another transport property because the repository has no such global
/// semantic yet. The caller must therefore provide the exact field id, entity
/// location and unit expected by the consuming model. This layer performs no
/// unit conversion.
struct ScalarConductivityFieldContract {
    EntityKind location;
    std::string id;
    std::string unit;
};

namespace conductivity_field_detail {

inline void require_text(
    const std::string& text,
    const char* message) {
    const bool nonblank =
        std::any_of(
            text.begin(),
            text.end(),
            [](unsigned char character) {
                return character > 32U;
            });
    if (!nonblank ||
        text.find('\0') !=
            std::string::npos) {
        throw std::invalid_argument(message);
    }
}

inline void validate_location(
    EntityKind location) {
    switch (location) {
    case EntityKind::vertex:
    case EntityKind::edge:
    case EntityKind::face:
    case EntityKind::cell:
        return;
    }
    throw std::invalid_argument(
        "mpmc::mesh::ScalarConductivityFieldContract: invalid field location");
}

} // namespace conductivity_field_detail

/// Resolve and validate a scalar conductivity field from a registry.
///
/// Values must be finite (already guaranteed by DenseFieldSnapshot) and
/// non-negative. Zero conductivity is permitted. Unit matching is exact and no
/// conversion or interpretation is performed.
[[nodiscard]] inline const DenseFieldSnapshot&
require_scalar_conductivity_field(
    const Topology& topology,
    const DenseFieldRegistry& registry,
    ScalarConductivityFieldContract contract) {
    conductivity_field_detail::validate_location(
        contract.location);
    conductivity_field_detail::require_text(
        contract.id,
        "mpmc::mesh::ScalarConductivityFieldContract: nonblank field id is required");
    conductivity_field_detail::require_text(
        contract.unit,
        "mpmc::mesh::ScalarConductivityFieldContract: explicit nonblank unit is required");

    const auto& field =
        registry.at(
            contract.location,
            contract.id);
    if (field.location() !=
            contract.location ||
        field.entity_count() !=
            topology.entity_count(
                contract.location) ||
        field.component_count() != 1U) {
        throw std::invalid_argument(
            "mpmc::mesh::require_scalar_conductivity_field: conductivity must be a scalar field aligned to the declared topology location");
    }
    if (field.metadata().unit !=
        contract.unit) {
        throw std::invalid_argument(
            "mpmc::mesh::require_scalar_conductivity_field: conductivity unit does not match the declared contract");
    }
    if (!std::all_of(
            field.values().begin(),
            field.values().end(),
            [](double value) {
                return value >= 0.0;
            })) {
        throw std::invalid_argument(
            "mpmc::mesh::require_scalar_conductivity_field: conductivity values must be non-negative");
    }
    return field;
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_CONDUCTIVITY_FIELD_HPP
