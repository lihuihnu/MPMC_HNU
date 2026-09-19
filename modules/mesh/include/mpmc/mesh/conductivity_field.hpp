#ifndef MPMC_MESH_CONDUCTIVITY_FIELD_HPP
#define MPMC_MESH_CONDUCTIVITY_FIELD_HPP

#include <mpmc/mesh/cartesian_symmetric_tensor_3d.hpp>
#include <mpmc/mesh/dense_field.hpp>
#include <mpmc/mesh/dense_field_registry.hpp>
#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

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

/// Contract for a full symmetric 3D conductivity tensor field.
///
/// The six components are interpreted in the mesh/world Cartesian basis using
/// the fixed order [xx, yy, zz, xy, xz, yz]. The physical interpretation
/// (thermal, electrical, or another conductivity) remains caller-defined
/// through the field id and exact unit string.
struct SymmetricConductivityFieldContract3D {
    EntityKind location;
    std::string id;
    std::string unit;
};

/// One full symmetric 3D conductivity tensor.
///
/// Values intentionally carry no hard-coded physical unit suffix because this
/// mesh layer does not decide whether conductivity is thermal, electrical, or
/// another transport coefficient. The associated typed field snapshot stores
/// and exposes the exact caller-declared unit.
struct CartesianSymmetricConductivityTensor3D {
    double xx;
    double yy;
    double zz;
    double xy;
    double xz;
    double yz;
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
        "mpmc::mesh::conductivity field contract: invalid field location");
}

inline void require_positive_semidefinite(
    const CartesianSymmetricConductivityTensor3D& tensor) {
    const auto status =
        classify_positive_semidefinite(
            CartesianSymmetricTensorComponents3D{
                tensor.xx,
                tensor.yy,
                tensor.zz,
                tensor.xy,
                tensor.xz,
                tensor.yz});
    switch (status) {
    case PositiveSemidefiniteStatus3D::valid:
        return;
    case PositiveSemidefiniteStatus3D::non_finite:
        throw std::invalid_argument(
            "mpmc::mesh::CartesianSymmetricConductivityField3D: tensor components must be finite");
    case PositiveSemidefiniteStatus3D::negative_diagonal:
        throw std::invalid_argument(
            "mpmc::mesh::CartesianSymmetricConductivityField3D: diagonal conductivity components must be non-negative");
    case PositiveSemidefiniteStatus3D::indefinite:
        throw std::invalid_argument(
            "mpmc::mesh::CartesianSymmetricConductivityField3D: conductivity tensor must be positive semidefinite");
    }
    throw std::logic_error(
        "mpmc::mesh::CartesianSymmetricConductivityField3D: invalid PSD classifier status");
}

} // namespace conductivity_field_detail

/// Resolve and validate a scalar conductivity field from a registry.
///
/// Values must be finite (already guaranteed by DenseFieldSnapshot) and
/// non-negative. Zero conductivity is permitted. Unit matching is exact and no
/// conversion or interpretation is performed.
[[nodiscard]] inline std::reference_wrapper<const DenseFieldSnapshot>
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
    return std::cref(field);
}

/// Immutable typed view of a location-aligned full symmetric 3D conductivity
/// field. The basis is always mesh/world Cartesian xyz; the unit is preserved
/// exactly as declared by the consuming contract.
class CartesianSymmetricConductivityField3D {
public:
    CartesianSymmetricConductivityField3D(
        EntityKind location,
        std::string unit,
        std::vector<CartesianSymmetricConductivityTensor3D> tensors)
        : location_(location),
          unit_(std::move(unit)),
          tensors_(std::move(tensors)) {
        conductivity_field_detail::validate_location(
            location_);
        conductivity_field_detail::require_text(
            unit_,
            "mpmc::mesh::CartesianSymmetricConductivityField3D: explicit nonblank unit is required");
        for (const auto& tensor : tensors_) {
            conductivity_field_detail::
                require_positive_semidefinite(
                    tensor);
        }
    }

    CartesianSymmetricConductivityField3D(
        const CartesianSymmetricConductivityField3D&) = default;
    CartesianSymmetricConductivityField3D(
        CartesianSymmetricConductivityField3D&&) noexcept = default;
    CartesianSymmetricConductivityField3D& operator=(
        const CartesianSymmetricConductivityField3D&) = delete;
    CartesianSymmetricConductivityField3D& operator=(
        CartesianSymmetricConductivityField3D&&) = delete;
    ~CartesianSymmetricConductivityField3D() = default;

    [[nodiscard]] static constexpr CartesianTensorBasis3D
    basis() noexcept {
        return CartesianTensorBasis3D::mesh_world_xyz;
    }

    [[nodiscard]] EntityKind location() const noexcept {
        return location_;
    }

    [[nodiscard]] std::string_view unit() const noexcept {
        return unit_;
    }

    [[nodiscard]] std::size_t entity_count() const noexcept {
        return tensors_.size();
    }

    [[nodiscard]] CartesianSymmetricConductivityTensor3D
    tensor(LocalIndex entity) const {
        return tensors_.at(
            static_cast<std::size_t>(
                entity.value()));
    }

private:
    EntityKind location_;
    std::string unit_;
    std::vector<CartesianSymmetricConductivityTensor3D> tensors_;
};

/// Resolve and materialize a full symmetric 3D conductivity tensor field.
///
/// Dense-field component order is [xx, yy, zz, xy, xz, yz]. Unit matching is
/// exact; no unit conversion or thermal/electrical interpretation is performed.
/// The tensor must be positive semidefinite at every entity. Zero principal
/// conductivity directions are allowed.
[[nodiscard]] inline CartesianSymmetricConductivityField3D
make_cartesian_symmetric_conductivity_field_3d(
    const Topology& topology,
    const DenseFieldRegistry& registry,
    SymmetricConductivityFieldContract3D contract) {
    conductivity_field_detail::validate_location(
        contract.location);
    conductivity_field_detail::require_text(
        contract.id,
        "mpmc::mesh::SymmetricConductivityFieldContract3D: nonblank field id is required");
    conductivity_field_detail::require_text(
        contract.unit,
        "mpmc::mesh::SymmetricConductivityFieldContract3D: explicit nonblank unit is required");

    const auto& field =
        registry.at(
            contract.location,
            contract.id);
    if (field.location() !=
            contract.location ||
        field.entity_count() !=
            topology.entity_count(
                contract.location) ||
        field.component_count() != 6U) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cartesian_symmetric_conductivity_field_3d: conductivity must be a six-component field aligned to the declared topology location");
    }
    if (field.metadata().unit !=
        contract.unit) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cartesian_symmetric_conductivity_field_3d: conductivity unit does not match the declared contract");
    }

    std::vector<CartesianSymmetricConductivityTensor3D>
        tensors;
    tensors.reserve(
        field.entity_count());
    for (std::size_t entity = 0U;
         entity < field.entity_count();
         ++entity) {
        const auto components =
            field.entity_values(
                LocalIndex{
                    static_cast<
                        LocalIndex::value_type>(
                        entity)});
        tensors.push_back(
            CartesianSymmetricConductivityTensor3D{
                components[0],
                components[1],
                components[2],
                components[3],
                components[4],
                components[5]});
    }

    return CartesianSymmetricConductivityField3D{
        contract.location,
        std::move(contract.unit),
        std::move(tensors)};
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_CONDUCTIVITY_FIELD_HPP
