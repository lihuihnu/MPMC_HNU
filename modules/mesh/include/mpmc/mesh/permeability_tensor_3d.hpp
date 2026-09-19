#ifndef MPMC_MESH_PERMEABILITY_TENSOR_3D_HPP
#define MPMC_MESH_PERMEABILITY_TENSOR_3D_HPP

#include <mpmc/mesh/dense_field.hpp>
#include <mpmc/mesh/dense_field_registry.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::mesh {

/// Explicit basis for material tensors represented by this mesh contract.
enum class CartesianTensorBasis3D : std::uint8_t {
    mesh_world_xyz = 0,
};

/// Cell permeability tensor for the existing diagonal baseline.
///
/// PERMX/PERMY/PERMZ are interpreted as Kxx/Kyy/Kzz in the mesh/world
/// Cartesian x/y/z basis.
struct CartesianDiagonalPermeabilityTensor3D {
    double kxx_m2;
    double kyy_m2;
    double kzz_m2;
};

/// Full symmetric 3D permeability tensor in the mesh/world Cartesian basis.
///
/// Component order used by the dense-field adapter is:
/// [Kxx, Kyy, Kzz, Kxy, Kxz, Kyz], all in square metres. Symmetry is explicit
/// in the type, so K_yx=K_xy, K_zx=K_xz and K_zy=K_yz are not duplicated.
struct CartesianSymmetricPermeabilityTensor3D {
    double kxx_m2;
    double kyy_m2;
    double kzz_m2;
    double kxy_m2;
    double kxz_m2;
    double kyz_m2;
};

/// Immutable cell-aligned diagonal permeability tensor snapshot.
class CellCartesianDiagonalPermeability3D {
public:
    explicit CellCartesianDiagonalPermeability3D(
        std::vector<CartesianDiagonalPermeabilityTensor3D> tensors)
        : tensors_(std::move(tensors)) {
        for (const auto tensor : tensors_) {
            if (!std::isfinite(tensor.kxx_m2) ||
                !std::isfinite(tensor.kyy_m2) ||
                !std::isfinite(tensor.kzz_m2) ||
                tensor.kxx_m2 < 0.0 ||
                tensor.kyy_m2 < 0.0 ||
                tensor.kzz_m2 < 0.0) {
                throw std::invalid_argument(
                    "mpmc::mesh::CellCartesianDiagonalPermeability3D: tensor components must be finite and non-negative");
            }
        }
    }

    CellCartesianDiagonalPermeability3D(
        const CellCartesianDiagonalPermeability3D&) = default;
    CellCartesianDiagonalPermeability3D(
        CellCartesianDiagonalPermeability3D&&) noexcept = default;
    CellCartesianDiagonalPermeability3D& operator=(
        const CellCartesianDiagonalPermeability3D&) = delete;
    CellCartesianDiagonalPermeability3D& operator=(
        CellCartesianDiagonalPermeability3D&&) = delete;
    ~CellCartesianDiagonalPermeability3D() = default;

    [[nodiscard]] static constexpr CartesianTensorBasis3D
    basis() noexcept {
        return CartesianTensorBasis3D::mesh_world_xyz;
    }

    [[nodiscard]] std::size_t cell_count() const noexcept {
        return tensors_.size();
    }

    [[nodiscard]] CartesianDiagonalPermeabilityTensor3D
    tensor(LocalIndex cell) const {
        return tensors_.at(
            static_cast<std::size_t>(cell.value()));
    }

private:
    std::vector<CartesianDiagonalPermeabilityTensor3D> tensors_;
};

namespace permeability_tensor_3d_detail {

inline void require_permeability_field(
    const DenseFieldSnapshot& field,
    const Topology& topology,
    std::string_view expected_id) {
    if (field.location() != EntityKind::cell ||
        field.component_count() != 1U ||
        field.entity_count() !=
            topology.entity_count(EntityKind::cell)) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cell_cartesian_diagonal_permeability_3d: permeability fields must be scalar cell fields aligned to topology");
    }
    if (field.metadata().id != expected_id ||
        field.metadata().unit != "m2") {
        throw std::invalid_argument(
            "mpmc::mesh::make_cell_cartesian_diagonal_permeability_3d: expected PERMX/PERMY/PERMZ with SI unit m2");
    }
}

inline void require_symmetric_permeability_field(
    const DenseFieldSnapshot& field,
    const Topology& topology) {
    if (field.location() != EntityKind::cell ||
        field.component_count() != 6U ||
        field.entity_count() !=
            topology.entity_count(EntityKind::cell)) {
        throw std::invalid_argument(
            "mpmc::mesh::make_cell_cartesian_symmetric_permeability_3d: field must be a six-component cell field aligned to topology");
    }
    if (field.metadata().unit != "m2") {
        throw std::invalid_argument(
            "mpmc::mesh::make_cell_cartesian_symmetric_permeability_3d: permeability tensor unit must be m2");
    }
}

inline void require_positive_semidefinite(
    const CartesianSymmetricPermeabilityTensor3D& tensor) {
    const double scale = std::max(
        {std::abs(tensor.kxx_m2),
         std::abs(tensor.kyy_m2),
         std::abs(tensor.kzz_m2),
         std::abs(tensor.kxy_m2),
         std::abs(tensor.kxz_m2),
         std::abs(tensor.kyz_m2)});

    if (!std::isfinite(scale)) {
        throw std::invalid_argument(
            "mpmc::mesh::CellCartesianSymmetricPermeability3D: tensor components must be finite");
    }
    if (tensor.kxx_m2 < 0.0 ||
        tensor.kyy_m2 < 0.0 ||
        tensor.kzz_m2 < 0.0) {
        throw std::invalid_argument(
            "mpmc::mesh::CellCartesianSymmetricPermeability3D: diagonal permeability components must be non-negative");
    }
    if (scale == 0.0) {
        return;
    }

    const double a = tensor.kxx_m2 / scale;
    const double b = tensor.kyy_m2 / scale;
    const double c = tensor.kzz_m2 / scale;
    const double d = tensor.kxy_m2 / scale;
    const double e = tensor.kxz_m2 / scale;
    const double f = tensor.kyz_m2 / scale;

    const double minor_xy = a * b - d * d;
    const double minor_xz = a * c - e * e;
    const double minor_yz = b * c - f * f;
    const double determinant =
        a * b * c +
        2.0 * d * e * f -
        a * f * f -
        b * e * e -
        c * d * d;

    constexpr double tolerance =
        4096.0 *
        std::numeric_limits<double>::epsilon();
    if (!std::isfinite(minor_xy) ||
        !std::isfinite(minor_xz) ||
        !std::isfinite(minor_yz) ||
        !std::isfinite(determinant) ||
        minor_xy < -tolerance ||
        minor_xz < -tolerance ||
        minor_yz < -tolerance ||
        determinant < -tolerance) {
        throw std::invalid_argument(
            "mpmc::mesh::CellCartesianSymmetricPermeability3D: permeability tensor must be positive semidefinite");
    }
}

} // namespace permeability_tensor_3d_detail

/// Immutable cell-aligned full symmetric permeability tensor snapshot.
class CellCartesianSymmetricPermeability3D {
public:
    explicit CellCartesianSymmetricPermeability3D(
        std::vector<CartesianSymmetricPermeabilityTensor3D> tensors)
        : tensors_(std::move(tensors)) {
        for (const auto& tensor : tensors_) {
            permeability_tensor_3d_detail::
                require_positive_semidefinite(tensor);
        }
    }

    CellCartesianSymmetricPermeability3D(
        const CellCartesianSymmetricPermeability3D&) = default;
    CellCartesianSymmetricPermeability3D(
        CellCartesianSymmetricPermeability3D&&) noexcept = default;
    CellCartesianSymmetricPermeability3D& operator=(
        const CellCartesianSymmetricPermeability3D&) = delete;
    CellCartesianSymmetricPermeability3D& operator=(
        CellCartesianSymmetricPermeability3D&&) = delete;
    ~CellCartesianSymmetricPermeability3D() = default;

    [[nodiscard]] static constexpr CartesianTensorBasis3D
    basis() noexcept {
        return CartesianTensorBasis3D::mesh_world_xyz;
    }

    [[nodiscard]] std::size_t cell_count() const noexcept {
        return tensors_.size();
    }

    [[nodiscard]] CartesianSymmetricPermeabilityTensor3D
    tensor(LocalIndex cell) const {
        return tensors_.at(
            static_cast<std::size_t>(cell.value()));
    }

private:
    std::vector<CartesianSymmetricPermeabilityTensor3D> tensors_;
};

[[nodiscard]] inline CellCartesianDiagonalPermeability3D
make_cell_cartesian_diagonal_permeability_3d(
    const Topology& topology,
    const DenseFieldSnapshot& permx,
    const DenseFieldSnapshot& permy,
    const DenseFieldSnapshot& permz) {
    permeability_tensor_3d_detail::require_permeability_field(
        permx, topology, "PERMX");
    permeability_tensor_3d_detail::require_permeability_field(
        permy, topology, "PERMY");
    permeability_tensor_3d_detail::require_permeability_field(
        permz, topology, "PERMZ");

    const std::size_t cell_count =
        topology.entity_count(EntityKind::cell);
    std::vector<CartesianDiagonalPermeabilityTensor3D>
        tensors;
    tensors.reserve(cell_count);
    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        const LocalIndex local{
            static_cast<LocalIndex::value_type>(cell)};
        tensors.push_back(
            CartesianDiagonalPermeabilityTensor3D{
                permx.value(local, 0U),
                permy.value(local, 0U),
                permz.value(local, 0U)});
    }

    return CellCartesianDiagonalPermeability3D{
        std::move(tensors)};
}

/// Materialize a full symmetric permeability tensor from one six-component
/// DenseFieldSnapshot. Component order is [Kxx,Kyy,Kzz,Kxy,Kxz,Kyz].
[[nodiscard]] inline CellCartesianSymmetricPermeability3D
make_cell_cartesian_symmetric_permeability_3d(
    const Topology& topology,
    const DenseFieldSnapshot& field) {
    permeability_tensor_3d_detail::
        require_symmetric_permeability_field(
            field,
            topology);

    const std::size_t cell_count =
        topology.entity_count(EntityKind::cell);
    std::vector<CartesianSymmetricPermeabilityTensor3D>
        tensors;
    tensors.reserve(cell_count);

    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        const LocalIndex local{
            static_cast<LocalIndex::value_type>(cell)};
        const auto components =
            field.entity_values(local);
        tensors.push_back(
            CartesianSymmetricPermeabilityTensor3D{
                components[0],
                components[1],
                components[2],
                components[3],
                components[4],
                components[5]});
    }

    return CellCartesianSymmetricPermeability3D{
        std::move(tensors)};
}

/// Resolve a named cell tensor from a location-aware registry and materialize
/// the typed symmetric permeability snapshot.
[[nodiscard]] inline CellCartesianSymmetricPermeability3D
make_cell_cartesian_symmetric_permeability_3d(
    const Topology& topology,
    const DenseFieldRegistry& registry,
    std::string_view field_id) {
    return make_cell_cartesian_symmetric_permeability_3d(
        topology,
        registry.at(EntityKind::cell, field_id));
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_PERMEABILITY_TENSOR_3D_HPP
