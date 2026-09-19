#ifndef MPMC_MESH_PERMEABILITY_TENSOR_3D_HPP
#define MPMC_MESH_PERMEABILITY_TENSOR_3D_HPP

#include <mpmc/mesh/dense_field.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::mesh {

/// Cell permeability tensor for the current baseline.
///
/// PERMX/PERMY/PERMZ are interpreted as Kxx/Kyy/Kzz in the existing mesh/world
/// Cartesian x/y/z basis. No local-axis rotation or off-diagonal permeability
/// components are inferred by this contract.
struct CartesianDiagonalPermeabilityTensor3D {
    double kxx_m2;
    double kyy_m2;
    double kzz_m2;
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

} // namespace permeability_tensor_3d_detail

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

} // namespace mpmc::mesh

#endif // MPMC_MESH_PERMEABILITY_TENSOR_3D_HPP
