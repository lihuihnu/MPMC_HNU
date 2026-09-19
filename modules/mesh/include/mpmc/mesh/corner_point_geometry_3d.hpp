#ifndef MPMC_MESH_CORNER_POINT_GEOMETRY_3D_HPP
#define MPMC_MESH_CORNER_POINT_GEOMETRY_3D_HPP

#include <mpmc/mesh/entity.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

struct Coordinate3D {
    double x_m;
    double y_m;
    double z_m;
};

/// Immutable cell-local corner-point geometry.
///
/// Each logical cell owns exactly eight corner vertices in the order
/// {LLL, HLL, LHL, HHL, LLH, HLH, LHH, HHH}, where I runs fastest.
/// Coordinates are SI metres. Cell volumes are cubic metres and may be zero
/// only when the importing contract explicitly retains a degenerate inactive
/// logical cell. No shared-face or shared-vertex identity is inferred here.
class CornerPointGeometry3D {
public:
    CornerPointGeometry3D(
        std::vector<Coordinate3D> vertex_coordinates_m,
        std::vector<double> cell_volumes_m3)
        : vertex_coordinates_m_(
              std::move(vertex_coordinates_m)),
          cell_volumes_m3_(
              std::move(cell_volumes_m3)) {
        validate();
    }

    CornerPointGeometry3D(
        const CornerPointGeometry3D&) = default;
    CornerPointGeometry3D(
        CornerPointGeometry3D&&) noexcept = default;
    CornerPointGeometry3D& operator=(
        const CornerPointGeometry3D&) = delete;
    CornerPointGeometry3D& operator=(
        CornerPointGeometry3D&&) = delete;
    ~CornerPointGeometry3D() = default;

    [[nodiscard]] std::size_t
    vertex_count() const noexcept {
        return vertex_coordinates_m_.size();
    }

    [[nodiscard]] std::size_t
    cell_count() const noexcept {
        return cell_volumes_m3_.size();
    }

    [[nodiscard]] std::span<const Coordinate3D>
    vertex_coordinates_m() const noexcept {
        return vertex_coordinates_m_;
    }

    [[nodiscard]] std::span<const double>
    cell_volumes_m3() const noexcept {
        return cell_volumes_m3_;
    }

    [[nodiscard]] Coordinate3D
    vertex_coordinate_m(LocalIndex vertex) const {
        return vertex_coordinates_m_.at(
            static_cast<std::size_t>(
                vertex.value()));
    }

    [[nodiscard]] double
    cell_volume_m3(LocalIndex cell) const {
        return cell_volumes_m3_.at(
            static_cast<std::size_t>(
                cell.value()));
    }

    [[nodiscard]] std::array<Coordinate3D, 8>
    cell_corners_m(LocalIndex cell) const {
        const std::size_t local =
            static_cast<std::size_t>(
                cell.value());
        if (local >= cell_count()) {
            throw std::out_of_range(
                "mpmc::mesh::CornerPointGeometry3D: cell index out of range");
        }
        const std::size_t begin =
            local * std::size_t{8U};
        return {
            vertex_coordinates_m_[begin],
            vertex_coordinates_m_[begin + 1U],
            vertex_coordinates_m_[begin + 2U],
            vertex_coordinates_m_[begin + 3U],
            vertex_coordinates_m_[begin + 4U],
            vertex_coordinates_m_[begin + 5U],
            vertex_coordinates_m_[begin + 6U],
            vertex_coordinates_m_[begin + 7U]};
    }

private:
    void validate() const {
        if (cell_volumes_m3_.size() >
            std::numeric_limits<std::size_t>::max() /
                std::size_t{8U}) {
            throw std::length_error(
                "mpmc::mesh::CornerPointGeometry3D: vertex count overflow");
        }
        if (vertex_coordinates_m_.size() !=
            cell_volumes_m3_.size() *
                std::size_t{8U}) {
            throw std::invalid_argument(
                "mpmc::mesh::CornerPointGeometry3D: exactly eight vertices per cell are required");
        }
        for (const auto coordinate :
             vertex_coordinates_m_) {
            if (!std::isfinite(coordinate.x_m) ||
                !std::isfinite(coordinate.y_m) ||
                !std::isfinite(coordinate.z_m)) {
                throw std::invalid_argument(
                    "mpmc::mesh::CornerPointGeometry3D: non-finite corner coordinate");
            }
        }
        for (const double volume :
             cell_volumes_m3_) {
            if (!std::isfinite(volume) ||
                volume < 0.0) {
                throw std::invalid_argument(
                    "mpmc::mesh::CornerPointGeometry3D: cell volume must be finite and nonnegative");
            }
        }
    }

    std::vector<Coordinate3D>
        vertex_coordinates_m_;
    std::vector<double> cell_volumes_m3_;
};

} // namespace mpmc::mesh

#endif // MPMC_MESH_CORNER_POINT_GEOMETRY_3D_HPP
