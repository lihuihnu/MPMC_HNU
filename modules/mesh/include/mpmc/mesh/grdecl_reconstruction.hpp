#ifndef MPMC_MESH_GRDECL_RECONSTRUCTION_HPP
#define MPMC_MESH_GRDECL_RECONSTRUCTION_HPP

#include <mpmc/mesh/mesh_exchange.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::mesh {

enum class GrdeclRepresentabilityDisposition : std::uint8_t {
    representable = 0,
    unsupported = 1,
};

struct GrdeclRepresentabilityIssue {
    std::string code;
    std::string message;
};

class GrdeclRepresentabilityReport {
public:
    [[nodiscard]] GrdeclRepresentabilityDisposition
    disposition() const noexcept {
        return disposition_;
    }

    [[nodiscard]] bool representable()
        const noexcept {
        return disposition_ ==
            GrdeclRepresentabilityDisposition::
                representable;
    }

    [[nodiscard]] const std::optional<
        std::array<std::size_t, 3>>&
    dimensions() const noexcept {
        return dimensions_;
    }

    [[nodiscard]] std::span<
        const GrdeclRepresentabilityIssue>
    issues() const noexcept {
        return issues_;
    }

    void set_dimensions(
        std::array<std::size_t, 3>
            dimensions) {
        dimensions_ = dimensions;
    }

    void note(
        std::string code,
        std::string message) {
        issues_.push_back(
            GrdeclRepresentabilityIssue{
                std::move(code),
                std::move(message)});
    }

    void reject(
        std::string code,
        std::string message) {
        disposition_ =
            GrdeclRepresentabilityDisposition::
                unsupported;
        issues_.push_back(
            GrdeclRepresentabilityIssue{
                std::move(code),
                std::move(message)});
    }

private:
    GrdeclRepresentabilityDisposition
        disposition_{
            GrdeclRepresentabilityDisposition::
                representable};
    std::optional<
        std::array<std::size_t, 3>>
        dimensions_;
    std::vector<
        GrdeclRepresentabilityIssue>
        issues_;
};

struct StructuredLogicalGridReconstruction3D {
    GrdeclRepresentabilityReport report;
    std::optional<LogicalCornerPointGrid3D>
        logical_grid;
    /// For each reconstructed I-fastest logical cell, source canonical local.
    std::vector<LocalIndex>
        source_cell_by_logical;
};

namespace grdecl_reconstruction_detail {

struct AxisClusters {
    std::vector<double> values;
    double tolerance;
};

[[nodiscard]] inline double
coordinate_tolerance(
    std::span<const Coordinate3D>
        coordinates) {
    double scale = 1.0;
    for (const auto value : coordinates) {
        scale = std::max(
            scale,
            std::max(
                {std::abs(value.x_m),
                 std::abs(value.y_m),
                 std::abs(value.z_m)}));
    }
    return 4096.0 *
           std::numeric_limits<double>::
               epsilon() *
           scale;
}

[[nodiscard]] inline AxisClusters
cluster_axis(
    std::vector<double> values,
    double tolerance) {
    std::sort(
        values.begin(),
        values.end());
    AxisClusters result{
        {},
        tolerance};
    for (const double value : values) {
        if (result.values.empty() ||
            std::abs(
                value -
                result.values.back()) >
                tolerance) {
            result.values.push_back(value);
        } else {
            // Use a deterministic representative that remains inside the
            // original floating-point cluster.
            result.values.back() =
                0.5 *
                (result.values.back() +
                 value);
        }
    }
    return result;
}

[[nodiscard]] inline std::optional<std::size_t>
axis_index(
    const AxisClusters& axis,
    double value) {
    const auto found =
        std::lower_bound(
            axis.values.begin(),
            axis.values.end(),
            value);
    std::optional<std::size_t> best;
    double best_distance =
        std::numeric_limits<double>::
            infinity();
    const auto consider =
        [&](auto iterator) {
            if (iterator ==
                axis.values.end()) {
                return;
            }
            const double distance =
                std::abs(
                    *iterator - value);
            if (distance <=
                    axis.tolerance &&
                distance <
                    best_distance) {
                best_distance =
                    distance;
                best =
                    static_cast<
                        std::size_t>(
                        std::distance(
                            axis.values.begin(),
                            iterator));
            }
        };
    consider(found);
    if (found !=
        axis.values.begin()) {
        consider(
            std::prev(found));
    }
    return best;
}

[[nodiscard]] inline std::size_t
checked_multiply(
    std::size_t left,
    std::size_t right,
    const char* message) {
    if (left != 0U &&
        right >
            std::numeric_limits<
                std::size_t>::max() /
                left) {
        throw std::length_error(message);
    }
    return left * right;
}

[[nodiscard]] inline std::size_t
logical_vertex(
    std::size_t i,
    std::size_t j,
    std::size_t k,
    std::size_t nx_vertices,
    std::size_t ny_vertices) {
    return i +
           nx_vertices *
               (j +
                ny_vertices * k);
}

[[nodiscard]] inline std::size_t
logical_cell(
    std::size_t i,
    std::size_t j,
    std::size_t k,
    std::size_t nx,
    std::size_t ny) {
    return i +
           nx *
               (j +
                ny * k);
}

[[nodiscard]] inline std::optional<
    std::reference_wrapper<
        const DenseFieldSnapshot>>
optional_scalar_field(
    const MeshExchangeDocument& document,
    std::string_view id,
    std::string_view unit,
    GrdeclRepresentabilityReport&
        report) {
    const auto* field =
        document.fields().find(
            EntityKind::cell,
            id);
    if (field == nullptr) {
        return std::nullopt;
    }
    if (field->component_count() != 1U ||
        field->metadata().unit != unit) {
        report.note(
            "grdecl_reconstruction.property_omitted",
            std::string{id} +
                " exists but does not match the scalar/unit contract and will not be projected");
        return std::nullopt;
    }
    return std::cref(*field);
}

[[nodiscard]] inline std::vector<double>
project_scalar_field(
    const DenseFieldSnapshot& field,
    std::span<const LocalIndex>
        source_cell_by_logical) {
    std::vector<double> values;
    values.reserve(
        source_cell_by_logical.size());
    for (const auto source :
         source_cell_by_logical) {
        values.push_back(
            field.value(
                source,
                0U));
    }
    return values;
}

} // namespace grdecl_reconstruction_detail

/// Detect whether a generic canonical mesh can be rebuilt as the current
/// GRDECL corner-point baseline without geometric approximation.
///
/// Baseline v1 intentionally requires a complete rectilinear tensor-product
/// hexahedral grid. Source vertex/cell ordering and stable IDs may be arbitrary,
/// but every coordinate must map to unique clustered x/y/z planes and every
/// logical cell must be present exactly once.
[[nodiscard]] inline StructuredLogicalGridReconstruction3D
reconstruct_structured_logical_grid_3d(
    const MeshExchangeDocument& document) {
    using namespace grdecl_reconstruction_detail;

    StructuredLogicalGridReconstruction3D
        result;

    if (document.logical_corner_point()
            .has_value()) {
        result.logical_grid =
            *document.logical_corner_point();
        result.report.set_dimensions(
            result.logical_grid->
                dimensions);
        return result;
    }

    if (document.dimension() != 3) {
        result.report.reject(
            "grdecl_reconstruction.requires_3d",
            "structured GRDECL reconstruction requires a 3D canonical mesh");
        return result;
    }

    const auto& topology =
        document.topology();
    if (!topology.has_relation(
            EntityKind::cell,
            EntityKind::vertex) ||
        topology.entity_count(
            EntityKind::cell) == 0U) {
        result.report.reject(
            "grdecl_reconstruction.cell_vertices_missing",
            "structured GRDECL reconstruction requires nonempty cell->vertex topology");
        return result;
    }

    const auto& cell_vertices =
        topology.relation(
            EntityKind::cell,
            EntityKind::vertex);
    for (std::size_t cell = 0U;
         cell < topology.entity_count(
             EntityKind::cell);
         ++cell) {
        if (cell_vertices.adjacent(
                LocalIndex{
                    static_cast<
                        LocalIndex::value_type>(
                        cell)})
                .size() != 8U) {
            result.report.reject(
                "grdecl_reconstruction.hexa_only",
                "structured GRDECL reconstruction currently requires every cell to be a linear hexahedron");
            return result;
        }
    }

    const auto coordinates =
        document.vertex_coordinates_m();
    if (coordinates.empty()) {
        result.report.reject(
            "grdecl_reconstruction.vertices_missing",
            "structured GRDECL reconstruction requires vertex coordinates");
        return result;
    }

    const double tolerance =
        coordinate_tolerance(
            coordinates);
    std::vector<double> x_values;
    std::vector<double> y_values;
    std::vector<double> z_values;
    x_values.reserve(
        coordinates.size());
    y_values.reserve(
        coordinates.size());
    z_values.reserve(
        coordinates.size());
    for (const auto coordinate :
         coordinates) {
        x_values.push_back(
            coordinate.x_m);
        y_values.push_back(
            coordinate.y_m);
        z_values.push_back(
            coordinate.z_m);
    }

    const auto x =
        cluster_axis(
            std::move(x_values),
            tolerance);
    const auto y =
        cluster_axis(
            std::move(y_values),
            tolerance);
    const auto z =
        cluster_axis(
            std::move(z_values),
            tolerance);

    if (x.values.size() < 2U ||
        y.values.size() < 2U ||
        z.values.size() < 2U) {
        result.report.reject(
            "grdecl_reconstruction.axis_cardinality",
            "structured GRDECL reconstruction requires at least two coordinate planes on each axis");
        return result;
    }

    const std::size_t nx =
        x.values.size() - 1U;
    const std::size_t ny =
        y.values.size() - 1U;
    const std::size_t nz =
        z.values.size() - 1U;
    const std::size_t expected_vertices =
        checked_multiply(
            checked_multiply(
                x.values.size(),
                y.values.size(),
                "grdecl reconstruction XY vertex count overflow"),
            z.values.size(),
            "grdecl reconstruction XYZ vertex count overflow");
    const std::size_t expected_cells =
        checked_multiply(
            checked_multiply(
                nx,
                ny,
                "grdecl reconstruction XY cell count overflow"),
            nz,
            "grdecl reconstruction XYZ cell count overflow");

    if (coordinates.size() !=
            expected_vertices ||
        topology.entity_count(
            EntityKind::cell) !=
            expected_cells) {
        result.report.reject(
            "grdecl_reconstruction.tensor_product_counts",
            "vertex/cell counts do not match a complete tensor-product structured grid");
        return result;
    }

    const std::size_t nx_vertices =
        x.values.size();
    const std::size_t ny_vertices =
        y.values.size();
    std::vector<std::optional<LocalIndex>>
        source_vertex_by_logical(
            expected_vertices);
    for (std::size_t source = 0U;
         source < coordinates.size();
         ++source) {
        const auto& coordinate =
            coordinates[source];
        const auto i =
            axis_index(
                x,
                coordinate.x_m);
        const auto j =
            axis_index(
                y,
                coordinate.y_m);
        const auto k =
            axis_index(
                z,
                coordinate.z_m);
        if (!i.has_value() ||
            !j.has_value() ||
            !k.has_value()) {
            result.report.reject(
                "grdecl_reconstruction.coordinate_plane_mapping",
                "one or more vertices cannot be mapped to reconstructed coordinate planes");
            return result;
        }
        const std::size_t logical =
            logical_vertex(
                *i,
                *j,
                *k,
                nx_vertices,
                ny_vertices);
        if (source_vertex_by_logical[
                logical]
                .has_value()) {
            result.report.reject(
                "grdecl_reconstruction.duplicate_tensor_vertex",
                "multiple source vertices map to one structured logical vertex");
            return result;
        }
        source_vertex_by_logical[
            logical] =
            LocalIndex{
                static_cast<
                    LocalIndex::value_type>(
                    source)};
    }
    if (std::any_of(
            source_vertex_by_logical.begin(),
            source_vertex_by_logical.end(),
            [](const auto& value) {
                return !value.has_value();
            })) {
        result.report.reject(
            "grdecl_reconstruction.tensor_vertex_hole",
            "the tensor-product vertex lattice contains missing vertices");
        return result;
    }

    result.source_cell_by_logical.assign(
        expected_cells,
        LocalIndex{
            std::numeric_limits<
                LocalIndex::value_type>::
                max()});
    std::vector<bool>
        logical_cell_seen(
            expected_cells,
            false);

    for (std::size_t source_cell = 0U;
         source_cell <
         topology.entity_count(
             EntityKind::cell);
         ++source_cell) {
        const auto vertices =
            cell_vertices.adjacent(
                LocalIndex{
                    static_cast<
                        LocalIndex::value_type>(
                        source_cell)});
        std::array<std::size_t, 8>
            is{};
        std::array<std::size_t, 8>
            js{};
        std::array<std::size_t, 8>
            ks{};

        for (std::size_t corner = 0U;
             corner < 8U;
             ++corner) {
            const auto coordinate =
                coordinates[
                    static_cast<std::size_t>(
                        vertices[corner].value())];
            const auto i =
                axis_index(
                    x,
                    coordinate.x_m);
            const auto j =
                axis_index(
                    y,
                    coordinate.y_m);
            const auto k =
                axis_index(
                    z,
                    coordinate.z_m);
            if (!i.has_value() ||
                !j.has_value() ||
                !k.has_value()) {
                result.report.reject(
                    "grdecl_reconstruction.cell_coordinate_mapping",
                    "cell vertex cannot be mapped to structured logical coordinates");
                return result;
            }
            is[corner] = *i;
            js[corner] = *j;
            ks[corner] = *k;
        }

        const auto [imin, imax] =
            std::minmax_element(
                is.begin(),
                is.end());
        const auto [jmin, jmax] =
            std::minmax_element(
                js.begin(),
                js.end());
        const auto [kmin, kmax] =
            std::minmax_element(
                ks.begin(),
                ks.end());
        if (*imax != *imin + 1U ||
            *jmax != *jmin + 1U ||
            *kmax != *kmin + 1U) {
            result.report.reject(
                "grdecl_reconstruction.non_adjacent_cell",
                "a hexahedron does not span exactly one adjacent logical interval on each axis");
            return result;
        }

        std::array<bool, 8>
            corners_seen{};
        for (std::size_t corner = 0U;
             corner < 8U;
             ++corner) {
            const std::size_t local_i =
                is[corner] - *imin;
            const std::size_t local_j =
                js[corner] - *jmin;
            const std::size_t local_k =
                ks[corner] - *kmin;
            const std::size_t slot =
                local_i +
                2U * local_j +
                4U * local_k;
            if (corners_seen[slot]) {
                result.report.reject(
                    "grdecl_reconstruction.cell_corner_duplicate",
                    "a hexahedron does not contain all eight logical corner combinations exactly once");
                return result;
            }
            corners_seen[slot] = true;
        }
        if (std::any_of(
                corners_seen.begin(),
                corners_seen.end(),
                [](bool seen) {
                    return !seen;
                })) {
            result.report.reject(
                "grdecl_reconstruction.cell_corner_hole",
                "a hexahedron is missing one or more structured logical corners");
            return result;
        }

        const std::size_t logical =
            logical_cell(
                *imin,
                *jmin,
                *kmin,
                nx,
                ny);
        if (logical_cell_seen[
                logical]) {
            result.report.reject(
                "grdecl_reconstruction.duplicate_logical_cell",
                "multiple source cells map to one structured logical cell");
            return result;
        }
        logical_cell_seen[
            logical] = true;
        result.source_cell_by_logical[
            logical] =
            LocalIndex{
                static_cast<
                    LocalIndex::value_type>(
                    source_cell)};
    }

    if (std::any_of(
            logical_cell_seen.begin(),
            logical_cell_seen.end(),
            [](bool seen) {
                return !seen;
            })) {
        result.report.reject(
            "grdecl_reconstruction.logical_cell_hole",
            "one or more structured logical cells are missing");
        result.source_cell_by_logical.clear();
        return result;
    }

    std::vector<double> coord_m;
    coord_m.reserve(
        checked_multiply(
            checked_multiply(
                nx_vertices,
                ny_vertices,
                "grdecl reconstruction pillar count overflow"),
            6U,
            "grdecl reconstruction COORD size overflow"));
    for (std::size_t j = 0U;
         j < ny_vertices;
         ++j) {
        for (std::size_t i = 0U;
             i < nx_vertices;
             ++i) {
            coord_m.push_back(
                x.values[i]);
            coord_m.push_back(
                y.values[j]);
            coord_m.push_back(
                z.values.front());
            coord_m.push_back(
                x.values[i]);
            coord_m.push_back(
                y.values[j]);
            coord_m.push_back(
                z.values.back());
        }
    }

    std::vector<double> zcorn_m(
        checked_multiply(
            expected_cells,
            8U,
            "grdecl reconstruction ZCORN size overflow"),
        0.0);
    const std::size_t doubled_nx =
        checked_multiply(
            nx,
            2U,
            "grdecl reconstruction 2*NX overflow");
    const std::size_t doubled_ny =
        checked_multiply(
            ny,
            2U,
            "grdecl reconstruction 2*NY overflow");
    for (std::size_t k = 0U;
         k < nz;
         ++k) {
        for (std::size_t j = 0U;
             j < ny;
             ++j) {
            for (std::size_t i = 0U;
                 i < nx;
                 ++i) {
                for (std::size_t local_k = 0U;
                     local_k < 2U;
                     ++local_k) {
                    for (std::size_t local_j = 0U;
                         local_j < 2U;
                         ++local_j) {
                        for (std::size_t local_i = 0U;
                             local_i < 2U;
                             ++local_i) {
                            const std::size_t z_index =
                                (2U * i + local_i) +
                                doubled_nx *
                                    ((2U * j + local_j) +
                                     doubled_ny *
                                         (2U * k + local_k));
                            zcorn_m[z_index] =
                                z.values[
                                    k + local_k];
                        }
                    }
                }
            }
        }
    }

    std::vector<std::uint8_t> active(
        expected_cells,
        std::uint8_t{1U});

    const auto poro =
        optional_scalar_field(
            document,
            "PORO",
            "1",
            result.report);
    const auto permx =
        optional_scalar_field(
            document,
            "PERMX",
            "m2",
            result.report);
    const auto permy =
        optional_scalar_field(
            document,
            "PERMY",
            "m2",
            result.report);
    const auto permz =
        optional_scalar_field(
            document,
            "PERMZ",
            "m2",
            result.report);

    result.logical_grid =
        LogicalCornerPointGrid3D{
            {nx, ny, nz},
            std::move(coord_m),
            std::move(zcorn_m),
            std::move(active),
            poro.has_value()
                ? project_scalar_field(
                      poro->get(),
                      result.source_cell_by_logical)
                : std::vector<double>{},
            permx.has_value()
                ? project_scalar_field(
                      permx->get(),
                      result.source_cell_by_logical)
                : std::vector<double>{},
            permy.has_value()
                ? project_scalar_field(
                      permy->get(),
                      result.source_cell_by_logical)
                : std::vector<double>{},
            permz.has_value()
                ? project_scalar_field(
                      permz->get(),
                      result.source_cell_by_logical)
                : std::vector<double>{},
            1.0,
            1.0};

    result.report.set_dimensions(
        {nx, ny, nz});
    return result;
}

[[nodiscard]] inline GrdeclRepresentabilityReport
detect_grdecl_representability(
    const MeshExchangeDocument& document) {
    return reconstruct_structured_logical_grid_3d(
               document)
        .report;
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_GRDECL_RECONSTRUCTION_HPP
