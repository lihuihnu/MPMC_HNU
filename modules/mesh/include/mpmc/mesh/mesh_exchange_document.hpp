#ifndef MPMC_MESH_MESH_EXCHANGE_DOCUMENT_HPP
#define MPMC_MESH_MESH_EXCHANGE_DOCUMENT_HPP

#include <mpmc/mesh/corner_point_geometry_3d.hpp>
#include <mpmc/mesh/dense_field_registry.hpp>
#include <mpmc/mesh/face_boundary.hpp>
#include <mpmc/mesh/topology.hpp>
#include <mpmc/mesh/mesh_exchange_group.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::mesh {

enum class MeshExchangeFormat : std::uint8_t {
    gmsh_4_1_ascii = 0,
    vtu_ascii = 1,
    grdecl = 2,
    face_based = 3,
};

enum class ConversionDisposition : std::uint8_t {
    lossless = 0,
    lossy = 1,
    unsupported = 2,
};

struct ConversionIssue {
    std::string code;
    std::string message;
};

class ConversionReport {
public:
    explicit ConversionReport(
        MeshExchangeFormat target_format)
        : target_format_(target_format) {}

    [[nodiscard]] MeshExchangeFormat
    target_format() const noexcept {
        return target_format_;
    }

    [[nodiscard]] ConversionDisposition
    disposition() const noexcept {
        return disposition_;
    }

    [[nodiscard]] bool lossless() const noexcept {
        return disposition_ ==
            ConversionDisposition::lossless;
    }

    [[nodiscard]] std::span<const ConversionIssue>
    issues() const noexcept {
        return issues_;
    }

    void note_lossy(
        std::string code,
        std::string message) {
        if (disposition_ ==
            ConversionDisposition::lossless) {
            disposition_ =
                ConversionDisposition::lossy;
        }
        issues_.push_back(
            ConversionIssue{
                std::move(code),
                std::move(message)});
    }

    void note_unsupported(
        std::string code,
        std::string message) {
        disposition_ =
            ConversionDisposition::unsupported;
        issues_.push_back(
            ConversionIssue{
                std::move(code),
                std::move(message)});
    }

private:
    MeshExchangeFormat target_format_;
    ConversionDisposition disposition_{
        ConversionDisposition::lossless};
    std::vector<ConversionIssue> issues_;
};



/// Source-preserving logical corner-point semantics.
///
/// COORD and ZCORN are stored in canonical SI metres while permeability is
/// stored in square metres. PORO/PERM arrays are optional: geometry-only GRDECL
/// documents keep them empty instead of inventing material properties.
/// source_*_scale_to_si retains the explicit import scale so a GRDECL writer
/// can reproduce the source numeric unit convention without guessing
/// FIELD/METRIC semantics.
struct LogicalCornerPointGrid3D {
    std::array<std::size_t, 3> dimensions;
    std::vector<double> coord_m;
    std::vector<double> zcorn_m;
    std::vector<std::uint8_t> active;
    std::vector<double> porosity;
    std::vector<double> permx_m2;
    std::vector<double> permy_m2;
    std::vector<double> permz_m2;
    double source_coordinate_scale_to_m;
    double source_permeability_scale_to_m2;
};


class MeshExchangeDocument {
public:
    MeshExchangeDocument(
        const MeshExchangeDocument&) = default;
    MeshExchangeDocument(
        MeshExchangeDocument&&) noexcept = default;
    MeshExchangeDocument& operator=(
        const MeshExchangeDocument&) = delete;
    MeshExchangeDocument& operator=(
        MeshExchangeDocument&&) = delete;
    ~MeshExchangeDocument() = default;

    [[nodiscard]] static MeshExchangeDocument
    create(
        MeshExchangeFormat source_format,
        int dimension,
        Topology topology,
        std::vector<Coordinate3D>
            vertex_coordinates_m,
        std::optional<FaceBoundarySnapshot>
            face_boundary,
        std::vector<DenseFieldSnapshot> fields,
        std::vector<MeshExchangeGroup> groups,
        std::optional<LogicalCornerPointGrid3D>
            logical_corner_point) {
        if (dimension != 2 &&
            dimension != 3) {
            throw std::invalid_argument(
                "mpmc::mesh::MeshExchangeDocument: dimension must be 2 or 3");
        }
        if (vertex_coordinates_m.size() !=
            topology.entity_count(
                EntityKind::vertex)) {
            throw std::invalid_argument(
                "mpmc::mesh::MeshExchangeDocument: vertex coordinate count does not match topology");
        }
        for (const auto coordinate :
             vertex_coordinates_m) {
            if (!std::isfinite(coordinate.x_m) ||
                !std::isfinite(coordinate.y_m) ||
                !std::isfinite(coordinate.z_m)) {
                throw std::invalid_argument(
                    "mpmc::mesh::MeshExchangeDocument: vertex coordinates must be finite");
            }
            if (dimension == 2 &&
                coordinate.z_m != 0.0) {
                throw std::invalid_argument(
                    "mpmc::mesh::MeshExchangeDocument: 2D canonical coordinates must lie in z=0");
            }
        }

        if (face_boundary.has_value() &&
            face_boundary->face_count() !=
                topology.entity_count(
                    EntityKind::face)) {
            throw std::invalid_argument(
                "mpmc::mesh::MeshExchangeDocument: face boundary count does not match topology");
        }

        validate_mesh_exchange_groups(topology, groups);
        if (logical_corner_point.has_value()) {
            validate_corner_point(
                *logical_corner_point);
            if (dimension != 3) {
                throw std::invalid_argument(
                    "mpmc::mesh::MeshExchangeDocument: logical corner-point semantics require dimension 3");
            }
        }

        auto registry =
            DenseFieldRegistry::create(
                topology,
                std::move(fields));

        return MeshExchangeDocument{
            source_format,
            dimension,
            std::move(topology),
            std::move(vertex_coordinates_m),
            std::move(face_boundary),
            std::move(registry),
            std::move(groups),
            std::move(logical_corner_point)};
    }

    [[nodiscard]] MeshExchangeFormat
    source_format() const noexcept {
        return source_format_;
    }

    [[nodiscard]] int dimension()
        const noexcept {
        return dimension_;
    }

    [[nodiscard]] const Topology&
    topology() const noexcept {
        return topology_;
    }

    [[nodiscard]] std::span<const Coordinate3D>
    vertex_coordinates_m() const noexcept {
        return vertex_coordinates_m_;
    }

    [[nodiscard]] const std::optional<
        FaceBoundarySnapshot>&
    face_boundary() const noexcept {
        return face_boundary_;
    }

    [[nodiscard]] const DenseFieldRegistry&
    fields() const noexcept {
        return fields_;
    }

    [[nodiscard]] std::span<
        const MeshExchangeGroup>
    groups() const noexcept {
        return groups_;
    }

    [[nodiscard]] const std::optional<
        LogicalCornerPointGrid3D>&
    logical_corner_point() const noexcept {
        return logical_corner_point_;
    }

private:
    MeshExchangeDocument(
        MeshExchangeFormat source_format,
        int dimension,
        Topology topology,
        std::vector<Coordinate3D>
            vertex_coordinates_m,
        std::optional<FaceBoundarySnapshot>
            face_boundary,
        DenseFieldRegistry fields,
        std::vector<MeshExchangeGroup> groups,
        std::optional<LogicalCornerPointGrid3D>
            logical_corner_point)
        : source_format_(source_format),
          dimension_(dimension),
          topology_(std::move(topology)),
          vertex_coordinates_m_(
              std::move(vertex_coordinates_m)),
          face_boundary_(
              std::move(face_boundary)),
          fields_(std::move(fields)),
          groups_(std::move(groups)),
          logical_corner_point_(
              std::move(logical_corner_point)) {}

    [[nodiscard]] static std::size_t
    checked_add(
        std::size_t left,
        std::size_t right,
        const char* message) {
        if (left >
            std::numeric_limits<std::size_t>::max() -
                right) {
            throw std::length_error(message);
        }
        return left + right;
    }

    [[nodiscard]] static std::size_t
    checked_multiply(
        std::size_t left,
        std::size_t right,
        const char* message) {
        if (left != 0U &&
            right >
                std::numeric_limits<std::size_t>::max() /
                    left) {
            throw std::length_error(message);
        }
        return left * right;
    }


    static void validate_corner_point(
        const LogicalCornerPointGrid3D& data) {
        const auto nx = data.dimensions[0];
        const auto ny = data.dimensions[1];
        const auto nz = data.dimensions[2];
        if (nx == 0U ||
            ny == 0U ||
            nz == 0U) {
            throw std::invalid_argument(
                "mpmc::mesh::MeshExchangeDocument: corner-point dimensions must be positive");
        }

        const std::size_t nxp1 =
            checked_add(
                nx,
                1U,
                "mpmc::mesh::MeshExchangeDocument: corner-point NX+1 overflow");
        const std::size_t nyp1 =
            checked_add(
                ny,
                1U,
                "mpmc::mesh::MeshExchangeDocument: corner-point NY+1 overflow");
        const std::size_t cell_count =
            checked_multiply(
                checked_multiply(
                    nx,
                    ny,
                    "mpmc::mesh::MeshExchangeDocument: corner-point XY count overflow"),
                nz,
                "mpmc::mesh::MeshExchangeDocument: corner-point cell count overflow");
        const std::size_t pillar_count =
            checked_multiply(
                nxp1,
                nyp1,
                "mpmc::mesh::MeshExchangeDocument: corner-point pillar count overflow");
        const std::size_t expected_coord =
            checked_multiply(
                pillar_count,
                6U,
                "mpmc::mesh::MeshExchangeDocument: corner-point COORD size overflow");
        const std::size_t expected_zcorn =
            checked_multiply(
                cell_count,
                8U,
                "mpmc::mesh::MeshExchangeDocument: corner-point ZCORN size overflow");

        const auto optional_cell_array =
            [cell_count](std::size_t size) {
                return size == 0U ||
                       size == cell_count;
            };
        if (data.coord_m.size() !=
                expected_coord ||
            data.zcorn_m.size() !=
                expected_zcorn ||
            data.active.size() !=
                cell_count ||
            !optional_cell_array(
                data.porosity.size()) ||
            !optional_cell_array(
                data.permx_m2.size()) ||
            !optional_cell_array(
                data.permy_m2.size()) ||
            !optional_cell_array(
                data.permz_m2.size())) {
            throw std::invalid_argument(
                "mpmc::mesh::MeshExchangeDocument: corner-point geometry arrays must match dimensions and optional property arrays must be empty or cell-aligned");
        }

        const auto finite =
            [](double value) {
                return std::isfinite(value);
            };
        if (!std::all_of(
                data.coord_m.begin(),
                data.coord_m.end(),
                finite) ||
            !std::all_of(
                data.zcorn_m.begin(),
                data.zcorn_m.end(),
                finite) ||
            !std::all_of(
                data.porosity.begin(),
                data.porosity.end(),
                finite) ||
            !std::all_of(
                data.permx_m2.begin(),
                data.permx_m2.end(),
                finite) ||
            !std::all_of(
                data.permy_m2.begin(),
                data.permy_m2.end(),
                finite) ||
            !std::all_of(
                data.permz_m2.begin(),
                data.permz_m2.end(),
                finite)) {
            throw std::invalid_argument(
                "mpmc::mesh::MeshExchangeDocument: corner-point numeric arrays must be finite");
        }

        for (const auto value : data.active) {
            if (value != std::uint8_t{0U} &&
                value != std::uint8_t{1U}) {
                throw std::invalid_argument(
                    "mpmc::mesh::MeshExchangeDocument: corner-point ACTNUM must be 0/1");
            }
        }
        for (const auto value :
             data.porosity) {
            if (value < 0.0 ||
                value > 1.0) {
                throw std::invalid_argument(
                    "mpmc::mesh::MeshExchangeDocument: corner-point porosity must lie in [0,1]");
            }
        }
        for (const auto* values :
             {&data.permx_m2,
              &data.permy_m2,
              &data.permz_m2}) {
            if (!std::all_of(
                    values->begin(),
                    values->end(),
                    [](double value) {
                        return value >= 0.0;
                    })) {
                throw std::invalid_argument(
                    "mpmc::mesh::MeshExchangeDocument: corner-point permeability must be non-negative");
            }
        }

        if (!std::isfinite(
                data.source_coordinate_scale_to_m) ||
            data.source_coordinate_scale_to_m <=
                0.0 ||
            !std::isfinite(
                data.source_permeability_scale_to_m2) ||
            data.source_permeability_scale_to_m2 <=
                0.0) {
            throw std::invalid_argument(
                "mpmc::mesh::MeshExchangeDocument: source unit scales must be finite and positive");
        }
    }

    MeshExchangeFormat source_format_;
    int dimension_;
    Topology topology_;
    std::vector<Coordinate3D>
        vertex_coordinates_m_;
    std::optional<FaceBoundarySnapshot>
        face_boundary_;
    DenseFieldRegistry fields_;
    std::vector<MeshExchangeGroup> groups_;
    std::optional<LogicalCornerPointGrid3D>
        logical_corner_point_;
};

} // namespace mpmc::mesh

#endif // MPMC_MESH_MESH_EXCHANGE_DOCUMENT_HPP
