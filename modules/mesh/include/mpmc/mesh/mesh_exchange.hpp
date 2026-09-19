#ifndef MPMC_MESH_MESH_EXCHANGE_HPP
#define MPMC_MESH_MESH_EXCHANGE_HPP

#include <mpmc/mesh/corner_point_geometry_3d.hpp>
#include <mpmc/mesh/dense_field_registry.hpp>
#include <mpmc/mesh/face_boundary.hpp>
#include <mpmc/mesh/topology.hpp>

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

struct MeshExchangeGroup {
    EntityKind location;
    std::uint32_t tag;
    std::string name;
    std::vector<GlobalEntityId> members;
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

        validate_groups(topology, groups);
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

    static void validate_groups(
        const Topology& topology,
        const std::vector<MeshExchangeGroup>&
            groups) {
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
            for (const auto member :
                 group.members) {
                if (std::find(
                        ids.begin(),
                        ids.end(),
                        member) == ids.end()) {
                    throw std::invalid_argument(
                        "mpmc::mesh::MeshExchangeDocument: group member is absent from topology");
                }
            }
        }
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

namespace mesh_exchange_detail {

[[nodiscard]] inline bool has_supported_generic_topology(
    const MeshExchangeDocument& document) {
    const auto& topology =
        document.topology();
    if (topology.entity_count(
            EntityKind::vertex) == 0U ||
        topology.entity_count(
            EntityKind::cell) == 0U ||
        topology.entity_count(
            EntityKind::edge) != 0U ||
        !topology.has_relation(
            EntityKind::cell,
            EntityKind::vertex)) {
        return false;
    }

    const auto& cell_vertices =
        topology.relation(
            EntityKind::cell,
            EntityKind::vertex);
    for (std::size_t cell = 0U;
         cell < topology.entity_count(
             EntityKind::cell);
         ++cell) {
        const auto vertices =
            cell_vertices.adjacent(
                LocalIndex{
                    static_cast<
                        LocalIndex::value_type>(
                        cell)});
        if (document.dimension() == 2) {
            if (vertices.size() != 3U &&
                vertices.size() != 4U) {
                return false;
            }
        } else if (
            vertices.size() != 4U &&
            vertices.size() != 5U &&
            vertices.size() != 6U &&
            vertices.size() != 8U) {
            return false;
        }
    }

    if (topology.entity_count(
            EntityKind::face) == 0U ||
        !topology.has_relation(
            EntityKind::cell,
            EntityKind::face) ||
        !topology.has_relation(
            EntityKind::face,
            EntityKind::vertex) ||
        !topology.has_relation(
            EntityKind::face,
            EntityKind::cell)) {
        return false;
    }

    const auto& face_vertices =
        topology.relation(
            EntityKind::face,
            EntityKind::vertex);
    for (std::size_t face = 0U;
         face < topology.entity_count(
             EntityKind::face);
         ++face) {
        const auto vertices =
            face_vertices.adjacent(
                LocalIndex{
                    static_cast<
                        LocalIndex::value_type>(
                        face)});
        if (document.dimension() == 2) {
            if (vertices.size() != 2U) {
                return false;
            }
        } else if (
            vertices.size() != 3U &&
            vertices.size() != 4U) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] inline bool has_tagged_faces(
    const MeshExchangeDocument& document) {
    if (!document.face_boundary()
             .has_value()) {
        return false;
    }
    const auto tags =
        document.face_boundary()->
            physical_tags();
    return std::any_of(
        tags.begin(),
        tags.end(),
        [](PhysicalTag tag) {
            return tag.is_tagged();
        });
}

[[nodiscard]] inline bool gmsh_ids_require_remap(
    const Topology& topology) {
    for (const auto id :
         topology.global_ids(
             EntityKind::vertex)) {
        if (id.value() == 0U) {
            return true;
        }
    }

    std::vector<
        GlobalEntityId::value_type>
        element_ids;
    const auto faces =
        topology.global_ids(
            EntityKind::face);
    const auto cells =
        topology.global_ids(
            EntityKind::cell);
    element_ids.reserve(
        faces.size() + cells.size());
    for (const auto id : faces) {
        if (id.value() == 0U) {
            return true;
        }
        element_ids.push_back(id.value());
    }
    for (const auto id : cells) {
        if (id.value() == 0U) {
            return true;
        }
        element_ids.push_back(id.value());
    }
    std::sort(
        element_ids.begin(),
        element_ids.end());
    return std::adjacent_find(
               element_ids.begin(),
               element_ids.end()) !=
           element_ids.end();
}

} // namespace mesh_exchange_detail

[[nodiscard]] inline ConversionReport
analyze_conversion(
    const MeshExchangeDocument& document,
    MeshExchangeFormat target_format) {
    using namespace mesh_exchange_detail;

    ConversionReport report{target_format};

    if (target_format ==
        MeshExchangeFormat::grdecl) {
        if (document.dimension() != 3) {
            report.note_unsupported(
                "grdecl.requires_3d",
                "GRDECL baseline requires a 3D document");
            return report;
        }
        if (!document.logical_corner_point()
                 .has_value()) {
            report.note_unsupported(
                "grdecl.logical_corner_point_missing",
                "GRDECL export requires preserved logical corner-point semantics; arbitrary linear topology is not reverse-engineered");
            return report;
        }

        if (!document.groups().empty()) {
            report.note_lossy(
                "grdecl.groups_not_representable",
                "minimal GRDECL baseline does not serialize generic named/physical groups");
        }
        if (has_tagged_faces(document)) {
            report.note_lossy(
                "grdecl.face_tags_not_representable",
                "minimal GRDECL baseline does not serialize generic face physical tags");
        }

        for (const auto& field :
             document.fields().fields()) {
            const bool supported =
                field.location() ==
                    EntityKind::cell &&
                (field.metadata().id == "PORO" ||
                 field.metadata().id == "PERMX" ||
                 field.metadata().id == "PERMY" ||
                 field.metadata().id == "PERMZ");
            if (!supported) {
                report.note_lossy(
                    "grdecl.field_not_representable",
                    "minimal GRDECL baseline serializes only cell PORO/PERMX/PERMY/PERMZ");
                break;
            }
        }
        return report;
    }

    if (!has_supported_generic_topology(
            document)) {
        report.note_unsupported(
            "generic.linear_topology_required",
            "Gmsh/VTU canonical writers require a supported linear 2D/3D topology with materialized faces");
        return report;
    }

    if (target_format ==
        MeshExchangeFormat::gmsh_4_1_ascii) {
        if (!document.fields().empty()) {
            report.note_lossy(
                "gmsh.fields_not_serialized",
                "current Gmsh writer bridge does not emit NodeData/ElementData; canonical fields are omitted");
        }
        if (document.logical_corner_point()
                .has_value()) {
            report.note_lossy(
                "gmsh.logical_corner_point_not_serialized",
                "Gmsh output contains the active generic mesh but not GRDECL logical corner-point/ACTNUM semantics");
        }
        for (const auto& group :
             document.groups()) {
            if (group.location !=
                    EntityKind::face &&
                group.location !=
                    EntityKind::cell) {
                report.note_lossy(
                    "gmsh.group_location_not_serialized",
                    "current Gmsh bridge serializes only face/cell physical groups");
                break;
            }
        }
        if (gmsh_ids_require_remap(
                document.topology())) {
            report.note_lossy(
                "gmsh.entity_ids_remapped",
                "Gmsh requires positive node tags and globally unique face/cell element tags; conflicting canonical IDs are deterministically remapped");
        }
        return report;
    }

    if (target_format ==
        MeshExchangeFormat::vtu_ascii) {
        if (!document.groups().empty()) {
            report.note_lossy(
                "vtu.groups_not_serialized",
                "current VTU bridge does not encode canonical named/physical groups");
        }
        if (has_tagged_faces(document)) {
            report.note_lossy(
                "vtu.face_tags_not_serialized",
                "current VTU bridge writes volume/surface topology fields but not boundary face PhysicalTag metadata");
        }
        if (document.logical_corner_point()
                .has_value()) {
            report.note_lossy(
                "vtu.logical_corner_point_not_serialized",
                "VTU output contains the active generic mesh but not GRDECL logical corner-point/ACTNUM semantics");
        }
        for (const auto& field :
             document.fields().fields()) {
            if (field.location() !=
                    EntityKind::vertex &&
                field.location() !=
                    EntityKind::cell) {
                report.note_lossy(
                    "vtu.non_point_cell_field_not_serialized",
                    "current VTU bridge serializes only vertex PointData and cell CellData fields");
                break;
            }
        }
        return report;
    }

    report.note_unsupported(
        "unknown_target_format",
        "unsupported canonical target format");
    return report;
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_MESH_EXCHANGE_HPP
