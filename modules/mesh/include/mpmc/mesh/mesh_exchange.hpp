#ifndef MPMC_MESH_MESH_EXCHANGE_HPP
#define MPMC_MESH_MESH_EXCHANGE_HPP

#include <mpmc/mesh/mesh_exchange_document.hpp>
#include <mpmc/mesh/grdecl_reconstruction.hpp>
#include <mpmc/mesh/grdecl.hpp>
#include <mpmc/mesh/active_corner_point.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mpmc::mesh {

namespace mesh_exchange_detail {

struct GrdeclIdentityMesh {
    Topology topology;
    std::vector<Coordinate3D> coordinates;
};

// Reuse the GRDECL geometry and active-grid implementation to predict the
// actual readback identities. Do not duplicate its vertex/face numbering rule.
[[nodiscard]] inline GrdeclIdentityMesh grdecl_identity_mesh(
    const LogicalCornerPointGrid3D& data) {
    const auto [nx, ny, nz] = data.dimensions;
    const auto count = data.active.size();
    const auto vertex_count = grdecl_detail::checked_multiply(
        count, 8U, "GRDECL identity corner count overflow");
    Topology::EntityIds ids;
    std::vector<Coordinate3D> coordinates;
    std::vector<double> volumes;
    std::vector<CsrAdjacency::Offset> offsets{0U};
    std::vector<LocalIndex> vertices;
    ids.cells.reserve(count);
    ids.vertices.reserve(vertex_count);
    coordinates.reserve(vertex_count);
    vertices.reserve(vertex_count);
    volumes.reserve(count);
    offsets.reserve(count + 1U);
    for (std::size_t k = 0U; k < nz; ++k) {
        for (std::size_t j = 0U; j < ny; ++j) {
            for (std::size_t i = 0U; i < nx; ++i) {
                const auto cell = i + nx * (j + ny * k);
                ids.cells.emplace_back(static_cast<std::uint64_t>(cell) + 1U);
                std::array<Coordinate3D, 8> corners{};
                for (std::size_t corner = 0U; corner < 8U; ++corner) {
                    const auto ci = corner % 2U;
                    const auto cj = (corner / 2U) % 2U;
                    const auto ck = corner / 4U;
                    const auto pillar = (i + ci + (nx + 1U) * (j + cj)) * 6U;
                    const auto z = (2U * i + ci) + 2U * nx *
                        ((2U * j + cj) + 2U * ny * (2U * k + ck));
                    corners[corner] = grdecl_detail::point_on_pillar(
                        grdecl_detail::Pillar{
                            {data.coord_m[pillar], data.coord_m[pillar + 1U], data.coord_m[pillar + 2U]},
                            {data.coord_m[pillar + 3U], data.coord_m[pillar + 4U], data.coord_m[pillar + 5U]}},
                        data.zcorn_m[z]);
                    const auto vertex = coordinates.size();
                    ids.vertices.emplace_back(static_cast<std::uint64_t>(vertex) + 1U);
                    vertices.push_back(grdecl_detail::local_index(vertex, "GRDECL identity vertex overflow"));
                    coordinates.push_back(corners[corner]);
                }
                volumes.push_back(grdecl_detail::check_cell_volume(corners, data.active[cell] != 0U).volume_m3);
                offsets.push_back(grdecl_detail::csr_offset(vertices.size(), "GRDECL identity CSR overflow"));
            }
        }
    }
    std::vector<CsrAdjacency> relations;
    relations.emplace_back(EntityKind::cell, EntityKind::vertex, vertex_count,
                          std::move(offsets), std::move(vertices));
    GrdeclImportResult raw{data.dimensions, Topology{std::move(ids), std::move(relations)},
        CornerPointGeometry3D{std::move(coordinates), std::move(volumes)}, data.active,
        {}, data.coord_m, data.zcorn_m, {1.0, 1.0}};
    if (raw.active_cell_count() == 0U) {
        const auto values = raw.geometry.vertex_coordinates_m();
        return {std::move(raw.topology), {values.begin(), values.end()}};
    }
    auto processed = process_active_corner_point_grid(raw);
    return {std::move(processed.topology), std::move(processed.vertex_coordinates_m)};
}

[[nodiscard]] inline bool grdecl_ids_change(
    const MeshExchangeDocument& document, const GrdeclIdentityMesh& target) {
    const auto coordinate_key = [](Coordinate3D coordinate) {
        return std::array<double, 3>{coordinate.x_m, coordinate.y_m, coordinate.z_m};
    };
    for (const auto kind : {EntityKind::vertex, EntityKind::edge, EntityKind::face, EntityKind::cell}) {
        const auto source_ids = document.topology().global_ids(kind);
        const auto target_ids = target.topology.global_ids(kind);
        if (source_ids.size() != target_ids.size()) return true;
        GroupEntityLookup lookup{target_ids};
        for (std::size_t local = 0U; local < source_ids.size(); ++local) {
            const auto found = lookup.find(source_ids[local]);
            if (!found) return true;
            if (kind == EntityKind::vertex) {
                if (coordinate_key(document.vertex_coordinates_m()[local]) !=
                    coordinate_key(target.coordinates[*found])) return true;
                continue;
            }
            if (!document.topology().has_relation(kind, EntityKind::vertex) ||
                !target.topology.has_relation(kind, EntityKind::vertex)) return true;
            std::vector<std::array<double, 3>> source_corners, target_corners;
            for (const auto vertex : document.topology().relation(kind, EntityKind::vertex).adjacent(
                     grdecl_detail::local_index(local, "source identity index overflow"))) {
                source_corners.push_back(coordinate_key(document.vertex_coordinates_m()[vertex.value()]));
            }
            for (const auto vertex : target.topology.relation(kind, EntityKind::vertex).adjacent(
                     grdecl_detail::local_index(*found, "target identity index overflow"))) {
                target_corners.push_back(coordinate_key(target.coordinates[vertex.value()]));
            }
            std::sort(source_corners.begin(), source_corners.end());
            std::sort(target_corners.begin(), target_corners.end());
            if (source_corners != target_corners) return true;
        }
    }
    return false;
}

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


// Both preflight and serialization use this preparation. A reconstructed
// logical grid is retained for the writer, so no second reconstruction occurs.
struct PreparedGrdeclConversion {
    ConversionReport report{MeshExchangeFormat::grdecl};
    std::optional<LogicalCornerPointGrid3D> reconstructed;
};

[[nodiscard]] inline PreparedGrdeclConversion prepare_grdecl_conversion(
    const MeshExchangeDocument& document) {
    PreparedGrdeclConversion prepared;
    auto& report = prepared.report;
    if (document.dimension() != 3) {
        report.note_unsupported("grdecl.requires_3d", "GRDECL baseline requires a 3D document");
        return prepared;
    }
    GrdeclRepresentabilityReport reconstruction_report;
    if (!document.logical_corner_point()) {
        auto reconstruction = reconstruct_structured_logical_grid_3d(document);
        reconstruction_report = std::move(reconstruction.report);
        if (!reconstruction_report.representable() || !reconstruction.logical_grid) {
            for (const auto& issue : reconstruction_report.issues()) {
                report.note_unsupported(issue.code, issue.message);
            }
            if (report.issues().empty()) {
                report.note_unsupported("grdecl_reconstruction.unsupported",
                    "generic canonical mesh is outside the structured GRDECL reconstruction baseline");
            }
            return prepared;
        }
        prepared.reconstructed.emplace(std::move(*reconstruction.logical_grid));
    }
    if (!document.groups().empty()) {
        report.note_lossy("grdecl.groups_not_representable",
            "minimal GRDECL baseline does not serialize generic named/physical groups");
    }
    if (has_tagged_faces(document)) {
        report.note_lossy("grdecl.face_tags_not_representable",
            "minimal GRDECL baseline does not serialize generic face physical tags");
    }
    for (const auto& field : document.fields().fields()) {
        if (field.location() != EntityKind::cell ||
            (field.metadata().id != "PORO" && field.metadata().id != "PERMX" &&
             field.metadata().id != "PERMY" && field.metadata().id != "PERMZ")) {
            report.note_lossy("grdecl.field_not_representable",
                "minimal GRDECL baseline serializes only cell PORO/PERMX/PERMY/PERMZ");
            break;
        }
    }
    for (const auto& issue : reconstruction_report.issues()) {
        report.note_lossy(issue.code, issue.message);
    }
    const auto& logical = prepared.reconstructed ? *prepared.reconstructed : *document.logical_corner_point();
    try {
        const auto target = grdecl_identity_mesh(logical);
        if (grdecl_ids_change(document, target)) {
            report.note_lossy("grdecl.entity_ids_remapped",
                "GRDECL does not store arbitrary entity IDs; readback changes vertex, face or cell identity bindings");
        }
    } catch (const std::invalid_argument& error) {
        report.note_unsupported("grdecl.computational_topology_unsupported", error.what());
    }
    return prepared;
}

} // namespace mesh_exchange_detail

[[nodiscard]] inline ConversionReport
analyze_conversion(
    const MeshExchangeDocument& document,
    MeshExchangeFormat target_format) {
    using namespace mesh_exchange_detail;

    ConversionReport report{target_format};

    if (target_format == MeshExchangeFormat::grdecl) {
        return prepare_grdecl_conversion(document).report;
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
        const auto& topology = document.topology();
        std::vector<std::uint32_t> face_tags(topology.entity_count(EntityKind::face), 0U);
        if (document.face_boundary()) {
            const auto tags = document.face_boundary()->physical_tags();
            for (std::size_t i = 0; i < tags.size(); ++i) face_tags[i] = tags[i].value();
        }
        GroupEntityLookup face_lookup{topology.global_ids(EntityKind::face)};
        for (const auto& group : document.groups()) {
            if (group.location != EntityKind::face) continue;
            for (const auto member : group.members) {
                const auto face = *face_lookup.find(member);
                const auto support = topology.relation(EntityKind::face, EntityKind::cell).adjacent(
                    LocalIndex{static_cast<LocalIndex::value_type>(face)});
                if (support.size() != 1U || (face_tags[face] != 0U && face_tags[face] != group.tag)) {
                    report.note_unsupported("gmsh.face_groups_not_representable",
                        "current Gmsh boundary bridge requires one physical tag per boundary face; overlapping, internal or conflicting face groups cannot be serialized");
                    return report;
                }
                face_tags[face] = group.tag;
            }
        }
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
