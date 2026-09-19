#ifndef MPMC_MESH_MESH_EXCHANGE_IO_HPP
#define MPMC_MESH_MESH_EXCHANGE_IO_HPP

#include <mpmc/mesh/active_corner_point.hpp>
#include <mpmc/mesh/gmsh_4_1.hpp>
#include <mpmc/mesh/gmsh_4_1_3d.hpp>
#include <mpmc/mesh/grdecl.hpp>
#include <mpmc/mesh/mesh_exchange.hpp>
#include <mpmc/mesh/vtu.hpp>
#include <mpmc/mesh/vtu_3d.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::mesh {

struct GrdeclExportOptions {
    /// Numeric output value multiplied by this scale yields SI metres.
    double coordinate_scale_to_m;
    /// Numeric output value multiplied by this scale yields SI square metres.
    double permeability_scale_to_m2;
};

struct MeshTextExportResult {
    std::optional<std::string> content;
    ConversionReport report;

    [[nodiscard]] bool exported()
        const noexcept {
        return content.has_value();
    }
};

namespace mesh_exchange_io_detail {

[[nodiscard]] inline EntityKind
gmsh_location(
    int mesh_dimension,
    int group_dimension) {
    if (group_dimension == 0) {
        return EntityKind::vertex;
    }
    if (mesh_dimension == 2) {
        if (group_dimension == 1) {
            return EntityKind::face;
        }
        if (group_dimension == 2) {
            return EntityKind::cell;
        }
    } else if (mesh_dimension == 3) {
        if (group_dimension == 1) {
            return EntityKind::edge;
        }
        if (group_dimension == 2) {
            return EntityKind::face;
        }
        if (group_dimension == 3) {
            return EntityKind::cell;
        }
    }
    throw std::invalid_argument(
        "mpmc::mesh::make_mesh_exchange_document: invalid Gmsh physical-group dimension");
}

using GroupKey =
    std::pair<EntityKind, std::uint32_t>;

[[nodiscard]] inline std::vector<MeshExchangeGroup>
gmsh_groups(
    int mesh_dimension,
    const Topology& topology,
    const FaceBoundarySnapshot& boundary,
    std::span<const GmshPhysicalName>
        physical_names,
    std::span<const GmshCellPhysicalGroups>
        cell_groups) {
    std::map<GroupKey, MeshExchangeGroup>
        groups;

    const auto ensure =
        [&](EntityKind location,
            std::uint32_t tag)
            -> MeshExchangeGroup& {
            const GroupKey key{
                location,
                tag};
            auto [found, inserted] =
                groups.emplace(
                    key,
                    MeshExchangeGroup{
                        location,
                        tag,
                        {},
                        {}});
            (void)inserted;
            return found->second;
        };

    for (const auto& name :
         physical_names) {
        auto& group =
            ensure(
                gmsh_location(
                    mesh_dimension,
                    name.dimension),
                name.tag);
        group.name = name.name;
    }

    const auto face_ids =
        topology.global_ids(
            EntityKind::face);
    for (std::size_t face = 0U;
         face < boundary.face_count();
         ++face) {
        const auto tag =
            boundary.physical_tag(
                LocalIndex{
                    static_cast<
                        LocalIndex::value_type>(
                        face)});
        if (!tag.is_tagged()) {
            continue;
        }
        ensure(
            EntityKind::face,
            tag.value())
            .members.push_back(
                face_ids[face]);
    }

    for (const auto& membership :
         cell_groups) {
        for (const auto tag :
             membership.physical_tags) {
            ensure(
                EntityKind::cell,
                tag)
                .members.push_back(
                    membership.cell_global_id);
        }
    }

    std::vector<MeshExchangeGroup> result;
    result.reserve(groups.size());
    for (auto& [key, group] : groups) {
        std::sort(
            group.members.begin(),
            group.members.end(),
            [](GlobalEntityId left,
               GlobalEntityId right) {
                return left.value() <
                    right.value();
            });
        group.members.erase(
            std::unique(
                group.members.begin(),
                group.members.end()),
            group.members.end());
        result.push_back(
            std::move(group));
        (void)key;
    }
    return result;
}

[[nodiscard]] inline std::vector<Coordinate3D>
coordinates_3d(
    const Geometry2D& geometry) {
    std::vector<Coordinate3D> result;
    result.reserve(
        geometry.vertex_count());
    for (const auto coordinate :
         geometry.vertex_coordinates_m()) {
        result.push_back(
            Coordinate3D{
                coordinate.x_m,
                coordinate.y_m,
                0.0});
    }
    return result;
}

[[nodiscard]] inline std::vector<DenseFieldSnapshot>
combine_fields(
    std::span<const DenseFieldSnapshot>
        first,
    std::span<const DenseFieldSnapshot>
        second) {
    std::vector<DenseFieldSnapshot> result;
    result.reserve(
        first.size() + second.size());
    for (const auto& field : first) {
        result.push_back(field);
    }
    for (const auto& field : second) {
        result.push_back(field);
    }
    return result;
}

[[nodiscard]] inline const DenseFieldSnapshot&
required_grdecl_field(
    const GrdeclImportResult& raw,
    std::string_view id,
    std::string_view unit) {
    const auto found =
        std::find_if(
            raw.cell_fields.begin(),
            raw.cell_fields.end(),
            [&](const auto& field) {
                return field.metadata().id ==
                    id;
            });
    if (found == raw.cell_fields.end() ||
        found->location() !=
            EntityKind::cell ||
        found->component_count() != 1U ||
        found->entity_count() !=
            raw.cell_count() ||
        found->metadata().unit != unit) {
        throw std::invalid_argument(
            "mpmc::mesh::make_mesh_exchange_document: required GRDECL field contract is absent");
    }
    return *found;
}

[[nodiscard]] inline std::vector<double>
scalar_values(
    const DenseFieldSnapshot& field) {
    return std::vector<double>{
        field.values().begin(),
        field.values().end()};
}

struct GmshTargetIds {
    std::vector<GlobalEntityId> vertices;
    std::vector<GlobalEntityId> faces;
    std::vector<GlobalEntityId> cells;
};

[[nodiscard]] inline std::vector<GlobalEntityId>
sequential_ids(
    std::size_t count,
    GlobalEntityId::value_type first) {
    if (count != 0U &&
        first >
            std::numeric_limits<
                GlobalEntityId::value_type>::max() -
                static_cast<
                    GlobalEntityId::value_type>(
                    count - 1U)) {
        throw std::length_error(
            "mpmc::mesh::canonical writer: stable ID allocation overflow");
    }
    std::vector<GlobalEntityId> ids;
    ids.reserve(count);
    for (std::size_t index = 0U;
         index < count;
         ++index) {
        ids.emplace_back(
            first +
            static_cast<
                GlobalEntityId::value_type>(
                index));
    }
    return ids;
}

[[nodiscard]] inline GmshTargetIds
gmsh_target_ids(
    const Topology& topology) {
    GmshTargetIds result{
        std::vector<GlobalEntityId>{
            topology.global_ids(
                EntityKind::vertex).begin(),
            topology.global_ids(
                EntityKind::vertex).end()},
        std::vector<GlobalEntityId>{
            topology.global_ids(
                EntityKind::face).begin(),
            topology.global_ids(
                EntityKind::face).end()},
        std::vector<GlobalEntityId>{
            topology.global_ids(
                EntityKind::cell).begin(),
            topology.global_ids(
                EntityKind::cell).end()}};

    if (std::any_of(
            result.vertices.begin(),
            result.vertices.end(),
            [](GlobalEntityId id) {
                return id.value() == 0U;
            })) {
        result.vertices =
            sequential_ids(
                result.vertices.size(),
                1U);
    }

    if (std::any_of(
            result.cells.begin(),
            result.cells.end(),
            [](GlobalEntityId id) {
                return id.value() == 0U;
            })) {
        result.cells =
            sequential_ids(
                result.cells.size(),
                1U);
    }

    std::vector<
        GlobalEntityId::value_type>
        element_ids;
    element_ids.reserve(
        result.faces.size() +
        result.cells.size());
    bool remap_faces =
        std::any_of(
            result.faces.begin(),
            result.faces.end(),
            [](GlobalEntityId id) {
                return id.value() == 0U;
            });
    for (const auto id : result.faces) {
        element_ids.push_back(id.value());
    }
    for (const auto id : result.cells) {
        element_ids.push_back(id.value());
    }
    std::sort(
        element_ids.begin(),
        element_ids.end());
    remap_faces =
        remap_faces ||
        std::adjacent_find(
            element_ids.begin(),
            element_ids.end()) !=
            element_ids.end();

    if (remap_faces) {
        GlobalEntityId::value_type maximum =
            0U;
        for (const auto id :
             result.cells) {
            maximum =
                std::max(
                    maximum,
                    id.value());
        }
        if (maximum ==
            std::numeric_limits<
                GlobalEntityId::value_type>::max()) {
            throw std::length_error(
                "mpmc::mesh::canonical Gmsh writer: cannot allocate face element tags");
        }
        result.faces =
            sequential_ids(
                result.faces.size(),
                maximum + 1U);
    }
    return result;
}

[[nodiscard]] inline std::vector<CsrAdjacency>
copy_relations(
    const Topology& topology) {
    constexpr std::array<EntityKind, 4>
        kinds{
            EntityKind::vertex,
            EntityKind::edge,
            EntityKind::face,
            EntityKind::cell};
    std::vector<CsrAdjacency> relations;
    relations.reserve(
        topology.relation_count());
    for (const auto source : kinds) {
        for (const auto target : kinds) {
            if (topology.has_relation(
                    source,
                    target)) {
                relations.push_back(
                    topology.relation(
                        source,
                        target));
            }
        }
    }
    return relations;
}

[[nodiscard]] inline Topology
topology_with_gmsh_ids(
    const Topology& source,
    const GmshTargetIds& ids) {
    Topology::EntityIds target_ids;
    target_ids.vertices = ids.vertices;
    target_ids.edges =
        std::vector<GlobalEntityId>{
            source.global_ids(
                EntityKind::edge).begin(),
            source.global_ids(
                EntityKind::edge).end()};
    target_ids.faces = ids.faces;
    target_ids.cells = ids.cells;
    return Topology{
        std::move(target_ids),
        copy_relations(source)};
}

[[nodiscard]] inline std::size_t
local_from_global(
    const Topology& topology,
    EntityKind location,
    GlobalEntityId id) {
    const auto ids =
        topology.global_ids(location);
    const auto found =
        std::find(
            ids.begin(),
            ids.end(),
            id);
    if (found == ids.end()) {
        throw std::invalid_argument(
            "mpmc::mesh::canonical writer: group member is absent from topology");
    }
    return static_cast<std::size_t>(
        std::distance(
            ids.begin(),
            found));
}

[[nodiscard]] inline FaceBoundarySnapshot
canonical_face_boundary(
    const MeshExchangeDocument& document,
    const Topology& target_topology) {
    const std::size_t face_count =
        target_topology.entity_count(
            EntityKind::face);
    std::vector<PhysicalTag> tags(
        face_count,
        PhysicalTag{0U});

    if (document.face_boundary()
            .has_value()) {
        const auto source_tags =
            document.face_boundary()->
                physical_tags();
        if (source_tags.size() !=
            face_count) {
            throw std::invalid_argument(
                "mpmc::mesh::canonical writer: face boundary count mismatch");
        }
        tags.assign(
            source_tags.begin(),
            source_tags.end());
    }

    for (const auto& group :
         document.groups()) {
        if (group.location !=
            EntityKind::face) {
            continue;
        }
        for (const auto member :
             group.members) {
            const std::size_t local =
                local_from_global(
                    document.topology(),
                    EntityKind::face,
                    member);
            if (tags[local].is_tagged() &&
                tags[local].value() !=
                    group.tag) {
                throw std::invalid_argument(
                    "mpmc::mesh::canonical Gmsh writer: one face cannot belong to multiple physical groups");
            }
            tags[local] =
                PhysicalTag{
                    group.tag};
        }
    }

    return make_face_boundary_snapshot(
        target_topology,
        tags);
}

[[nodiscard]] inline std::vector<GmshPhysicalName>
canonical_physical_names(
    const MeshExchangeDocument& document) {
    std::vector<GmshPhysicalName> names;
    for (const auto& group :
         document.groups()) {
        if (group.name.empty() ||
            (group.location !=
                 EntityKind::face &&
             group.location !=
                 EntityKind::cell)) {
            continue;
        }
        int dimension = 0;
        if (document.dimension() == 2) {
            dimension =
                group.location ==
                        EntityKind::face
                    ? 1
                    : 2;
        } else {
            dimension =
                group.location ==
                        EntityKind::face
                    ? 2
                    : 3;
        }
        names.push_back(
            GmshPhysicalName{
                dimension,
                group.tag,
                group.name});
    }
    return names;
}

[[nodiscard]] inline std::vector<
    GmshCellPhysicalGroups>
canonical_cell_groups(
    const MeshExchangeDocument& document,
    const Topology& target_topology) {
    const std::size_t cell_count =
        target_topology.entity_count(
            EntityKind::cell);
    std::vector<std::vector<std::uint32_t>>
        tags(cell_count);

    for (const auto& group :
         document.groups()) {
        if (group.location !=
            EntityKind::cell) {
            continue;
        }
        for (const auto member :
             group.members) {
            const std::size_t local =
                local_from_global(
                    document.topology(),
                    EntityKind::cell,
                    member);
            tags[local].push_back(
                group.tag);
        }
    }

    std::vector<GmshCellPhysicalGroups>
        result;
    result.reserve(cell_count);
    const auto target_ids =
        target_topology.global_ids(
            EntityKind::cell);
    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        auto& cell_tags = tags[cell];
        std::sort(
            cell_tags.begin(),
            cell_tags.end());
        cell_tags.erase(
            std::unique(
                cell_tags.begin(),
                cell_tags.end()),
            cell_tags.end());
        result.push_back(
            GmshCellPhysicalGroups{
                target_ids[cell],
                std::move(cell_tags)});
    }
    return result;
}

[[nodiscard]] inline Geometry2D
canonical_geometry_2d(
    const MeshExchangeDocument& document,
    const Topology& topology) {
    if (document.dimension() != 2) {
        throw std::invalid_argument(
            "mpmc::mesh::canonical 2D writer: dimension mismatch");
    }

    std::vector<Coordinate2D> coordinates;
    coordinates.reserve(
        document.vertex_coordinates_m()
            .size());
    for (const auto coordinate :
         document.vertex_coordinates_m()) {
        coordinates.push_back(
            Coordinate2D{
                coordinate.x_m,
                coordinate.y_m});
    }

    const auto& cell_vertices =
        topology.relation(
            EntityKind::cell,
            EntityKind::vertex);
    std::vector<Coordinate2D>
        cell_centroids;
    std::vector<double> cell_areas;
    cell_centroids.reserve(
        topology.entity_count(
            EntityKind::cell));
    cell_areas.reserve(
        topology.entity_count(
            EntityKind::cell));
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
        std::vector<Coordinate2D>
            polygon;
        polygon.reserve(vertices.size());
        for (const auto vertex :
             vertices) {
            polygon.push_back(
                coordinates[
                    static_cast<std::size_t>(
                        vertex.value())]);
        }
        const auto metric =
            vtu_detail::polygon_metric(
                polygon);
        cell_centroids.push_back(
            metric.centroid);
        cell_areas.push_back(
            metric.area);
    }

    const auto& face_vertices =
        topology.relation(
            EntityKind::face,
            EntityKind::vertex);
    const auto& face_cells =
        topology.relation(
            EntityKind::face,
            EntityKind::cell);
    std::vector<Coordinate2D>
        face_centroids;
    std::vector<double> face_lengths;
    std::vector<LocalIndex> face_owners;
    std::vector<UnitVector2D>
        face_normals;
    const std::size_t face_count =
        topology.entity_count(
            EntityKind::face);
    face_centroids.reserve(face_count);
    face_lengths.reserve(face_count);
    face_owners.reserve(face_count);
    face_normals.reserve(face_count);

    for (std::size_t face = 0U;
         face < face_count;
         ++face) {
        const LocalIndex local{
            static_cast<
                LocalIndex::value_type>(
                face)};
        const auto vertices =
            face_vertices.adjacent(local);
        if (vertices.size() != 2U) {
            throw std::invalid_argument(
                "mpmc::mesh::canonical 2D writer: every face must contain two vertices");
        }
        const auto a =
            coordinates[
                static_cast<std::size_t>(
                    vertices[0].value())];
        const auto b =
            coordinates[
                static_cast<std::size_t>(
                    vertices[1].value())];
        const double dx =
            b.x_m - a.x_m;
        const double dy =
            b.y_m - a.y_m;
        const double length =
            std::hypot(dx, dy);
        if (!std::isfinite(length) ||
            length <= 0.0) {
            throw std::invalid_argument(
                "mpmc::mesh::canonical 2D writer: invalid face length");
        }
        const Coordinate2D centroid{
            0.5 * (a.x_m + b.x_m),
            0.5 * (a.y_m + b.y_m)};
        const auto adjacent =
            face_cells.adjacent(local);
        if (adjacent.empty() ||
            adjacent.size() > 2U) {
            throw std::invalid_argument(
                "mpmc::mesh::canonical 2D writer: face support must be one or two cells");
        }
        const LocalIndex owner =
            adjacent.front();
        const auto owner_centroid =
            cell_centroids[
                static_cast<std::size_t>(
                    owner.value())];
        UnitVector2D normal{
            -dy / length,
            dx / length};
        const double dot =
            normal.x *
                (centroid.x_m -
                 owner_centroid.x_m) +
            normal.y *
                (centroid.y_m -
                 owner_centroid.y_m);
        const double tolerance =
            256.0 *
            std::numeric_limits<double>::
                epsilon() *
            std::max(1.0, length);
        if (!std::isfinite(dot) ||
            std::abs(dot) <= tolerance) {
            throw std::invalid_argument(
                "mpmc::mesh::canonical 2D writer: cannot orient face normal");
        }
        if (dot < 0.0) {
            normal.x = -normal.x;
            normal.y = -normal.y;
        }
        face_centroids.push_back(
            centroid);
        face_lengths.push_back(
            length);
        face_owners.push_back(
            owner);
        face_normals.push_back(
            normal);
    }

    return Geometry2D{
        std::move(coordinates),
        std::move(cell_centroids),
        std::move(cell_areas),
        std::move(face_centroids),
        std::move(face_lengths),
        std::move(face_owners),
        std::move(face_normals)};
}

[[nodiscard]] inline LinearCellType3D
linear_cell_type(
    std::size_t vertex_count) {
    if (vertex_count == 4U) {
        return LinearCellType3D::tetrahedron;
    }
    if (vertex_count == 8U) {
        return LinearCellType3D::hexahedron;
    }
    if (vertex_count == 6U) {
        return LinearCellType3D::wedge;
    }
    if (vertex_count == 5U) {
        return LinearCellType3D::pyramid;
    }
    throw std::invalid_argument(
        "mpmc::mesh::canonical 3D writer: unsupported linear cell width");
}

[[nodiscard]] inline LinearMesh3D
canonical_linear_mesh_3d(
    const MeshExchangeDocument& document,
    std::span<const GlobalEntityId>
        vertex_ids,
    std::span<const GlobalEntityId>
        face_ids,
    std::span<const GlobalEntityId>
        cell_ids,
    const FaceBoundarySnapshot& boundary) {
    if (document.dimension() != 3) {
        throw std::invalid_argument(
            "mpmc::mesh::canonical 3D writer: dimension mismatch");
    }
    const auto& topology =
        document.topology();
    const auto& cell_vertices =
        topology.relation(
            EntityKind::cell,
            EntityKind::vertex);
    const auto& face_vertices =
        topology.relation(
            EntityKind::face,
            EntityKind::vertex);

    std::vector<LinearCell3D> cells;
    cells.reserve(cell_ids.size());
    for (std::size_t cell = 0U;
         cell < cell_ids.size();
         ++cell) {
        const auto vertices =
            cell_vertices.adjacent(
                LocalIndex{
                    static_cast<
                        LocalIndex::value_type>(
                        cell)});
        cells.push_back(
            LinearCell3D{
                cell_ids[cell],
                linear_cell_type(
                    vertices.size()),
                std::vector<LocalIndex>{
                    vertices.begin(),
                    vertices.end()}});
    }

    std::vector<LinearFaceAnnotation3D>
        annotations;
    annotations.reserve(face_ids.size());
    for (std::size_t face = 0U;
         face < face_ids.size();
         ++face) {
        const LocalIndex local{
            static_cast<
                LocalIndex::value_type>(
                face)};
        const auto vertices =
            face_vertices.adjacent(local);
        annotations.push_back(
            LinearFaceAnnotation3D{
                std::vector<LocalIndex>{
                    vertices.begin(),
                    vertices.end()},
                face_ids[face],
                boundary.physical_tag(
                    local)});
    }

    GlobalEntityId::value_type floor =
        0U;
    for (const auto id : face_ids) {
        floor =
            std::max(
                floor,
                id.value());
    }
    for (const auto id : cell_ids) {
        floor =
            std::max(
                floor,
                id.value());
    }

    return make_linear_mesh_3d(
        vertex_ids,
        document.vertex_coordinates_m(),
        cells,
        annotations,
        floor);
}

inline void split_vtu_fields(
    const MeshExchangeDocument& document,
    std::vector<DenseFieldSnapshot>&
        point_fields,
    std::vector<DenseFieldSnapshot>&
        cell_fields) {
    for (const auto& field :
         document.fields().fields()) {
        if (field.location() ==
            EntityKind::vertex) {
            point_fields.push_back(field);
        } else if (
            field.location() ==
            EntityKind::cell) {
            cell_fields.push_back(field);
        }
    }
}

inline void require_export_scale(
    double value,
    const char* message) {
    if (!std::isfinite(value) ||
        value <= 0.0) {
        throw std::invalid_argument(message);
    }
}

inline void write_double_record(
    std::ostringstream& output,
    std::string_view keyword,
    std::span<const double> values,
    double scale_to_si,
    std::size_t values_per_line) {
    output << keyword << '\n';
    for (std::size_t index = 0U;
         index < values.size();
         ++index) {
        const double external =
            values[index] / scale_to_si;
        if (!std::isfinite(external)) {
            throw std::invalid_argument(
                "mpmc::mesh::export_grdecl_ascii: scaled output is non-finite");
        }
        output << external;
        if ((index + 1U) %
                values_per_line ==
                0U ||
            index + 1U == values.size()) {
            output << '\n';
        } else {
            output << ' ';
        }
    }
    output << "/\n\n";
}

inline void write_actnum(
    std::ostringstream& output,
    std::span<const std::uint8_t> active) {
    output << "ACTNUM\n";
    for (std::size_t index = 0U;
         index < active.size();
         ++index) {
        output <<
            (active[index] ==
                     std::uint8_t{0U}
                 ? 0
                 : 1);
        if ((index + 1U) % 16U == 0U ||
            index + 1U == active.size()) {
            output << '\n';
        } else {
            output << ' ';
        }
    }
    output << "/\n\n";
}

} // namespace mesh_exchange_io_detail

[[nodiscard]] inline MeshExchangeDocument
make_mesh_exchange_document(
    const Gmsh41ImportResult& source) {
    return MeshExchangeDocument::create(
        MeshExchangeFormat::gmsh_4_1_ascii,
        2,
        source.topology,
        mesh_exchange_io_detail::
            coordinates_3d(source.geometry),
        source.face_boundary,
        {},
        mesh_exchange_io_detail::gmsh_groups(
            2,
            source.topology,
            source.face_boundary,
            source.physical_names,
            source.cell_physical_groups),
        std::nullopt);
}

[[nodiscard]] inline MeshExchangeDocument
make_mesh_exchange_document(
    const Gmsh41ImportResult3D& source) {
    return MeshExchangeDocument::create(
        MeshExchangeFormat::gmsh_4_1_ascii,
        3,
        source.topology,
        source.vertex_coordinates_m,
        source.face_boundary,
        {},
        mesh_exchange_io_detail::gmsh_groups(
            3,
            source.topology,
            source.face_boundary,
            source.physical_names,
            source.cell_physical_groups),
        std::nullopt);
}

[[nodiscard]] inline MeshExchangeDocument
make_mesh_exchange_document(
    const VtuImportResult& source) {
    return MeshExchangeDocument::create(
        MeshExchangeFormat::vtu_ascii,
        2,
        source.topology,
        mesh_exchange_io_detail::
            coordinates_3d(source.geometry),
        make_face_boundary_snapshot(
            source.topology),
        mesh_exchange_io_detail::
            combine_fields(
                source.point_fields,
                source.cell_fields),
        {},
        std::nullopt);
}

[[nodiscard]] inline MeshExchangeDocument
make_mesh_exchange_document(
    const VtuImportResult3D& source) {
    return MeshExchangeDocument::create(
        MeshExchangeFormat::vtu_ascii,
        3,
        source.topology,
        source.vertex_coordinates_m,
        source.face_boundary,
        mesh_exchange_io_detail::
            combine_fields(
                source.point_fields,
                source.cell_fields),
        {},
        std::nullopt);
}

[[nodiscard]] inline MeshExchangeDocument
make_mesh_exchange_document(
    const GrdeclImportResult& source) {
    const auto& poro =
        mesh_exchange_io_detail::
            required_grdecl_field(
                source,
                "PORO",
                "1");
    const auto& permx =
        mesh_exchange_io_detail::
            required_grdecl_field(
                source,
                "PERMX",
                "m2");
    const auto& permy =
        mesh_exchange_io_detail::
            required_grdecl_field(
                source,
                "PERMY",
                "m2");
    const auto& permz =
        mesh_exchange_io_detail::
            required_grdecl_field(
                source,
                "PERMZ",
                "m2");

    LogicalCornerPointGrid3D corner_point{
        source.dimensions,
        source.coord_m,
        source.zcorn_m,
        source.active,
        mesh_exchange_io_detail::
            scalar_values(poro),
        mesh_exchange_io_detail::
            scalar_values(permx),
        mesh_exchange_io_detail::
            scalar_values(permy),
        mesh_exchange_io_detail::
            scalar_values(permz),
        source.source_options
            .coordinate_scale_to_m,
        source.source_options
            .permeability_scale_to_m2};

    if (source.active_cell_count() != 0U) {
        const auto processed =
            process_active_corner_point_grid(
                source);
        return MeshExchangeDocument::create(
            MeshExchangeFormat::grdecl,
            3,
            processed.topology,
            processed.vertex_coordinates_m,
            make_face_boundary_snapshot(
                processed.topology),
            processed.cell_fields,
            {},
            std::move(corner_point));
    }

    return MeshExchangeDocument::create(
        MeshExchangeFormat::grdecl,
        3,
        source.topology,
        std::vector<Coordinate3D>{
            source.geometry
                .vertex_coordinates_m()
                .begin(),
            source.geometry
                .vertex_coordinates_m()
                .end()},
        std::nullopt,
        source.cell_fields,
        {},
        std::move(corner_point));
}

[[nodiscard]] inline MeshTextExportResult
export_gmsh_4_1_ascii(
    const MeshExchangeDocument& document) {
    using namespace mesh_exchange_io_detail;

    auto report =
        analyze_conversion(
            document,
            MeshExchangeFormat::gmsh_4_1_ascii);
    if (report.disposition() ==
        ConversionDisposition::unsupported) {
        return MeshTextExportResult{
            std::nullopt,
            std::move(report)};
    }

    const auto target_ids =
        gmsh_target_ids(
            document.topology());

    if (document.dimension() == 2) {
        auto topology =
            topology_with_gmsh_ids(
                document.topology(),
                target_ids);
        auto boundary =
            canonical_face_boundary(
                document,
                topology);
        auto geometry =
            canonical_geometry_2d(
                document,
                topology);
        auto physical_names =
            canonical_physical_names(
                document);
        auto cell_groups =
            canonical_cell_groups(
                document,
                topology);

        Gmsh41ImportResult target{
            std::move(topology),
            std::move(geometry),
            std::move(boundary),
            std::move(physical_names),
            std::move(cell_groups)};
        return MeshTextExportResult{
            export_gmsh_4_1_ascii(
                target),
            std::move(report)};
    }

    auto target_topology =
        topology_with_gmsh_ids(
            document.topology(),
            target_ids);
    auto target_boundary =
        canonical_face_boundary(
            document,
            target_topology);
    auto linear =
        canonical_linear_mesh_3d(
            document,
            target_ids.vertices,
            target_ids.faces,
            target_ids.cells,
            target_boundary);
    auto physical_names =
        canonical_physical_names(
            document);
    auto cell_groups =
        canonical_cell_groups(
            document,
            linear.topology);

    Gmsh41ImportResult3D target{
        std::move(linear.topology),
        std::move(
            linear.vertex_coordinates_m),
        std::move(
            linear.cell_volumes_m3),
        std::move(
            linear.face_geometry),
        std::move(
            linear.face_boundary),
        std::move(physical_names),
        std::move(cell_groups)};
    return MeshTextExportResult{
        export_gmsh_4_1_ascii_3d(
            target),
        std::move(report)};
}

[[nodiscard]] inline MeshTextExportResult
export_vtu_ascii(
    const MeshExchangeDocument& document) {
    using namespace mesh_exchange_io_detail;

    auto report =
        analyze_conversion(
            document,
            MeshExchangeFormat::vtu_ascii);
    if (report.disposition() ==
        ConversionDisposition::unsupported) {
        return MeshTextExportResult{
            std::nullopt,
            std::move(report)};
    }

    std::vector<DenseFieldSnapshot>
        point_fields;
    std::vector<DenseFieldSnapshot>
        cell_fields;
    split_vtu_fields(
        document,
        point_fields,
        cell_fields);

    if (document.dimension() == 2) {
        auto geometry =
            canonical_geometry_2d(
                document,
                document.topology());
        VtuImportResult target{
            document.topology(),
            std::move(geometry),
            std::move(point_fields),
            std::move(cell_fields)};
        return MeshTextExportResult{
            export_vtu_ascii(target),
            std::move(report)};
    }

    const auto ids =
        document.topology().global_ids(
            EntityKind::vertex);
    const auto face_ids =
        document.topology().global_ids(
            EntityKind::face);
    const auto cell_ids =
        document.topology().global_ids(
            EntityKind::cell);
    const auto boundary =
        canonical_face_boundary(
            document,
            document.topology());
    auto linear =
        canonical_linear_mesh_3d(
            document,
            ids,
            face_ids,
            cell_ids,
            boundary);

    VtuImportResult3D target{
        std::move(linear.topology),
        std::move(
            linear.vertex_coordinates_m),
        std::move(
            linear.cell_volumes_m3),
        std::move(
            linear.face_geometry),
        std::move(
            linear.face_boundary),
        std::move(point_fields),
        std::move(cell_fields)};
    return MeshTextExportResult{
        export_vtu_ascii_3d(target),
        std::move(report)};
}

[[nodiscard]] inline MeshTextExportResult
export_grdecl_ascii(
    const MeshExchangeDocument& document,
    GrdeclExportOptions options) {
    using namespace mesh_exchange_io_detail;

    require_export_scale(
        options.coordinate_scale_to_m,
        "mpmc::mesh::export_grdecl_ascii: coordinate_scale_to_m must be finite and positive");
    require_export_scale(
        options.permeability_scale_to_m2,
        "mpmc::mesh::export_grdecl_ascii: permeability_scale_to_m2 must be finite and positive");

    auto report =
        analyze_conversion(
            document,
            MeshExchangeFormat::grdecl);
    if (!report.lossless()) {
        return MeshTextExportResult{
            std::nullopt,
            std::move(report)};
    }

    const auto& data =
        *document.logical_corner_point();

    std::ostringstream output;
    output <<
        std::setprecision(
            std::numeric_limits<double>::max_digits10);
    output
        << "-- MPMC minimal GRDECL export\n"
        << "-- coordinate_scale_to_m "
        << options.coordinate_scale_to_m
        << "\n"
        << "-- permeability_scale_to_m2 "
        << options.permeability_scale_to_m2
        << "\n\n";

    output << "SPECGRID\n"
           << data.dimensions[0] << ' '
           << data.dimensions[1] << ' '
           << data.dimensions[2]
           << " 1 F /\n\n";

    write_double_record(
        output,
        "COORD",
        data.coord_m,
        options.coordinate_scale_to_m,
        6U);
    write_double_record(
        output,
        "ZCORN",
        data.zcorn_m,
        options.coordinate_scale_to_m,
        8U);
    write_actnum(
        output,
        data.active);
    write_double_record(
        output,
        "PORO",
        data.porosity,
        1.0,
        8U);
    write_double_record(
        output,
        "PERMX",
        data.permx_m2,
        options.permeability_scale_to_m2,
        8U);
    write_double_record(
        output,
        "PERMY",
        data.permy_m2,
        options.permeability_scale_to_m2,
        8U);
    write_double_record(
        output,
        "PERMZ",
        data.permz_m2,
        options.permeability_scale_to_m2,
        8U);

    return MeshTextExportResult{
        output.str(),
        std::move(report)};
}

[[nodiscard]] inline MeshTextExportResult
export_grdecl_ascii(
    const MeshExchangeDocument& document) {
    const auto report =
        analyze_conversion(
            document,
            MeshExchangeFormat::grdecl);
    if (!report.lossless()) {
        return MeshTextExportResult{
            std::nullopt,
            report};
    }
    const auto& data =
        *document.logical_corner_point();
    return export_grdecl_ascii(
        document,
        GrdeclExportOptions{
            data.source_coordinate_scale_to_m,
            data.source_permeability_scale_to_m2});
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_MESH_EXCHANGE_IO_HPP
