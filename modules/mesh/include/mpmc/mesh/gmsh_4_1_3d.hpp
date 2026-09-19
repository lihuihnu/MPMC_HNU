#ifndef MPMC_MESH_GMSH_4_1_3D_HPP
#define MPMC_MESH_GMSH_4_1_3D_HPP

#include <mpmc/mesh/gmsh_4_1.hpp>
#include <mpmc/mesh/linear_cell_mesh_3d.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace mpmc::mesh {

struct Gmsh41ImportResult3D {
    Topology topology;
    std::vector<Coordinate3D> vertex_coordinates_m;
    std::vector<double> cell_volumes_m3;
    FaceGeometry3D face_geometry;
    FaceBoundarySnapshot face_boundary;
    std::vector<GmshPhysicalName> physical_names;
    std::vector<GmshCellPhysicalGroups> cell_physical_groups;
};

namespace gmsh41_3d_detail {

struct ParsedDocument {
    bool entities_seen{false};
    std::map<gmsh41_detail::EntityKey, std::string>
        physical_names;
    std::map<gmsh41_detail::EntityKey, gmsh41_detail::EntityInfo>
        entities;
    std::vector<gmsh41_detail::NodeRecord> nodes;
    std::vector<gmsh41_detail::ElementRecord> elements;
};

[[nodiscard]] inline ParsedDocument parse_document(
    std::string_view content) {
    using namespace gmsh41_detail;

    std::istringstream input{std::string(content)};
    bool mesh_format_seen = false;
    bool physical_names_seen = false;
    bool entities_seen = false;
    bool nodes_seen = false;
    bool elements_seen = false;

    ParsedDocument parsed;
    std::string section;
    while (input >> section) {
        if (section == "$MeshFormat") {
            if (mesh_format_seen) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii_3d: duplicate $MeshFormat section");
            }
            double version = 0.0;
            int file_type = -1;
            int data_size = 0;
            read_value(
                input, version,
                "mpmc::mesh::import_gmsh_4_1_ascii_3d: invalid MSH version");
            read_value(
                input, file_type,
                "mpmc::mesh::import_gmsh_4_1_ascii_3d: invalid MSH file type");
            read_value(
                input, data_size,
                "mpmc::mesh::import_gmsh_4_1_ascii_3d: invalid MSH data size");
            if (std::abs(version - 4.1) >
                8.0 *
                    std::numeric_limits<double>::epsilon()) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii_3d: only MSH version 4.1 is supported");
            }
            if (file_type != 0) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii_3d: binary MSH is unsupported; ASCII is required");
            }
            if (data_size != 8) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii_3d: only 8-byte MSH data-size files are supported");
            }
            require_token(
                input, "$EndMeshFormat",
                "mpmc::mesh::import_gmsh_4_1_ascii_3d: missing $EndMeshFormat");
            mesh_format_seen = true;
        } else if (section == "$PhysicalNames") {
            if (physical_names_seen) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii_3d: duplicate $PhysicalNames section");
            }
            parse_physical_names(
                input, parsed.physical_names);
            physical_names_seen = true;
        } else if (section == "$Entities") {
            if (entities_seen) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii_3d: duplicate $Entities section");
            }
            parse_entities(
                input, parsed.entities);
            entities_seen = true;
        } else if (section == "$Nodes") {
            if (nodes_seen) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii_3d: duplicate $Nodes section");
            }
            parse_nodes(
                input, parsed.nodes);
            nodes_seen = true;
        } else if (section == "$Elements") {
            if (elements_seen) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii_3d: duplicate $Elements section");
            }
            parse_elements(
                input, parsed.elements);
            elements_seen = true;
        } else {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii_3d: unsupported or unexpected MSH section");
        }
    }

    if (!mesh_format_seen ||
        !nodes_seen ||
        !elements_seen) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii_3d: $MeshFormat, $Nodes and $Elements are required");
    }
    if (physical_names_seen &&
        !entities_seen) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii_3d: $Entities is required when Physical Groups are declared");
    }
    parsed.entities_seen = entities_seen;
    return parsed;
}

[[nodiscard]] inline LocalIndex lookup_node(
    const std::map<std::uint64_t, LocalIndex>& local_by_tag,
    std::uint64_t tag) {
    const auto found =
        local_by_tag.find(tag);
    if (found == local_by_tag.end()) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii_3d: element references missing node");
    }
    return found->second;
}

[[nodiscard]] inline std::array<double, 6> bounds(
    std::span<const LocalIndex> vertices,
    std::span<const Coordinate3D> coordinates) {
    std::array<double, 6> result{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity()};
    for (const auto vertex : vertices) {
        const auto point =
            coordinates[
                static_cast<std::size_t>(
                    vertex.value())];
        result[0] = std::min(result[0], point.x_m);
        result[1] = std::min(result[1], point.y_m);
        result[2] = std::min(result[2], point.z_m);
        result[3] = std::max(result[3], point.x_m);
        result[4] = std::max(result[4], point.y_m);
        result[5] = std::max(result[5], point.z_m);
    }
    return result;
}

inline void require_export_name(
    const GmshPhysicalName& name) {
    if ((name.dimension != 2 &&
         name.dimension != 3) ||
        name.tag == 0U ||
        name.name.empty() ||
        name.name.find('\0') !=
            std::string::npos ||
        name.name.find('\n') !=
            std::string::npos ||
        name.name.find('\r') !=
            std::string::npos) {
        throw std::invalid_argument(
            "mpmc::mesh::export_gmsh_4_1_ascii_3d: physical names must describe represented 2D/3D groups with nonempty single-line text");
    }
}

} // namespace gmsh41_3d_detail

[[nodiscard]] inline Gmsh41ImportResult3D
import_gmsh_4_1_ascii_3d(
    std::string_view content,
    double coordinate_scale_to_m) {
    using namespace gmsh41_detail;
    using namespace gmsh41_3d_detail;

    if (!std::isfinite(
            coordinate_scale_to_m) ||
        coordinate_scale_to_m <= 0.0) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii_3d: coordinate_scale_to_m must be finite and positive");
    }

    auto parsed =
        parse_document(content);
    std::sort(
        parsed.nodes.begin(),
        parsed.nodes.end(),
        [](const NodeRecord& left,
           const NodeRecord& right) {
            return left.tag < right.tag;
        });

    std::vector<GlobalEntityId>
        vertex_ids;
    std::vector<Coordinate3D>
        coordinates;
    std::map<std::uint64_t, LocalIndex>
        node_local_by_tag;
    vertex_ids.reserve(parsed.nodes.size());
    coordinates.reserve(parsed.nodes.size());

    for (std::size_t index = 0U;
         index < parsed.nodes.size();
         ++index) {
        const auto& node =
            parsed.nodes[index];
        const Coordinate3D coordinate{
            node.x * coordinate_scale_to_m,
            node.y * coordinate_scale_to_m,
            node.z * coordinate_scale_to_m};
        if (!std::isfinite(coordinate.x_m) ||
            !std::isfinite(coordinate.y_m) ||
            !std::isfinite(coordinate.z_m)) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii_3d: scaled node coordinate is non-finite");
        }
        vertex_ids.emplace_back(node.tag);
        coordinates.push_back(coordinate);
        node_local_by_tag.emplace(
            node.tag,
            gmsh41_detail::checked_local(
                index,
                "mpmc::mesh::import_gmsh_4_1_ascii_3d: node local index overflow"));
    }

    std::vector<ElementRecord>
        volume_elements;
    std::vector<ElementRecord>
        surface_elements;
    std::uint64_t maximum_element_tag = 0U;
    for (const auto& element :
         parsed.elements) {
        maximum_element_tag =
            std::max(
                maximum_element_tag,
                element.tag);
        if (parsed.entities_seen &&
            find_entity(
                parsed.entities,
                element.dimension,
                element.entity_tag) ==
                nullptr) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii_3d: element references entity absent from $Entities");
        }
        if (element.dimension == 3) {
            volume_elements.push_back(element);
        } else if (element.dimension == 2) {
            surface_elements.push_back(element);
        }
    }
    if (volume_elements.empty()) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii_3d: no supported linear 3D volume cells found");
    }

    std::sort(
        volume_elements.begin(),
        volume_elements.end(),
        [](const ElementRecord& left,
           const ElementRecord& right) {
            return left.tag < right.tag;
        });
    std::sort(
        surface_elements.begin(),
        surface_elements.end(),
        [](const ElementRecord& left,
           const ElementRecord& right) {
            return left.tag < right.tag;
        });

    std::vector<LinearCell3D> cells;
    std::vector<GmshCellPhysicalGroups>
        cell_groups;
    cells.reserve(volume_elements.size());
    cell_groups.reserve(volume_elements.size());

    std::set<std::uint32_t>
        represented_volume_tags;
    for (const auto& element :
         volume_elements) {
        LinearCellType3D type =
            LinearCellType3D::tetrahedron;
        if (element.element_type == 4 &&
            element.node_tags.size() == 4U) {
            type =
                LinearCellType3D::tetrahedron;
        } else if (
            element.element_type == 5 &&
            element.node_tags.size() == 8U) {
            type =
                LinearCellType3D::hexahedron;
        } else if (
            element.element_type == 6 &&
            element.node_tags.size() == 6U) {
            type =
                LinearCellType3D::wedge;
        } else if (
            element.element_type == 7 &&
            element.node_tags.size() == 5U) {
            type =
                LinearCellType3D::pyramid;
        } else {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii_3d: volume baseline supports only linear tetrahedron, hexahedron, prism/wedge and pyramid");
        }

        std::vector<LocalIndex> vertices;
        vertices.reserve(
            element.node_tags.size());
        for (const auto tag :
             element.node_tags) {
            vertices.push_back(
                lookup_node(
                    node_local_by_tag,
                    tag));
        }
        cells.push_back(
            LinearCell3D{
                GlobalEntityId{
                    element.tag},
                type,
                std::move(vertices)});

        std::vector<std::uint32_t> groups;
        if (parsed.entities_seen) {
            const auto* entity =
                find_entity(
                    parsed.entities,
                    3,
                    element.entity_tag);
            if (entity == nullptr) {
                throw std::logic_error(
                    "mpmc::mesh::import_gmsh_4_1_ascii_3d: volume entity lookup drift");
            }
            groups =
                entity->physical_tags;
            represented_volume_tags.insert(
                groups.begin(),
                groups.end());
        }
        cell_groups.push_back(
            GmshCellPhysicalGroups{
                GlobalEntityId{
                    element.tag},
                std::move(groups)});
    }

    std::vector<LinearFaceAnnotation3D>
        annotations;
    annotations.reserve(
        surface_elements.size());
    std::set<std::uint32_t>
        represented_surface_tags;
    for (const auto& element :
         surface_elements) {
        if ((element.element_type != 2 ||
             element.node_tags.size() != 3U) &&
            (element.element_type != 3 ||
             element.node_tags.size() != 4U)) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii_3d: surface baseline supports only linear triangle and quadrilateral faces");
        }

        std::uint32_t physical_tag = 0U;
        if (parsed.entities_seen) {
            const auto* entity =
                find_entity(
                    parsed.entities,
                    2,
                    element.entity_tag);
            if (entity == nullptr) {
                throw std::logic_error(
                    "mpmc::mesh::import_gmsh_4_1_ascii_3d: surface entity lookup drift");
            }
            if (entity->physical_tags.size() >
                1U) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii_3d: one face cannot map multiple Physical Groups in FaceBoundarySnapshot");
            }
            if (!entity->physical_tags.empty()) {
                physical_tag =
                    entity->physical_tags.front();
                represented_surface_tags.insert(
                    physical_tag);
            }
        }

        std::vector<LocalIndex> vertices;
        vertices.reserve(
            element.node_tags.size());
        for (const auto tag :
             element.node_tags) {
            vertices.push_back(
                lookup_node(
                    node_local_by_tag,
                    tag));
        }
        annotations.push_back(
            LinearFaceAnnotation3D{
                std::move(vertices),
                GlobalEntityId{
                    element.tag},
                PhysicalTag{
                    physical_tag}});
    }

    auto mesh =
        make_linear_mesh_3d(
            vertex_ids,
            coordinates,
            cells,
            annotations,
            maximum_element_tag);

    std::vector<GmshPhysicalName>
        physical_names;
    for (const auto& [key, name] :
         parsed.physical_names) {
        if (key.dimension == 2) {
            const auto tag =
                gmsh41_detail::checked_physical_tag(
                    key.tag,
                    "mpmc::mesh::import_gmsh_4_1_ascii_3d: invalid surface PhysicalName tag");
            if (represented_surface_tags.contains(
                    tag)) {
                physical_names.push_back(
                    GmshPhysicalName{
                        2, tag, name});
            }
        } else if (key.dimension == 3) {
            const auto tag =
                gmsh41_detail::checked_physical_tag(
                    key.tag,
                    "mpmc::mesh::import_gmsh_4_1_ascii_3d: invalid volume PhysicalName tag");
            if (represented_volume_tags.contains(
                    tag)) {
                physical_names.push_back(
                    GmshPhysicalName{
                        3, tag, name});
            }
        }
    }
    std::sort(
        physical_names.begin(),
        physical_names.end(),
        [](const GmshPhysicalName& left,
           const GmshPhysicalName& right) {
            return std::tie(
                       left.dimension,
                       left.tag) <
                   std::tie(
                       right.dimension,
                       right.tag);
        });

    return Gmsh41ImportResult3D{
        std::move(mesh.topology),
        std::move(
            mesh.vertex_coordinates_m),
        std::move(mesh.cell_volumes_m3),
        std::move(mesh.face_geometry),
        std::move(mesh.face_boundary),
        std::move(physical_names),
        std::move(cell_groups)};
}

[[nodiscard]] inline std::string
export_gmsh_4_1_ascii_3d(
    const Gmsh41ImportResult3D& mesh) {
    using namespace gmsh41_3d_detail;

    const auto& topology = mesh.topology;
    const std::size_t vertex_count =
        topology.entity_count(
            EntityKind::vertex);
    const std::size_t face_count =
        topology.entity_count(
            EntityKind::face);
    const std::size_t cell_count =
        topology.entity_count(
            EntityKind::cell);
    if (vertex_count == 0U ||
        face_count == 0U ||
        cell_count == 0U ||
        topology.entity_count(
            EntityKind::edge) != 0U) {
        throw std::invalid_argument(
            "mpmc::mesh::export_gmsh_4_1_ascii_3d: nonempty vertex/face/cell topology without materialized edges is required");
    }
    if (mesh.vertex_coordinates_m.size() !=
            vertex_count ||
        mesh.cell_volumes_m3.size() !=
            cell_count ||
        mesh.face_geometry.face_count() !=
            face_count ||
        mesh.face_boundary.face_count() !=
            face_count) {
        throw std::invalid_argument(
            "mpmc::mesh::export_gmsh_4_1_ascii_3d: topology and 3D geometry snapshots are not aligned");
    }
    for (const auto& relation :
         {std::pair{EntityKind::cell,
                    EntityKind::vertex},
          std::pair{EntityKind::cell,
                    EntityKind::face},
          std::pair{EntityKind::face,
                    EntityKind::vertex},
          std::pair{EntityKind::face,
                    EntityKind::cell}}) {
        if (!topology.has_relation(
                relation.first,
                relation.second)) {
            throw std::invalid_argument(
                "mpmc::mesh::export_gmsh_4_1_ascii_3d: required 3D topology relation is missing");
        }
    }

    const auto& cell_vertices =
        topology.relation(
            EntityKind::cell,
            EntityKind::vertex);
    const auto& face_vertices =
        topology.relation(
            EntityKind::face,
            EntityKind::vertex);
    const auto vertex_ids =
        topology.global_ids(
            EntityKind::vertex);
    const auto face_ids =
        topology.global_ids(
            EntityKind::face);
    const auto cell_ids =
        topology.global_ids(
            EntityKind::cell);

    std::set<std::uint64_t>
        element_ids;
    for (const auto id : face_ids) {
        if (id.value() == 0U ||
            !element_ids.insert(
                id.value()).second) {
            throw std::invalid_argument(
                "mpmc::mesh::export_gmsh_4_1_ascii_3d: face/cell element tags must be unique and positive");
        }
    }
    for (const auto id : cell_ids) {
        if (id.value() == 0U ||
            !element_ids.insert(
                id.value()).second) {
            throw std::invalid_argument(
                "mpmc::mesh::export_gmsh_4_1_ascii_3d: face/cell element tags must be globally unique and positive");
        }
    }
    for (const auto id : vertex_ids) {
        if (id.value() == 0U) {
            throw std::invalid_argument(
                "mpmc::mesh::export_gmsh_4_1_ascii_3d: node tags must be positive");
        }
    }

    std::map<std::uint64_t,
             std::vector<std::uint32_t>>
        cell_groups_by_id;
    for (const auto& record :
         mesh.cell_physical_groups) {
        auto groups =
            record.physical_tags;
        std::sort(
            groups.begin(),
            groups.end());
        if (std::adjacent_find(
                groups.begin(),
                groups.end()) !=
            groups.end() ||
            std::find(
                groups.begin(),
                groups.end(),
                0U) != groups.end() ||
            !cell_groups_by_id.emplace(
                record.cell_global_id.value(),
                std::move(groups))
                 .second) {
            throw std::invalid_argument(
                "mpmc::mesh::export_gmsh_4_1_ascii_3d: invalid or duplicate cell Physical Group record");
        }
    }
    for (const auto id : cell_ids) {
        if (!cell_groups_by_id.contains(
                id.value())) {
            cell_groups_by_id.emplace(
                id.value(),
                std::vector<std::uint32_t>{});
        }
    }

    std::set<std::pair<int, std::uint32_t>>
        physical_name_keys;
    std::set<std::uint32_t>
        represented_surface_tags;
    std::set<std::uint32_t>
        represented_volume_tags;
    for (std::size_t face = 0U;
         face < face_count;
         ++face) {
        const auto tag =
            mesh.face_boundary.physical_tag(
                gmsh41_detail::checked_local(
                    face,
                    "mpmc::mesh::export_gmsh_4_1_ascii_3d: face local index overflow"));
        if (tag.is_tagged()) {
            represented_surface_tags.insert(
                tag.value());
        }
    }
    for (const auto& [id, groups] :
         cell_groups_by_id) {
        (void)id;
        represented_volume_tags.insert(
            groups.begin(),
            groups.end());
    }
    auto physical_names =
        mesh.physical_names;
    std::sort(
        physical_names.begin(),
        physical_names.end(),
        [](const GmshPhysicalName& left,
           const GmshPhysicalName& right) {
            return std::tie(
                       left.dimension,
                       left.tag) <
                   std::tie(
                       right.dimension,
                       right.tag);
        });
    for (const auto& name :
         physical_names) {
        require_export_name(name);
        if (!physical_name_keys
                 .insert(
                     {name.dimension,
                      name.tag})
                 .second) {
            throw std::invalid_argument(
                "mpmc::mesh::export_gmsh_4_1_ascii_3d: duplicate PhysicalName key");
        }
        const bool represented =
            name.dimension == 2
                ? represented_surface_tags
                      .contains(name.tag)
                : represented_volume_tags
                      .contains(name.tag);
        if (!represented) {
            throw std::invalid_argument(
                "mpmc::mesh::export_gmsh_4_1_ascii_3d: PhysicalName has no represented group membership");
        }
    }

    std::ostringstream output;
    output << std::setprecision(
        std::numeric_limits<double>::max_digits10);
    output << "$MeshFormat\n"
           << "4.1 0 8\n"
           << "$EndMeshFormat\n";

    if (!physical_names.empty()) {
        output << "$PhysicalNames\n"
               << physical_names.size()
               << '\n';
        for (const auto& name :
             physical_names) {
            output << name.dimension << ' '
                   << name.tag << ' '
                   << std::quoted(name.name)
                   << '\n';
        }
        output << "$EndPhysicalNames\n";
    }

    output << "$Entities\n"
           << vertex_count << " 0 "
           << face_count << ' '
           << cell_count << '\n';

    for (std::size_t vertex = 0U;
         vertex < vertex_count;
         ++vertex) {
        const auto point =
            mesh.vertex_coordinates_m[vertex];
        output << (vertex + 1U) << ' '
               << point.x_m << ' '
               << point.y_m << ' '
               << point.z_m
               << " 0\n";
    }

    for (std::size_t face = 0U;
         face < face_count;
         ++face) {
        const auto local =
            gmsh41_detail::checked_local(
                face,
                "mpmc::mesh::export_gmsh_4_1_ascii_3d: face local index overflow");
        const auto vertices =
            face_vertices.adjacent(local);
        if (vertices.size() != 3U &&
            vertices.size() != 4U) {
            throw std::invalid_argument(
                "mpmc::mesh::export_gmsh_4_1_ascii_3d: faces must be linear triangles or quadrilaterals");
        }
        const auto box =
            bounds(
                vertices,
                mesh.vertex_coordinates_m);
        const auto physical_tag =
            mesh.face_boundary
                .physical_tag(local);
        output << (face + 1U);
        for (const auto value : box) {
            output << ' ' << value;
        }
        if (physical_tag.is_tagged()) {
            output << " 1 "
                   << physical_tag.value();
        } else {
            output << " 0";
        }
        output << " 0\n";
    }

    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        const auto local =
            gmsh41_detail::checked_local(
                cell,
                "mpmc::mesh::export_gmsh_4_1_ascii_3d: cell local index overflow");
        const auto vertices =
            cell_vertices.adjacent(local);
        if (vertices.size() != 4U &&
            vertices.size() != 5U &&
            vertices.size() != 6U &&
            vertices.size() != 8U) {
            throw std::invalid_argument(
                "mpmc::mesh::export_gmsh_4_1_ascii_3d: cells must be linear tetrahedra, pyramids, prisms/wedges or hexahedra");
        }
        const auto box =
            bounds(
                vertices,
                mesh.vertex_coordinates_m);
        const auto& groups =
            cell_groups_by_id.at(
                cell_ids[cell].value());
        output << (cell + 1U);
        for (const auto value : box) {
            output << ' ' << value;
        }
        output << ' ' << groups.size();
        for (const auto tag : groups) {
            output << ' ' << tag;
        }
        output << " 0\n";
    }
    output << "$EndEntities\n";

    const auto [node_minimum, node_maximum] =
        std::minmax_element(
            vertex_ids.begin(),
            vertex_ids.end(),
            [](GlobalEntityId left,
               GlobalEntityId right) {
                return left.value() <
                       right.value();
            });
    output << "$Nodes\n"
           << vertex_count << ' '
           << vertex_count << ' '
           << node_minimum->value() << ' '
           << node_maximum->value()
           << '\n';
    for (std::size_t vertex = 0U;
         vertex < vertex_count;
         ++vertex) {
        const auto point =
            mesh.vertex_coordinates_m[vertex];
        output << "0 " << (vertex + 1U)
               << " 0 1\n"
               << vertex_ids[vertex].value()
               << '\n'
               << point.x_m << ' '
               << point.y_m << ' '
               << point.z_m << '\n';
    }
    output << "$EndNodes\n";

    const auto element_minimum =
        *element_ids.begin();
    const auto element_maximum =
        *element_ids.rbegin();
    output << "$Elements\n"
           << (face_count + cell_count)
           << ' '
           << (face_count + cell_count)
           << ' '
           << element_minimum << ' '
           << element_maximum << '\n';

    for (std::size_t face = 0U;
         face < face_count;
         ++face) {
        const auto local =
            gmsh41_detail::checked_local(
                face,
                "mpmc::mesh::export_gmsh_4_1_ascii_3d: face local index overflow");
        const auto vertices =
            face_vertices.adjacent(local);
        const int element_type =
            vertices.size() == 3U
                ? 2
                : 3;
        output << "2 " << (face + 1U)
               << ' ' << element_type
               << " 1\n"
               << face_ids[face].value();
        for (const auto vertex :
             vertices) {
            output << ' '
                   << vertex_ids[
                          static_cast<std::size_t>(
                              vertex.value())]
                          .value();
        }
        output << '\n';
    }

    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        const auto local =
            gmsh41_detail::checked_local(
                cell,
                "mpmc::mesh::export_gmsh_4_1_ascii_3d: cell local index overflow");
        const auto vertices =
            cell_vertices.adjacent(local);
        int element_type = 0;
        if (vertices.size() == 4U) {
            element_type = 4;
        } else if (vertices.size() == 8U) {
            element_type = 5;
        } else if (vertices.size() == 6U) {
            element_type = 6;
        } else if (vertices.size() == 5U) {
            element_type = 7;
        } else {
            throw std::logic_error(
                "mpmc::mesh::export_gmsh_4_1_ascii_3d: validated cell width drift");
        }
        output << "3 " << (cell + 1U)
               << ' ' << element_type
               << " 1\n"
               << cell_ids[cell].value();
        for (const auto vertex :
             vertices) {
            output << ' '
                   << vertex_ids[
                          static_cast<std::size_t>(
                              vertex.value())]
                          .value();
        }
        output << '\n';
    }
    output << "$EndElements\n";

    return output.str();
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_GMSH_4_1_3D_HPP
