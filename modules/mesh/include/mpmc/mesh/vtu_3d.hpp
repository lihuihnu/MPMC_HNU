#ifndef MPMC_MESH_VTU_3D_HPP
#define MPMC_MESH_VTU_3D_HPP

#include <mpmc/mesh/linear_cell_mesh_3d.hpp>
#include <mpmc/mesh/vtu.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::mesh {

struct VtuImportResult3D {
    Topology topology;
    std::vector<Coordinate3D> vertex_coordinates_m;
    std::vector<double> cell_volumes_m3;
    FaceGeometry3D face_geometry;
    FaceBoundarySnapshot face_boundary;
    std::vector<DenseFieldSnapshot> point_fields;
    std::vector<DenseFieldSnapshot> cell_fields;
};

[[nodiscard]] inline VtuImportResult3D
import_vtu_ascii_3d(
    std::string_view content) {
    using namespace vtu_detail;

    if (content.find("<AppendedData") !=
            std::string_view::npos ||
        content.find("format=\"binary\"") !=
            std::string_view::npos ||
        content.find("format='binary'") !=
            std::string_view::npos ||
        content.find("format=\"appended\"") !=
            std::string_view::npos ||
        content.find("format='appended'") !=
            std::string_view::npos) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii_3d: binary/appended VTU is unsupported");
    }

    const auto vtk_file =
        require_unique_element(
            content,
            "VTKFile",
            "mpmc::mesh::import_vtu_ascii_3d: VTKFile root is required");
    const auto& vtk_type =
        require_attribute(
            vtk_file,
            "type",
            "mpmc::mesh::import_vtu_ascii_3d: VTKFile type is required");
    if (vtk_type != "UnstructuredGrid") {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii_3d: only VTKFile type=UnstructuredGrid is supported");
    }
    if (vtk_file.attributes.contains(
            "compressor")) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii_3d: compressed VTU is unsupported");
    }

    const auto grid =
        require_unique_element(
            vtk_file.body,
            "UnstructuredGrid",
            "mpmc::mesh::import_vtu_ascii_3d: UnstructuredGrid element is required");
    const auto piece =
        require_unique_element(
            grid.body,
            "Piece",
            "mpmc::mesh::import_vtu_ascii_3d: exactly one Piece is required");
    const std::size_t point_count =
        parse_size(
            require_attribute(
                piece,
                "NumberOfPoints",
                "mpmc::mesh::import_vtu_ascii_3d: Piece NumberOfPoints is required"),
            "mpmc::mesh::import_vtu_ascii_3d: invalid NumberOfPoints");
    const std::size_t cell_count =
        parse_size(
            require_attribute(
                piece,
                "NumberOfCells",
                "mpmc::mesh::import_vtu_ascii_3d: Piece NumberOfCells is required"),
            "mpmc::mesh::import_vtu_ascii_3d: invalid NumberOfCells");
    if (point_count == 0U ||
        cell_count == 0U) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii_3d: Piece must contain points and cells");
    }

    const auto points_section =
        require_unique_element(
            piece.body,
            "Points",
            "mpmc::mesh::import_vtu_ascii_3d: Points section is required");
    const auto point_arrays =
        data_arrays(
            points_section.body);
    if (point_arrays.size() != 1U) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii_3d: Points must contain exactly one DataArray");
    }
    const auto& point_array =
        point_arrays.front();
    require_ascii_data_array(point_array);
    const auto& point_type =
        require_attribute(
            point_array,
            "type",
            "mpmc::mesh::import_vtu_ascii_3d: Points DataArray type is required");
    if (point_type != "Float32" &&
        point_type != "Float64") {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii_3d: Points must use Float32 or Float64");
    }
    if (component_count(
            point_array) != 3U) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii_3d: Points must have NumberOfComponents=3");
    }
    const auto point_values =
        parse_double_values(
            point_array.body,
            "mpmc::mesh::import_vtu_ascii_3d: invalid point coordinates");
    if (point_values.size() !=
        point_count * 3U) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii_3d: point coordinate count mismatch");
    }

    std::vector<Coordinate3D>
        coordinates;
    std::vector<GlobalEntityId>
        vertex_ids;
    coordinates.reserve(point_count);
    vertex_ids.reserve(point_count);
    for (std::size_t point = 0U;
         point < point_count;
         ++point) {
        coordinates.push_back(
            Coordinate3D{
                point_values[
                    point * 3U],
                point_values[
                    point * 3U + 1U],
                point_values[
                    point * 3U + 2U]});
        vertex_ids.emplace_back(
            static_cast<std::uint64_t>(
                point) +
            1U);
    }

    const auto cells_section =
        require_unique_element(
            piece.body,
            "Cells",
            "mpmc::mesh::import_vtu_ascii_3d: Cells section is required");
    const auto cell_arrays =
        data_arrays(
            cells_section.body);
    if (cell_arrays.size() != 3U) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii_3d: Cells must contain connectivity, offsets and types DataArrays");
    }

    std::optional<XmlElement>
        connectivity_array;
    std::optional<XmlElement>
        offsets_array;
    std::optional<XmlElement>
        types_array;
    for (const auto& array :
         cell_arrays) {
        require_ascii_data_array(array);
        const auto& name =
            require_attribute(
                array,
                "Name",
                "mpmc::mesh::import_vtu_ascii_3d: Cells DataArray Name is required");
        if (name == "connectivity") {
            if (connectivity_array.has_value()) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_vtu_ascii_3d: duplicate connectivity array");
            }
            connectivity_array = array;
        } else if (name == "offsets") {
            if (offsets_array.has_value()) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_vtu_ascii_3d: duplicate offsets array");
            }
            offsets_array = array;
        } else if (name == "types") {
            if (types_array.has_value()) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_vtu_ascii_3d: duplicate types array");
            }
            types_array = array;
        } else {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii_3d: unsupported Cells DataArray");
        }
    }
    if (!connectivity_array.has_value() ||
        !offsets_array.has_value() ||
        !types_array.has_value()) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii_3d: incomplete Cells arrays");
    }

    const auto& connectivity_type =
        require_attribute(
            *connectivity_array,
            "type",
            "mpmc::mesh::import_vtu_ascii_3d: connectivity type is required");
    const auto& offsets_type =
        require_attribute(
            *offsets_array,
            "type",
            "mpmc::mesh::import_vtu_ascii_3d: offsets type is required");
    const auto& types_type =
        require_attribute(
            *types_array,
            "type",
            "mpmc::mesh::import_vtu_ascii_3d: types type is required");
    if ((connectivity_type != "Int32" &&
         connectivity_type != "Int64") ||
        (offsets_type != "Int32" &&
         offsets_type != "Int64") ||
        types_type != "UInt8") {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii_3d: Cells requires Int32/Int64 connectivity+offsets and UInt8 types");
    }

    const auto connectivity =
        parse_integer_values<std::int64_t>(
            connectivity_array->body,
            "mpmc::mesh::import_vtu_ascii_3d: invalid connectivity");
    const auto offsets =
        parse_integer_values<std::int64_t>(
            offsets_array->body,
            "mpmc::mesh::import_vtu_ascii_3d: invalid offsets");
    const auto types =
        parse_integer_values<std::int64_t>(
            types_array->body,
            "mpmc::mesh::import_vtu_ascii_3d: invalid types");
    if (offsets.size() != cell_count ||
        types.size() != cell_count) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii_3d: offsets/types count mismatch");
    }

    std::vector<
        std::vector<LocalIndex>>
        cell_vertices;
    std::vector<LinearCellType3D>
        cell_types;
    cell_vertices.reserve(cell_count);
    cell_types.reserve(cell_count);
    std::size_t previous_offset = 0U;
    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        if (offsets[cell] < 0) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii_3d: cell offset cannot be negative");
        }
        const auto raw_end =
            static_cast<std::uint64_t>(
                offsets[cell]);
        if (raw_end >
            static_cast<std::uint64_t>(
                std::numeric_limits<
                    std::size_t>::max())) {
            throw std::length_error(
                "mpmc::mesh::import_vtu_ascii_3d: cell offset exceeds size_t");
        }
        const std::size_t end =
            static_cast<std::size_t>(
                raw_end);
        if (end <= previous_offset ||
            end > connectivity.size()) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii_3d: invalid monotone cell offsets");
        }

        const std::size_t width =
            end - previous_offset;
        LinearCellType3D type =
            LinearCellType3D::tetrahedron;
        std::size_t expected = 0U;
        if (types[cell] == 10) {
            type =
                LinearCellType3D::tetrahedron;
            expected = 4U;
        } else if (types[cell] == 12) {
            type =
                LinearCellType3D::hexahedron;
            expected = 8U;
        }
        if (expected == 0U ||
            width != expected) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii_3d: only VTK_TETRA(10) and VTK_HEXAHEDRON(12) are supported with matching connectivity widths");
        }

        std::vector<LocalIndex> vertices;
        vertices.reserve(width);
        std::set<std::size_t> unique;
        for (std::size_t index =
                 previous_offset;
             index < end;
             ++index) {
            if (connectivity[index] < 0 ||
                static_cast<std::uint64_t>(
                    connectivity[index]) >=
                    static_cast<std::uint64_t>(
                        point_count)) {
                throw std::out_of_range(
                    "mpmc::mesh::import_vtu_ascii_3d: connectivity point index out of range");
            }
            const auto point =
                static_cast<std::size_t>(
                    connectivity[index]);
            if (!unique.insert(
                    point).second) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_vtu_ascii_3d: cell contains repeated point index");
            }
            vertices.push_back(
                local_index(
                    point,
                    "mpmc::mesh::import_vtu_ascii_3d: point local index overflow"));
        }
        cell_vertices.push_back(
            std::move(vertices));
        cell_types.push_back(type);
        previous_offset = end;
    }
    if (previous_offset !=
        connectivity.size()) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii_3d: final cell offset must consume all connectivity");
    }

    std::vector<std::uint64_t>
        cell_global_ids;
    std::vector<ParsedField>
        parsed_cell_fields;
    const auto cell_data =
        find_element(
            piece.body,
            "CellData");
    if (cell_data.has_value()) {
        std::set<std::string>
            field_names;
        for (const auto& array :
             data_arrays(
                 cell_data->body)) {
            require_ascii_data_array(array);
            const auto& name =
                require_attribute(
                    array,
                    "Name",
                    "mpmc::mesh::import_vtu_ascii_3d: CellData DataArray Name is required");
            if (!field_names.insert(
                    name).second) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_vtu_ascii_3d: duplicate CellData Name");
            }

            if (name ==
                global_cell_id_name) {
                const auto& type =
                    require_attribute(
                        array,
                        "type",
                        "mpmc::mesh::import_vtu_ascii_3d: global cell ID type is required");
                if (type != "UInt64" &&
                    type != "Int64") {
                    throw std::invalid_argument(
                        "mpmc::mesh::import_vtu_ascii_3d: mpmc_global_cell_id must use UInt64 or Int64");
                }
                std::set<std::uint64_t>
                    unique;
                if (type == "UInt64") {
                    const auto values =
                        parse_integer_values<
                            std::uint64_t>(
                            array.body,
                            "mpmc::mesh::import_vtu_ascii_3d: invalid UInt64 global cell IDs");
                    if (values.size() !=
                        cell_count) {
                        throw std::invalid_argument(
                            "mpmc::mesh::import_vtu_ascii_3d: global cell ID count mismatch");
                    }
                    for (const auto id :
                         values) {
                        if (!unique.insert(
                                id)
                                 .second) {
                            throw std::invalid_argument(
                                "mpmc::mesh::import_vtu_ascii_3d: global cell IDs must be unique");
                        }
                        cell_global_ids.push_back(
                            id);
                    }
                } else {
                    const auto values =
                        parse_integer_values<
                            std::int64_t>(
                            array.body,
                            "mpmc::mesh::import_vtu_ascii_3d: invalid Int64 global cell IDs");
                    if (values.size() !=
                        cell_count) {
                        throw std::invalid_argument(
                            "mpmc::mesh::import_vtu_ascii_3d: global cell ID count mismatch");
                    }
                    for (const auto value :
                         values) {
                        if (value < 0) {
                            throw std::invalid_argument(
                                "mpmc::mesh::import_vtu_ascii_3d: global cell ID cannot be negative");
                        }
                        const auto id =
                            static_cast<
                                std::uint64_t>(
                                value);
                        if (!unique.insert(
                                id)
                                 .second) {
                            throw std::invalid_argument(
                                "mpmc::mesh::import_vtu_ascii_3d: global cell IDs must be unique");
                        }
                        cell_global_ids.push_back(
                            id);
                    }
                }
            } else {
                parsed_cell_fields.push_back(
                    parse_field(
                        array,
                        EntityKind::cell,
                        cell_count));
            }
        }
    }
    if (cell_global_ids.empty()) {
        cell_global_ids.reserve(
            cell_count);
        for (std::size_t cell = 0U;
             cell < cell_count;
             ++cell) {
            cell_global_ids.push_back(
                static_cast<std::uint64_t>(
                    cell) +
                1U);
        }
    }

    std::vector<ParsedField>
        parsed_point_fields;
    const auto point_data =
        find_element(
            piece.body,
            "PointData");
    if (point_data.has_value()) {
        std::set<std::string>
            field_names;
        for (const auto& array :
             data_arrays(
                 point_data->body)) {
            const auto& name =
                require_attribute(
                    array,
                    "Name",
                    "mpmc::mesh::import_vtu_ascii_3d: PointData DataArray Name is required");
            if (!field_names.insert(
                    name).second) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_vtu_ascii_3d: duplicate PointData Name");
            }
            parsed_point_fields.push_back(
                parse_field(
                    array,
                    EntityKind::vertex,
                    point_count));
        }
    }

    std::vector<LinearCell3D> cells;
    cells.reserve(cell_count);
    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        cells.push_back(
            LinearCell3D{
                GlobalEntityId{
                    cell_global_ids[cell]},
                cell_types[cell],
                std::move(
                    cell_vertices[cell])});
    }

    auto mesh =
        make_linear_mesh_3d(
            vertex_ids,
            coordinates,
            cells);

    std::vector<DenseFieldSnapshot>
        point_fields;
    point_fields.reserve(
        parsed_point_fields.size());
    for (auto& field :
         parsed_point_fields) {
        point_fields.push_back(
            DenseFieldSnapshot::create(
                mesh.topology,
                EntityKind::vertex,
                field.component_count,
                std::move(field.values),
                std::move(field.metadata)));
    }

    std::vector<DenseFieldSnapshot>
        cell_fields;
    cell_fields.reserve(
        parsed_cell_fields.size());
    for (auto& field :
         parsed_cell_fields) {
        cell_fields.push_back(
            DenseFieldSnapshot::create(
                mesh.topology,
                EntityKind::cell,
                field.component_count,
                std::move(field.values),
                std::move(field.metadata)));
    }

    return VtuImportResult3D{
        std::move(mesh.topology),
        std::move(
            mesh.vertex_coordinates_m),
        std::move(mesh.cell_volumes_m3),
        std::move(mesh.face_geometry),
        std::move(mesh.face_boundary),
        std::move(point_fields),
        std::move(cell_fields)};
}

[[nodiscard]] inline std::string
export_vtu_ascii_3d(
    const VtuImportResult3D& mesh) {
    using namespace vtu_detail;

    const auto& topology = mesh.topology;
    const std::size_t point_count =
        topology.entity_count(
            EntityKind::vertex);
    const std::size_t cell_count =
        topology.entity_count(
            EntityKind::cell);
    if (point_count == 0U ||
        cell_count == 0U ||
        topology.entity_count(
            EntityKind::edge) != 0U) {
        throw std::invalid_argument(
            "mpmc::mesh::export_vtu_ascii_3d: nonempty 3D vertex/face/cell mesh without materialized edges is required");
    }
    if (mesh.vertex_coordinates_m.size() !=
            point_count ||
        mesh.cell_volumes_m3.size() !=
            cell_count ||
        mesh.face_geometry.face_count() !=
            topology.entity_count(
                EntityKind::face) ||
        mesh.face_boundary.face_count() !=
            topology.entity_count(
                EntityKind::face)) {
        throw std::invalid_argument(
            "mpmc::mesh::export_vtu_ascii_3d: topology and geometry counts are not aligned");
    }
    if (!topology.has_relation(
            EntityKind::cell,
            EntityKind::vertex)) {
        throw std::invalid_argument(
            "mpmc::mesh::export_vtu_ascii_3d: cell->vertex relation is required");
    }

    validate_export_fields(
        mesh.point_fields,
        EntityKind::vertex,
        point_count);
    validate_export_fields(
        mesh.cell_fields,
        EntityKind::cell,
        cell_count);

    const auto& cell_vertices =
        topology.relation(
            EntityKind::cell,
            EntityKind::vertex);
    std::vector<std::size_t> offsets;
    std::vector<int> types;
    offsets.reserve(cell_count);
    types.reserve(cell_count);
    std::size_t connectivity_count = 0U;
    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        const auto vertices =
            cell_vertices.adjacent(
                local_index(
                    cell,
                    "mpmc::mesh::export_vtu_ascii_3d: cell local index overflow"));
        int type = 0;
        if (vertices.size() == 4U) {
            type = 10;
        } else if (
            vertices.size() == 8U) {
            type = 12;
        } else {
            throw std::invalid_argument(
                "mpmc::mesh::export_vtu_ascii_3d: only linear tetrahedron/hexahedron cells are supported");
        }
        if (connectivity_count >
            std::numeric_limits<
                std::size_t>::max() -
                vertices.size()) {
            throw std::length_error(
                "mpmc::mesh::export_vtu_ascii_3d: connectivity size overflow");
        }
        connectivity_count +=
            vertices.size();
        offsets.push_back(
            connectivity_count);
        types.push_back(type);
    }

    std::ostringstream output;
    output << std::setprecision(
        std::numeric_limits<double>::max_digits10);
    output << "<?xml version=\"1.0\"?>\n"
           << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
           << "  <UnstructuredGrid>\n"
           << "    <Piece NumberOfPoints=\""
           << point_count
           << "\" NumberOfCells=\""
           << cell_count
           << "\">\n";

    output << "      <PointData>\n";
    for (const auto& field :
         mesh.point_fields) {
        write_field(
            output,
            field,
            8U);
    }
    output << "      </PointData>\n";

    output << "      <CellData>\n"
           << "        <DataArray type=\"UInt64\" Name=\""
           << global_cell_id_name
           << "\" format=\"ascii\">\n"
           << "          ";
    const auto cell_ids =
        topology.global_ids(
            EntityKind::cell);
    for (std::size_t cell = 0U;
         cell < cell_ids.size();
         ++cell) {
        if (cell != 0U) {
            output << ' ';
        }
        output << cell_ids[cell].value();
    }
    output << "\n        </DataArray>\n";
    for (const auto& field :
         mesh.cell_fields) {
        write_field(
            output,
            field,
            8U);
    }
    output << "      </CellData>\n";

    output << "      <Points>\n"
           << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n"
           << "          ";
    for (std::size_t point = 0U;
         point < point_count;
         ++point) {
        if (point != 0U) {
            output << ' ';
        }
        const auto coordinate =
            mesh.vertex_coordinates_m[
                point];
        output << coordinate.x_m
               << ' '
               << coordinate.y_m
               << ' '
               << coordinate.z_m;
    }
    output << "\n        </DataArray>\n"
           << "      </Points>\n";

    output << "      <Cells>\n"
           << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n"
           << "          ";
    bool first_connectivity = true;
    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        const auto vertices =
            cell_vertices.adjacent(
                local_index(
                    cell,
                    "mpmc::mesh::export_vtu_ascii_3d: cell local index overflow"));
        for (const auto vertex :
             vertices) {
            if (!first_connectivity) {
                output << ' ';
            }
            first_connectivity = false;
            output << vertex.value();
        }
    }
    output << "\n        </DataArray>\n"
           << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n"
           << "          ";
    for (std::size_t index = 0U;
         index < offsets.size();
         ++index) {
        if (index != 0U) {
            output << ' ';
        }
        output << offsets[index];
    }
    output << "\n        </DataArray>\n"
           << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n"
           << "          ";
    for (std::size_t index = 0U;
         index < types.size();
         ++index) {
        if (index != 0U) {
            output << ' ';
        }
        output << types[index];
    }
    output << "\n        </DataArray>\n"
           << "      </Cells>\n"
           << "    </Piece>\n"
           << "  </UnstructuredGrid>\n"
           << "</VTKFile>\n";

    return output.str();
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_VTU_3D_HPP
