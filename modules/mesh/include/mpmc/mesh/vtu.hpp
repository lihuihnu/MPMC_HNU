#ifndef MPMC_MESH_VTU_HPP
#define MPMC_MESH_VTU_HPP

#include <mpmc/mesh/csr_adjacency.hpp>
#include <mpmc/mesh/dense_field.hpp>
#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/geometry_2d.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace mpmc::mesh {

struct VtuImportResult {
    Topology topology;
    Geometry2D geometry;
    std::vector<DenseFieldSnapshot> point_fields;
    std::vector<DenseFieldSnapshot> cell_fields;
};

namespace vtu_detail {

inline constexpr std::string_view global_cell_id_name =
    "mpmc_global_cell_id";

struct XmlElement {
    std::map<std::string, std::string> attributes;
    std::string_view body;
    std::size_t next_offset;
    bool self_closing;
};

struct EdgeKey {
    std::size_t first;
    std::size_t second;

    [[nodiscard]] friend bool operator<(const EdgeKey& left,
                                        const EdgeKey& right) noexcept {
        return std::tie(left.first, left.second) <
               std::tie(right.first, right.second);
    }
};

struct FaceBuild {
    EdgeKey key;
    std::vector<std::size_t> adjacent_cells;
};

struct CellMetric {
    Coordinate2D centroid;
    double area;
};

struct ParsedField {
    std::string name;
    std::size_t component_count;
    std::vector<double> values;
    DenseFieldMetadata metadata;
};

[[nodiscard]] inline bool is_xml_space(char value) noexcept {
    return value == ' ' || value == '\t' ||
           value == '\n' || value == '\r';
}

[[nodiscard]] inline std::string xml_unescape(std::string_view value) {
    std::string result;
    result.reserve(value.size());
    for (std::size_t i = 0U; i < value.size(); ++i) {
        if (value[i] != '&') {
            result.push_back(value[i]);
            continue;
        }
        const auto semicolon = value.find(';', i + 1U);
        if (semicolon == std::string_view::npos) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: unterminated XML entity");
        }
        const auto entity =
            value.substr(i, semicolon - i + 1U);
        if (entity == "&amp;") result.push_back('&');
        else if (entity == "&lt;") result.push_back('<');
        else if (entity == "&gt;") result.push_back('>');
        else if (entity == "&quot;") result.push_back('"');
        else if (entity == "&apos;") result.push_back('\'');
        else {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: unsupported XML entity");
        }
        i = semicolon;
    }
    return result;
}

[[nodiscard]] inline std::string xml_escape(std::string_view value) {
    std::string result;
    result.reserve(value.size());
    for (const char character : value) {
        switch (character) {
        case '&': result += "&amp;"; break;
        case '<': result += "&lt;"; break;
        case '>': result += "&gt;"; break;
        case '"': result += "&quot;"; break;
        case '\'': result += "&apos;"; break;
        case '\0':
            throw std::invalid_argument(
                "mpmc::mesh::export_vtu_ascii: XML attribute text cannot contain NUL");
        default:
            result.push_back(character);
            break;
        }
    }
    return result;
}

[[nodiscard]] inline std::map<std::string, std::string>
parse_attributes(std::string_view open_tag,
                 std::string_view element_name) {
    const std::size_t name_start =
        open_tag.find(element_name);
    if (name_start == std::string_view::npos) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: malformed XML opening tag");
    }
    std::size_t cursor =
        name_start + element_name.size();
    std::map<std::string, std::string> attributes;

    while (cursor < open_tag.size()) {
        while (cursor < open_tag.size() &&
               is_xml_space(open_tag[cursor])) {
            ++cursor;
        }
        if (cursor >= open_tag.size() ||
            open_tag[cursor] == '>' ||
            open_tag[cursor] == '/') {
            break;
        }

        const std::size_t key_begin = cursor;
        while (cursor < open_tag.size() &&
               !is_xml_space(open_tag[cursor]) &&
               open_tag[cursor] != '=' &&
               open_tag[cursor] != '>' &&
               open_tag[cursor] != '/') {
            ++cursor;
        }
        if (cursor == key_begin) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: malformed XML attribute name");
        }
        const std::string key{
            open_tag.substr(
                key_begin, cursor - key_begin)};

        while (cursor < open_tag.size() &&
               is_xml_space(open_tag[cursor])) {
            ++cursor;
        }
        if (cursor >= open_tag.size() ||
            open_tag[cursor] != '=') {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: XML attribute is missing '='");
        }
        ++cursor;
        while (cursor < open_tag.size() &&
               is_xml_space(open_tag[cursor])) {
            ++cursor;
        }
        if (cursor >= open_tag.size() ||
            (open_tag[cursor] != '"' &&
             open_tag[cursor] != '\'')) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: XML attribute value must be quoted");
        }
        const char quote = open_tag[cursor++];
        const std::size_t value_begin = cursor;
        const auto value_end =
            open_tag.find(quote, cursor);
        if (value_end == std::string_view::npos) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: unterminated XML attribute value");
        }
        const std::string value =
            xml_unescape(
                open_tag.substr(
                    value_begin,
                    value_end - value_begin));
        cursor = value_end + 1U;

        if (!attributes.emplace(key, value).second) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: duplicate XML attribute");
        }
    }
    return attributes;
}

[[nodiscard]] inline bool tag_name_matches(
    std::string_view xml,
    std::size_t offset,
    std::string_view name) {
    if (offset + 1U + name.size() >
        xml.size() ||
        xml[offset] != '<' ||
        xml.substr(offset + 1U, name.size()) !=
            name) {
        return false;
    }
    const std::size_t after =
        offset + 1U + name.size();
    return after < xml.size() &&
           (is_xml_space(xml[after]) ||
            xml[after] == '>' ||
            xml[after] == '/');
}

[[nodiscard]] inline std::optional<XmlElement>
find_element(std::string_view xml,
             std::string_view name,
             std::size_t start = 0U) {
    std::size_t open = start;
    while (true) {
        open = xml.find('<', open);
        if (open == std::string_view::npos) {
            return std::nullopt;
        }
        if (tag_name_matches(xml, open, name)) {
            break;
        }
        ++open;
    }

    const auto open_end =
        xml.find('>', open + 1U);
    if (open_end == std::string_view::npos) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: unterminated XML opening tag");
    }
    const auto open_tag =
        xml.substr(
            open, open_end - open + 1U);
    std::size_t tail = open_end;
    while (tail > open &&
           is_xml_space(xml[tail - 1U])) {
        --tail;
    }
    const bool self_closing =
        tail > open && xml[tail - 1U] == '/';
    const auto attributes =
        parse_attributes(open_tag, name);

    if (self_closing) {
        return XmlElement{
            attributes, {}, open_end + 1U, true};
    }

    const std::string closing =
        "</" + std::string(name) + ">";
    const auto close =
        xml.find(closing, open_end + 1U);
    if (close == std::string_view::npos) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: missing XML closing tag");
    }
    return XmlElement{
        attributes,
        xml.substr(
            open_end + 1U,
            close - (open_end + 1U)),
        close + closing.size(),
        false};
}

[[nodiscard]] inline XmlElement require_unique_element(
    std::string_view xml,
    std::string_view name,
    const char* message) {
    const auto first =
        find_element(xml, name);
    if (!first.has_value()) {
        throw std::invalid_argument(message);
    }
    if (find_element(
            xml, name, first->next_offset)
            .has_value()) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: duplicate XML element where one is required");
    }
    return *first;
}

[[nodiscard]] inline const std::string& require_attribute(
    const XmlElement& element,
    std::string_view name,
    const char* message) {
    const auto found =
        element.attributes.find(
            std::string(name));
    if (found == element.attributes.end()) {
        throw std::invalid_argument(message);
    }
    return found->second;
}

template <typename Integer>
[[nodiscard]] inline Integer parse_integer(
    std::string_view text,
    const char* message) {
    Integer value{};
    const char* first = text.data();
    const char* last = text.data() + text.size();
    const auto [end, error] =
        std::from_chars(first, last, value);
    if (error != std::errc{} || end != last) {
        throw std::invalid_argument(message);
    }
    return value;
}

[[nodiscard]] inline std::size_t parse_size(
    std::string_view text,
    const char* message) {
    const auto value =
        parse_integer<std::uint64_t>(
            text, message);
    if (value >
        static_cast<std::uint64_t>(
            std::numeric_limits<std::size_t>::max())) {
        throw std::length_error(message);
    }
    return static_cast<std::size_t>(value);
}

template <typename Integer>
[[nodiscard]] inline std::vector<Integer>
parse_integer_values(std::string_view body,
                     const char* message) {
    std::istringstream input{std::string(body)};
    std::vector<Integer> values;
    std::string token;
    while (input >> token) {
        values.push_back(
            parse_integer<Integer>(
                token, message));
    }
    return values;
}

[[nodiscard]] inline std::vector<double>
parse_double_values(std::string_view body,
                    const char* message) {
    std::istringstream input{std::string(body)};
    std::vector<double> values;
    double value = 0.0;
    while (input >> value) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(message);
        }
        values.push_back(value);
    }
    if (!input.eof()) {
        throw std::invalid_argument(message);
    }
    return values;
}

[[nodiscard]] inline std::vector<XmlElement>
data_arrays(std::string_view section) {
    std::vector<XmlElement> arrays;
    std::size_t cursor = 0U;
    while (true) {
        const auto element =
            find_element(
                section, "DataArray", cursor);
        if (!element.has_value()) break;
        if (element->self_closing) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: ASCII DataArray cannot be self-closing");
        }
        if (element->body.find('<') !=
            std::string_view::npos) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: nested DataArray content is unsupported in the minimal ASCII baseline");
        }
        arrays.push_back(*element);
        cursor = element->next_offset;
    }
    return arrays;
}

inline void require_ascii_data_array(
    const XmlElement& array) {
    const auto& format =
        require_attribute(
            array, "format",
            "mpmc::mesh::import_vtu_ascii: DataArray format is required");
    if (format != "ascii") {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: only DataArray format=\"ascii\" is supported");
    }
}

[[nodiscard]] inline std::size_t component_count(
    const XmlElement& array) {
    const auto found =
        array.attributes.find(
            "NumberOfComponents");
    if (found == array.attributes.end()) {
        return 1U;
    }
    const std::size_t count =
        parse_size(
            found->second,
            "mpmc::mesh::import_vtu_ascii: invalid NumberOfComponents");
    if (count == 0U) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: NumberOfComponents must be positive");
    }
    return count;
}

[[nodiscard]] inline FieldSourceKind
parse_source_kind(std::string_view value) {
    if (value == "file_import") {
        return FieldSourceKind::file_import;
    }
    if (value == "user_supplied") {
        return FieldSourceKind::user_supplied;
    }
    if (value == "generated") {
        return FieldSourceKind::generated;
    }
    if (value == "synthetic_test") {
        return FieldSourceKind::synthetic_test;
    }
    throw std::invalid_argument(
        "mpmc::mesh::import_vtu_ascii: invalid mpmc_source_kind");
}

[[nodiscard]] inline std::string source_kind_text(
    FieldSourceKind kind) {
    switch (kind) {
    case FieldSourceKind::file_import:
        return "file_import";
    case FieldSourceKind::user_supplied:
        return "user_supplied";
    case FieldSourceKind::generated:
        return "generated";
    case FieldSourceKind::synthetic_test:
        return "synthetic_test";
    case FieldSourceKind::unspecified:
        break;
    }
    throw std::invalid_argument(
        "mpmc::mesh::export_vtu_ascii: invalid field source kind");
}

[[nodiscard]] inline DenseFieldMetadata field_metadata(
    const XmlElement& array,
    std::string_view name,
    EntityKind location) {
    const auto lookup =
        [&](std::string_view key)
            -> std::optional<std::string> {
        const auto found =
            array.attributes.find(
                std::string(key));
        if (found == array.attributes.end()) {
            return std::nullopt;
        }
        return found->second;
    };

    DenseFieldMetadata metadata;
    metadata.id = std::string(name);
    metadata.unit =
        lookup("mpmc_unit").value_or("1");
    const auto source_kind =
        lookup("mpmc_source_kind");
    metadata.source.kind =
        source_kind.has_value()
            ? parse_source_kind(*source_kind)
            : FieldSourceKind::file_import;
    metadata.source.reference =
        lookup("mpmc_source_reference")
            .value_or(
                "VTK XML UnstructuredGrid");
    metadata.source.revision =
        lookup("mpmc_source_revision")
            .value_or("VTU ASCII baseline");
    metadata.source.locator =
        lookup("mpmc_source_locator")
            .value_or(
                location == EntityKind::vertex
                    ? "PointData/" +
                          std::string(name)
                    : "CellData/" +
                          std::string(name));
    return metadata;
}

[[nodiscard]] inline ParsedField parse_field(
    const XmlElement& array,
    EntityKind location,
    std::size_t entity_count) {
    require_ascii_data_array(array);
    const auto& type =
        require_attribute(
            array, "type",
            "mpmc::mesh::import_vtu_ascii: field DataArray type is required");
    if (type != "Float32" &&
        type != "Float64") {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: point/cell fields must use Float32 or Float64 in the minimal baseline");
    }
    const auto& name =
        require_attribute(
            array, "Name",
            "mpmc::mesh::import_vtu_ascii: point/cell field Name is required");
    if (name.empty()) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: field Name cannot be empty");
    }
    const std::size_t components =
        component_count(array);
    if (entity_count != 0U &&
        components >
            std::numeric_limits<std::size_t>::max() /
                entity_count) {
        throw std::length_error(
            "mpmc::mesh::import_vtu_ascii: field value count overflow");
    }
    auto values =
        parse_double_values(
            array.body,
            "mpmc::mesh::import_vtu_ascii: invalid field values");
    if (values.size() !=
        entity_count * components) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: field value count does not match Piece entity count and NumberOfComponents");
    }
    return ParsedField{
        name,
        components,
        std::move(values),
        field_metadata(
            array, name, location)};
}

[[nodiscard]] inline EdgeKey edge_key(
    std::size_t first,
    std::size_t second) {
    if (first == second) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: cell edge has repeated point index");
    }
    return first < second
               ? EdgeKey{first, second}
               : EdgeKey{second, first};
}

[[nodiscard]] inline double cross(
    Coordinate2D a,
    Coordinate2D b,
    Coordinate2D c) {
    return (b.x_m - a.x_m) *
               (c.y_m - a.y_m) -
           (b.y_m - a.y_m) *
               (c.x_m - a.x_m);
}

[[nodiscard]] inline bool segments_intersect_strictly(
    Coordinate2D a,
    Coordinate2D b,
    Coordinate2D c,
    Coordinate2D d,
    double tolerance) {
    const double c1 = cross(a, b, c);
    const double c2 = cross(a, b, d);
    const double c3 = cross(c, d, a);
    const double c4 = cross(c, d, b);
    if (std::abs(c1) <= tolerance ||
        std::abs(c2) <= tolerance ||
        std::abs(c3) <= tolerance ||
        std::abs(c4) <= tolerance) {
        return false;
    }
    return (c1 > 0.0) != (c2 > 0.0) &&
           (c3 > 0.0) != (c4 > 0.0);
}

[[nodiscard]] inline CellMetric polygon_metric(
    const std::vector<Coordinate2D>& points) {
    if (points.size() != 3U &&
        points.size() != 4U) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: only linear triangle and quad cells are supported");
    }

    double scale = 1.0;
    for (const auto point : points) {
        scale = std::max(
            scale,
            std::max(
                std::abs(point.x_m),
                std::abs(point.y_m)));
    }
    const double tolerance =
        256.0 *
        std::numeric_limits<double>::epsilon() *
        scale * scale;

    if (points.size() == 4U) {
        if (segments_intersect_strictly(
                points[0], points[1],
                points[2], points[3],
                tolerance) ||
            segments_intersect_strictly(
                points[1], points[2],
                points[3], points[0],
                tolerance)) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: self-intersecting quad is unsupported");
        }
        double sign = 0.0;
        for (std::size_t i = 0U; i < 4U; ++i) {
            const double turn =
                cross(
                    points[i],
                    points[(i + 1U) % 4U],
                    points[(i + 2U) % 4U]);
            if (std::abs(turn) <= tolerance) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_vtu_ascii: degenerate quad corner");
            }
            if (sign == 0.0) sign = turn;
            else if ((turn > 0.0) !=
                     (sign > 0.0)) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_vtu_ascii: concave quad is unsupported by the minimal baseline");
            }
        }
    }

    double twice_area = 0.0;
    double centroid_x = 0.0;
    double centroid_y = 0.0;
    for (std::size_t i = 0U;
         i < points.size();
         ++i) {
        const auto& a = points[i];
        const auto& b =
            points[(i + 1U) %
                   points.size()];
        const double term =
            a.x_m * b.y_m -
            b.x_m * a.y_m;
        twice_area += term;
        centroid_x +=
            (a.x_m + b.x_m) * term;
        centroid_y +=
            (a.y_m + b.y_m) * term;
    }
    if (!std::isfinite(twice_area) ||
        std::abs(twice_area) <= tolerance) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: degenerate cell area");
    }
    const Coordinate2D centroid{
        centroid_x / (3.0 * twice_area),
        centroid_y / (3.0 * twice_area)};
    const double area =
        0.5 * std::abs(twice_area);
    if (!std::isfinite(centroid.x_m) ||
        !std::isfinite(centroid.y_m) ||
        !std::isfinite(area) ||
        area <= 0.0) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: invalid cell metric");
    }
    return CellMetric{centroid, area};
}

[[nodiscard]] inline LocalIndex local_index(
    std::size_t value,
    const char* message) {
    if (value >
        static_cast<std::size_t>(
            std::numeric_limits<
                LocalIndex::value_type>::max())) {
        throw std::length_error(message);
    }
    return LocalIndex{
        static_cast<
            LocalIndex::value_type>(value)};
}

[[nodiscard]] inline CsrAdjacency::Offset offset_value(
    std::size_t value,
    const char* message) {
    if (value >
        static_cast<std::size_t>(
            std::numeric_limits<
                CsrAdjacency::Offset>::max())) {
        throw std::length_error(message);
    }
    return static_cast<
        CsrAdjacency::Offset>(value);
}

[[nodiscard]] inline std::string attribute(
    const DenseFieldMetadata& metadata,
    std::string_view key) {
    if (key == "mpmc_unit") {
        return xml_escape(metadata.unit);
    }
    if (key == "mpmc_source_kind") {
        return xml_escape(
            source_kind_text(
                metadata.source.kind));
    }
    if (key == "mpmc_source_reference") {
        return xml_escape(
            metadata.source.reference);
    }
    if (key == "mpmc_source_revision") {
        return xml_escape(
            metadata.source.revision);
    }
    if (key == "mpmc_source_locator") {
        return xml_escape(
            metadata.source.locator);
    }
    throw std::logic_error(
        "mpmc::mesh::export_vtu_ascii: unknown metadata attribute");
}

inline void validate_export_fields(
    std::span<const DenseFieldSnapshot> fields,
    EntityKind expected_location,
    std::size_t expected_count) {
    std::set<std::string> names;
    for (const auto& field : fields) {
        if (field.location() != expected_location ||
            field.entity_count() != expected_count) {
            throw std::invalid_argument(
                "mpmc::mesh::export_vtu_ascii: field location/entity count does not match association");
        }
        if (!names.insert(
                field.metadata().id).second) {
            throw std::invalid_argument(
                "mpmc::mesh::export_vtu_ascii: duplicate field ID in one association");
        }
        if (expected_location ==
                EntityKind::cell &&
            field.metadata().id ==
                global_cell_id_name) {
            throw std::invalid_argument(
                "mpmc::mesh::export_vtu_ascii: cell field name is reserved for stable cell identity");
        }
    }
}

inline void write_field(
    std::ostringstream& output,
    const DenseFieldSnapshot& field,
    std::size_t indent) {
    const std::string spaces(indent, ' ');
    output << spaces
           << "<DataArray type=\"Float64\" Name=\""
           << xml_escape(field.metadata().id)
           << "\" NumberOfComponents=\""
           << field.component_count()
           << "\" format=\"ascii\""
           << " mpmc_unit=\""
           << attribute(
                  field.metadata(),
                  "mpmc_unit")
           << "\" mpmc_source_kind=\""
           << attribute(
                  field.metadata(),
                  "mpmc_source_kind")
           << "\" mpmc_source_reference=\""
           << attribute(
                  field.metadata(),
                  "mpmc_source_reference")
           << "\" mpmc_source_revision=\""
           << attribute(
                  field.metadata(),
                  "mpmc_source_revision")
           << "\" mpmc_source_locator=\""
           << attribute(
                  field.metadata(),
                  "mpmc_source_locator")
           << "\">\n";
    output << spaces << "  ";
    for (std::size_t i = 0U;
         i < field.values().size();
         ++i) {
        if (i != 0U) output << ' ';
        output << field.values()[i];
    }
    output << '\n'
           << spaces << "</DataArray>\n";
}

} // namespace vtu_detail

[[nodiscard]] inline VtuImportResult
import_vtu_ascii(std::string_view content) {
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
            "mpmc::mesh::import_vtu_ascii: binary/appended VTU is unsupported");
    }

    const auto vtk_file =
        require_unique_element(
            content, "VTKFile",
            "mpmc::mesh::import_vtu_ascii: VTKFile root is required");
    const auto& vtk_type =
        require_attribute(
            vtk_file, "type",
            "mpmc::mesh::import_vtu_ascii: VTKFile type is required");
    if (vtk_type != "UnstructuredGrid") {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: only VTKFile type=UnstructuredGrid is supported");
    }
    if (vtk_file.attributes.contains(
            "compressor")) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: compressed VTU is unsupported");
    }

    const auto grid =
        require_unique_element(
            vtk_file.body,
            "UnstructuredGrid",
            "mpmc::mesh::import_vtu_ascii: UnstructuredGrid element is required");
    const auto piece =
        require_unique_element(
            grid.body, "Piece",
            "mpmc::mesh::import_vtu_ascii: exactly one Piece is required");
    const std::size_t point_count =
        parse_size(
            require_attribute(
                piece,
                "NumberOfPoints",
                "mpmc::mesh::import_vtu_ascii: Piece NumberOfPoints is required"),
            "mpmc::mesh::import_vtu_ascii: invalid NumberOfPoints");
    const std::size_t cell_count =
        parse_size(
            require_attribute(
                piece,
                "NumberOfCells",
                "mpmc::mesh::import_vtu_ascii: Piece NumberOfCells is required"),
            "mpmc::mesh::import_vtu_ascii: invalid NumberOfCells");
    if (point_count == 0U ||
        cell_count == 0U) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: Piece must contain points and cells");
    }

    const auto points_section =
        require_unique_element(
            piece.body, "Points",
            "mpmc::mesh::import_vtu_ascii: Points section is required");
    const auto point_arrays =
        data_arrays(points_section.body);
    if (point_arrays.size() != 1U) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: Points must contain exactly one DataArray");
    }
    const auto& point_array =
        point_arrays.front();
    require_ascii_data_array(point_array);
    const auto& point_type =
        require_attribute(
            point_array, "type",
            "mpmc::mesh::import_vtu_ascii: Points DataArray type is required");
    if (point_type != "Float32" &&
        point_type != "Float64") {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: Points must use Float32 or Float64");
    }
    if (component_count(point_array) != 3U) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: Points must have NumberOfComponents=3");
    }
    const auto point_values =
        parse_double_values(
            point_array.body,
            "mpmc::mesh::import_vtu_ascii: invalid point coordinates");
    if (point_values.size() !=
        point_count * 3U) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: point coordinate count mismatch");
    }

    std::vector<Coordinate2D> coordinates;
    coordinates.reserve(point_count);
    for (std::size_t point = 0U;
         point < point_count;
         ++point) {
        const double x =
            point_values[point * 3U];
        const double y =
            point_values[point * 3U + 1U];
        const double z =
            point_values[point * 3U + 2U];
        const double scale =
            std::max(
                1.0,
                std::max(
                    std::abs(x),
                    std::abs(y)));
        const double tolerance =
            256.0 *
            std::numeric_limits<double>::epsilon() *
            scale;
        if (std::abs(z) > tolerance) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: minimal Geometry2D baseline requires all points in the XY plane");
        }
        coordinates.push_back(
            Coordinate2D{x, y});
    }

    const auto cells_section =
        require_unique_element(
            piece.body, "Cells",
            "mpmc::mesh::import_vtu_ascii: Cells section is required");
    const auto cell_arrays =
        data_arrays(cells_section.body);
    if (cell_arrays.size() != 3U) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: Cells must contain connectivity, offsets and types DataArrays");
    }
    std::optional<XmlElement> connectivity_array;
    std::optional<XmlElement> offsets_array;
    std::optional<XmlElement> types_array;
    for (const auto& array : cell_arrays) {
        require_ascii_data_array(array);
        const auto& name =
            require_attribute(
                array, "Name",
                "mpmc::mesh::import_vtu_ascii: Cells DataArray Name is required");
        if (name == "connectivity") {
            if (connectivity_array.has_value()) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_vtu_ascii: duplicate connectivity array");
            }
            connectivity_array = array;
        } else if (name == "offsets") {
            if (offsets_array.has_value()) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_vtu_ascii: duplicate offsets array");
            }
            offsets_array = array;
        } else if (name == "types") {
            if (types_array.has_value()) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_vtu_ascii: duplicate types array");
            }
            types_array = array;
        } else {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: unsupported Cells DataArray");
        }
    }
    if (!connectivity_array.has_value() ||
        !offsets_array.has_value() ||
        !types_array.has_value()) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: incomplete Cells arrays");
    }

    const auto connectivity =
        parse_integer_values<std::int64_t>(
            connectivity_array->body,
            "mpmc::mesh::import_vtu_ascii: invalid connectivity");
    const auto offsets =
        parse_integer_values<std::int64_t>(
            offsets_array->body,
            "mpmc::mesh::import_vtu_ascii: invalid offsets");
    const auto types =
        parse_integer_values<std::int64_t>(
            types_array->body,
            "mpmc::mesh::import_vtu_ascii: invalid types");
    if (offsets.size() != cell_count ||
        types.size() != cell_count) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: offsets/types count mismatch");
    }

    std::vector<std::vector<std::size_t>>
        cell_point_indices;
    cell_point_indices.reserve(cell_count);
    std::size_t previous_offset = 0U;
    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        if (offsets[cell] < 0) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: cell offset cannot be negative");
        }
        const std::size_t end =
            static_cast<std::size_t>(
                offsets[cell]);
        if (end <= previous_offset ||
            end > connectivity.size()) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: invalid monotone cell offsets");
        }
        const std::size_t width =
            end - previous_offset;
        const std::size_t expected =
            types[cell] == 5
                ? 3U
                : types[cell] == 9
                      ? 4U
                      : 0U;
        if (expected == 0U ||
            width != expected) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: only VTK_TRIANGLE(5) and VTK_QUAD(9) are supported with matching connectivity widths");
        }

        std::vector<std::size_t> points;
        points.reserve(width);
        std::set<std::size_t> unique;
        for (std::size_t i = previous_offset;
             i < end;
             ++i) {
            if (connectivity[i] < 0 ||
                static_cast<std::uint64_t>(
                    connectivity[i]) >=
                    static_cast<std::uint64_t>(
                        point_count)) {
                throw std::out_of_range(
                    "mpmc::mesh::import_vtu_ascii: connectivity point index out of range");
            }
            const auto point =
                static_cast<std::size_t>(
                    connectivity[i]);
            if (!unique.insert(point).second) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_vtu_ascii: cell contains repeated point index");
            }
            points.push_back(point);
        }
        cell_point_indices.push_back(
            std::move(points));
        previous_offset = end;
    }
    if (previous_offset != connectivity.size()) {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: final cell offset must consume all connectivity");
    }

    std::vector<std::uint64_t>
        cell_global_ids;
    std::vector<ParsedField>
        parsed_cell_fields;
    const auto cell_data =
        find_element(piece.body, "CellData");
    if (cell_data.has_value()) {
        std::set<std::string> field_names;
        for (const auto& array :
             data_arrays(cell_data->body)) {
            require_ascii_data_array(array);
            const auto& name =
                require_attribute(
                    array, "Name",
                    "mpmc::mesh::import_vtu_ascii: CellData DataArray Name is required");
            if (!field_names.insert(name).second) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_vtu_ascii: duplicate CellData Name");
            }
            if (name == global_cell_id_name) {
                const auto& type =
                    require_attribute(
                        array, "type",
                        "mpmc::mesh::import_vtu_ascii: global cell ID type is required");
                if (type != "UInt64" &&
                    type != "Int64") {
                    throw std::invalid_argument(
                        "mpmc::mesh::import_vtu_ascii: mpmc_global_cell_id must use UInt64 or Int64");
                }
                const auto values =
                    parse_integer_values<std::int64_t>(
                        array.body,
                        "mpmc::mesh::import_vtu_ascii: invalid global cell IDs");
                if (values.size() != cell_count) {
                    throw std::invalid_argument(
                        "mpmc::mesh::import_vtu_ascii: global cell ID count mismatch");
                }
                cell_global_ids.clear();
                cell_global_ids.reserve(cell_count);
                std::set<std::uint64_t> unique;
                for (const auto value : values) {
                    if (value < 0) {
                        throw std::invalid_argument(
                            "mpmc::mesh::import_vtu_ascii: global cell ID cannot be negative");
                    }
                    const auto id =
                        static_cast<std::uint64_t>(
                            value);
                    if (!unique.insert(id).second) {
                        throw std::invalid_argument(
                            "mpmc::mesh::import_vtu_ascii: global cell IDs must be unique");
                    }
                    cell_global_ids.push_back(id);
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
        cell_global_ids.reserve(cell_count);
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
        find_element(piece.body, "PointData");
    if (point_data.has_value()) {
        std::set<std::string> field_names;
        for (const auto& array :
             data_arrays(point_data->body)) {
            const auto& name =
                require_attribute(
                    array, "Name",
                    "mpmc::mesh::import_vtu_ascii: PointData DataArray Name is required");
            if (!field_names.insert(name).second) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_vtu_ascii: duplicate PointData Name");
            }
            parsed_point_fields.push_back(
                parse_field(
                    array,
                    EntityKind::vertex,
                    point_count));
        }
    }

    std::vector<CellMetric> cell_metrics;
    cell_metrics.reserve(cell_count);
    std::map<EdgeKey, FaceBuild> faces_by_edge;
    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        std::vector<Coordinate2D> polygon;
        polygon.reserve(
            cell_point_indices[cell].size());
        for (const auto point :
             cell_point_indices[cell]) {
            polygon.push_back(
                coordinates[point]);
        }
        cell_metrics.push_back(
            polygon_metric(polygon));

        const auto& points =
            cell_point_indices[cell];
        for (std::size_t edge = 0U;
             edge < points.size();
             ++edge) {
            const auto key =
                edge_key(
                    points[edge],
                    points[
                        (edge + 1U) %
                        points.size()]);
            auto [found, inserted] =
                faces_by_edge.emplace(
                    key,
                    FaceBuild{key, {}});
            found->second.adjacent_cells
                .push_back(cell);
            if (found->second
                    .adjacent_cells.size() >
                2U) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_vtu_ascii: non-manifold edge belongs to more than two cells");
            }
            (void)inserted;
        }
    }

    std::vector<FaceBuild> faces;
    faces.reserve(faces_by_edge.size());
    for (auto& [key, face] :
         faces_by_edge) {
        std::sort(
            face.adjacent_cells.begin(),
            face.adjacent_cells.end());
        faces.push_back(std::move(face));
        (void)key;
    }

    Topology::EntityIds ids;
    ids.vertices.reserve(point_count);
    for (std::size_t point = 0U;
         point < point_count;
         ++point) {
        ids.vertices.emplace_back(
            static_cast<std::uint64_t>(
                point) +
            1U);
    }
    ids.faces.reserve(faces.size());
    for (std::size_t face = 0U;
         face < faces.size();
         ++face) {
        ids.faces.emplace_back(
            static_cast<std::uint64_t>(
                face) +
            1U);
    }
    ids.cells.reserve(cell_count);
    for (const auto id : cell_global_ids) {
        ids.cells.emplace_back(id);
    }

    std::map<EdgeKey, std::size_t>
        face_local_by_edge;
    for (std::size_t face = 0U;
         face < faces.size();
         ++face) {
        face_local_by_edge.emplace(
            faces[face].key, face);
    }

    std::vector<CsrAdjacency::Offset>
        cell_vertex_offsets{0U};
    std::vector<LocalIndex> cell_vertices;
    std::vector<CsrAdjacency::Offset>
        cell_face_offsets{0U};
    std::vector<LocalIndex> cell_faces;
    cell_vertex_offsets.reserve(
        cell_count + 1U);
    cell_face_offsets.reserve(
        cell_count + 1U);

    for (const auto& points :
         cell_point_indices) {
        for (const auto point : points) {
            cell_vertices.push_back(
                local_index(
                    point,
                    "mpmc::mesh::import_vtu_ascii: vertex local index overflow"));
        }
        for (std::size_t edge = 0U;
             edge < points.size();
             ++edge) {
            const auto key =
                edge_key(
                    points[edge],
                    points[
                        (edge + 1U) %
                        points.size()]);
            const auto found =
                face_local_by_edge.find(key);
            if (found ==
                face_local_by_edge.end()) {
                throw std::logic_error(
                    "mpmc::mesh::import_vtu_ascii: internal face lookup failed");
            }
            cell_faces.push_back(
                local_index(
                    found->second,
                    "mpmc::mesh::import_vtu_ascii: face local index overflow"));
        }
        cell_vertex_offsets.push_back(
            offset_value(
                cell_vertices.size(),
                "mpmc::mesh::import_vtu_ascii: cell->vertex CSR overflow"));
        cell_face_offsets.push_back(
            offset_value(
                cell_faces.size(),
                "mpmc::mesh::import_vtu_ascii: cell->face CSR overflow"));
    }

    std::vector<CsrAdjacency::Offset>
        face_vertex_offsets{0U};
    std::vector<LocalIndex> face_vertices;
    std::vector<CsrAdjacency::Offset>
        face_cell_offsets{0U};
    std::vector<LocalIndex> face_cells;
    face_vertex_offsets.reserve(
        faces.size() + 1U);
    face_cell_offsets.reserve(
        faces.size() + 1U);

    for (const auto& face : faces) {
        face_vertices.push_back(
            local_index(
                face.key.first,
                "mpmc::mesh::import_vtu_ascii: face vertex local index overflow"));
        face_vertices.push_back(
            local_index(
                face.key.second,
                "mpmc::mesh::import_vtu_ascii: face vertex local index overflow"));
        for (const auto cell :
             face.adjacent_cells) {
            face_cells.push_back(
                local_index(
                    cell,
                    "mpmc::mesh::import_vtu_ascii: face cell local index overflow"));
        }
        face_vertex_offsets.push_back(
            offset_value(
                face_vertices.size(),
                "mpmc::mesh::import_vtu_ascii: face->vertex CSR overflow"));
        face_cell_offsets.push_back(
            offset_value(
                face_cells.size(),
                "mpmc::mesh::import_vtu_ascii: face->cell CSR overflow"));
    }

    std::vector<CsrAdjacency> relations;
    relations.reserve(4U);
    relations.emplace_back(
        EntityKind::cell,
        EntityKind::vertex,
        point_count,
        std::move(cell_vertex_offsets),
        std::move(cell_vertices));
    relations.emplace_back(
        EntityKind::cell,
        EntityKind::face,
        faces.size(),
        std::move(cell_face_offsets),
        std::move(cell_faces));
    relations.emplace_back(
        EntityKind::face,
        EntityKind::vertex,
        point_count,
        std::move(face_vertex_offsets),
        std::move(face_vertices));
    relations.emplace_back(
        EntityKind::face,
        EntityKind::cell,
        cell_count,
        std::move(face_cell_offsets),
        std::move(face_cells));

    Topology topology{
        std::move(ids),
        std::move(relations)};

    std::vector<Coordinate2D> cell_centroids;
    std::vector<double> cell_areas;
    cell_centroids.reserve(cell_count);
    cell_areas.reserve(cell_count);
    for (const auto metric : cell_metrics) {
        cell_centroids.push_back(
            metric.centroid);
        cell_areas.push_back(
            metric.area);
    }

    std::vector<Coordinate2D> face_centroids;
    std::vector<double> face_lengths;
    std::vector<LocalIndex> face_owners;
    std::vector<UnitVector2D> face_normals;
    face_centroids.reserve(faces.size());
    face_lengths.reserve(faces.size());
    face_owners.reserve(faces.size());
    face_normals.reserve(faces.size());

    for (const auto& face : faces) {
        const auto a =
            coordinates[face.key.first];
        const auto b =
            coordinates[face.key.second];
        const double dx = b.x_m - a.x_m;
        const double dy = b.y_m - a.y_m;
        const double length =
            std::hypot(dx, dy);
        if (!std::isfinite(length) ||
            length <= 0.0) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: face length must be finite and positive");
        }
        const Coordinate2D centroid{
            std::midpoint(a.x_m, b.x_m),
            std::midpoint(a.y_m, b.y_m)};
        const std::size_t owner =
            face.adjacent_cells.front();
        const auto owner_centroid =
            cell_metrics[owner].centroid;
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
        if (dot == 0.0) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: cannot orient face normal");
        }
        if (dot < 0.0) {
            normal.x = -normal.x;
            normal.y = -normal.y;
        }
        face_centroids.push_back(centroid);
        face_lengths.push_back(length);
        face_owners.push_back(
            local_index(
                owner,
                "mpmc::mesh::import_vtu_ascii: face owner local index overflow"));
        face_normals.push_back(normal);
    }

    Geometry2D geometry{
        std::move(coordinates),
        std::move(cell_centroids),
        std::move(cell_areas),
        std::move(face_centroids),
        std::move(face_lengths),
        std::move(face_owners),
        std::move(face_normals)};

    std::vector<DenseFieldSnapshot> point_fields;
    point_fields.reserve(
        parsed_point_fields.size());
    for (auto& field :
         parsed_point_fields) {
        point_fields.push_back(
            DenseFieldSnapshot::create(
                topology,
                EntityKind::vertex,
                field.component_count,
                std::move(field.values),
                std::move(field.metadata)));
    }
    std::vector<DenseFieldSnapshot> cell_fields;
    cell_fields.reserve(
        parsed_cell_fields.size());
    for (auto& field :
         parsed_cell_fields) {
        cell_fields.push_back(
            DenseFieldSnapshot::create(
                topology,
                EntityKind::cell,
                field.component_count,
                std::move(field.values),
                std::move(field.metadata)));
    }

    return VtuImportResult{
        std::move(topology),
        std::move(geometry),
        std::move(point_fields),
        std::move(cell_fields)};
}

[[nodiscard]] inline std::string
export_vtu_ascii(const VtuImportResult& mesh) {
    using namespace vtu_detail;

    const auto& topology = mesh.topology;
    const auto& geometry = mesh.geometry;
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
            "mpmc::mesh::export_vtu_ascii: only nonempty 2D vertex/face/cell meshes are supported");
    }
    if (geometry.vertex_count() !=
            point_count ||
        geometry.cell_count() !=
            cell_count) {
        throw std::invalid_argument(
            "mpmc::mesh::export_vtu_ascii: topology and geometry counts are not aligned");
    }
    if (!topology.has_relation(
            EntityKind::cell,
            EntityKind::vertex)) {
        throw std::invalid_argument(
            "mpmc::mesh::export_vtu_ascii: cell->vertex relation is required");
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
    std::size_t connectivity_count = 0U;
    offsets.reserve(cell_count);
    types.reserve(cell_count);
    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        const auto vertices =
            cell_vertices.adjacent(
                local_index(
                    cell,
                    "mpmc::mesh::export_vtu_ascii: cell local index overflow"));
        if (vertices.size() != 3U &&
            vertices.size() != 4U) {
            throw std::invalid_argument(
                "mpmc::mesh::export_vtu_ascii: only linear triangle/quad cells are supported");
        }
        if (connectivity_count >
            std::numeric_limits<std::size_t>::max() -
                vertices.size()) {
            throw std::length_error(
                "mpmc::mesh::export_vtu_ascii: connectivity size overflow");
        }
        connectivity_count +=
            vertices.size();
        offsets.push_back(
            connectivity_count);
        types.push_back(
            vertices.size() == 3U
                ? 5
                : 9);
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
        write_field(output, field, 8U);
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
    for (std::size_t i = 0U;
         i < cell_ids.size();
         ++i) {
        if (i != 0U) output << ' ';
        output << cell_ids[i].value();
    }
    output << "\n        </DataArray>\n";
    for (const auto& field :
         mesh.cell_fields) {
        write_field(output, field, 8U);
    }
    output << "      </CellData>\n";

    output << "      <Points>\n"
           << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n"
           << "          ";
    for (std::size_t point = 0U;
         point < point_count;
         ++point) {
        if (point != 0U) output << ' ';
        const auto coordinate =
            geometry.vertex_coordinate_m(
                local_index(
                    point,
                    "mpmc::mesh::export_vtu_ascii: point local index overflow"));
        output << coordinate.x_m << ' '
               << coordinate.y_m
               << " 0";
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
                    "mpmc::mesh::export_vtu_ascii: cell local index overflow"));
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
    for (std::size_t i = 0U;
         i < offsets.size();
         ++i) {
        if (i != 0U) output << ' ';
        output << offsets[i];
    }
    output << "\n        </DataArray>\n"
           << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n"
           << "          ";
    for (std::size_t i = 0U;
         i < types.size();
         ++i) {
        if (i != 0U) output << ' ';
        output << types[i];
    }
    output << "\n        </DataArray>\n"
           << "      </Cells>\n"
           << "    </Piece>\n"
           << "  </UnstructuredGrid>\n"
           << "</VTKFile>\n";

    return output.str();
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_VTU_HPP
