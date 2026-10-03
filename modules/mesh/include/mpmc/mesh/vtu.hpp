#ifndef MPMC_MESH_VTU_HPP
#define MPMC_MESH_VTU_HPP

#include <mpmc/mesh/csr_adjacency.hpp>
#include <mpmc/mesh/dense_field.hpp>
#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/geometry_2d.hpp>
#include <mpmc/mesh/linear_cell_mesh_2d.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <initializer_list>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <span>
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
inline constexpr std::string_view global_vertex_id_name =
    "mpmc_global_vertex_id";

struct XmlElement {
    std::map<std::string, std::string> attributes;
    std::string_view body;
    std::size_t next_offset;
    bool self_closing;
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
        auto element =
            find_element(
                section, "DataArray", cursor);
        if (!element.has_value()) break;
        if (element->self_closing) {
            throw std::invalid_argument(
                "mpmc::mesh::import_vtu_ascii: ASCII DataArray cannot be self-closing");
        }
        const auto nested = element->body.find('<');
        if (nested != std::string_view::npos) {
            // Official VTK ASCII writers append this derived magnitude-range
            // cache, including to Points. It is not array data or a physical
            // unit. Accept only its documented shape; do not silently discard
            // arbitrary InformationKey entries (which can carry semantics).
            constexpr const char* invalid =
                "mpmc::mesh::import_vtu_ascii: unsupported or malformed DataArray InformationKey";
            const auto metadata = element->body.substr(nested);
            if (!tag_name_matches(metadata, 0U, "InformationKey")) {
                throw std::invalid_argument(invalid);
            }
            const auto key = find_element(metadata, "InformationKey");
            if (!key || key->self_closing || key->attributes.size() != 3U ||
                require_attribute(*key, "name", invalid) != "L2_NORM_RANGE" ||
                require_attribute(*key, "location", invalid) != "vtkDataArray" ||
                require_attribute(*key, "length", invalid) != "2" ||
                !std::all_of(metadata.begin() + static_cast<std::ptrdiff_t>(key->next_offset),
                             metadata.end(), is_xml_space)) {
                throw std::invalid_argument(invalid);
            }
            std::size_t position = 0U;
            for (const auto index : {"0", "1"}) {
                while (position < key->body.size() && is_xml_space(key->body[position])) ++position;
                if (!tag_name_matches(key->body, position, "Value")) {
                    throw std::invalid_argument(invalid);
                }
                const auto value = find_element(key->body, "Value", position);
                if (!value || value->self_closing || value->attributes.size() != 1U ||
                    require_attribute(*value, "index", invalid) != index ||
                    parse_double_values(value->body, invalid).size() != 1U) {
                    throw std::invalid_argument(invalid);
                }
                position = value->next_offset;
            }
            if (!std::all_of(key->body.begin() + static_cast<std::ptrdiff_t>(position),
                             key->body.end(), is_xml_space)) {
                throw std::invalid_argument(invalid);
            }
            element->body = element->body.substr(0U, nested);
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

// Identity is integer metadata, never a floating-point scientific field.
[[nodiscard]] inline std::vector<GlobalEntityId> parse_vertex_ids(
    const XmlElement& array, std::size_t count) {
    require_ascii_data_array(array);
    if (component_count(array) != 1U) {
        throw std::invalid_argument("mpmc::mesh::import_vtu_ascii: vertex IDs require one component");
    }
    const auto tuples = array.attributes.find("NumberOfTuples");
    if (tuples != array.attributes.end() &&
        parse_size(tuples->second, "invalid vertex ID NumberOfTuples") != count) {
        throw std::invalid_argument("mpmc::mesh::import_vtu_ascii: vertex ID tuple count mismatch");
    }
    const auto& type = require_attribute(array, "type", "vertex ID type is required");
    std::vector<std::uint64_t> values;
    if (type == "UInt64") {
        values = parse_integer_values<std::uint64_t>(array.body, "invalid UInt64 vertex IDs");
    } else if (type == "Int64") {
        for (const auto value : parse_integer_values<std::int64_t>(array.body, "invalid Int64 vertex IDs")) {
            if (value < 0) throw std::invalid_argument("vertex ID cannot be negative");
            values.push_back(static_cast<std::uint64_t>(value));
        }
    } else {
        throw std::invalid_argument("vertex IDs must use UInt64 or Int64");
    }
    if (values.size() != count) throw std::invalid_argument("vertex ID count mismatch");
    std::set<std::uint64_t> unique;
    std::vector<GlobalEntityId> result;
    result.reserve(count);
    for (const auto value : values) {
        if (!unique.insert(value).second) throw std::invalid_argument("vertex IDs must be unique");
        result.emplace_back(value);
    }
    return result;
}

inline void write_vertex_ids(std::ostringstream& output, const Topology& topology) {
    output << "        <DataArray type=\"UInt64\" Name=\"" << global_vertex_id_name
           << "\" NumberOfComponents=\"1\" format=\"ascii\">\n          ";
    const auto ids = topology.global_ids(EntityKind::vertex);
    for (std::size_t i = 0; i < ids.size(); ++i) {
        if (i != 0U) output << ' ';
        output << ids[i].value();
    }
    output << "\n        </DataArray>\n";
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
    if (name == global_vertex_id_name || name == global_cell_id_name) {
        throw std::invalid_argument("identity array name is reserved for its entity association");
    }
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
        if (field.metadata().id == global_cell_id_name ||
            field.metadata().id == global_vertex_id_name) {
            throw std::invalid_argument(
                "mpmc::mesh::export_vtu_ascii: field name is reserved for stable entity identity");
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

    const auto& connectivity_type =
        require_attribute(
            *connectivity_array, "type",
            "mpmc::mesh::import_vtu_ascii: connectivity type is required");
    const auto& offsets_type =
        require_attribute(
            *offsets_array, "type",
            "mpmc::mesh::import_vtu_ascii: offsets type is required");
    const auto& types_type =
        require_attribute(
            *types_array, "type",
            "mpmc::mesh::import_vtu_ascii: types type is required");
    if ((connectivity_type != "Int32" &&
         connectivity_type != "Int64") ||
        (offsets_type != "Int32" &&
         offsets_type != "Int64") ||
        types_type != "UInt8") {
        throw std::invalid_argument(
            "mpmc::mesh::import_vtu_ascii: Cells requires Int32/Int64 connectivity+offsets and UInt8 types");
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
                cell_global_ids.clear();
                cell_global_ids.reserve(cell_count);
                std::set<std::uint64_t> unique;
                if (type == "UInt64") {
                    const auto values =
                        parse_integer_values<std::uint64_t>(
                            array.body,
                            "mpmc::mesh::import_vtu_ascii: invalid UInt64 global cell IDs");
                    if (values.size() != cell_count) {
                        throw std::invalid_argument(
                            "mpmc::mesh::import_vtu_ascii: global cell ID count mismatch");
                    }
                    for (const auto id : values) {
                        if (!unique.insert(id).second) {
                            throw std::invalid_argument(
                                "mpmc::mesh::import_vtu_ascii: global cell IDs must be unique");
                        }
                        cell_global_ids.push_back(id);
                    }
                } else {
                    const auto values =
                        parse_integer_values<std::int64_t>(
                            array.body,
                            "mpmc::mesh::import_vtu_ascii: invalid Int64 global cell IDs");
                    if (values.size() != cell_count) {
                        throw std::invalid_argument(
                            "mpmc::mesh::import_vtu_ascii: global cell ID count mismatch");
                    }
                    for (const auto value : values) {
                        if (value < 0) {
                            throw std::invalid_argument(
                                "mpmc::mesh::import_vtu_ascii: global cell ID cannot be negative");
                        }
                        const auto id =
                            static_cast<std::uint64_t>(value);
                        if (!unique.insert(id).second) {
                            throw std::invalid_argument(
                                "mpmc::mesh::import_vtu_ascii: global cell IDs must be unique");
                        }
                        cell_global_ids.push_back(id);
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

    std::vector<GlobalEntityId> vertex_ids;
    vertex_ids.reserve(point_count);
    for (std::size_t point = 0; point < point_count; ++point) {
        vertex_ids.emplace_back(static_cast<std::uint64_t>(point) + 1U);
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
            if (name == global_vertex_id_name) {
                vertex_ids = parse_vertex_ids(array, point_count);
                continue;
            }
            parsed_point_fields.push_back(
                parse_field(
                    array,
                    EntityKind::vertex,
                    point_count));
        }
    }

    std::vector<LinearCell2D> cells;
    cells.reserve(cell_count);
    for (std::size_t i = 0; i < cell_count; ++i) {
        LinearCell2D cell{GlobalEntityId{cell_global_ids[i]},
            types[i] == 5 ? LinearCellType2D::triangle : LinearCellType2D::quadrilateral, {}};
        for (const auto point : cell_point_indices[i]) {
            cell.vertices.push_back(local_index(point,
                "mpmc::mesh::import_vtu_ascii: vertex local index overflow"));
        }
        cells.push_back(std::move(cell));
    }
    auto mesh = make_linear_mesh_2d(std::move(vertex_ids), std::move(coordinates), cells);
    auto& topology = mesh.topology;

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
        std::move(mesh.geometry),
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
    write_vertex_ids(output, topology);
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
