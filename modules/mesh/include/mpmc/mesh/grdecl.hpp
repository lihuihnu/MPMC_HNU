#ifndef MPMC_MESH_GRDECL_HPP
#define MPMC_MESH_GRDECL_HPP

#include <mpmc/mesh/corner_point_geometry_3d.hpp>
#include <mpmc/mesh/csr_adjacency.hpp>
#include <mpmc/mesh/dense_field.hpp>
#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::mesh {

struct GrdeclImportOptions {
    double coordinate_scale_to_m;
    double permeability_scale_to_m2;
};

struct GrdeclImportResult {
    std::array<std::size_t, 3> dimensions;
    Topology topology;
    CornerPointGeometry3D geometry;
    std::vector<std::uint8_t> active;
    std::vector<DenseFieldSnapshot> cell_fields;
    std::vector<double> coord_m;
    std::vector<double> zcorn_m;
    GrdeclImportOptions source_options;

    [[nodiscard]] std::size_t cell_count() const noexcept {
        return active.size();
    }

    [[nodiscard]] bool is_active(LocalIndex cell) const {
        const std::size_t local =
            static_cast<std::size_t>(cell.value());
        if (local >= active.size()) {
            throw std::out_of_range(
                "mpmc::mesh::GrdeclImportResult: cell index out of range");
        }
        return active[local] != std::uint8_t{0U};
    }

    [[nodiscard]] std::size_t active_cell_count() const noexcept {
        return static_cast<std::size_t>(
            std::count(
                active.begin(),
                active.end(),
                std::uint8_t{1U}));
    }
};

namespace grdecl_detail {

using Records =
    std::map<std::string, std::vector<std::string>>;

[[nodiscard]] inline std::string uppercase(
    std::string_view text) {
    std::string result{text};
    for (char& character : result) {
        if (character >= 'a' &&
            character <= 'z') {
            character = static_cast<char>(
                character - 'a' + 'A');
        }
    }
    return result;
}

[[nodiscard]] inline std::vector<std::string>
tokenize(std::string_view content) {
    std::string cleaned;
    cleaned.reserve(content.size());

    std::size_t line_begin = 0U;
    while (line_begin < content.size()) {
        const auto line_end =
            content.find('\n', line_begin);
        const std::size_t end =
            line_end == std::string_view::npos
                ? content.size()
                : line_end;
        const auto line =
            content.substr(
                line_begin,
                end - line_begin);
        const auto comment =
            line.find("--");
        const auto visible =
            comment == std::string_view::npos
                ? line
                : line.substr(0U, comment);
        cleaned.append(
            visible.begin(),
            visible.end());
        cleaned.push_back('\n');

        if (line_end == std::string_view::npos) {
            break;
        }
        line_begin = line_end + 1U;
    }

    std::vector<std::string> tokens;
    std::string current;
    const auto flush = [&]() {
        if (!current.empty()) {
            tokens.push_back(
                std::move(current));
            current.clear();
        }
    };

    for (const char character : cleaned) {
        if (character == '/') {
            flush();
            tokens.emplace_back("/");
        } else if (character == ' ' ||
                   character == '\t' ||
                   character == '\n' ||
                   character == '\r') {
            flush();
        } else {
            current.push_back(character);
        }
    }
    flush();
    return tokens;
}

[[nodiscard]] inline Records parse_records(
    std::string_view content) {
    const auto tokens = tokenize(content);
    if (tokens.empty()) {
        throw std::invalid_argument(
            "mpmc::mesh::import_grdecl: input is empty");
    }

    const std::set<std::string> supported{
        "SPECGRID",
        "COORD",
        "ZCORN",
        "ACTNUM",
        "PORO",
        "PERMX",
        "PERMY",
        "PERMZ"};

    Records records;
    std::size_t cursor = 0U;
    while (cursor < tokens.size()) {
        if (tokens[cursor] == "/") {
            throw std::invalid_argument(
                "mpmc::mesh::import_grdecl: unexpected '/' outside keyword record");
        }
        const std::string keyword =
            uppercase(tokens[cursor++]);
        if (!supported.contains(keyword)) {
            throw std::invalid_argument(
                "mpmc::mesh::import_grdecl: unsupported keyword '" +
                keyword + "'");
        }
        if (records.contains(keyword)) {
            throw std::invalid_argument(
                "mpmc::mesh::import_grdecl: duplicate keyword '" +
                keyword + "'");
        }

        std::vector<std::string> values;
        while (cursor < tokens.size() &&
               tokens[cursor] != "/") {
            values.push_back(
                tokens[cursor++]);
        }
        if (cursor >= tokens.size()) {
            throw std::invalid_argument(
                "mpmc::mesh::import_grdecl: keyword '" +
                keyword + "' is missing '/' terminator");
        }
        ++cursor;
        records.emplace(
            keyword,
            std::move(values));
    }
    return records;
}

[[nodiscard]] inline const std::vector<std::string>&
require_record(const Records& records,
               std::string_view keyword) {
    const auto found =
        records.find(std::string(keyword));
    if (found == records.end()) {
        throw std::invalid_argument(
            "mpmc::mesh::import_grdecl: required keyword '" +
            std::string(keyword) +
            "' is missing");
    }
    return found->second;
}

template <typename Integer>
[[nodiscard]] inline Integer parse_integer(
    std::string_view token,
    const char* message) {
    Integer value{};
    const char* begin = token.data();
    const char* end =
        token.data() + token.size();
    const auto [next, error] =
        std::from_chars(
            begin, end, value);
    if (error != std::errc{} ||
        next != end) {
        throw std::invalid_argument(message);
    }
    return value;
}

[[nodiscard]] inline double parse_double(
    std::string_view token,
    const char* message) {
    std::string normalized{token};
    for (char& character : normalized) {
        if (character == 'd') character = 'e';
        if (character == 'D') character = 'E';
    }
    std::size_t consumed = 0U;
    double value = 0.0;
    try {
        value = std::stod(
            normalized,
            &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument(message);
    }
    if (consumed != normalized.size() ||
        !std::isfinite(value)) {
        throw std::invalid_argument(message);
    }
    return value;
}

template <typename ParseValue>
[[nodiscard]] inline auto expand_numeric_record(
    const std::vector<std::string>& tokens,
    ParseValue parse_value,
    const char* message) {
    using Value =
        decltype(parse_value(
            std::string_view{}, message));
    std::vector<Value> values;

    for (const auto& token : tokens) {
        const auto star =
            token.find('*');
        if (star == std::string::npos) {
            values.push_back(
                parse_value(token, message));
            continue;
        }
        if (token.find('*', star + 1U) !=
            std::string::npos ||
            star == 0U ||
            star + 1U == token.size()) {
            throw std::invalid_argument(
                "mpmc::mesh::import_grdecl: repeat syntax must be N*value; default-only repeats are unsupported for arrays");
        }

        const auto raw_count =
            parse_integer<std::uint64_t>(
                std::string_view{token}.substr(
                    0U, star),
                "mpmc::mesh::import_grdecl: invalid repeat count");
        if (raw_count == 0U) {
            throw std::invalid_argument(
                "mpmc::mesh::import_grdecl: repeat count must be positive");
        }
        if (raw_count >
            static_cast<std::uint64_t>(
                std::numeric_limits<std::size_t>::max() -
                values.size())) {
            throw std::length_error(
                "mpmc::mesh::import_grdecl: expanded array size overflow");
        }

        const Value value =
            parse_value(
                std::string_view{token}.substr(
                    star + 1U),
                message);
        values.insert(
            values.end(),
            static_cast<std::size_t>(
                raw_count),
            value);
    }
    return values;
}

[[nodiscard]] inline std::size_t checked_multiply(
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

[[nodiscard]] inline std::size_t checked_add(
    std::size_t left,
    std::size_t right,
    const char* message) {
    if (right >
        std::numeric_limits<std::size_t>::max() -
            left) {
        throw std::length_error(message);
    }
    return left + right;
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

[[nodiscard]] inline CsrAdjacency::Offset csr_offset(
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

struct Specgrid {
    std::size_t nx;
    std::size_t ny;
    std::size_t nz;
};

[[nodiscard]] inline Specgrid parse_specgrid(
    const std::vector<std::string>& record) {
    if (record.size() < 3U ||
        record.size() > 5U) {
        throw std::invalid_argument(
            "mpmc::mesh::import_grdecl: SPECGRID requires 3..5 items");
    }

    const auto parse_dimension =
        [](const std::string& token) {
            const auto value =
                parse_integer<std::uint64_t>(
                    token,
                    "mpmc::mesh::import_grdecl: invalid SPECGRID dimension");
            if (value == 0U ||
                value >
                    static_cast<std::uint64_t>(
                        std::numeric_limits<std::size_t>::max())) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_grdecl: SPECGRID dimensions must be positive and representable");
            }
            return static_cast<std::size_t>(
                value);
        };

    const std::size_t nx =
        parse_dimension(record[0]);
    const std::size_t ny =
        parse_dimension(record[1]);
    const std::size_t nz =
        parse_dimension(record[2]);

    std::uint64_t numres = 1U;
    if (record.size() >= 4U &&
        record[3] != "*" &&
        record[3] != "1*") {
        numres =
            parse_integer<std::uint64_t>(
                record[3],
                "mpmc::mesh::import_grdecl: invalid SPECGRID NUMRES");
    }
    if (numres != 1U) {
        throw std::invalid_argument(
            "mpmc::mesh::import_grdecl: minimal baseline requires SPECGRID NUMRES=1");
    }

    std::string type = "F";
    if (record.size() >= 5U &&
        record[4] != "*" &&
        record[4] != "1*") {
        type = uppercase(record[4]);
    }
    if (type != "F") {
        throw std::invalid_argument(
            "mpmc::mesh::import_grdecl: radial SPECGRID type is unsupported; Cartesian F is required");
    }

    return Specgrid{nx, ny, nz};
}

struct Pillar {
    Coordinate3D top;
    Coordinate3D bottom;
};

[[nodiscard]] inline Coordinate3D point_on_pillar(
    const Pillar& pillar,
    double z_m) {
    const double dz =
        pillar.bottom.z_m -
        pillar.top.z_m;
    const double scale =
        std::max(
            1.0,
            std::max(
                std::abs(pillar.top.z_m),
                std::abs(pillar.bottom.z_m)));
    const double tolerance =
        128.0 *
        std::numeric_limits<double>::epsilon() *
        scale;
    if (std::abs(dz) <= tolerance) {
        throw std::invalid_argument(
            "mpmc::mesh::import_grdecl: pillar has zero z extent and cannot be interpolated");
    }

    const double t =
        (z_m - pillar.top.z_m) /
        dz;
    const Coordinate3D point{
        pillar.top.x_m +
            t * (pillar.bottom.x_m -
                 pillar.top.x_m),
        pillar.top.y_m +
            t * (pillar.bottom.y_m -
                 pillar.top.y_m),
        z_m};
    if (!std::isfinite(point.x_m) ||
        !std::isfinite(point.y_m) ||
        !std::isfinite(point.z_m)) {
        throw std::invalid_argument(
            "mpmc::mesh::import_grdecl: non-finite interpolated corner point");
    }
    return point;
}

[[nodiscard]] inline double determinant(
    Coordinate3D a,
    Coordinate3D b,
    Coordinate3D c) {
    return
        a.x_m *
            (b.y_m * c.z_m -
             b.z_m * c.y_m) -
        a.y_m *
            (b.x_m * c.z_m -
             b.z_m * c.x_m) +
        a.z_m *
            (b.x_m * c.y_m -
             b.y_m * c.x_m);
}

[[nodiscard]] inline Coordinate3D subtract(
    Coordinate3D left,
    Coordinate3D right) {
    return Coordinate3D{
        left.x_m - right.x_m,
        left.y_m - right.y_m,
        left.z_m - right.z_m};
}

[[nodiscard]] inline double signed_tetra_volume(
    Coordinate3D a,
    Coordinate3D b,
    Coordinate3D c,
    Coordinate3D d) {
    return determinant(
               subtract(b, a),
               subtract(c, a),
               subtract(d, a)) /
           6.0;
}

struct CellVolumeCheck {
    double volume_m3;
    bool degenerate;
};

[[nodiscard]] inline CellVolumeCheck
check_cell_volume(
    const std::array<Coordinate3D, 8>& corners,
    bool active) {
    double min_x = corners.front().x_m;
    double max_x = corners.front().x_m;
    double min_y = corners.front().y_m;
    double max_y = corners.front().y_m;
    double min_z = corners.front().z_m;
    double max_z = corners.front().z_m;
    for (const auto point : corners) {
        min_x = std::min(min_x, point.x_m);
        max_x = std::max(max_x, point.x_m);
        min_y = std::min(min_y, point.y_m);
        max_y = std::max(max_y, point.y_m);
        min_z = std::min(min_z, point.z_m);
        max_z = std::max(max_z, point.z_m);
    }
    const double extent =
        std::max(
            {max_x - min_x,
             max_y - min_y,
             max_z - min_z,
             1.0});
    const double tolerance =
        4096.0 *
        std::numeric_limits<double>::epsilon() *
        extent * extent * extent;

    static constexpr std::array<
        std::array<std::size_t, 4>, 5>
        tetrahedra{{
            {{0U, 1U, 2U, 4U}},
            {{1U, 3U, 2U, 7U}},
            {{1U, 2U, 4U, 7U}},
            {{1U, 4U, 5U, 7U}},
            {{2U, 4U, 7U, 6U}}
        }};

    double volume = 0.0;
    bool degenerate = false;
    for (const auto& tetrahedron :
         tetrahedra) {
        const double tetra_volume =
            signed_tetra_volume(
                corners[tetrahedron[0]],
                corners[tetrahedron[1]],
                corners[tetrahedron[2]],
                corners[tetrahedron[3]]);
        if (!std::isfinite(tetra_volume)) {
            throw std::invalid_argument(
                "mpmc::mesh::import_grdecl: non-finite cell tetrahedron volume");
        }
        if (tetra_volume < -tolerance) {
            throw std::invalid_argument(
                "mpmc::mesh::import_grdecl: flipped/inverted corner-point cell");
        }
        if (tetra_volume <= tolerance) {
            degenerate = true;
        }
        volume +=
            std::max(0.0, tetra_volume);
    }

    if (volume <= tolerance) {
        degenerate = true;
        volume = 0.0;
    }
    if (active && degenerate) {
        throw std::invalid_argument(
            "mpmc::mesh::import_grdecl: active corner-point cell is degenerate");
    }
    return CellVolumeCheck{
        volume, degenerate};
}

[[nodiscard]] inline DenseFieldMetadata
field_metadata(std::string id,
               std::string unit) {
    return DenseFieldMetadata{
        id,
        unit,
        FieldSourceMetadata{
            FieldSourceKind::file_import,
            "GRDECL input",
            "SPECGRID/COORD/ZCORN minimal baseline",
            id}};
}

} // namespace grdecl_detail

[[nodiscard]] inline GrdeclImportResult
import_grdecl(
    std::string_view content,
    GrdeclImportOptions options) {
    using namespace grdecl_detail;

    if (!std::isfinite(
            options.coordinate_scale_to_m) ||
        options.coordinate_scale_to_m <= 0.0) {
        throw std::invalid_argument(
            "mpmc::mesh::import_grdecl: coordinate_scale_to_m must be finite and positive");
    }
    if (!std::isfinite(
            options.permeability_scale_to_m2) ||
        options.permeability_scale_to_m2 <= 0.0) {
        throw std::invalid_argument(
            "mpmc::mesh::import_grdecl: permeability_scale_to_m2 must be finite and positive");
    }

    const Records records =
        parse_records(content);
    const Specgrid spec =
        parse_specgrid(
            require_record(
                records, "SPECGRID"));

    const std::size_t nxp1 =
        checked_add(
            spec.nx, 1U,
            "mpmc::mesh::import_grdecl: NX+1 overflow");
    const std::size_t nyp1 =
        checked_add(
            spec.ny, 1U,
            "mpmc::mesh::import_grdecl: NY+1 overflow");
    const std::size_t pillar_count =
        checked_multiply(
            nxp1, nyp1,
            "mpmc::mesh::import_grdecl: pillar count overflow");
    const std::size_t cell_count =
        checked_multiply(
            checked_multiply(
                spec.nx, spec.ny,
                "mpmc::mesh::import_grdecl: XY cell count overflow"),
            spec.nz,
            "mpmc::mesh::import_grdecl: XYZ cell count overflow");

    const auto coord =
        expand_numeric_record(
            require_record(records, "COORD"),
            [](std::string_view token,
               const char* message) {
                return parse_double(
                    token, message);
            },
            "mpmc::mesh::import_grdecl: invalid COORD value");
    const auto zcorn =
        expand_numeric_record(
            require_record(records, "ZCORN"),
            [](std::string_view token,
               const char* message) {
                return parse_double(
                    token, message);
            },
            "mpmc::mesh::import_grdecl: invalid ZCORN value");
    const auto actnum =
        expand_numeric_record(
            require_record(records, "ACTNUM"),
            [](std::string_view token,
               const char* message) {
                return parse_integer<std::int64_t>(
                    token, message);
            },
            "mpmc::mesh::import_grdecl: invalid ACTNUM value");
    const auto poro =
        expand_numeric_record(
            require_record(records, "PORO"),
            [](std::string_view token,
               const char* message) {
                return parse_double(
                    token, message);
            },
            "mpmc::mesh::import_grdecl: invalid PORO value");
    const auto permx =
        expand_numeric_record(
            require_record(records, "PERMX"),
            [](std::string_view token,
               const char* message) {
                return parse_double(
                    token, message);
            },
            "mpmc::mesh::import_grdecl: invalid PERMX value");
    const auto permy =
        expand_numeric_record(
            require_record(records, "PERMY"),
            [](std::string_view token,
               const char* message) {
                return parse_double(
                    token, message);
            },
            "mpmc::mesh::import_grdecl: invalid PERMY value");
    const auto permz =
        expand_numeric_record(
            require_record(records, "PERMZ"),
            [](std::string_view token,
               const char* message) {
                return parse_double(
                    token, message);
            },
            "mpmc::mesh::import_grdecl: invalid PERMZ value");

    const std::size_t expected_coord =
        checked_multiply(
            pillar_count, 6U,
            "mpmc::mesh::import_grdecl: COORD size overflow");
    const std::size_t expected_zcorn =
        checked_multiply(
            cell_count, 8U,
            "mpmc::mesh::import_grdecl: ZCORN size overflow");
    if (coord.size() != expected_coord) {
        throw std::invalid_argument(
            "mpmc::mesh::import_grdecl: COORD must contain 6*(NX+1)*(NY+1) values");
    }
    if (zcorn.size() != expected_zcorn) {
        throw std::invalid_argument(
            "mpmc::mesh::import_grdecl: ZCORN must contain 8*NX*NY*NZ values");
    }
    for (const auto size :
         {actnum.size(),
          poro.size(),
          permx.size(),
          permy.size(),
          permz.size()}) {
        if (size != cell_count) {
            throw std::invalid_argument(
                "mpmc::mesh::import_grdecl: ACTNUM/PORO/PERM arrays must contain NX*NY*NZ values");
        }
    }

    std::vector<double> coord_m;
    coord_m.reserve(coord.size());
    for (const double value : coord) {
        const double scaled =
            value * options.coordinate_scale_to_m;
        if (!std::isfinite(scaled)) {
            throw std::invalid_argument(
                "mpmc::mesh::import_grdecl: scaled COORD value is non-finite");
        }
        coord_m.push_back(scaled);
    }
    std::vector<double> zcorn_m;
    zcorn_m.reserve(zcorn.size());
    for (const double value : zcorn) {
        const double scaled =
            value * options.coordinate_scale_to_m;
        if (!std::isfinite(scaled)) {
            throw std::invalid_argument(
                "mpmc::mesh::import_grdecl: scaled ZCORN value is non-finite");
        }
        zcorn_m.push_back(scaled);
    }

    std::vector<Pillar> pillars;
    pillars.reserve(pillar_count);
    for (std::size_t pillar = 0U;
         pillar < pillar_count;
         ++pillar) {
        const std::size_t base =
            pillar * 6U;
        pillars.push_back(
            Pillar{
                Coordinate3D{
                    coord_m[base],
                    coord_m[base + 1U],
                    coord_m[base + 2U]},
                Coordinate3D{
                    coord_m[base + 3U],
                    coord_m[base + 4U],
                    coord_m[base + 5U]}});
    }

    std::vector<std::uint8_t> active;
    active.reserve(cell_count);
    for (const auto value : actnum) {
        if (value != 0 &&
            value != 1) {
            throw std::invalid_argument(
                "mpmc::mesh::import_grdecl: ACTNUM values must be 0 or 1");
        }
        active.push_back(
            value == 0
                ? std::uint8_t{0U}
                : std::uint8_t{1U});
    }

    for (const double value : poro) {
        if (value < 0.0 || value > 1.0) {
            throw std::invalid_argument(
                "mpmc::mesh::import_grdecl: PORO must lie in [0,1]");
        }
    }
    for (const auto* values :
         {&permx, &permy, &permz}) {
        for (const double value : *values) {
            if (value < 0.0) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_grdecl: permeability values must be nonnegative");
            }
        }
    }

    if (cell_count >
        std::numeric_limits<std::size_t>::max() /
            std::size_t{8U}) {
        throw std::length_error(
            "mpmc::mesh::import_grdecl: corner vertex count overflow");
    }
    const std::size_t vertex_count =
        cell_count * std::size_t{8U};
    if (vertex_count >
        static_cast<std::size_t>(
            std::numeric_limits<
                LocalIndex::value_type>::max()) +
            std::size_t{1U}) {
        throw std::length_error(
            "mpmc::mesh::import_grdecl: corner vertex count exceeds LocalIndex capacity");
    }

    Topology::EntityIds ids;
    ids.cells.reserve(cell_count);
    ids.vertices.reserve(vertex_count);
    for (std::size_t cell = 0U;
         cell < cell_count;
         ++cell) {
        ids.cells.emplace_back(
            static_cast<std::uint64_t>(
                cell) +
            1U);
        for (std::size_t corner = 0U;
             corner < 8U;
             ++corner) {
            ids.vertices.emplace_back(
                static_cast<std::uint64_t>(
                    cell * 8U + corner) +
                1U);
        }
    }

    std::vector<CsrAdjacency::Offset>
        cell_vertex_offsets;
    std::vector<LocalIndex>
        cell_vertices;
    cell_vertex_offsets.reserve(
        cell_count + 1U);
    cell_vertices.reserve(vertex_count);
    cell_vertex_offsets.push_back(
        CsrAdjacency::Offset{0U});

    std::vector<Coordinate3D>
        vertex_coordinates;
    std::vector<double> cell_volumes;
    vertex_coordinates.reserve(vertex_count);
    cell_volumes.reserve(cell_count);

    const std::size_t doubled_nx =
        checked_multiply(
            spec.nx, 2U,
            "mpmc::mesh::import_grdecl: 2*NX overflow");
    const std::size_t doubled_ny =
        checked_multiply(
            spec.ny, 2U,
            "mpmc::mesh::import_grdecl: 2*NY overflow");

    for (std::size_t k = 0U;
         k < spec.nz;
         ++k) {
        for (std::size_t j = 0U;
             j < spec.ny;
             ++j) {
            for (std::size_t i = 0U;
                 i < spec.nx;
                 ++i) {
                const std::size_t cell =
                    i +
                    spec.nx *
                        (j + spec.ny * k);

                const std::array<std::size_t, 4>
                    pillar_indices{
                        i + nxp1 * j,
                        i + 1U + nxp1 * j,
                        i + nxp1 * (j + 1U),
                        i + 1U +
                            nxp1 * (j + 1U)};

                std::array<Coordinate3D, 8>
                    corners{};
                for (std::size_t corner = 0U;
                     corner < 8U;
                     ++corner) {
                    const std::size_t local_i =
                        corner % 2U;
                    const std::size_t local_j =
                        (corner / 2U) % 2U;
                    const std::size_t local_k =
                        corner / 4U;
                    const std::size_t z_index =
                        (2U * i + local_i) +
                        doubled_nx *
                            ((2U * j + local_j) +
                             doubled_ny *
                                 (2U * k + local_k));
                    if (z_index >= zcorn.size()) {
                        throw std::logic_error(
                            "mpmc::mesh::import_grdecl: internal ZCORN index overflow");
                    }
                    const double z_m =
                        zcorn_m[z_index];
                    const std::size_t pillar_slot =
                        local_i +
                        2U * local_j;
                    corners[corner] =
                        point_on_pillar(
                            pillars[
                                pillar_indices[
                                    pillar_slot]],
                            z_m);
                }

                const auto volume =
                    check_cell_volume(
                        corners,
                        active[cell] !=
                            std::uint8_t{0U});
                cell_volumes.push_back(
                    volume.volume_m3);

                for (std::size_t corner = 0U;
                     corner < 8U;
                     ++corner) {
                    const std::size_t vertex =
                        cell * 8U + corner;
                    vertex_coordinates.push_back(
                        corners[corner]);
                    cell_vertices.push_back(
                        local_index(
                            vertex,
                            "mpmc::mesh::import_grdecl: vertex local index overflow"));
                }
                cell_vertex_offsets.push_back(
                    csr_offset(
                        cell_vertices.size(),
                        "mpmc::mesh::import_grdecl: cell->vertex CSR overflow"));
            }
        }
    }

    std::vector<CsrAdjacency> relations;
    relations.emplace_back(
        EntityKind::cell,
        EntityKind::vertex,
        vertex_count,
        std::move(cell_vertex_offsets),
        std::move(cell_vertices));

    Topology topology{
        std::move(ids),
        std::move(relations)};
    CornerPointGeometry3D geometry{
        std::move(vertex_coordinates),
        std::move(cell_volumes)};

    const auto make_field =
        [&](std::string id,
            std::string unit,
            const std::vector<double>& values,
            double scale) {
        std::vector<double> converted;
        converted.reserve(values.size());
        for (const double value : values) {
            const double scaled =
                value * scale;
            if (!std::isfinite(scaled)) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_grdecl: scaled property is non-finite");
            }
            converted.push_back(scaled);
        }
        return DenseFieldSnapshot::create(
            topology,
            EntityKind::cell,
            1U,
            std::move(converted),
            field_metadata(
                std::move(id),
                std::move(unit)));
    };

    std::vector<DenseFieldSnapshot>
        fields;
    fields.reserve(4U);
    fields.push_back(
        make_field(
            "PORO", "1", poro, 1.0));
    fields.push_back(
        make_field(
            "PERMX", "m2", permx,
            options.permeability_scale_to_m2));
    fields.push_back(
        make_field(
            "PERMY", "m2", permy,
            options.permeability_scale_to_m2));
    fields.push_back(
        make_field(
            "PERMZ", "m2", permz,
            options.permeability_scale_to_m2));

    return GrdeclImportResult{
        {spec.nx, spec.ny, spec.nz},
        std::move(topology),
        std::move(geometry),
        std::move(active),
        std::move(fields),
        std::move(coord_m),
        std::move(zcorn_m),
        options};
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_GRDECL_HPP
