#ifndef MPMC_MESH_MESH_EXCHANGE_IO_HPP
#define MPMC_MESH_MESH_EXCHANGE_IO_HPP

#include <mpmc/mesh/gmsh_4_1.hpp>
#include <mpmc/mesh/gmsh_4_1_3d.hpp>
#include <mpmc/mesh/grdecl.hpp>
#include <mpmc/mesh/mesh_exchange.hpp>
#include <mpmc/mesh/vtu.hpp>
#include <mpmc/mesh/vtu_3d.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
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
