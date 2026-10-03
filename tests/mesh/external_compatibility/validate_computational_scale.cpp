#include <mpmc/mesh/cartesian_3d.hpp>
#include <mpmc/mesh/computational_mesh.hpp>
#include <mpmc/mesh/dof_numbering.hpp>
#include <mpmc/mesh/mesh_exchange_io.hpp>

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {
namespace mesh = mpmc::mesh;
using Clock = std::chrono::steady_clock;

void require_scale(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

double seconds_since(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

mesh::MeshTextExportResult export_scale(
    const mesh::MeshExchangeDocument& document, std::string_view format) {
    if (format == "vtu") return mesh::export_vtu_ascii(document);
    if (format == "gmsh") return mesh::export_gmsh_4_1_ascii(document);
    if (format == "grdecl") return mesh::export_grdecl_ascii(document);
    throw std::invalid_argument("scale format must be gmsh, vtu or grdecl");
}

mesh::MeshExchangeDocument import_scale(std::string_view content, std::string_view format) {
    if (format == "vtu") return mesh::make_mesh_exchange_document(mesh::import_vtu_ascii_3d(content));
    if (format == "gmsh") return mesh::make_mesh_exchange_document(mesh::import_gmsh_4_1_ascii_3d(content, 1.0));
    return mesh::make_mesh_exchange_document(mesh::import_grdecl(content, {1.0, 1.0}));
}
} // namespace

// Synthetic orthogonal grid: each cell is 1 x 2 x 3 metres. This validates
// file exchange -> checked geometry -> existing serial DoF numbering and a
// cell-centred two-point adjacency footprint. It is not a PDE solver or a
// general stencil policy. The analytic linear field has zero Laplacian.
void validate_computational_scale(const char* format_arg, const char* extent_arg,
                                  const char* output_stem) {
    const std::string_view format{format_arg};
    const std::string_view extent{extent_arg};
    std::size_t n = 0;
    const auto parsed = std::from_chars(extent.data(), extent.data() + extent.size(), n);
    if (parsed.ec != std::errc{} || parsed.ptr != extent.data() + extent.size() || n < 1U || n > 128U) {
        throw std::invalid_argument("scale extent must be an integer in [1,128]");
    }
    const auto total_start = Clock::now();
    std::vector<double> x(n + 1U), y(n + 1U), z(n + 1U);
    for (std::size_t i = 0; i <= n; ++i) {
        x[i] = static_cast<double>(i);
        y[i] = 2.0 * static_cast<double>(i);
        z[i] = 3.0 * static_cast<double>(i);
    }
    auto grid = mesh::make_cartesian_mesh_3d(x, y, z);
    const auto source = mesh::MeshExchangeDocument::create(
        mesh::MeshExchangeFormat::vtu_ascii, 3, std::move(grid.topology),
        std::move(grid.vertex_coordinates_m), std::move(grid.face_boundary), {}, {}, std::nullopt);
    const double generation_seconds = seconds_since(total_start);
    auto start = Clock::now();
    const auto exported = export_scale(source, format);
    require_scale(exported.exported(), "scale source export unsupported");
    const double export_seconds = seconds_since(start);
    const std::string stem{output_stem};
    const std::string suffix = format == "gmsh" ? ".msh" : "." + std::string{format};
    {
        std::ofstream file(stem + suffix, std::ios::binary);
        file << *exported.content;
        require_scale(static_cast<bool>(file), "scale mesh write failed");
    }
    start = Clock::now();
    const auto document = import_scale(*exported.content, format);
    const double import_seconds = seconds_since(start);
    start = Clock::now();
    const auto prepared = mesh::prepare_linear_mesh_3d(document);
    const auto geometry = mesh::make_cell_face_geometric_operator_3d(prepared);
    const double preparation_seconds = seconds_since(start);
    start = Clock::now();
    const auto& topology = prepared.topology;
    const std::size_t cells = topology.entity_count(mesh::EntityKind::cell);
    const std::size_t faces = topology.entity_count(mesh::EntityKind::face);
    require_scale(cells == n * n * n, "scale cell count mismatch");
    require_scale(faces == 3U * n * n * (n + 1U), "scale shared face count mismatch");
    const auto layout = mesh::DofLayout::create(topology, {{"synthetic.scalar", mesh::EntityKind::cell, 1U}});
    const auto numbering = mesh::DofNumberingSnapshot::create_serial(
        layout, mesh::make_serial_partition_snapshot(topology));
    require_scale(numbering.global_dof_count() == cells && numbering.ghost_dof_count() == 0U,
                  "scale DoF count mismatch");
    std::vector<double> linear_balance(cells, 0.0);
    std::vector<std::size_t> row_width(cells, 1U);
    std::size_t internal_faces = 0;
    for (std::size_t f = 0; f < faces; ++f) {
        const auto face = mesh::LocalIndex{static_cast<mesh::LocalIndex::value_type>(f)};
        const auto connection = geometry.face_connection_geometry(face);
        const auto owner = connection.owner.value();
        const auto normal = connection.owner_unit_normal;
        // F = grad(x + 2*y - 3*z); divergence F = 0, in physical metres.
        const double flux = connection.area_m2 * (normal.x + 2.0 * normal.y - 3.0 * normal.z);
        linear_balance[owner] += flux;
        if (connection.neighbour.has_value()) {
            ++internal_faces;
            const auto neighbour = connection.neighbour->value();
            linear_balance[neighbour] -= flux;
            ++row_width[owner];
            ++row_width[neighbour];
            const auto a = geometry.cell_reference_point_m(connection.owner);
            const auto b = geometry.cell_reference_point_m(*connection.neighbour);
            const double difference = (b.x_m - a.x_m) + 2.0 * (b.y_m - a.y_m) - 3.0 * (b.z_m - a.z_m);
            const double distance = connection.owner_normal_distance_m + *connection.neighbour_normal_distance_m;
            require_scale(std::abs(connection.area_m2 * difference / distance - flux) < 1e-11,
                          "scale affine two-point flux mismatch");
        }
    }
    std::size_t structural_entries = 0;
    double max_balance = 0.0;
    for (std::size_t c = 0; c < cells; ++c) {
        require_scale(std::abs(prepared.cell_volumes_m3[c] - 6.0) < 1e-11,
                      "scale analytic volume mismatch");
        max_balance = std::max(max_balance, std::abs(linear_balance[c]));
        require_scale(row_width[c] <= 7U, "scale duplicate/malformed cell coupling");
        require_scale(numbering.local_scalar(numbering.global_index(c)) == c, "scale DoF mapping mismatch");
        structural_entries += row_width[c];
    }
    require_scale(internal_faces == 3U * n * n * (n - 1U), "scale internal face mismatch");
    require_scale(structural_entries == cells + 2U * internal_faces, "scale structural count mismatch");
    require_scale(max_balance < 1e-11, "scale affine cell balance mismatch");
    const double assembly_preparation_seconds = seconds_since(start);
    std::ofstream report(stem + ".metrics.txt");
    report << std::setprecision(17)
           << "schema=mpmc.mesh.computational_scale.v1\nformat=" << format
           << "\nextent=" << n << "\ncells=" << cells << "\nfaces=" << faces
           << "\ninput_bytes=" << exported.content->size()
           << "\nstructural_entries=" << structural_entries
           << "\nmax_affine_balance=" << max_balance
           << "\ngeneration_seconds=" << generation_seconds
           << "\nexport_seconds=" << export_seconds << "\nimport_seconds=" << import_seconds
           << "\npreparation_seconds=" << preparation_seconds
           << "\nassembly_preparation_seconds=" << assembly_preparation_seconds
           << "\ntotal_seconds=" << seconds_since(total_start) << '\n';
    require_scale(static_cast<bool>(report), "scale metrics write failed");
    std::cout << "[PASS] computational.mesh.scale format=" << format << " cells=" << cells << '\n';
}
