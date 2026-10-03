#include <mpmc/mesh/cartesian_3d.hpp>
#include <mpmc/mesh/dense_field.hpp>
#include <mpmc/mesh/face_boundary.hpp>
#include <mpmc/mesh/face_geometry_3d.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#if defined(__linux__)
#include <sys/resource.h>
#elif defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <psapi.h>
#endif

namespace mesh = mpmc::mesh;

namespace {

struct Workload {
    std::size_t nx{64U};
    std::size_t ny{64U};
    std::size_t nz{64U};
};
constexpr std::size_t field_components = 4U;
constexpr std::size_t sample_count = 5U;
constexpr std::size_t topology_repetitions = 8U;
constexpr std::size_t field_repetitions = 64U;
constexpr double gib = 1024.0 * 1024.0 * 1024.0;

[[nodiscard]] std::size_t checked_add(
    std::size_t left,
    std::size_t right) {
    if (left >
        std::numeric_limits<std::size_t>::max() -
            right) {
        throw std::length_error(
            "mesh baseline byte count overflow");
    }
    return left + right;
}

[[nodiscard]] std::size_t checked_multiply(
    std::size_t left,
    std::size_t right) {
    if (left != 0U &&
        right >
            std::numeric_limits<std::size_t>::max() /
                left) {
        throw std::length_error(
            "mesh baseline byte count overflow");
    }
    return left * right;
}

template <typename T>
[[nodiscard]] std::size_t span_bytes(
    std::span<const T> values) {
    return checked_multiply(
        values.size(),
        sizeof(T));
}

struct ExpectedCounts {
    std::size_t cells;
    std::size_t vertices;
    std::size_t faces;
    std::size_t boundary_faces;
};

[[nodiscard]] ExpectedCounts expected_counts(Workload size) {
    const auto product3 = [](std::size_t a, std::size_t b, std::size_t c) {
        return checked_multiply(checked_multiply(a, b), c);
    };
    const std::size_t cells = product3(size.nx, size.ny, size.nz);
    const std::size_t vertices = product3(
        checked_add(size.nx, 1U), checked_add(size.ny, 1U), checked_add(size.nz, 1U));
    const std::size_t faces = checked_add(
        checked_add(product3(checked_add(size.nx, 1U), size.ny, size.nz),
                    product3(size.nx, checked_add(size.ny, 1U), size.nz)),
        product3(size.nx, size.ny, checked_add(size.nz, 1U)));
    const std::size_t boundary = checked_multiply(2U, checked_add(
        checked_add(checked_multiply(size.nx, size.ny), checked_multiply(size.nx, size.nz)),
        checked_multiply(size.ny, size.nz)));
    const std::uint64_t capacity =
        static_cast<std::uint64_t>(std::numeric_limits<mesh::LocalIndex::value_type>::max()) + 1ULL;
    // Check all dimensions/counts before axis allocation. The Cartesian builder
    // also reserves up to two face-to-cell entries per face.
    if (cells == 0U || cells > capacity || vertices > capacity || faces > capacity ||
        checked_multiply(cells, 8U) > std::numeric_limits<mesh::CsrAdjacency::Offset>::max() ||
        checked_multiply(faces, 4U) > std::numeric_limits<mesh::CsrAdjacency::Offset>::max()) {
        throw std::length_error("mesh benchmark workload exceeds local/CSR capacity");
    }
    return {cells, vertices, faces, boundary};
}

[[nodiscard]] std::size_t parse_dimension(std::string_view text) {
    std::size_t value{};
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (text.empty() || parsed.ec != std::errc{} ||
        parsed.ptr != text.data() + text.size() || value == 0U) {
        throw std::invalid_argument("mesh benchmark dimensions must be positive decimal integers");
    }
    return value;
}

[[nodiscard]] Workload parse_workload(std::span<const std::string_view> args) {
    Workload result;
    if (!args.empty()) {
        if (args.size() != 4U || args[0] != "--grid") {
            throw std::invalid_argument("usage: mpmc_mesh_baseline_benchmark [--grid NX NY NZ]");
        }
        result = {parse_dimension(args[1]), parse_dimension(args[2]), parse_dimension(args[3])};
    }
    static_cast<void>(expected_counts(result));
    return result;
}

// This benchmark is its own existing CI owner. Exercise the new configuration
// contract on every invocation, including the unchanged no-argument cloud run.
void verify_configuration_contract() {
    const auto default_size = parse_workload({});
    const auto default_counts = expected_counts(default_size);
    const std::array<std::string_view, 4U> asymmetric{"--grid", "2", "3", "4"};
    const auto counts = expected_counts(parse_workload(asymmetric));
    if (default_size.nx != 64U || default_size.ny != 64U || default_size.nz != 64U ||
        default_counts.cells != 262144U || counts.cells != 24U ||
        counts.vertices != 60U || counts.faces != 98U || counts.boundary_faces != 52U) {
        throw std::logic_error("mesh benchmark configuration regression");
    }
    for (const auto invalid : {"", "0", "-1", "+1", "1.5", "2x", " 2", "2 ",
                               "18446744073709551616"}) {
        bool rejected = false;
        try {
            static_cast<void>(parse_dimension(invalid));
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        if (!rejected) {
            throw std::logic_error("invalid benchmark dimension accepted");
        }
    }
    const std::array<std::string_view, 3U> missing{"--grid", "2", "3"};
    bool rejected = false;
    try {
        static_cast<void>(parse_workload(missing));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    if (!rejected) {
        throw std::logic_error("incomplete benchmark workload accepted");
    }
    for (const auto oversized : {Workload{std::numeric_limits<std::size_t>::max(), 1U, 1U},
                                 Workload{65536U, 65536U, 1U}}) {
        rejected = false;
        try {
            static_cast<void>(expected_counts(oversized));
        } catch (const std::length_error&) {
            rejected = true;
        }
        if (!rejected) {
            throw std::logic_error("oversized benchmark workload accepted");
        }
    }
}

[[nodiscard]] std::vector<double> axis(
    std::size_t cell_count) {
    std::vector<double> values;
    values.reserve(cell_count + 1U);
    for (std::size_t index = 0U;
         index <= cell_count;
         ++index) {
        values.push_back(
            static_cast<double>(index));
    }
    return values;
}

[[nodiscard]] mesh::DenseFieldMetadata
field_metadata(Workload size) {
    return mesh::DenseFieldMetadata{
        "benchmark.cell.field4",
        "1",
        mesh::FieldSourceMetadata{
            mesh::FieldSourceKind::synthetic_test,
            "benchmark://mesh/cartesian/" + std::to_string(size.nx) + "x" +
                std::to_string(size.ny) + "x" + std::to_string(size.nz),
            "baseline-v1",
            "mesh_baseline_benchmark.cpp"}};
}

[[nodiscard]] std::size_t topology_payload_bytes(
    const mesh::Topology& topology) {
    std::size_t bytes = 0U;
    for (const auto kind :
         {mesh::EntityKind::vertex,
          mesh::EntityKind::edge,
          mesh::EntityKind::face,
          mesh::EntityKind::cell}) {
        bytes = checked_add(
            bytes,
            span_bytes(
                topology.global_ids(kind)));
    }

    const std::array<
        std::pair<
            mesh::EntityKind,
            mesh::EntityKind>,
        4>
        relations{{
            {mesh::EntityKind::cell,
             mesh::EntityKind::vertex},
            {mesh::EntityKind::cell,
             mesh::EntityKind::face},
            {mesh::EntityKind::face,
             mesh::EntityKind::vertex},
            {mesh::EntityKind::face,
             mesh::EntityKind::cell},
        }};
    for (const auto& relation :
         relations) {
        const auto& csr =
            topology.relation(
                relation.first,
                relation.second);
        bytes = checked_add(
            bytes,
            span_bytes(csr.offsets()));
        bytes = checked_add(
            bytes,
            span_bytes(csr.indices()));
    }
    return bytes;
}

[[nodiscard]] std::size_t geometry_payload_bytes(
    const mesh::LinearMesh3D& grid) {
    std::size_t bytes = 0U;
    bytes = checked_add(
        bytes,
        span_bytes(
            std::span<const mesh::Coordinate3D>{
                grid.vertex_coordinates_m}));
    bytes = checked_add(
        bytes,
        span_bytes(
            std::span<const double>{
                grid.cell_volumes_m3}));
    bytes = checked_add(
        bytes,
        span_bytes(
            grid.face_geometry.face_centroids_m()));
    bytes = checked_add(
        bytes,
        span_bytes(
            grid.face_geometry.face_areas_m2()));
    bytes = checked_add(
        bytes,
        span_bytes(
            grid.face_geometry.face_owners()));
    bytes = checked_add(
        bytes,
        span_bytes(
            grid.face_geometry
                .face_owner_unit_normals()));
    return bytes;
}

[[nodiscard]] std::size_t boundary_payload_bytes(
    const mesh::FaceBoundarySnapshot& boundary) {
    return checked_add(
        span_bytes(boundary.classifications()),
        span_bytes(boundary.physical_tags()));
}

#if defined(_MSC_VER)
#define MPMC_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define MPMC_NOINLINE __attribute__((noinline))
#else
#define MPMC_NOINLINE
#endif

MPMC_NOINLINE
[[nodiscard]] std::uint64_t
traverse_topology_once(
    const mesh::Topology& topology) {
    std::uint64_t checksum = 1469598103934665603ULL;
    for (const auto kind :
         {mesh::EntityKind::vertex,
          mesh::EntityKind::edge,
          mesh::EntityKind::face,
          mesh::EntityKind::cell}) {
        for (const auto id :
             topology.global_ids(kind)) {
            checksum ^=
                id.value() +
                0x9e3779b97f4a7c15ULL;
            checksum *=
                1099511628211ULL;
        }
    }

    const std::array<
        std::pair<
            mesh::EntityKind,
            mesh::EntityKind>,
        4>
        relations{{
            {mesh::EntityKind::cell,
             mesh::EntityKind::vertex},
            {mesh::EntityKind::cell,
             mesh::EntityKind::face},
            {mesh::EntityKind::face,
             mesh::EntityKind::vertex},
            {mesh::EntityKind::face,
             mesh::EntityKind::cell},
        }};
    for (const auto& relation :
         relations) {
        const auto& csr =
            topology.relation(
                relation.first,
                relation.second);
        for (const auto offset :
             csr.offsets()) {
            checksum ^=
                static_cast<std::uint64_t>(
                    offset) +
                0x517cc1b727220a95ULL;
            checksum *=
                1099511628211ULL;
        }
        for (const auto index :
             csr.indices()) {
            checksum ^=
                static_cast<std::uint64_t>(
                    index.value()) +
                0x94d049bb133111ebULL;
            checksum *=
                1099511628211ULL;
        }
    }
    return checksum;
}

MPMC_NOINLINE
[[nodiscard]] std::uint64_t
traverse_field_once(
    const mesh::DenseFieldSnapshot& field) {
    double sum = 0.0;
    double weighted = 0.0;
    std::size_t component = 0U;
    for (const double value :
         field.values()) {
        sum += value;
        weighted +=
            value *
            static_cast<double>(
                component + 1U);
        component =
            (component + 1U) %
            field.component_count();
    }
    return std::bit_cast<std::uint64_t>(sum) ^
           std::bit_cast<std::uint64_t>(weighted);
}

template <typename Function>
[[nodiscard]] double measure_seconds(
    Function&& function,
    std::size_t repetitions,
    std::uint64_t& checksum) {
    const auto begin =
        std::chrono::steady_clock::now();
    for (std::size_t repeat = 0U;
         repeat < repetitions;
         ++repeat) {
        checksum ^=
            function() +
            static_cast<std::uint64_t>(
                repeat + 1U) *
                0x9e3779b97f4a7c15ULL;
    }
    const auto end =
        std::chrono::steady_clock::now();
    const double seconds =
        std::chrono::duration<double>(
            end - begin)
            .count();
    if (!std::isfinite(seconds) ||
        seconds <= 0.0) {
        throw std::runtime_error(
            "non-positive benchmark duration");
    }
    return seconds;
}

[[nodiscard]] double median(
    std::array<double, sample_count> values) {
    std::sort(
        values.begin(),
        values.end());
    return values[sample_count / 2U];
}

[[nodiscard]] std::uint64_t
peak_rss_kib() {
#if defined(__linux__)
    rusage usage{};
    if (getrusage(
            RUSAGE_SELF,
            &usage) != 0 ||
        usage.ru_maxrss < 0) {
        throw std::runtime_error(
            "getrusage failed");
    }
    return static_cast<std::uint64_t>(
        usage.ru_maxrss);
#elif defined(_WIN32)
    PROCESS_MEMORY_COUNTERS counters{};
    counters.cb = static_cast<DWORD>(sizeof(counters));
    if (GetProcessMemoryInfo(GetCurrentProcess(), &counters, counters.cb) == 0) {
        throw std::runtime_error("GetProcessMemoryInfo failed");
    }
    return static_cast<std::uint64_t>(counters.PeakWorkingSetSize) / 1024ULL;
#else
    return 0U;
#endif
}

[[nodiscard]] const char* peak_rss_source() {
#if defined(__linux__)
    return "linux_getrusage_ru_maxrss";
#elif defined(_WIN32)
    return "windows_peak_working_set";
#else
    return "unavailable";
#endif
}

void verify_grid(const mesh::LinearMesh3D& grid, const ExpectedCounts& expected) {
    const auto& topology = grid.topology;
    const auto boundaries = grid.face_boundary.classifications();
    if (topology.entity_count(mesh::EntityKind::cell) != expected.cells ||
        topology.entity_count(mesh::EntityKind::vertex) != expected.vertices ||
        topology.entity_count(mesh::EntityKind::face) != expected.faces ||
        topology.relation(mesh::EntityKind::cell, mesh::EntityKind::vertex).entry_count() != 8U * expected.cells ||
        topology.relation(mesh::EntityKind::cell, mesh::EntityKind::face).entry_count() != 6U * expected.cells ||
        topology.relation(mesh::EntityKind::face, mesh::EntityKind::vertex).entry_count() != 4U * expected.faces ||
        topology.relation(mesh::EntityKind::face, mesh::EntityKind::cell).entry_count() != 6U * expected.cells ||
        static_cast<std::size_t>(std::count(boundaries.begin(), boundaries.end(),
                                          mesh::FaceClassification::boundary)) != expected.boundary_faces ||
        grid.cell_volumes_m3.size() != expected.cells ||
        !std::all_of(grid.cell_volumes_m3.begin(), grid.cell_volumes_m3.end(),
                     [](double volume) { return volume == 1.0; })) {
        throw std::runtime_error("Cartesian benchmark analytic invariant failure");
    }
}

void print_samples(
    std::string_view key,
    const std::array<double, sample_count>& values) {
    std::cout << key << '=';
    for (std::size_t index = 0U;
         index < values.size();
         ++index) {
        if (index != 0U) {
            std::cout << ',';
        }
        std::cout << values[index];
    }
    std::cout << '\n';
}

} // namespace

int main(int argc, char** argv) {
    try {
        verify_configuration_contract();
        if (argc == 2 && std::string_view{argv[1]} == "--help") {
            std::cout << "usage: mpmc_mesh_baseline_benchmark [--grid NX NY NZ]\n";
            return 0;
        }
        std::vector<std::string_view> arguments;
        for (int index = 1; index < argc; ++index) {
            arguments.emplace_back(argv[index]);
        }
        const auto workload = parse_workload(arguments);
        const auto [nx, ny, nz] = workload;
        const auto expected = expected_counts(workload);
        const auto construction_begin = std::chrono::steady_clock::now();
        const auto x = axis(nx);
        const auto y = axis(ny);
        const auto z = axis(nz);
        const auto grid =
            mesh::make_cartesian_mesh_3d(
                x,
                y,
                z);

        const double mesh_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - construction_begin).count();
        const auto rss_after_mesh = peak_rss_kib();
        verify_grid(grid, expected);

        const std::size_t cell_count =
            grid.topology.entity_count(
                mesh::EntityKind::cell);
        const auto field_begin = std::chrono::steady_clock::now();
        std::vector<double> values;
        values.resize(
            checked_multiply(
                cell_count,
                field_components));
        for (std::size_t cell = 0U;
             cell < cell_count;
             ++cell) {
            for (std::size_t component = 0U;
                 component < field_components;
                 ++component) {
                values[
                    cell * field_components +
                    component] =
                    static_cast<double>(
                        (cell % 1024U) +
                        component + 1U) /
                    1024.0;
            }
        }
        const auto field =
            mesh::DenseFieldSnapshot::create(
                grid.topology,
                mesh::EntityKind::cell,
                field_components,
                std::move(values),
                field_metadata(workload));
        const double field_construction_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - field_begin).count();
        const auto rss_after_field = peak_rss_kib();

        const std::size_t topology_bytes =
            topology_payload_bytes(
                grid.topology);
        const std::size_t geometry_bytes =
            geometry_payload_bytes(grid);
        const std::size_t boundary_bytes =
            boundary_payload_bytes(
                grid.face_boundary);
        const std::size_t field_bytes =
            span_bytes(field.values());
        const std::size_t logical_payload =
            checked_add(
                checked_add(
                    topology_bytes,
                    geometry_bytes),
                checked_add(
                    boundary_bytes,
                    field_bytes));

        std::uint64_t checksum = 0U;
        checksum ^=
            traverse_topology_once(
                grid.topology);
        checksum ^=
            traverse_field_once(field);

        std::array<double, sample_count>
            topology_samples{};
        std::array<double, sample_count>
            field_samples{};

        for (std::size_t sample = 0U;
             sample < sample_count;
             ++sample) {
            const double topology_seconds =
                measure_seconds(
                    [&] {
                        return traverse_topology_once(
                            grid.topology);
                    },
                    topology_repetitions,
                    checksum);
            topology_samples[sample] =
                (static_cast<double>(
                     topology_bytes) *
                 static_cast<double>(
                     topology_repetitions)) /
                topology_seconds /
                gib;

            const double field_seconds =
                measure_seconds(
                    [&] {
                        return traverse_field_once(
                            field);
                    },
                    field_repetitions,
                    checksum);
            field_samples[sample] =
                (static_cast<double>(
                     field_bytes) *
                 static_cast<double>(
                     field_repetitions)) /
                field_seconds /
                gib;
        }

        const std::uint64_t rss_kib =
            peak_rss_kib();
        if (checksum == 0U ||
            logical_payload == 0U ||
            topology_bytes == 0U ||
            field_bytes == 0U) {
            throw std::runtime_error(
                "benchmark checksum/payload sanity failure");
        }

        std::cout
            << std::fixed
            << std::setprecision(6);
        std::cout
            << "benchmark_schema=mpmc.mesh.baseline.v1\n"
            << "workload=cartesian_3d_" << nx << 'x' << ny << 'x' << nz << "_field4\n"
            << "configuration_contract=passed\n"
            << "analytic_grid_invariants=passed\n"
            << std::setprecision(9)
            << "mesh_construction_seconds=" << mesh_seconds << '\n'
            << "field_construction_seconds=" << field_construction_seconds << '\n'
            << std::setprecision(6)
            << "peak_rss_source=" << peak_rss_source() << '\n'
            << "peak_rss_after_mesh_kib=" << rss_after_mesh << '\n'
            << "peak_rss_after_field_kib=" << rss_after_field << '\n'
            << "nx=" << nx << '\n'
            << "ny=" << ny << '\n'
            << "nz=" << nz << '\n'
            << "cells=" << cell_count << '\n'
            << "vertices="
            << grid.topology.entity_count(
                   mesh::EntityKind::vertex)
            << '\n'
            << "faces="
            << grid.topology.entity_count(
                   mesh::EntityKind::face)
            << '\n'
            << "field_components="
            << field_components << '\n'
            << "logical_topology_bytes="
            << topology_bytes << '\n'
            << "logical_geometry_bytes="
            << geometry_bytes << '\n'
            << "logical_boundary_bytes="
            << boundary_bytes << '\n'
            << "logical_field_bytes="
            << field_bytes << '\n'
            << "logical_total_bytes="
            << logical_payload << '\n'
            << "peak_rss_kib="
            << rss_kib << '\n'
            << "topology_repetitions="
            << topology_repetitions << '\n'
            << "field_repetitions="
            << field_repetitions << '\n';
        print_samples(
            "topology_samples_gib_s",
            topology_samples);
        std::cout
            << "topology_median_gib_s="
            << median(topology_samples)
            << '\n';
        print_samples(
            "field_samples_gib_s",
            field_samples);
        std::cout
            << "field_median_gib_s="
            << median(field_samples)
            << '\n'
            << "checksum="
            << checksum << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr
            << "[FAIL] mesh baseline benchmark: "
            << error.what()
            << '\n';
        return 1;
    }
}
