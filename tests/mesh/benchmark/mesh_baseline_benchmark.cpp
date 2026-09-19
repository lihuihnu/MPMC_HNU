#include <mpmc/mesh/cartesian_3d.hpp>
#include <mpmc/mesh/dense_field.hpp>
#include <mpmc/mesh/face_boundary.hpp>
#include <mpmc/mesh/face_geometry_3d.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <array>
#include <bit>
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
#include <utility>
#include <vector>

#if defined(__linux__)
#include <sys/resource.h>
#endif

namespace mesh = mpmc::mesh;

namespace {

constexpr std::size_t nx = 64U;
constexpr std::size_t ny = 64U;
constexpr std::size_t nz = 64U;
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
field_metadata() {
    return mesh::DenseFieldMetadata{
        "benchmark.cell.field4",
        "1",
        mesh::FieldSourceMetadata{
            mesh::FieldSourceKind::synthetic_test,
            "benchmark://mesh/fixed-64-cubed",
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
#else
    return 0U;
#endif
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

int main() {
    try {
        const auto x = axis(nx);
        const auto y = axis(ny);
        const auto z = axis(nz);
        const auto grid =
            mesh::make_cartesian_mesh_3d(
                x,
                y,
                z);

        const std::size_t cell_count =
            grid.topology.entity_count(
                mesh::EntityKind::cell);
        if (cell_count !=
            nx * ny * nz) {
            throw std::runtime_error(
                "fixed workload cell count drift");
        }

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
                field_metadata());

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
            << "workload=cartesian_3d_64x64x64_field4\n"
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
