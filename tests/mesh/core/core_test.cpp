#include <mpmc/mesh/csr_adjacency.hpp>
#include <mpmc/mesh/entity.hpp>

#include <array>
#include <cstdint>
#include <exception>
#include <iostream>
#include <source_location>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace {
namespace mesh = mpmc::mesh;

void require(bool condition,
             std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " + std::string(message));
    }
}

template <typename Exception, typename Function>
void expect_throw(Function&& function) {
    bool caught = false;
    try {
        std::forward<Function>(function)();
    } catch (const Exception&) {
        caught = true;
    }
    require(caught, "expected exception was not thrown");
}

mesh::CsrAdjacency two_cell_quad_cell_to_vertex() {
    //  v3 ---- v4 ---- v5
    //   |  c0  |  c1  |
    //  v0 ---- v1 ---- v2
    //
    // Connectivity is synthetic software-test data, not a physical mesh dataset.
    std::vector<mesh::CsrAdjacency::Offset> offsets{0U, 4U, 8U};
    std::vector<mesh::LocalIndex> vertices{
        mesh::LocalIndex{0U}, mesh::LocalIndex{1U}, mesh::LocalIndex{4U},
        mesh::LocalIndex{3U}, mesh::LocalIndex{1U}, mesh::LocalIndex{2U},
        mesh::LocalIndex{5U}, mesh::LocalIndex{4U},
    };
    return mesh::CsrAdjacency{mesh::EntityKind::cell, mesh::EntityKind::vertex, 6U,
                              std::move(offsets), std::move(vertices)};
}

void strong_indices() {
    static_assert(sizeof(mesh::LocalIndex) == sizeof(std::uint32_t));
    static_assert(sizeof(mesh::GlobalEntityId) == sizeof(std::uint64_t));
    static_assert(std::is_trivially_copyable_v<mesh::LocalIndex>);
    static_assert(std::is_trivially_copyable_v<mesh::GlobalEntityId>);
    static_assert(!std::is_default_constructible_v<mesh::LocalIndex>);
    static_assert(!std::is_default_constructible_v<mesh::GlobalEntityId>);
    static_assert(!std::is_convertible_v<mesh::LocalIndex, std::uint32_t>);
    static_assert(!std::is_convertible_v<mesh::GlobalEntityId, std::uint64_t>);
    static_assert(!std::is_same_v<mesh::LocalIndex, mesh::GlobalEntityId>);

    constexpr mesh::LocalIndex local{17U};
    constexpr mesh::GlobalEntityId global{9000000001ULL};
    static_assert(local.value() == 17U);
    static_assert(global.value() == 9000000001ULL);
    require(mesh::EntityKind::cell != mesh::EntityKind::face, "entity kinds collapsed");
}

void two_cell_quad() {
    const auto adjacency = two_cell_quad_cell_to_vertex();
    require(adjacency.source_kind() == mesh::EntityKind::cell, "source kind lost");
    require(adjacency.target_kind() == mesh::EntityKind::vertex, "target kind lost");
    require(adjacency.source_count() == 2U, "wrong source count");
    require(adjacency.target_count() == 6U, "wrong target count");
    require(adjacency.entry_count() == 8U, "wrong entry count");

    const std::array<std::uint32_t, 4> expected0{0U, 1U, 4U, 3U};
    const std::array<std::uint32_t, 4> expected1{1U, 2U, 5U, 4U};
    const auto row0 = adjacency.adjacent(mesh::LocalIndex{0U});
    const auto row1 = adjacency.adjacent(mesh::LocalIndex{1U});
    require(row0.size() == expected0.size() && row1.size() == expected1.size(), "row width");
    for (std::size_t i = 0; i < expected0.size(); ++i) {
        require(row0[i].value() == expected0[i], "cell 0 connectivity mismatch");
        require(row1[i].value() == expected1[i], "cell 1 connectivity mismatch");
    }

    const auto offsets = adjacency.offsets();
    const auto indices = adjacency.indices();
    require(offsets.size() == 3U && offsets[0] == 0U && offsets[1] == 4U &&
                offsets[2] == 8U,
            "CSR offsets mismatch");
    require(indices.size() == 8U, "CSR values mismatch");
    require(row0.data() == indices.data(),
            "row view must alias contiguous value storage");
    require(row1.data() == indices.data() + 4,
            "second row must be contiguous without per-row storage");

    expect_throw<std::out_of_range>(
        [&] { (void)adjacency.adjacent(mesh::LocalIndex{2U}); });
}

void invalid_csr() {
    using Offset = mesh::CsrAdjacency::Offset;
    const auto make = [](std::size_t target_count,
                         std::vector<Offset> offsets,
                         std::vector<mesh::LocalIndex> indices) {
        return mesh::CsrAdjacency{mesh::EntityKind::cell, mesh::EntityKind::vertex,
                                  target_count, std::move(offsets), std::move(indices)};
    };

    expect_throw<std::invalid_argument>([&] { (void)make(0U, {}, {}); });
    expect_throw<std::invalid_argument>([&] { (void)make(2U, {1U}, {}); });
    expect_throw<std::invalid_argument>([&] {
        (void)make(3U, {0U, 2U, 1U},
                   {mesh::LocalIndex{0U}, mesh::LocalIndex{1U}});
    });
    expect_throw<std::invalid_argument>([&] {
        (void)make(3U, {0U, 1U},
                   {mesh::LocalIndex{0U}, mesh::LocalIndex{1U}});
    });
    expect_throw<std::out_of_range>([&] {
        (void)make(3U, {0U, 1U}, {mesh::LocalIndex{3U}});
    });

    const auto empty = make(0U, {0U}, {});
    require(empty.source_count() == 0U && empty.entry_count() == 0U,
            "empty relation rejected");
}

void headers() {
    static_assert(
        std::is_same_v<decltype(std::declval<const mesh::CsrAdjacency&>().indices()),
                       std::span<const mesh::LocalIndex>>);
    static_assert(
        std::is_same_v<decltype(std::declval<const mesh::CsrAdjacency&>().offsets()),
                       std::span<const mesh::CsrAdjacency::Offset>>);
    require(true, "headers");
}

} // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "provide one named mesh core test");
        const std::string_view name{argv[1]};
        if (name == "strong_indices") { strong_indices(); }
        else if (name == "two_cell_quad") { two_cell_quad(); }
        else if (name == "invalid_csr") { invalid_csr(); }
        else if (name == "headers") { headers(); }
        else { throw std::invalid_argument("unknown mesh core test"); }
        std::cout << "[PASS] " << name << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
