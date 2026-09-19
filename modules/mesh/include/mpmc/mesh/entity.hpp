#ifndef MPMC_MESH_ENTITY_HPP
#define MPMC_MESH_ENTITY_HPP

#include <compare>
#include <cstdint>

namespace mpmc::mesh {

/// Topological entity classes used by mesh connectivity and later field/DoF layouts.
enum class EntityKind : std::uint8_t {
    vertex = 0,
    edge = 1,
    face = 2,
    cell = 3,
};

/// Stable identity intended to survive local reordering and partitioning.
class GlobalEntityId {
public:
    using value_type = std::uint64_t;

    explicit constexpr GlobalEntityId(value_type value) noexcept : value_(value) {}

    [[nodiscard]] constexpr value_type value() const noexcept { return value_; }

    [[nodiscard]] friend constexpr auto operator<=>(const GlobalEntityId&,
                                                     const GlobalEntityId&) noexcept = default;

private:
    value_type value_;
};

/// Compact process-local index for hot-path topology access.
class LocalIndex {
public:
    using value_type = std::uint32_t;

    explicit constexpr LocalIndex(value_type value) noexcept : value_(value) {}

    [[nodiscard]] constexpr value_type value() const noexcept { return value_; }

    [[nodiscard]] friend constexpr auto operator<=>(const LocalIndex&,
                                                     const LocalIndex&) noexcept = default;

private:
    value_type value_;
};

} // namespace mpmc::mesh

#endif // MPMC_MESH_ENTITY_HPP
