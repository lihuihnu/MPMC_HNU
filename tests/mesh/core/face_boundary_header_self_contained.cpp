#include <mpmc/mesh/face_boundary.hpp>

#include <cstdint>
#include <span>
#include <type_traits>
#include <utility>

static_assert(sizeof(mpmc::mesh::FaceClassification) == sizeof(std::uint8_t));
static_assert(sizeof(mpmc::mesh::PhysicalTag) == sizeof(std::uint32_t));
static_assert(std::is_trivially_copyable_v<mpmc::mesh::PhysicalTag>);
static_assert(!std::is_default_constructible_v<mpmc::mesh::FaceBoundarySnapshot>);
static_assert(!std::is_copy_assignable_v<mpmc::mesh::FaceBoundarySnapshot>);
static_assert(std::is_same_v<
    decltype(std::declval<const mpmc::mesh::FaceBoundarySnapshot&>().physical_tags()),
    std::span<const mpmc::mesh::PhysicalTag>>);
