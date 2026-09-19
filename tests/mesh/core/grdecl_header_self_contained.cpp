#include <mpmc/mesh/grdecl.hpp>

#include <type_traits>

static_assert(std::is_move_constructible_v<mpmc::mesh::GrdeclImportResult>);
static_assert(!std::is_default_constructible_v<mpmc::mesh::GrdeclImportResult>);
