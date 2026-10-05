#include <mpmc/mesh/mesh_exchange_group.hpp>

static_assert(sizeof(mpmc::mesh::MeshExchangeGroup) > 0U);

bool mesh_exchange_group_header_self_contained() {
    const mpmc::mesh::MeshExchangeGroup group{mpmc::mesh::EntityKind::cell, 1U, "probe", {}};
    return group.tag == 1U;
}
