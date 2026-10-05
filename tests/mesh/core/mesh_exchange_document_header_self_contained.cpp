#include <mpmc/mesh/mesh_exchange_document.hpp>

static_assert(sizeof(mpmc::mesh::MeshExchangeDocument) > 0U);

bool mesh_exchange_document_header_self_contained() {
    const mpmc::mesh::ConversionReport report{mpmc::mesh::MeshExchangeFormat::vtu_ascii};
    return report.target_format() == mpmc::mesh::MeshExchangeFormat::vtu_ascii;
}
