#include <mpmc/flow_discretization/cell_source.hpp>

#include <string_view>

static_assert(
    mpmc::flow_discretization::
        CellSourceLinearization3D::convention ==
    mpmc::flow_discretization::
        cell_source_convention);

bool cell_source_header() {
    return
        mpmc::flow_discretization::
            cell_source_convention ==
        std::string_view{
            "flow_discretization/cell-source-linearization/v1"};
}
