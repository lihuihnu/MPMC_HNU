#include <mpmc/well_discretization/cell_source_adapter.hpp>

#include <string_view>

static_assert(
    mpmc::well_discretization::
        WellConnectionCellSourceAdapterResult3P::
            convention ==
    mpmc::well_discretization::
        well_connection_cell_source_adapter_convention);

int main() {
    return
        mpmc::well_discretization::
            well_connection_cell_source_adapter_convention ==
        std::string_view{
            "well-discretization/connection-to-cell-source/sign-bridge/v1"}
        ? 0
        : 1;
}
