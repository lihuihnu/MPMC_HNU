#include <mpmc/flow_discretization/tpfa_component_molar_flux.hpp>

#include <string_view>
#include <type_traits>
#include <utility>

static_assert(
    mpmc::flow_discretization::
        MaterializedTpfaInternalFaceComponentMolarFluxLinearization3D::
            convention ==
    mpmc::flow_discretization::
        tpfa_component_molar_flux_convention);

static_assert(
    std::is_same_v<
        decltype(
            (std::declval<
                const mpmc::flow_discretization::
                    MaterializedTpfaInternalFaceComponentMolarFluxLinearization3D&>()
                 .total_molar_flux_mol_per_s)),
        const double&>);

bool tpfa_component_molar_flux_header() {
    return true;
}
