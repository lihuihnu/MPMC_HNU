#include <mpmc/flow_discretization/tpfa_phase_darcy_flux.hpp>

#include <string_view>
#include <type_traits>
#include <utility>

static_assert(
    mpmc::flow_discretization::
        MaterializedTpfaInternalFacePhaseDarcyFluxLinearization3D::
            convention ==
    mpmc::flow_discretization::
        materialized_tpfa_phase_flux_convention);

static_assert(
    std::is_same_v<
        decltype(
            (std::declval<
                const mpmc::flow_discretization::
                    MaterializedTpfaPhaseDarcyFluxEntry3D&>()
                 .volumetric_flux_m3_per_s)),
        const double&>);

bool tpfa_phase_darcy_flux_header() {
    return true;
}
