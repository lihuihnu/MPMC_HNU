#include <mpmc/flow/phase_transport.hpp>

#include <string_view>
#include <type_traits>
#include <utility>

static_assert(
    mpmc::flow::
        PhaseTransportPropertyNaturalVariableLinearization3P::
            convention ==
    mpmc::flow::
        phase_transport_linearization_convention);

static_assert(
    mpmc::flow::
        LocalPhaseMobilityLinearization3P::
            convention ==
    mpmc::flow::
        local_phase_mobility_convention);

static_assert(
    std::is_same_v<
        decltype(
            (std::declval<
                const mpmc::flow::
                    LocalPhaseMobilityLinearization3P&>()
                 .mobility_per_pa_s)),
        const std::array<double, 3>&>);

bool phase_transport_header() {
    return true;
}
