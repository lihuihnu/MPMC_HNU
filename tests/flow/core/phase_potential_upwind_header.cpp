#include <mpmc/flow/phase_potential_upwind.hpp>

#include <string_view>
#include <type_traits>
#include <utility>

static_assert(
    mpmc::flow::
        TwoCellPhasePotentialUpwindLinearization3P::
            convention ==
    mpmc::flow::
        two_cell_phase_potential_convention);

static_assert(
    std::is_same_v<
        decltype(
            (std::declval<
                const mpmc::flow::
                    TwoCellPhasePotentialEntry3P&>()
                 .phase_potential_difference_pa)),
        const double&>);

bool phase_potential_upwind_header() {
    return true;
}
