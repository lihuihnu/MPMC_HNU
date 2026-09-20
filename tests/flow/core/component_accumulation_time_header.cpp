#include <mpmc/flow/component_accumulation_time.hpp>

#include <cstddef>
#include <string_view>
#include <type_traits>
#include <utility>

static_assert(
    mpmc::flow::
        BackwardEulerComponentAccumulationResidual3P::
            convention ==
    mpmc::flow::
        backward_euler_component_accumulation_convention);

static_assert(
    std::is_same_v<
        decltype(
            std::declval<
                const mpmc::flow::
                    BackwardEulerComponentAccumulationResidual3P&>()
                .component_row_identity(
                    std::size_t{})),
        std::string_view>);

bool component_accumulation_time_header() {
    return true;
}
