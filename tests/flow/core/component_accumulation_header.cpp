#include <mpmc/flow/component_accumulation.hpp>

#include <cstddef>
#include <string_view>
#include <type_traits>
#include <utility>

static_assert(
    mpmc::flow::
        PoreVolumeComponentAccumulationSnapshot3P::
            convention ==
    mpmc::flow::
        pore_volume_component_accumulation_convention);

static_assert(
    std::is_same_v<
        decltype(
            std::declval<
                const mpmc::flow::
                    PoreVolumeComponentAccumulationLinearization3P&>()
                .d_component(
                    std::size_t{},
                    std::size_t{})),
        double>);

bool component_accumulation_header() {
    return true;
}
