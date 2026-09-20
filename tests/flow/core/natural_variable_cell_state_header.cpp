#include <mpmc/flow/natural_variable_cell_state.hpp>

#include <cstddef>
#include <optional>
#include <type_traits>
#include <utility>

static_assert(
    std::is_same_v<
        decltype(
            std::declval<
                const mpmc::flow::
                    NaturalVariableLayout3P&>()
                .independent_saturation_unknown_index(
                    mpmc::flow::PhaseSlot3::phase0)),
        std::optional<std::size_t>>);

static_assert(
    std::is_same_v<
        decltype(
            std::declval<
                const mpmc::flow::
                    NaturalVariableLayout3P&>()
                .composition_unknown_identity(
                    std::size_t{})),
        std::optional<
            mpmc::flow::
                NaturalVariableCompositionUnknownIdentity3P>>);

bool natural_variable_cell_state_header() {
    return true;
}
