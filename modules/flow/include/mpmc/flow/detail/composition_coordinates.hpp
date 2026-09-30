#ifndef MPMC_FLOW_DETAIL_COMPOSITION_COORDINATES_HPP
#define MPMC_FLOW_DETAIL_COMPOSITION_COORDINATES_HPP

#include <mpmc/flow/natural_variable_cell_state.hpp>

#include <cstddef>

namespace mpmc::flow::composition_coordinate_detail {

/// One composition block omits its frozen dependent component. Independent
/// columns follow canonical component order. Thus dx_i/dq is +1 for its own
/// column, -1 for the reconstructed component, and zero elsewhere.
///
/// The caller owns index validation and diagnostics. This allocation-free
/// lookup accepts a validated, representable block (component_count >= 2).
[[nodiscard]] inline double block_derivative(
    std::size_t component_count,
    std::size_t dependent_component,
    std::size_t first_column,
    std::size_t component,
    std::size_t column) noexcept {
    if (column < first_column ||
        column - first_column >= component_count - 1U) {
        return 0.0;
    }
    const std::size_t rank = column - first_column;
    const std::size_t independent_component =
        rank < dependent_component ? rank : rank + 1U;
    if (component == independent_component) {
        return 1.0;
    }
    return component == dependent_component ? -1.0 : 0.0;
}

/// Descriptor and 3P layouts already expose an inverse column map. Keep that
/// map's checked-column behavior, including its original exception messages.
template <class Layout>
[[nodiscard]] inline double derivative(
    const Layout& layout,
    std::size_t phase,
    std::size_t component,
    std::size_t column) {
    const auto identity = layout.composition_unknown_identity(column);
    if (!identity || static_cast<std::size_t>(identity->phase) != phase) {
        return 0.0;
    }
    if (component == identity->component) {
        return 1.0;
    }
    return component == layout.dependent_composition_component(identity->phase)
        ? -1.0 : 0.0;
}

} // namespace mpmc::flow::composition_coordinate_detail

#endif // MPMC_FLOW_DETAIL_COMPOSITION_COORDINATES_HPP
