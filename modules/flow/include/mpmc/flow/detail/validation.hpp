#ifndef MPMC_FLOW_DETAIL_VALIDATION_HPP
#define MPMC_FLOW_DETAIL_VALIDATION_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <string>

namespace mpmc::flow::validation_detail {

// Share the comparison algorithm, not a universal scientific tolerance.
// Each caller keeps its existing epsilon multiplier and optional closure scale.
// This policy has a unit scale floor; rate comparisons with a different floor
// deliberately retain their own implementation.
template <std::size_t EpsilonMultiplier>
[[nodiscard]] inline bool near_roundoff(
    double first, double second, double extra_scale = 0.0) {
    static_assert(EpsilonMultiplier > 0U);
    if (!std::isfinite(first) || !std::isfinite(second) ||
        !std::isfinite(extra_scale) || extra_scale < 0.0) {
        return false;
    }
    const double scale =
        std::max({1.0, std::abs(first), std::abs(second), extra_scale});
    return std::abs(first - second) <=
        static_cast<double>(EpsilonMultiplier) *
            std::numeric_limits<double>::epsilon() * scale;
}

// Compare borrowed identity sequences without constructing a vector of strings.
// Order and length are part of the identity; this is not set membership.
[[nodiscard]] inline bool same_component_ids(
    std::span<const std::string> first, std::span<const std::string> second) {
    return std::equal(first.begin(), first.end(), second.begin(), second.end());
}

// Layout descriptors include the active phase count and dependent-component
// pivots. Matching only the number of unknowns would accept a different chart.
template <typename Layout>
[[nodiscard]] inline bool same_layout(const Layout& first, const Layout& second) {
    return first.component_count() == second.component_count() &&
        first.phase_count() == second.phase_count() &&
        first.unknown_count() == second.unknown_count() &&
        first.composition_pivot().dependent_components() ==
            second.composition_pivot().dependent_components();
}

// Identity comparison is distinct from input validation: callers retain their
// positive-support checks and diagnostics. Inactive sidecars must remain exactly
// zero/empty even when active values are compared within a roundoff allowance.
template <typename Identity, typename Near>
[[nodiscard]] inline bool same_state_identity(
    const Identity& first, const Identity& second, Near near) {
    if (!same_layout(first.layout, second.layout) ||
        !same_component_ids(first.component_ids, second.component_ids) ||
        !near(first.reference_pressure_pa, second.reference_pressure_pa) ||
        !near(first.temperature_k, second.temperature_k)) {
        return false;
    }
    for (std::size_t phase = 0U; phase < first.layout.phase_count(); ++phase) {
        if (!near(first.saturation[phase], second.saturation[phase]) ||
            first.phase_composition[phase].size() !=
                second.phase_composition[phase].size()) {
            return false;
        }
        for (std::size_t component = 0U;
             component < first.phase_composition[phase].size(); ++component) {
            if (!near(first.phase_composition[phase][component],
                      second.phase_composition[phase][component])) {
                return false;
            }
        }
    }
    for (std::size_t phase = first.layout.phase_count();
         phase < first.saturation.size(); ++phase) {
        if (first.saturation[phase] != 0.0 || second.saturation[phase] != 0.0 ||
            !first.phase_composition[phase].empty() ||
            !second.phase_composition[phase].empty()) {
            return false;
        }
    }
    return true;
}

} // namespace mpmc::flow::validation_detail

#endif // MPMC_FLOW_DETAIL_VALIDATION_HPP
