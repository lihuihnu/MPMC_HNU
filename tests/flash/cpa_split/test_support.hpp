#ifndef MPMC_TEST_CPA_SPLIT_SUPPORT_HPP
#define MPMC_TEST_CPA_SPLIT_SUPPORT_HPP

#include "../cpa_stability/test_support.hpp"

#include <array>
#include <vector>

namespace cpa_split_test {
using Vec = std::vector<double>;

inline constexpr double pressure_pa = 1.0e6;
inline constexpr double temperature_k = 160.0;

// Independent offline chemical-potential equality solve of the explicitly
// synthetic non-associating CPA/SRK-limit binary from cpa_stability/test_support.hpp.
// These are structural numerical anchors only, not experimental phase data.
inline const std::array<Vec, 2>& coexistence_phases() {
    static const std::array<Vec, 2> values{{
        {0.9374784092011380, 0.0625215907988620}, // upper-density branch
        {0.5071007718530185, 0.4928992281469815}  // lower-density branch
    }};
    return values;
}

inline Vec feed() {
    const auto& phases = coexistence_phases();
    return {
        0.5 * (phases[0][0] + phases[1][0]),
        0.5 * (phases[0][1] + phases[1][1])};
}

inline std::vector<Vec> starts(bool swapped = false) {
    std::vector<Vec> values;
    values.reserve(coexistence_phases().size());
    for (const auto& phase : coexistence_phases()) {
        if (swapped) {
            values.push_back({phase[1], phase[0]});
        } else {
            values.push_back(phase);
        }
    }
    return values;
}

inline Vec feed(bool swapped) {
    const auto value = feed();
    return swapped ? Vec{value[1], value[0]} : value;
}

} // namespace cpa_split_test

#endif // MPMC_TEST_CPA_SPLIT_SUPPORT_HPP
