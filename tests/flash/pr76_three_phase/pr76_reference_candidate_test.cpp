#include <mpmc/flash/pr76_three_phase.hpp>

#include "synthetic_fixture.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace {
namespace fl = mpmc::flash;
using pr76_max3_test::Vec;

Vec log_ratio(const Vec& numerator, const Vec& denominator) {
    Vec value(numerator.size(), 0.0);
    for (std::size_t i = 0; i < value.size(); ++i) {
        value[i] = std::log(numerator[i]) - std::log(denominator[i]);
    }
    return value;
}

} // namespace

int main() {
    try {
        const auto model = pr76_max3_test::model();
        const auto& phases = pr76_max3_test::reference_phases();
        fl::Pr76ThreePhaseEvaluator evaluator(
            model,
            {fl::Pr76RootSide::lower_admissible,
             fl::Pr76RootSide::upper_admissible,
             fl::Pr76RootSide::lower_admissible});
        const auto result = fl::iterate_pt_three_phase(
            1.0e6, 250.0, pr76_max3_test::equal_feed(),
            log_ratio(phases[1], phases[0]),
            log_ratio(phases[2], phases[0]),
            {0.30, 0.30}, evaluator);
        if (result.status != fl::PtThreePhaseStatus::converged_candidate ||
            result.candidate() == nullptr ||
            result.candidate()->chemical_potential_norm >
                result.options.chemical_potential_tolerance) {
            throw std::runtime_error(
                "exact PR76 structural phase compositions did not close in fixed-three-phase primitive");
        }
        return 0;
    } catch (const std::exception&) {
        return 1;
    }
}
