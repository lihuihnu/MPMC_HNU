#include <mpmc/flash/pr76_three_phase.hpp>

#include "synthetic_fixture.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
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

        const auto property0 = evaluator(1.0e6, 250.0, phases[0], 0U);
        const auto property1 = evaluator(1.0e6, 250.0, phases[1], 1U);
        const auto property2 = evaluator(1.0e6, 250.0, phases[2], 2U);
        double anchor_mu_norm = 0.0;
        for (std::size_t i = 0; i < phases[0].size(); ++i) {
            const double mu0 = std::log(phases[0][i]) + property0.activity.ln_phi[i];
            const double mu1 = std::log(phases[1][i]) + property1.activity.ln_phi[i];
            const double mu2 = std::log(phases[2][i]) + property2.activity.ln_phi[i];
            anchor_mu_norm = std::max({
                anchor_mu_norm, std::abs(mu0 - mu1), std::abs(mu0 - mu2)});
        }

        const auto result = fl::iterate_pt_three_phase(
            1.0e6, 250.0, pr76_max3_test::equal_feed(),
            log_ratio(phases[1], phases[0]),
            log_ratio(phases[2], phases[0]),
            {0.30, 0.30}, evaluator);
        if (result.status != fl::PtThreePhaseStatus::converged_candidate ||
            result.candidate() == nullptr ||
            result.candidate()->chemical_potential_norm >
                result.options.chemical_potential_tolerance) {
            std::cerr << "anchor_mu_norm=" << anchor_mu_norm
                      << " root_branches=" << property0.activity.branch << ','
                      << property1.activity.branch << ','
                      << property2.activity.branch
                      << " z=" << property0.z << ',' << property1.z << ','
                      << property2.z
                      << " status=" << static_cast<int>(result.status)
                      << " iterations=" << result.iterations
                      << " evaluations=" << result.evaluations
                      << " final_norm="
                      << (result.point ? result.point->chemical_potential_norm : -1.0)
                      << " diagnostic=" << result.diagnostic << '\n';
            throw std::runtime_error(
                "exact PR76 structural phase compositions did not close in fixed-three-phase primitive");
        }
        if (anchor_mu_norm > 5.0e-8) {
            std::cerr << "anchor_mu_norm=" << anchor_mu_norm << '\n';
            throw std::runtime_error(
                "independent structural anchor does not satisfy PR76 chemical-potential equality");
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
