#include <mpmc/ad/math.hpp>
#include <mpmc/ad/runtime_differentiate.hpp>
#include <mpmc/thermodynamics/cpa_residual_helmholtz.hpp>

#include "test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace ad = mpmc::ad;
namespace th = mpmc::thermodynamics;

void require(bool condition, const std::string& message) {
    if (!condition) { throw std::runtime_error(message); }
}

void require_abs(double actual, double expected, double tolerance,
                 const std::string& message) {
    if (!(std::abs(actual - expected) <= tolerance)) {
        std::ostringstream stream;
        stream << std::setprecision(17) << message
               << " actual=" << actual
               << " expected=" << expected
               << " delta=" << std::abs(actual - expected)
               << " tolerance=" << tolerance;
        throw std::runtime_error(stream.str());
    }
}

double roundoff(double scale) {
    constexpr double multiplier = 4096.0;
    return multiplier * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(scale));
}

struct HelmholtzValues {
    double cubic{};
    double association{};
    double total{};
};

struct TemperatureDerivatives {
    double cubic{};
    double association{};
    double total{};
};

double direct_association_value(
    std::span<const double> mole_numbers,
    const th::CpaAssociationResult& association) {
    double value = 0.0;
    for (const auto& site : association.sites) {
        const double x = site.unbonded_fraction;
        value += mole_numbers[site.component_index] *
            static_cast<double>(site.multiplicity) *
            (std::log(x) - 0.5 * x + 0.5);
    }
    return value;
}

std::vector<double> normalized_composition(
    std::span<const double> mole_numbers) {
    const double total_moles = std::accumulate(
        mole_numbers.begin(), mole_numbers.end(), 0.0);
    require(total_moles > 0.0 && std::isfinite(total_moles),
            "temperature-derivative regression requires positive total moles");

    std::vector<double> composition;
    composition.reserve(mole_numbers.size());
    for (const double value : mole_numbers) {
        require(value >= 0.0 && std::isfinite(value),
                "temperature-derivative regression requires finite nonnegative moles");
        composition.push_back(value / total_moles);
    }
    return composition;
}

HelmholtzValues direct_values_at_temperature(
    double temperature_k,
    double volume_m3,
    std::span<const double> mole_numbers,
    const th::CpaParameterSet& parameters) {
    require(temperature_k > 0.0 && std::isfinite(temperature_k),
            "temperature-derivative direct value requires positive temperature");
    require(volume_m3 > 0.0 && std::isfinite(volume_m3),
            "temperature-derivative direct value requires positive volume");
    require(mole_numbers.size() == parameters.size(),
            "temperature-derivative direct value dimension mismatch");

    const double total_moles = std::accumulate(
        mole_numbers.begin(), mole_numbers.end(), 0.0);
    require(total_moles > 0.0 && std::isfinite(total_moles),
            "temperature-derivative direct value requires positive total moles");
    const auto composition = normalized_composition(mole_numbers);
    const double rho = total_moles / volume_m3;
    const auto association = th::solve_cpa_association(
        temperature_k, rho, composition, parameters);
    require(association.converged(),
            "temperature-derivative finite-difference association did not converge");

    const double cubic = th::cpa_cubic_residual_helmholtz_reduced(
        temperature_k, volume_m3, mole_numbers, parameters);
    const double association_value = direct_association_value(
        mole_numbers, association);
    return {cubic, association_value, cubic + association_value};
}

TemperatureDerivatives ad_temperature_derivatives(
    double temperature_k,
    double volume_m3,
    std::span<const double> mole_numbers,
    const th::CpaParameterSet& parameters,
    const th::CpaAssociationResult& primal_association) {
    require(mole_numbers.size() == parameters.size(),
            "temperature-derivative AD dimension mismatch");

    const std::array<double, 1> inputs{{temperature_k}};
    using Number = ad::Dual<double, 4>;
    ad::RuntimeJacobianWorkspace<double, 4> workspace;
    const auto callback = [&](std::span<const Number> variables,
                              std::span<Number> outputs) {
        const Number& active_temperature = variables.front();
        const Number volume{volume_m3};
        std::vector<Number> active_moles;
        active_moles.reserve(mole_numbers.size());
        for (const double value : mole_numbers) {
            active_moles.emplace_back(value);
        }
        outputs[0] = th::cpa_cubic_residual_helmholtz_reduced(
            active_temperature, volume, std::span<const Number>{active_moles},
            parameters);
        outputs[1] = th::cpa_association_q_reduced(
            active_temperature, volume, std::span<const Number>{active_moles},
            parameters, primal_association);
        outputs[2] = th::cpa_residual_helmholtz_reduced(
            active_temperature, volume, std::span<const Number>{active_moles},
            parameters, primal_association);
    };

    const auto result = ad::value_and_jacobian_runtime<4>(
        callback, std::span<const double>{inputs}, 3U, workspace,
        {64U, 8U, 512U});
    require(result.input_count == 1U && result.output_count == 3U,
            "temperature-derivative AD returned unexpected Jacobian shape");

    return {
        result.jacobian[0U],
        result.jacobian[1U],
        result.jacobian[2U]};
}

TemperatureDerivatives central_difference(
    const HelmholtzValues& plus,
    const HelmholtzValues& minus,
    double step_k) {
    const double denominator = 2.0 * step_k;
    return {
        (plus.cubic - minus.cubic) / denominator,
        (plus.association - minus.association) / denominator,
        (plus.total - minus.total) / denominator};
}

void run_temperature_derivative_regression() {
    constexpr double temperature_k = cpa_physical_test::temperature_k;
    constexpr double rho = 200.0;
    constexpr double methanol_fraction = 0.35;
    constexpr double step_k = 0.005;
    constexpr double cubic_tolerance_per_k = 1.0e-10;
    constexpr double association_tolerance_per_k = 1.0e-9;
    constexpr double total_tolerance_per_k = 1.0e-9;
    constexpr double order_tolerance_per_k = 1.0e-12;
    constexpr std::array<bool, 2> swapped_orders{{false, true}};

    std::array<TemperatureDerivatives, 2> ad_by_order{};
    std::array<TemperatureDerivatives, 2> fd_by_order{};
    double max_abs_association_value = 0.0;
    double max_abs_cubic_delta = 0.0;
    double max_abs_association_delta = 0.0;
    double max_abs_total_delta = 0.0;

    for (std::size_t order_index = 0U;
         order_index < swapped_orders.size(); ++order_index) {
        const bool swapped = swapped_orders[order_index];
        const auto parameters = cpa_physical_test::parameters(swapped);
        const auto composition = cpa_physical_test::composition(
            methanol_fraction, swapped);
        const std::vector<double> mole_numbers(
            composition.begin(), composition.end());
        const double volume_m3 = 1.0 / rho;

        const auto base_association = th::solve_cpa_association(
            temperature_k, rho, composition, parameters);
        require(base_association.converged(),
                "temperature-derivative base association did not converge");
        max_abs_association_value = std::max(
            max_abs_association_value,
            std::abs(direct_association_value(
                mole_numbers, base_association)));

        const auto ad_values = ad_temperature_derivatives(
            temperature_k, volume_m3, mole_numbers,
            parameters, base_association);
        const auto plus = direct_values_at_temperature(
            temperature_k + step_k, volume_m3, mole_numbers, parameters);
        const auto minus = direct_values_at_temperature(
            temperature_k - step_k, volume_m3, mole_numbers, parameters);
        const auto fd_values = central_difference(plus, minus, step_k);

        ad_by_order[order_index] = ad_values;
        fd_by_order[order_index] = fd_values;

        const double cubic_delta = std::abs(
            ad_values.cubic - fd_values.cubic);
        const double association_delta = std::abs(
            ad_values.association - fd_values.association);
        const double total_delta = std::abs(
            ad_values.total - fd_values.total);
        max_abs_cubic_delta = std::max(
            max_abs_cubic_delta, cubic_delta);
        max_abs_association_delta = std::max(
            max_abs_association_delta, association_delta);
        max_abs_total_delta = std::max(
            max_abs_total_delta, total_delta);

        require_abs(
            ad_values.cubic, fd_values.cubic,
            cubic_tolerance_per_k + roundoff(ad_values.cubic) +
                roundoff(fd_values.cubic),
            "CPA Helmholtz cubic temperature derivative mismatch");
        require_abs(
            ad_values.association, fd_values.association,
            association_tolerance_per_k + roundoff(ad_values.association) +
                roundoff(fd_values.association),
            "CPA Helmholtz association temperature derivative mismatch");
        require_abs(
            ad_values.total, fd_values.total,
            total_tolerance_per_k + roundoff(ad_values.total) +
                roundoff(fd_values.total),
            "CPA Helmholtz total temperature derivative mismatch");
        require_abs(
            ad_values.total, ad_values.cubic + ad_values.association,
            roundoff(ad_values.total) +
                roundoff(ad_values.cubic + ad_values.association),
            "CPA Helmholtz AD temperature derivative decomposition mismatch");
        require_abs(
            fd_values.total, fd_values.cubic + fd_values.association,
            roundoff(fd_values.total) +
                roundoff(fd_values.cubic + fd_values.association),
            "CPA Helmholtz finite-difference derivative decomposition mismatch");
    }

    require(max_abs_association_value > 1.0e-8,
            "temperature-derivative regression did not exercise association");

    double max_abs_order_delta = 0.0;
    for (const auto* values : {&ad_by_order, &fd_by_order}) {
        for (const auto pair : std::array<std::pair<double, double>, 3>{{
                 {(*values)[0].cubic, (*values)[1].cubic},
                 {(*values)[0].association, (*values)[1].association},
                 {(*values)[0].total, (*values)[1].total}}}) {
            const double delta = std::abs(pair.first - pair.second);
            max_abs_order_delta = std::max(max_abs_order_delta, delta);
            require_abs(
                pair.first, pair.second,
                order_tolerance_per_k + roundoff(pair.first) +
                    roundoff(pair.second),
                "CPA Helmholtz temperature derivative changed under component permutation");
        }
    }

    std::cout << std::setprecision(17)
              << "CPA_HELMHOLTZ_TEMPERATURE_DERIVATIVE_OK"
              << " component_orders=2"
              << " rho=" << rho
              << " step_k=" << step_k
              << " max_abs_assoc_F=" << max_abs_association_value
              << " dF_dT_ad=" << ad_by_order[0].total
              << " dF_dT_fd=" << fd_by_order[0].total
              << " max_abs_d_cubic=" << max_abs_cubic_delta
              << " max_abs_d_assoc=" << max_abs_association_delta
              << " max_abs_d_total=" << max_abs_total_delta
              << " max_abs_order_delta=" << max_abs_order_delta
              << '\n';
}

} // namespace

int main() {
    try {
        run_temperature_derivative_regression();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
