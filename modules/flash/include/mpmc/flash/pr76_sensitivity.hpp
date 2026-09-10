#ifndef MPMC_FLASH_PR76_SENSITIVITY_HPP
#define MPMC_FLASH_PR76_SENSITIVITY_HPP

#include <mpmc/flash/detail/pr76_sensitivity_local.hpp>
#include <mpmc/flash/pr76_split.hpp>
#include <mpmc/flash/pt_sensitivity.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace mpmc::flash {

struct Pr76PtVleSensitivityResult {
    PtVleSensitivityResult sensitivity;
    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string model_profile{thermodynamics::pr76_profile};
    std::string phase_convention{thermodynamics::pr76_pt_convention};
    thermodynamics::Pr76RootOptions root_options;
};

namespace detail {

inline void pt_sensitivity_set_failure(
    PtVleSensitivityResult& result, PtSensitivityStatus status,
    const std::string& diagnostic) {
    result.status = status;
    result.diagnostic = diagnostic;
}

[[nodiscard]] inline bool pt_sensitivity_vector_finite(
    const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

} // namespace detail

/// Differentiate one already accepted PR76 VLE solution by the local implicit
/// equilibrium equations. No SSI, RR bisection, line-search or TPD iteration is
/// differentiated. Columns are q=(p_Pa,T_K,z_0,...,z_{N-2}), with the final feed
/// component dependent. Success is local to the accepted phase set and stored
/// liquid/vapor PR76 simple-root branches.
[[nodiscard]] inline Pr76PtVleSensitivityResult differentiate_pr76_pt_vle(
    const Pr76PtSplitResult& split, const Pr76VleEvaluator& evaluator,
    PtSensitivityOptions options = {}) {
    if (!std::isfinite(options.minimum_derivative_phase_fraction) ||
        !(options.minimum_derivative_phase_fraction > 0.0) ||
        !(options.minimum_derivative_phase_fraction < 0.5) ||
        !std::isfinite(options.minimum_reciprocal_condition) ||
        options.minimum_reciprocal_condition < 0.0 ||
        options.minimum_reciprocal_condition >= 1.0 ||
        options.max_components < 2) {
        throw std::invalid_argument("PT sensitivity: invalid options");
    }

    const auto& model = evaluator.model();
    const auto& accepted_feed = split.solution.initial_stability.feed;
    const double pressure_pa = split.solution.initial_stability.pressure_pa;
    const double temperature_k = split.solution.initial_stability.temperature_k;
    const std::size_t n = model.size();
    if (n < 2 || n > options.max_components || accepted_feed.size() != n) {
        throw std::length_error(
            "PT sensitivity: component count unsupported or split/model mismatch");
    }
    if (!std::isfinite(pressure_pa) || !(pressure_pa > 0.0) ||
        !std::isfinite(temperature_k) || !(temperature_k > 0.0)) {
        throw std::domain_error("PT sensitivity: accepted split has invalid p or T");
    }
    (void)detail::stability_check_composition(accepted_feed);

    Pr76PtVleSensitivityResult result;
    result.dataset_id = split.dataset_id;
    result.revision = split.revision;
    result.component_ids = split.component_ids;
    result.model_profile = split.model_profile;
    result.phase_convention = split.phase_convention;
    result.root_options = split.root_options;
    auto& sensitivity = result.sensitivity;
    sensitivity.component_count = n;
    sensitivity.input_count = n + 1;
    sensitivity.pressure_pa = pressure_pa;
    sensitivity.temperature_k = temperature_k;
    sensitivity.feed = accepted_feed;

    const auto& parameters = model.parameters();
    if (split.dataset_id != parameters.dataset_id() ||
        split.revision != parameters.revision() ||
        split.component_ids.size() != n ||
        split.model_profile != thermodynamics::pr76_profile ||
        split.phase_convention != thermodynamics::pr76_pt_convention ||
        split.root_options.max_iterations != evaluator.root_options().max_iterations) {
        throw std::invalid_argument(
            "PT sensitivity: split result and evaluator metadata do not match");
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (split.component_ids[i] != parameters.components().items()[i].id) {
            throw std::invalid_argument(
                "PT sensitivity: ordered component snapshot does not match");
        }
    }

    if (split.solution.status != PtSplitStatus::two_phase_no_instability_found ||
        !split.solution.candidate() || !split.solution.final_stability ||
        split.solution.final_stability->status != StabilityStatus::no_instability_found) {
        detail::pt_sensitivity_set_failure(
            sensitivity, PtSensitivityStatus::solution_not_accepted,
            "PT sensitivity: requires an accepted two-phase flash result");
        return result;
    }
    for (double value : accepted_feed) {
        if (!(value > 0.0)) {
            detail::pt_sensitivity_set_failure(
                sensitivity, PtSensitivityStatus::unsupported_feed_support,
                "PT sensitivity: zero-feed support is not implemented in v1");
            return result;
        }
    }

    const auto& candidate = *split.solution.candidate();
    if (candidate.log_k.size() != n || candidate.fractions.liquid.size() != n ||
        candidate.fractions.vapor.size() != n ||
        candidate.liquid.activity.ln_phi.size() != n ||
        candidate.vapor.activity.ln_phi.size() != n) {
        throw std::invalid_argument("PT sensitivity: split candidate dimension mismatch");
    }

    double log_k_contrast = 0.0;
    for (double value : candidate.log_k) {
        log_k_contrast = std::max(log_k_contrast, std::abs(value));
    }
    if (!(log_k_contrast > split.solution.options.iteration.log_k_separation) ||
        !(candidate.vapor.z - candidate.liquid.z >
          split.solution.options.iteration.relative_z_separation *
              std::max(candidate.vapor.z, candidate.liquid.z))) {
        detail::pt_sensitivity_set_failure(
            sensitivity, PtSensitivityStatus::phase_boundary,
            "PT sensitivity: phases are inside the local separation guard");
        return result;
    }

    const double beta = candidate.fractions.vapor_fraction;
    const double phase_guard = std::max(
        options.minimum_derivative_phase_fraction,
        split.solution.options.iteration.minimum_phase_fraction);
    if (!(beta > phase_guard) || !(1.0 - beta > phase_guard)) {
        detail::pt_sensitivity_set_failure(
            sensitivity, PtSensitivityStatus::phase_boundary,
            "PT sensitivity: phase fraction is inside the derivative boundary guard");
        return result;
    }
    if (!(candidate.fugacity_norm <=
          split.solution.options.iteration.fugacity_tolerance)) {
        detail::pt_sensitivity_set_failure(
            sensitivity, PtSensitivityStatus::solution_not_accepted,
            "PT sensitivity: candidate no longer satisfies fugacity tolerance");
        return result;
    }

    for (std::size_t i = 0; i < n; ++i) {
        if (!(candidate.fractions.liquid[i] > 0.0) ||
            !(candidate.fractions.vapor[i] > 0.0) ||
            !std::isfinite(candidate.log_k[i])) {
            detail::pt_sensitivity_set_failure(
                sensitivity, PtSensitivityStatus::unsupported_feed_support,
                "PT sensitivity: positive phase compositions are required");
            return result;
        }
        const double recovered =
            std::fma(beta, candidate.fractions.vapor[i],
                     (1.0 - beta) * candidate.fractions.liquid[i]);
        const double error = std::abs(recovered - accepted_feed[i]);
        if (error > split.solution.options.iteration.mass_absolute_tolerance ||
            error / accepted_feed[i] >
                split.solution.options.iteration.mass_relative_tolerance) {
            throw std::invalid_argument(
                "PT sensitivity: accepted feed does not match the split candidate");
        }
    }

    detail::Pr76SensitivityLocalSystem local_system;
    try {
        local_system =
            detail::pr76_sensitivity_local_system(split, evaluator, candidate);
    } catch (const detail::PtSensitivityFailure& error) {
        detail::pt_sensitivity_set_failure(sensitivity, error.status(), error.what());
        return result;
    } catch (const thermodynamics::Pr76PhaseError& error) {
        detail::pt_sensitivity_set_failure(
            sensitivity, PtSensitivityStatus::property_failure, error.what());
        return result;
    } catch (const std::range_error& error) {
        detail::pt_sensitivity_set_failure(
            sensitivity, PtSensitivityStatus::property_failure, error.what());
        return result;
    } catch (const std::domain_error& error) {
        detail::pt_sensitivity_set_failure(
            sensitivity, PtSensitivityStatus::arithmetic_failure, error.what());
        return result;
    } catch (const std::runtime_error& error) {
        detail::pt_sensitivity_set_failure(
            sensitivity, PtSensitivityStatus::arithmetic_failure, error.what());
        return result;
    }

    const auto& local = local_system.values;
    sensitivity.equilibrium_residual_norm = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        sensitivity.equilibrium_residual_norm = std::max(
            sensitivity.equilibrium_residual_norm, std::abs(local.values[i]));
    }
    const double arithmetic_scale =
        8192.0 * std::numeric_limits<double>::epsilon() *
        static_cast<double>(n);
    const double residual_limit =
        split.solution.options.iteration.fugacity_tolerance + arithmetic_scale;
    if (sensitivity.equilibrium_residual_norm > residual_limit) {
        detail::pt_sensitivity_set_failure(
            sensitivity, PtSensitivityStatus::solution_not_accepted,
            "PT sensitivity: stored split data do not reproduce equilibrium");
        return result;
    }

    const std::size_t q_count = n + 1;
    const std::size_t variable_count = local_system.variable_count;
    const auto q_variable = [n](std::size_t column) { return n + column; };
    std::vector<double> a(n * n);
    std::vector<double> b(n * q_count);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t column = 0; column < n; ++column) {
            a[row * n + column] =
                local.jacobian[row * variable_count + column];
        }
        for (std::size_t column = 0; column < q_count; ++column) {
            b[row * q_count + column] =
                local.jacobian[row * variable_count + q_variable(column)];
        }
    }

    detail::PtSensitivityLu factor;
    if (!detail::pt_sensitivity_factor(a, n, factor)) {
        detail::pt_sensitivity_set_failure(
            sensitivity, PtSensitivityStatus::ill_conditioned_equilibrium,
            "PT sensitivity: dF/dlogK cannot be stably factorized");
        return result;
    }
    sensitivity.equilibrium_jacobian_rcond =
        detail::pt_sensitivity_rcond(factor);
    const double required_rcond = std::max(
        256.0 * std::numeric_limits<double>::epsilon() *
            static_cast<double>(n),
        options.minimum_reciprocal_condition);
    if (!std::isfinite(sensitivity.equilibrium_jacobian_rcond) ||
        !(sensitivity.equilibrium_jacobian_rcond > required_rcond)) {
        detail::pt_sensitivity_set_failure(
            sensitivity, PtSensitivityStatus::ill_conditioned_equilibrium,
            "PT sensitivity: local equilibrium Jacobian is ill-conditioned");
        return result;
    }

    std::vector<double> d_log_k(n * q_count);
    for (std::size_t column = 0; column < q_count; ++column) {
        std::vector<double> rhs(n);
        for (std::size_t row = 0; row < n; ++row) {
            rhs[row] = -b[row * q_count + column];
        }
        if (!detail::pt_sensitivity_solve(factor, rhs)) {
            detail::pt_sensitivity_set_failure(
                sensitivity, PtSensitivityStatus::ill_conditioned_equilibrium,
                "PT sensitivity: implicit linear solve failed");
            return result;
        }
        for (std::size_t row = 0; row < n; ++row) {
            d_log_k[row * q_count + column] = rhs[row];
        }
    }

    double backward_error = 0.0;
    for (std::size_t column = 0; column < q_count; ++column) {
        for (std::size_t row = 0; row < n; ++row) {
            double residual = b[row * q_count + column];
            double scale = std::abs(residual);
            for (std::size_t j = 0; j < n; ++j) {
                const double term =
                    a[row * n + j] * d_log_k[j * q_count + column];
                residual += term;
                scale += std::abs(term);
            }
            backward_error = std::max(
                backward_error,
                std::abs(residual) / std::max(1.0, scale));
        }
    }
    sensitivity.linear_solve_backward_error = backward_error;
    if (!std::isfinite(backward_error) || backward_error > arithmetic_scale) {
        detail::pt_sensitivity_set_failure(
            sensitivity, PtSensitivityStatus::arithmetic_failure,
            "PT sensitivity: implicit solve failed backward-error check");
        return result;
    }

    sensitivity.log_k_jacobian = d_log_k;
    sensitivity.vapor_fraction_gradient.assign(q_count, 0.0);
    sensitivity.liquid_jacobian.assign(n * q_count, 0.0);
    sensitivity.vapor_jacobian.assign(n * q_count, 0.0);
    const auto total_derivative =
        [&](std::size_t output, std::size_t q_column) {
            double value = local.jacobian[
                output * variable_count + q_variable(q_column)];
            for (std::size_t j = 0; j < n; ++j) {
                value += local.jacobian[output * variable_count + j] *
                         d_log_k[j * q_count + q_column];
            }
            return value;
        };

    for (std::size_t column = 0; column < q_count; ++column) {
        sensitivity.vapor_fraction_gradient[column] =
            total_derivative(n, column);
        for (std::size_t i = 0; i < n; ++i) {
            sensitivity.liquid_jacobian[i * q_count + column] =
                total_derivative(n + 1 + i, column);
            sensitivity.vapor_jacobian[i * q_count + column] =
                total_derivative(2 * n + 1 + i, column);
        }
    }

    if (!detail::pt_sensitivity_vector_finite(sensitivity.log_k_jacobian) ||
        !detail::pt_sensitivity_vector_finite(sensitivity.liquid_jacobian) ||
        !detail::pt_sensitivity_vector_finite(sensitivity.vapor_jacobian) ||
        !detail::pt_sensitivity_vector_finite(
            sensitivity.vapor_fraction_gradient)) {
        detail::pt_sensitivity_set_failure(
            sensitivity, PtSensitivityStatus::arithmetic_failure,
            "PT sensitivity: nonfinite total derivative");
        return result;
    }

    const double beta_base = local.values[n];
    for (std::size_t column = 0; column < q_count; ++column) {
        double liquid_sum = 0.0;
        double vapor_sum = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            liquid_sum += sensitivity.d_liquid(i, column);
            vapor_sum += sensitivity.d_vapor(i, column);
        }
        if (std::abs(liquid_sum) > arithmetic_scale ||
            std::abs(vapor_sum) > arithmetic_scale) {
            detail::pt_sensitivity_set_failure(
                sensitivity, PtSensitivityStatus::arithmetic_failure,
                "PT sensitivity: differentiated phase normalization failed");
            return result;
        }
        for (std::size_t i = 0; i < n; ++i) {
            const double dz =
                column < 2
                    ? 0.0
                    : (i == column - 2 ? 1.0 : (i + 1 == n ? -1.0 : 0.0));
            const double x = local.values[n + 1 + i];
            const double y = local.values[2 * n + 1 + i];
            const double derivative =
                (1.0 - beta_base) * sensitivity.d_liquid(i, column) +
                beta_base * sensitivity.d_vapor(i, column) +
                (y - x) * sensitivity.d_vapor_fraction(column) - dz;
            const double scale = std::max(
                1.0,
                std::abs(sensitivity.d_liquid(i, column)) +
                    std::abs(sensitivity.d_vapor(i, column)) +
                    std::abs(sensitivity.d_vapor_fraction(column)) +
                    std::abs(dz));
            if (std::abs(derivative) > arithmetic_scale * scale) {
                detail::pt_sensitivity_set_failure(
                    sensitivity, PtSensitivityStatus::arithmetic_failure,
                    "PT sensitivity: differentiated material balance failed");
                return result;
            }
        }
    }

    sensitivity.status = PtSensitivityStatus::success;
    sensitivity.diagnostic = "PT sensitivity: local implicit derivative available";
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_PR76_SENSITIVITY_HPP
