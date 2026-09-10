#ifndef MPMC_FLASH_DETAIL_PR76_SENSITIVITY_LOCAL_HPP
#define MPMC_FLASH_DETAIL_PR76_SENSITIVITY_LOCAL_HPP

#include <mpmc/ad/math.hpp>
#include <mpmc/ad/runtime_differentiate.hpp>
#include <mpmc/flash/detail/pt_sensitivity_detail.hpp>
#include <mpmc/flash/pr76_split.hpp>

#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace mpmc::flash::detail {

struct Pr76SensitivityLocalSystem {
    std::size_t variable_count{};
    std::size_t output_count{};
    ad::RuntimeValueAndJacobian<double> values;
};

// Variables are [logK_0..logK_N-1, p, T, z_0..z_N-2]. Outputs are
// [F_0..F_N-1, beta, x_0..x_N-1, y_0..y_N-1]. The callback is a local
// differential model at one accepted state: its RR primal comes from that state
// and only the algebraic IFT slopes depend on logK/feed seeds. The PR76 phase
// kernel performs its own existing simple-root IFT, so no iterative algorithm is
// differentiated anywhere in this path.
[[nodiscard]] inline Pr76SensitivityLocalSystem pr76_sensitivity_local_system(
    const Pr76PtSplitResult& split, const Pr76VleEvaluator& evaluator,
    const PtSplitState& candidate) {
    constexpr std::size_t direction_width = 4;
    using Number = ad::Dual<double, direction_width>;

    const auto& model = evaluator.model();
    const auto& feed = split.solution.initial_stability.feed;
    const std::size_t n = model.size();
    const std::size_t variable_count = 2 * n + 1;
    const std::size_t output_count = 3 * n + 1;

    std::vector<double> inputs(variable_count);
    for (std::size_t i = 0; i < n; ++i) {
        inputs[i] = candidate.log_k[i];
    }
    inputs[n] = split.solution.initial_stability.pressure_pa;
    inputs[n + 1] = split.solution.initial_stability.temperature_k;
    for (std::size_t i = 0; i + 1 < n; ++i) {
        inputs[n + 2 + i] = feed[i];
    }

    ad::RuntimeJacobianWorkspace<double, direction_width> jacobian_workspace;
    thermodynamics::Pr76PhaseWorkspace<Number> liquid_workspace;
    thermodynamics::Pr76PhaseWorkspace<Number> vapor_workspace;
    const auto equations =
        [&](std::span<const Number> variables, std::span<Number> outputs) {
            std::vector<Number> log_k(n);
            std::vector<Number> z(n);
            for (std::size_t i = 0; i < n; ++i) {
                log_k[i] = variables[i];
            }
            Number feed_sum{0.0};
            for (std::size_t i = 0; i + 1 < n; ++i) {
                z[i] = variables[n + 2 + i];
                feed_sum += z[i];
            }
            z[n - 1] = Number{1.0} - feed_sum;

            const auto rr = pt_sensitivity_rr_linearization(
                std::span<const Number>{z.data(), z.size()},
                std::span<const Number>{log_k.data(), log_k.size()}, candidate);
            const auto liquid = model.evaluate_reduced(
                variables[n], variables[n + 1],
                std::span<const Number>{rr.liquid.data(), n - 1},
                candidate.liquid.activity.branch, liquid_workspace,
                split.root_options);
            const auto vapor = model.evaluate_reduced(
                variables[n], variables[n + 1],
                std::span<const Number>{rr.vapor.data(), n - 1},
                candidate.vapor.activity.branch, vapor_workspace,
                split.root_options);

            using std::log;
            for (std::size_t i = 0; i < n; ++i) {
                outputs[i] =
                    log(rr.liquid[i]) - log(rr.vapor[i]) +
                    liquid.ln_phi[i] - vapor.ln_phi[i];
            }
            outputs[n] = rr.vapor_fraction;
            for (std::size_t i = 0; i < n; ++i) {
                outputs[n + 1 + i] = rr.liquid[i];
                outputs[2 * n + 1 + i] = rr.vapor[i];
            }
        };

    ad::RuntimeJacobianLimits limits;
    limits.max_inputs = variable_count;
    limits.max_outputs = output_count;
    if (variable_count >
        std::numeric_limits<std::size_t>::max() / output_count) {
        throw std::length_error("PT sensitivity: Jacobian shape overflow");
    }
    limits.max_jacobian_entries = variable_count * output_count;
    return {
        variable_count,
        output_count,
        ad::value_and_jacobian_runtime<direction_width>(
            equations, inputs, output_count, jacobian_workspace, limits)
    };
}

} // namespace mpmc::flash::detail

#endif // MPMC_FLASH_DETAIL_PR76_SENSITIVITY_LOCAL_HPP
