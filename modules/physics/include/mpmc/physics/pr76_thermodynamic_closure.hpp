#ifndef MPMC_PHYSICS_PR76_THERMODYNAMIC_CLOSURE_HPP
#define MPMC_PHYSICS_PR76_THERMODYNAMIC_CLOSURE_HPP

#include <mpmc/ad/math.hpp>
#include <mpmc/ad/runtime_differentiate.hpp>
#include <mpmc/flash/pr76_sensitivity.hpp>
#include <mpmc/physics/thermodynamic_closure.hpp>
#include <mpmc/thermodynamics/pr76_pure.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace mpmc::physics {

struct Pr76PtVleClosureOptions {
    mpmc::flash::PtSensitivityOptions sensitivity;
};

namespace detail {

struct Pr76PhaseZPartials {
    double value{};
    // Derivatives with respect to (p, T, w_0, ..., w_{N-2}).
    std::vector<double> gradient;
};

inline void pr76_closure_validate_snapshot(
    const mpmc::flash::Pr76PtSplitResult& split,
    const mpmc::flash::Pr76VleEvaluator& evaluator) {
    namespace th = mpmc::thermodynamics;
    const auto& model = evaluator.model();
    const auto& parameters = model.parameters();
    const std::size_t n = model.size();
    if (n < 2 || split.component_ids.size() != n ||
        split.solution.initial_stability.feed.size() != n ||
        split.dataset_id != parameters.dataset_id() ||
        split.revision != parameters.revision() ||
        split.model_profile != th::pr76_profile ||
        split.phase_convention != th::pr76_pt_convention ||
        split.root_options.max_iterations != evaluator.root_options().max_iterations) {
        throw std::invalid_argument(
            "PR76 closure: split result and evaluator metadata do not match");
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (split.component_ids[i] != parameters.components().items()[i].id) {
            throw std::invalid_argument(
                "PR76 closure: ordered component snapshot does not match");
        }
    }
    const double p = split.solution.initial_stability.pressure_pa;
    const double t = split.solution.initial_stability.temperature_k;
    if (!std::isfinite(p) || !(p > 0.0) || !std::isfinite(t) || !(t > 0.0)) {
        throw std::domain_error("PR76 closure: invalid accepted p/T base point");
    }
}

[[nodiscard]] inline bool closure_vector_finite(
    const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

[[nodiscard]] inline Pr76PhaseZPartials pr76_phase_z_partials(
    const mpmc::thermodynamics::Pr76Phase<double>& model,
    double pressure_pa, double temperature_k,
    std::span<const double> composition, std::size_t root_index,
    mpmc::thermodynamics::Pr76RootOptions root_options) {
    constexpr std::size_t direction_width = 4;
    using Number = mpmc::ad::Dual<double, direction_width>;
    const std::size_t n = model.size();
    if (n < 2 || composition.size() != n) {
        throw std::invalid_argument("PR76 closure: invalid phase composition shape");
    }

    std::vector<double> inputs(n + 1);
    inputs[0] = pressure_pa;
    inputs[1] = temperature_k;
    for (std::size_t i = 0; i + 1 < n; ++i) {
        inputs[2 + i] = composition[i];
    }

    mpmc::ad::RuntimeJacobianWorkspace<double, direction_width> jacobian_workspace;
    mpmc::thermodynamics::Pr76PhaseWorkspace<Number> phase_workspace;
    const auto equation =
        [&](std::span<const Number> variables, std::span<Number> outputs) {
            const std::span<const Number> reduced{
                variables.data() + 2, n - 1};
            outputs[0] = model.evaluate_reduced(
                variables[0], variables[1], reduced, root_index,
                phase_workspace, root_options).z;
        };

    mpmc::ad::RuntimeJacobianLimits limits;
    limits.max_inputs = inputs.size();
    limits.max_outputs = 1;
    limits.max_jacobian_entries = inputs.size();
    const auto values = mpmc::ad::value_and_jacobian_runtime<direction_width>(
        equation, inputs, 1, jacobian_workspace, limits);
    return {values.values[0], values.jacobian};
}

[[nodiscard]] inline ThermodynamicClosureLinearizationReason map_sensitivity_failure(
    mpmc::flash::PtSensitivityStatus status) {
    using Source = mpmc::flash::PtSensitivityStatus;
    switch (status) {
    case Source::success:
        return ThermodynamicClosureLinearizationReason::none;
    case Source::phase_boundary:
        return ThermodynamicClosureLinearizationReason::phase_boundary;
    case Source::ill_conditioned_equilibrium:
        return ThermodynamicClosureLinearizationReason::ill_conditioned_equilibrium;
    case Source::unsupported_feed_support:
        return ThermodynamicClosureLinearizationReason::unsupported_feed_support;
    case Source::property_failure:
        return ThermodynamicClosureLinearizationReason::property_failure;
    case Source::arithmetic_failure:
        return ThermodynamicClosureLinearizationReason::arithmetic_failure;
    case Source::solution_not_accepted:
        return ThermodynamicClosureLinearizationReason::solution_not_accepted;
    }
    return ThermodynamicClosureLinearizationReason::arithmetic_failure;
}

} // namespace detail

// Consume an already solved PR76 PT flash result. This adapter does not solve a
// conservation law and does not own Newton globalization. It publishes a valid
// two-phase primal whenever the accepted flash state can be reproduced, and an
// optional local linearization only when all implicit-sensitivity/property checks
// succeed on the same atomic snapshot.
[[nodiscard]] inline ThermodynamicClosureSnapshot build_pr76_pt_vle_closure(
    const mpmc::flash::Pr76PtSplitResult& split,
    const mpmc::flash::Pr76VleEvaluator& evaluator,
    Pr76PtVleClosureOptions options = {}) {
    namespace fl = mpmc::flash;
    namespace th = mpmc::thermodynamics;

    detail::pr76_closure_validate_snapshot(split, evaluator);

    ThermodynamicClosureSnapshot result;
    result.pressure_pa = split.solution.initial_stability.pressure_pa;
    result.temperature_k = split.solution.initial_stability.temperature_k;
    result.feed = split.solution.initial_stability.feed;
    result.component_ids = split.component_ids;
    result.thermodynamic_model = split.model_profile;
    result.dataset_id = split.dataset_id;
    result.revision = split.revision;

    if (split.solution.status == fl::PtSplitStatus::single_phase_no_instability_found) {
        result.primal_status =
            ThermodynamicClosurePrimalStatus::phase_regime_not_implemented;
        result.diagnostic =
            "PR76 closure: single-phase closure is outside the v1 PT-VLE adapter";
        return result;
    }
    if (split.solution.status == fl::PtSplitStatus::phase_set_unstable ||
        split.solution.status == fl::PtSplitStatus::indeterminate) {
        result.primal_status = ThermodynamicClosurePrimalStatus::indeterminate;
        result.diagnostic =
            "PR76 closure: flash phase set is not an accepted thermodynamic primal";
        return result;
    }

    const auto* candidate = split.solution.candidate();
    if (candidate == nullptr || !split.solution.final_stability ||
        split.solution.final_stability->status != fl::StabilityStatus::no_instability_found) {
        result.primal_status = ThermodynamicClosurePrimalStatus::indeterminate;
        result.diagnostic =
            "PR76 closure: accepted two-phase status lacks its accepted phase set";
        return result;
    }

    const std::size_t n = evaluator.model().size();
    const double beta = candidate->fractions.vapor_fraction;
    if (!std::isfinite(beta) || !(beta > 0.0) || !(beta < 1.0) ||
        candidate->fractions.liquid.size() != n ||
        candidate->fractions.vapor.size() != n) {
        result.primal_status = ThermodynamicClosurePrimalStatus::indeterminate;
        result.diagnostic = "PR76 closure: invalid accepted two-phase primal shape";
        return result;
    }

    try {
        th::Pr76PhaseWorkspace<double> liquid_workspace;
        th::Pr76PhaseWorkspace<double> vapor_workspace;
        const std::span<const double> liquid_reduced{
            candidate->fractions.liquid.data(), n - 1};
        const std::span<const double> vapor_reduced{
            candidate->fractions.vapor.data(), n - 1};
        const auto liquid_values = evaluator.model().evaluate_reduced(
            result.pressure_pa, result.temperature_k, liquid_reduced,
            candidate->liquid.activity.branch, liquid_workspace,
            split.root_options);
        const auto vapor_values = evaluator.model().evaluate_reduced(
            result.pressure_pa, result.temperature_k, vapor_reduced,
            candidate->vapor.activity.branch, vapor_workspace,
            split.root_options);

        const double r = th::Pr76Pure<double>::gas_constant();
        const double liquid_density = result.pressure_pa /
            (liquid_values.z * r * result.temperature_k);
        const double vapor_density = result.pressure_pa /
            (vapor_values.z * r * result.temperature_k);
        if (!std::isfinite(liquid_values.z) || !(liquid_values.z > 0.0) ||
            !std::isfinite(vapor_values.z) || !(vapor_values.z > 0.0) ||
            !std::isfinite(liquid_density) || !(liquid_density > 0.0) ||
            !std::isfinite(vapor_density) || !(vapor_density > 0.0)) {
            throw std::range_error("PR76 closure: nonrepresentable phase volume state");
        }

        PtVleThermodynamicState primal;
        primal.liquid.mole_phase_fraction = 1.0 - beta;
        primal.liquid.composition = candidate->fractions.liquid;
        primal.liquid.compressibility_factor = liquid_values.z;
        primal.liquid.molar_density_mol_per_m3 = liquid_density;
        primal.vapor.mole_phase_fraction = beta;
        primal.vapor.composition = candidate->fractions.vapor;
        primal.vapor.compressibility_factor = vapor_values.z;
        primal.vapor.molar_density_mol_per_m3 = vapor_density;
        result.primal = std::move(primal);
        result.primal_status = ThermodynamicClosurePrimalStatus::valid;
    } catch (const std::exception& error) {
        result.primal_status = ThermodynamicClosurePrimalStatus::indeterminate;
        result.primal.reset();
        result.diagnostic = std::string("PR76 closure: primal reproduction failed: ") +
                            error.what();
        return result;
    }

    const auto sensitivity =
        fl::differentiate_pr76_pt_vle(split, evaluator, options.sensitivity);
    if (sensitivity.sensitivity.status != fl::PtSensitivityStatus::success) {
        result.linearization_reason = detail::map_sensitivity_failure(
            sensitivity.sensitivity.status);
        result.diagnostic = sensitivity.sensitivity.diagnostic;
        if (sensitivity.sensitivity.status ==
            fl::PtSensitivityStatus::solution_not_accepted) {
            // The independent sensitivity path could not reproduce equilibrium;
            // do not publish this inconsistent state even for a residual.
            result.primal_status = ThermodynamicClosurePrimalStatus::indeterminate;
            result.primal.reset();
        }
        return result;
    }

    const auto& s = sensitivity.sensitivity;
    const std::size_t q_count = s.input_count;
    if (q_count != n + 1 || !result.primal) {
        result.primal_status = ThermodynamicClosurePrimalStatus::indeterminate;
        result.primal.reset();
        result.linearization_reason =
            ThermodynamicClosureLinearizationReason::solution_not_accepted;
        result.diagnostic = "PR76 closure: inconsistent sensitivity shape";
        return result;
    }

    try {
        const auto liquid_z = detail::pr76_phase_z_partials(
            evaluator.model(), result.pressure_pa, result.temperature_k,
            result.primal->liquid.composition, candidate->liquid.activity.branch,
            split.root_options);
        const auto vapor_z = detail::pr76_phase_z_partials(
            evaluator.model(), result.pressure_pa, result.temperature_k,
            result.primal->vapor.composition, candidate->vapor.activity.branch,
            split.root_options);
        const double z_guard =
            1024.0 * std::numeric_limits<double>::epsilon();
        const auto same_z = [z_guard](double a, double b) {
            return std::abs(a - b) <= z_guard * std::max(1.0, std::abs(b));
        };
        if (!same_z(liquid_z.value, result.primal->liquid.compressibility_factor) ||
            !same_z(vapor_z.value, result.primal->vapor.compressibility_factor) ||
            liquid_z.gradient.size() != n + 1 ||
            vapor_z.gradient.size() != n + 1) {
            throw std::runtime_error(
                "PR76 closure: phase derivative path changed the accepted primal");
        }

        PtVleThermodynamicLinearization linearization;
        linearization.component_count = n;
        linearization.input_count = q_count;
        linearization.vapor_fraction_gradient = s.vapor_fraction_gradient;
        linearization.liquid_composition_jacobian = s.liquid_jacobian;
        linearization.vapor_composition_jacobian = s.vapor_jacobian;
        linearization.liquid_compressibility_gradient.assign(q_count, 0.0);
        linearization.vapor_compressibility_gradient.assign(q_count, 0.0);
        linearization.liquid_molar_density_gradient.assign(q_count, 0.0);
        linearization.vapor_molar_density_gradient.assign(q_count, 0.0);

        for (std::size_t column = 0; column < q_count; ++column) {
            double dz_liquid = column == 0 ? liquid_z.gradient[0] :
                               (column == 1 ? liquid_z.gradient[1] : 0.0);
            double dz_vapor = column == 0 ? vapor_z.gradient[0] :
                              (column == 1 ? vapor_z.gradient[1] : 0.0);
            for (std::size_t component = 0; component + 1 < n; ++component) {
                dz_liquid += liquid_z.gradient[2 + component] *
                    s.d_liquid(component, column);
                dz_vapor += vapor_z.gradient[2 + component] *
                    s.d_vapor(component, column);
            }
            linearization.liquid_compressibility_gradient[column] = dz_liquid;
            linearization.vapor_compressibility_gradient[column] = dz_vapor;

            const double direct_pressure =
                column == 0 ? 1.0 / result.pressure_pa : 0.0;
            const double direct_temperature =
                column == 1 ? 1.0 / result.temperature_k : 0.0;
            const auto density_derivative =
                [&](const ThermodynamicPhaseState& phase, double dz) {
                    return phase.molar_density_mol_per_m3 *
                        (direct_pressure - dz / phase.compressibility_factor -
                         direct_temperature);
                };
            linearization.liquid_molar_density_gradient[column] =
                density_derivative(result.primal->liquid, dz_liquid);
            linearization.vapor_molar_density_gradient[column] =
                density_derivative(result.primal->vapor, dz_vapor);
        }

        if (!detail::closure_vector_finite(
                linearization.vapor_fraction_gradient) ||
            !detail::closure_vector_finite(
                linearization.liquid_composition_jacobian) ||
            !detail::closure_vector_finite(
                linearization.vapor_composition_jacobian) ||
            !detail::closure_vector_finite(
                linearization.liquid_compressibility_gradient) ||
            !detail::closure_vector_finite(
                linearization.vapor_compressibility_gradient) ||
            !detail::closure_vector_finite(
                linearization.liquid_molar_density_gradient) ||
            !detail::closure_vector_finite(
                linearization.vapor_molar_density_gradient)) {
            throw std::runtime_error("PR76 closure: nonfinite total derivative");
        }

        result.linearization = std::move(linearization);
        result.linearization_status =
            ThermodynamicClosureLinearizationStatus::available;
        result.linearization_reason =
            ThermodynamicClosureLinearizationReason::none;
        result.diagnostic =
            "PR76 closure: two-phase primal and local linearization available";
        return result;
    } catch (const th::Pr76PhaseError& error) {
        result.linearization_reason =
            ThermodynamicClosureLinearizationReason::property_failure;
        result.diagnostic = error.what();
    } catch (const std::range_error& error) {
        result.linearization_reason =
            ThermodynamicClosureLinearizationReason::property_failure;
        result.diagnostic = error.what();
    } catch (const std::exception& error) {
        result.linearization_reason =
            ThermodynamicClosureLinearizationReason::arithmetic_failure;
        result.diagnostic = error.what();
    }
    result.linearization.reset();
    result.linearization_status =
        ThermodynamicClosureLinearizationStatus::unavailable;
    return result;
}

} // namespace mpmc::physics

#endif // MPMC_PHYSICS_PR76_THERMODYNAMIC_CLOSURE_HPP
