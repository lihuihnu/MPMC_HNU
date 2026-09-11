#ifndef MPMC_FLASH_DETAIL_SW92_PROFILE_C_SENSITIVITY_LOCAL_HPP
#define MPMC_FLASH_DETAIL_SW92_PROFILE_C_SENSITIVITY_LOCAL_HPP

#include <mpmc/ad/math.hpp>
#include <mpmc/ad/runtime_differentiate.hpp>
#include <mpmc/flash/sw92_profile_c_phase_set.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace mpmc::flash::detail {

struct Sw92ProfileCSensitivityLocalSystem {
    std::size_t component_count{};
    std::size_t phase_count{};
    std::size_t unknown_count{};
    std::size_t input_count{};
    std::size_t variable_count{};
    std::size_t residual_count{};
    std::size_t phase_fraction_output{};
    std::size_t composition_output{};
    std::size_t compressibility_output{};
    std::size_t density_output{};
    std::size_t output_count{};
    // Internal chart only. The reference component/phase is the largest accepted
    // entry, so all base log ratios are <= 0 and softmax reconstruction is
    // well-scaled. External q remains (p,T,z_0,...,z_{N-2}).
    std::vector<std::size_t> phase_reference_components;
    std::size_t phase_fraction_reference{};
    ad::RuntimeValueAndJacobian<double> values;
};

// Local fixed-phase-set coordinates:
//
// For each accepted phase alpha choose reference component r_alpha at the
// largest accepted x. Internal composition unknowns are
//     eta_{alpha,i} = log(x_{alpha,i}/x_{alpha,r}), i != r.
// Choose the largest accepted phase fraction as beta reference and use
//     tau_alpha = log(beta_alpha/beta_r), alpha != r.
// Softmax reconstruction enforces positivity and exact simplex closure while
// avoiding both trace-component subtraction and the 1/x scaling of raw mole-
// fraction coordinates. These are internal IFT coordinates only.
//
// q = [p, T, z_0..z_{N-2}], with the last overall-feed component dependent.
// The residual has exactly PN-1 equations for PN-1 unknowns:
//   - N common-chemical-potential equations for each phase a>0 versus phase 0;
//   - N-1 component material balances.
//
// SW92 phase properties are differentiated by forward AD through the existing
// corrected-original pure/mixing/fugacity kernel. The selected PR cubic root is
// differentiated by its simple-root IFT inside Sw92Phase; root iteration is not
// differentiated.
[[nodiscard]] inline Sw92ProfileCSensitivityLocalSystem
sw92_profile_c_sensitivity_local_system(
    const Sw92ProfileCPtPhaseSetResult& source,
    const thermodynamics::Sw92Phase<double>& model,
    thermodynamics::Sw92RootOptions root_options = {}) {
    constexpr std::size_t direction_width = 4;
    using Number = ad::Dual<double, direction_width>;

    const auto* accepted = source.solution.accepted_phase_set();
    if (accepted == nullptr || accepted->phases.empty() ||
        accepted->phases.size() != source.phase_metadata.size()) {
        throw std::invalid_argument(
            "SW92 sensitivity local system: accepted phase-set shape mismatch");
    }
    const std::size_t n = source.component_ids.size();
    const std::size_t pcount = accepted->phases.size();
    if (n < 2 || source.solution.feed.size() != n || pcount > 3U) {
        throw std::invalid_argument(
            "SW92 sensitivity local system: unsupported component/phase shape");
    }

    const std::size_t unknown_count = pcount * (n - 1U) + (pcount - 1U);
    const std::size_t input_count = n + 1U;
    if (unknown_count > std::numeric_limits<std::size_t>::max() - input_count) {
        throw std::length_error("SW92 sensitivity local system: variable shape overflow");
    }
    const std::size_t variable_count = unknown_count + input_count;
    const std::size_t residual_count = unknown_count;
    const std::size_t phase_fraction_output = residual_count;
    const std::size_t composition_output = phase_fraction_output + pcount;
    const std::size_t compressibility_output = composition_output + pcount * n;
    const std::size_t density_output = compressibility_output + pcount;
    const std::size_t output_count = density_output + pcount;

    std::vector<std::size_t> reference_component(pcount, 0U);
    std::vector<double> inputs(variable_count, 0.0);
    for (std::size_t phase = 0; phase < pcount; ++phase) {
        const auto& composition = accepted->phases[phase].composition;
        if (composition.size() != n) {
            throw std::invalid_argument(
                "SW92 sensitivity local system: phase composition dimension mismatch");
        }
        reference_component[phase] = static_cast<std::size_t>(
            std::distance(
                composition.begin(),
                std::max_element(composition.begin(), composition.end())));
        const double reference = composition[reference_component[phase]];
        if (!std::isfinite(reference) || !(reference > 0.0)) {
            throw std::domain_error(
                "SW92 sensitivity local system: invalid composition reference");
        }
        std::size_t local_column = 0U;
        for (std::size_t i = 0; i < n; ++i) {
            if (i == reference_component[phase]) { continue; }
            if (!std::isfinite(composition[i]) || !(composition[i] > 0.0)) {
                throw std::domain_error(
                    "SW92 sensitivity local system: log-ratio chart requires positive composition");
            }
            inputs[phase * (n - 1U) + local_column] =
                std::log(composition[i] / reference);
            ++local_column;
        }
    }

    const std::size_t beta_offset = pcount * (n - 1U);
    std::size_t beta_reference = 0U;
    for (std::size_t phase = 1U; phase < pcount; ++phase) {
        if (accepted->phases[phase].mole_phase_fraction >
            accepted->phases[beta_reference].mole_phase_fraction) {
            beta_reference = phase;
        }
    }
    const double beta_reference_value =
        accepted->phases[beta_reference].mole_phase_fraction;
    if (!std::isfinite(beta_reference_value) || !(beta_reference_value > 0.0)) {
        throw std::domain_error(
            "SW92 sensitivity local system: invalid phase-fraction reference");
    }
    std::size_t beta_column = 0U;
    for (std::size_t phase = 0; phase < pcount; ++phase) {
        if (phase == beta_reference) { continue; }
        const double value = accepted->phases[phase].mole_phase_fraction;
        if (!std::isfinite(value) || !(value > 0.0)) {
            throw std::domain_error(
                "SW92 sensitivity local system: log-ratio chart requires positive phase fractions");
        }
        inputs[beta_offset + beta_column] =
            std::log(value / beta_reference_value);
        ++beta_column;
    }

    const std::size_t q_offset = unknown_count;
    inputs[q_offset] = source.solution.pressure_pa;
    inputs[q_offset + 1U] = source.solution.temperature_k;
    for (std::size_t i = 0; i + 1U < n; ++i) {
        inputs[q_offset + 2U + i] = source.solution.feed[i];
    }

    ad::RuntimeJacobianWorkspace<double, direction_width> jacobian_workspace;
    thermodynamics::Sw92PhaseWorkspace<Number> phase_workspace;
    const auto equations =
        [&](std::span<const Number> variables, std::span<Number> outputs) {
            const Number& pressure = variables[q_offset];
            const Number& temperature = variables[q_offset + 1U];

            std::vector<Number> feed(n);
            Number feed_sum{0.0};
            for (std::size_t i = 0; i + 1U < n; ++i) {
                feed[i] = variables[q_offset + 2U + i];
                feed_sum += feed[i];
            }
            feed[n - 1U] = Number{1.0} - feed_sum;

            using std::exp;
            std::vector<std::vector<Number>> composition(
                pcount, std::vector<Number>(n));
            for (std::size_t phase = 0; phase < pcount; ++phase) {
                Number denominator{1.0};
                std::size_t local_column = 0U;
                for (std::size_t i = 0; i < n; ++i) {
                    if (i == reference_component[phase]) { continue; }
                    composition[phase][i] =
                        exp(variables[phase * (n - 1U) + local_column]);
                    denominator += composition[phase][i];
                    ++local_column;
                }
                composition[phase][reference_component[phase]] = Number{1.0};
                for (std::size_t i = 0; i < n; ++i) {
                    composition[phase][i] /= denominator;
                }
            }

            std::vector<Number> beta(pcount);
            Number beta_denominator{1.0};
            beta_column = 0U;
            for (std::size_t phase = 0; phase < pcount; ++phase) {
                if (phase == beta_reference) { continue; }
                beta[phase] = exp(variables[beta_offset + beta_column]);
                beta_denominator += beta[phase];
                ++beta_column;
            }
            beta[beta_reference] = Number{1.0};
            for (auto& value : beta) { value /= beta_denominator; }

            std::vector<thermodynamics::Sw92PhaseValues<Number, double>> values;
            values.reserve(pcount);
            for (std::size_t phase = 0; phase < pcount; ++phase) {
                values.push_back(model.evaluate_full(
                    pressure, temperature,
                    std::span<const Number>{
                        composition[phase].data(), composition[phase].size()},
                    source.nacl_molality_mol_per_kg_water,
                    source.phase_metadata[phase].thermodynamic_family,
                    accepted->phases[phase].activity.branch,
                    phase_workspace, root_options));
            }

            using std::log;
            std::size_t row = 0;
            for (std::size_t phase = 1U; phase < pcount; ++phase) {
                for (std::size_t i = 0; i < n; ++i) {
                    outputs[row++] =
                        log(composition[0][i]) + values[0].ln_phi[i] -
                        log(composition[phase][i]) - values[phase].ln_phi[i];
                }
            }
            for (std::size_t i = 0; i + 1U < n; ++i) {
                Number recovered{0.0};
                for (std::size_t phase = 0; phase < pcount; ++phase) {
                    recovered += beta[phase] * composition[phase][i];
                }
                outputs[row++] = recovered - feed[i];
            }
            if (row != residual_count) {
                throw std::logic_error(
                    "SW92 sensitivity local system: residual shape mismatch");
            }

            for (std::size_t phase = 0; phase < pcount; ++phase) {
                outputs[phase_fraction_output + phase] = beta[phase];
                for (std::size_t i = 0; i < n; ++i) {
                    outputs[composition_output + phase * n + i] =
                        composition[phase][i];
                }
                outputs[compressibility_output + phase] = values[phase].z;
                outputs[density_output + phase] =
                    pressure /
                    (values[phase].z *
                     thermodynamics::Sw92Pure<double>::gas_constant() *
                     temperature);
            }
        };

    ad::RuntimeJacobianLimits limits;
    limits.max_inputs = variable_count;
    limits.max_outputs = output_count;
    if (variable_count >
        std::numeric_limits<std::size_t>::max() / output_count) {
        throw std::length_error("SW92 sensitivity local system: Jacobian overflow");
    }
    limits.max_jacobian_entries = variable_count * output_count;
    auto differentiated = ad::value_and_jacobian_runtime<direction_width>(
        equations, inputs, output_count, jacobian_workspace, limits);

    return {
        n,
        pcount,
        unknown_count,
        input_count,
        variable_count,
        residual_count,
        phase_fraction_output,
        composition_output,
        compressibility_output,
        density_output,
        output_count,
        std::move(reference_component),
        beta_reference,
        std::move(differentiated)};
}

} // namespace mpmc::flash::detail

#endif // MPMC_FLASH_DETAIL_SW92_PROFILE_C_SENSITIVITY_LOCAL_HPP
