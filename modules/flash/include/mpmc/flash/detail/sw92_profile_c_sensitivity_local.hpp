#ifndef MPMC_FLASH_DETAIL_SW92_PROFILE_C_SENSITIVITY_LOCAL_HPP
#define MPMC_FLASH_DETAIL_SW92_PROFILE_C_SENSITIVITY_LOCAL_HPP

#include <mpmc/ad/math.hpp>
#include <mpmc/ad/runtime_differentiate.hpp>
#include <mpmc/flash/sw92_profile_c_phase_set.hpp>

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
    ad::RuntimeValueAndJacobian<double> values;
};

// Local fixed-phase-set coordinates:
//
// u = [x^0_0..x^0_{N-2}, ..., x^{P-1}_0..x^{P-1}_{N-2},
//      beta_0..beta_{P-2}]
// q = [p, T, z_0..z_{N-2}]
//
// The last composition in every phase, last phase fraction and last feed
// component are dependent simplex coordinates. The residual has exactly PN-1
// equations for PN-1 unknowns:
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

    std::vector<double> inputs(variable_count, 0.0);
    for (std::size_t phase = 0; phase < pcount; ++phase) {
        const auto& composition = accepted->phases[phase].composition;
        if (composition.size() != n) {
            throw std::invalid_argument(
                "SW92 sensitivity local system: phase composition dimension mismatch");
        }
        for (std::size_t i = 0; i + 1U < n; ++i) {
            inputs[phase * (n - 1U) + i] = composition[i];
        }
    }
    const std::size_t beta_offset = pcount * (n - 1U);
    for (std::size_t phase = 0; phase + 1U < pcount; ++phase) {
        inputs[beta_offset + phase] =
            accepted->phases[phase].mole_phase_fraction;
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

            std::vector<std::vector<Number>> composition(
                pcount, std::vector<Number>(n));
            std::vector<Number> beta(pcount);
            for (std::size_t phase = 0; phase < pcount; ++phase) {
                Number sum{0.0};
                for (std::size_t i = 0; i + 1U < n; ++i) {
                    composition[phase][i] =
                        variables[phase * (n - 1U) + i];
                    sum += composition[phase][i];
                }
                composition[phase][n - 1U] = Number{1.0} - sum;
            }
            Number beta_sum{0.0};
            for (std::size_t phase = 0; phase + 1U < pcount; ++phase) {
                beta[phase] = variables[beta_offset + phase];
                beta_sum += beta[phase];
            }
            beta[pcount - 1U] = Number{1.0} - beta_sum;

            std::vector<thermodynamics::Sw92PhaseValues<Number, double>> values;
            values.reserve(pcount);
            for (std::size_t phase = 0; phase < pcount; ++phase) {
                const auto independent = std::span<const Number>{
                    composition[phase].data(), n - 1U};
                values.push_back(model.evaluate_reduced(
                    pressure, temperature, independent,
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
        std::move(differentiated)};
}

} // namespace mpmc::flash::detail

#endif // MPMC_FLASH_DETAIL_SW92_PROFILE_C_SENSITIVITY_LOCAL_HPP
