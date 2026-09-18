#include <mpmc/ad/math.hpp>
#include <mpmc/ad/runtime_differentiate.hpp>
#include <mpmc/thermodynamics/cpa_pt_phase.hpp>
#include <mpmc/thermodynamics/cpa_residual_helmholtz.hpp>

#include "cpa_clapeyron_gate_e_generated.hpp"
#include "cpa_thermopack_parameter_snapshot.hpp"
#include "test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace ad = mpmc::ad;
namespace th = mpmc::thermodynamics;
namespace external = cpa_clapeyron_gate_e;

constexpr double max_abs_f_res = 1.0e-10;
constexpr double max_abs_pressure_pa = 5.0e-6;
constexpr double max_abs_mu_res_over_rt = 1.0e-10;
constexpr double max_abs_ln_phi = 1.0e-10;

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

struct Derived {
    double f_res{};
    double pressure_pa{};
    std::array<double, 2> mu_res_over_rt{};
    std::array<double, 2> ln_phi{};
};

Derived derive(
    const external::ClapeyronState& reference,
    const th::CpaParameterSet& parameters) {
    const auto composition = cpa_physical_test::composition(
        reference.composition_methanol, false);
    const auto primal = th::evaluate_cpa_phase_at_density(
        reference.temperature_k,
        reference.molar_density_mol_per_m3,
        composition,
        parameters);
    require(
        primal.association.converged(),
        "Gate-E MPMC association solve did not converge");

    const double volume_m3 = 1.0 / reference.molar_density_mol_per_m3;
    std::vector<double> inputs;
    inputs.reserve(3U);
    inputs.push_back(volume_m3);
    inputs.insert(inputs.end(), composition.begin(), composition.end());

    using Number = ad::Dual<double, 4>;
    ad::RuntimeJacobianWorkspace<double, 4> workspace;
    const auto callback = [&](std::span<const Number> variables,
                              std::span<Number> outputs) {
        const Number temperature{reference.temperature_k};
        const Number& volume = variables.front();
        const auto mole_numbers = variables.subspan(1U);
        outputs[0] = th::cpa_residual_helmholtz_reduced(
            temperature,
            volume,
            mole_numbers,
            parameters,
            primal.association);
    };

    const auto result = ad::value_and_jacobian_runtime<4>(
        callback,
        std::span<const double>{inputs},
        1U,
        workspace,
        {1024U, 8U, 8192U});
    require(
        result.input_count == inputs.size() && result.output_count == 1U,
        "Gate-E Helmholtz AD returned unexpected Jacobian shape");

    const double total_moles = std::accumulate(
        composition.begin(), composition.end(), 0.0);
    const double rt =
        th::cpa_gas_constant_j_per_mol_k * reference.temperature_k;
    const double ideal_pressure_pa = total_moles * rt / volume_m3;

    Derived derived;
    derived.f_res = result.values[0];
    derived.pressure_pa =
        ideal_pressure_pa - rt * result.jacobian[0U];

    const double z = derived.pressure_pa * volume_m3 /
        (total_moles * rt);
    require(
        z > 0.0 && std::isfinite(z),
        "Gate-E MPMC Helmholtz derivative produced invalid Z");
    const double log_z = std::log(z);

    for (std::size_t i = 0U; i < 2U; ++i) {
        derived.mu_res_over_rt[i] = result.jacobian[i + 1U];
        derived.ln_phi[i] = derived.mu_res_over_rt[i] - log_z;
    }
    return derived;
}

void run_gate_e() {
    static_assert(external::states.size() == 10U);
    require(
        external::gas_constant_j_per_mol_k ==
            th::cpa_gas_constant_j_per_mol_k,
        "Gate-E Clapeyron/MPMC gas constant mismatch");
    const auto parameters = cpa_thermopack_snapshot::parameters(false);

    double max_f_delta = 0.0;
    double max_pressure_delta = 0.0;
    double max_mu_delta = 0.0;
    double max_lnphi_delta = 0.0;
    double max_external_assoc_residual = 0.0;
    bool gate_pass = true;

    for (std::size_t state_index = 0U;
         state_index < external::states.size();
         ++state_index) {
        const auto& reference = external::states[state_index];
        const auto derived = derive(reference, parameters);

        const double f_delta =
            std::abs(derived.f_res - reference.f_res);
        const double pressure_delta =
            std::abs(derived.pressure_pa - reference.pressure_pa);
        max_f_delta = std::max(max_f_delta, f_delta);
        max_pressure_delta = std::max(
            max_pressure_delta, pressure_delta);
        max_external_assoc_residual = std::max(
            max_external_assoc_residual,
            reference.association_equation_max_residual);

        gate_pass = gate_pass &&
            f_delta <= max_abs_f_res &&
            pressure_delta <= max_abs_pressure_pa;

        const std::array<double, 2> ref_mu{{
            reference.mu_res_over_rt_methanol,
            reference.mu_res_over_rt_water}};
        const std::array<double, 2> ref_lnphi{{
            reference.ln_phi_methanol,
            reference.ln_phi_water}};
        for (std::size_t component = 0U; component < 2U; ++component) {
            const double mu_delta = std::abs(
                derived.mu_res_over_rt[component] - ref_mu[component]);
            const double lnphi_delta = std::abs(
                derived.ln_phi[component] - ref_lnphi[component]);
            max_mu_delta = std::max(max_mu_delta, mu_delta);
            max_lnphi_delta = std::max(max_lnphi_delta, lnphi_delta);

            gate_pass = gate_pass &&
                mu_delta <= max_abs_mu_res_over_rt &&
                lnphi_delta <= max_abs_ln_phi;
        }
    }

    std::cout << std::setprecision(17)
              << "CPA_HELMHOLTZ_CLAPEYRON_GATE_E_AUDIT"
              << " source_commit=" << external::source_commit
              << " states=" << external::states.size()
              << " max_abs_dF_res=" << max_f_delta
              << " max_abs_dP_pa=" << max_pressure_delta
              << " max_abs_dmu_over_rt=" << max_mu_delta
              << " max_abs_dlnphi=" << max_lnphi_delta
              << " max_external_assoc_equation_residual="
              << max_external_assoc_residual
              << " threshold_dF_res=" << max_abs_f_res
              << " threshold_dP_pa=" << max_abs_pressure_pa
              << " threshold_dmu_over_rt=" << max_abs_mu_res_over_rt
              << " threshold_dlnphi=" << max_abs_ln_phi
              << " result=" << (gate_pass ? "PASS" : "FAIL")
              << '\n';

    require(
        gate_pass,
        "Clapeyron Gate-E frozen acceptance contract exceeded; "
        "see CPA_HELMHOLTZ_CLAPEYRON_GATE_E_AUDIT summary");
}

} // namespace

int main() {
    try {
        run_gate_e();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
