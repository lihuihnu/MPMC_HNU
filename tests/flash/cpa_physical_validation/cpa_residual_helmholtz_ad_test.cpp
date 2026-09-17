#include <mpmc/ad/math.hpp>
#include <mpmc/ad/runtime_differentiate.hpp>
#include <mpmc/thermodynamics/cpa_pt_phase.hpp>
#include <mpmc/thermodynamics/cpa_residual_helmholtz.hpp>

#include "cpa_thermopack_parameter_snapshot.hpp"
#include "cpa_thermopack_parity_thresholds.hpp"
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
#include <string_view>
#include <vector>

namespace {
namespace ad = mpmc::ad;
namespace th = mpmc::thermodynamics;
namespace external = cpa_thermopack_phase_kernel;
namespace thresholds = cpa_thermopack_parity_thresholds;

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

struct Derived {
    double cubic_value{};
    double association_q_value{};
    double total_value{};
    double pressure_physical_pa{};
    double pressure_association_pa{};
    double pressure_total_pa{};
    std::vector<double> mu_cubic;
    std::vector<double> mu_association;
    std::vector<double> mu_total;
    std::vector<double> ln_phi;
};

Derived derive(
    double target_pressure_pa,
    double temperature_k,
    double rho,
    std::span<const double> composition,
    const th::CpaParameterSet& parameters,
    const th::CpaAssociationResult& association) {
    require(composition.size() == parameters.size(),
            "Helmholtz AD adapter composition dimension mismatch");
    require(target_pressure_pa > 0.0 && std::isfinite(target_pressure_pa),
            "Helmholtz AD adapter requires positive target pressure");

    std::vector<double> inputs;
    inputs.reserve(composition.size() + 1U);
    inputs.push_back(1.0 / rho);
    inputs.insert(inputs.end(), composition.begin(), composition.end());

    using Number = ad::Dual<double, 4>;
    ad::RuntimeJacobianWorkspace<double, 4> workspace;
    const auto callback = [&](std::span<const Number> variables,
                              std::span<Number> outputs) {
        const Number temperature{temperature_k};
        const Number& volume = variables.front();
        const auto mole_numbers = variables.subspan(1U);
        outputs[0] = th::cpa_cubic_residual_helmholtz_reduced(
            temperature, volume, mole_numbers, parameters);
        outputs[1] = th::cpa_association_q_reduced(
            temperature, volume, mole_numbers, parameters, association);
        outputs[2] = th::cpa_residual_helmholtz_reduced(
            temperature, volume, mole_numbers, parameters, association);
    };

    const auto result = ad::value_and_jacobian_runtime<4>(
        callback, std::span<const double>{inputs}, 3U, workspace,
        {1024U, 8U, 8192U});
    require(result.input_count == inputs.size() && result.output_count == 3U,
            "Helmholtz AD adapter returned unexpected Jacobian shape");

    const std::size_t input_count = result.input_count;
    const auto d = [&](std::size_t output, std::size_t input) {
        return result.jacobian[output * input_count + input];
    };
    const double total_moles = std::accumulate(
        composition.begin(), composition.end(), 0.0);
    const double rt = th::cpa_gas_constant_j_per_mol_k * temperature_k;
    const double ideal_pressure = total_moles * rt / inputs.front();
    const double z = target_pressure_pa / (rho * rt);
    require(z > 0.0 && std::isfinite(z),
            "Helmholtz AD adapter produced invalid target Z");

    Derived derived;
    derived.cubic_value = result.values[0];
    derived.association_q_value = result.values[1];
    derived.total_value = result.values[2];
    derived.pressure_physical_pa = ideal_pressure - rt * d(0U, 0U);
    derived.pressure_association_pa = -rt * d(1U, 0U);
    derived.pressure_total_pa = ideal_pressure - rt * d(2U, 0U);
    derived.mu_cubic.resize(parameters.size());
    derived.mu_association.resize(parameters.size());
    derived.mu_total.resize(parameters.size());
    derived.ln_phi.resize(parameters.size());
    const double log_z = std::log(z);
    for (std::size_t i = 0U; i < parameters.size(); ++i) {
        derived.mu_cubic[i] = d(0U, i + 1U);
        derived.mu_association[i] = d(1U, i + 1U);
        derived.mu_total[i] = d(2U, i + 1U);
        derived.ln_phi[i] = derived.mu_total[i] - log_z;
    }
    return derived;
}

double direct_assoc_value(
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

std::vector<double> analytic_assoc_mu(
    double rho,
    std::span<const double> composition,
    const th::CpaParameterSet& parameters,
    const th::CpaAssociationResult& association) {
    std::vector<double> log_x(parameters.size(), 0.0);
    double association_sum = 0.0;
    for (const auto& site : association.sites) {
        const double multiplicity = static_cast<double>(site.multiplicity);
        log_x[site.component_index] += multiplicity * std::log(site.unbonded_fraction);
        association_sum += composition[site.component_index] * multiplicity *
            (1.0 - site.unbonded_fraction);
    }
    std::vector<double> result(parameters.size(), 0.0);
    for (std::size_t i = 0U; i < parameters.size(); ++i) {
        result[i] = log_x[i] - (1.9 / 8.0) * rho *
            parameters.pure(i).b_m3_per_mol * association.radial_distribution *
            association_sum;
    }
    return result;
}

void check_internal_state(
    double temperature_k,
    double rho,
    std::span<const double> composition,
    const th::CpaParameterSet& parameters,
    std::string_view label,
    double& max_q_delta,
    double& max_pressure_delta,
    double& max_mu_delta,
    double& max_lnphi_delta) {
    const auto state = th::evaluate_cpa_phase_at_density(
        temperature_k, rho, composition, parameters);
    require(state.pressure_pa > 0.0,
            std::string(label) + ": internal state pressure must be positive");

    th::CpaPtRoot analytic_root;
    analytic_root.molar_density_mol_per_m3 = rho;
    th::cpa_detail::cpa_fill_ln_phi(
        state.pressure_pa, temperature_k, composition,
        parameters, state, analytic_root);

    const auto derived = derive(
        state.pressure_pa, temperature_k, rho,
        composition, parameters, state.association);
    const double direct_assoc = direct_assoc_value(composition, state.association);
    const auto assoc_mu = analytic_assoc_mu(
        rho, composition, parameters, state.association);
    const double rt = th::cpa_gas_constant_j_per_mol_k * temperature_k;
    const double z = state.pressure_pa / (rho * rt);

    const double q_delta = std::abs(derived.association_q_value - direct_assoc);
    max_q_delta = std::max(max_q_delta, q_delta);
    require_abs(
        derived.association_q_value, direct_assoc,
        1.0e-10 + 1.0e-12 * std::max(
            std::abs(derived.association_q_value), std::abs(direct_assoc)),
        std::string(label) + ": stationary Q value mismatch");
    require_abs(
        derived.total_value,
        derived.cubic_value + derived.association_q_value,
        roundoff(derived.cubic_value) + roundoff(derived.association_q_value),
        std::string(label) + ": total Helmholtz value mismatch");

    const double p_physical_delta = std::abs(
        derived.pressure_physical_pa - state.pressure_physical_pa);
    const double p_assoc_delta = std::abs(
        derived.pressure_association_pa - state.pressure_association_pa);
    const double p_total_delta = std::abs(
        derived.pressure_total_pa - state.pressure_pa);
    max_pressure_delta = std::max(
        {max_pressure_delta, p_physical_delta, p_assoc_delta, p_total_delta});
    require_abs(
        derived.pressure_physical_pa, state.pressure_physical_pa,
        roundoff(state.pressure_physical_pa),
        std::string(label) + ": AD physical pressure mismatch");
    require_abs(
        derived.pressure_association_pa, state.pressure_association_pa,
        roundoff(state.pressure_association_pa),
        std::string(label) + ": AD association pressure mismatch");
    const double pressure_sum_tolerance =
        roundoff(state.pressure_physical_pa) +
        roundoff(state.pressure_association_pa) +
        roundoff(state.pressure_pa);
    require_abs(
        derived.pressure_total_pa, state.pressure_pa,
        pressure_sum_tolerance,
        std::string(label) + ": AD total pressure mismatch after componentwise roundoff propagation");
    require_abs(
        derived.pressure_total_pa,
        derived.pressure_physical_pa + derived.pressure_association_pa,
        pressure_sum_tolerance,
        std::string(label) + ": total derivative disagrees with derived pressure decomposition");

    for (std::size_t i = 0U; i < parameters.size(); ++i) {
        const double total_mu = analytic_root.ln_phi[i] + std::log(z);
        const double cubic_mu = total_mu - assoc_mu[i];
        const double cubic_delta = std::abs(derived.mu_cubic[i] - cubic_mu);
        const double assoc_delta = std::abs(derived.mu_association[i] - assoc_mu[i]);
        const double total_delta = std::abs(derived.mu_total[i] - total_mu);
        const double lnphi_delta = std::abs(derived.ln_phi[i] - analytic_root.ln_phi[i]);
        max_mu_delta = std::max(
            {max_mu_delta, cubic_delta, assoc_delta, total_delta});
        max_lnphi_delta = std::max(max_lnphi_delta, lnphi_delta);

        require_abs(
            derived.mu_cubic[i], cubic_mu, roundoff(cubic_mu),
            std::string(label) + ": AD cubic residual chemical potential mismatch");
        require_abs(
            derived.mu_association[i], assoc_mu[i], roundoff(assoc_mu[i]),
            std::string(label) + ": AD association residual chemical potential mismatch");
        const double mu_sum_tolerance =
            roundoff(cubic_mu) + roundoff(assoc_mu[i]) + roundoff(total_mu);
        require_abs(
            derived.mu_total[i], total_mu, mu_sum_tolerance,
            std::string(label) + ": AD total residual chemical potential mismatch");
        require_abs(
            derived.mu_total[i], derived.mu_cubic[i] + derived.mu_association[i],
            mu_sum_tolerance,
            std::string(label) + ": total derivative disagrees with derived chemical-potential decomposition");
        require_abs(
            derived.ln_phi[i], analytic_root.ln_phi[i],
            mu_sum_tolerance + roundoff(std::log(z)),
            std::string(label) + ": AD ln(phi) mismatch");
    }
}

void run_internal_equivalence() {
    double max_q_delta = 0.0;
    double max_pressure_delta = 0.0;
    double max_mu_delta = 0.0;
    double max_lnphi_delta = 0.0;

    for (const bool swapped : {false, true}) {
        const auto literature = cpa_physical_test::parameters(swapped);
        const auto parity = cpa_thermopack_snapshot::parameters(swapped);
        for (std::size_t state_index = 0U;
             state_index < external::states.size(); ++state_index) {
            const auto& state = external::states[state_index];
            for (const auto* phase : {&state.liquid, &state.vapor}) {
                const auto composition = cpa_physical_test::composition(
                    phase->composition_methanol, swapped);
                const std::string suffix =
                    " state=" + std::to_string(state_index) +
                    (phase == &state.liquid ? " phase=liquid" : " phase=vapor") +
                    (swapped ? " order=swapped" : " order=normal");
                check_internal_state(
                    state.temperature_k, phase->molar_density_mol_per_m3,
                    composition, literature, "literature" + suffix,
                    max_q_delta, max_pressure_delta, max_mu_delta, max_lnphi_delta);
                check_internal_state(
                    state.temperature_k, phase->molar_density_mol_per_m3,
                    composition, parity, "parity" + suffix,
                    max_q_delta, max_pressure_delta, max_mu_delta, max_lnphi_delta);
            }
        }
    }

    std::cout << std::setprecision(17)
              << "CPA_HELMHOLTZ_INTERNAL_OK"
              << " max_abs_dQ_value=" << max_q_delta
              << " max_abs_dP_pa=" << max_pressure_delta
              << " max_abs_dmu_over_rt=" << max_mu_delta
              << " max_abs_dlnphi=" << max_lnphi_delta
              << " snapshots=literature,thermopack-parity"
              << " component_orders=normal,swapped"
              << '\n';
}

void check_external_phase(
    const external::ThermoPackPhaseKernelState& state,
    const external::PhaseReference& phase,
    const th::CpaParameterSet& parameters,
    double& max_p_physical,
    double& max_p_assoc,
    double& max_mu_cubic,
    double& max_mu_assoc,
    double& max_lnphi) {
    const auto composition = cpa_physical_test::composition(
        phase.composition_methanol, false);
    const auto primal = th::evaluate_cpa_phase_at_density(
        state.temperature_k, phase.molar_density_mol_per_m3,
        composition, parameters);
    const auto derived = derive(
        state.pressure_pa, state.temperature_k,
        phase.molar_density_mol_per_m3,
        composition, parameters, primal.association);

    const double d_p_physical = std::abs(
        derived.pressure_physical_pa - phase.pressure_physical_pa);
    const double d_p_assoc = std::abs(
        derived.pressure_association_pa - phase.pressure_association_pa);
    max_p_physical = std::max(max_p_physical, d_p_physical);
    max_p_assoc = std::max(max_p_assoc, d_p_assoc);
    require(d_p_physical <= thresholds::phase_max_abs_pressure_physical_pa,
            "Helmholtz AD physical pressure exceeded frozen ThermoPack parity gate");
    require(d_p_assoc <= thresholds::phase_max_abs_pressure_association_pa,
            "Helmholtz AD association pressure exceeded frozen ThermoPack parity gate");

    const std::array<double, 2> ref_mu_cubic{{
        phase.mu_cubic_over_rt_methanol,
        phase.mu_cubic_over_rt_water}};
    const std::array<double, 2> ref_mu_assoc{{
        phase.mu_association_over_rt_methanol,
        phase.mu_association_over_rt_water}};
    const std::array<double, 2> ref_lnphi{{
        phase.ln_phi_methanol,
        phase.ln_phi_water}};
    for (std::size_t i = 0U; i < 2U; ++i) {
        const double d_mu_cubic = std::abs(derived.mu_cubic[i] - ref_mu_cubic[i]);
        const double d_mu_assoc = std::abs(derived.mu_association[i] - ref_mu_assoc[i]);
        const double d_lnphi = std::abs(derived.ln_phi[i] - ref_lnphi[i]);
        max_mu_cubic = std::max(max_mu_cubic, d_mu_cubic);
        max_mu_assoc = std::max(max_mu_assoc, d_mu_assoc);
        max_lnphi = std::max(max_lnphi, d_lnphi);
        require(d_mu_cubic <= thresholds::phase_max_abs_mu_cubic_over_rt,
                "Helmholtz AD cubic mu exceeded frozen ThermoPack parity gate");
        require(d_mu_assoc <= thresholds::phase_max_abs_mu_association_over_rt,
                "Helmholtz AD association mu exceeded frozen ThermoPack parity gate");
        require(d_lnphi <= thresholds::phase_max_abs_ln_phi,
                "Helmholtz AD ln(phi) exceeded frozen ThermoPack parity gate");
    }
}

void run_external_parity() {
    const auto parameters = cpa_thermopack_snapshot::parameters(false);
    double max_p_physical = 0.0;
    double max_p_assoc = 0.0;
    double max_mu_cubic = 0.0;
    double max_mu_assoc = 0.0;
    double max_lnphi = 0.0;

    for (const auto& state : external::states) {
        check_external_phase(
            state, state.liquid, parameters,
            max_p_physical, max_p_assoc,
            max_mu_cubic, max_mu_assoc, max_lnphi);
        check_external_phase(
            state, state.vapor, parameters,
            max_p_physical, max_p_assoc,
            max_mu_cubic, max_mu_assoc, max_lnphi);
    }

    std::cout << std::setprecision(17)
              << "CPA_HELMHOLTZ_THERMOPACK_OK"
              << " threshold_contract=" << thresholds::contract
              << " max_abs_dPphysical_pa=" << max_p_physical
              << " max_abs_dPassociation_pa=" << max_p_assoc
              << " max_abs_dmu_cubic_over_rt=" << max_mu_cubic
              << " max_abs_dmu_association_over_rt=" << max_mu_assoc
              << " max_abs_dlnphi=" << max_lnphi
              << " states=10"
              << '\n';
}

} // namespace

int main() {
    try {
        run_internal_equivalence();
        run_external_parity();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
