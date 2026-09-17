#include <mpmc/ad/math.hpp>
#include <mpmc/ad/runtime_differentiate.hpp>
#include <mpmc/thermodynamics/cpa_pt_phase.hpp>
#include <mpmc/thermodynamics/cpa_residual_helmholtz.hpp>

#include "cpa_thermopack_phase_kernel_generated.hpp"
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
            "Helmholtz limit adapter composition dimension mismatch");
    require(target_pressure_pa > 0.0 && std::isfinite(target_pressure_pa),
            "Helmholtz limit adapter requires positive target pressure");

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
    const std::size_t input_count = result.input_count;
    const auto derivative = [&](std::size_t output, std::size_t input) {
        return result.jacobian[output * input_count + input];
    };

    const double total_moles = std::accumulate(
        composition.begin(), composition.end(), 0.0);
    const double rt = th::cpa_gas_constant_j_per_mol_k * temperature_k;
    const double ideal_pressure = total_moles * rt / inputs.front();
    const double z = target_pressure_pa / (rho * rt);
    require(z > 0.0 && std::isfinite(z),
            "Helmholtz limit adapter produced invalid target Z");

    Derived derived;
    derived.cubic_value = result.values[0];
    derived.association_q_value = result.values[1];
    derived.total_value = result.values[2];
    derived.pressure_physical_pa = ideal_pressure - rt * derivative(0U, 0U);
    derived.pressure_association_pa = -rt * derivative(1U, 0U);
    derived.pressure_total_pa = ideal_pressure - rt * derivative(2U, 0U);
    derived.mu_cubic.resize(parameters.size());
    derived.mu_association.resize(parameters.size());
    derived.mu_total.resize(parameters.size());
    derived.ln_phi.resize(parameters.size());
    const double log_z = std::log(z);
    for (std::size_t i = 0U; i < parameters.size(); ++i) {
        derived.mu_cubic[i] = derivative(0U, i + 1U);
        derived.mu_association[i] = derivative(1U, i + 1U);
        derived.mu_total[i] = derivative(2U, i + 1U);
        derived.ln_phi[i] = derived.mu_total[i] - log_z;
    }
    return derived;
}

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

th::CpaParameterSet pure_parameters(std::string_view component_id) {
    const auto baseline = cpa_physical_test::parameters(false);
    std::size_t selected = baseline.size();
    for (std::size_t i = 0U; i < baseline.size(); ++i) {
        if (baseline.components().at(i).id == component_id) {
            selected = i;
            break;
        }
    }
    require(selected < baseline.size(), "unknown pure-component fixture id");

    th::CpaParameterInput input;
    input.dataset_id = baseline.dataset_id() + "__pure-limit__" +
                       std::string(component_id);
    input.revision = baseline.revision() + "__pure-limit-projection";
    input.applicability = baseline.applicability();
    input.pure.push_back(baseline.pure_records()[selected]);
    for (const auto& pair : baseline.association_records()) {
        if (pair.first_component_id == component_id &&
            pair.second_component_id == component_id) {
            input.association_pairs.push_back(pair);
        }
    }

    const std::vector<th::Component> catalog{
        baseline.components().at(selected)};
    const std::vector<std::string> order{std::string(component_id)};
    return th::CpaParameterSet::create(catalog, order, input);
}

enum class AssociationLimitVariant {
    no_sites,
    no_pairs
};

th::CpaParameterSet association_limit_parameters(
    AssociationLimitVariant variant) {
    const auto baseline = cpa_physical_test::parameters(false);
    th::CpaParameterInput input;
    input.dataset_id = baseline.dataset_id() +
        (variant == AssociationLimitVariant::no_sites
             ? "__nonassociating-limit"
             : "__zero-association-limit");
    input.revision = baseline.revision() + "__association-limit-projection";
    input.applicability = baseline.applicability();
    input.pure.assign(
        baseline.pure_records().begin(), baseline.pure_records().end());
    input.binary.assign(
        baseline.binary_records().begin(), baseline.binary_records().end());
    input.association_pairs.clear();
    if (variant == AssociationLimitVariant::no_sites) {
        for (auto& pure : input.pure) { pure.sites.clear(); }
    }

    const auto selected_components = baseline.components().items();
    const std::vector<th::Component> catalog(
        selected_components.begin(), selected_components.end());
    std::vector<std::string> order;
    order.reserve(catalog.size());
    for (const auto& component : catalog) { order.push_back(component.id); }
    return th::CpaParameterSet::create(catalog, order, input);
}

void run_pure_component_limit() {
    constexpr double temperature_k = cpa_physical_test::temperature_k;
    constexpr std::array<double, 2> densities{{20.0, 200.0}};
    double max_pressure_delta = 0.0;
    double max_mu_delta = 0.0;
    double max_lnphi_delta = 0.0;
    double max_q_value_delta = 0.0;
    double max_abs_association_q = 0.0;

    for (const std::string_view component_id : {std::string_view{"METHANOL"},
                                                 std::string_view{"WATER"}}) {
        const auto parameters = pure_parameters(component_id);
        require(parameters.size() == 1U,
                "pure-component projection did not create one-component model");
        const std::array<double, 1> composition{{1.0}};
        for (const double rho : densities) {
            const auto state = th::evaluate_cpa_phase_at_density(
                temperature_k, rho, composition, parameters);
            require(state.pressure_pa > 0.0,
                    "pure-component limit state must have positive pressure");
            require(state.association.converged(),
                    "pure-component association did not converge");

            const auto derived = derive(
                state.pressure_pa, temperature_k, rho,
                composition, parameters, state.association);
            th::CpaPtRoot analytic_root;
            analytic_root.molar_density_mol_per_m3 = rho;
            th::cpa_detail::cpa_fill_ln_phi(
                state.pressure_pa, temperature_k, composition,
                parameters, state, analytic_root);

            const double direct_q = direct_association_value(
                composition, state.association);
            max_abs_association_q = std::max(
                max_abs_association_q, std::abs(direct_q));
            max_q_value_delta = std::max(
                max_q_value_delta,
                std::abs(derived.association_q_value - direct_q));
            require_abs(
                derived.association_q_value, direct_q,
                1.0e-10 + roundoff(direct_q),
                "pure-component stationary Q mismatch");

            const double p_delta = std::abs(
                derived.pressure_total_pa - state.pressure_pa);
            max_pressure_delta = std::max(max_pressure_delta, p_delta);
            const double pressure_tolerance =
                roundoff(state.pressure_physical_pa) +
                roundoff(state.pressure_association_pa) +
                roundoff(state.pressure_pa);
            require_abs(
                derived.pressure_total_pa, state.pressure_pa,
                pressure_tolerance,
                "pure-component Helmholtz pressure mismatch");

            const double rt = th::cpa_gas_constant_j_per_mol_k * temperature_k;
            const double z = state.pressure_pa / (rho * rt);
            const double analytic_mu = analytic_root.ln_phi.front() + std::log(z);
            const double mu_delta = std::abs(
                derived.mu_total.front() - analytic_mu);
            const double lnphi_delta = std::abs(
                derived.ln_phi.front() - analytic_root.ln_phi.front());
            max_mu_delta = std::max(max_mu_delta, mu_delta);
            max_lnphi_delta = std::max(max_lnphi_delta, lnphi_delta);
            require_abs(
                derived.mu_total.front(), analytic_mu,
                roundoff(analytic_mu) + roundoff(derived.mu_total.front()),
                "pure-component residual chemical potential mismatch");
            require_abs(
                derived.ln_phi.front(), analytic_root.ln_phi.front(),
                roundoff(analytic_root.ln_phi.front()) + roundoff(derived.ln_phi.front()),
                "pure-component ln(phi) mismatch");
        }
    }

    require(max_abs_association_q > 1.0e-8,
            "pure-component regression did not exercise nonzero association");
    std::cout << std::setprecision(17)
              << "CPA_HELMHOLTZ_PURE_LIMIT_OK"
              << " components=2"
              << " densities=2"
              << " max_abs_assoc_Q=" << max_abs_association_q
              << " max_abs_dQ=" << max_q_value_delta
              << " max_abs_dP_pa=" << max_pressure_delta
              << " max_abs_dmu_over_rt=" << max_mu_delta
              << " max_abs_dlnphi=" << max_lnphi_delta
              << '\n';
}

void run_zero_association_limits() {
    constexpr double temperature_k = cpa_physical_test::temperature_k;
    const auto& phase = external::states[2].vapor;
    const double rho = phase.molar_density_mol_per_m3;
    const std::array<double, 2> composition{{0.35, 0.65}};

    const auto no_sites = association_limit_parameters(
        AssociationLimitVariant::no_sites);
    const auto no_pairs = association_limit_parameters(
        AssociationLimitVariant::no_pairs);
    const auto no_sites_state = th::evaluate_cpa_phase_at_density(
        temperature_k, rho, composition, no_sites);
    const auto no_pairs_state = th::evaluate_cpa_phase_at_density(
        temperature_k, rho, composition, no_pairs);

    require(no_sites_state.association.status ==
                th::CpaAssociationStatus::no_associating_sites,
            "non-associating limit did not report no-associating-sites");
    require(no_sites_state.association.sites.empty(),
            "non-associating limit unexpectedly retained site states");
    require(no_pairs_state.association.status == th::CpaAssociationStatus::success,
            "zero-association pair limit did not solve explicit site state");
    require(!no_pairs_state.association.sites.empty(),
            "zero-association pair limit did not retain site classes");
    for (const auto& site : no_pairs_state.association.sites) {
        require_abs(site.unbonded_fraction, 1.0, 0.0,
                    "zero-association limit did not give X=1");
    }

    const auto no_sites_derived = derive(
        no_sites_state.pressure_pa, temperature_k, rho,
        composition, no_sites, no_sites_state.association);
    const auto no_pairs_derived = derive(
        no_pairs_state.pressure_pa, temperature_k, rho,
        composition, no_pairs, no_pairs_state.association);

    for (const auto* derived : {&no_sites_derived, &no_pairs_derived}) {
        require_abs(derived->association_q_value, 0.0, 0.0,
                    "zero-association limit has nonzero association Helmholtz term");
        require_abs(derived->pressure_association_pa, 0.0, 0.0,
                    "zero-association limit has nonzero association pressure");
        require_abs(
            derived->total_value, derived->cubic_value,
            roundoff(derived->cubic_value),
            "zero-association total Helmholtz did not reduce to SRK");
        require_abs(
            derived->pressure_total_pa, derived->pressure_physical_pa,
            roundoff(derived->pressure_physical_pa),
            "zero-association total pressure did not reduce to SRK");
        for (std::size_t i = 0U; i < derived->mu_total.size(); ++i) {
            require_abs(derived->mu_association[i], 0.0, 0.0,
                        "zero-association limit has nonzero association mu");
            require_abs(
                derived->mu_total[i], derived->mu_cubic[i],
                roundoff(derived->mu_cubic[i]),
                "zero-association total mu did not reduce to SRK");
        }
    }

    require_abs(
        no_sites_derived.total_value, no_pairs_derived.total_value,
        roundoff(no_sites_derived.total_value) + roundoff(no_pairs_derived.total_value),
        "no-sites and zero-pairs SRK Helmholtz limits disagree");
    require_abs(
        no_sites_derived.pressure_total_pa, no_pairs_derived.pressure_total_pa,
        roundoff(no_sites_derived.pressure_total_pa) +
            roundoff(no_pairs_derived.pressure_total_pa),
        "no-sites and zero-pairs SRK pressure limits disagree");
    for (std::size_t i = 0U; i < composition.size(); ++i) {
        require_abs(
            no_sites_derived.mu_total[i], no_pairs_derived.mu_total[i],
            roundoff(no_sites_derived.mu_total[i]) +
                roundoff(no_pairs_derived.mu_total[i]),
            "no-sites and zero-pairs SRK chemical-potential limits disagree");
        require_abs(
            no_sites_derived.ln_phi[i], no_pairs_derived.ln_phi[i],
            roundoff(no_sites_derived.ln_phi[i]) +
                roundoff(no_pairs_derived.ln_phi[i]),
            "no-sites and zero-pairs SRK fugacity limits disagree");
    }

    std::cout << std::setprecision(17)
              << "CPA_HELMHOLTZ_ZERO_ASSOC_LIMIT_OK"
              << " rho=" << rho
              << " no_sites_status="
              << static_cast<int>(no_sites_state.association.status)
              << " zero_pairs_status="
              << static_cast<int>(no_pairs_state.association.status)
              << " site_classes_with_zero_pairs="
              << no_pairs_state.association.sites.size()
              << '\n';
}

struct DiluteObservables {
    double density{};
    double abs_f_res_per_mol{};
    double abs_z_minus_one{};
    double max_abs_lnphi{};
};

void run_dilute_gas_limit() {
    constexpr double temperature_k = cpa_physical_test::temperature_k;
    constexpr std::array<double, 4> densities{{1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4}};
    const std::array<double, 2> composition{{0.4, 0.6}};
    const auto parameters = cpa_physical_test::parameters(false);
    std::array<DiluteObservables, densities.size()> observations{};

    for (std::size_t k = 0U; k < densities.size(); ++k) {
        const double rho = densities[k];
        const auto state = th::evaluate_cpa_phase_at_density(
            temperature_k, rho, composition, parameters);
        require(state.pressure_pa > 0.0,
                "dilute-gas limit produced nonpositive pressure");
        const auto derived = derive(
            state.pressure_pa, temperature_k, rho,
            composition, parameters, state.association);
        const double rt = th::cpa_gas_constant_j_per_mol_k * temperature_k;
        const double z = derived.pressure_total_pa / (rho * rt);
        double max_abs_lnphi = 0.0;
        for (const double value : derived.ln_phi) {
            max_abs_lnphi = std::max(max_abs_lnphi, std::abs(value));
        }
        observations[k] = {
            rho,
            std::abs(derived.total_value),
            std::abs(z - 1.0),
            max_abs_lnphi};
    }

    constexpr double decade_ratio_limit = 0.15;
    for (std::size_t k = 1U; k < observations.size(); ++k) {
        const auto& previous = observations[k - 1U];
        const auto& current = observations[k];
        require(
            current.abs_f_res_per_mol <=
                decade_ratio_limit * previous.abs_f_res_per_mol,
            "dilute-gas F_res did not show first-order decade convergence");
        require(
            current.abs_z_minus_one <=
                decade_ratio_limit * previous.abs_z_minus_one,
            "dilute-gas Z did not show first-order decade convergence");
        require(
            current.max_abs_lnphi <=
                decade_ratio_limit * previous.max_abs_lnphi,
            "dilute-gas ln(phi) did not show first-order decade convergence");
    }

    const auto& terminal = observations.back();
    constexpr double terminal_limit = 5.0e-7;
    require(terminal.abs_f_res_per_mol <= terminal_limit,
            "dilute-gas terminal residual Helmholtz is not near zero");
    require(terminal.abs_z_minus_one <= terminal_limit,
            "dilute-gas terminal Z is not near one");
    require(terminal.max_abs_lnphi <= terminal_limit,
            "dilute-gas terminal ln(phi) is not near zero");

    std::cout << std::setprecision(17)
              << "CPA_HELMHOLTZ_DILUTE_LIMIT_OK"
              << " rho_first=" << observations.front().density
              << " rho_last=" << terminal.density
              << " final_abs_Fres=" << terminal.abs_f_res_per_mol
              << " final_abs_Z_minus_1=" << terminal.abs_z_minus_one
              << " final_max_abs_lnphi=" << terminal.max_abs_lnphi
              << " decade_ratio_limit=" << decade_ratio_limit
              << '\n';
}

void compare_permuted_state(
    const external::ThermoPackPhaseKernelState& reference_state,
    const external::PhaseReference& phase,
    double& max_value_delta,
    double& max_pressure_delta,
    double& max_mu_delta,
    double& max_lnphi_delta) {
    const auto normal_parameters = cpa_physical_test::parameters(false);
    const auto swapped_parameters = cpa_physical_test::parameters(true);
    const auto normal_composition = cpa_physical_test::composition(
        phase.composition_methanol, false);
    const auto swapped_composition = cpa_physical_test::composition(
        phase.composition_methanol, true);

    const auto normal_state = th::evaluate_cpa_phase_at_density(
        reference_state.temperature_k, phase.molar_density_mol_per_m3,
        normal_composition, normal_parameters);
    const auto swapped_state = th::evaluate_cpa_phase_at_density(
        reference_state.temperature_k, phase.molar_density_mol_per_m3,
        swapped_composition, swapped_parameters);
    const auto normal = derive(
        reference_state.pressure_pa, reference_state.temperature_k,
        phase.molar_density_mol_per_m3, normal_composition,
        normal_parameters, normal_state.association);
    const auto swapped = derive(
        reference_state.pressure_pa, reference_state.temperature_k,
        phase.molar_density_mol_per_m3, swapped_composition,
        swapped_parameters, swapped_state.association);

    for (const auto values : std::array<std::pair<double, double>, 3>{{
             {normal.cubic_value, swapped.cubic_value},
             {normal.association_q_value, swapped.association_q_value},
             {normal.total_value, swapped.total_value}}}) {
        const double delta = std::abs(values.first - values.second);
        max_value_delta = std::max(max_value_delta, delta);
        require_abs(
            values.first, values.second,
            roundoff(values.first) + roundoff(values.second),
            "component permutation changed residual Helmholtz value");
    }

    for (const auto values : std::array<std::pair<double, double>, 3>{{
             {normal.pressure_physical_pa, swapped.pressure_physical_pa},
             {normal.pressure_association_pa, swapped.pressure_association_pa},
             {normal.pressure_total_pa, swapped.pressure_total_pa}}}) {
        const double delta = std::abs(values.first - values.second);
        max_pressure_delta = std::max(max_pressure_delta, delta);
        require_abs(
            values.first, values.second,
            roundoff(values.first) + roundoff(values.second),
            "component permutation changed Helmholtz-derived pressure");
    }

    for (std::size_t normal_index = 0U; normal_index < 2U; ++normal_index) {
        const std::size_t swapped_index = 1U - normal_index;
        for (const auto values : std::array<std::pair<double, double>, 3>{{
                 {normal.mu_cubic[normal_index], swapped.mu_cubic[swapped_index]},
                 {normal.mu_association[normal_index],
                  swapped.mu_association[swapped_index]},
                 {normal.mu_total[normal_index], swapped.mu_total[swapped_index]}}}) {
            const double delta = std::abs(values.first - values.second);
            max_mu_delta = std::max(max_mu_delta, delta);
            require_abs(
                values.first, values.second,
                roundoff(values.first) + roundoff(values.second),
                "component permutation changed residual chemical potential");
        }
        const double lnphi_delta = std::abs(
            normal.ln_phi[normal_index] - swapped.ln_phi[swapped_index]);
        max_lnphi_delta = std::max(max_lnphi_delta, lnphi_delta);
        require_abs(
            normal.ln_phi[normal_index], swapped.ln_phi[swapped_index],
            roundoff(normal.ln_phi[normal_index]) +
                roundoff(swapped.ln_phi[swapped_index]),
            "component permutation changed ln(phi)");
    }
}

void run_component_permutation_limit() {
    double max_value_delta = 0.0;
    double max_pressure_delta = 0.0;
    double max_mu_delta = 0.0;
    double max_lnphi_delta = 0.0;

    compare_permuted_state(
        external::states[2], external::states[2].liquid,
        max_value_delta, max_pressure_delta, max_mu_delta, max_lnphi_delta);
    compare_permuted_state(
        external::states[4], external::states[4].vapor,
        max_value_delta, max_pressure_delta, max_mu_delta, max_lnphi_delta);

    std::cout << std::setprecision(17)
              << "CPA_HELMHOLTZ_PERMUTATION_OK"
              << " phase_states=2"
              << " max_abs_dF=" << max_value_delta
              << " max_abs_dP_pa=" << max_pressure_delta
              << " max_abs_dmu_over_rt=" << max_mu_delta
              << " max_abs_dlnphi=" << max_lnphi_delta
              << '\n';
}

} // namespace

int main() {
    try {
        run_pure_component_limit();
        run_zero_association_limits();
        run_dilute_gas_limit();
        run_component_permutation_limit();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
