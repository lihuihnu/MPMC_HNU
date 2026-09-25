#include <mpmc/thermodynamics/cpa_pt_phase.hpp>
#include <mpmc/thermodynamics/selected_phase_fugacity.hpp>
#include <mpmc/thermodynamics/selected_phase_density.hpp>

#include "test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

double pure_a(double temperature_k, const th::CpaPureParameters& pure) {
    const double base = 1.0 + pure.c1_dimensionless *
        (1.0 - std::sqrt(temperature_k / pure.critical_temperature_k));
    return pure.a0_pa_m6_per_mol2 * base * base;
}

struct ReferenceMix {
    double a{};
    double b{};
};

ReferenceMix reference_mix(double temperature_k, const Vec& x,
                           const th::CpaParameterSet& parameters) {
    Vec a_i(parameters.size(), 0.0);
    for (std::size_t i = 0; i < parameters.size(); ++i) {
        a_i[i] = pure_a(temperature_k, parameters.pure(i));
    }
    ReferenceMix mixed;
    for (std::size_t i = 0; i < parameters.size(); ++i) {
        mixed.b += x[i] * parameters.pure(i).b_m3_per_mol;
        for (std::size_t j = 0; j < parameters.size(); ++j) {
            mixed.a += x[i] * x[j] * std::sqrt(a_i[i] * a_i[j]) *
                       (1.0 - parameters.kij(i, j));
        }
    }
    return mixed;
}

double analytic_fixture_site_fraction(double temperature_k, double density,
                                      const Vec& x,
                                      const th::CpaParameterSet& parameters) {
    const std::size_t a_index = parameters.components().index_of("A");
    const auto mixed = reference_mix(temperature_k, x, parameters);
    const double eta = mixed.b * density / 4.0;
    const double g = 1.0 / (1.0 - 1.9 * eta);
    const auto* pair = parameters.association_pair(
        a_index, "H", a_index, "H");
    require(pair != nullptr, "binary fixture lost A-H self association");
    const double delta = g *
        std::expm1(pair->epsilon_j_per_mol /
                   (th::cpa_gas_constant_j_per_mol_k * temperature_k)) *
        parameters.pure(a_index).b_m3_per_mol * pair->beta_dimensionless;
    const double q = density * x[a_index] * delta;
    return 2.0 / (1.0 + std::sqrt(1.0 + 4.0 * q));
}

double reference_residual_helmholtz_rt(
    double temperature_k, double density, const Vec& x,
    const th::CpaParameterSet& parameters) {
    const auto mixed = reference_mix(temperature_k, x, parameters);
    const double b_rho = mixed.b * density;
    const double cubic = -std::log1p(-b_rho) -
        mixed.a / (mixed.b * th::cpa_gas_constant_j_per_mol_k * temperature_k) *
        std::log1p(b_rho);
    const std::size_t a_index = parameters.components().index_of("A");
    const double site_x = analytic_fixture_site_fraction(
        temperature_k, density, x, parameters);
    const double association = x[a_index] *
        (std::log(site_x) - 0.5 * site_x + 0.5);
    return cubic + association;
}

double reference_pressure(double temperature_k, double density, const Vec& x,
                          const th::CpaParameterSet& parameters) {
    const auto mixed = reference_mix(temperature_k, x, parameters);
    const double b_rho = mixed.b * density;
    const double rt = th::cpa_gas_constant_j_per_mol_k * temperature_k;
    const double physical = rt * density / (1.0 - b_rho) -
        mixed.a * density * density / (1.0 + b_rho);
    const std::size_t a_index = parameters.components().index_of("A");
    const double site_x = analytic_fixture_site_fraction(
        temperature_k, density, x, parameters);
    const double eta = b_rho / 4.0;
    const double rho_dln_g = (1.9 * eta) / (1.0 - 1.9 * eta);
    const double association = -0.5 * rt * density *
        (1.0 + rho_dln_g) * x[a_index] * (1.0 - site_x);
    return physical + association;
}

const th::CpaPtRoot& nearest_root(const th::CpaPtRootSet& roots, double density) {
    require(!roots.roots.empty(), "CPA PT root set unexpectedly empty");
    const auto found = std::min_element(
        roots.roots.begin(), roots.roots.end(),
        [density](const auto& first, const auto& second) {
            return std::abs(first.molar_density_mol_per_m3 - density) <
                   std::abs(second.molar_density_mol_per_m3 - density);
        });
    return *found;
}

void nonassociating_srk_three_roots() {
    const auto parameters = cpa_pt_test::nonassociating_pure();
    const auto phase = th::CpaPtPhase::from_parameters(parameters);
    const Vec x{1.0};
    constexpr double pressure = 1.0e5;
    constexpr double temperature = 200.0;
    const auto roots = phase.roots(pressure, temperature, x);
    require(roots.status == th::CpaPtRootStatus::success,
            "non-associating SRK root scan did not resolve simple roots");
    require(roots.roots.size() == 3U,
            "CPA non-associating limit did not recover all three SRK roots");
    require(roots.roots[0].pressure_slope_sign == 1 &&
                roots.roots[1].pressure_slope_sign == -1 &&
                roots.roots[2].pressure_slope_sign == 1,
            "CPA SRK-limit roots lost stable/unstable crossing topology");

    const double a = pure_a(temperature, parameters.pure(0));
    const double b = parameters.pure(0).b_m3_per_mol;
    const double rt = th::cpa_gas_constant_j_per_mol_k * temperature;
    const double A = a * pressure / (rt * rt);
    const double B = b * pressure / rt;
    for (const auto& root : roots.roots) {
        const double z = root.compressibility_factor;
        const double polynomial =
            ((z - 1.0) * z + (A - B - B * B)) * z - A * B;
        require(std::abs(polynomial) < 2.0e-9,
                "CPA SRK-limit density root disagrees with independent cubic polynomial");
        const double expected_ln_phi =
            (z - 1.0) - std::log(z - B) -
            (A / B) * std::log1p(B / z);
        require(root.ln_phi.size() == 1U &&
                    std::abs(root.ln_phi[0] - expected_ln_phi) < 2.0e-8,
                "CPA non-associating ln(phi) disagrees with independent SRK expression");
    }
}

void associating_helmholtz_chemical_potential() {
    const auto parameters = cpa_pt_test::associating_binary(false);
    const auto phase = th::CpaPtPhase::from_parameters(parameters);
    constexpr double temperature = 330.0;
    constexpr double density = 5000.0;
    const Vec x{0.7, 0.3};
    const double pressure = reference_pressure(temperature, density, x, parameters);
    require(std::isfinite(pressure) && pressure > 0.0,
            "independent CPA fixture pressure is invalid");

    const auto roots = phase.roots(pressure, temperature, x);
    require(roots.status == th::CpaPtRootStatus::success,
            "associating CPA finite root scan did not resolve simple roots");
    const auto& root = nearest_root(roots, density);
    require(std::abs(root.molar_density_mol_per_m3 - density) <=
                2.0e-8 * density,
            "CPA associating root left the independent density anchor");

    constexpr double total_moles = 1.0;
    const double volume = total_moles / density;
    const double z = pressure /
        (density * th::cpa_gas_constant_j_per_mol_k * temperature);
    constexpr double h = 1.0e-6;
    for (std::size_t component = 0; component < x.size(); ++component) {
        Vec n_plus = x;
        Vec n_minus = x;
        n_plus[component] += h;
        n_minus[component] -= h;
        const double total_plus = 1.0 + h;
        const double total_minus = 1.0 - h;
        Vec x_plus(x.size(), 0.0);
        Vec x_minus(x.size(), 0.0);
        for (std::size_t i = 0; i < x.size(); ++i) {
            x_plus[i] = n_plus[i] / total_plus;
            x_minus[i] = n_minus[i] / total_minus;
        }
        const double rho_plus = total_plus / volume;
        const double rho_minus = total_minus / volume;
        const double extensive_plus = total_plus *
            reference_residual_helmholtz_rt(
                temperature, rho_plus, x_plus, parameters);
        const double extensive_minus = total_minus *
            reference_residual_helmholtz_rt(
                temperature, rho_minus, x_minus, parameters);
        const double mu_residual =
            (extensive_plus - extensive_minus) / (2.0 * h);
        const double expected_ln_phi = mu_residual - std::log(z);
        require(root.ln_phi.size() == x.size() &&
                    std::abs(root.ln_phi[component] - expected_ln_phi) < 3.0e-7,
                "CPA ln(phi) disagrees with independent fixed-T,V Helmholtz derivative");
    }
}

void associating_component_permutation() {
    const auto first_parameters = cpa_pt_test::associating_binary(false);
    const auto second_parameters = cpa_pt_test::associating_binary(true);
    const auto first_phase = th::CpaPtPhase::from_parameters(first_parameters);
    const auto second_phase = th::CpaPtPhase::from_parameters(second_parameters);
    constexpr double temperature = 330.0;
    constexpr double density = 5000.0;
    const Vec first_x{0.7, 0.3};
    const Vec second_x{0.3, 0.7};
    const double pressure = reference_pressure(
        temperature, density, first_x, first_parameters);
    const auto first = first_phase.roots(pressure, temperature, first_x);
    const auto second = second_phase.roots(pressure, temperature, second_x);
    require(first.status == th::CpaPtRootStatus::success &&
                second.status == th::CpaPtRootStatus::success,
            "CPA component permutation changed root-search status");
    const auto& a = nearest_root(first, density);
    const auto& b = nearest_root(second, density);
    require(std::abs(a.molar_density_mol_per_m3 - b.molar_density_mol_per_m3) <=
                2.0e-10 * density,
            "CPA density root changed under component permutation");
    require(a.ln_phi.size() == 2U && b.ln_phi.size() == 2U &&
                std::abs(a.ln_phi[0] - b.ln_phi[1]) < 2.0e-10 &&
                std::abs(a.ln_phi[1] - b.ln_phi[0]) < 2.0e-10,
            "CPA fugacity coefficients changed under component permutation");
}


void selected_phase_fugacity_contract() {
    const auto parameters = cpa_pt_test::nonassociating_pure();
    const auto phase = th::CpaPtPhase::from_parameters(parameters);
    const Vec x{1.0};
    constexpr double pressure = 1.0e5;
    constexpr double temperature = 200.0;
    const auto roots = phase.roots(pressure, temperature, x);
    require(roots.status == th::CpaPtRootStatus::success &&
                roots.roots.size() == 3U,
            "CPA selected-phase fugacity root fixture");
    constexpr std::size_t root_index = 2U;
    const auto wrapped = th::evaluate_selected_phase_fugacity(
        phase, pressure, temperature, std::span<const double>{x},
        th::CpaSelectedPhase{root_index, {}});
    require(wrapped.ln_phi.size() == roots.roots[root_index].ln_phi.size(),
            "CPA selected-phase fugacity size changed");
    for (std::size_t i = 0U; i < wrapped.ln_phi.size(); ++i) {
        require(std::abs(wrapped.ln_phi[i] - roots.roots[root_index].ln_phi[i]) < 1.0e-14,
                "CPA selected-phase fugacity changed selected root value");
    }
    bool caught = false;
    try {
        (void)th::evaluate_selected_phase_fugacity(
            phase, pressure, temperature, std::span<const double>{x},
            th::CpaSelectedPhase{3U, {}});
    } catch (const std::out_of_range&) {
        caught = true;
    }
    require(caught, "CPA selected-phase fugacity silently changed invalid root index");
    static_assert(
        th::SelectedPhaseFugacityCapabilities<th::CpaPtPhase>::derivative_support ==
        th::SelectedPhaseFugacityDerivativeSupport::scalar_generic_first_order);
}


void selected_phase_density_contract() {
    const auto parameters = cpa_pt_test::associating_binary(false);
    const auto phase = th::CpaPtPhase::from_parameters(parameters);
    constexpr double temperature = 330.0;
    constexpr double density = 5000.0;
    constexpr double x0 = 0.7;
    const Vec composition{x0, 1.0 - x0};
    const double pressure = reference_pressure(
        temperature, density, composition, parameters);
    const auto roots = phase.roots(
        pressure, temperature, composition);
    require(
        roots.status == th::CpaPtRootStatus::success,
        "CPA selected density fixture has no root");
    const auto found = std::min_element(
        roots.roots.begin(),
        roots.roots.end(),
        [](const auto& first, const auto& second) {
            return std::abs(
                       first.molar_density_mol_per_m3 -
                       5000.0) <
                   std::abs(
                       second.molar_density_mol_per_m3 -
                       5000.0);
        });
    const std::size_t root_index =
        static_cast<std::size_t>(
            std::distance(
                roots.roots.begin(),
                found));

    const auto primal =
        th::evaluate_selected_phase_molar_density(
            phase,
            pressure,
            temperature,
            std::span<const double>{composition},
            th::CpaSelectedPhase{
                root_index, {}});
    require(
        std::abs(
            primal.molar_density_mol_per_m3 -
            roots.roots[root_index]
                .molar_density_mol_per_m3) <
            1.0e-12 *
                roots.roots[root_index]
                    .molar_density_mol_per_m3,
        "CPA selected density changed selected root primal");

    using D = mpmc::ad::Dual<double, 3U>;
    const std::array<D, 2> composition_ad{
        D::variable(x0, 2U),
        D{
            1.0 - x0,
            D::Gradient{
                0.0, 0.0, -1.0}}};
    const auto differentiated =
        th::evaluate_selected_phase_molar_density(
            phase,
            D::variable(pressure, 0U),
            D::variable(temperature, 1U),
            std::span<const D>{
                composition_ad},
            th::CpaSelectedPhase{
                root_index, {}});
    require(
        std::abs(
            differentiated
                .molar_density_mol_per_m3
                .value() -
            primal.molar_density_mol_per_m3) <
            1.0e-12 *
                primal.molar_density_mol_per_m3,
        "CPA selected density derivative changed primal");
    for (std::size_t lane = 0U;
         lane < 3U;
         ++lane) {
        require(
            std::isfinite(
                differentiated
                    .molar_density_mol_per_m3
                    .derivative(lane)),
            "CPA selected density IFT derivative is non-finite");
    }
}

std::size_t nearest_root_index(
    const th::CpaPtRootSet& roots,
    double density) {
    require(!roots.roots.empty(),
            "CPA PT root set unexpectedly empty");
    const auto found = std::min_element(
        roots.roots.begin(),
        roots.roots.end(),
        [density](const auto& first, const auto& second) {
            return std::abs(first.molar_density_mol_per_m3 - density) <
                   std::abs(second.molar_density_mol_per_m3 - density);
        });
    return static_cast<std::size_t>(
        std::distance(roots.roots.begin(), found));
}

void associating_selected_phase_derivatives() {
    const auto parameters = cpa_pt_test::associating_binary(false);
    const auto phase = th::CpaPtPhase::from_parameters(parameters);
    constexpr double temperature = 330.0;
    constexpr double density = 5000.0;
    constexpr double x0 = 0.7;
    const Vec composition{x0, 1.0 - x0};
    const double pressure = reference_pressure(
        temperature, density, composition, parameters);
    const auto roots = phase.roots(
        pressure, temperature, composition);
    require(roots.status == th::CpaPtRootStatus::success,
            "CPA derivative fixture has no resolved PT roots");
    const std::size_t root_index =
        nearest_root_index(roots, density);

    using D = mpmc::ad::Dual<double, 3U>;
    const D pressure_ad = D::variable(pressure, 0U);
    const D temperature_ad = D::variable(temperature, 1U);
    const std::array<D, 2> composition_ad{
        D::variable(x0, 2U),
        D{1.0 - x0, D::Gradient{0.0, 0.0, -1.0}}};

    const auto evaluated = th::evaluate_selected_phase_fugacity(
        phase,
        pressure_ad,
        temperature_ad,
        std::span<const D>{composition_ad},
        th::CpaSelectedPhase{root_index, {}});

    require(evaluated.ln_phi.size() == 2U,
            "CPA selected derivative fugacity size changed");
    for (std::size_t component = 0U; component < 2U; ++component) {
        require(std::abs(evaluated.ln_phi[component].value() -
                         roots.roots[root_index].ln_phi[component]) < 1.0e-14,
                "CPA selected derivative changed primal ln(phi)");
    }

    const auto fresh_ln_phi =
        [&](double p, double temp, double first_fraction,
            std::size_t component) {
            const Vec x{first_fraction, 1.0 - first_fraction};
            const auto perturbed = phase.roots(p, temp, x);
            require(perturbed.status == th::CpaPtRootStatus::success,
                    "CPA fresh derivative perturbation lost root set");
            require(root_index < perturbed.roots.size(),
                    "CPA fresh derivative perturbation changed selected root count");
            return perturbed.roots[root_index].ln_phi[component];
        };

    const std::array<double, 3> steps{
        std::max(100.0, 1.0e-5 * pressure),
        1.0e-3,
        1.0e-5};

    for (std::size_t component = 0U; component < 2U; ++component) {
        const double fd_pressure =
            (fresh_ln_phi(pressure + steps[0], temperature, x0, component) -
             fresh_ln_phi(pressure - steps[0], temperature, x0, component)) /
            (2.0 * steps[0]);
        const double fd_temperature =
            (fresh_ln_phi(pressure, temperature + steps[1], x0, component) -
             fresh_ln_phi(pressure, temperature - steps[1], x0, component)) /
            (2.0 * steps[1]);
        const double fd_composition =
            (fresh_ln_phi(pressure, temperature, x0 + steps[2], component) -
             fresh_ln_phi(pressure, temperature, x0 - steps[2], component)) /
            (2.0 * steps[2]);

        const std::array<double, 3> expected{
            fd_pressure,
            fd_temperature,
            fd_composition};
        const std::array<double, 3> absolute_tolerances{
            5.0e-11,
            5.0e-6,
            5.0e-5};
        for (std::size_t lane = 0U; lane < expected.size(); ++lane) {
            const double actual =
                evaluated.ln_phi[component].derivative(lane);
            const double tolerance =
                absolute_tolerances[lane] +
                2.0e-3 * std::abs(expected[lane]);
            require(std::isfinite(actual) &&
                        std::abs(actual - expected[lane]) <= tolerance,
                    "CPA selected-phase analytic/IFT derivative disagrees with fresh PT re-solves");
        }
    }
}

void evaluation_budget_is_explicit() {
    const auto parameters = cpa_pt_test::nonassociating_pure();
    const auto phase = th::CpaPtPhase::from_parameters(parameters);
    const Vec x{1.0};
    th::CpaPtOptions options;
    options.scan_intervals = 16U;
    options.max_evaluations = 16U;
    const auto roots = phase.roots(1.0e5, 200.0, x, options);
    require(roots.status == th::CpaPtRootStatus::evaluation_limit,
            "CPA PT root solver hid an exhausted property-evaluation budget");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"nonassociating_three_roots", nonassociating_srk_three_roots},
    {"associating_helmholtz", associating_helmholtz_chemical_potential},
    {"component_permutation", associating_component_permutation},
    {"selected_phase_fugacity", selected_phase_fugacity_contract},
    {"selected_phase_density", selected_phase_density_contract},
    {"selected_phase_derivatives", associating_selected_phase_derivatives},
    {"evaluation_budget", evaluation_budget_is_explicit}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) { throw std::invalid_argument("one test name required"); }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                std::cout << "[PASS] " << name << '\n';
                return 0;
            }
        }
        throw std::invalid_argument("unknown test name");
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
