#include <mpmc/thermodynamics/cpa_phase.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {
namespace th = mpmc::thermodynamics;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

th::Provenance source(std::string locator) {
    return {th::SourceKind::synthetic_test,
            "MPMC_HNU CPA structural regression",
            "v1", std::move(locator),
            "Synthetic values exercise CPA contracts/equations only; not physical validation",
            "tests/thermodynamics/cpa_baseline/cpa_baseline_test.cpp",
            "Repository structural regression"};
}

th::CpaSourcedValue value(double number, std::string locator) {
    return {number, source(std::move(locator))};
}

th::Component component(std::string id) {
    const auto provenance = source("component-" + id);
    return {id, id, th::ComponentKind::pure, provenance, std::nullopt};
}

th::Applicability applicability() {
    return {std::nullopt, std::nullopt, source("dataset-applicability")};
}

th::CpaParameterSet one_component(bool associating) {
    const std::vector<th::Component> catalog{component("A")};
    const std::vector<std::string> order{"A"};
    th::CpaParameterInput input;
    input.dataset_id = "synthetic-cpa-one";
    input.revision = "v1";
    input.applicability = applicability();
    th::CpaPureParameterInput pure;
    pure.component_id = "A";
    pure.critical_temperature_k = value(500.0, "Tc-A");
    pure.a0_pa_m6_per_mol2 = value(0.25, "a0-A");
    pure.b_m3_per_mol = value(4.0e-5, "b-A");
    pure.c1_dimensionless = value(0.7, "c1-A");
    if (associating) { pure.sites.push_back({"H", 1U}); }
    input.pure.push_back(pure);
    if (associating) {
        input.association_pairs.push_back({
            "A", "H", "A", "H",
            value(12000.0, "epsilon-A-H-A-H"),
            value(0.03, "beta-A-H-A-H")});
    }
    return th::CpaParameterSet::create(
        catalog, order, input, th::DataPolicy::allow_synthetic_tests);
}

th::CpaParameterSet binary(bool swapped) {
    const std::vector<th::Component> catalog{component("A"), component("B")};
    const std::vector<std::string> order = swapped
        ? std::vector<std::string>{"B", "A"}
        : std::vector<std::string>{"A", "B"};
    th::CpaParameterInput input;
    input.dataset_id = "synthetic-cpa-binary";
    input.revision = "v1";
    input.applicability = applicability();
    input.pure.push_back({
        "A", value(500.0, "Tc-A"), value(0.25, "a0-A"),
        value(4.0e-5, "b-A"), value(0.7, "c1-A"), {{"H", 1U}}});
    input.pure.push_back({
        "B", value(400.0, "Tc-B"), value(0.18, "a0-B"),
        value(5.0e-5, "b-B"), value(0.5, "c1-B"), {}});
    input.binary.push_back({"A", "B", value(0.04, "kij-A-B")});
    input.association_pairs.push_back({
        "A", "H", "A", "H",
        value(10000.0, "epsilon-A-H-A-H"),
        value(0.02, "beta-A-H-A-H")});
    return th::CpaParameterSet::create(
        catalog, order, input, th::DataPolicy::allow_synthetic_tests);
}

void parameter_reordering() {
    const auto first = binary(false);
    const auto second = binary(true);
    require(first.components().at(0).id == "A" && second.components().at(1).id == "A",
            "CPA ordered component snapshot did not reorder");
    require(std::abs(first.pure(0).a0_pa_m6_per_mol2 -
                     second.pure(1).a0_pa_m6_per_mol2) < 1e-15,
            "CPA pure parameters lost component alignment under reordering");
    require(std::abs(first.kij(0, 1) - second.kij(0, 1)) < 1e-15,
            "CPA binary interaction changed under reordering");
    require(first.association_pair(0, "H", 0, "H") != nullptr &&
            second.association_pair(1, "H", 1, "H") != nullptr,
            "CPA site-pair interaction lost component identity under reordering");
}

void runtime_subset_rebuild() {
    const std::vector<th::Component> catalog{
        component("A"), component("B"), component("C")};
    th::CpaParameterInput input;
    input.dataset_id = "synthetic-cpa-superset";
    input.revision = "v1";
    input.applicability = applicability();
    input.pure.push_back({
        "A", value(500.0, "Tc-A"), value(0.25, "a0-A"),
        value(4.0e-5, "b-A"), value(0.7, "c1-A"), {{"H", 1U}}});
    input.pure.push_back({
        "B", value(400.0, "Tc-B"), value(0.18, "a0-B"),
        value(5.0e-5, "b-B"), value(0.5, "c1-B"), {}});
    input.pure.push_back({
        "C", value(450.0, "Tc-C"), value(0.21, "a0-C"),
        value(4.5e-5, "b-C"), value(0.6, "c1-C"), {{"Q", 1U}}});
    input.binary.push_back({"A", "B", value(0.04, "kij-A-B")});
    input.binary.push_back({"A", "C", value(0.02, "kij-A-C")});
    input.binary.push_back({"B", "C", value(0.03, "kij-B-C")});
    input.association_pairs.push_back({
        "A", "H", "A", "H",
        value(10000.0, "epsilon-A-H-A-H"),
        value(0.02, "beta-A-H-A-H")});
    input.association_pairs.push_back({
        "C", "Q", "C", "Q",
        value(9000.0, "epsilon-C-Q-C-Q"),
        value(0.015, "beta-C-Q-C-Q")});

    const std::vector<std::string> order_ab{"A", "B"};
    const std::vector<std::string> order_bc{"B", "C"};
    const auto ab = th::CpaParameterSet::create(
        catalog, order_ab, input, th::DataPolicy::allow_synthetic_tests);
    const auto bc = th::CpaParameterSet::create(
        catalog, order_bc, input, th::DataPolicy::allow_synthetic_tests);

    require(ab.components().at(0).id == "A" && ab.components().at(1).id == "B" &&
                ab.association_records().size() == 1U &&
                ab.association_records()[0].first_component_id == "A" &&
                ab.association_pair(0, "H", 0, "H") != nullptr,
            "CPA A/B subset did not select only its association records");
    require(bc.components().at(0).id == "B" && bc.components().at(1).id == "C" &&
                bc.association_records().size() == 1U &&
                bc.association_records()[0].first_component_id == "C" &&
                bc.association_pair(1, "Q", 1, "Q") != nullptr,
            "CPA B/C subset did not select only its association records");
    require(std::abs(ab.kij(0, 1) - 0.04) < 1e-15 &&
                std::abs(bc.kij(0, 1) - 0.03) < 1e-15,
            "CPA runtime subset rebuild lost selected binary alignment");
}

void provenance_retention() {
    const auto first = binary(false);
    const auto second = binary(true);
    require(first.pure_records().size() == 2U && second.pure_records().size() == 2U,
            "CPA sourced pure records were not retained");
    require(first.pure_records()[0].component_id == "A" &&
                first.pure_records()[0].a0_pa_m6_per_mol2.source.locator == "a0-A" &&
                second.pure_records()[1].component_id == "A" &&
                second.pure_records()[1].a0_pa_m6_per_mol2.source.locator == "a0-A",
            "CPA sourced pure provenance lost ordered component identity");
    require(first.binary_records().size() == 1U &&
                first.binary_records()[0].kij_dimensionless.source.locator == "kij-A-B",
            "CPA binary provenance was not retained");
    require(first.association_records().size() == 1U &&
                first.association_records()[0].epsilon_j_per_mol.source.locator ==
                    "epsilon-A-H-A-H" &&
                first.association_records()[0].beta_dimensionless.source.locator ==
                    "beta-A-H-A-H",
            "CPA association provenance was not retained");
}

void analytic_one_site_association() {
    const auto parameters = one_component(true);
    const std::vector<double> x{1.0};
    constexpr double temperature = 320.0;
    constexpr double density = 8000.0;
    const auto result = th::solve_cpa_association(
        temperature, density, x, parameters);
    require(result.status == th::CpaAssociationStatus::success &&
                result.sites.size() == 1U,
            "one-site CPA association did not converge");

    const auto* pair = parameters.association_pair(0, "H", 0, "H");
    require(pair != nullptr, "one-site association pair missing");
    const double delta = result.radial_distribution *
        std::expm1(pair->epsilon_j_per_mol /
                   (th::cpa_gas_constant_j_per_mol_k * temperature)) *
        parameters.pure(0).b_m3_per_mol * pair->beta_dimensionless;
    const double q = density * delta;
    const double exact = 2.0 / (1.0 + std::sqrt(1.0 + 4.0 * q));
    require(std::abs(result.sites[0].unbonded_fraction - exact) < 3e-12,
            "CPA one-site fixed point disagrees with independent quadratic solution");
}

double association_helmholtz_rt(
    double temperature, double density,
    const std::vector<double>& x,
    const th::CpaParameterSet& parameters) {
    const auto association = th::solve_cpa_association(
        temperature, density, x, parameters);
    require(association.status == th::CpaAssociationStatus::success,
            "association Helmholtz reference did not converge");
    double value_rt = 0.0;
    for (const auto& site : association.sites) {
        const double site_x = site.unbonded_fraction;
        value_rt += x[site.component_index] *
            static_cast<double>(site.multiplicity) *
            (std::log(site_x) - 0.5 * site_x + 0.5);
    }
    return value_rt;
}

void association_pressure_helmholtz_crosscheck() {
    const auto parameters = one_component(true);
    const std::vector<double> x{1.0};
    constexpr double temperature = 320.0;
    constexpr double density = 8000.0;
    const auto state = th::evaluate_cpa_phase_at_density(
        temperature, density, x, parameters);
    const double h = density * 1.0e-5;
    const double plus = association_helmholtz_rt(
        temperature, density + h, x, parameters);
    const double minus = association_helmholtz_rt(
        temperature, density - h, x, parameters);
    const double pressure_reference =
        density * density * th::cpa_gas_constant_j_per_mol_k * temperature *
        (plus - minus) / (2.0 * h);
    require(std::abs(state.pressure_association_pa - pressure_reference) <=
                2.0e-8 * std::max(1.0, std::abs(pressure_reference)),
            "CPA association pressure disagrees with independent Helmholtz density derivative");
}

void nonassociating_srk_limit() {
    const auto parameters = one_component(false);
    const std::vector<double> x{1.0};
    constexpr double temperature = 350.0;
    constexpr double density = 1500.0;
    const auto state = th::evaluate_cpa_phase_at_density(
        temperature, density, x, parameters);
    require(state.association.status == th::CpaAssociationStatus::no_associating_sites,
            "non-associating CPA component did not expose SRK limit");
    require(state.pressure_association_pa == 0.0,
            "non-associating CPA state has nonzero association pressure");

    const double tr = temperature / parameters.pure(0).critical_temperature_k;
    const double base = 1.0 + parameters.pure(0).c1_dimensionless *
        (1.0 - std::sqrt(tr));
    const double a = parameters.pure(0).a0_pa_m6_per_mol2 * base * base;
    const double b = parameters.pure(0).b_m3_per_mol;
    const double brho = b * density;
    const double expected = th::cpa_gas_constant_j_per_mol_k * temperature * density /
        (1.0 - brho) - a * density * density / (1.0 + brho);
    require(std::abs(state.pressure_pa - expected) <=
                2e-13 * std::max(1.0, std::abs(expected)),
            "CPA non-associating pressure does not reduce to configured SRK physical term");
}

void association_pressure_and_permutation() {
    const auto first = binary(false);
    const auto second = binary(true);
    const std::vector<double> x_first{0.7, 0.3};
    const std::vector<double> x_second{0.3, 0.7};
    const auto a = th::evaluate_cpa_phase_at_density(
        330.0, 5000.0, x_first, first);
    const auto b = th::evaluate_cpa_phase_at_density(
        330.0, 5000.0, x_second, second);
    require(a.association.status == th::CpaAssociationStatus::success &&
                a.pressure_association_pa < 0.0,
            "associating CPA mixture did not produce a converged negative association pressure contribution");
    require(std::abs(a.pressure_pa - b.pressure_pa) <=
                2e-12 * std::max(1.0, std::abs(a.pressure_pa)),
            "CPA phase pressure changed under runtime component permutation");
}

void missing_binary_rejected() {
    const std::vector<th::Component> catalog{component("A"), component("B")};
    const std::vector<std::string> order{"A", "B"};
    th::CpaParameterInput input;
    input.dataset_id = "missing-binary";
    input.revision = "v1";
    input.applicability = applicability();
    input.pure.push_back({
        "A", value(500.0, "Tc-A"), value(0.25, "a0-A"),
        value(4e-5, "b-A"), value(0.7, "c1-A"), {}});
    input.pure.push_back({
        "B", value(400.0, "Tc-B"), value(0.18, "a0-B"),
        value(5e-5, "b-B"), value(0.5, "c1-B"), {}});
    bool caught = false;
    try {
        (void)th::CpaParameterSet::create(
            catalog, order, input, th::DataPolicy::allow_synthetic_tests);
    } catch (const th::ContractError& error) {
        caught = error.code() == th::ContractErrorCode::missing_parameter;
    }
    require(caught, "CPA parameter contract accepted a missing binary interaction");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"parameter_reordering", parameter_reordering},
    {"runtime_subset", runtime_subset_rebuild},
    {"provenance_retention", provenance_retention},
    {"analytic_one_site", analytic_one_site_association},
    {"association_pressure", association_pressure_helmholtz_crosscheck},
    {"nonassociating_srk", nonassociating_srk_limit},
    {"association_permutation", association_pressure_and_permutation},
    {"missing_binary", missing_binary_rejected}};

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
