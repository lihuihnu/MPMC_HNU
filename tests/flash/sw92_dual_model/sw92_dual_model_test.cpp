#include <mpmc/flash/sw92_dual_model.hpp>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <optional>
#include <source_location>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool sw92_dual_model_header();

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " + std::string(message));
    }
}

void near(double actual, long double expected, long double relative = 8e-9L,
          long double absolute = 8e-12L,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(static_cast<long double>(actual) - expected) >
            absolute + relative * std::abs(expected)) {
        require(false, "reference mismatch", where);
    }
}

th::Provenance paper(std::string locator) {
    return {th::SourceKind::literature,
            "doi:10.1016/0378-3812(92)85105-H",
            "Fluid Phase Equilibria 77 (1992) 217-240; authors' errata in supplied PDF",
            std::move(locator),
            "SW92 corrected-original dual-model regression; not experimental validation",
            "User PDF sha256:cb5b1d5034d78d934e887449ce0c89692d43431835d1371da5606d95b2c6bf58",
            "Bibliographic/formula facts only; source PDF is not redistributed"};
}

th::Provenance synthetic(std::string locator) {
    return {th::SourceKind::synthetic_test,
            "MPMC_HNU SW92 dual-model structural fixture", "v1", std::move(locator),
            "Artificial gas/gas BIP only for multicomponent phase-label contract testing",
            "tests/flash/sw92_dual_model/sw92_dual_model_test.cpp",
            "Repository test data; no physical-validation claim"};
}

th::SourcedScalar sourced(double value, th::Unit unit, const th::Provenance& source,
                          std::string original = "SI or dimensionless",
                          std::string conversion = "identity") {
    return {value, unit, source, std::move(original), std::move(conversion)};
}

struct GasSpec {
    const char* id;
    const char* display;
    th::Sw92Species species;
    double tc;
    double pc_bar;
    double omega;
    double nonaqueous_kij;
};

constexpr GasSpec co2{"carbon-dioxide", "Carbon dioxide", th::Sw92Species::carbon_dioxide,
                      304.2, 73.8, 0.2273, 0.1896};
constexpr GasSpec methane{"methane", "Methane", th::Sw92Species::hydrocarbon,
                          190.6, 46.0, 0.0108, 0.4850};
constexpr GasSpec water_spec{"water", "Water", th::Sw92Species::water,
                             647.3, 221.2, 0.3434, 0.0};

void add_pure(th::Sw92ParameterInput& input, std::vector<th::Component>& catalog,
              const GasSpec& spec) {
    const auto identity = paper("Table 3 component identities");
    const auto properties = paper("Table 3 physical properties");
    catalog.push_back({spec.id, spec.display, th::ComponentKind::pure, identity, std::nullopt});
    input.pure.push_back({spec.id, spec.species,
        sourced(spec.tc, th::Unit::kelvin, properties, "K", "identity"),
        sourced(spec.pc_bar * 100000.0, th::Unit::pascal, properties,
                "bar", "bar * 100000 -> Pa"),
        sourced(spec.omega, th::Unit::dimensionless, properties)});
}

void add_water_pair(th::Sw92ParameterInput& input, const GasSpec& gas) {
    input.water_binary.push_back({gas.id, th::Sw92NonAqueousWaterRule::constant,
        sourced(gas.nonaqueous_kij, th::Unit::dimensionless,
                paper("Table 5 non-aqueous water BIP"))});
}

th::Sw92Phase<double> binary_model(const GasSpec& gas, bool reverse = false) {
    std::vector<th::Component> catalog;
    th::Sw92ParameterInput input;
    input.model_id = std::string(th::sw92_corrected_profile);
    input.dataset_id = std::string("SW92-Table3-Table5-") + gas.id + "-dual-model";
    input.revision = "user-PDF-sha256-cb5b1d50";
    input.applicability = {{std::nullopt, std::nullopt,
                            paper("Eqs. (9)-(17), Tables 2-5; no global T/p bounds inferred")},
                           std::nullopt};
    add_pure(input, catalog, gas);
    add_pure(input, catalog, water_spec);
    add_water_pair(input, gas);
    const std::vector<std::string> order = reverse
        ? std::vector<std::string>{water_spec.id, gas.id}
        : std::vector<std::string>{gas.id, water_spec.id};
    return th::Sw92Phase<double>::from_parameters(
        th::Sw92ParameterSet::create(catalog, order, input));
}

th::Sw92Phase<double> ternary_structural_model() {
    std::vector<th::Component> catalog;
    th::Sw92ParameterInput input;
    input.model_id = std::string(th::sw92_corrected_profile);
    input.dataset_id = "SW92-dual-model-ternary-structural";
    input.revision = "synthetic-nonwater-bip-v1";
    input.applicability = {{std::nullopt, std::nullopt,
                            paper("SW92 equations; structural routing test only")},
                           std::nullopt};
    add_pure(input, catalog, co2);
    add_pure(input, catalog, methane);
    add_pure(input, catalog, water_spec);
    add_water_pair(input, co2);
    add_water_pair(input, methane);
    const auto source = synthetic("CO2/CH4 family-specific non-water pair");
    input.nonwater_binary.push_back({co2.id, methane.id,
        sourced(0.0, th::Unit::dimensionless, source),
        sourced(0.0, th::Unit::dimensionless, source)});
    const std::vector<std::string> order{co2.id, methane.id, water_spec.id};
    return th::Sw92Phase<double>::from_parameters(
        th::Sw92ParameterSet::create(catalog, order, input,
                                     th::DataPolicy::allow_synthetic_tests));
}

struct DualGolden {
    const GasSpec* gas;
    double pressure_pa;
    double temperature_k;
    double molality;
    double feed_gas;
    long double aqueous_liquid_gas;
    long double aqueous_vapor_gas;
    long double aqueous_vapor_fraction;
    long double nonaqueous_liquid_gas;
    long double nonaqueous_vapor_gas;
    long double nonaqueous_vapor_fraction;
    long double cross_ratio_gas;
    long double cross_ratio_water;
};

constexpr DualGolden golden[] = {
    {&co2, 3.0e6, 340.0, 0.0, 0.7,
     0.005936172024213170056428932801071659387598344017016019L,
     0.98741239776990554564891282081338673931115242042601548L,
     0.70716315868829186790522005951676476191789314445481349L,
     0.00028414840081531502922309549898233567384253874407420669L,
     0.98872581451967009667119228772574961667799726909104273L,
     0.70789797272168778243120726503517591720856121160092583L,
     166.55949498881378773966141761680876885103599239182512L,
     0.011341510638494450131621101780021805693331383211128272L},
    {&methane, 1.0e7, 350.0, 1.0, 0.5,
     0.00090687277385112663805494943543515398677630890790446910L,
     0.99043450340463085644642841949477174162636913432131321L,
     0.50437512988697371210319463020514765271259334999788643L,
     0.0000061874447993037728438291445940435604668979585134072584L,
     0.99467372174749769026766292177496815969371532754437181L,
     0.50267430604912253510194884751912280998048583603410677L,
     1096.8172718688151563916858388328206926482921644799898L,
     0.0053311128936398784351512017990406043727637758682054346L}
};

void check_dual(const fl::Sw92WhitsonDualModelResult& result,
                const DualGolden& reference, bool reverse = false) {
    require(result.status == fl::Sw92DualModelStatus::observables_available &&
                result.compatibility_observables_available(),
            "dual-model observables were not available");
    require(result.equilibrium_algorithm == fl::sw92_whitson_dual_model_algorithm &&
                result.model_profile == th::sw92_corrected_profile &&
                result.phase_convention == th::sw92_pt_convention &&
                result.nacl_molality_mol_per_kg_water == reference.molality,
            "outer dual-model metadata changed");
    require(result.aqueous_run.equilibrium_algorithm == fl::sw92_family_vle_algorithm &&
                result.nonaqueous_run.equilibrium_algorithm == fl::sw92_family_vle_algorithm &&
                result.aqueous_run.family == th::SwPhaseFamily::aqueous &&
                result.nonaqueous_run.family == th::SwPhaseFamily::nonaqueous,
            "nested fixed-family identities changed");
    require(result.aqueous_run.dataset_id == result.nonaqueous_run.dataset_id &&
                result.aqueous_run.revision == result.nonaqueous_run.revision &&
                result.aqueous_run.component_ids == result.nonaqueous_run.component_ids,
            "family runs did not preserve the same ordered snapshot");
    require(result.aqueous_run.solution.status == fl::PtSplitStatus::two_phase_no_instability_found &&
                result.nonaqueous_run.solution.status == fl::PtSplitStatus::two_phase_no_instability_found,
            "one of the independent family runs was not accepted");

    const auto* aq = result.aqueous_run.solution.candidate();
    const auto* na = result.nonaqueous_run.solution.candidate();
    require(aq != nullptr && na != nullptr, "accepted family run lacks a candidate");
    const std::size_t gas_index = reverse ? 1U : 0U;
    const std::size_t water_index = reverse ? 0U : 1U;
    near(aq->fractions.liquid[gas_index], reference.aqueous_liquid_gas);
    near(aq->fractions.vapor[gas_index], reference.aqueous_vapor_gas);
    near(aq->fractions.vapor_fraction, reference.aqueous_vapor_fraction);
    near(na->fractions.liquid[gas_index], reference.nonaqueous_liquid_gas);
    near(na->fractions.vapor[gas_index], reference.nonaqueous_vapor_gas);
    near(na->fractions.vapor_fraction, reference.nonaqueous_vapor_fraction);

    const auto& observables = *result.observables;
    require(observables.water_index == water_index &&
                observables.aqueous.source_role == fl::PtPhaseRole::liquid_candidate &&
                observables.nonaqueous.source_role == fl::PtPhaseRole::vapor_candidate,
            "binary target-phase routing changed");
    near(observables.aqueous.composition[gas_index], reference.aqueous_liquid_gas);
    near(observables.nonaqueous.composition[gas_index], reference.nonaqueous_vapor_gas);
    near(observables.cross_model_equilibrium_ratio[gas_index], reference.cross_ratio_gas);
    near(observables.cross_model_equilibrium_ratio[water_index], reference.cross_ratio_water);

    // The two complete model-pass phase fractions are intentionally independent.
    require(std::abs(aq->fractions.vapor_fraction - na->fractions.vapor_fraction) > 1e-8,
            "dual-model runs were accidentally collapsed into one phase fraction");
}

void co2_dual_model() {
    const auto& reference = golden[0];
    const auto model = binary_model(*reference.gas);
    const Vec feed{reference.feed_gas, 1.0 - reference.feed_gas};
    const auto result = fl::solve_sw92_whitson_dual_model_observables(
        reference.pressure_pa, reference.temperature_k, feed,
        model, reference.molality);
    check_dual(result, reference);
}

void methane_dual_model() {
    const auto& reference = golden[1];
    const auto model = binary_model(*reference.gas);
    const Vec feed{reference.feed_gas, 1.0 - reference.feed_gas};
    const auto result = fl::solve_sw92_whitson_dual_model_observables(
        reference.pressure_pa, reference.temperature_k, feed,
        model, reference.molality);
    check_dual(result, reference);
}

void permutation() {
    const auto& reference = golden[0];
    const auto normal_model = binary_model(*reference.gas, false);
    const auto reverse_model = binary_model(*reference.gas, true);
    const auto first = fl::solve_sw92_whitson_dual_model_observables(
        reference.pressure_pa, reference.temperature_k,
        Vec{reference.feed_gas, 1.0 - reference.feed_gas}, normal_model, reference.molality);
    const auto second = fl::solve_sw92_whitson_dual_model_observables(
        reference.pressure_pa, reference.temperature_k,
        Vec{1.0 - reference.feed_gas, reference.feed_gas}, reverse_model, reference.molality);
    check_dual(first, reference, false);
    check_dual(second, reference, true);
    require(first.component_ids[0] == second.component_ids[1] &&
                first.component_ids[1] == second.component_ids[0],
            "component permutation metadata did not follow the snapshot");
}

void family_failure_isolated() {
    const auto& reference = golden[0];
    const auto model = binary_model(*reference.gas);
    fl::Sw92DualModelOptions options;
    options.aqueous.root_options.max_iterations = 1;
    const auto result = fl::solve_sw92_whitson_dual_model_observables(
        reference.pressure_pa, reference.temperature_k,
        Vec{reference.feed_gas, 1.0 - reference.feed_gas}, model, reference.molality, options);
    require(result.status == fl::Sw92DualModelStatus::aqueous_observable_unavailable,
            "AQ-only resource failure did not remain family-specific");
    require(result.aqueous_run.solution.status == fl::PtSplitStatus::indeterminate &&
                result.aqueous_run.solution.initial_stability.reference_issue ==
                    fl::StabilityPropertyIssue::root_iteration_limit,
            "AQ root-budget failure semantics changed");
    require(result.nonaqueous_run.solution.status ==
                fl::PtSplitStatus::two_phase_no_instability_found,
            "NA run was lost when AQ failed");
    require(!result.observables.has_value(),
            "cross-model observables were published with one unavailable family target");
}

void multicomponent_label_guard() {
    const auto model = ternary_structural_model();
    const auto result = fl::solve_sw92_whitson_dual_model_observables(
        3.0e6, 340.0, Vec{0.5, 0.2, 0.3}, model, 0.0);
    require(result.status == fl::Sw92DualModelStatus::multicomponent_phase_label_not_implemented,
            "multicomponent call silently acquired a physical water-cutoff label");
    require(!result.observables.has_value() &&
                result.aqueous_run.family == th::SwPhaseFamily::aqueous &&
                result.nonaqueous_run.family == th::SwPhaseFamily::nonaqueous &&
                result.aqueous_run.component_ids.size() == 3 &&
                result.nonaqueous_run.component_ids.size() == 3,
            "multicomponent family runs were not independently retained");
}

void headers() {
    require(sw92_dual_model_header(), "standalone SW92 dual-model public header");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test cases[] = {
    {"co2_dual_model", co2_dual_model},
    {"methane_dual_model", methane_dual_model},
    {"permutation", permutation},
    {"family_failure_isolated", family_failure_isolated},
    {"multicomponent_label_guard", multicomponent_label_guard},
    {"headers", headers}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument("exactly one test name required");
        }
        for (const auto& [name, run] : cases) {
            if (name == argv[1]) {
                run();
                std::cout << "[PASS] " << name << '\n';
                return 0;
            }
        }
        throw std::invalid_argument("unknown test");
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
