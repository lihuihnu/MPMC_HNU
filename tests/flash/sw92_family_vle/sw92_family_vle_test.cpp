#include <mpmc/flash/sw92_split.hpp>

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

bool sw92_split_header();

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

void near(double actual, long double expected, long double relative = 5e-9L,
          long double absolute = 5e-12L,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(static_cast<long double>(actual) - expected) >
            absolute + relative * std::abs(expected)) {
        require(false, "reference mismatch", where);
    }
}

template <typename Error, typename Function>
void expect_error(Function&& function) {
    bool caught = false;
    try {
        function();
    } catch (const Error&) {
        caught = true;
    }
    require(caught, "expected exception missing");
}

th::Provenance sw92_source(std::string locator) {
    return {th::SourceKind::literature,
            "doi:10.1016/0378-3812(92)85105-H",
            "Fluid Phase Equilibria 77 (1992) 217-240; authors' errata in supplied PDF",
            std::move(locator),
            "SW92 corrected-original one-family VLE regression; not experimental validation",
            "User PDF sha256:cb5b1d5034d78d934e887449ce0c89692d43431835d1371da5606d95b2c6bf58",
            "Bibliographic/formula facts only; source PDF is not redistributed"};
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
constexpr GasSpec water{"water", "Water", th::Sw92Species::water,
                        647.3, 221.2, 0.3434, 0.0};

th::Sw92Phase<double> binary_model(const GasSpec& gas, bool reverse = false) {
    const auto identity = sw92_source("Table 3 component identities");
    const auto properties = sw92_source("Table 3 physical properties");
    const auto na_pair = sw92_source("Table 5 non-aqueous water BIP");

    std::vector<th::Component> catalog{
        {gas.id, gas.display, th::ComponentKind::pure, identity, std::nullopt},
        {water.id, water.display, th::ComponentKind::pure, identity, std::nullopt}};

    th::Sw92ParameterInput input;
    input.model_id = std::string(th::sw92_corrected_profile);
    input.dataset_id = std::string("SW92-Table3-Table5-") + gas.id + "-family-vle";
    input.revision = "user-PDF-sha256-cb5b1d50";
    input.applicability = {{std::nullopt, std::nullopt,
                            sw92_source("Eqs. (9)-(17), Tables 2-5; no global T/p bounds inferred")},
                           std::nullopt};
    input.pure.push_back({gas.id, gas.species,
        sourced(gas.tc, th::Unit::kelvin, properties, "K", "identity"),
        sourced(gas.pc_bar * 100000.0, th::Unit::pascal, properties,
                "bar", "bar * 100000 -> Pa"),
        sourced(gas.omega, th::Unit::dimensionless, properties)});
    input.pure.push_back({water.id, water.species,
        sourced(water.tc, th::Unit::kelvin, properties, "K", "identity"),
        sourced(water.pc_bar * 100000.0, th::Unit::pascal, properties,
                "bar", "bar * 100000 -> Pa"),
        sourced(water.omega, th::Unit::dimensionless, properties)});
    input.water_binary.push_back({gas.id, th::Sw92NonAqueousWaterRule::constant,
                                  sourced(gas.nonaqueous_kij, th::Unit::dimensionless, na_pair)});

    const std::vector<std::string> order = reverse
        ? std::vector<std::string>{water.id, gas.id}
        : std::vector<std::string>{gas.id, water.id};
    return th::Sw92Phase<double>::from_parameters(
        th::Sw92ParameterSet::create(catalog, order, input));
}

struct VleGolden {
    const GasSpec* gas;
    th::SwPhaseFamily family;
    double pressure_pa;
    double temperature_k;
    double molality;
    double feed_gas;
    long double liquid_gas;
    long double vapor_gas;
    long double vapor_fraction;
    long double liquid_z;
    long double vapor_z;
};

constexpr VleGolden golden[] = {
    {&co2, th::SwPhaseFamily::nonaqueous, 3.0e6, 340.0, 0.0, 0.7,
     0.000284148400815315029223095498982335673842538744L,
     0.988725814519670096671192287725749616677997269L,
     0.707897972721687782431207265035175917208561212L,
     0.0231907380927525026099073080035L,
     0.886568486878578795928863543918L},
    {&methane, th::SwPhaseFamily::aqueous, 1.0e7, 350.0, 1.0, 0.5,
     0.000906872773851126638054949435435153986776308908L,
     0.990434503404630856446428419494771741626369134L,
     0.504375129886973712103194630205147652712593350L,
     0.0755670459371378192055074638222L,
     0.904713093689994850424138680252L}
};

void check_accepted(const fl::Sw92FamilyPtSplitResult& result,
                    const VleGolden& reference, bool reverse = false) {
    const auto& solution = result.solution;
    require(solution.status == fl::PtSplitStatus::two_phase_no_instability_found,
            "one-family VLE was not accepted");
    require(solution.candidate() != nullptr && solution.final_stability.has_value(),
            "accepted one-family VLE is missing candidate/final stability");
    require(solution.final_stability->status == fl::StabilityStatus::no_instability_found,
            "final same-family stability did not pass");
    require(!solution.global_stability_proven,
            "finite one-family search cannot prove global stability");

    const auto& candidate = *solution.candidate();
    const std::size_t gas_index = reverse ? 1U : 0U;
    const std::size_t water_index = reverse ? 0U : 1U;
    near(candidate.fractions.vapor_fraction, reference.vapor_fraction);
    near(candidate.fractions.liquid[gas_index], reference.liquid_gas);
    near(candidate.fractions.vapor[gas_index], reference.vapor_gas);
    near(candidate.fractions.liquid[water_index], 1.0L - reference.liquid_gas);
    near(candidate.fractions.vapor[water_index], 1.0L - reference.vapor_gas);
    near(candidate.liquid.z, reference.liquid_z);
    near(candidate.vapor.z, reference.vapor_z);
    require(candidate.liquid.z < candidate.vapor.z,
            "candidate density ordering changed");
    require(candidate.fugacity_norm <= solution.options.iteration.fugacity_tolerance,
            "accepted state exceeds fugacity tolerance");
    require(candidate.fractions.mass_absolute <= solution.options.iteration.mass_absolute_tolerance &&
                candidate.fractions.mass_relative <= solution.options.iteration.mass_relative_tolerance,
            "accepted state exceeds material-balance tolerance");

    require(result.family == reference.family &&
                result.nacl_molality_mol_per_kg_water == reference.molality &&
                result.model_profile == th::sw92_corrected_profile &&
                result.phase_convention == th::sw92_pt_convention &&
                result.component_ids.size() == 2,
            "one-family VLE metadata changed");
}

void provider_roles() {
    const auto model = binary_model(co2);
    fl::Sw92FamilyVleEvaluator evaluator(
        model, 0.0, th::SwPhaseFamily::nonaqueous);
    const Vec composition{0.7, 0.3};
    const auto liquid = evaluator(
        3.0e6, 340.0, composition, fl::PtPhaseRole::liquid_candidate);
    const auto vapor = evaluator(
        3.0e6, 340.0, composition, fl::PtPhaseRole::vapor_candidate);
    require(liquid.activity.branch == 0U && vapor.activity.branch == 2U,
            "requested same-family root roles changed");
    require(liquid.z < vapor.z,
            "same-family liquid/vapor candidate ordering changed");
    require(evaluator.family() == th::SwPhaseFamily::nonaqueous &&
                evaluator.nacl_molality_mol_per_kg_water() == 0.0,
            "provider lost fixed family/molality state");
}

void co2_nonaqueous_vle() {
    const auto& reference = golden[0];
    const auto model = binary_model(*reference.gas);
    fl::Sw92FamilyVleEvaluator evaluator(
        model, reference.molality, reference.family);
    const Vec feed{reference.feed_gas, 1.0 - reference.feed_gas};
    const auto result = fl::solve_sw92_pt_family_vle(
        reference.pressure_pa, reference.temperature_k, feed, evaluator);
    check_accepted(result, reference);
}

void methane_aqueous_vle() {
    const auto& reference = golden[1];
    const auto model = binary_model(*reference.gas);
    fl::Sw92FamilyVleEvaluator evaluator(
        model, reference.molality, reference.family);
    const Vec feed{reference.feed_gas, 1.0 - reference.feed_gas};
    const auto result = fl::solve_sw92_pt_family_vle(
        reference.pressure_pa, reference.temperature_k, feed, evaluator);
    check_accepted(result, reference);
}

void permutation() {
    const auto& reference = golden[0];
    const auto normal_model = binary_model(*reference.gas, false);
    const auto reverse_model = binary_model(*reference.gas, true);
    fl::Sw92FamilyVleEvaluator normal(
        normal_model, reference.molality, reference.family);
    fl::Sw92FamilyVleEvaluator reversed(
        reverse_model, reference.molality, reference.family);
    const auto first = fl::solve_sw92_pt_family_vle(
        reference.pressure_pa, reference.temperature_k,
        Vec{reference.feed_gas, 1.0 - reference.feed_gas}, normal);
    const auto second = fl::solve_sw92_pt_family_vle(
        reference.pressure_pa, reference.temperature_k,
        Vec{1.0 - reference.feed_gas, reference.feed_gas}, reversed);
    check_accepted(first, reference, false);
    check_accepted(second, reference, true);
    near(first.solution.candidate()->fractions.vapor_fraction,
         second.solution.candidate()->fractions.vapor_fraction, 1e-10L, 1e-12L);
}

void contracts_and_failures() {
    const auto model = binary_model(co2);
    fl::Sw92FamilyVleEvaluator evaluator(
        model, 0.0, th::SwPhaseFamily::nonaqueous);
    expect_error<std::invalid_argument>([&] {
        (void)fl::solve_sw92_pt_family_vle(
            3.0e6, 340.0, Vec{1.0}, evaluator);
    });
    expect_error<std::invalid_argument>([&] {
        (void)evaluator(3.0e6, 340.0, Vec{0.7, 0.3},
                        static_cast<fl::PtPhaseRole>(99));
    });

    fl::Sw92FamilyVleEvaluator limited(
        model, 0.0, th::SwPhaseFamily::nonaqueous, th::Sw92RootOptions{1});
    const auto failed = fl::solve_sw92_pt_family_vle(
        3.0e6, 340.0, Vec{0.7, 0.3}, limited);
    require(failed.solution.status == fl::PtSplitStatus::indeterminate &&
                failed.solution.initial_stability.reference_issue ==
                    fl::StabilityPropertyIssue::root_iteration_limit,
            "root-iteration failure was not preserved");

    expect_error<std::domain_error>([&] {
        fl::Sw92FamilyVleEvaluator bad(
            model, -1.0, th::SwPhaseFamily::aqueous);
    });
}

void headers() {
    require(sw92_split_header(), "standalone SW92 split public header");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test cases[] = {
    {"provider_roles", provider_roles},
    {"co2_nonaqueous_vle", co2_nonaqueous_vle},
    {"methane_aqueous_vle", methane_aqueous_vle},
    {"permutation", permutation},
    {"contracts_and_failures", contracts_and_failures},
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
