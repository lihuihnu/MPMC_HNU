#include <mpmc/flash/sw92_stability.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <source_location>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool stability_sw92_header();

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

void near(double actual, long double expected, double relative = 2e-11,
          double absolute = 2e-13,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(static_cast<long double>(actual) - expected) >
            absolute + relative * std::abs(expected)) {
        std::ostringstream out;
        out << std::setprecision(24) << "actual=" << actual << " expected=" << expected;
        require(false, out.str(), where);
    }
}

template <typename Error, typename Function>
void expect_error(Function&& fn) {
    bool caught = false;
    try {
        fn();
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
            "SW92 corrected-original fixed-family stability regression; not experimental validation",
            "User PDF sha256:cb5b1d5034d78d934e887449ce0c89692d43431835d1371da5606d95b2c6bf58",
            "Bibliographic/formula facts only; source PDF is not redistributed"};
}

th::SourcedScalar sourced(double value, th::Unit unit, const th::Provenance& source,
                          std::string original = "SI or dimensionless",
                          std::string conversion = "identity") {
    return {value, unit, source, std::move(original), std::move(conversion)};
}

th::Sw92Phase<double> co2_water_model(
    bool reverse = false,
    std::optional<th::ClosedInterval> molality_bounds = std::nullopt) {
    const auto identity = sw92_source("Table 3 component identities");
    const auto properties = sw92_source("Table 3 physical properties");
    const auto na_pair = sw92_source("Table 5 CO2/water non-aqueous BIP");

    std::vector<th::Component> catalog{
        {"carbon-dioxide", "Carbon dioxide", th::ComponentKind::pure, identity, std::nullopt},
        {"water", "Water", th::ComponentKind::pure, identity, std::nullopt}};

    th::Sw92ParameterInput input;
    input.model_id = std::string(th::sw92_corrected_profile);
    input.dataset_id = "SW92-Table3-Table5-corrected-stability";
    input.revision = "user-PDF-sha256-cb5b1d50";
    input.applicability = {{std::nullopt, std::nullopt,
                            sw92_source("Eqs. (9)-(17), Tables 2-5; no global T/p bounds inferred")},
                           molality_bounds};
    input.pure.push_back({"carbon-dioxide", th::Sw92Species::carbon_dioxide,
        sourced(304.2, th::Unit::kelvin, properties, "K", "identity"),
        sourced(73.8e5, th::Unit::pascal, properties, "bar", "bar * 100000 -> Pa"),
        sourced(0.2273, th::Unit::dimensionless, properties)});
    input.pure.push_back({"water", th::Sw92Species::water,
        sourced(647.3, th::Unit::kelvin, properties, "K", "identity"),
        sourced(221.2e5, th::Unit::pascal, properties, "bar", "bar * 100000 -> Pa"),
        sourced(0.3434, th::Unit::dimensionless, properties)});
    input.water_binary.push_back({"carbon-dioxide", th::Sw92NonAqueousWaterRule::constant,
                                  sourced(0.1896, th::Unit::dimensionless, na_pair)});

    const std::vector<std::string> order = reverse
        ? std::vector<std::string>{"water", "carbon-dioxide"}
        : std::vector<std::string>{"carbon-dioxide", "water"};
    const auto parameters = th::Sw92ParameterSet::create(catalog, order, input);
    return th::Sw92Phase<double>::from_parameters(parameters);
}

void family_root_selection() {
    const auto model = co2_water_model();
    fl::Sw92FamilyStabilityEvaluator evaluator(
        model, 0.0, th::SwPhaseFamily::nonaqueous);
    th::Sw92PhaseWorkspace<double> workspace;
    const Vec composition{0.7, 0.3};
    const auto roots = model.roots(3e6, 340.0, composition, 0.0,
                                   th::SwPhaseFamily::nonaqueous, workspace);
    require(roots.status == th::Sw92RootStatus::success && roots.count == 3,
            "CO2/water NA fixture must retain three algebraic roots");
    require(roots.roots[0].slope_sign > 0 && roots.roots[1].slope_sign < 0 &&
                roots.roots[2].slope_sign > 0,
            "mechanical root topology changed");

    const auto selected = evaluator(3e6, 340.0, composition);
    require(selected.branch == 2 && selected.smooth,
            "minimum-Gibbs nonaqueous root selection changed");
    for (std::size_t root_index : {std::size_t{0}, std::size_t{2}}) {
        const auto candidate = model.evaluate(3e6, 340.0, composition, 0.0,
            th::SwPhaseFamily::nonaqueous, root_index, workspace);
        double difference = 0.0;
        for (std::size_t i = 0; i < composition.size(); ++i) {
            difference += composition[i] * (selected.ln_phi[i] - candidate.ln_phi[i]);
        }
        require(difference <= 3e-13, "selected NA root is not the Gibbs minimum");
    }
}

void aqueous_route() {
    const auto model = co2_water_model();
    fl::Sw92FamilyStabilityEvaluator evaluator(
        model, 1.0, th::SwPhaseFamily::aqueous);
    const Vec composition{0.02, 0.98};
    const auto selected = evaluator(20e6, 423.15, composition);
    require(selected.branch == 0 && selected.smooth,
            "aqueous single-root route changed");
    require(evaluator.family() == th::SwPhaseFamily::aqueous &&
                evaluator.nacl_molality_mol_per_kg_water() == 1.0,
            "fixed family/molality state changed");
    th::Sw92PhaseWorkspace<double> workspace;
    const auto direct = model.evaluate(20e6, 423.15, composition, 1.0,
        th::SwPhaseFamily::aqueous, 0, workspace);
    for (std::size_t i = 0; i < composition.size(); ++i) {
        near(selected.ln_phi[i], direct.ln_phi[i], 2e-13, 2e-13);
    }
}

void generic_bridge() {
    const auto model = co2_water_model();
    fl::Sw92FamilyStabilityEvaluator evaluator(
        model, 0.0, th::SwPhaseFamily::nonaqueous);
    fl::StabilityOptions options;
    options.automatic_starts = false;
    const Vec feed{0.7, 0.3};
    const std::vector<Vec> starts{{0.5, 0.5}};
    const auto result = fl::test_sw92_pt_family_stability(
        3e6, 340.0, feed, evaluator, options, starts);

    require(result.search.status == fl::StabilityStatus::unstable,
            "fixed-family negative TPD witness was lost");
    require(result.search.reference && result.search.reference->branch == 2,
            "feed must use the same-family minimum-Gibbs root");
    require(result.search.trials.size() == 1 && result.search.trials[0].point &&
                result.search.trials[0].status == fl::StabilityTrialStatus::negative_tpd &&
                result.search.trials[0].point->branch == 0,
            "trial must select its own same-family minimum-Gibbs root");

    constexpr long double sw92_tpd_golden[] = {
        -0.1645820525617525951779453326142646500459L
    };
    near(result.search.trials[0].point->value, sw92_tpd_golden[0]);
    require(!result.search.global_stability_proven,
            "finite fixed-family search cannot prove global stability");
    require(result.dataset_id == "SW92-Table3-Table5-corrected-stability" &&
                result.component_ids == std::vector<std::string>{"carbon-dioxide", "water"} &&
                result.family == th::SwPhaseFamily::nonaqueous &&
                result.nacl_molality_mol_per_kg_water == 0.0 &&
                result.model_profile == th::sw92_corrected_profile &&
                result.phase_convention == th::sw92_pt_convention,
            "fixed-family metadata changed");
}

void contracts_and_failures() {
    const auto model = co2_water_model();
    expect_error<std::domain_error>([&] {
        fl::Sw92FamilyStabilityEvaluator bad(model, -1.0,
                                             th::SwPhaseFamily::aqueous);
    });
    expect_error<std::domain_error>([&] {
        fl::Sw92FamilyStabilityEvaluator bad(
            model, std::numeric_limits<double>::quiet_NaN(),
            th::SwPhaseFamily::aqueous);
    });
    expect_error<std::invalid_argument>([&] {
        fl::Sw92FamilyStabilityEvaluator bad(
            model, 0.0, static_cast<th::SwPhaseFamily>(99));
    });
    expect_error<std::invalid_argument>([&] {
        fl::Sw92FamilyStabilityEvaluator bad(
            model, 0.0, th::SwPhaseFamily::aqueous, th::Sw92RootOptions{0});
    });

    const auto bounded = co2_water_model(false, th::ClosedInterval{0.0, 2.0});
    expect_error<std::domain_error>([&] {
        fl::Sw92FamilyStabilityEvaluator bad(
            bounded, 3.0, th::SwPhaseFamily::aqueous);
    });

    fl::Sw92FamilyStabilityEvaluator evaluator(
        model, 0.0, th::SwPhaseFamily::nonaqueous);
    expect_error<std::invalid_argument>([&] {
        (void)fl::test_sw92_pt_family_stability(
            3e6, 340.0, Vec{1.0}, evaluator);
    });

    fl::Sw92FamilyStabilityEvaluator limited(
        model, 0.0, th::SwPhaseFamily::nonaqueous, th::Sw92RootOptions{1});
    const auto failed = fl::test_sw92_pt_family_stability(
        3e6, 340.0, Vec{0.7, 0.3}, limited);
    require(failed.search.status == fl::StabilityStatus::indeterminate &&
                failed.search.reference_issue ==
                    fl::StabilityPropertyIssue::root_iteration_limit,
            "SW92 root-iteration failure mapping changed");
}

void permutation() {
    const auto normal = co2_water_model(false);
    const auto reversed = co2_water_model(true);
    fl::Sw92FamilyStabilityEvaluator a(
        normal, 0.0, th::SwPhaseFamily::nonaqueous);
    fl::Sw92FamilyStabilityEvaluator b(
        reversed, 0.0, th::SwPhaseFamily::nonaqueous);
    const auto first = a(3e6, 340.0, Vec{0.7, 0.3});
    const auto second = b(3e6, 340.0, Vec{0.3, 0.7});
    require(first.branch == second.branch && first.smooth == second.smooth,
            "root diagnostic changed under component permutation");
    near(first.ln_phi[0], second.ln_phi[1], 2e-13, 2e-13);
    near(first.ln_phi[1], second.ln_phi[0], 2e-13, 2e-13);
}

void headers() {
    require(stability_sw92_header(), "standalone SW92 stability public header");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test cases[] = {
    {"family_root_selection", family_root_selection},
    {"aqueous_route", aqueous_route},
    {"generic_bridge", generic_bridge},
    {"contracts_and_failures", contracts_and_failures},
    {"permutation", permutation},
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
