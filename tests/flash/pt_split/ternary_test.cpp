#include <mpmc/flash/pr76_split.hpp>

#include "ternary_references.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
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

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

constexpr double pressure_pa = 4.0e6;
constexpr double temperature_k = 220.0;
constexpr double reference_tolerance = 2e-9;
constexpr std::string_view dataset_id =
    "DeitersBell-aic16730-PengRobinson1976-ternary";

void require(bool value, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!value) {
        throw std::runtime_error(
            std::string(where.file_name()) + ":" + std::to_string(where.line()) +
            ": " + std::string(message));
    }
}

void near(double actual, long double expected,
          double tolerance = reference_tolerance) {
    require(
        std::isfinite(actual) &&
            std::abs(static_cast<long double>(actual) - expected) <=
                static_cast<long double>(tolerance) * (1 + std::abs(expected)),
        "ternary independent-reference mismatch");
}

struct TernaryFixture {
    std::vector<th::Component> catalog;
    th::PrParameterInput input;
};

th::Provenance pure_source() {
    return {
        th::SourceKind::literature,
        "https://doi.org/10.1002/aic.16730",
        "AIChE Journal 65(11), e16730; first published 26 July 2019",
        "Table 1, methane/ethane/propane Peng-Robinson parameters",
        "Tc/Pc/omega only; this regression does not use that paper's fitted "
        "methane-propane binary interaction parameter",
        "Read open Wiley full text and Table 1",
        "Limited attributed factual parameters; no article text reproduced"
    };
}

th::Provenance zero_kij_source() {
    return {
        th::SourceKind::literature,
        "https://doi.org/10.1021/i160057a011",
        "Ind. Eng. Chem. Fundam. 15(1), 1976",
        "journal p.62, methane-ethane-propane ternary discussion and Figure 4",
        "Original ternary example states that no interaction coefficients were "
        "used; all three off-diagonal kij are therefore explicit zero records",
        "Read project-supplied Peng-Robinson 1976 PDF",
        "Limited attributed factual model statement; no article text reproduced"
    };
}

th::SourcedScalar scalar(
    double value, th::Unit unit, const th::Provenance& source,
    std::string original_unit, std::string conversion) {
    return {value, unit, source, std::move(original_unit), std::move(conversion)};
}

TernaryFixture make_fixture() {
    const auto pure = pure_source();
    const auto pairs = zero_kij_source();
    constexpr std::array<std::string_view, 3> ids{
        "methane", "ethane", "propane"
    };
    constexpr std::array<std::array<double, 3>, 3> specs{{
        {190.555, 4.595e6, 0.0},
        {305.4, 4.88e6, 0.099},
        {369.825, 4.248e6, 0.15308}
    }};

    TernaryFixture fixture;
    fixture.input.model_id = std::string(th::pr76_profile);
    fixture.input.dataset_id = std::string(dataset_id);
    fixture.input.revision = "ternary-regression-v1";
    fixture.input.applicability = {std::nullopt, std::nullopt, pairs};

    for (std::size_t i = 0; i < ids.size(); ++i) {
        fixture.catalog.push_back({
            std::string(ids[i]), std::string(ids[i]), th::ComponentKind::pure,
            pure, {}
        });
        fixture.input.pure.push_back({
            std::string(ids[i]),
            scalar(specs[i][0], th::Unit::kelvin, pure, "K", "identity"),
            scalar(
                specs[i][1], th::Unit::pascal, pure, "MPa",
                "MPa * 1e6 -> Pa"),
            scalar(
                specs[i][2], th::Unit::dimensionless, pure, "dimensionless",
                "identity")
        });
    }

    const auto zero = scalar(
        0.0, th::Unit::dimensionless, pairs, "dimensionless", "identity");
    fixture.input.binary.push_back({"methane", "ethane", zero});
    fixture.input.binary.push_back({"methane", "propane", zero});
    fixture.input.binary.push_back({"ethane", "propane", zero});
    return fixture;
}

th::Pr76Phase<double> make_model(std::span<const std::string> order) {
    const auto fixture = make_fixture();
    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(fixture.catalog, order, fixture.input));
}

void check_default_thresholds(const fl::PtSplitResult& solution) {
    require(
        solution.options.iteration.fugacity_tolerance == 1e-11 &&
            solution.options.iteration.mass_absolute_tolerance == 1e-12 &&
            solution.options.iteration.mass_relative_tolerance == 1e-10 &&
            solution.options.iteration.minimum_phase_fraction == 1e-10,
        "existing production acceptance thresholds changed");
}

void check_reference(
    const fl::Pr76PtSplitResult& result, const th::Pr76Phase<double>& model,
    const ternary_reference::TernaryReference& reference,
    std::span<const std::size_t> canonical_index) {
    const auto& solution = result.solution;
    check_default_thresholds(solution);
    require(
        solution.status == fl::PtSplitStatus::two_phase_no_instability_found,
        "ternary PT flash did not finish with the existing two-phase success state");
    require(
        solution.equations_converged() && solution.candidate(),
        "ternary PT flash lost its converged candidate");
    require(
        solution.initial_stability.status == fl::StabilityStatus::unstable,
        "ternary feed did not preserve initial instability evidence");
    require(
        solution.final_stability &&
            solution.final_stability->status ==
                fl::StabilityStatus::no_instability_found,
        "ternary final common-tangent review did not complete");
    require(!solution.global_stability_proven, "finite TPD search became a global proof");

    const auto& candidate = *solution.candidate();
    require(
        candidate.fugacity_norm <= solution.options.iteration.fugacity_tolerance,
        "stored ternary fugacity residual exceeds the unchanged threshold");
    require(
        candidate.fractions.mass_absolute <=
                solution.options.iteration.mass_absolute_tolerance &&
            candidate.fractions.mass_relative <=
                solution.options.iteration.mass_relative_tolerance,
        "stored ternary material balance exceeds the unchanged threshold");
    require(
        candidate.fractions.vapor_fraction >
                solution.options.iteration.minimum_phase_fraction &&
            1 - candidate.fractions.vapor_fraction >
                solution.options.iteration.minimum_phase_fraction,
        "reference state unexpectedly exercises the disappearance gate");

    for (std::size_t i = 0; i < canonical_index.size(); ++i) {
        const std::size_t canonical = canonical_index[i];
        near(candidate.fractions.liquid[i], reference.liquid[canonical]);
        near(candidate.fractions.vapor[i], reference.vapor[canonical]);
    }
    near(candidate.fractions.vapor_fraction, reference.vapor_fraction);
    near(candidate.liquid.z, reference.liquid_z);
    near(candidate.vapor.z, reference.vapor_z);
    require(solution.gibbs_change < 0, "accepted ternary split does not lower Gibbs energy");

    th::Pr76PhaseWorkspace<double> workspace;
    const auto liquid = model.evaluate_full(
        pressure_pa, temperature_k, candidate.fractions.liquid,
        candidate.liquid.activity.branch, workspace);
    const auto vapor = model.evaluate_full(
        pressure_pa, temperature_k, candidate.fractions.vapor,
        candidate.vapor.activity.branch, workspace);
    near(liquid.z, reference.liquid_z);
    near(vapor.z, reference.vapor_z);

    double recomputed_fugacity = 0.0;
    for (std::size_t i = 0; i < canonical_index.size(); ++i) {
        const double residual =
            std::log(candidate.fractions.liquid[i]) + liquid.ln_phi[i] -
            std::log(candidate.fractions.vapor[i]) - vapor.ln_phi[i];
        recomputed_fugacity = std::max(recomputed_fugacity, std::abs(residual));
    }
    require(
        recomputed_fugacity <=
            solution.options.iteration.fugacity_tolerance +
                32 * std::numeric_limits<double>::epsilon(),
        "fresh PR76 fugacity recomputation exceeds the production threshold");
}

void parameters_runtime() {
    const auto fixture = make_fixture();
    const std::vector<std::string> canonical{"methane", "ethane", "propane"};
    const auto parameters =
        th::PrParameterSet::create(fixture.catalog, canonical, fixture.input);
    require(parameters.components().size() == 3, "ternary snapshot size");
    require(parameters.binary_records().size() == 3, "explicit ternary pair count");
    for (const auto& pair : parameters.binary_records()) {
        require(pair.kij && pair.kij->value == 0.0, "zero kij is not explicit");
        require(
            pair.kij->source.reference == "https://doi.org/10.1021/i160057a011",
            "zero kij provenance lost");
    }

    auto incomplete = fixture.input;
    incomplete.binary.pop_back();
    bool caught = false;
    try {
        (void)th::PrParameterSet::create(fixture.catalog, canonical, incomplete);
    } catch (const th::ContractError& error) {
        caught = error.code() == th::ContractErrorCode::missing_parameter;
    }
    require(caught, "missing ternary pair was silently interpreted as zero");

    const std::array<std::vector<std::string>, 4> orders{{
        {"methane"},
        {"methane", "ethane", "propane"},
        {"propane", "methane"},
        {"ethane", "propane", "methane"}
    }};
    th::Pr76PhaseWorkspace<double> workspace;
    for (const auto& order : orders) {
        const auto model = make_model(order);
        const Vec composition(
            order.size(), 1.0 / static_cast<double>(order.size()));
        const auto roots =
            model.roots_full(pressure_pa, temperature_k, composition, workspace);
        require(
            roots.status == th::Pr76RootStatus::success && roots.count > 0,
            "runtime component-shape PR76 evaluation failed");
    }
}

void flash_reference() {
    const std::vector<std::string> order{"methane", "ethane", "propane"};
    const auto model = make_model(order);
    fl::Pr76VleEvaluator evaluator(model);
    const std::array<std::size_t, 3> identity{0, 1, 2};

    for (const auto& reference : ternary_reference::ternary_references) {
        Vec feed;
        feed.reserve(3);
        for (long double value : reference.feed) {
            feed.push_back(static_cast<double>(value));
        }
        const auto result = fl::solve_pr76_pt_vle(
            pressure_pa, temperature_k, feed, evaluator);
        check_reference(result, model, reference, identity);
        require(
            std::string_view(result.dataset_id) == dataset_id &&
                result.revision == "ternary-regression-v1" &&
                result.component_ids == order,
            "ternary result metadata does not match the ordered snapshot");
        std::cout
            << "beta=" << result.solution.candidate()->fractions.vapor_fraction
            << " fugacity=" << result.solution.candidate()->fugacity_norm
            << " evaluations=" << result.solution.split_evaluations << '\n';
    }
}

void permutations() {
    constexpr std::array<std::array<std::size_t, 3>, 6> permutations{{
        {0, 1, 2}, {0, 2, 1}, {1, 0, 2},
        {1, 2, 0}, {2, 0, 1}, {2, 1, 0}
    }};
    constexpr std::array<std::string_view, 3> ids{
        "methane", "ethane", "propane"
    };
    const auto& reference = ternary_reference::ternary_references[1];

    for (const auto& permutation : permutations) {
        std::vector<std::string> order;
        Vec feed;
        order.reserve(3);
        feed.reserve(3);
        for (const std::size_t canonical : permutation) {
            order.emplace_back(ids[canonical]);
            feed.push_back(static_cast<double>(reference.feed[canonical]));
        }
        const auto model = make_model(order);
        fl::Pr76VleEvaluator evaluator(model);
        const auto result = fl::solve_pr76_pt_vle(
            pressure_pa, temperature_k, feed, evaluator);
        check_reference(result, model, reference, permutation);
        require(result.component_ids == order, "permuted component IDs lost");
    }
}

void final_indeterminate() {
    const std::vector<std::string> order{"methane", "ethane", "propane"};
    const auto model = make_model(order);
    fl::Pr76VleEvaluator evaluator(model);
    const auto& reference = ternary_reference::ternary_references[1];
    Vec feed;
    feed.reserve(3);
    for (long double value : reference.feed) {
        feed.push_back(static_cast<double>(value));
    }
    fl::PtSplitOptions options;
    options.final_stability.max_evaluations = 1;

    const auto result = fl::solve_pr76_pt_vle(
        pressure_pa, temperature_k, feed, evaluator, options);
    const auto& solution = result.solution;
    check_default_thresholds(solution);
    require(
        solution.status == fl::PtSplitStatus::indeterminate,
        "final TPD budget exhaustion did not stay indeterminate");
    require(
        solution.equations_converged() && solution.candidate(),
        "final TPD indeterminate state discarded a converged candidate");
    require(
        solution.final_stability &&
            solution.final_stability->status == fl::StabilityStatus::indeterminate,
        "final TPD budget exhaustion did not remain explicitly unresolved");
    near(
        solution.candidate()->fractions.vapor_fraction,
        reference.vapor_fraction);
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"parameters_runtime", parameters_runtime},
    {"flash_reference", flash_reference},
    {"permutations", permutations},
    {"final_indeterminate", final_indeterminate}
};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument("one test name required");
        }
        for (const auto& [name, run] : tests) {
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
