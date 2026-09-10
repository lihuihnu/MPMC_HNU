#include <mpmc/flash/pr76_sensitivity.hpp>

#include "sensitivity_references.hpp"
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

void require(bool value, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!value) {
        throw std::runtime_error(
            std::string(where.file_name()) + ":" + std::to_string(where.line()) +
            ": " + std::string(message));
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

// Derivatives span Pa, K and dimensionless feed coordinates. This independent
// reference budget does not modify flash or sensitivity acceptance thresholds.
void near_derivative(double actual, long double expected) {
    const long double scale = std::max(1.0e-7L, std::abs(expected));
    const long double tolerance = 2.0e-8L * scale;
    require(
        std::isfinite(actual) &&
            std::abs(static_cast<long double>(actual) - expected) <= tolerance,
        "implicit sensitivity differs from independent Decimal reference");
}

th::SourcedScalar sourced_scalar(double value, th::Unit unit,
                                  const th::Provenance& source,
                                  std::string original_unit,
                                  std::string conversion) {
    return {value, unit, source, std::move(original_unit), std::move(conversion)};
}

th::Pr76Phase<double> binary_model() {
    const th::Provenance source{
        th::SourceKind::literature,
        "https://academicweb.nd.edu/~markst/zm97a.pdf",
        "revised October 1997",
        "Section 4.3 and Table 4",
        "PR numerical regression, not experimental validation",
        "Read author-hosted PDF pages",
        "Limited attributed factual parameters; no article text or code reproduced"
    };
    const std::array<std::string, 2> ids{"nitrogen", "ethane"};
    const std::array<std::array<double, 3>, 2> specs{{
        {126.2, 3.39e6, 0.04},
        {305.4, 4.88e6, 0.098}
    }};
    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "Hua1997-nitrogen-ethane";
    input.revision = "split-regression-v1";
    input.applicability = {std::nullopt, std::nullopt, source};
    for (std::size_t i = 0; i < ids.size(); ++i) {
        catalog.push_back({ids[i], ids[i], th::ComponentKind::pure, source, {}});
        input.pure.push_back({
            ids[i],
            sourced_scalar(specs[i][0], th::Unit::kelvin, source, "K", "identity"),
            sourced_scalar(specs[i][1], th::Unit::pascal, source, "bar",
                           "bar * 100000 -> Pa"),
            sourced_scalar(specs[i][2], th::Unit::dimensionless, source,
                           "dimensionless", "identity")
        });
    }
    input.binary.push_back({
        ids[0], ids[1],
        sourced_scalar(0.08, th::Unit::dimensionless, source,
                       "dimensionless", "identity")
    });
    const std::vector<std::string> order{ids.begin(), ids.end()};
    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(catalog, order, input));
}

th::Pr76Phase<double> ternary_model() {
    const th::Provenance pure{
        th::SourceKind::literature,
        "https://doi.org/10.1002/aic.16730",
        "AIChE Journal 65(11), e16730; first published 26 July 2019",
        "Table 1, methane/ethane/propane Peng-Robinson parameters",
        "Tc/Pc/omega only; fitted methane-propane binary kij is not used here",
        "Read open Wiley full text and Table 1",
        "Limited attributed factual parameters; no article text reproduced"
    };
    const th::Provenance zero_kij{
        th::SourceKind::literature,
        "https://doi.org/10.1021/i160057a011",
        "Ind. Eng. Chem. Fundam. 15(1), 1976",
        "journal p.62, methane-ethane-propane ternary discussion and Figure 4",
        "Original ternary example used no interaction coefficients; explicit zero kij",
        "Read project-supplied Peng-Robinson 1976 PDF",
        "Limited attributed factual model statement; no article text reproduced"
    };
    constexpr std::array<std::string_view, 3> ids{
        "methane", "ethane", "propane"
    };
    constexpr std::array<std::array<double, 3>, 3> specs{{
        {190.555, 4.595e6, 0.0},
        {305.4, 4.88e6, 0.099},
        {369.825, 4.248e6, 0.15308}
    }};
    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "DeitersBell-aic16730-PengRobinson1976-ternary";
    input.revision = "ternary-regression-v1";
    input.applicability = {std::nullopt, std::nullopt, zero_kij};
    for (std::size_t i = 0; i < ids.size(); ++i) {
        catalog.push_back({
            std::string(ids[i]), std::string(ids[i]), th::ComponentKind::pure,
            pure, {}
        });
        input.pure.push_back({
            std::string(ids[i]),
            sourced_scalar(specs[i][0], th::Unit::kelvin, pure, "K", "identity"),
            sourced_scalar(specs[i][1], th::Unit::pascal, pure, "MPa",
                           "MPa * 1e6 -> Pa"),
            sourced_scalar(specs[i][2], th::Unit::dimensionless, pure,
                           "dimensionless", "identity")
        });
    }
    const auto zero = sourced_scalar(
        0.0, th::Unit::dimensionless, zero_kij, "dimensionless", "identity");
    input.binary.push_back({"methane", "ethane", zero});
    input.binary.push_back({"methane", "propane", zero});
    input.binary.push_back({"ethane", "propane", zero});
    const std::vector<std::string> order{"methane", "ethane", "propane"};
    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(catalog, order, input));
}

void check_common(const fl::Pr76PtSplitResult& split,
                  const fl::Pr76PtVleSensitivityResult& result,
                  std::size_t expected_components) {
    const auto& sensitivity = result.sensitivity;
    require(split.solution.status == fl::PtSplitStatus::two_phase_no_instability_found,
            "reference flash state is no longer accepted two phase");
    require(sensitivity.status == fl::PtSensitivityStatus::success,
            sensitivity.diagnostic);
    require(sensitivity.component_count == expected_components,
            "sensitivity component count");
    require(sensitivity.input_count == expected_components + 1,
            "sensitivity reduced-input count");
    require(sensitivity.pressure_pa == split.solution.initial_stability.pressure_pa &&
                sensitivity.temperature_k == split.solution.initial_stability.temperature_k &&
                sensitivity.feed == split.solution.initial_stability.feed,
            "sensitivity base point metadata");
    require(sensitivity.pressure_column() == 0 &&
                sensitivity.temperature_column() == 1 &&
                sensitivity.dependent_feed_component() == expected_components - 1,
            "sensitivity coordinate helpers");
    require(sensitivity.log_k_jacobian.size() ==
                expected_components * sensitivity.input_count &&
                sensitivity.liquid_jacobian.size() ==
                expected_components * sensitivity.input_count &&
                sensitivity.vapor_jacobian.size() ==
                expected_components * sensitivity.input_count &&
                sensitivity.vapor_fraction_gradient.size() == sensitivity.input_count,
            "sensitivity matrix shape");
    require(std::isfinite(sensitivity.equilibrium_jacobian_rcond) &&
                sensitivity.equilibrium_jacobian_rcond > 0.1,
            "reference equilibrium Jacobian unexpectedly ill-conditioned");
    require(std::isfinite(sensitivity.linear_solve_backward_error) &&
                sensitivity.linear_solve_backward_error < 1e-12,
            "implicit solve backward error");
    require(sensitivity.equilibrium_residual_norm <=
                split.solution.options.iteration.fugacity_tolerance + 1e-11,
            "sensitivity local equation residual");
}

template <std::size_t N, std::size_t Q>
void compare_matrix(const fl::PtVleSensitivityResult& result,
                    const std::array<std::array<long double, Q>, N>& expected,
                    bool liquid, bool vapor, bool log_k) {
    static_assert(Q == N + 1);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t column = 0; column < Q; ++column) {
            double actual = 0.0;
            if (liquid) {
                actual = result.d_liquid(i, column);
            } else if (vapor) {
                actual = result.d_vapor(i, column);
            } else if (log_k) {
                actual = result.d_log_k(i, column);
            } else {
                throw std::logic_error("invalid test matrix selector");
            }
            near_derivative(actual, expected[i][column]);
        }
    }
}

void binary_reference() {
    const auto model = binary_model();
    fl::Pr76VleEvaluator evaluator(model);
    const Vec feed{0.30, 0.70};
    const auto split = fl::solve_pr76_pt_vle(7.6e6, 270.0, feed, evaluator);
    const auto result = fl::differentiate_pr76_pt_vle(split, evaluator);
    check_common(split, result, 2);
    const auto& sensitivity = result.sensitivity;
    require(sensitivity.feed_column(0) == 2, "binary feed coordinate");
    for (std::size_t column = 0; column < 3; ++column) {
        near_derivative(sensitivity.d_vapor_fraction(column),
                        sensitivity_reference::binary_vapor_fraction[column]);
    }
    compare_matrix(sensitivity, sensitivity_reference::binary_log_k,
                   false, false, true);
    compare_matrix(sensitivity, sensitivity_reference::binary_liquid,
                   true, false, false);
    compare_matrix(sensitivity, sensitivity_reference::binary_vapor,
                   false, true, false);
}

void ternary_reference() {
    const auto model = ternary_model();
    fl::Pr76VleEvaluator evaluator(model);
    const auto& anchor = ternary_reference::ternary_references[1];
    Vec feed;
    for (const long double value : anchor.feed) {
        feed.push_back(static_cast<double>(value));
    }
    const auto split = fl::solve_pr76_pt_vle(4.0e6, 220.0, feed, evaluator);
    const auto result = fl::differentiate_pr76_pt_vle(split, evaluator);
    check_common(split, result, 3);
    const auto& sensitivity = result.sensitivity;
    require(sensitivity.feed_column(0) == 2 && sensitivity.feed_column(1) == 3,
            "ternary feed coordinates");
    for (std::size_t column = 0; column < 4; ++column) {
        near_derivative(sensitivity.d_vapor_fraction(column),
                        sensitivity_reference::ternary_vapor_fraction[column]);
    }
    compare_matrix(sensitivity, sensitivity_reference::ternary_log_k,
                   false, false, true);
    compare_matrix(sensitivity, sensitivity_reference::ternary_liquid,
                   true, false, false);
    compare_matrix(sensitivity, sensitivity_reference::ternary_vapor,
                   false, true, false);
}

void guards() {
    const auto model = binary_model();
    fl::Pr76VleEvaluator evaluator(model);
    const Vec feed{0.30, 0.70};
    const auto split = fl::solve_pr76_pt_vle(7.6e6, 270.0, feed, evaluator);
    require(split.solution.status == fl::PtSplitStatus::two_phase_no_instability_found,
            "guard fixture no longer accepted");

    auto unaccepted = split;
    unaccepted.solution.status = fl::PtSplitStatus::indeterminate;
    require(fl::differentiate_pr76_pt_vle(unaccepted, evaluator).sensitivity.status ==
                fl::PtSensitivityStatus::solution_not_accepted,
            "unaccepted phase set produced a derivative");

    fl::PtSensitivityOptions boundary;
    boundary.minimum_derivative_phase_fraction = 0.45;
    require(fl::differentiate_pr76_pt_vle(split, evaluator, boundary).sensitivity.status ==
                fl::PtSensitivityStatus::phase_boundary,
            "derivative phase-boundary guard ignored");

    fl::PtSensitivityOptions conditioning;
    conditioning.minimum_reciprocal_condition = 0.90;
    require(fl::differentiate_pr76_pt_vle(split, evaluator, conditioning).sensitivity.status ==
                fl::PtSensitivityStatus::ill_conditioned_equilibrium,
            "caller conditioning requirement ignored");

    auto zero_support = split;
    zero_support.solution.initial_stability.feed = {0.0, 1.0};
    require(fl::differentiate_pr76_pt_vle(zero_support, evaluator).sensitivity.status ==
                fl::PtSensitivityStatus::unsupported_feed_support,
            "zero-feed support silently differentiated");

    auto wrong_model = split;
    wrong_model.component_ids[0] = "not-nitrogen";
    expect_error<std::invalid_argument>([&] {
        (void)fl::differentiate_pr76_pt_vle(wrong_model, evaluator);
    });

    fl::PtSensitivityOptions invalid;
    invalid.minimum_derivative_phase_fraction = 0.0;
    expect_error<std::invalid_argument>([&] {
        (void)fl::differentiate_pr76_pt_vle(split, evaluator, invalid);
    });
}

void invariants() {
    const auto model = ternary_model();
    fl::Pr76VleEvaluator evaluator(model);
    const auto& anchor = ternary_reference::ternary_references[1];
    Vec feed;
    for (const long double value : anchor.feed) {
        feed.push_back(static_cast<double>(value));
    }
    const auto split = fl::solve_pr76_pt_vle(4.0e6, 220.0, feed, evaluator);
    const auto result = fl::differentiate_pr76_pt_vle(split, evaluator);
    check_common(split, result, 3);
    const auto& sensitivity = result.sensitivity;
    const auto& point = *split.solution.candidate();

    for (std::size_t column = 0; column < sensitivity.input_count; ++column) {
        double dx_sum = 0.0;
        double dy_sum = 0.0;
        for (std::size_t i = 0; i < 3; ++i) {
            dx_sum += sensitivity.d_liquid(i, column);
            dy_sum += sensitivity.d_vapor(i, column);
        }
        require(std::abs(dx_sum) < 1e-11 && std::abs(dy_sum) < 1e-11,
                "differentiated phase normalization");

        for (std::size_t i = 0; i < 3; ++i) {
            const double dz =
                column < 2 ? 0.0
                           : (i == column - 2 ? 1.0 : (i == 2 ? -1.0 : 0.0));
            const double derivative =
                (1.0 - point.fractions.vapor_fraction) *
                    sensitivity.d_liquid(i, column) +
                point.fractions.vapor_fraction *
                    sensitivity.d_vapor(i, column) +
                (point.fractions.vapor[i] - point.fractions.liquid[i]) *
                    sensitivity.d_vapor_fraction(column) -
                dz;
            require(std::abs(derivative) < 2e-10,
                    "differentiated material balance");
        }
    }

    expect_error<std::out_of_range>([&] {
        (void)sensitivity.feed_column(2);
    });
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"binary_reference", binary_reference},
    {"ternary_reference", ternary_reference},
    {"guards", guards},
    {"invariants", invariants}
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
