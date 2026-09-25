#include <mpmc/flow/fugacity_equilibrium_linearization.hpp>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool fugacity_equilibrium_linearization_header();

namespace {

namespace flow = mpmc::flow;

void require(
    bool condition,
    std::string_view message,
    std::source_location where =
        std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(
            std::string{where.file_name()} +
            ":" +
            std::to_string(where.line()) +
            ": " +
            std::string{message});
    }
}

void near(
    double actual,
    double expected,
    std::source_location where =
        std::source_location::current()) {
    if (!std::isfinite(actual) ||
        !std::isfinite(expected) ||
        std::abs(actual - expected) >
            1.0e-14) {
        std::cerr
            << "actual=" << actual
            << " expected=" << expected
            << '\n';
        require(false, "numeric mismatch", where);
    }
}

template <class Function>
void expect_invalid(
    Function&& function,
    std::string_view fragment) {
    try {
        function();
    } catch (const std::invalid_argument& error) {
        require(
            std::string_view{error.what()}.find(fragment) !=
                std::string_view::npos,
            "invalid_argument diagnostic changed");
        return;
    }
    throw std::runtime_error(
        "expected std::invalid_argument");
}

flow::NaturalVariableLayout3P layout() {
    return flow::NaturalVariableLayout3P{
        flow::NaturalVariableCompositionPivot3P::
            from_dependent_components(
                3U,
                {1U, 0U, 2U})};
}

std::vector<double> jacobian(
    std::size_t rows,
    std::size_t columns) {
    std::vector<double> values(
        rows * columns,
        0.0);
    for (std::size_t row = 0U;
         row < rows;
         ++row) {
        for (std::size_t column = 0U;
             column < columns;
             ++column) {
            values[
                row * columns +
                column] =
                0.001 *
                static_cast<double>(
                    1U +
                    row * 100U +
                    column);
        }
    }
    return values;
}

void carrier() {
    const auto chart = layout();
    const std::size_t q =
        chart.unknown_count();
    const std::vector<double>
        values{0.1, 0.2, 0.3, -0.4, -0.5, -0.6};
    const auto derivatives =
        jacobian(
            values.size(),
            q);

    const flow::
        FugacityEquilibriumResidualLinearization3P
        linearization{
            chart,
            {"A", "B", "C"},
            flow::FugacityEquilibriumResidual3P<double>{
                3U,
                values},
            q,
            derivatives};

    require(
        linearization.component_count() == 3U &&
            linearization.input_count() == q &&
            linearization.residual_count() == 6U &&
            linearization.values().size() == 6U &&
            linearization.jacobian().size() ==
                6U * q,
        "fugacity linearization shape changed");

    for (std::size_t phase = 1U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<flow::PhaseSlot3>(
                phase);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            const std::size_t local =
                (phase - 1U) * 3U +
                component;
            require(
                linearization.local_row_index(
                    slot,
                    component) ==
                    local &&
                linearization.equation_index(
                    slot,
                    component) ==
                    4U + local,
                "fugacity linearization row mapping changed");
            near(
                linearization.residual(
                    slot,
                    component),
                values[local]);

            for (std::size_t column = 0U;
                 column < q;
                 ++column) {
                near(
                    linearization.d_residual(
                        slot,
                        component,
                        column),
                    derivatives[
                        local * q +
                        column]);
            }
        }
    }
}

void invalid_inputs() {
    const auto chart = layout();
    const std::size_t q =
        chart.unknown_count();
    const std::vector<double>
        values{0.1, 0.2, 0.3, -0.4, -0.5, -0.6};
    const auto good =
        jacobian(
            values.size(),
            q);

    expect_invalid(
        [&] {
            (void)flow::
                FugacityEquilibriumResidualLinearization3P{
                    chart,
                    {"A", "B", "C"},
                    flow::FugacityEquilibriumResidual3P<double>{
                        3U,
                        values},
                    q - 1U,
                    good};
        },
        "shape/layout");

    expect_invalid(
        [&] {
            auto wrong = good;
            wrong.pop_back();
            (void)flow::
                FugacityEquilibriumResidualLinearization3P{
                    chart,
                    {"A", "B", "C"},
                    flow::FugacityEquilibriumResidual3P<double>{
                        3U,
                        values},
                    q,
                    std::move(wrong)};
        },
        "shape/layout");

    expect_invalid(
        [&] {
            auto nonfinite = good;
            nonfinite[0] =
                std::numeric_limits<double>::infinity();
            (void)flow::
                FugacityEquilibriumResidualLinearization3P{
                    chart,
                    {"A", "B", "C"},
                    flow::FugacityEquilibriumResidual3P<double>{
                        3U,
                        values},
                    q,
                    std::move(nonfinite)};
        },
        "non-finite");

    expect_invalid(
        [&] {
            auto bad_values = values;
            bad_values[2] =
                std::numeric_limits<double>::quiet_NaN();
            (void)flow::
                FugacityEquilibriumResidualLinearization3P{
                    chart,
                    {"A", "B", "C"},
                    flow::FugacityEquilibriumResidual3P<double>{
                        3U,
                        std::move(bad_values)},
                    q,
                    good};
        },
        "non-finite");

    expect_invalid(
        [&] {
            (void)flow::
                FugacityEquilibriumResidualLinearization3P{
                    chart,
                    {"A", "A", "C"},
                    flow::FugacityEquilibriumResidual3P<double>{
                        3U,
                        values},
                    q,
                    good};
        },
        "unique");
}

void headers() {
    require(
        fugacity_equilibrium_linearization_header(),
        "fugacity equilibrium linearization header probe failed");
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"carrier", carrier},
    {"invalid_inputs", invalid_inputs},
    {"headers", headers}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument(
                "one test name required");
        }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                std::cout
                    << "[PASS] "
                    << name
                    << '\n';
                return 0;
            }
        }
        throw std::invalid_argument(
            "unknown test");
    } catch (const std::exception& error) {
        std::cerr
            << "[FAIL] "
            << error.what()
            << '\n';
        return 1;
    }
}
