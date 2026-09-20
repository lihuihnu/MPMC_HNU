#include <mpmc/flow/component_accumulation_time.hpp>

#include <algorithm>
#include <array>
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

bool component_accumulation_time_header();

namespace {

namespace fl = mpmc::flow;

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
    double relative = 2.0e-13,
    double absolute = 2.0e-13,
    std::source_location where =
        std::source_location::current()) {
    if (!std::isfinite(actual) ||
        !std::isfinite(expected) ||
        std::abs(actual - expected) >
            absolute +
                relative *
                    std::max(
                        std::abs(actual),
                        std::abs(expected))) {
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

fl::PoreVolumeComponentAccumulationSnapshot3P
snapshot(
    double porosity,
    std::vector<double> components) {
    double total = 0.0;
    for (double value : components) {
        total += value;
    }
    return {
        porosity,
        {"A", "B", "C"},
        std::move(components),
        total};
}

fl::PoreVolumeComponentAccumulationLinearization3P
linearization(double porosity) {
    const auto pivot =
        fl::NaturalVariableCompositionPivot3P::
            from_dependent_components(
                3U,
                {1U, 0U, 2U});
    const fl::NaturalVariableLayout3P layout{
        pivot};
    const std::size_t q =
        layout.unknown_count();

    std::vector<double> jacobian(
        3U * q,
        0.0);
    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        for (std::size_t column = 0U;
             column < q;
             ++column) {
            jacobian[
                component * q + column] =
                static_cast<double>(
                    (component + 1U) *
                    (column + 2U)) *
                0.125;
        }
    }

    std::vector<double> total_gradient(
        q,
        0.0);
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            total_gradient[column] +=
                jacobian[
                    component * q + column];
        }
    }

    return {
        layout,
        porosity,
        {"A", "B", "C"},
        q,
        std::move(jacobian),
        std::move(total_gradient)};
}

void backward_euler() {
    constexpr double porosity = 0.25;
    constexpr double dt = 10.0;

    const auto pair =
        fl::make_pore_volume_component_accumulation_pair(
            snapshot(
                porosity,
                {100.0, 200.0, 300.0}),
            snapshot(
                porosity,
                {90.0, 220.0, 270.0}));
    const auto current_linearization =
        linearization(porosity);

    const auto result =
        fl::
            build_backward_euler_component_accumulation_residual(
                pair,
                current_linearization,
                dt);

    require(
        result.component_ids ==
                std::vector<std::string>{"A", "B", "C"} &&
            result.component_count() == 3U &&
            result.input_count ==
                current_linearization.input_count &&
            result.current_layout
                    .composition_pivot()
                    .dependent_components() ==
                current_linearization.layout
                    .composition_pivot()
                    .dependent_components(),
        "backward-Euler identity/chart changed");
    near(result.porosity, porosity);
    near(result.time_step_seconds, dt);

    near(result.residual(0U), 1.0);
    near(result.residual(1U), -2.0);
    near(result.residual(2U), 3.0);
    near(
        result.total_residual_mol_per_bulk_m3_s,
        2.0);
    require(
        result.component_row_identity(0U) == "A" &&
            result.component_row_identity(1U) == "B" &&
            result.component_row_identity(2U) == "C",
        "backward-Euler component row identity changed");

    for (std::size_t column = 0U;
         column < result.input_count;
         ++column) {
        double sum = 0.0;
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            const double expected =
                current_linearization
                    .d_component(
                        component,
                        column) /
                dt;
            near(
                result.d_residual(
                    component,
                    column),
                expected);
            sum += expected;
        }
        near(
            result.d_total(column),
            current_linearization
                    .d_total(column) /
                dt);
        near(
            sum,
            result.d_total(column));
    }
}

void previous_is_frozen() {
    constexpr double porosity = 0.25;
    constexpr double dt = 5.0;
    const auto current =
        snapshot(
            porosity,
            {120.0, 180.0, 300.0});
    const auto current_linearization =
        linearization(porosity);

    const auto first =
        fl::
            build_backward_euler_component_accumulation_residual(
                fl::
                    make_pore_volume_component_accumulation_pair(
                        current,
                        snapshot(
                            porosity,
                            {100.0, 170.0, 280.0})),
                current_linearization,
                dt);
    const auto second =
        fl::
            build_backward_euler_component_accumulation_residual(
                fl::
                    make_pore_volume_component_accumulation_pair(
                        current,
                        snapshot(
                            porosity,
                            {110.0, 160.0, 290.0})),
                current_linearization,
                dt);

    require(
        first.component_residual_mol_per_bulk_m3_s !=
            second.component_residual_mol_per_bulk_m3_s,
        "changing previous snapshot did not change residual");
    require(
        first.component_jacobian ==
                second.component_jacobian &&
            first.total_residual_gradient ==
                second.total_residual_gradient,
        "previous snapshot incorrectly entered current-state Jacobian");
}


void near_steady_roundoff() {
    constexpr double porosity = 0.25;

    auto current =
        snapshot(
            porosity,
            {1.0e12, 2.0e12, 3.0e12});
    auto previous =
        snapshot(
            porosity,
            {1.0e12 - 1.0,
             2.0e12 + 2.0,
             3.0e12 - 1.0});

    // Each snapshot remains within the upstream machine-roundoff closure
    // contract at this inventory scale, even though differencing magnifies the
    // relative effect of the stored total's last-bit-scale discrepancy.
    current.total_accumulation_mol_per_bulk_m3 +=
        0.25;
    previous.total_accumulation_mol_per_bulk_m3 -=
        0.25;

    const auto pair =
        fl::make_pore_volume_component_accumulation_pair(
            current,
            previous);
    const auto result =
        fl::
            build_backward_euler_component_accumulation_residual(
                pair,
                linearization(porosity),
                1.0);

    require(
        std::isfinite(
            result.total_residual_mol_per_bulk_m3_s),
        "near-steady backward-Euler residual became non-finite");
    near(
        result.residual(0U) +
            result.residual(1U) +
            result.residual(2U),
        0.0,
        0.0,
        0.0);
}

void invalid_inputs() {
    constexpr double porosity = 0.25;
    const auto pair =
        fl::make_pore_volume_component_accumulation_pair(
            snapshot(
                porosity,
                {100.0, 200.0, 300.0}),
            snapshot(
                porosity,
                {90.0, 220.0, 270.0}));
    const auto valid =
        linearization(porosity);

    expect_invalid(
        [&] {
            (void)fl::
                build_backward_euler_component_accumulation_residual(
                    pair,
                    valid,
                    0.0);
        },
        "time step");

    expect_invalid(
        [&] {
            (void)fl::
                build_backward_euler_component_accumulation_residual(
                    pair,
                    valid,
                    std::numeric_limits<double>::quiet_NaN());
        },
        "time step");

    auto wrong_ids = valid;
    wrong_ids.component_ids[2] = "X";
    expect_invalid(
        [&] {
            (void)fl::
                build_backward_euler_component_accumulation_residual(
                    pair,
                    wrong_ids,
                    1.0);
        },
        "identity/porosity");

    auto wrong_porosity = valid;
    wrong_porosity.porosity = 0.24;
    expect_invalid(
        [&] {
            (void)fl::
                build_backward_euler_component_accumulation_residual(
                    pair,
                    wrong_porosity,
                    1.0);
        },
        "identity/porosity");

    auto wrong_shape = valid;
    wrong_shape.component_jacobian.pop_back();
    expect_invalid(
        [&] {
            (void)fl::
                build_backward_euler_component_accumulation_residual(
                    pair,
                    wrong_shape,
                    1.0);
        },
        "shape");

    auto nonfinite = valid;
    nonfinite.component_jacobian[0] =
        std::numeric_limits<double>::infinity();
    expect_invalid(
        [&] {
            (void)fl::
                build_backward_euler_component_accumulation_residual(
                    pair,
                    nonfinite,
                    1.0);
        },
        "non-finite");

    auto inconsistent = valid;
    inconsistent.component_jacobian[0] +=
        1.0;
    expect_invalid(
        [&] {
            (void)fl::
                build_backward_euler_component_accumulation_residual(
                    pair,
                    inconsistent,
                    1.0);
        },
        "closure");
}

void headers() {
    require(
        component_accumulation_time_header(),
        "component accumulation time header probe failed");
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"backward_euler", backward_euler},
    {"previous_is_frozen", previous_is_frozen},
    {"near_steady_roundoff", near_steady_roundoff},
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
