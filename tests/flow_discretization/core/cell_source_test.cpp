#include <mpmc/flow_discretization/cell_source.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

bool cell_source_header();

namespace {
namespace fd = mpmc::flow_discretization;

void require(
    bool condition,
    std::string_view message,
    std::source_location where =
        std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(
            std::string{where.file_name()} + ":" +
            std::to_string(where.line()) + ": " +
            std::string{message});
    }
}

void near(
    double actual,
    double expected,
    std::source_location where =
        std::source_location::current()) {
    require(
        std::isfinite(actual) &&
            std::abs(actual - expected) <=
                1.0e-13 *
                    std::max(
                        {1.0,
                         std::abs(actual),
                         std::abs(expected)}),
        "numeric mismatch",
        where);
}

fd::CellSourceLinearization3D source() {
    return {
        "test/cell-source/v1",
        {"A", "B"},
        3U,
        {2.0, -1.0},
        {
            0.3, -0.6, 0.9,
            -0.2, 0.4, -0.8},
        50.0,
        {3.0, -4.0, 5.0}};
}

void normalized_sign_units_and_jacobian() {
    const auto value =
        fd::normalize_cell_source_by_bulk_volume(
            source(),
            2.0);
    require(
        value.component_ids ==
                std::vector<std::string>{"A", "B"} &&
            value.input_count == 3U &&
            value.bulk_volume_m3 == 2.0,
        "cell source metadata changed");
    near(value.component_residual_mol_per_bulk_m3_s[0], -1.0);
    near(value.component_residual_mol_per_bulk_m3_s[1], 0.5);
    near(value.d_component_residual(0U, 0U), -0.15);
    near(value.d_component_residual(1U, 2U), 0.4);
    near(value.energy_residual_w_per_bulk_m3, -25.0);
    near(value.d_energy_residual(1U), 2.0);
}

void invalid_inputs() {
    auto bad = source();
    bad.component_molar_rate_jacobian_mol_per_s.pop_back();
    bool rejected = false;
    try {
        (void)fd::normalize_cell_source_by_bulk_volume(
            bad,
            2.0);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "malformed source was accepted");

    rejected = false;
    try {
        (void)fd::normalize_cell_source_by_bulk_volume(
            source(),
            0.0);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "zero bulk volume was accepted");
}

void headers() {
    require(
        cell_source_header(),
        "cell source header probe failed");
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument(
                "one test name required");
        }
        const std::string_view name{argv[1]};
        if (name ==
            "normalized_sign_units_and_jacobian") {
            normalized_sign_units_and_jacobian();
        } else if (name == "invalid_inputs") {
            invalid_inputs();
        } else if (name == "headers") {
            headers();
        } else {
            throw std::invalid_argument(
                "unknown test");
        }
        std::cout << "[PASS] " << name << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
