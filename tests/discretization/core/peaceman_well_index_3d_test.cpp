#include <mpmc/discretization/peaceman_well_index_3d.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace {

namespace discretization = mpmc::discretization;
namespace mesh = mpmc::mesh;

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

void require_close(
    double actual,
    double expected,
    double relative,
    double absolute,
    std::string_view message,
    std::source_location where =
        std::source_location::current()) {
    const double scale =
        std::max(
            {1.0,
             std::abs(actual),
             std::abs(expected)});
    if (!std::isfinite(actual) ||
        !std::isfinite(expected) ||
        std::abs(actual - expected) >
            absolute + relative * scale) {
        throw std::runtime_error(
            std::string{where.file_name()} +
            ":" +
            std::to_string(where.line()) +
            ": " +
            std::string{message});
    }
}

template <class Function>
void expect_invalid(Function&& function) {
    bool caught = false;
    try {
        std::forward<Function>(function)();
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    require(caught, "expected invalid_argument");
}

void isotropic_square_limit() {
    constexpr double k = 100.0e-15;
    const auto result =
        discretization::make_peaceman_well_index_3d(
            {10.0, 10.0, 5.0},
            {k, k, 7.0e-15},
            discretization::
                AxisAlignedWellDirection3D::z,
            0.10,
            0.0);

    const double expected_r0 =
        0.14 *
        std::sqrt(10.0 * 10.0 +
                  10.0 * 10.0);
    const double expected_wi =
        2.0 *
        std::numbers::pi_v<double> *
        k *
        5.0 /
        std::log(expected_r0 / 0.10);

    require_close(
        result.equivalent_radius_m,
        expected_r0,
        2.0e-15,
        1.0e-15,
        "isotropic square Peaceman radius");
    require_close(
        result.effective_radial_permeability_m2,
        k,
        2.0e-15,
        1.0e-28,
        "isotropic square effective permeability");
    require_close(
        result.completion_length_m,
        5.0,
        0.0,
        0.0,
        "Z-well completion length");
    require_close(
        result.well_index_m3,
        expected_wi,
        3.0e-15,
        1.0e-27,
        "isotropic square Peaceman WI");
}

void anisotropic_rectangular_z() {
    const auto result =
        discretization::make_peaceman_well_index_3d(
            {20.0, 10.0, 5.0},
            {4.0e-13, 1.0e-13, 8.0e-14},
            discretization::
                AxisAlignedWellDirection3D::z,
            0.10,
            0.0);

    // For Ky/Kx=1/4:
    // numerator=sqrt(20^2*1/2 + 10^2*2)=20,
    // denominator=1/sqrt(2)+sqrt(2)=3/sqrt(2).
    const double expected_r0 =
        5.6 *
        std::sqrt(2.0) /
        3.0;
    const double expected_kh =
        2.0e-13;
    const double expected_wi =
        2.0 *
        std::numbers::pi_v<double> *
        expected_kh *
        5.0 /
        std::log(expected_r0 / 0.10);

    require_close(
        result.equivalent_radius_m,
        expected_r0,
        3.0e-15,
        1.0e-15,
        "anisotropic rectangular equivalent radius");
    require_close(
        result.effective_radial_permeability_m2,
        expected_kh,
        2.0e-15,
        1.0e-28,
        "anisotropic rectangular kh");
    require_close(
        result.well_index_m3,
        expected_wi,
        4.0e-15,
        1.0e-27,
        "anisotropic rectangular WI");
}

void axis_permutation_and_cell_binding() {
    const mesh::CellCartesianDiagonalPermeability3D
        permeability{
            {
                mesh::CartesianDiagonalPermeabilityTensor3D{
                    8.0e-14,
                    4.0e-13,
                    1.0e-13}
            }};

    const auto x_connection =
        discretization::make_cell_peaceman_well_index_3d(
            permeability,
            mesh::LocalIndex{0U},
            {5.0, 20.0, 10.0},
            discretization::
                AxisAlignedWellDirection3D::x,
            0.10,
            0.0);

    const auto z_reference =
        discretization::make_peaceman_well_index_3d(
            {20.0, 10.0, 5.0},
            {4.0e-13, 1.0e-13, 8.0e-14},
            discretization::
                AxisAlignedWellDirection3D::z,
            0.10,
            0.0);

    require(
        x_connection.cell ==
            mesh::LocalIndex{0U},
        "cell Peaceman binding lost LocalIndex");
    require_close(
        x_connection.connection.equivalent_radius_m,
        z_reference.equivalent_radius_m,
        2.0e-15,
        1.0e-15,
        "axis permutation equivalent radius");
    require_close(
        x_connection.connection
            .effective_radial_permeability_m2,
        z_reference
            .effective_radial_permeability_m2,
        2.0e-15,
        1.0e-28,
        "axis permutation kh");
    require_close(
        x_connection.connection.completion_length_m,
        z_reference.completion_length_m,
        0.0,
        0.0,
        "axis permutation completion length");
    require_close(
        x_connection.connection.well_index_m3,
        z_reference.well_index_m3,
        3.0e-15,
        1.0e-27,
        "axis permutation WI");
}

void skin_factor() {
    const auto zero_skin =
        discretization::make_peaceman_well_index_3d(
            {20.0, 10.0, 5.0},
            {4.0e-13, 1.0e-13, 8.0e-14},
            discretization::
                AxisAlignedWellDirection3D::z,
            0.10,
            0.0);
    const auto positive_skin =
        discretization::make_peaceman_well_index_3d(
            {20.0, 10.0, 5.0},
            {4.0e-13, 1.0e-13, 8.0e-14},
            discretization::
                AxisAlignedWellDirection3D::z,
            0.10,
            2.0);

    require_close(
        positive_skin.logarithmic_denominator,
        zero_skin.logarithmic_denominator + 2.0,
        2.0e-15,
        1.0e-15,
        "skin denominator");
    require(
        positive_skin.well_index_m3 <
            zero_skin.well_index_m3,
        "positive skin must reduce WI");
}

void invalid_inputs() {
    const auto valid_k =
        mesh::CartesianDiagonalPermeabilityTensor3D{
            4.0e-13,
            1.0e-13,
            8.0e-14};

    expect_invalid(
        [&] {
            (void)discretization::
                make_peaceman_well_index_3d(
                    {0.0, 10.0, 5.0},
                    valid_k,
                    discretization::
                        AxisAlignedWellDirection3D::z,
                    0.1);
        });
    expect_invalid(
        [&] {
            (void)discretization::
                make_peaceman_well_index_3d(
                    {20.0, 10.0, 5.0},
                    {4.0e-13, 0.0, 8.0e-14},
                    discretization::
                        AxisAlignedWellDirection3D::z,
                    0.1);
        });
    expect_invalid(
        [&] {
            (void)discretization::
                make_peaceman_well_index_3d(
                    {20.0, 10.0, 5.0},
                    {-4.0e-13, 1.0e-13, 8.0e-14},
                    discretization::
                        AxisAlignedWellDirection3D::z,
                    0.1);
        });
    expect_invalid(
        [&] {
            (void)discretization::
                make_peaceman_well_index_3d(
                    {20.0, 10.0, 5.0},
                    valid_k,
                    discretization::
                        AxisAlignedWellDirection3D::z,
                    0.0);
        });
    expect_invalid(
        [&] {
            (void)discretization::
                make_peaceman_well_index_3d(
                    {10.0, 10.0, 5.0},
                    {
                        100.0e-15,
                        100.0e-15,
                        10.0e-15},
                    discretization::
                        AxisAlignedWellDirection3D::z,
                    10.0,
                    0.0);
        });
    expect_invalid(
        [&] {
            (void)discretization::
                make_peaceman_well_index_3d(
                    {20.0, 10.0, 5.0},
                    valid_k,
                    static_cast<
                        discretization::
                            AxisAlignedWellDirection3D>(
                                99),
                    0.1);
        });
    expect_invalid(
        [&] {
            (void)discretization::
                make_peaceman_well_index_3d(
                    {
                        std::numeric_limits<double>::
                            quiet_NaN(),
                        10.0,
                        5.0},
                    valid_k,
                    discretization::
                        AxisAlignedWellDirection3D::z,
                    0.1);
        });
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"isotropic_square_limit",
     isotropic_square_limit},
    {"anisotropic_rectangular_z",
     anisotropic_rectangular_z},
    {"axis_permutation_and_cell_binding",
     axis_permutation_and_cell_binding},
    {"skin_factor", skin_factor},
    {"invalid_inputs", invalid_inputs}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument(
                "one test name required");
        }
        const std::string_view name{
            argv[1]};
        for (const auto& [test_name, run] :
             tests) {
            if (name == test_name) {
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
