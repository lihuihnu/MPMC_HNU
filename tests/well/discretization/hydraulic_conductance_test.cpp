#include <mpmc/well_discretization/hydraulic_conductance.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace flow = mpmc::flow;
namespace well = mpmc::well;
namespace wd = mpmc::well_discretization;

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
    double relative = 3.0e-13,
    double absolute = 3.0e-25,
    std::source_location where =
        std::source_location::current()) {
    const double scale =
        std::max(
            {1.0,
             std::abs(actual),
             std::abs(expected)});
    require(
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                absolute + relative * scale,
        "numeric mismatch",
        where);
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

well::PeacemanWellIndex3D connection() {
    return well::make_peaceman_well_index_3d(
        {20.0, 10.0, 5.0},
        {4.0e-13, 1.0e-13, 8.0e-14},
        well::AxisAlignedWellDirection3D::z,
        0.10,
        0.0);
}

flow::LocalPhaseMobilityLinearization3P
mobility_fixture() {
    flow::LocalPhaseMobilityLinearization3P
        result{};
    result.state_identity.layout =
        flow::NaturalVariableLayoutDescriptor{
            2U,
            3U,
            {1U, 1U, 1U}};
    result.state_identity.component_ids =
        {"A", "B"};
    result.state_identity.reference_pressure_pa =
        1.0e6;
    result.state_identity.temperature_k =
        350.0;
    result.state_identity.saturation =
        {0.2, 0.3, 0.5};
    result.state_identity.phase_composition = {
        std::vector<double>{0.7, 0.3},
        std::vector<double>{0.5, 0.5},
        std::vector<double>{0.2, 0.8}};

    const std::size_t q =
        result.state_identity.layout
            .unknown_count();
    result.mobility_per_pa_s =
        {1.0e-9, 2.0e-9, 3.0e-9};
    for (std::size_t phase = 0U;
         phase < flow::fixed_three_phase_count;
         ++phase) {
        result.mobility_gradient[phase]
            .resize(q);
        for (std::size_t column = 0U;
             column < q;
             ++column) {
            result.mobility_gradient[phase][column] =
                static_cast<double>(
                    (phase + 1U) *
                    (column + 1U)) *
                1.0e-12;
        }
    }
    return result;
}

void product_and_analytic_gradient() {
    const auto wi =
        connection();
    const auto mobility =
        mobility_fixture();
    const auto result =
        wd::make_hydraulic_conductance_linearization_3p(
            wi,
            mobility,
            flow::PhaseSlot3::phase1);

    require(
        result.phase ==
                flow::PhaseSlot3::phase1 &&
            result.state_identity.component_ids ==
                std::vector<std::string>{"A", "B"} &&
            result.hydraulic_conductance_gradient.size() ==
                7U,
        "hydraulic-conductance metadata changed");
    near(
        result.well_index_m3,
        wi.well_index_m3);
    near(
        result.phase_mobility_per_pa_s,
        2.0e-9);
    near(
        result.hydraulic_conductance_m3_per_pa_s,
        wi.well_index_m3 * 2.0e-9);
    for (std::size_t column = 0U;
         column < 7U;
         ++column) {
        near(
            result.d_hydraulic_conductance(
                column),
            wi.well_index_m3 *
                mobility
                    .mobility_gradient[1][column]);
    }
}

void phase_selection() {
    const auto wi =
        connection();
    const auto mobility =
        mobility_fixture();
    const auto phase0 =
        wd::make_hydraulic_conductance_linearization_3p(
            wi,
            mobility,
            flow::PhaseSlot3::phase0);
    const auto phase2 =
        wd::make_hydraulic_conductance_linearization_3p(
            wi,
            mobility,
            flow::PhaseSlot3::phase2);

    near(
        phase0.hydraulic_conductance_m3_per_pa_s,
        wi.well_index_m3 * 1.0e-9);
    near(
        phase2.hydraulic_conductance_m3_per_pa_s,
        wi.well_index_m3 * 3.0e-9);
    near(
        phase0.d_hydraulic_conductance(3U),
        wi.well_index_m3 * 4.0e-12);
    near(
        phase2.d_hydraulic_conductance(3U),
        wi.well_index_m3 * 12.0e-12);
}

void zero_mobility_nonzero_gradient() {
    const auto wi =
        connection();
    auto mobility =
        mobility_fixture();
    mobility.mobility_per_pa_s[2] =
        0.0;
    mobility.mobility_gradient[2][0] =
        4.0e-12;

    const auto result =
        wd::make_hydraulic_conductance_linearization_3p(
            wi,
            mobility,
            flow::PhaseSlot3::phase2);

    require(
        result.hydraulic_conductance_m3_per_pa_s ==
            0.0,
        "zero mobility did not preserve exact zero conductance");
    near(
        result.d_hydraulic_conductance(0U),
        wi.well_index_m3 * 4.0e-12);
    require(
        result.d_hydraulic_conductance(0U) !=
            0.0,
        "zero mobility incorrectly erased analytic derivative");
}

void invalid_inputs() {
    const auto wi =
        connection();

    {
        auto mobility =
            mobility_fixture();
        mobility.mobility_per_pa_s[1] =
            -1.0;
        expect_invalid(
            [&] {
                (void)wd::
                    make_hydraulic_conductance_linearization_3p(
                        wi,
                        mobility,
                        flow::PhaseSlot3::phase1);
            });
    }
    {
        auto mobility =
            mobility_fixture();
        mobility.mobility_gradient[1]
            .pop_back();
        expect_invalid(
            [&] {
                (void)wd::
                    make_hydraulic_conductance_linearization_3p(
                        wi,
                        mobility,
                        flow::PhaseSlot3::phase1);
            });
    }
    {
        auto mobility =
            mobility_fixture();
        mobility.mobility_gradient[1][2] =
            std::numeric_limits<double>::
                quiet_NaN();
        expect_invalid(
            [&] {
                (void)wd::
                    make_hydraulic_conductance_linearization_3p(
                        wi,
                        mobility,
                        flow::PhaseSlot3::phase1);
            });
    }
    {
        auto invalid_wi = wi;
        invalid_wi.well_index_m3 =
            0.0;
        expect_invalid(
            [&] {
                (void)wd::
                    make_hydraulic_conductance_linearization_3p(
                        invalid_wi,
                        mobility_fixture(),
                        flow::PhaseSlot3::phase1);
            });
    }
    {
        expect_invalid(
            [&] {
                (void)wd::
                    make_hydraulic_conductance_linearization_3p(
                        wi,
                        mobility_fixture(),
                        static_cast<
                            flow::PhaseSlot3>(
                                99));
            });
    }
    {
        auto mobility =
            mobility_fixture();
        mobility.state_identity.layout =
            flow::NaturalVariableLayoutDescriptor{
                2U,
                2U,
                {1U, 1U}};
        expect_invalid(
            [&] {
                (void)wd::
                    make_hydraulic_conductance_linearization_3p(
                        wi,
                        mobility,
                        flow::PhaseSlot3::phase1);
            });
    }
}

using Test =
    std::pair<
        std::string_view,
        void (*)()>;

constexpr Test tests[]{
    {"product_and_analytic_gradient",
     product_and_analytic_gradient},
    {"phase_selection",
     phase_selection},
    {"zero_mobility_nonzero_gradient",
     zero_mobility_nonzero_gradient},
    {"invalid_inputs",
     invalid_inputs}};

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
