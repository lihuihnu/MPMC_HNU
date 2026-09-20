#include <mpmc/ad/dual.hpp>
#include <mpmc/flow/saturation_constitutive.hpp>

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

bool saturation_constitutive_header();

namespace {

namespace ad = mpmc::ad;
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
    double relative = 3.0e-13,
    double absolute = 3.0e-13,
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

struct SyntheticRelativePermeability3P {
    template <typename Number>
    fl::RelativePermeabilityEvaluation3P<Number>
    operator()(
        const fl::ThreePhaseSaturationState3P<Number>&
            state) const {
        const auto& s = state.saturation;
        return {
            std::array<Number, 3>{
                s[0] * s[0],
                s[1] * s[1],
                s[2] * s[2]}};
    }
};

struct SyntheticCapillaryPressure3P {
    template <typename Number>
    fl::CapillaryPressureOffsetsEvaluation3P<Number>
    operator()(
        const fl::ThreePhaseSaturationState3P<Number>&
            state) const {
        const auto& s = state.saturation;
        return {
            std::array<Number, 2>{
                Number{2.0e5} *
                    (s[1] - s[0]),
                Number{3.0e5} *
                    (s[2] - s[0])}};
    }
};

static_assert(
    fl::RelativePermeabilityEvaluator3P<
        SyntheticRelativePermeability3P,
        double>);
static_assert(
    fl::CapillaryPressureEvaluator3P<
        SyntheticCapillaryPressure3P,
        double>);
static_assert(
    fl::CapillaryPressureEvaluator3P<
        fl::NoCapillaryPressure3P,
        double>);

void primal() {
    const auto evaluated =
        fl::evaluate_three_phase_saturation_constitutive(
            10.0e6,
            0.20,
            0.30,
            SyntheticRelativePermeability3P{},
            SyntheticCapillaryPressure3P{});

    near(
        evaluated.saturation_state.saturation[0],
        0.20);
    near(
        evaluated.saturation_state.saturation[1],
        0.30);
    near(
        evaluated.saturation_state.saturation[2],
        0.50);

    near(
        evaluated.relative_permeability[0],
        0.04);
    near(
        evaluated.relative_permeability[1],
        0.09);
    near(
        evaluated.relative_permeability[2],
        0.25);

    near(
        evaluated.capillary_pressure_offset_pa[0],
        0.0);
    near(
        evaluated.capillary_pressure_offset_pa[1],
        2.0e4);
    near(
        evaluated.capillary_pressure_offset_pa[2],
        3.0e4);

    near(
        evaluated.phase_pressure_pa[0],
        10.0e6);
    near(
        evaluated.phase_pressure_pa[1],
        10.02e6);
    near(
        evaluated.phase_pressure_pa[2],
        10.03e6);
}

void differentiable_chart() {
    using D = ad::Dual<double, 2U>;

    const D s0 =
        D::variable(0.20, 0U);
    const D s1 =
        D::variable(0.30, 1U);

    const auto evaluated =
        fl::evaluate_three_phase_saturation_constitutive(
            D{10.0e6},
            s0,
            s1,
            SyntheticRelativePermeability3P{},
            SyntheticCapillaryPressure3P{});

    const auto& s2 =
        evaluated.saturation_state.saturation[2];
    near(s2.value(), 0.50);
    near(s2.derivative(0U), -1.0);
    near(s2.derivative(1U), -1.0);

    near(
        evaluated.relative_permeability[0]
            .derivative(0U),
        0.40);
    near(
        evaluated.relative_permeability[0]
            .derivative(1U),
        0.0);
    near(
        evaluated.relative_permeability[1]
            .derivative(0U),
        0.0);
    near(
        evaluated.relative_permeability[1]
            .derivative(1U),
        0.60);
    near(
        evaluated.relative_permeability[2]
            .derivative(0U),
        -1.0);
    near(
        evaluated.relative_permeability[2]
            .derivative(1U),
        -1.0);

    near(
        evaluated.capillary_pressure_offset_pa[0]
            .derivative(0U),
        0.0);
    near(
        evaluated.capillary_pressure_offset_pa[0]
            .derivative(1U),
        0.0);

    near(
        evaluated.capillary_pressure_offset_pa[1]
            .derivative(0U),
        -2.0e5);
    near(
        evaluated.capillary_pressure_offset_pa[1]
            .derivative(1U),
        2.0e5);

    near(
        evaluated.capillary_pressure_offset_pa[2]
            .derivative(0U),
        -6.0e5);
    near(
        evaluated.capillary_pressure_offset_pa[2]
            .derivative(1U),
        -3.0e5);

    near(
        evaluated.phase_pressure_pa[0]
            .derivative(0U),
        0.0);
    near(
        evaluated.phase_pressure_pa[1]
            .derivative(0U),
        -2.0e5);
    near(
        evaluated.phase_pressure_pa[1]
            .derivative(1U),
        2.0e5);
    near(
        evaluated.phase_pressure_pa[2]
            .derivative(0U),
        -6.0e5);
    near(
        evaluated.phase_pressure_pa[2]
            .derivative(1U),
        -3.0e5);
}

void no_capillary() {
    using D = ad::Dual<double, 2U>;
    const auto evaluated =
        fl::evaluate_three_phase_saturation_constitutive(
            D{8.0e6},
            D::variable(0.25, 0U),
            D::variable(0.35, 1U),
            SyntheticRelativePermeability3P{},
            fl::NoCapillaryPressure3P{});

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        near(
            evaluated
                .capillary_pressure_offset_pa[
                    phase]
                .value(),
            0.0);
        near(
            evaluated
                .capillary_pressure_offset_pa[
                    phase]
                .derivative(0U),
            0.0);
        near(
            evaluated
                .capillary_pressure_offset_pa[
                    phase]
                .derivative(1U),
            0.0);
        near(
            evaluated.phase_pressure_pa[phase]
                .value(),
            8.0e6);
        near(
            evaluated.phase_pressure_pa[phase]
                .derivative(0U),
            0.0);
        near(
            evaluated.phase_pressure_pa[phase]
                .derivative(1U),
            0.0);
    }
}

struct NegativeRelativePermeability3P {
    template <typename Number>
    fl::RelativePermeabilityEvaluation3P<Number>
    operator()(
        const fl::ThreePhaseSaturationState3P<Number>&) const {
        return {
            std::array<Number, 3>{
                Number{0.1},
                Number{-0.1},
                Number{0.2}}};
    }
};

struct NonFiniteCapillaryPressure3P {
    fl::CapillaryPressureOffsetsEvaluation3P<double>
    operator()(
        const fl::ThreePhaseSaturationState3P<double>&) const {
        return {
            std::array<double, 2>{
                0.0,
                std::numeric_limits<double>::
                    infinity()}};
    }
};

struct NonPositiveResolvedPressure3P {
    fl::CapillaryPressureOffsetsEvaluation3P<double>
    operator()(
        const fl::ThreePhaseSaturationState3P<double>&) const {
        return {
            std::array<double, 2>{
                -11.0e6,
                0.0}};
    }
};

void invalid_inputs() {
    expect_invalid(
        [] {
            (void)fl::
                evaluate_three_phase_saturation_constitutive(
                    10.0e6,
                    0.0,
                    0.30,
                    SyntheticRelativePermeability3P{},
                    fl::NoCapillaryPressure3P{});
        },
        "phase0 saturation");

    expect_invalid(
        [] {
            (void)fl::
                evaluate_three_phase_saturation_constitutive(
                    10.0e6,
                    0.70,
                    0.30,
                    SyntheticRelativePermeability3P{},
                    fl::NoCapillaryPressure3P{});
        },
        "phase2 saturation");

    expect_invalid(
        [] {
            (void)fl::
                evaluate_three_phase_saturation_constitutive(
                    0.0,
                    0.20,
                    0.30,
                    SyntheticRelativePermeability3P{},
                    fl::NoCapillaryPressure3P{});
        },
        "reference phase pressure");

    expect_invalid(
        [] {
            (void)fl::
                evaluate_three_phase_saturation_constitutive(
                    10.0e6,
                    0.20,
                    0.30,
                    NegativeRelativePermeability3P{},
                    fl::NoCapillaryPressure3P{});
        },
        "relative permeability");

    expect_invalid(
        [] {
            (void)fl::
                evaluate_three_phase_saturation_constitutive(
                    10.0e6,
                    0.20,
                    0.30,
                    SyntheticRelativePermeability3P{},
                    NonFiniteCapillaryPressure3P{});
        },
        "capillary pressure");

    expect_invalid(
        [] {
            (void)fl::
                evaluate_three_phase_saturation_constitutive(
                    10.0e6,
                    0.20,
                    0.30,
                    SyntheticRelativePermeability3P{},
                    NonPositiveResolvedPressure3P{});
        },
        "resolved phase pressure");
}

void headers() {
    require(
        saturation_constitutive_header(),
        "saturation constitutive header probe failed");
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"primal", primal},
    {"differentiable_chart", differentiable_chart},
    {"no_capillary", no_capillary},
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
