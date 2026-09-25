#include <mpmc/ad/dual.hpp>
#include <mpmc/flow/saturation_constitutive.hpp>

#include "spe3_kenyon_behie_saturation_reference.hpp"

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
namespace spe3 = mpmc::flow::test_reference::spe3;

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

template <typename Number>
[[nodiscard]] double
reference_primal(const Number& value) {
    if constexpr (requires { value.value(); }) {
        return static_cast<double>(value.value());
    } else {
        return static_cast<double>(value);
    }
}

template <typename Number>
[[nodiscard]] Number
reference_linear_segment(
    const Number& coordinate,
    double x0,
    double y0,
    double x1,
    double y1,
    std::string_view quantity) {
    const double x =
        reference_primal(coordinate);
    if (!std::isfinite(x) ||
        x < x0 ||
        x > x1) {
        throw std::invalid_argument(
            std::string{quantity} +
            " outside SPE3 reference bracket");
    }
    const double slope =
        (y1 - y0) / (x1 - x0);
    return Number{y0} +
        (coordinate - Number{x0}) *
            slope;
}

/// Test-only table evaluator.  It is intentionally restricted to the source
/// brackets required for the registered SPE3 point and is not a production
/// interpolation/extrapolation model.
struct Spe3ReferenceRelativePermeability3P {
    template <typename Number>
    [[nodiscard]]
    fl::RelativePermeabilityEvaluation3P<Number>
    operator()(
        const fl::ThreePhaseSaturationState3P<Number>&
            state) const {
        const auto& s = state.saturation;

        const double so0 =
            1.0 -
            spe3::oil_low_saturation_source
                .water_saturation;
        const double so1 =
            1.0 -
            spe3::oil_high_saturation_source
                .water_saturation;

        const Number kro =
            reference_linear_segment(
                s[0],
                so0,
                spe3::oil_low_saturation_source
                    .oil_relative_permeability,
                so1,
                spe3::oil_high_saturation_source
                    .oil_relative_permeability,
                "oil saturation");
        const Number krw =
            reference_linear_segment(
                s[1],
                spe3::water_lower.water_saturation,
                spe3::water_lower
                    .water_relative_permeability,
                spe3::water_upper.water_saturation,
                spe3::water_upper
                    .water_relative_permeability,
                "water saturation");
        const Number krg =
            reference_linear_segment(
                s[2],
                spe3::gas_lower.gas_saturation,
                spe3::gas_lower
                    .gas_relative_permeability,
                spe3::gas_upper.gas_saturation,
                spe3::gas_upper
                    .gas_relative_permeability,
                "gas saturation");

        return {
            std::array<Number, 3>{
                kro,
                krw,
                krg}};
    }
};

struct Spe3ReferenceCapillaryPressure3P {
    template <typename Number>
    [[nodiscard]]
    fl::CapillaryPressureOffsetsEvaluation3P<Number>
    operator()(
        const fl::ThreePhaseSaturationState3P<Number>&
            state) const {
        const Number pcow_psi =
            reference_linear_segment(
                state.saturation[1],
                spe3::water_lower.water_saturation,
                spe3::water_lower.pcow_psi,
                spe3::water_upper.water_saturation,
                spe3::water_upper.pcow_psi,
                "water saturation");
        const Number pcgo_psi =
            reference_linear_segment(
                state.saturation[2],
                spe3::gas_lower.gas_saturation,
                spe3::gas_lower.pcgo_psi,
                spe3::gas_upper.gas_saturation,
                spe3::gas_upper.pcgo_psi,
                "gas saturation");

        // The repository contract stores offsets relative to phase0=oil.
        // SPE/OPM signs are Pcow=Po-Pw and Pcgo=Pg-Po.
        return {
            std::array<Number, 2>{
                -pcow_psi *
                    spe3::psi_to_pa,
                pcgo_psi *
                    spe3::psi_to_pa}};
    }
};

static_assert(
    fl::RelativePermeabilityEvaluator3P<
        Spe3ReferenceRelativePermeability3P,
        double>);
static_assert(
    fl::CapillaryPressureEvaluator3P<
        Spe3ReferenceCapillaryPressure3P,
        double>);

void reference_backed_spe3() {
    constexpr double reference_pressure_pa =
        20.0e6;

    const auto primal =
        fl::evaluate_three_phase_saturation_constitutive(
            reference_pressure_pa,
            spe3::oil_saturation,
            spe3::water_saturation,
            Spe3ReferenceRelativePermeability3P{},
            Spe3ReferenceCapillaryPressure3P{});

    near(
        primal.saturation_state.saturation[2],
        spe3::gas_saturation);
    near(
        primal.relative_permeability[0],
        spe3::oil_relative_permeability);
    near(
        primal.relative_permeability[1],
        spe3::water_relative_permeability);
    near(
        primal.relative_permeability[2],
        spe3::gas_relative_permeability);
    near(
        primal.capillary_pressure_offset_pa[0],
        0.0);
    near(
        primal.capillary_pressure_offset_pa[1],
        spe3::water_minus_oil_pressure_pa,
        3.0e-13,
        1.0e-8);
    near(
        primal.capillary_pressure_offset_pa[2],
        spe3::gas_minus_oil_pressure_pa);
    near(
        primal.phase_pressure_pa[1],
        reference_pressure_pa +
            spe3::water_minus_oil_pressure_pa,
        3.0e-13,
        1.0e-8);
    near(
        primal.phase_pressure_pa[2],
        reference_pressure_pa);

    using D = ad::Dual<double, 2U>;
    const auto differentiated =
        fl::evaluate_three_phase_saturation_constitutive(
            D{reference_pressure_pa},
            D::variable(
                spe3::oil_saturation,
                0U),
            D::variable(
                spe3::water_saturation,
                1U),
            Spe3ReferenceRelativePermeability3P{},
            Spe3ReferenceCapillaryPressure3P{});

    near(
        differentiated
            .relative_permeability[0]
            .derivative(0U),
        spe3::d_kro_d_so);
    near(
        differentiated
            .relative_permeability[0]
            .derivative(1U),
        0.0);
    near(
        differentiated
            .relative_permeability[1]
            .derivative(0U),
        0.0);
    near(
        differentiated
            .relative_permeability[1]
            .derivative(1U),
        spe3::d_krw_d_sw);
    near(
        differentiated
            .relative_permeability[2]
            .derivative(0U),
        -spe3::d_krg_d_sg);
    near(
        differentiated
            .relative_permeability[2]
            .derivative(1U),
        -spe3::d_krg_d_sg);
    near(
        differentiated
            .capillary_pressure_offset_pa[1]
            .derivative(0U),
        0.0);
    near(
        differentiated
            .capillary_pressure_offset_pa[1]
            .derivative(1U),
        spe3::
            d_water_minus_oil_pressure_d_sw_pa,
        3.0e-13,
        1.0e-7);
    near(
        differentiated
            .capillary_pressure_offset_pa[2]
            .derivative(0U),
        0.0);
    near(
        differentiated
            .capillary_pressure_offset_pa[2]
            .derivative(1U),
        0.0);

    // Both bracket-endpoint states lie exactly on source table rows and must
    // remain admissible; no extrapolation is needed for this regression.
    const auto endpoint_a =
        fl::evaluate_three_phase_saturation_constitutive(
            reference_pressure_pa,
            0.56,
            0.28,
            Spe3ReferenceRelativePermeability3P{},
            Spe3ReferenceCapillaryPressure3P{});
    near(
        endpoint_a.relative_permeability[0],
        0.150);
    near(
        endpoint_a.relative_permeability[1],
        0.020);
    near(
        endpoint_a.relative_permeability[2],
        0.040);
    near(
        endpoint_a.capillary_pressure_offset_pa[1],
        -15.5 * spe3::psi_to_pa,
        3.0e-13,
        1.0e-8);

    const auto endpoint_b =
        fl::evaluate_three_phase_saturation_constitutive(
            reference_pressure_pa,
            0.52,
            0.32,
            Spe3ReferenceRelativePermeability3P{},
            Spe3ReferenceCapillaryPressure3P{});
    near(
        endpoint_b.relative_permeability[0],
        0.112);
    near(
        endpoint_b.relative_permeability[1],
        0.033);
    near(
        endpoint_b.relative_permeability[2],
        0.040);
    near(
        endpoint_b.capillary_pressure_offset_pa[1],
        -12.0 * spe3::psi_to_pa,
        3.0e-13,
        1.0e-8);

    expect_invalid(
        [&] {
            (void)fl::
                evaluate_three_phase_saturation_constitutive(
                    reference_pressure_pa,
                    0.57,
                    0.28,
                    Spe3ReferenceRelativePermeability3P{},
                    Spe3ReferenceCapillaryPressure3P{});
        },
        "oil saturation outside SPE3 reference bracket");
}

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
        9.0e4);

    near(
        evaluated.phase_pressure_pa[0],
        10.0e6);
    near(
        evaluated.phase_pressure_pa[1],
        10.02e6);
    near(
        evaluated.phase_pressure_pa[2],
        10.09e6);
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
                    0.80,
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
    {"reference_backed_spe3", reference_backed_spe3},
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
