#include <mpmc/ad/dual.hpp>
#include <mpmc/ad/math.hpp>
#include <mpmc/flow/fugacity_equilibrium_residual.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <span>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool fugacity_equilibrium_residual_header();
bool thermodynamics_fugacity_adapters_header();

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
    double relative = 2.0e-12,
    double absolute = 2.0e-12,
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
        require(
            false,
            "numeric mismatch",
            where);
    }
}

template <class Function>
void expect_invalid(
    Function&& function,
    std::string_view fragment) {
    try {
        function();
    } catch (
        const std::invalid_argument& error) {
        require(
            std::string_view{
                error.what()}
                    .find(fragment) !=
                std::string_view::npos,
            "invalid_argument diagnostic changed");
        return;
    }
    throw std::runtime_error(
        "expected std::invalid_argument");
}

template <class Function>
void expect_range(
    Function&& function,
    std::string_view fragment) {
    try {
        function();
    } catch (
        const std::range_error& error) {
        require(
            std::string_view{
                error.what()}
                    .find(fragment) !=
                std::string_view::npos,
            "range_error diagnostic changed");
        return;
    }
    throw std::runtime_error(
        "expected std::range_error");
}

template <typename Number>
fl::PhaseLnFugacityCoefficientEvaluation<
    Number>
synthetic_ln_phi(
    fl::PhaseSlot3 slot,
    const Number& pressure_pa,
    const Number& temperature_k,
    std::span<const Number> composition) {
    const std::size_t phase =
        static_cast<std::size_t>(
            slot) +
        1U;

    std::vector<Number> values;
    values.reserve(composition.size());
    for (std::size_t component = 0U;
         component < composition.size();
         ++component) {
        const double p_coefficient =
            0.01 *
            static_cast<double>(
                phase *
                (component + 1U));
        const double t_coefficient =
            0.002 *
            static_cast<double>(
                component + 1U);
        const double x_coefficient =
            0.2 *
            static_cast<double>(
                phase);
        values.push_back(
            Number{p_coefficient} *
                (pressure_pa /
                 Number{1.0e6}) +
            Number{t_coefficient} *
                (temperature_k /
                 Number{300.0}) +
            Number{x_coefficient} *
                composition[component]);
    }
    return {std::move(values)};
}

struct SyntheticEvaluator {
    template <typename Number>
    fl::PhaseLnFugacityCoefficientEvaluation<
        Number>
    operator()(
        fl::PhaseSlot3 slot,
        const Number& pressure_pa,
        const Number& temperature_k,
        std::span<const Number>
            composition) const {
        return synthetic_ln_phi(
            slot,
            pressure_pa,
            temperature_k,
            composition);
    }
};

static_assert(
    fl::PhaseLnFugacityCoefficientEvaluator3P<
        SyntheticEvaluator,
        double>);

using Dual2 = ad::Dual<double, 2U>;

static_assert(
    fl::PhaseLnFugacityCoefficientEvaluator3P<
        SyntheticEvaluator,
        Dual2>);

struct DoubleStateStorage {
    std::array<
        std::vector<double>,
        3>
        compositions{
            std::vector<double>{
                0.20,
                0.30,
                0.50},
            std::vector<double>{
                0.25,
                0.35,
                0.40},
            std::vector<double>{
                0.10,
                0.55,
                0.35}};

    [[nodiscard]]
    fl::FugacityEquilibriumStateView3P<
        double>
    view(
        std::array<double, 3> pressures =
            {10.0e6, 11.0e6, 9.0e6},
        double temperature = 350.0) const {
        return {
            pressures,
            temperature,
            {
                std::span<const double>{
                    compositions[0]},
                std::span<const double>{
                    compositions[1]},
                std::span<const double>{
                    compositions[2]}}};
    }
};

double expected_ln_phi(
    std::size_t phase,
    std::size_t component,
    double pressure_pa,
    double temperature_k,
    double composition) {
    const double phase_factor =
        static_cast<double>(
            phase + 1U);
    const double component_factor =
        static_cast<double>(
            component + 1U);
    return
        0.01 *
            phase_factor *
            component_factor *
            (pressure_pa / 1.0e6) +
        0.002 *
            component_factor *
            (temperature_k / 300.0) +
        0.2 *
            phase_factor *
            composition;
}

void residual_values() {
    const DoubleStateStorage storage;
    const auto state = storage.view();
    const auto residual =
        fl::
            evaluate_fugacity_equilibrium_residual_3p(
                state,
                SyntheticEvaluator{});

    require(
        residual.component_count() == 3U &&
            residual.residual_count() == 6U &&
            residual.values().size() == 6U,
        "2*Nc residual shape changed");

    for (std::size_t phase = 1U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<fl::PhaseSlot3>(
                phase);
        for (std::size_t component = 0U;
             component < 3U;
             ++component) {
            const double expected =
                std::log(
                    storage
                        .compositions[0]
                                    [component] /
                    storage
                        .compositions[phase]
                                    [component]) +
                expected_ln_phi(
                    0U,
                    component,
                    state.phase_pressures_pa[
                        0U],
                    state.temperature_k,
                    storage
                        .compositions[0]
                                    [component]) -
                expected_ln_phi(
                    phase,
                    component,
                    state.phase_pressures_pa[
                        phase],
                    state.temperature_k,
                    storage
                        .compositions[phase]
                                    [component]) +
                std::log(
                    state.phase_pressures_pa[
                        0U] /
                    state.phase_pressures_pa[
                        phase]);

            near(
                residual.residual(
                    slot,
                    component),
                expected);

            const std::size_t
                local_row =
                    (phase - 1U) *
                        3U +
                    component;
            require(
                residual.local_row_index(
                    slot,
                    component) ==
                    local_row,
                "local fugacity row order changed");

            const fl::NaturalVariableLayout3P
                layout{3U};
            require(
                layout
                    .fugacity_equilibrium_equation_index(
                        slot,
                        component) ==
                    4U + local_row,
                "local/global fugacity row mapping changed");
        }
    }
}

struct RecordingEvaluator {
    std::array<double, 3>* pressures;

    fl::PhaseLnFugacityCoefficientEvaluation<
        double>
    operator()(
        fl::PhaseSlot3 slot,
        const double& pressure_pa,
        const double&,
        std::span<const double>
            composition) const {
        const std::size_t phase =
            static_cast<std::size_t>(
                slot);
        pressures->at(phase) =
            pressure_pa;
        return {
            std::vector<double>(
                composition.size(),
                0.0)};
    }
};

void actual_phase_pressure() {
    const DoubleStateStorage storage;
    const std::array<double, 3>
        pressures{
            7.5e6,
            8.25e6,
            6.75e6};
    const auto state =
        storage.view(
            pressures,
            340.0);

    std::array<double, 3> seen{
        0.0,
        0.0,
        0.0};
    const auto residual =
        fl::
            evaluate_fugacity_equilibrium_residual_3p(
                state,
                RecordingEvaluator{
                    &seen});

    require(
        seen == pressures,
        "thermodynamic evaluator did not receive actual phase pressures");

    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        near(
            residual.residual(
                fl::PhaseSlot3::phase1,
                component),
            std::log(
                storage.compositions[0]
                                    [component] /
                storage.compositions[1]
                                    [component]) +
                std::log(
                    pressures[0] /
                    pressures[1]));
        near(
            residual.residual(
                fl::PhaseSlot3::phase2,
                component),
            std::log(
                storage.compositions[0]
                                    [component] /
                storage.compositions[2]
                                    [component]) +
                std::log(
                    pressures[0] /
                    pressures[2]));
    }
}

void differentiable_interface() {
    const std::array<
        std::vector<Dual2>,
        3>
        compositions{
            std::vector<Dual2>{
                Dual2{0.20},
                Dual2{0.30},
                Dual2{0.50}},
            std::vector<Dual2>{
                Dual2::variable(
                    0.25,
                    1U),
                Dual2{0.35},
                Dual2{
                    0.40,
                    Dual2::Gradient{
                        0.0,
                        -1.0}}},
            std::vector<Dual2>{
                Dual2{0.10},
                Dual2{0.55},
                Dual2{0.35}}};

    const fl::
        FugacityEquilibriumStateView3P<
            Dual2>
        state{
            {
                Dual2{10.0e6},
                Dual2::variable(
                    11.0e6,
                    0U),
                Dual2{9.0e6}},
            Dual2{350.0},
            {
                std::span<const Dual2>{
                    compositions[0]},
                std::span<const Dual2>{
                    compositions[1]},
                std::span<const Dual2>{
                    compositions[2]}}};

    const auto residual =
        fl::
            evaluate_fugacity_equilibrium_residual_3p(
                state,
                SyntheticEvaluator{});

    const auto& phase1_component0 =
        residual.residual(
            fl::PhaseSlot3::phase1,
            0U);

    const double expected_dp =
        -(0.02 / 1.0e6 +
          1.0 / 11.0e6);
    const double expected_dx =
        -(1.0 / 0.25 +
          0.4);

    near(
        phase1_component0.derivative(
            0U),
        expected_dp,
        2.0e-12,
        2.0e-14);
    near(
        phase1_component0.derivative(
            1U),
        expected_dx);

    const auto& phase1_component2 =
        residual.residual(
            fl::PhaseSlot3::phase1,
            2U);
    near(
        phase1_component2.derivative(
            1U),
        1.0 / 0.40 + 0.4);

    const auto& phase2_component0 =
        residual.residual(
            fl::PhaseSlot3::phase2,
            0U);
    near(
        phase2_component0.derivative(
            0U),
        0.0,
        0.0,
        1.0e-15);
    near(
        phase2_component0.derivative(
            1U),
        0.0,
        0.0,
        1.0e-15);
}

struct WrongSizeEvaluator {
    template <typename Number>
    fl::PhaseLnFugacityCoefficientEvaluation<
        Number>
    operator()(
        fl::PhaseSlot3,
        const Number&,
        const Number&,
        std::span<const Number>) const {
        return {
            std::vector<Number>{
                Number{0.0}}};
    }
};

struct NonFiniteEvaluator {
    fl::PhaseLnFugacityCoefficientEvaluation<
        double>
    operator()(
        fl::PhaseSlot3,
        const double&,
        const double&,
        std::span<const double>
            composition) const {
        std::vector<double> values(
            composition.size(),
            0.0);
        values[0] =
            std::numeric_limits<double>::
                infinity();
        return {std::move(values)};
    }
};

struct ThrowingEvaluator {
    fl::PhaseLnFugacityCoefficientEvaluation<
        double>
    operator()(
        fl::PhaseSlot3,
        const double&,
        const double&,
        std::span<const double>) const {
        throw std::domain_error(
            "synthetic selected-root failure");
    }
};

void invalid_inputs() {
    const DoubleStateStorage storage;

    {
        auto state = storage.view();
        state.phase_pressures_pa[1] =
            0.0;
        expect_invalid(
            [&] {
                (void)fl::
                    evaluate_fugacity_equilibrium_residual_3p(
                        state,
                        SyntheticEvaluator{});
            },
            "actual phase pressure");
    }
    {
        auto state = storage.view();
        state.temperature_k = 0.0;
        expect_invalid(
            [&] {
                (void)fl::
                    evaluate_fugacity_equilibrium_residual_3p(
                        state,
                        SyntheticEvaluator{});
            },
            "temperature");
    }
    {
        auto bad = storage;
        bad.compositions[1][0] =
            0.0;
        const auto state = bad.view();
        expect_invalid(
            [&] {
                (void)fl::
                    evaluate_fugacity_equilibrium_residual_3p(
                        state,
                        SyntheticEvaluator{});
            },
            "strictly positive");
    }
    {
        auto bad = storage;
        bad.compositions[2][0] +=
            1.0e-3;
        const auto state = bad.view();
        expect_invalid(
            [&] {
                (void)fl::
                    evaluate_fugacity_equilibrium_residual_3p(
                        state,
                        SyntheticEvaluator{});
            },
            "must be normalized");
    }
    {
        const auto state = storage.view();
        expect_invalid(
            [&] {
                (void)fl::
                    evaluate_fugacity_equilibrium_residual_3p(
                        state,
                        WrongSizeEvaluator{});
            },
            "wrong component count");
    }
    {
        const auto state = storage.view();
        expect_range(
            [&] {
                (void)fl::
                    evaluate_fugacity_equilibrium_residual_3p(
                        state,
                        NonFiniteEvaluator{});
            },
            "non-finite");
    }
    {
        const auto state = storage.view();
        try {
            (void)fl::
                evaluate_fugacity_equilibrium_residual_3p(
                    state,
                    ThrowingEvaluator{});
        } catch (
            const std::domain_error& error) {
            require(
                std::string_view{
                    error.what()} ==
                    "synthetic selected-root failure",
                "thermodynamic exception was rewritten");
            return;
        }
        throw std::runtime_error(
            "thermodynamic evaluator failure was not propagated");
    }
}

void headers() {
    require(
        fugacity_equilibrium_residual_header(),
        "fugacity equilibrium public-header probe failed");
    require(
        thermodynamics_fugacity_adapters_header(),
        "thermodynamics fugacity adapter public-header probe failed");
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"residual_values", residual_values},
    {"actual_phase_pressure", actual_phase_pressure},
    {"differentiable_interface", differentiable_interface},
    {"invalid_inputs", invalid_inputs},
    {"headers", headers}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument(
                "one test name required");
        }
        for (const auto& [name, run] :
             tests) {
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
