#include <mpmc/ad/dual.hpp>
#include <mpmc/flow/phase_transport.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <source_location>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool phase_transport_header();

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
    double relative,
    double absolute,
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
        const fl::ThreePhaseSaturationState3P<Number>& state) const {
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
        const fl::ThreePhaseSaturationState3P<Number>& state) const {
        const auto& s = state.saturation;
        return {
            std::array<Number, 2>{
                Number{2.0e5} *
                    (s[1] - s[0]),
                Number{3.0e5} *
                    (s[2] - s[0])}};
    }
};

struct SyntheticTransportLaw3P {
    std::array<double, 3> rho_base{
        700.0, 760.0, 820.0};
    std::array<double, 3> rho_p{
        1.2e-5, 0.9e-5, 1.1e-5};
    std::array<double, 3> rho_t{
        0.20, 0.25, 0.18};
    std::array<double, 3> rho_s0{
        15.0, 20.0, 12.0};
    std::array<double, 3> rho_s1{
        10.0, 14.0, 18.0};

    std::array<double, 3> mu_base{
        1.0e-3, 1.4e-3, 1.8e-3};
    std::array<double, 3> mu_p{
        1.0e-10, 0.8e-10, 1.2e-10};
    std::array<double, 3> mu_t{
        2.0e-6, 1.5e-6, 2.5e-6};
    std::array<double, 3> mu_s0{
        4.0e-4, 3.0e-4, 2.0e-4};
    std::array<double, 3> mu_s1{
        2.0e-4, 4.0e-4, 3.0e-4};

    std::array<std::array<double, 3>, 3>
        rho_x{{
            {30.0, 20.0, 10.0},
            {15.0, 35.0, 25.0},
            {22.0, 12.0, 32.0}}};
    std::array<std::array<double, 3>, 3>
        mu_x{{
            {3.0e-4, 2.0e-4, 1.0e-4},
            {1.0e-4, 4.0e-4, 2.0e-4},
            {2.5e-4, 1.5e-4, 3.5e-4}}};

    double rho(
        std::size_t phase,
        double pressure,
        double temperature,
        double s0,
        double s1,
        std::span<const double> x) const {
        double value =
            rho_base[phase] +
            rho_p[phase] * pressure +
            rho_t[phase] * temperature +
            rho_s0[phase] * s0 +
            rho_s1[phase] * s1;
        for (std::size_t component = 0U;
             component < x.size();
             ++component) {
            value +=
                rho_x[phase][component] *
                x[component];
        }
        return value;
    }

    double mu(
        std::size_t phase,
        double pressure,
        double temperature,
        double s0,
        double s1,
        std::span<const double> x) const {
        double value =
            mu_base[phase] +
            mu_p[phase] * pressure +
            mu_t[phase] * temperature +
            mu_s0[phase] * s0 +
            mu_s1[phase] * s1;
        for (std::size_t component = 0U;
             component < x.size();
             ++component) {
            value +=
                mu_x[phase][component] *
                x[component];
        }
        return value;
    }
};

const std::array<std::vector<double>, 3>&
phase_compositions() {
    static const std::array<std::vector<double>, 3>
        values{
            std::vector<double>{0.10, 0.70, 0.20},
            std::vector<double>{0.60, 0.20, 0.20},
            std::vector<double>{0.20, 0.30, 0.50}};
    return values;
}

std::vector<double> independent_values(
    std::span<const double> full,
    std::size_t dependent) {
    std::vector<double> result;
    result.reserve(full.size() - 1U);
    for (std::size_t component = 0U;
         component < full.size();
         ++component) {
        if (component != dependent) {
            result.push_back(full[component]);
        }
    }
    return result;
}

fl::NaturalVariableCellState3P make_state(
    const SyntheticTransportLaw3P& law) {
    constexpr double pressure = 10.0e6;
    constexpr double temperature = 350.0;
    constexpr double s0 = 0.20;
    constexpr double s1 = 0.30;

    const auto& compositions =
        phase_compositions();
    const auto pivot =
        fl::NaturalVariableCompositionPivot3P::
            select(compositions);

    fl::NaturalVariableCellStateInput3P input;
    input.component_ids = {"A", "B", "C"};
    input.reference_pressure_pa = pressure;
    input.temperature_k = temperature;
    input.independent_saturations = {s0, s1};
    input.composition_pivot = pivot;

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        input.independent_phase_compositions[phase] =
            independent_values(
                compositions[phase],
                pivot.dependent_component(
                    static_cast<fl::PhaseSlot3>(
                        phase)));
        const double rho =
            law.rho(
                phase,
                pressure,
                temperature,
                s0,
                s1,
                compositions[phase]);
        const double mu =
            law.mu(
                phase,
                pressure,
                temperature,
                s0,
                s1,
                compositions[phase]);

        input.phase_properties[phase] = {
            3000.0 +
                100.0 *
                    static_cast<double>(phase),
            rho,
            mu,
            2.0e5,
            1.5e5};
    }

    return fl::NaturalVariableCellState3P::create(
        std::move(input));
}

std::vector<double> natural_inputs(
    const fl::NaturalVariableCellState3P& state) {
    const auto& layout = state.layout();
    std::vector<double> inputs(
        layout.unknown_count(),
        0.0);
    inputs[layout.pressure_unknown_index()] =
        state.reference_pressure_pa();
    inputs[layout.temperature_unknown_index()] =
        state.temperature_k();
    inputs[*layout.independent_saturation_unknown_index(
        fl::PhaseSlot3::phase0)] =
        state.phase_saturation(
            fl::PhaseSlot3::phase0);
    inputs[*layout.independent_saturation_unknown_index(
        fl::PhaseSlot3::phase1)] =
        state.phase_saturation(
            fl::PhaseSlot3::phase1);

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<fl::PhaseSlot3>(
                phase);
        const auto x =
            state.phase_composition(slot);
        for (std::size_t component = 0U;
             component < x.size();
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        slot,
                        component);
            if (column) {
                inputs[*column] =
                    x[component];
            }
        }
    }
    return inputs;
}

std::array<std::vector<double>, 3>
reconstruct_compositions(
    const fl::NaturalVariableLayout3P& layout,
    std::span<const double> inputs) {
    std::array<std::vector<double>, 3> result;
    const std::size_t n =
        layout.component_count();

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto slot =
            static_cast<fl::PhaseSlot3>(
                phase);
        const std::size_t dependent =
            layout.dependent_composition_component(
                slot);
        result[phase].assign(n, 0.0);
        double sum = 0.0;
        for (std::size_t component = 0U;
             component < n;
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        slot,
                        component);
            if (!column) {
                continue;
            }
            result[phase][component] =
                inputs[*column];
            sum += inputs[*column];
        }
        result[phase][dependent] =
            1.0 - sum;
    }
    return result;
}

struct FreshProperties {
    std::array<double, 3> pressure{};
    std::array<double, 3> rho{};
    std::array<double, 3> mu{};
    std::array<double, 3> kr{};
    std::array<double, 3> mobility{};
};

FreshProperties evaluate_fresh(
    const fl::NaturalVariableLayout3P& layout,
    std::span<const double> inputs,
    const SyntheticTransportLaw3P& law) {
    const auto compositions =
        reconstruct_compositions(
            layout,
            inputs);
    const double p =
        inputs[layout.pressure_unknown_index()];
    const double t =
        inputs[layout.temperature_unknown_index()];
    const double s0 =
        inputs[*layout.independent_saturation_unknown_index(
            fl::PhaseSlot3::phase0)];
    const double s1 =
        inputs[*layout.independent_saturation_unknown_index(
            fl::PhaseSlot3::phase1)];
    const double s2 =
        1.0 - s0 - s1;
    const std::array<double, 3>
        saturation{s0, s1, s2};

    FreshProperties result;
    result.pressure = {
        p,
        p + 2.0e5 * (s1 - s0),
        p + 3.0e5 * (s2 - s0)};
    result.kr = {
        s0 * s0,
        s1 * s1,
        s2 * s2};

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        result.rho[phase] =
            law.rho(
                phase,
                p,
                t,
                s0,
                s1,
                compositions[phase]);
        result.mu[phase] =
            law.mu(
                phase,
                p,
                t,
                s0,
                s1,
                compositions[phase]);
        result.mobility[phase] =
            result.kr[phase] /
            result.mu[phase];
    }
    return result;
}

std::array<std::vector<double>, 3>
transport_gradient(
    const fl::NaturalVariableCellState3P& state,
    const SyntheticTransportLaw3P& law,
    bool viscosity) {
    const auto& layout =
        state.layout();
    const std::size_t q =
        layout.unknown_count();
    std::array<std::vector<double>, 3>
        gradient;

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        gradient[phase].assign(q, 0.0);

        gradient[phase][
            layout.pressure_unknown_index()] =
            viscosity
                ? law.mu_p[phase]
                : law.rho_p[phase];
        gradient[phase][
            layout.temperature_unknown_index()] =
            viscosity
                ? law.mu_t[phase]
                : law.rho_t[phase];

        gradient[phase][
            *layout.independent_saturation_unknown_index(
                fl::PhaseSlot3::phase0)] =
            viscosity
                ? law.mu_s0[phase]
                : law.rho_s0[phase];
        gradient[phase][
            *layout.independent_saturation_unknown_index(
                fl::PhaseSlot3::phase1)] =
            viscosity
                ? law.mu_s1[phase]
                : law.rho_s1[phase];

        const auto slot =
            static_cast<fl::PhaseSlot3>(
                phase);
        const std::size_t dependent =
            layout
                .dependent_composition_component(
                    slot);
        for (std::size_t component = 0U;
             component <
             layout.component_count();
             ++component) {
            const auto column =
                layout
                    .independent_composition_unknown_index(
                        slot,
                        component);
            if (!column) {
                continue;
            }
            const double component_coeff =
                viscosity
                    ? law.mu_x[phase][component]
                    : law.rho_x[phase][component];
            const double dependent_coeff =
                viscosity
                    ? law.mu_x[phase][dependent]
                    : law.rho_x[phase][dependent];
            gradient[phase][*column] =
                component_coeff -
                dependent_coeff;
        }
    }
    return gradient;
}

fl::ThreePhaseSaturationCoordinateDerivatives3P
saturation_derivatives() {
    using D = ad::Dual<double, 2U>;
    const auto evaluated =
        fl::evaluate_three_phase_saturation_constitutive(
            D{10.0e6},
            D::variable(0.20, 0U),
            D::variable(0.30, 1U),
            SyntheticRelativePermeability3P{},
            SyntheticCapillaryPressure3P{});

    fl::ThreePhaseSaturationCoordinateDerivatives3P
        result{};
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        for (std::size_t direction = 0U;
             direction < 2U;
             ++direction) {
            result.relative_permeability
                [phase][direction] =
                evaluated
                    .relative_permeability[phase]
                    .derivative(direction);
            result.capillary_pressure_offset_pa
                [phase][direction] =
                evaluated
                    .capillary_pressure_offset_pa[
                        phase]
                    .derivative(direction);
        }
    }
    return result;
}

void local_mobility() {
    const SyntheticTransportLaw3P law;
    const auto state =
        make_state(law);

    const auto transport =
        fl::make_phase_transport_property_linearization(
            state,
            transport_gradient(
                state,
                law,
                false),
            transport_gradient(
                state,
                law,
                true),
            {"synthetic-rho",
             "structural-fixture",
             "v1"},
            {"synthetic-mu",
             "structural-fixture",
             "v1"});

    const auto constitutive_primal =
        fl::evaluate_three_phase_saturation_constitutive(
            state,
            SyntheticRelativePermeability3P{},
            SyntheticCapillaryPressure3P{});
    const auto constitutive =
        fl::make_saturation_constitutive_natural_variable_linearization(
            state,
            constitutive_primal,
            saturation_derivatives());

    const auto mobility =
        fl::build_local_phase_mobility_linearization(
            state,
            transport,
            constitutive);

    require(
        mobility.state_identity.component_ids ==
                std::vector<std::string>{"A", "B", "C"} &&
            mobility.state_identity.layout
                    .composition_pivot()
                    .dependent_components() ==
                state.layout()
                    .composition_pivot()
                    .dependent_components() &&
            mobility.mass_density_provenance.model ==
                "synthetic-rho" &&
            mobility.viscosity_provenance.model ==
                "synthetic-mu",
        "local mobility lost exact-state identity/provenance");

    const auto inputs =
        natural_inputs(state);
    const auto fresh =
        evaluate_fresh(
            state.layout(),
            inputs,
            law);

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        near(
            mobility.phase_pressure_pa[phase],
            fresh.pressure[phase],
            1.0e-13,
            1.0e-8);
        near(
            mobility.mass_density_kg_per_m3[phase],
            fresh.rho[phase],
            1.0e-13,
            1.0e-10);
        near(
            mobility.dynamic_viscosity_pa_s[phase],
            fresh.mu[phase],
            1.0e-13,
            1.0e-14);
        near(
            mobility.relative_permeability[phase],
            fresh.kr[phase],
            1.0e-13,
            1.0e-14);
        near(
            mobility.mobility_per_pa_s[phase],
            fresh.mobility[phase],
            1.0e-13,
            1.0e-10);
    }

    for (std::size_t column = 0U;
         column < inputs.size();
         ++column) {
        double step = 1.0e-6;
        if (column ==
            state.layout().pressure_unknown_index()) {
            step = 10.0;
        } else if (
            column ==
            state.layout().temperature_unknown_index()) {
            step = 1.0e-4;
        }

        auto plus = inputs;
        auto minus = inputs;
        plus[column] += step;
        minus[column] -= step;

        const auto plus_value =
            evaluate_fresh(
                state.layout(),
                plus,
                law);
        const auto minus_value =
            evaluate_fresh(
                state.layout(),
                minus,
                law);

        for (std::size_t phase = 0U;
             phase < 3U;
             ++phase) {
            const auto fd =
                [&](double first,
                    double second) {
                    return (first - second) /
                        (2.0 * step);
                };

            near(
                mobility.phase_pressure_gradient
                    [phase][column],
                fd(
                    plus_value.pressure[phase],
                    minus_value.pressure[phase]),
                2.0e-7,
                2.0e-5);
            near(
                mobility.mass_density_gradient
                    [phase][column],
                fd(
                    plus_value.rho[phase],
                    minus_value.rho[phase]),
                2.0e-6,
                2.0e-6);
            near(
                mobility.dynamic_viscosity_gradient
                    [phase][column],
                fd(
                    plus_value.mu[phase],
                    minus_value.mu[phase]),
                2.0e-6,
                2.0e-10);
            near(
                mobility.relative_permeability_gradient
                    [phase][column],
                fd(
                    plus_value.kr[phase],
                    minus_value.kr[phase]),
                2.0e-6,
                2.0e-8);
            near(
                mobility.mobility_gradient
                    [phase][column],
                fd(
                    plus_value.mobility[phase],
                    minus_value.mobility[phase]),
                5.0e-6,
                5.0e-5);
        }
    }
}

void no_capillary_pressure_gradient() {
    const SyntheticTransportLaw3P law;
    const auto state =
        make_state(law);
    const auto primal =
        fl::evaluate_three_phase_saturation_constitutive(
            state,
            SyntheticRelativePermeability3P{},
            fl::NoCapillaryPressure3P{});

    using D = ad::Dual<double, 2U>;
    const auto dual =
        fl::evaluate_three_phase_saturation_constitutive(
            D{state.reference_pressure_pa()},
            D::variable(
                state.phase_saturation(
                    fl::PhaseSlot3::phase0),
                0U),
            D::variable(
                state.phase_saturation(
                    fl::PhaseSlot3::phase1),
                1U),
            SyntheticRelativePermeability3P{},
            fl::NoCapillaryPressure3P{});

    fl::ThreePhaseSaturationCoordinateDerivatives3P
        derivatives{};
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        for (std::size_t direction = 0U;
             direction < 2U;
             ++direction) {
            derivatives.relative_permeability
                [phase][direction] =
                dual.relative_permeability[phase]
                    .derivative(direction);
            derivatives.capillary_pressure_offset_pa
                [phase][direction] =
                dual.capillary_pressure_offset_pa[phase]
                    .derivative(direction);
        }
    }

    const auto linearized =
        fl::make_saturation_constitutive_natural_variable_linearization(
            state,
            primal,
            derivatives);

    const auto& layout =
        state.layout();
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        near(
            linearized.phase_pressure_gradient
                [phase][layout.pressure_unknown_index()],
            1.0,
            0.0,
            0.0);
        for (std::size_t column = 1U;
             column < layout.unknown_count();
             ++column) {
            near(
                linearized.phase_pressure_gradient
                    [phase][column],
                0.0,
                0.0,
                0.0);
        }
    }
}

void invalid_inputs() {
    const SyntheticTransportLaw3P law;
    const auto state =
        make_state(law);
    const auto rho_gradient =
        transport_gradient(
            state,
            law,
            false);
    const auto mu_gradient =
        transport_gradient(
            state,
            law,
            true);

    expect_invalid(
        [&] {
            (void)fl::
                make_phase_transport_property_linearization(
                    state,
                    rho_gradient,
                    mu_gradient,
                    {"",
                     "structural-fixture",
                     "v1"},
                    {"synthetic-mu",
                     "structural-fixture",
                     "v1"});
        },
        "provenance");

    auto nonfinite = mu_gradient;
    nonfinite[1][0] =
        std::numeric_limits<double>::infinity();
    expect_invalid(
        [&] {
            (void)fl::
                make_phase_transport_property_linearization(
                    state,
                    rho_gradient,
                    nonfinite,
                    {"synthetic-rho",
                     "structural-fixture",
                     "v1"},
                    {"synthetic-mu",
                     "structural-fixture",
                     "v1"});
        },
        "non-finite");

    const auto transport =
        fl::make_phase_transport_property_linearization(
            state,
            rho_gradient,
            mu_gradient,
            {"synthetic-rho",
             "structural-fixture",
             "v1"},
            {"synthetic-mu",
             "structural-fixture",
             "v1"});

    const auto primal =
        fl::evaluate_three_phase_saturation_constitutive(
            state,
            SyntheticRelativePermeability3P{},
            SyntheticCapillaryPressure3P{});

    auto bad_derivatives =
        saturation_derivatives();
    bad_derivatives
        .capillary_pressure_offset_pa[0][0] =
        1.0;
    expect_invalid(
        [&] {
            (void)fl::
                make_saturation_constitutive_natural_variable_linearization(
                    state,
                    primal,
                    bad_derivatives);
        },
        "reference-phase");

    auto inconsistent_primal =
        primal;
    inconsistent_primal.phase_pressure_pa[1] +=
        1000.0;
    expect_invalid(
        [&] {
            (void)fl::
                make_saturation_constitutive_natural_variable_linearization(
                    state,
                    inconsistent_primal,
                    saturation_derivatives());
        },
        "inconsistent");

    const auto constitutive =
        fl::make_saturation_constitutive_natural_variable_linearization(
            state,
            primal,
            saturation_derivatives());

    auto wrong_transport =
        transport;
    wrong_transport.dynamic_viscosity_pa_s[0] *=
        2.0;
    expect_invalid(
        [&] {
            (void)fl::
                build_local_phase_mobility_linearization(
                    state,
                    wrong_transport,
                    constitutive);
        },
        "viscosity primal");

    auto wrong_identity =
        transport;
    wrong_identity.state_identity.temperature_k +=
        1.0;
    expect_invalid(
        [&] {
            (void)fl::
                build_local_phase_mobility_linearization(
                    state,
                    wrong_identity,
                    constitutive);
        },
        "state identity");
}

void headers() {
    require(
        phase_transport_header(),
        "phase transport header probe failed");
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"local_mobility", local_mobility},
    {"no_capillary_pressure_gradient", no_capillary_pressure_gradient},
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
