#include <mpmc/well_discretization/energy_rate.hpp>

#include <algorithm>
#include <cmath>
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

void require_energy(
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

void near_energy(
    double actual,
    double expected,
    double relative = 5.0e-12,
    double absolute = 5.0e-20,
    std::source_location where =
        std::source_location::current()) {
    const double scale =
        std::max(
            {1.0e-30,
             std::abs(actual),
             std::abs(expected)});
    require_energy(
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                absolute + relative * scale,
        "energy-rate numeric mismatch",
        where);
}

template <class Function>
void expect_energy_invalid(
    Function&& function) {
    bool caught = false;
    try {
        std::forward<Function>(function)();
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    require_energy(
        caught,
        "expected invalid_argument");
}

well::PeacemanWellIndex3D
energy_connection() {
    return well::make_peaceman_well_index_3d(
        {20.0, 10.0, 5.0},
        {4.0e-13, 1.0e-13, 8.0e-14},
        well::AxisAlignedWellDirection3D::z,
        0.10,
        0.0);
}

flow::LocalPhaseMobilityLinearization3P
energy_mobility() {
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
    result.phase_pressure_pa =
        {1.20e6, 1.10e6, 0.90e6};
    result.mobility_per_pa_s =
        {1.0e-9, 2.0e-9, 3.0e-9};

    for (std::size_t phase = 0U;
         phase < flow::fixed_three_phase_count;
         ++phase) {
        result.phase_pressure_gradient[phase]
            .assign(q, 0.0);
        result.phase_pressure_gradient[phase][0] =
            1.0;
        result.phase_pressure_gradient[phase][2] =
            100.0 *
            static_cast<double>(
                phase + 1U);
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

flow::PhaseTransportPropertyNaturalVariableLinearization3P
energy_transport(
    const flow::NaturalVariableStateIdentity3P&
        identity) {
    flow::PhaseTransportPropertyNaturalVariableLinearization3P
        result{};
    result.state_identity = identity;
    result.mass_density_provenance =
        {"test-density", "energy-rate", "v1"};
    result.mass_density_kg_per_m3 =
        {700.0, 800.0, 900.0};

    const std::size_t q =
        identity.layout.unknown_count();
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        result.mass_density_gradient[phase]
            .resize(q);
        for (std::size_t column = 0U;
             column < q;
             ++column) {
            result.mass_density_gradient[phase][column] =
                0.01 *
                static_cast<double>(
                    (phase + 1U) *
                    (column + 1U));
        }
    }
    return result;
}

flow::PhaseCaloricPropertyNaturalVariableLinearization3P
energy_caloric(
    const flow::NaturalVariableStateIdentity3P&
        identity) {
    flow::PhaseCaloricPropertyNaturalVariableLinearization3P
        result{};
    result.state_identity = identity;
    result.enthalpy_provenance =
        {"test-enthalpy", "energy-rate", "v1"};
    result.specific_enthalpy_j_per_kg =
        {1000.0, 2000.0, 3000.0};

    const std::size_t q =
        identity.layout.unknown_count();
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        result.specific_enthalpy_gradient[phase]
            .resize(q);
        for (std::size_t column = 0U;
             column < q;
             ++column) {
            result.specific_enthalpy_gradient[phase][column] =
                0.5 *
                static_cast<double>(
                    (phase + 1U) *
                    (column + 1U));
        }
    }
    return result;
}

wd::InjectionPhaseSpecificEnthalpy3P
injection_enthalpy() {
    return {
        "test/injection-enthalpy/v1",
        {7000.0, 8000.0, 9000.0}};
}

void production_and_injection_enthalpy_selection() {
    const auto wi =
        energy_connection();
    const auto mobility =
        energy_mobility();
    const auto transport =
        energy_transport(
            mobility.state_identity);
    const auto caloric =
        energy_caloric(
            mobility.state_identity);
    const auto injection =
        injection_enthalpy();

    const auto result =
        wd::make_connection_advective_energy_rate_linearization_3p(
            wi,
            mobility,
            transport,
            caloric,
            injection,
            1.00e6);

    require_energy(
        result.phase_volumetric_rate_m3_per_s[0] >
                0.0 &&
            result.phase_volumetric_rate_m3_per_s[1] >
                0.0 &&
            result.phase_volumetric_rate_m3_per_s[2] <
                0.0,
        "mixed production/injection signs were not preserved");
    require_energy(
        result.enthalpy_source[0] ==
                wd::WellConnectionEnergyEnthalpySource3P::
                    reservoir_production_or_zero_tie &&
            result.enthalpy_source[1] ==
                wd::WellConnectionEnergyEnthalpySource3P::
                    reservoir_production_or_zero_tie &&
            result.enthalpy_source[2] ==
                wd::WellConnectionEnergyEnthalpySource3P::
                    explicit_well_injection,
        "energy enthalpy direction selection mismatch");

    near_energy(
        result.selected_specific_enthalpy_j_per_kg[0],
        1000.0);
    near_energy(
        result.selected_specific_enthalpy_j_per_kg[1],
        2000.0);
    near_energy(
        result.selected_specific_enthalpy_j_per_kg[2],
        9000.0);

    double expected = 0.0;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        expected +=
            result.phase_volumetric_rate_m3_per_s[
                phase] *
            transport.mass_density_kg_per_m3[
                phase] *
            result.selected_specific_enthalpy_j_per_kg[
                phase];
    }
    near_energy(
        result.total_advective_energy_rate_w,
        expected);
}

void production_reservoir_enthalpy_product_rule() {
    const auto wi =
        energy_connection();
    const auto mobility =
        energy_mobility();
    const auto transport =
        energy_transport(
            mobility.state_identity);
    const auto caloric =
        energy_caloric(
            mobility.state_identity);
    const auto result =
        wd::make_connection_advective_energy_rate_linearization_3p(
            wi,
            mobility,
            transport,
            caloric,
            injection_enthalpy(),
            1.00e6);

    constexpr std::size_t column = 2U;
    double expected = 0.0;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto rate =
            wd::make_pressure_drawdown_phase_rate_linearization_3p(
                wi,
                mobility,
                static_cast<flow::PhaseSlot3>(
                    phase),
                1.00e6);
        const double q =
            rate.phase_volumetric_rate_m3_per_s;
        const double rho =
            transport.mass_density_kg_per_m3[
                phase];
        const double h =
            q >= 0.0
                ? caloric.specific_enthalpy_j_per_kg[
                      phase]
                : injection_enthalpy()
                      .specific_enthalpy_j_per_kg[
                          phase];
        const double dh =
            q >= 0.0
                ? caloric.specific_enthalpy_gradient[
                      phase][column]
                : 0.0;
        expected +=
            rate.d_phase_rate_d_reservoir_unknown(
                column) *
                rho *
                h +
            q *
                transport.mass_density_gradient[
                    phase][column] *
                h +
            q *
                rho *
                dh;
    }
    near_energy(
        result.d_energy_rate_d_reservoir_unknown(
            column),
        expected);
}

void injection_does_not_use_reservoir_enthalpy_derivative() {
    const auto wi =
        energy_connection();
    auto mobility =
        energy_mobility();
    const auto transport =
        energy_transport(
            mobility.state_identity);
    auto caloric =
        energy_caloric(
            mobility.state_identity);

    // Isolate phase2 injection. Give its reservoir enthalpy derivative a huge
    // value: the injection branch must ignore it.
    mobility.mobility_per_pa_s[0] =
        0.0;
    mobility.mobility_per_pa_s[1] =
        0.0;
    for (std::size_t column = 0U;
         column < 7U;
         ++column) {
        mobility.mobility_gradient[0][column] =
            0.0;
        mobility.mobility_gradient[1][column] =
            0.0;
    }
    caloric.specific_enthalpy_gradient[2][0] =
        1.0e30;

    const auto injection =
        injection_enthalpy();
    const auto result =
        wd::make_connection_advective_energy_rate_linearization_3p(
            wi,
            mobility,
            transport,
            caloric,
            injection,
            1.00e6);

    const auto rate =
        wd::make_pressure_drawdown_phase_rate_linearization_3p(
            wi,
            mobility,
            flow::PhaseSlot3::phase2,
            1.00e6);
    const double q =
        rate.phase_volumetric_rate_m3_per_s;
    require_energy(
        q < 0.0,
        "phase2 fixture must be injection");

    const double h_inj =
        injection.specific_enthalpy_j_per_kg[2];
    const double expected =
        rate.d_phase_rate_d_reservoir_unknown(
            0U) *
            transport.mass_density_kg_per_m3[2] *
            h_inj +
        q *
            transport.mass_density_gradient[2][0] *
            h_inj;
    near_energy(
        result.d_energy_rate_d_reservoir_unknown(
            0U),
        expected);
}

void bhp_derivative_uses_selected_enthalpy() {
    const auto wi =
        energy_connection();
    const auto mobility =
        energy_mobility();
    const auto transport =
        energy_transport(
            mobility.state_identity);
    const auto caloric =
        energy_caloric(
            mobility.state_identity);
    const auto injection =
        injection_enthalpy();

    const auto result =
        wd::make_connection_advective_energy_rate_linearization_3p(
            wi,
            mobility,
            transport,
            caloric,
            injection,
            1.00e6);

    double expected = 0.0;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto rate =
            wd::make_pressure_drawdown_phase_rate_linearization_3p(
                wi,
                mobility,
                static_cast<flow::PhaseSlot3>(
                    phase),
                1.00e6);
        expected +=
            rate.d_phase_rate_d_bottom_hole_pressure() *
            transport.mass_density_kg_per_m3[
                phase] *
            result.selected_specific_enthalpy_j_per_kg[
                phase];
    }
    near_energy(
        result.d_energy_rate_d_bottom_hole_pressure(),
        expected);
}

void zero_rate_tie_uses_reservoir_enthalpy() {
    const auto wi =
        energy_connection();
    auto mobility =
        energy_mobility();
    const auto transport =
        energy_transport(
            mobility.state_identity);
    const auto caloric =
        energy_caloric(
            mobility.state_identity);
    auto injection =
        injection_enthalpy();

    mobility.mobility_per_pa_s[0] =
        0.0;
    mobility.mobility_per_pa_s[2] =
        0.0;
    for (std::size_t column = 0U;
         column < 7U;
         ++column) {
        mobility.mobility_gradient[0][column] =
            0.0;
        mobility.mobility_gradient[2][column] =
            0.0;
    }
    injection.specific_enthalpy_j_per_kg[1] =
        9.0e9;

    const double bhp =
        mobility.phase_pressure_pa[1];
    const auto result =
        wd::make_connection_advective_energy_rate_linearization_3p(
            wi,
            mobility,
            transport,
            caloric,
            injection,
            bhp);

    require_energy(
        result.phase_volumetric_rate_m3_per_s[1] ==
                0.0 &&
            result.phase_advective_energy_rate_w[1] ==
                0.0 &&
            result.enthalpy_source[1] ==
                wd::WellConnectionEnergyEnthalpySource3P::
                    reservoir_production_or_zero_tie,
        "exact-zero rate did not use deterministic reservoir enthalpy tie branch");
    near_energy(
        result.selected_specific_enthalpy_j_per_kg[1],
        caloric.specific_enthalpy_j_per_kg[1]);

    const auto rate =
        wd::make_pressure_drawdown_phase_rate_linearization_3p(
            wi,
            mobility,
            flow::PhaseSlot3::phase1,
            bhp);
    const double expected =
        rate.d_phase_rate_d_reservoir_unknown(
            0U) *
        transport.mass_density_kg_per_m3[1] *
        caloric.specific_enthalpy_j_per_kg[1];
    near_energy(
        result.d_energy_rate_d_reservoir_unknown(
            0U),
        expected);
}

void energy_invalid_inputs() {
    const auto wi =
        energy_connection();
    const auto mobility =
        energy_mobility();

    {
        auto injection =
            injection_enthalpy();
        injection.provenance.clear();
        expect_energy_invalid(
            [&] {
                (void)wd::
                    make_connection_advective_energy_rate_linearization_3p(
                        wi,
                        mobility,
                        energy_transport(
                            mobility.state_identity),
                        energy_caloric(
                            mobility.state_identity),
                        injection,
                        1.0e6);
            });
    }
    {
        auto injection =
            injection_enthalpy();
        injection.specific_enthalpy_j_per_kg[2] =
            std::numeric_limits<double>::
                quiet_NaN();
        expect_energy_invalid(
            [&] {
                (void)wd::
                    make_connection_advective_energy_rate_linearization_3p(
                        wi,
                        mobility,
                        energy_transport(
                            mobility.state_identity),
                        energy_caloric(
                            mobility.state_identity),
                        injection,
                        1.0e6);
            });
    }
    {
        auto transport =
            energy_transport(
                mobility.state_identity);
        transport.mass_density_gradient[0]
            .pop_back();
        expect_energy_invalid(
            [&] {
                (void)wd::
                    make_connection_advective_energy_rate_linearization_3p(
                        wi,
                        mobility,
                        transport,
                        energy_caloric(
                            mobility.state_identity),
                        injection_enthalpy(),
                        1.0e6);
            });
    }
    {
        auto caloric =
            energy_caloric(
                mobility.state_identity);
        caloric.specific_enthalpy_gradient[0][0] =
            std::numeric_limits<double>::
                quiet_NaN();
        expect_energy_invalid(
            [&] {
                (void)wd::
                    make_connection_advective_energy_rate_linearization_3p(
                        wi,
                        mobility,
                        energy_transport(
                            mobility.state_identity),
                        caloric,
                        injection_enthalpy(),
                        1.0e6);
            });
    }
    {
        auto transport =
            energy_transport(
                mobility.state_identity);
        transport.state_identity.temperature_k +=
            1.0;
        expect_energy_invalid(
            [&] {
                (void)wd::
                    make_connection_advective_energy_rate_linearization_3p(
                        wi,
                        mobility,
                        transport,
                        energy_caloric(
                            mobility.state_identity),
                        injection_enthalpy(),
                        1.0e6);
            });
    }
}

} // namespace

bool run_energy_rate_case(
    std::string_view name) {
    if (name ==
        "energy_rate_production_and_injection_enthalpy_selection") {
        production_and_injection_enthalpy_selection();
    } else if (
        name ==
        "energy_rate_production_reservoir_enthalpy_product_rule") {
        production_reservoir_enthalpy_product_rule();
    } else if (
        name ==
        "energy_rate_injection_ignores_reservoir_enthalpy_derivative") {
        injection_does_not_use_reservoir_enthalpy_derivative();
    } else if (
        name ==
        "energy_rate_bhp_derivative_selected_enthalpy") {
        bhp_derivative_uses_selected_enthalpy();
    } else if (
        name ==
        "energy_rate_zero_tie_reservoir_enthalpy") {
        zero_rate_tie_uses_reservoir_enthalpy();
    } else if (
        name ==
        "energy_rate_invalid_inputs") {
        energy_invalid_inputs();
    } else {
        return false;
    }
    return true;
}
