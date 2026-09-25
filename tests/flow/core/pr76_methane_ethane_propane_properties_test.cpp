#include <mpmc/flow/pr76_methane_ethane_propane_properties.hpp>

#include <mpmc/ad/dual.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace ad = mpmc::ad;
namespace flow = mpmc::flow;
namespace th = mpmc::thermodynamics;

void require_real(
    bool condition,
    std::string_view message) {
    if (!condition) {
        throw std::runtime_error(
            std::string{message});
    }
}

void near_real(
    double actual,
    double expected,
    double relative = 2.0e-9,
    double absolute = 2.0e-11) {
    const double scale =
        std::max(
            {1.0,
             std::abs(actual),
             std::abs(expected)});
    require_real(
        std::isfinite(actual) &&
            std::isfinite(expected) &&
            std::abs(actual - expected) <=
                absolute +
                    relative * scale,
        "real PR76 property numeric reference changed");
}

th::Provenance source_real(
    std::string locator) {
    return {
        th::SourceKind::synthetic_test,
        "test://pr76-real-property-provider",
        "v1",
        std::move(locator),
        "Software regression fixture using the production curated numeric identity",
        "Constructed in the flow-core test",
        "No third-party source text redistributed"};
}

th::SourcedScalar datum_real(
    double value,
    th::Unit unit,
    std::string locator) {
    return {
        value,
        unit,
        source_real(
            std::move(locator)),
        "specified SI",
        "identity"};
}

[[nodiscard]] th::PrParameterSet
curated_compatible_parameters(
    std::string dataset_id =
        "DeitersBell-aic16730-PengRobinson1976-ternary") {
    constexpr std::array<std::string_view, 3>
        ids{
            "methane",
            "ethane",
            "propane"};
    constexpr std::array<double, 3>
        tc{
            190.555,
            305.4,
            369.825};
    constexpr std::array<double, 3>
        pc{
            4.595e6,
            4.88e6,
            4.248e6};
    constexpr std::array<double, 3>
        omega{
            0.0,
            0.099,
            0.15308};
    constexpr std::array<double, 3>
        mw{
            0.0160425,
            0.0300690,
            0.0440956};

    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    input.model_id =
        std::string{th::pr76_profile};
    input.dataset_id =
        std::move(dataset_id);
    input.revision =
        "curated-compatible-test-v1";
    input.applicability = {
        std::nullopt,
        std::nullopt,
        source_real(
            "unknown applicability")};

    for (std::size_t i = 0U;
         i < ids.size();
         ++i) {
        const std::string id{ids[i]};
        catalog.push_back({
            id,
            id,
            th::ComponentKind::pure,
            source_real(
                "component-" +
                std::to_string(i)),
            datum_real(
                mw[i],
                th::Unit::kilogram_per_mole,
                "molar-mass-" +
                    std::to_string(i))});
        input.pure.push_back({
            id,
            datum_real(
                tc[i],
                th::Unit::kelvin,
                "Tc-" +
                    std::to_string(i)),
            datum_real(
                pc[i],
                th::Unit::pascal,
                "Pc-" +
                    std::to_string(i)),
            datum_real(
                omega[i],
                th::Unit::dimensionless,
                "omega-" +
                    std::to_string(i))});
    }

    const auto zero =
        datum_real(
            0.0,
            th::Unit::dimensionless,
            "zero-kij");
    input.binary.push_back(
        {"methane", "ethane", zero});
    input.binary.push_back(
        {"methane", "propane", zero});
    input.binary.push_back(
        {"ethane", "propane", zero});

    const std::array<std::string, 3> order{
        "methane",
        "ethane",
        "propane"};
    return th::PrParameterSet::create(
        catalog,
        order,
        input,
        th::DataPolicy::allow_synthetic_tests);
}

template <typename Function>
void expect_invalid_real(Function&& function) {
    bool caught = false;
    try {
        std::forward<Function>(
            function)();
    } catch (const std::invalid_argument&) {
        caught = true;
    } catch (const std::out_of_range&) {
        caught = true;
    }
    require_real(
        caught,
        "expected real-property provider error was not reported");
}

void nist_ideal_enthalpy_contract(
    const flow::
        Pr76MethaneEthanePropanePropertyProvider&
            provider) {
    constexpr std::array<double, 3>
        composition{
            0.25,
            0.50,
            0.25};

    near_real(
        provider
            .ideal_gas_molar_enthalpy_j_per_mol(
                298.15,
                std::span<const double>{
                    composition}),
        0.0,
        0.0,
        1.0e-12);

    // Independent trapezoidal integration of the published NIST SRD 69
    // recommended Cp knots from 298.15 K to 450 K.
    constexpr double expected_450 =
        9581.981625;
    near_real(
        provider
            .ideal_gas_molar_enthalpy_j_per_mol(
                450.0,
                std::span<const double>{
                    composition}),
        expected_450);

    using D = ad::Dual<double, 1>;
    const std::array<D, 3> x{
        D{composition[0]},
        D{composition[1]},
        D{composition[2]}};
    const auto active =
        provider
            .ideal_gas_molar_enthalpy_j_per_mol(
                D::variable(
                    450.0,
                    0U),
                std::span<const D>{x});

    // At 450 K piecewise-linear Cp gives the midpoint of each 400/500 K pair.
    constexpr double expected_cp_450 =
        72.5825;
    near_real(
        active.derivative(0U),
        expected_cp_450);

    expect_invalid_real([&] {
        (void)provider
            .ideal_gas_molar_enthalpy_j_per_mol(
                298.14,
                std::span<const double>{
                    composition});
    });
    expect_invalid_real([&] {
        (void)provider
            .ideal_gas_molar_enthalpy_j_per_mol(
                1500.01,
                std::span<const double>{
                    composition});
    });
}

void lbc_contract(
    const flow::
        Pr76MethaneEthanePropanePropertyProvider&
            provider) {
    constexpr std::array<double, 3>
        composition{
            0.25,
            0.50,
            0.25};
    constexpr double mixture_molar_mass =
        0.030069025;
    constexpr double molar_density =
        5000.0;

    // Independent field-unit LBC/Stiel-Thodos/Herning-Zipperer evaluation
    // using Tc/Pc from the curated PR snapshot and NIST critical volumes.
    constexpr double expected_pa_s =
        2.0590343064952366e-5;
    const double actual =
        provider.lbc_dynamic_viscosity_pa_s(
            450.0,
            std::span<const double>{
                composition},
            mixture_molar_mass,
            molar_density);
    near_real(
        actual,
        expected_pa_s,
        2.0e-10,
        1.0e-14);
}

void residual_enthalpy_contract(
    const th::Pr76Phase<double>& model,
    const flow::
        Pr76MethaneEthanePropanePropertyProvider&
            provider) {
    constexpr std::array<double, 3>
        composition{
            0.25,
            0.50,
            0.25};
    th::Pr76PhaseWorkspace<double>
        workspace;
    const auto roots =
        model.roots_full(
            1.0e6,
            450.0,
            std::span<const double>{
                composition},
            workspace);
    require_real(
        roots.status ==
                th::Pr76RootStatus::success &&
            roots.count == 1U,
        "real PR76 enthalpy fixture root topology changed");

    const th::Pr76SelectedPhase
        selection{0U, {}};
    const double residual =
        provider
            .pr76_residual_molar_enthalpy_j_per_mol(
                selection,
                1.0e6,
                450.0,
                std::span<const double>{
                    composition});

    // Independent direct-Z PR departure evaluation at this fixture.
    constexpr double expected_residual =
        -296.87363274649726;
    near_real(
        residual,
        expected_residual,
        5.0e-9,
        1.0e-8);
}

void closure_ad_contract(
    const th::Pr76Phase<double>& model) {
    constexpr std::array<double, 3>
        composition{
            0.25,
            0.50,
            0.25};
    constexpr double mixture_molar_mass =
        0.030069025;

    auto closure =
        flow::
            make_pr76_methane_ethane_propane_property_closure(
                model,
                std::vector<
                    th::Pr76SelectedPhase>{
                    {0U, {}}});

    using D = ad::Dual<double, 3>;
    const std::array<D, 3> active_x{
        D{
            composition[0],
            {0.0, 0.0, 1.0}},
        D{composition[1]},
        D{
            composition[2],
            {0.0, 0.0, -1.0}}};

    const auto active =
        closure.evaluate(
            0U,
            D::variable(
                1.0e6,
                0U),
            D::variable(
                450.0,
                1U),
            std::span<const D>{
                active_x});

    near_real(
        active
            .mixture_molar_mass_kg_per_mol
            .value(),
        mixture_molar_mass);
    near_real(
        active
            .specific_enthalpy_j_per_kg
            .value(),
        308793.1182422278,
        5.0e-9,
        1.0e-6);

    const std::array<double, 3> steps{
        10.0,
        1.0e-3,
        1.0e-6};

    for (std::size_t lane = 0U;
         lane < steps.size();
         ++lane) {
        double plus_pressure = 1.0e6;
        double minus_pressure = 1.0e6;
        double plus_temperature = 450.0;
        double minus_temperature = 450.0;
        auto plus_x = composition;
        auto minus_x = composition;

        if (lane == 0U) {
            plus_pressure += steps[lane];
            minus_pressure -= steps[lane];
        } else if (lane == 1U) {
            plus_temperature += steps[lane];
            minus_temperature -= steps[lane];
        } else {
            plus_x[0] += steps[lane];
            plus_x[2] -= steps[lane];
            minus_x[0] -= steps[lane];
            minus_x[2] += steps[lane];
        }

        const auto plus =
            closure.evaluate(
                0U,
                plus_pressure,
                plus_temperature,
                std::span<const double>{
                    plus_x});
        const auto minus =
            closure.evaluate(
                0U,
                minus_pressure,
                minus_temperature,
                std::span<const double>{
                    minus_x});

        const auto central =
            [&](double positive,
                double negative) {
                return
                    (positive - negative) /
                    (2.0 * steps[lane]);
            };
        const auto compare =
            [&](double automatic,
                double finite_difference,
                const char* message) {
                const double scale =
                    std::max(
                        {1.0,
                         std::abs(automatic),
                         std::abs(
                             finite_difference)});
                require_real(
                    std::isfinite(automatic) &&
                        std::isfinite(
                            finite_difference) &&
                        std::abs(
                            automatic -
                            finite_difference) <=
                            4.0e-4 * scale,
                    message);
            };

        compare(
            active
                .dynamic_viscosity_pa_s
                .derivative(lane),
            central(
                plus
                    .dynamic_viscosity_pa_s,
                minus
                    .dynamic_viscosity_pa_s),
            "LBC AD derivative disagrees with fresh perturbation");

        compare(
            active
                .specific_enthalpy_j_per_kg
                .derivative(lane),
            central(
                plus
                    .specific_enthalpy_j_per_kg,
                minus
                    .specific_enthalpy_j_per_kg),
            "PR/NIST enthalpy AD derivative disagrees with fresh perturbation");

        compare(
            active
                .specific_internal_energy_j_per_kg
                .derivative(lane),
            central(
                plus
                    .specific_internal_energy_j_per_kg,
                minus
                    .specific_internal_energy_j_per_kg),
            "PR/NIST internal-energy AD derivative disagrees with fresh perturbation");
    }
}

} // namespace

void pr76_methane_ethane_propane_provider_regression() {
    auto parameters =
        curated_compatible_parameters();
    auto model =
        th::Pr76Phase<double>::
            from_parameters(parameters);
    const flow::
        Pr76MethaneEthanePropanePropertyProvider
            provider{model};

    nist_ideal_enthalpy_contract(
        provider);
    lbc_contract(
        provider);
    residual_enthalpy_contract(
        model,
        provider);
    closure_ad_contract(
        model);

    expect_invalid_real([&] {
        auto wrong_parameters =
            curated_compatible_parameters(
                "not-the-curated-dataset");
        auto wrong_model =
            th::Pr76Phase<double>::
                from_parameters(
                    wrong_parameters);
        (void)flow::
            Pr76MethaneEthanePropanePropertyProvider{
                wrong_model};
    });
}
