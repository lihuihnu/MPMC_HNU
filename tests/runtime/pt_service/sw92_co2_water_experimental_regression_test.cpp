#include <mpmc/flash/sw92_profile_c_pt_flash_backend.hpp>
#include <mpmc/runtime/pt_service.hpp>

#include "../../thermodynamics/sw92/test_support.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {
namespace fl = mpmc::flash;
namespace rt = mpmc::runtime;
namespace th = mpmc::thermodynamics;

// Independent experimental validation data, not production-model golden values.
// Original measurements:
//   A. Zawisza and B. Malesinska,
//   "Solubility of carbon dioxide in liquid water and of water in gaseous
//    carbon dioxide in the range 0.2-5 MPa and at temperatures up to 473 K",
//   J. Chem. Eng. Data 26 (1981) 388-391.
//   DOI: https://doi.org/10.1021/je00026a012
// Machine-readable compilation used for the exact rows below:
//   IUPAC-NIST Solubility Database, SRD 106, system 62_57,
//   https://srdata.nist.gov/solubility/sol_detail.aspx?sysID=62_57
// NIST identifies component 1 as CO2 and reports the tabulated x1 values below.
// The NIST page notes that these mole fractions are compiler-calculated from the
// original measurements. Therefore the acceptance budgets below are MODEL
// regression tolerances; they are not presented as experimental uncertainty.
struct ExperimentalPoint {
    double temperature_k;
    double pressure_bar;
    double aqueous_co2_mole_fraction;
};

constexpr std::array<ExperimentalPoint, 9> experimental_points{{
    {323.15, 8.92, 0.0030},
    {323.15, 14.41, 0.0046},
    {323.15, 25.05, 0.0076},
    {348.15, 2.37, 0.0005},
    {348.15, 21.00, 0.0046},
    {348.15, 35.91, 0.0076},
    {373.15, 3.56, 0.0005},
    {373.15, 27.10, 0.0046},
    {373.15, 45.60, 0.0076},
}};

// These are deliberately absolute mole-fraction budgets. They are wide enough
// to test the published corrected-original SW92 model rather than refit it to
// this independent 1981 dataset, but tight enough to catch coefficient/routing,
// unit-conversion, phase-selection, or order-of-magnitude regressions.
constexpr double max_absolute_error_budget = 0.0025;
constexpr double mean_absolute_error_budget = 0.0015;

void require(bool condition, std::string_view message) {
    if (!condition) { throw std::runtime_error(std::string(message)); }
}

double component_fraction(const rt::PtServicePhase& phase,
                          std::string_view component_id) {
    for (const auto& component : phase.components) {
        if (component.component_id == component_id) {
            return component.mole_fraction;
        }
    }
    throw std::runtime_error("experimental regression lost component ID");
}

double predicted_aqueous_co2(rt::PtService& service,
                             const ExperimentalPoint& point) {
    rt::PtServiceRequest request;
    request.configured_backend_id = "sw92-exp";
    request.pressure_pa = point.pressure_bar * 100000.0;
    request.temperature_k = point.temperature_k;
    // An equimolar overall feed is only a phase-split carrier. At accepted
    // two-phase equilibrium the water-richer phase composition is an intensive
    // equilibrium result and does not depend on assigning a public phase label.
    request.feed = {{sw92_test::water.id, 0.5},
                    {sw92_test::co2.id, 0.5}};

    const auto response = service.solve(request);
    require(response.structurally_valid(),
            "experimental regression returned structurally invalid service response");
    require(response.outcome == rt::PtServiceOutcome::accepted &&
                response.result.has_value(),
            "experimental regression state was not accepted");
    require(response.result->phases.size() == 2U,
            "CO2-water experimental anchor must close as two phases");

    const rt::PtServicePhase* water_richer = nullptr;
    double greatest_water_fraction = -std::numeric_limits<double>::infinity();
    for (const auto& phase : response.result->phases) {
        const double water_fraction =
            component_fraction(phase, sw92_test::water.id);
        if (water_fraction > greatest_water_fraction) {
            greatest_water_fraction = water_fraction;
            water_richer = &phase;
        }
    }
    require(water_richer != nullptr,
            "experimental regression could not identify water-richer composition coordinate");
    return component_fraction(*water_richer, sw92_test::co2.id);
}

} // namespace

int main() {
    try {
        const auto sw_model = th::Sw92Phase<double>::from_parameters(
            sw92_test::binary_parameters(sw92_test::co2));
        fl::Sw92ProfileCPtFlashBackend sw_backend(sw_model);
        const std::array<rt::PtServiceBackendRegistration, 1> registrations{{
            {"sw92-exp", &sw_backend}}};
        rt::PtService service(registrations);

        double max_absolute_error = 0.0;
        double absolute_error_sum = 0.0;
        std::cout << std::setprecision(17);
        for (const auto& point : experimental_points) {
            const double predicted = predicted_aqueous_co2(service, point);
            const double absolute_error =
                std::abs(predicted - point.aqueous_co2_mole_fraction);
            max_absolute_error = std::max(max_absolute_error, absolute_error);
            absolute_error_sum += absolute_error;
            std::cout << "SW92_CO2_WATER_EXPERIMENT"
                      << " t_k=" << point.temperature_k
                      << " p_bar=" << point.pressure_bar
                      << " x_co2_exp=" << point.aqueous_co2_mole_fraction
                      << " x_co2_sw92=" << predicted
                      << " abs_error=" << absolute_error << '\n';
        }
        const double mean_absolute_error =
            absolute_error_sum / static_cast<double>(experimental_points.size());
        std::cout << "SW92_CO2_WATER_EXPERIMENT_SUMMARY"
                  << " points=" << experimental_points.size()
                  << " max_abs_error=" << max_absolute_error
                  << " mean_abs_error=" << mean_absolute_error
                  << " max_budget=" << max_absolute_error_budget
                  << " mean_budget=" << mean_absolute_error_budget << '\n';

        require(max_absolute_error <= max_absolute_error_budget,
                "SW92 CO2-water maximum absolute experimental error exceeded budget");
        require(mean_absolute_error <= mean_absolute_error_budget,
                "SW92 CO2-water mean absolute experimental error exceeded budget");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SW92_CO2_WATER_EXPERIMENT_FAIL " << error.what() << '\n';
        return 1;
    }
}
