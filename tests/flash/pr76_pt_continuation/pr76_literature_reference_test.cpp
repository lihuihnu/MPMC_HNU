#include "literature_kpoint_reference.hpp"

#include <cmath>
#include <stdexcept>
#include <string_view>

namespace {
namespace ref = pr76_literature_reference;
namespace th = mpmc::thermodynamics;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

void near(double actual, double expected, double tolerance, const char* message) {
    require(std::isfinite(actual) && std::abs(actual - expected) <= tolerance,
            message);
}

void check_source(const th::Provenance& source) {
    require(source.kind == th::SourceKind::literature,
            "literature reference lost source kind");
    require(!source.reference.empty() && !source.revision.empty() &&
                !source.locator.empty() && !source.acquisition.empty() &&
                !source.usage_terms.empty(),
            "literature provenance is incomplete");
}

} // namespace

int main() {
    try {
        const auto anchor = ref::experimental_kpoint();
        require(anchor.experimental && anchor.critical_endpoint &&
                    !anchor.suitable_as_distinct_three_phase_flash_oracle,
                "experimental K-point semantics were weakened");
        check_source(anchor.source);
        near(anchor.temperature_k, 317.43, 0.0,
             "Smits K-point temperature changed");
        near(anchor.pressure_pa, 9.22e6, 0.0,
             "Smits K-point pressure changed");
        const auto feed = anchor.ordered_feed();
        near(feed[0], 0.9501, 1.0e-15,
             "Smits K-point CO2 loading changed");
        near(feed[1], 0.00697103, 1.0e-12,
             "derived 1-pentanol loading changed");
        near(feed[2], 0.04292897, 1.0e-12,
             "derived n-tridecane loading changed");
        near(feed[0] + feed[1] + feed[2], 1.0, 2.0e-15,
             "derived literature loading is not normalized");

        const auto model = ref::illustrative_pr76_model();
        const auto& parameters = model.parameters();
        require(parameters.dataset_id() ==
                    "Mushrif-2004-Fig3.4-PR76-illustrative-kpoint-surface" &&
                    parameters.revision() == "literature-transcription-v1",
                "Mushrif reference identity changed");
        require(parameters.components().size() == 3U &&
                    parameters.components().at(0).id == "CO2" &&
                    parameters.components().at(1).id == "1-pentanol" &&
                    parameters.components().at(2).id == "n-tridecane",
                "literature PR76 component order changed");

        const auto tc = parameters.critical_temperatures_k();
        const auto pc = parameters.critical_pressures_pa();
        const auto omega = parameters.acentric_factors();
        near(tc[0], 304.1, 0.0, "CO2 Tc changed");
        near(pc[0], 7.38e6, 0.0, "CO2 Pc changed");
        near(omega[0], 0.2390, 0.0, "CO2 omega changed");
        near(tc[1], 588.2, 0.0, "1-pentanol Tc changed");
        near(pc[1], 3.91e6, 0.0, "1-pentanol Pc changed");
        near(omega[1], 0.5784, 0.0, "1-pentanol omega changed");
        near(tc[2], 676.0, 0.0, "n-tridecane Tc changed");
        near(pc[2], 1.72e6, 0.0, "n-tridecane Pc changed");
        near(omega[2], 0.6203, 0.0, "n-tridecane omega changed");
        near(parameters.kij(0U, 1U), 0.17, 0.0,
             "Mushrif Figure 3.4 CO2/pentanol BIP changed");
        near(parameters.kij(0U, 2U), 0.15, 0.0,
             "Mushrif Figure 3.4 CO2/tridecane BIP changed");
        near(parameters.kij(1U, 2U), 0.10, 0.0,
             "Mushrif Figure 3.4 pentanol/tridecane BIP changed");

        for (const auto& component : parameters.components().items()) {
            check_source(component.definition);
        }
        for (const auto& record : parameters.pure_records()) {
            check_source(record.critical_temperature->source);
            check_source(record.critical_pressure->source);
            check_source(record.acentric_factor->source);
        }
        for (const auto& record : parameters.binary_records()) {
            check_source(record.kij->source);
        }
        check_source(parameters.applicability().declaration);
        return 0;
    } catch (const std::exception&) {
        return 1;
    }
}
