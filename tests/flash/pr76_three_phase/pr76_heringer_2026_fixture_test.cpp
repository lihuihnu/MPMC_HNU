#include "heringer_2026_sour_gas_fixture.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace hg = pr76_heringer_2026_sour_gas_test;

namespace {

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

void require_close(double actual, double expected, double tolerance, const char* message) {
    if (!(std::abs(actual - expected) <= tolerance)) {
        throw std::runtime_error(message);
    }
}

} // namespace

int main() {
    try {
        const auto model = hg::model();
        const auto& parameters = model.parameters();
        require(parameters.dataset_id() == "Heringer-et-al-2026-sour-gas-PR76",
                "Heringer fixture lost dataset identity");
        require(parameters.revision() == "appendix-b-table-b1-v1",
                "Heringer fixture lost dataset revision");
        require(parameters.critical_temperatures_k().size() == hg::canonical_ids.size(),
                "Heringer fixture changed component count");
        require(parameters.pure_records().size() == hg::canonical_ids.size(),
                "Heringer fixture lost pure records");
        require(parameters.binary_records().size() == 15U,
                "Heringer fixture did not publish a complete six-component pair set");

        constexpr std::string_view doi =
            "https://doi.org/10.1016/j.fluid.2025.114653";
        for (std::size_t i = 0; i < hg::canonical_ids.size(); ++i) {
            require(parameters.critical_temperatures_k()[i] == hg::critical_temperatures_k[i],
                    "Heringer fixture changed a Table B1 critical temperature");
            require(parameters.critical_pressures_pa()[i] == hg::critical_pressures_pa[i],
                    "Heringer fixture changed a Table B1 critical pressure");
            require(parameters.acentric_factors()[i] == hg::acentric_factors[i],
                    "Heringer fixture changed a Table B1 acentric factor");
            const auto& record = parameters.pure_records()[i];
            require(record.critical_temperature && record.critical_pressure &&
                        record.acentric_factor,
                    "Heringer fixture lost a sourced pure-component field");
            require(record.critical_temperature->source.reference == doi &&
                        record.critical_pressure->source.reference == doi &&
                        record.acentric_factor->source.reference == doi,
                    "Heringer fixture lost Table B1 DOI provenance");
        }

        for (std::size_t i = 0; i < hg::canonical_ids.size(); ++i) {
            require(parameters.kij(i, i) == 0.0,
                    "Heringer fixture changed structural diagonal kij");
            for (std::size_t j = i + 1U; j < hg::canonical_ids.size(); ++j) {
                require(parameters.kij(i, j) == hg::kij(i, j) &&
                            parameters.kij(j, i) == hg::kij(i, j),
                        "Heringer fixture changed a Table B1 symmetric kij");
            }
        }

        double raw_total = 0.0;
        for (const double value : hg::table_b1_mole_percent) { raw_total += value; }
        require_close(raw_total, 99.98, 1.0e-12,
                      "Heringer Table B1 rounded Mol% total changed");
        const auto feed = hg::table_b1_feed();
        double normalized_total = 0.0;
        for (const double value : feed) { normalized_total += value; }
        require_close(normalized_total, 1.0, 2.0e-15,
                      "Heringer Table B1 feed normalization changed");
        require_close(feed[0], 70.59 / 99.98, 2.0e-15,
                      "Heringer Table B1 CO2 normalization changed");

        const auto reverse_order = hg::reversed_order();
        const auto reverse_model = hg::model(reverse_order);
        const auto& reverse_parameters = reverse_model.parameters();
        const auto reverse_feed = hg::table_b1_feed(reverse_order);
        for (std::size_t i = 0; i < reverse_order.size(); ++i) {
            const std::size_t canonical = hg::canonical_index(reverse_order[i]);
            require(reverse_parameters.critical_temperatures_k()[i] ==
                        hg::critical_temperatures_k[canonical],
                    "Heringer fixture lost pure-parameter alignment after reordering");
            require_close(reverse_feed[i], feed[canonical], 0.0,
                          "Heringer fixture lost feed alignment after reordering");
        }

        std::cout << "PASS heringer_2026_sour_gas_fixture\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
