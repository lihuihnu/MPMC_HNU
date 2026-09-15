#include "test_support.hpp"

#include <iostream>
#include <stdexcept>
#include <string_view>

namespace {
using namespace model_test;

void expect_rejected(mc::Pr76ExecutableModel& model, const fl::PtFlashRequest& request,
                     std::string_view field) {
    try { (void)model.solve(request); }
    catch (const mc::Pr76SolveRequestError& error) {
        require(error.field() == field, "wrong applicability rejection field");
        const auto* domain = dynamic_cast<const std::domain_error*>(&error);
        require(domain != nullptr, "applicability rejection changed standard category");
        return;
    }
    throw std::runtime_error("out-of-applicability request was accepted");
}

void run() {
    const auto native = pr76_max3_test::model();

    auto lower_only = definition(native);
    lower_only.applicability.temperature_lower_k = 250.0;
    auto lower_model = create(lower_only, preset());
    require(lower_model->parameter_snapshot().definition().applicability.assess_temperature(249.0) ==
                mc::ModelApplicabilityAssessment::outside_declared_bounds,
            "public one-sided lower assessment changed");
    require(!lower_model->parameter_snapshot().parameters().applicability().temperature_k,
            "one-sided lower fabricated native interval");
    expect_rejected(*lower_model, {1.0e6, 249.0, single.feed}, "temperature_k");
    compare(lower_model->solve(single), direct(native, single));

    auto upper_open = definition(native);
    upper_open.applicability.pressure_upper_pa = 1.0e6;
    upper_open.applicability.pressure_upper_exclusive = true;
    auto upper_model = create(upper_open, preset());
    require(upper_model->parameter_snapshot().definition().applicability.assess_pressure(1.0e6) ==
                mc::ModelApplicabilityAssessment::outside_declared_bounds,
            "public open upper endpoint equality changed");
    expect_rejected(*upper_model, single, "pressure_pa");
    const fl::PtFlashRequest below_open{0.9e6, 250.0, single.feed};
    compare(upper_model->solve(below_open), direct(native, below_open));

    auto mixed = definition(native);
    mixed.applicability.temperature_lower_k = 240.0;
    mixed.applicability.temperature_upper_k = 260.0;
    mixed.applicability.temperature_upper_exclusive = true;
    mixed.applicability.pressure_lower_pa = 0.9e6;
    mixed.applicability.pressure_lower_exclusive = true;
    mixed.applicability.pressure_upper_pa = 1.1e6;
    auto mixed_model = create(mixed, preset());
    require(mixed_model->parameter_snapshot().parameters().applicability().temperature_k.has_value() &&
            mixed_model->parameter_snapshot().parameters().applicability().pressure_pa.has_value(),
            "complete interval did not retain native closed envelope");
    expect_rejected(*mixed_model, {1.0e6, 260.0, single.feed}, "temperature_k");
    expect_rejected(*mixed_model, {0.9e6, 250.0, single.feed}, "pressure_pa");
    const fl::PtFlashRequest closed_edges{1.1e6, 240.0, single.feed};
    compare(mixed_model->solve(closed_edges), direct(native, closed_edges));

    // Public hints are orthogonal to applicability. The same public state gate
    // applies before a hinted native solve and does not persist or mutate hints.
    auto hints = mc::make_pt_solve_hints_v1();
    hints.initial_stability_starts = {single.feed};
    try { (void)upper_model->solve(single, hints); }
    catch (const mc::Pr76SolveRequestError& error) {
        require(error.field() == "pressure_pa", "hinted applicability field changed");
        compare(upper_model->solve(below_open), direct(native, below_open));
        return;
    }
    throw std::runtime_error("hinted open-endpoint request was accepted");
}
} // namespace

int main() {
    try {
        run();
        std::cout << "PASS executable_applicability_endpoints\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
