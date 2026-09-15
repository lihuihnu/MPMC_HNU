#include <mpmc/model_configuration/pr76_parameters.hpp>

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {
namespace mc = mpmc::model_configuration;
namespace th = mpmc::thermodynamics;
using Assessment = mc::ModelApplicabilityAssessment;
using Code = mc::ModelConfigurationErrorCode;

void require(bool value, std::string_view message) {
    if (!value) { throw std::runtime_error(std::string(message)); }
}

mc::ModelProvenance source() {
    return {mc::SourceKind::synthetic_test, "test://applicability", "v1", "fixture",
            "software contract only", "constructed in test", "no third-party data"};
}

mc::ModelScalar scalar(double value, std::string unit) {
    return {value, source(), std::move(unit), "identity"};
}

mc::ThermodynamicModelDefinition fixture() {
    mc::ThermodynamicModelDefinition d;
    d.version = mc::model_parameter_definition_v1;
    d.family = mc::ThermodynamicModelFamily::peng_robinson_1976;
    d.display_name = "Synthetic applicability fixture";
    d.dataset_id = "applicability-fixture";
    d.revision = "v1";
    d.provenance = source();
    d.applicability.provenance = source();
    d.components.push_back({"test:0", "Synthetic component", mc::ComponentKind::pure,
                            source(), scalar(0.02, "kg/mol")});
    mc::Pr76ParameterDefinition pr;
    pr.pure.push_back({"test:0", scalar(300.0, "K"), scalar(5.0e6, "Pa"),
                       scalar(0.1, "1")});
    d.parameters = std::move(pr);
    return d;
}

mc::Pr76ModelParameters prepare(const mc::ThermodynamicModelDefinition& d) {
    return mc::Pr76ModelParameters::create(d, {}, mc::ModelDataPolicy::allow_synthetic_tests);
}

template <class F>
void expect_error(Code code, std::string_view field, F&& action) {
    try { action(); }
    catch (const mc::ModelConfigurationError& error) {
        require(error.code() == code, "wrong error code");
        require(error.field() == field, error.field());
        return;
    }
    throw std::runtime_error("invalid applicability accepted");
}

void run() {
    mc::ModelApplicability a;
    require(a.assess(300.0, 1.0e6) == Assessment::unknown, "absent bounds must be unknown");

    a.temperature_lower_k = 250.0;
    require(a.assess_temperature(249.0) == Assessment::outside_declared_bounds,
            "known one-sided lower violation not outside");
    require(a.assess_temperature(250.0) == Assessment::unknown,
            "closed one-sided endpoint should remain unknown, not inside");
    require(a.assess_temperature(400.0) == Assessment::unknown,
            "missing upper endpoint was treated as infinity");
    a.temperature_lower_exclusive = true;
    require(a.assess_temperature(250.0) == Assessment::outside_declared_bounds,
            "open lower endpoint equality must be outside");

    a.temperature_upper_k = 400.0;
    require(a.assess_temperature(251.0) == Assessment::inside_declared_bounds,
            "complete open/closed interval interior not inside");
    require(a.assess_temperature(400.0) == Assessment::inside_declared_bounds,
            "closed upper endpoint equality not inside");
    a.temperature_upper_exclusive = true;
    require(a.assess_temperature(400.0) == Assessment::outside_declared_bounds,
            "open upper endpoint equality must be outside");

    a.pressure_upper_pa = 2.0e6;
    require(a.assess(300.0, 2.1e6) == Assessment::outside_declared_bounds,
            "one definite axis violation must dominate unknown axis");
    require(a.assess(300.0, 1.0e6) == Assessment::inside_declared_bounds,
            "complete temperature and pressure bounds should be inside");
    require(a.assess(std::numeric_limits<double>::quiet_NaN(), 1.0e6) == Assessment::unknown,
            "non-finite query assessment must not claim validity");

    auto d = fixture();
    d.applicability.temperature_lower_k = 250.0;
    const auto one_sided = prepare(d);
    require(!one_sided.parameters().applicability().temperature_k,
            "one-sided public bound fabricated a native opposite endpoint");
    require(one_sided.definition().applicability.temperature_lower_k == 250.0,
            "public one-sided endpoint lost from immutable snapshot");

    d.applicability.temperature_upper_k = 400.0;
    d.applicability.temperature_lower_exclusive = true;
    const auto complete = prepare(d);
    require(complete.parameters().applicability().temperature_k.has_value(),
            "complete interval did not map to native envelope");
    require(complete.parameters().applicability().temperature_k->lower == 250.0 &&
            complete.parameters().applicability().temperature_k->upper == 400.0,
            "complete interval numeric envelope changed");
    require(complete.definition().applicability.temperature_lower_exclusive,
            "exclusive marker lost from public snapshot");

    d = fixture();
    d.applicability.temperature_lower_exclusive = true;
    expect_error(Code::invalid_range, "applicability.temperature_lower_exclusive",
                 [&] { (void)prepare(d); });

    d = fixture();
    d.applicability.temperature_lower_k = 400.0;
    d.applicability.temperature_upper_k = 300.0;
    expect_error(Code::invalid_range, "applicability.temperature_k",
                 [&] { (void)prepare(d); });

    d = fixture();
    d.applicability.temperature_lower_k = 300.0;
    d.applicability.temperature_upper_k = 300.0;
    (void)prepare(d); // Closed singleton is a valid explicit declaration.
    d.applicability.temperature_upper_exclusive = true;
    expect_error(Code::invalid_range, "applicability.temperature_k",
                 [&] { (void)prepare(d); });

    for (double bad : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::quiet_NaN()}) {
        d = fixture();
        d.applicability.pressure_upper_pa = bad;
        expect_error(Code::invalid_range, "applicability.pressure_upper_pa",
                     [&] { (void)prepare(d); });
    }
}
} // namespace

int main() {
    try {
        run();
        std::cout << "PASS applicability_endpoint_contract\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
