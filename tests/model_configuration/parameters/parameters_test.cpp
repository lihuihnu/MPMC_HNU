#include <mpmc/model_configuration/pr76_parameters.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

bool public_header_compiles();
bool adapter_header_compiles();

namespace {
namespace mc = mpmc::model_configuration;
namespace th = mpmc::thermodynamics;
using Code = mc::ModelConfigurationErrorCode;

void require(bool condition, std::string_view reason,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " + std::string(reason));
    }
}
template <class F>
void expect_error(Code code, F&& action, std::string_view field = {}) {
    try { action(); }
    catch (const mc::ModelConfigurationError& error) {
        require(error.code() == code, error.what());
        require(!error.field().empty(), "missing field diagnostic");
        if (!field.empty()) { require(error.field() == field, error.field()); }
        return;
    }
    throw std::runtime_error("invalid configuration accepted");
}

// Artificial software-contract values ONLY. No real fluids, experimental
// validation, literature fit, EOS outputs or independent flash reference here.
mc::ModelProvenance source(std::string locator = "definition") {
    return {mc::SourceKind::synthetic_test, "test://model-parameters", "fixture-v1",
            std::move(locator), "Artificial contract fixture; not physical validation",
            "Constructed in this test", "No third-party data copied"};
}
mc::ModelScalar scalar(double value, std::string unit) {
    return {value, source("scalar"), std::move(unit), "identity"};
}
mc::ThermodynamicModelDefinition fixture(std::size_t n = 3) {
    mc::ThermodynamicModelDefinition d;
    d.version = mc::model_parameter_definition_v1;
    d.family = mc::ThermodynamicModelFamily::peng_robinson_1976;
    d.display_name = "Synthetic model";
    d.dataset_id = "synthetic-contract";
    d.revision = "r1";
    d.provenance = source();
    d.applicability.provenance = source("unknown bounds");
    mc::Pr76ParameterDefinition pr;
    for (std::size_t i = 0; i < n; ++i) {
        const auto id = "test:" + std::to_string(i);
        const auto x = static_cast<double>(i);
        d.components.push_back({id, "Synthetic component", mc::ComponentKind::pure,
                                source(id), scalar(0.02 + 0.01 * x, "kg/mol")});
        pr.pure.push_back({id, scalar(200 + 100 * x, "K"),
                           scalar(1e6 + 1e6 * x, "Pa"), scalar(-0.1 + 0.1 * x, "1")});
        for (std::size_t j = 0; j < i; ++j) {
            pr.binary.push_back({"test:" + std::to_string(j), id,
                                 scalar(0.01 * static_cast<double>(i + j), "1")});
        }
    }
    d.parameters = std::move(pr);
    return d;
}
mc::Pr76ParameterDefinition& pr(mc::ThermodynamicModelDefinition& d) {
    return std::get<mc::Pr76ParameterDefinition>(d.parameters);
}
mc::Pr76ModelParameters prepare(const mc::ThermodynamicModelDefinition& d,
                               mc::ModelConfigurationLimits limits = {}) {
    return mc::Pr76ModelParameters::create(d, limits, mc::ModelDataPolicy::allow_synthetic_tests);
}

void complete_models() {
    for (const std::size_t n : {1U, 2U, 3U}) {
        auto d = fixture(n);
        d.components.back().kind = mc::ComponentKind::pseudo;
        const auto prepared = prepare(d);
        const auto& p = prepared.parameters();
        require(p.components().size() == n && p.pure_records().size() == n, "component count");
        require(p.binary_records().size() == n * (n - 1) / 2, "pair count");
        require(p.components().at(n - 1).kind == th::ComponentKind::pseudo, "pseudo lost");
        for (std::size_t i = 0; i < n; ++i) {
            const double x = static_cast<double>(i);
            require(p.critical_temperatures_k()[i] == 200 + 100 * x, "Tc mapping");
            require(p.critical_pressures_pa()[i] == 1e6 + 1e6 * x, "Pc mapping");
            require(p.acentric_factors()[i] == -0.1 + 0.1 * x, "omega mapping");
            require(p.kij(i, i) == 0, "structural diagonal");
            for (std::size_t j = 0; j < i; ++j) {
                require(p.kij(i, j) == 0.01 * static_cast<double>(i + j) &&
                        p.kij(i, j) == p.kij(j, i), "symmetric pair mapping");
            }
        }
        require(p.dataset_id() == d.dataset_id && p.revision() == d.revision, "dataset metadata");
    }
}
void order_mutations() {
    auto d = fixture();
    const auto original = prepare(d);
    std::array<std::size_t, 3> order{0, 1, 2};
    do {
        auto reordered = d;
        for (std::size_t i = 0; i < 3; ++i) { reordered.components[i] = d.components[order[i]]; }
        // Parameter storage order is independent of selected component order.
        std::reverse(pr(reordered).pure.begin(), pr(reordered).pure.end());
        const auto result = prepare(reordered);
        for (std::size_t i = 0; i < 3; ++i) {
            require(result.parameters().critical_temperatures_k()[i] ==
                    original.parameters().critical_temperatures_k()[order[i]], "reordered Tc");
            for (std::size_t j = 0; j < 3; ++j) {
                require(result.parameters().kij(i, j) == original.parameters().kij(order[i], order[j]),
                        "reordered matrix");
            }
        }
    } while (std::next_permutation(order.begin(), order.end()));
    // Remove, then add a complete component using the same runtime interface.
    d.components.pop_back();
    pr(d).pure.pop_back();
    pr(d).binary.resize(1);
    const auto removed = prepare(d);
    require(removed.parameters().components().size() == 2, "removal");
    d = fixture(4);
    const auto added = prepare(d);
    require(added.parameters().components().size() == 4, "addition");
    d.components[0].component_id = "replacement";
    expect_error(Code::unknown_component, [&] { (void)prepare(d); });
    pr(d).pure[0].component_id = "replacement";
    for (auto& pair : pr(d).binary) {
        if (pair.first_component_id == "test:0") { pair.first_component_id = "replacement"; }
    }
    const auto replaced = prepare(d);
    require(replaced.parameters().components().at(0).id == "replacement", "replacement");
}
void missing_parameters() {
    for (int field = 0; field < 4; ++field) {
        auto d = fixture(2);
        if (field == 0) { pr(d).pure[0].critical_temperature_k.reset(); }
        if (field == 1) { pr(d).pure[0].critical_pressure_pa.reset(); }
        if (field == 2) { pr(d).pure[0].acentric_factor.reset(); }
        if (field == 3) { pr(d).binary[0].kij.reset(); }
        expect_error(Code::missing_parameter, [&] { (void)prepare(d); });
    }
    auto d = fixture(2);
    pr(d).binary.clear();
    expect_error(Code::missing_parameter, [&] { (void)prepare(d); }, "parameters.binary");
    d = fixture(2);
    pr(d).pure.pop_back();
    expect_error(Code::missing_parameter, [&] { (void)prepare(d); });
    d = fixture(2);
    pr(d).binary[0].kij->value = 0;
    const auto zero = prepare(d);
    require(zero.parameters().kij(0, 1) == 0, "explicit zero must remain legal");
    pr(d).pure[0].critical_temperature_k = mc::ModelScalar{};
    expect_error(Code::invalid_value, [&] { (void)prepare(d); });
}
void identity_pairs() {
    auto d = fixture();
    d.components[1].component_id = d.components[0].component_id;
    expect_error(Code::duplicate_identifier, [&] { (void)prepare(d); });
    d = fixture();
    d.components[0].component_id = "";
    expect_error(Code::invalid_identifier, [&] { (void)prepare(d); });
    d = fixture();
    pr(d).pure.push_back(pr(d).pure[0]);
    expect_error(Code::duplicate_parameter, [&] { (void)prepare(d); });
    for (bool reversed : {false, true}) {
        d = fixture();
        auto pair = pr(d).binary[0];
        if (reversed) { std::swap(pair.first_component_id, pair.second_component_id); }
        pr(d).binary.push_back(pair);
        expect_error(Code::duplicate_parameter, [&] { (void)prepare(d); });
    }
    d = fixture();
    pr(d).binary[0].second_component_id = pr(d).binary[0].first_component_id;
    expect_error(Code::invalid_pair, [&] { (void)prepare(d); });
    d = fixture();
    pr(d).binary[0].second_component_id = "absent";
    expect_error(Code::unknown_component, [&] { (void)prepare(d); });
    d = fixture();
    d.components[0].kind = static_cast<mc::ComponentKind>(99);
    expect_error(Code::invalid_value, [&] { (void)prepare(d); });
}
void numeric_domains() {
    for (double bad : {std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::infinity(),
                       -std::numeric_limits<double>::infinity()}) {
        for (int field = 0; field < 5; ++field) {
            auto d = fixture();
            if (field == 0) { pr(d).pure[0].critical_temperature_k->value = bad; }
            if (field == 1) { pr(d).pure[0].critical_pressure_pa->value = bad; }
            if (field == 2) { pr(d).pure[0].acentric_factor->value = bad; }
            if (field == 3) { pr(d).binary[0].kij->value = bad; }
            if (field == 4) { d.components[0].molar_mass_kg_per_mol->value = bad; }
            expect_error(Code::invalid_value, [&] { (void)prepare(d); });
        }
    }
    for (double bad : {0.0, -1.0}) {
        auto d = fixture();
        pr(d).pure[0].critical_temperature_k->value = bad;
        expect_error(Code::invalid_value, [&] { (void)prepare(d); },
                     "parameters.pure[test:0].critical_temperature_k");
        d = fixture();
        pr(d).pure[0].critical_pressure_pa->value = bad;
        expect_error(Code::invalid_value, [&] { (void)prepare(d); });
    }
    auto d = fixture();
    pr(d).pure[0].acentric_factor->value = -0.5;
    pr(d).binary[0].kij->value = -1.2;
    const auto unchanged = prepare(d);
    require(unchanged.parameters().kij(0, 1) == -1.2, "no invented kij range or clipping");
}
void provenance() {
    auto d = fixture();
    expect_error(Code::invalid_source, [&] { (void)mc::Pr76ModelParameters::create(d); });
    auto& pc = *pr(d).pure[0].critical_pressure_pa;
    pc.original_unit = "MPa";
    pc.conversion = "MPa * 1e6 -> Pa";
    pc.provenance.locator = "canonical pressure already converted";
    const auto p = prepare(d);
    const auto& mapped = *p.parameters().pure_records()[0].critical_pressure;
    require(mapped.value == pc.value && mapped.unit == th::Unit::pascal, "no duplicate conversion");
    require(mapped.original_unit == pc.original_unit && mapped.conversion == pc.conversion &&
            mapped.source.locator == pc.provenance.locator, "source metadata lost");
    require(p.definition().provenance.note == d.provenance.note, "model source lost");
    require(p.parameters().components().at(0).molar_mass->value == 0.02, "molar mass lost");
    for (int field = 0; field < 5; ++field) {
        auto invalid = d;
        if (field == 0) { invalid.provenance.reference.clear(); }
        if (field == 1) { invalid.components[0].provenance.revision.clear(); }
        if (field == 2) { pr(invalid).pure[0].critical_pressure_pa->provenance.locator.clear(); }
        if (field == 3) { pr(invalid).binary[0].kij->conversion.clear(); }
        if (field == 4) { invalid.applicability.provenance.acquisition.clear(); }
        expect_error(Code::missing_field, [&] { (void)prepare(invalid); });
    }
    d.provenance.kind = static_cast<mc::SourceKind>(99);
    expect_error(Code::invalid_source, [&] { (void)prepare(d); });
    d = fixture();
    // Test source-kind transport only; these values remain a synthetic fixture.
    d.provenance.kind = mc::SourceKind::user_supplied;
    const auto user_source = prepare(d);
    require(user_source.definition().provenance.kind == mc::SourceKind::user_supplied,
            "user source must not be upgraded to literature");
}
void applicability() {
    auto d = fixture();
    const auto unknown = prepare(d);
    require(unknown.parameters().applicability().assess(300, 1e6) == th::RangeAssessment::unknown,
            "absent bounds are unknown");
    d.applicability.temperature_lower_k = 200;
    expect_error(Code::invalid_range, [&] { (void)prepare(d); });
    d.applicability.temperature_upper_k = 400;
    const auto temperature_only = prepare(d);
    require(temperature_only.parameters().applicability().assess(300, 1e6) == th::RangeAssessment::unknown,
            "missing pressure is unknown");
    require(temperature_only.parameters().applicability().assess(401, 1e6) ==
            th::RangeAssessment::outside_declared_bounds, "declared bound must still reject");
    d.applicability.pressure_lower_pa = 1e5;
    d.applicability.pressure_upper_pa = 2e6;
    const auto bounded = prepare(d);
    require(bounded.parameters().applicability().assess(200, 2e6) ==
            th::RangeAssessment::inside_declared_bounds, "closed endpoints");
    for (double bad : {0.0, 500.0, std::numeric_limits<double>::infinity()}) {
        d.applicability.temperature_lower_k = bad;
        expect_error(Code::invalid_range, [&] { (void)prepare(d); });
    }
}
void immutability() {
    static_assert(!std::is_copy_assignable_v<mc::Pr76ModelParameters>);
    static_assert(!std::is_move_assignable_v<mc::Pr76ModelParameters>);
    auto draft = fixture();
    const auto a = prepare(draft);
    draft.revision = "r2";
    draft.display_name = "Changed";
    pr(draft).binary[0].kij->value = 0.25;
    const auto b = prepare(draft);
    draft = {};
    require(a.parameters().kij(0, 1) == 0.01 && b.parameters().kij(0, 1) == 0.25,
            "parameter snapshots alias");
    require(a.definition().revision == "r1" && b.definition().revision == "r2", "revision alias");
    require(a.definition().display_name == "Synthetic model", "label alias");
    const auto copy = a;
    require(copy.parameters().components().at(0).id == "test:0", "owning snapshot copy");
}
void quotas() {
    const auto d = fixture();
    for (int field = 0; field < 7; ++field) {
        mc::ModelConfigurationLimits limits;
        if (field == 0) { limits.max_components = 2; }
        if (field == 1) { limits.max_pair_records = 2; }
        if (field == 2) { limits.max_matrix_entries = 8; }
        if (field == 3) { limits.max_identifier_bytes = 2; }
        if (field == 4) { limits.max_display_name_bytes = 2; }
        if (field == 5) { limits.max_provenance_field_bytes = 2; }
        if (field == 6) { limits.max_total_text_bytes = 2; }
        expect_error(Code::resource_limit, [&] { (void)prepare(d, limits); });
    }
    mc::ModelConfigurationLimits exact;
    exact.max_components = 3;
    exact.max_pair_records = 3;
    exact.max_matrix_entries = 9;
    const auto within = prepare(d, exact);
    require(within.parameters().components().size() == 3, "exact shape quota");
    exact.max_total_text_bytes = 2048;
    expect_error(Code::resource_limit, [&] { (void)prepare(d, exact); });
    auto oversized = d;
    pr(oversized).pure.resize(257);
    expect_error(Code::resource_limit, [&] { (void)prepare(oversized); });
    oversized = d;
    pr(oversized).binary[0].second_component_id.assign(129, 'x');
    expect_error(Code::resource_limit, [&] { (void)prepare(oversized); },
                 "parameters.binary[0].second_component_id");
    oversized = d;
    pr(oversized).binary[0].kij->provenance.note.assign(4097, 'x');
    expect_error(Code::resource_limit, [&] { (void)prepare(oversized); },
                 "parameters.binary[0].kij.provenance.note");
}
void versions() {
    for (auto family : {mc::ThermodynamicModelFamily::unspecified,
                        mc::ThermodynamicModelFamily::soreide_whitson_1992,
                        mc::ThermodynamicModelFamily::cubic_plus_association,
                        static_cast<mc::ThermodynamicModelFamily>(99)}) {
        auto d = fixture();
        d.family = family;
        expect_error(Code::unsupported_family, [&] { (void)prepare(d); }, "family");
    }
    auto d = fixture();
    d.version.clear();
    expect_error(Code::unsupported_version, [&] { (void)prepare(d); });
    d = fixture();
    d.parameters = std::monostate{};
    expect_error(Code::missing_field, [&] { (void)prepare(d); });
    d = fixture();
    expect_error(Code::invalid_source, [&] {
        (void)mc::Pr76ModelParameters::create(d, {}, static_cast<mc::ModelDataPolicy>(99));
    });
}
} // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "one test case required");
        const std::string_view name(argv[1]);
        if (name == "complete_models") { complete_models(); }
        else if (name == "order_mutations") { order_mutations(); }
        else if (name == "missing_parameters") { missing_parameters(); }
        else if (name == "identity_pairs") { identity_pairs(); }
        else if (name == "numeric_domains") { numeric_domains(); }
        else if (name == "provenance") { provenance(); }
        else if (name == "applicability") { applicability(); }
        else if (name == "immutability") { immutability(); }
        else if (name == "quotas") { quotas(); }
        else if (name == "versions") { versions(); }
        else if (name == "headers") { require(public_header_compiles() && adapter_header_compiles(), "headers"); }
        else { throw std::runtime_error("unknown test case"); }
        std::cout << "PASS " << name << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
