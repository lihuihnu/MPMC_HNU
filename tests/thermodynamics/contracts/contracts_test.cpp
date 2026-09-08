#include <mpmc/thermodynamics/components.hpp>
#include <mpmc/thermodynamics/pr_parameters.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <exception>
#include <iostream>
#include <limits>
#include <optional>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

bool components_header_compiles();
std::string_view pr_profile_from_separate_translation_unit();

namespace {
namespace th = mpmc::thermodynamics;
constexpr auto testing = th::DataPolicy::allow_synthetic_tests;
using Code = th::ContractErrorCode;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " + std::string(message));
    }
}

template <typename Function>
void expect_error(Code expected, Function&& function) {
    bool caught = false;
    try {
        function();
    } catch (const th::ContractError& error) {
        require(error.code() == expected, "wrong structured error code");
        require(!error.field().empty(), "missing diagnostic field context");
        caught = true;
    }
    require(caught, "invalid contract was accepted");
}

// Every numerical value below is an artificial software fixture. None represents
// measured/estimated data for a real chemical, and none validates an EOS model.
th::Provenance source(std::string locator) {
    return {th::SourceKind::synthetic_test, "test://thermodynamics/contracts", "fixture-v1",
            std::move(locator), "Artificial test only; not experimental or fitted physical data",
            "Constructed in this test source", "Synthetic fixture; no third-party data copied"};
}
th::SourcedScalar datum(double value, th::Unit unit, std::string locator) {
    std::string unit_name;
    switch (unit) {
    case th::Unit::kelvin: unit_name = "K"; break;
    case th::Unit::pascal: unit_name = "Pa"; break;
    case th::Unit::kilogram_per_mole: unit_name = "kg/mol"; break;
    default: unit_name = "1"; break;
    }
    return {value, unit, source(std::move(locator)), std::move(unit_name), "identity"};
}
struct Fixture {
    std::vector<th::Component> catalog;
    std::vector<std::string> order;
    th::PrParameterInput data;

    [[nodiscard]] th::PrParameterSet build(th::ContractLimits limits = {}) const {
        return th::PrParameterSet::create(catalog, order, data, testing, limits);
    }
};
Fixture fixture() {
    Fixture f;
    f.order = {"test:a", "test:b", "test:c", "test:d"};
    const std::array<double, 4> tc{200, 300, 400, 500};
    const std::array<double, 4> pc{1e6, 2e6, 3e6, 4e6};
    const std::array<double, 4> omega{-0.1, 0.0, 0.2, 0.6};
    for (std::size_t i = 0; i < f.order.size(); ++i) {
        const auto& id = f.order[i];
        f.catalog.push_back({id, "Synthetic label", i == 3 ? th::ComponentKind::pseudo
                                                          : th::ComponentKind::pure,
                             source(id + ".definition"), std::nullopt});
        f.data.pure.push_back({id, datum(tc[i], th::Unit::kelvin, id + ".Tc"),
                              datum(pc[i], th::Unit::pascal, id + ".Pc"),
                              datum(omega[i], th::Unit::dimensionless, id + ".omega")});
    }
    f.catalog[0].molar_mass = datum(0.01, th::Unit::kilogram_per_mole, "test:a.mass");
    f.data.model_id = th::pr76_profile;
    f.data.dataset_id = "synthetic-four-component-set";
    f.data.revision = "test-v1";
    f.data.applicability = {th::ClosedInterval{250, 450}, th::ClosedInterval{1e5, 5e6},
                            source("declared-synthetic-domain")};
    f.data.binary = {
        {"test:b", "test:a", datum(0.125, th::Unit::dimensionless, "AB")},
        {"test:a", "test:c", datum(-0.25, th::Unit::dimensionless, "AC")},
        {"test:d", "test:a", datum(1.25, th::Unit::dimensionless, "AD")},
        {"test:c", "test:b", datum(0.5, th::Unit::dimensionless, "BC")},
        {"test:b", "test:d", datum(-0.125, th::Unit::dimensionless, "BD")},
        {"test:d", "test:c", datum(0.75, th::Unit::dimensionless, "CD")}};
    std::rotate(f.data.pure.begin(), f.data.pure.begin() + 2, f.data.pure.end());
    return f;
}

// Independent identity-indexed expectation, not extracted from the built snapshot.
void check_binding(const th::PrParameterSet& result, const std::vector<std::string>& order) {
    const std::array<std::string, 4> ids{"test:a", "test:b", "test:c", "test:d"};
    const std::array<double, 4> tc{200, 300, 400, 500};
    const std::array<double, 4> pc{1e6, 2e6, 3e6, 4e6};
    const std::array<double, 4> omega{-0.1, 0.0, 0.2, 0.6};
    const std::array<std::array<double, 4>, 4> kij{{
        {0, 0.125, -0.25, 1.25}, {0.125, 0, 0.5, -0.125},
        {-0.25, 0.5, 0, 0.75}, {1.25, -0.125, 0.75, 0}}};
    const auto original_index = [&](const std::string& id) {
        const auto it = std::find(ids.begin(), ids.end(), id);
        require(it != ids.end(), "unexpected fixture ID");
        return static_cast<std::size_t>(it - ids.begin());
    };
    require(result.components().size() == order.size(), "incorrect active count");
    require(result.kij_matrix().size() == order.size() * order.size(), "incorrect matrix shape");
    for (std::size_t i = 0; i < order.size(); ++i) {
        const auto k = original_index(order[i]);
        require(result.components().at(i).id == order[i], "component order mismatch");
        require(result.components().index_of(order[i]) == i, "index map mismatch");
        require(result.pure_records()[i].component_id == order[i], "pure record order mismatch");
        require(result.critical_temperatures_k()[i] == tc[k], "Tc bound by position, not ID");
        require(result.critical_pressures_pa()[i] == pc[k], "Pc bound by position, not ID");
        require(result.acentric_factors()[i] == omega[k], "omega bound by position, not ID");
        for (std::size_t j = 0; j < order.size(); ++j) {
            require(result.kij(i, j) == kij[k][original_index(order[j])],
                    "pair permutation mismatch");
        }
    }
    std::size_t pair = 0;
    for (std::size_t i = 0; i < order.size(); ++i) {
        for (std::size_t j = i + 1; j < order.size(); ++j) {
            const auto& record = result.binary_records()[pair++];
            require(record.first_id == order[i] && record.second_id == order[j],
                    "pair record order");
            require(record.kij->value == result.kij(i, j), "pair provenance/numeric alignment");
            require(record.kij->source.kind == th::SourceKind::synthetic_test, "source label lost");
        }
    }
    require(result.binary_records().size() == pair, "unexpected diagonal/extra pair records");
}

void ordered_identity() {
    auto f = fixture();
    const auto ordered = th::OrderedComponents::select(f.catalog, f.order, testing);
    require(ordered.at(0).display_name == ordered.at(1).display_name, "labels need not be unique");
    require(ordered.at(3).kind == th::ComponentKind::pseudo, "pseudo identity lost");
    require(ordered.at(0).molar_mass->value == 0.01 && !ordered.at(1).molar_mass, "optional mass");
    f.order[0] = "test:A";
    expect_error(Code::unknown_component, [&] { (void)f.build(); });
}
void order_mutations() {
    auto f = fixture();
    const std::vector<std::vector<std::string>> orders{
        {"test:a", "test:b"}, {"test:c", "test:a", "test:b"}, {"test:b"},
        {"test:d", "test:b"}, {"test:b", "test:d"}, {"test:a", "test:b", "test:c", "test:d"}};
    for (const auto& order : orders) {
        f.order = order;
        const auto result = f.build();
        check_binding(result, order);
    }
}
void all_permutations() {
    auto f = fixture();
    std::size_t count = 0;
    do {
        const auto result = f.build();
        check_binding(result, f.order);
        ++count;
    } while (std::next_permutation(f.order.begin(), f.order.end()));
    require(count == 24, "all four-component permutations must be checked");
}
void identity_errors() {
    auto f = fixture();
    f.catalog.push_back(f.catalog.front());
    expect_error(Code::duplicate_identifier, [&] { (void)f.build(); });
    f = fixture(); f.order[1] = f.order[0];
    expect_error(Code::duplicate_identifier, [&] { (void)f.build(); });
    f = fixture(); f.order[1] = "absent";
    expect_error(Code::unknown_component, [&] { (void)f.build(); });
    f = fixture(); f.catalog[0].id = "bad id";
    expect_error(Code::invalid_identifier, [&] { (void)f.build(); });
    f = fixture(); f.catalog[0].display_name = " \t";
    expect_error(Code::missing_field, [&] { (void)f.build(); });
    f = fixture(); f.catalog[0].kind = th::ComponentKind::unspecified;
    expect_error(Code::invalid_value, [&] { (void)f.build(); });
    f = fixture(); f.catalog[3].definition.locator.clear();
    expect_error(Code::missing_field, [&] { (void)f.build(); });
}
void missing_parameters() {
    auto f = fixture(); f.data.pure.pop_back();
    expect_error(Code::missing_parameter, [&] { (void)f.build(); });
    for (int field = 0; field < 3; ++field) {
        f = fixture();
        auto& pure = f.data.pure[0];
        if (field == 0) { pure.critical_temperature.reset(); }
        if (field == 1) { pure.critical_pressure.reset(); }
        if (field == 2) { pure.acentric_factor.reset(); }
        expect_error(Code::missing_parameter, [&] { (void)f.build(); });
    }
    f = fixture(); f.data.binary.pop_back();
    expect_error(Code::missing_parameter, [&] { (void)f.build(); });
    f = fixture(); f.data.binary[0].kij.reset();
    expect_error(Code::missing_parameter, [&] { (void)f.build(); });
    // A database may lack data for unselected components, but present bad records
    // are not silently ignored. A one-component PR snapshot needs no fitted pair.
    f = fixture(); f.order = {"test:c"}; f.data.pure.resize(1); f.data.binary.clear();
    const auto single = f.build();
    require(single.kij(0, 0) == 0 && single.binary_records().empty(), "pure structural diagonal");
    f.order.push_back("test:b");
    expect_error(Code::missing_parameter, [&] { (void)f.build(); });
}
void duplicate_and_pair_errors() {
    auto f = fixture(); f.data.pure.push_back(f.data.pure[0]);
    expect_error(Code::duplicate_parameter, [&] { (void)f.build(); });
    f = fixture(); auto reversed = f.data.binary[0];
    std::swap(reversed.first_id, reversed.second_id); f.data.binary.push_back(reversed);
    expect_error(Code::duplicate_parameter, [&] { (void)f.build(); });
    f = fixture(); f.data.binary[0].second_id = f.data.binary[0].first_id;
    expect_error(Code::invalid_pair, [&] { (void)f.build(); });
    f = fixture(); f.data.binary[0].first_id = "unknown";
    expect_error(Code::unknown_component, [&] { (void)f.build(); });
    f = fixture(); f.data.pure[0].component_id = "unknown";
    expect_error(Code::unknown_component, [&] { (void)f.build(); });
}
void values_and_units() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (double invalid : {0.0, -1.0, nan, inf, -inf}) {
        auto f = fixture(); f.data.pure[0].critical_temperature->value = invalid;
        expect_error(Code::invalid_value, [&] { (void)f.build(); });
        f = fixture(); f.data.pure[0].critical_pressure->value = invalid;
        expect_error(Code::invalid_value, [&] { (void)f.build(); });
        f = fixture(); f.catalog[0].molar_mass->value = invalid;
        expect_error(Code::invalid_value, [&] { (void)f.build(); });
    }
    for (double invalid : {nan, inf, -inf}) {
        auto f = fixture(); f.data.pure[0].acentric_factor->value = invalid;
        expect_error(Code::invalid_value, [&] { (void)f.build(); });
        f = fixture(); f.data.binary[0].kij->value = invalid;
        expect_error(Code::invalid_value, [&] { (void)f.build(); });
    }
    auto f = fixture(); f.data.pure[0].critical_pressure->unit = th::Unit::kelvin;
    expect_error(Code::invalid_unit, [&] { (void)f.build(); });
    f = fixture(); f.data.pure[0].acentric_factor->unit = th::Unit::pascal;
    expect_error(Code::invalid_unit, [&] { (void)f.build(); });
    f = fixture(); f.data.binary[0].kij->unit = th::Unit::unspecified;
    expect_error(Code::invalid_unit, [&] { (void)f.build(); });
    f = fixture(); f.catalog[0].molar_mass->unit = th::Unit::dimensionless;
    expect_error(Code::invalid_unit, [&] { (void)f.build(); });
    f = fixture(); const auto valid = f.build();
    require(valid.acentric_factors()[0] == -0.1 && valid.kij(0, 3) == 1.25,
            "must not impose invented universal omega/kij bounds");
}
void provenance_policy() {
    auto f = fixture();
    expect_error(Code::invalid_source, [&] {
        (void)th::PrParameterSet::create(f.catalog, f.order, f.data); // No test-data opt-in.
    });
    f.data.pure[0].critical_temperature->source.reference.clear();
    expect_error(Code::missing_field, [&] { (void)f.build(); });
    f = fixture(); f.data.binary[0].kij->source.revision.clear();
    expect_error(Code::missing_field, [&] { (void)f.build(); });
    f = fixture(); f.data.pure[0].critical_pressure->original_unit.clear();
    expect_error(Code::missing_field, [&] { (void)f.build(); });
    f = fixture(); f.data.pure[0].critical_pressure->conversion.clear();
    expect_error(Code::missing_field, [&] { (void)f.build(); });
    f = fixture(); f.data.pure[0].critical_pressure->source.usage_terms.clear();
    expect_error(Code::missing_field, [&] { (void)f.build(); });
    f = fixture(); f.catalog[0].definition.acquisition.clear();
    expect_error(Code::missing_field, [&] { (void)f.build(); });
    f = fixture(); f.data.binary[0].kij->source.kind = th::SourceKind::unspecified;
    expect_error(Code::invalid_source, [&] { (void)f.build(); });
    f = fixture(); f.data.binary[0].kij->source.kind = th::SourceKind::assumption;
    f.data.binary[0].kij->source.note.clear();
    expect_error(Code::missing_field, [&] { (void)f.build(); });
    f.data.binary[0].kij->source.note = "Explicit zero assumption in this synthetic software test";
    f.data.binary[0].kij->value = 0;
    const auto explicit_zero = f.build();
    require(explicit_zero.kij(0, 1) == 0,
            "explicit zero assumption must not be treated as missing");
    require(explicit_zero.binary_records()[0].kij->source.kind == th::SourceKind::assumption,
            "assumption label discarded");
    // Validation of non-synthetic source kinds uses metadata only, not chemical data.
    for (auto kind : {th::SourceKind::literature, th::SourceKind::database,
                      th::SourceKind::user_supplied}) {
        auto metadata = source("metadata-schema-probe");
        metadata.kind = kind;
        th::detail::validate_source(metadata, "test-metadata", th::DataPolicy::ordinary);
    }
}
void model_revision() {
    for (std::string model : {"", "PR", "PR78", "SW", "CPA"}) {
        auto f = fixture(); f.data.model_id = std::move(model);
        expect_error(Code::unsupported_model, [&] { (void)f.build(); });
    }
    auto f = fixture(); f.data.dataset_id.clear();
    expect_error(Code::missing_field, [&] { (void)f.build(); });
    f = fixture(); f.data.revision.clear();
    expect_error(Code::missing_field, [&] { (void)f.build(); });
    f = fixture(); const auto old = f.build();
    f.data.revision = "test-v2";
    const auto next = f.build();
    require(old.revision() == "test-v1" && next.revision() == "test-v2", "revision identity");
    require(old.model_id() == th::pr76_profile, "model profile identity");
}
void applicability() {
    auto f = fixture();
    const auto data = f.build();
    using Assessment = th::RangeAssessment;
    require(data.applicability().assess(250, 1e5) == Assessment::inside_declared_bounds,
            "lower bound");
    require(data.applicability().assess(450, 5e6) == Assessment::inside_declared_bounds,
            "upper bound");
    require(data.applicability().assess(249, 1e5) == Assessment::outside_declared_bounds,
            "outside T");
    require(data.applicability().assess(350, 6e6) == Assessment::outside_declared_bounds,
            "outside p");
    f.data.applicability.pressure_pa.reset();
    const auto partial = f.build();
    require(partial.applicability().assess(350, 1e7) == Assessment::unknown,
            "unknown is not unlimited");
    require(partial.applicability().assess(249, 1e7) == Assessment::outside_declared_bounds,
            "known outside overrides missing other bound");
    f.data.applicability.temperature_k.reset();
    const auto unknown = f.build();
    require(unknown.applicability().assess(350, 1e5) == Assessment::unknown, "missing bounds");
    for (const auto range : {th::ClosedInterval{400, 300}, th::ClosedInterval{0, 300},
                            th::ClosedInterval{250, std::numeric_limits<double>::infinity()}}) {
        f.data.applicability.temperature_k = range;
        expect_error(Code::invalid_range, [&] { (void)f.build(); });
    }
    expect_error(Code::invalid_value, [&] { (void)data.applicability().assess(0, 1e5); });
}
void ownership_and_recovery() {
    auto f = fixture();
    const auto initial = f.build();
    f.catalog[0].display_name = "changed label";
    f.data.pure[0].critical_temperature->value = -1;
    f.order = {"test:d", "test:a"};
    expect_error(Code::invalid_value, [&] { (void)f.build(); });
    require(initial.components().at(0).display_name == "Synthetic label", "catalog alias");
    require(initial.critical_temperatures_k()[2] == 400, "source input alias");
    const auto copy = initial;
    f = fixture();
    const auto restored = f.build();
    check_binding(initial, f.order); check_binding(copy, f.order); check_binding(restored, f.order);
    bool caught = false;
    try { (void)restored.kij(4, 0); } catch (const std::out_of_range&) { caught = true; }
    require(caught, "row bounds check");
    caught = false;
    try { (void)restored.kij(0, 4); } catch (const std::out_of_range&) { caught = true; }
    require(caught, "column bounds check");
    // Returning from the builder must not retain spans/pointers into a temporary dataset.
    const auto temporary_input = fixture().build();
    check_binding(temporary_input, f.order);
}
void resource_limits() {
    auto f = fixture();
    expect_error(Code::size_limit, [&] { (void)f.build({3, 6, 16}); });
    expect_error(Code::size_limit, [&] { (void)f.build({4, 5, 16}); });
    expect_error(Code::size_limit, [&] { (void)f.build({4, 6, 15}); });
    const auto exact = f.build({4, 6, 16});
    require(exact.kij_matrix().size() == 16, "inclusive quotas");
    f.order.clear();
    expect_error(Code::missing_field, [&] { (void)f.build(); });
}

template <typename T>
concept TemporaryItems = requires(T value) { std::move(value).items(); };
template <typename T>
concept TemporaryMatrix = requires(T value) { std::move(value).kij_matrix(); };
static_assert(!TemporaryItems<th::OrderedComponents>);
static_assert(!TemporaryMatrix<th::PrParameterSet>);
static_assert(!std::is_default_constructible_v<th::PrParameterSet>);
static_assert(!std::is_copy_assignable_v<th::PrParameterSet>);
static_assert(std::is_copy_constructible_v<th::PrParameterSet>);
static_assert(std::is_same_v<decltype(std::declval<const th::PrParameterSet&>().kij_matrix()),
                             std::span<const double>>);
void headers() {
    require(components_header_compiles(), "common header standalone integration");
    require(pr_profile_from_separate_translation_unit() == th::pr76_profile, "PR header ODR");
}
} // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "provide one named contract test");
        const std::string_view name{argv[1]};
        if (name == "ordered_identity") { ordered_identity(); }
        else if (name == "order_mutations") { order_mutations(); }
        else if (name == "all_permutations") { all_permutations(); }
        else if (name == "identity_errors") { identity_errors(); }
        else if (name == "missing_parameters") { missing_parameters(); }
        else if (name == "duplicate_and_pair_errors") { duplicate_and_pair_errors(); }
        else if (name == "values_and_units") { values_and_units(); }
        else if (name == "provenance_policy") { provenance_policy(); }
        else if (name == "model_revision") { model_revision(); }
        else if (name == "applicability") { applicability(); }
        else if (name == "ownership_and_recovery") { ownership_and_recovery(); }
        else if (name == "resource_limits") { resource_limits(); }
        else if (name == "headers") { headers(); }
        else { throw std::invalid_argument("unknown contract test"); }
        std::cout << "[PASS] " << name << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
