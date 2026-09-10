#include <mpmc/flash/sw92_asymmetric_max2.hpp>

#include "test_support.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " +
                                 std::string(message));
    }
}

void near(double actual, long double expected, long double relative = 2e-8L,
          long double absolute = 2e-11L,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) ||
        std::abs(static_cast<long double>(actual) - expected) >
            absolute + relative * std::abs(expected)) {
        std::cerr << "actual=" << actual << " expected=" << expected << '\n';
        require(false, "reference mismatch", where);
    }
}

template<class Error, class Function>
void expect_error(Function&& function) {
    bool caught = false;
    try {
        function();
    } catch (const Error&) {
        caught = true;
    }
    require(caught, "expected exception missing");
}

th::Provenance fateen_ch4_co2_source() {
    return {
        th::SourceKind::literature,
        "doi:10.1016/j.jare.2012.03.004",
        "Journal of Advanced Research 4 (2013) 137-145",
        "Table 1 row 34 (Methane/carbon dioxide), constant kij=0.0919; Fig. 4",
        "External PR + classical van-der-Waals gas/gas BIP supplied identically to AQ and NA slots; not an SW92 water-pair parameter",
        "Open-access PMCID:PMC4195455; datum transcribed from the published table",
        "CC BY-NC-ND article; numeric datum cited for validation fixture only; article content is not redistributed"};
}

sw92_test::Prepared physical_ternary_input(bool reverse = false) {
    sw92_test::Prepared prepared;
    prepared.input.model_id = std::string(th::sw92_corrected_profile);
    prepared.input.dataset_id = "SW92-corrected+Fateen2013-CH4-CO2-kij0919";
    prepared.input.revision = "SW92-pdf-cb5b1d50+Fateen2013-Table1-row34";
    prepared.input.applicability = {
        {std::nullopt, std::nullopt,
         sw92_test::paper("SW92 Eqs. (9)-(17), Tables 2-5; hybrid gas/gas BIP is separately sourced")},
        std::nullopt};

    sw92_test::add_pure(prepared, sw92_test::methane);
    sw92_test::add_pure(prepared, sw92_test::co2);
    sw92_test::add_pure(prepared, sw92_test::water);

    prepared.input.water_binary =
        sw92_test::binary_input(sw92_test::methane).input.water_binary;
    prepared.input.water_binary.push_back(
        sw92_test::binary_input(sw92_test::co2).input.water_binary.front());

    const auto gas_pair_source = fateen_ch4_co2_source();
    prepared.input.nonwater_binary.push_back({
        sw92_test::methane.id,
        sw92_test::co2.id,
        sw92_test::scalar(0.0919, th::Unit::dimensionless, gas_pair_source),
        sw92_test::scalar(0.0919, th::Unit::dimensionless, gas_pair_source)});

    prepared.order = reverse
        ? std::vector<std::string>{sw92_test::water.id,
                                   sw92_test::co2.id,
                                   sw92_test::methane.id}
        : std::vector<std::string>{sw92_test::methane.id,
                                   sw92_test::co2.id,
                                   sw92_test::water.id};
    return prepared;
}

th::Sw92ParameterSet physical_ternary_parameters(bool reverse = false) {
    auto prepared = physical_ternary_input(reverse);
    return th::Sw92ParameterSet::create(
        prepared.catalog, prepared.order, prepared.input);
}

th::Sw92Phase<double> physical_ternary_model(bool reverse = false) {
    return th::Sw92Phase<double>::from_parameters(
        physical_ternary_parameters(reverse));
}

struct TernaryGolden {
    long double water_rich_ch4;
    long double water_rich_co2;
    long double water_rich_water;
    long double gas_rich_ch4;
    long double gas_rich_co2;
    long double gas_rich_water;
    long double gas_rich_fraction;
    long double water_rich_z;
    long double gas_rich_z;
    long double feed_aq_gibbs;
    long double feed_na_gibbs;
    long double feed_aq_minus_na_gibbs;
    long double pair_gibbs;
    long double pair_minus_feed;
    long double common_ch4;
    long double common_co2;
    long double common_water;
    long double water_rich_aq_minus_na_gibbs;
    long double gas_rich_aq_minus_na_gibbs;
};

constexpr TernaryGolden golden{
    0.000926595565135305376744149839868598720L,
    0.00410222546705065601534780964593278797L,
    0.994971178967814038607908040514198613L,
    0.694825377331009477598127996194254247L,
    0.294122280631855616926095936126350533L,
    0.0110523420371349054757760676793952200L,
    0.503060984696532100395670072167533754L,
    0.0758151675933028182659305846173684704L,
    0.860286072014983298778271289543848281L,
    -2.28631878230116056890700781726502666L,
    -1.59857965464664726953666024786125066L,
    -0.687739127654513299370347569403776002L,
    -3.09595922989048202539397969359502534L,
    -0.809640447589321456486971876329998672L,
    -0.462255240302127022411452093427274697L,
    -1.50007205423714639252586253581859718L,
    -5.41831817529833121734218416104537923L,
    -0.0158605870907434694828533831379943984L,
    -0.00607565267345640010750747291986237355L,
};

struct OrientedPair {
    const fl::Sw92AsymmetricFixedPairPhase* water_rich{};
    const fl::Sw92AsymmetricFixedPairPhase* gas_rich{};
    const fl::Sw92AsymmetricFamilyAssignmentCheck* water_rich_assignment{};
    const fl::Sw92AsymmetricFamilyAssignmentCheck* gas_rich_assignment{};
};

OrientedPair orient_pair(const fl::Sw92AsymmetricFixedPairState& pair,
                         std::size_t water_index) {
    if (pair.phase0.composition.at(water_index) >=
        pair.phase1.composition.at(water_index)) {
        return {&pair.phase0, &pair.phase1,
                &pair.phase0_assignment, &pair.phase1_assignment};
    }
    return {&pair.phase1, &pair.phase0,
            &pair.phase1_assignment, &pair.phase0_assignment};
}

void traceable_parameter_snapshot() {
    const auto parameters = physical_ternary_parameters();
    require(parameters.components().size() == 3 && parameters.water_index() == 2,
            "traceable ternary component/water indexing changed");
    require(parameters.components().items()[0].id == sw92_test::methane.id &&
                parameters.components().items()[1].id == sw92_test::co2.id &&
                parameters.components().items()[2].id == sw92_test::water.id,
            "traceable ternary ordered identities changed");

    const auto records = parameters.nonwater_binary_records();
    require(records.size() == 1 && records.front().aqueous_kij &&
                records.front().nonaqueous_kij,
            "traceable CH4/CO2 non-water pair missing");
    const auto& record = records.front();
    near(record.aqueous_kij->value, 0.0919L, 0.0L, 1e-15L);
    near(record.nonaqueous_kij->value, 0.0919L, 0.0L, 1e-15L);
    for (const auto* datum : {&*record.aqueous_kij, &*record.nonaqueous_kij}) {
        require(datum->source.kind == th::SourceKind::literature &&
                    datum->source.reference == "doi:10.1016/j.jare.2012.03.004" &&
                    datum->source.locator.find("Table 1 row 34") != std::string::npos,
                "external CH4/CO2 BIP provenance was not retained");
    }
    require(record.aqueous_kij->source.reference !=
                sw92_test::paper("Table 5 non-aqueous water BIP").reference,
            "external gas/gas BIP was incorrectly labeled as an SW92 datum");
}

void ternary_max2_reference() {
    const auto model = physical_ternary_model();
    const Vec feed{0.35, 0.15, 0.50};
    const auto result = fl::solve_sw92_xu_asymmetric_max2(
        1.0e7, 350.0, feed, model, 0.0);

    require(result.status ==
                fl::Sw92AsymmetricMax2Status::two_phase_no_instability_found &&
                result.accepted_phase_set() != nullptr &&
                result.accepted_phase_count() == 2 &&
                result.final_stability &&
                result.final_stability->status == fl::StabilityStatus::no_instability_found &&
                result.selection.initial_stability.reference_family ==
                    th::SwPhaseFamily::aqueous &&
                !result.global_stability_proven,
            "traceable ternary max2 state was not accepted by the finite two-family review");

    const auto* pair =
        result.selection.selected_pair_candidate_pending_final_stability();
    require(pair != nullptr, "traceable ternary selected pair missing");
    const auto oriented = orient_pair(*pair, 2U);
    const auto& water_rich = *oriented.water_rich;
    const auto& gas_rich = *oriented.gas_rich;
    require(water_rich.family == th::SwPhaseFamily::aqueous &&
                gas_rich.family == th::SwPhaseFamily::aqueous &&
                oriented.water_rich_assignment->status ==
                    fl::Sw92AsymmetricFamilyAssignmentStatus::assigned_lower &&
                oriented.gas_rich_assignment->status ==
                    fl::Sw92AsymmetricFamilyAssignmentStatus::assigned_lower,
            "traceable ternary lower-envelope family assignment changed");

    near(water_rich.composition[0], golden.water_rich_ch4);
    near(water_rich.composition[1], golden.water_rich_co2);
    near(water_rich.composition[2], golden.water_rich_water);
    near(gas_rich.composition[0], golden.gas_rich_ch4);
    near(gas_rich.composition[1], golden.gas_rich_co2);
    near(gas_rich.composition[2], golden.gas_rich_water);
    near(gas_rich.mole_phase_fraction, golden.gas_rich_fraction);
    near(water_rich.compressibility_factor, golden.water_rich_z);
    near(gas_rich.compressibility_factor, golden.gas_rich_z);

    near(result.lower_feed_reduced_gibbs, golden.feed_aq_gibbs);
    near(result.selection.initial_stability.nonaqueous.feed_reduced_gibbs,
         golden.feed_na_gibbs);
    near(result.selection.initial_stability.aqueous_minus_nonaqueous_feed_gibbs,
         golden.feed_aq_minus_na_gibbs);
    near(result.selected_pair_reduced_gibbs, golden.pair_gibbs);
    near(result.pair_minus_lower_feed_reduced_gibbs, golden.pair_minus_feed);
    require(result.pair_minus_lower_feed_reduced_gibbs <
                -result.pair_feed_gibbs_combined_guard,
            "traceable ternary pair is not resolved below the feed Gibbs state");

    require(result.final_stability->common_log_activity.size() == 3,
            "traceable ternary final common tangent dimension changed");
    near(result.final_stability->common_log_activity[0], golden.common_ch4);
    near(result.final_stability->common_log_activity[1], golden.common_co2);
    near(result.final_stability->common_log_activity[2], golden.common_water);
    near(oriented.water_rich_assignment->aqueous_minus_nonaqueous_gibbs,
         golden.water_rich_aq_minus_na_gibbs);
    near(oriented.gas_rich_assignment->aqueous_minus_nonaqueous_gibbs,
         golden.gas_rich_aq_minus_na_gibbs);

    require(pair->chemical_potential_residual.size() == 3,
            "traceable ternary chemical-potential residual dimension changed");
    for (std::size_t i = 0; i < feed.size(); ++i) {
        const double recovered =
            pair->phase0.mole_phase_fraction * pair->phase0.composition[i] +
            pair->phase1.mole_phase_fraction * pair->phase1.composition[i];
        near(recovered, feed[i], 0.0L, 2e-10L);
        require(std::abs(pair->chemical_potential_residual[i]) <=
                    result.options.selection.fixed_pair.chemical_potential_tolerance,
                "traceable ternary common chemical potential tolerance failed");
    }
}

void component_permutation() {
    const auto normal_model = physical_ternary_model(false);
    const auto reverse_model = physical_ternary_model(true);
    const auto normal = fl::solve_sw92_xu_asymmetric_max2(
        1.0e7, 350.0, Vec{0.35, 0.15, 0.50}, normal_model, 0.0);
    const auto reverse = fl::solve_sw92_xu_asymmetric_max2(
        1.0e7, 350.0, Vec{0.50, 0.15, 0.35}, reverse_model, 0.0);
    require(normal.status == fl::Sw92AsymmetricMax2Status::two_phase_no_instability_found &&
                reverse.status == normal.status &&
                normal.accepted_phase_count() == 2 && reverse.accepted_phase_count() == 2,
            "traceable ternary permutation changed max2 acceptance");

    const auto* first_pair =
        normal.selection.selected_pair_candidate_pending_final_stability();
    const auto* second_pair =
        reverse.selection.selected_pair_candidate_pending_final_stability();
    require(first_pair != nullptr && second_pair != nullptr,
            "traceable ternary permuted pair missing");
    const auto first = orient_pair(*first_pair, 2U);
    const auto second = orient_pair(*second_pair, 0U);

    near(first.water_rich->mole_phase_fraction,
         second.water_rich->mole_phase_fraction, 2e-8L, 2e-11L);
    near(first.gas_rich->mole_phase_fraction,
         second.gas_rich->mole_phase_fraction, 2e-8L, 2e-11L);
    near(first.water_rich->composition[0], second.water_rich->composition[2]);
    near(first.water_rich->composition[1], second.water_rich->composition[1]);
    near(first.water_rich->composition[2], second.water_rich->composition[0]);
    near(first.gas_rich->composition[0], second.gas_rich->composition[2]);
    near(first.gas_rich->composition[1], second.gas_rich->composition[1]);
    near(first.gas_rich->composition[2], second.gas_rich->composition[0]);
    near(normal.selected_pair_reduced_gibbs, reverse.selected_pair_reduced_gibbs,
         2e-8L, 2e-11L);
    require(normal.component_ids.size() == 3 && reverse.component_ids.size() == 3 &&
                normal.component_ids[0] == reverse.component_ids[2] &&
                normal.component_ids[1] == reverse.component_ids[1] &&
                normal.component_ids[2] == reverse.component_ids[0],
            "traceable ternary component metadata did not follow the requested order");
}

void missing_nonwater_pair_rejected() {
    auto missing = physical_ternary_input();
    missing.input.nonwater_binary.clear();
    expect_error<th::ContractError>([&] {
        (void)th::Sw92ParameterSet::create(
            missing.catalog, missing.order, missing.input);
    });

    auto missing_family_value = physical_ternary_input();
    missing_family_value.input.nonwater_binary.front().nonaqueous_kij.reset();
    expect_error<th::ContractError>([&] {
        (void)th::Sw92ParameterSet::create(
            missing_family_value.catalog, missing_family_value.order,
            missing_family_value.input);
    });
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test cases[] = {
    {"traceable_parameter_snapshot", traceable_parameter_snapshot},
    {"ternary_max2_reference", ternary_max2_reference},
    {"component_permutation", component_permutation},
    {"missing_nonwater_pair_rejected", missing_nonwater_pair_rejected},
};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument("exactly one test name required");
        }
        for (const auto& [name, run] : cases) {
            if (name == argv[1]) {
                run();
                std::cout << "[PASS] " << name << '\n';
                return 0;
            }
        }
        throw std::invalid_argument("unknown test name");
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
