#include <mpmc/flash/sw92_asymmetric_max2.hpp>

#include "test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
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

void near(double actual, long double expected, long double relative = 8e-9L,
          long double absolute = 8e-12L,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) ||
        std::abs(static_cast<long double>(actual) - expected) >
            absolute + relative * std::abs(expected)) {
        std::cerr << "actual=" << actual << " expected=" << expected << '\n';
        require(false, "reference mismatch", where);
    }
}

struct BinaryGolden {
    long double low_gas;
    long double high_gas;
    long double high_fraction;
    long double low_z;
    long double high_z;
    long double lower_feed_gibbs;
    long double pair_gibbs;
    long double pair_minus_feed;
    long double common_d_gas;
    long double common_d_water;
    long double low_aq_minus_na_gibbs;
    long double high_aq_minus_na_gibbs;
    long double aq_na_low_gas;
    long double aq_na_high_gas;
};

constexpr BinaryGolden golden[] = {
    {
        0.0012488516673043369772690695862810869L,
        0.98976716870138545359588622097568157L,
        0.50454416447146141705145186641410852L,
        0.075641724396787743253196769889981955L,
        0.90395573410374091463758399061390327L,
        -1.7966801075584243682333600397004459L,
        -2.7666079361625361393408003803537260L,
        -0.96992782860411177110744034065328018L,
        -0.11855782477689377854625233874433222L,
        -5.4146580475481785001353484219631199L,
        -0.0065287747292147990593988720606446854L,
        -0.0062801265280233769480425833604655190L,
        0.0012548893978171366881230600626601386L,
        0.99448254517749047477840934275427544L,
    },
    {
        0.00057927917850318821367216301214118614L,
        0.98563843743740236519943874461280656L,
        0.50699566278255549394787817509225852L,
        0.14290115194351380537689004320025247L,
        0.92414076037139241686440971190588936L,
        -1.7674270686613408746172242606525397L,
        -2.6750517768359819185769968256602970L,
        -0.90762470817464104395977256500775723L,
        -0.14172250487679300817885858011795890L,
        -5.2083810487951708289751350712026350L,
        -0.0017962814315828905309606338801866296L,
        -0.0090032749539800564092525982378769837L,
        0.00058324135361841879751362867136907871L,
        0.99242009503855451670813326717549818L,
    },
};

th::Sw92Phase<double> methane_model(bool reverse = false) {
    return th::Sw92Phase<double>::from_parameters(
        sw92_test::binary_parameters(sw92_test::methane, reverse));
}

th::Sw92Phase<double> co2_model() {
    return th::Sw92Phase<double>::from_parameters(
        sw92_test::binary_parameters(sw92_test::co2));
}

th::Sw92Phase<double> methane_family_control_model(double nonaqueous_kij) {
    auto prepared = sw92_test::binary_input(sw92_test::methane);
    const auto source = sw92_test::synthetic(
        "controlled CH4/water NA BIP for family-envelope switch/tie routing only");
    prepared.input.dataset_id = "SW92-family-switch-structural";
    prepared.input.revision = "synthetic-family-envelope-v1";
    prepared.input.water_binary.front().nonaqueous_kij = sw92_test::scalar(
        nonaqueous_kij, th::Unit::dimensionless, source);
    return th::Sw92Phase<double>::from_parameters(
        th::Sw92ParameterSet::create(
            prepared.catalog, prepared.order, prepared.input,
            th::DataPolicy::allow_synthetic_tests));
}

th::Sw92Phase<double> ternary_structural_model(bool reverse = false) {
    sw92_test::Prepared prepared;
    prepared.input.model_id = std::string(th::sw92_corrected_profile);
    prepared.input.dataset_id = "SW92-runtime-ternary-structural";
    prepared.input.revision = "synthetic-nonwater-bip-v1";
    prepared.input.applicability = {
        {std::nullopt, std::nullopt,
         sw92_test::paper("SW92 model basis; synthetic CH4/N2 non-water pair only")},
        std::nullopt};
    sw92_test::add_pure(prepared, sw92_test::methane);
    sw92_test::add_pure(prepared, sw92_test::nitrogen);
    sw92_test::add_pure(prepared, sw92_test::water);
    prepared.input.water_binary =
        sw92_test::binary_input(sw92_test::methane).input.water_binary;
    prepared.input.water_binary.push_back(
        sw92_test::binary_input(sw92_test::nitrogen).input.water_binary.front());
    const auto synthetic = sw92_test::synthetic(
        "CH4/N2 kij=0 in both families; runtime-dimension structure only, not physical validation");
    prepared.input.nonwater_binary.push_back({
        sw92_test::methane.id, sw92_test::nitrogen.id,
        sw92_test::scalar(0.0, th::Unit::dimensionless, synthetic),
        sw92_test::scalar(0.0, th::Unit::dimensionless, synthetic)});
    prepared.order = reverse
        ? std::vector<std::string>{sw92_test::water.id,
                                   sw92_test::nitrogen.id,
                                   sw92_test::methane.id}
        : std::vector<std::string>{sw92_test::methane.id,
                                   sw92_test::nitrogen.id,
                                   sw92_test::water.id};
    return th::Sw92Phase<double>::from_parameters(
        th::Sw92ParameterSet::create(
            prepared.catalog, prepared.order, prepared.input,
            th::DataPolicy::allow_synthetic_tests));
}

const fl::Sw92AsymmetricFixedPairState& selected_pair(
    const fl::Sw92AsymmetricMax2Result& result) {
    const auto* pair = result.selection.selected_pair_candidate_pending_final_stability();
    require(pair != nullptr, "selected Gate-3B.2 pair missing");
    return *pair;
}

const fl::Sw92AsymmetricFixedPairPhase& low_gas_phase(
    const fl::Sw92AsymmetricFixedPairState& pair, std::size_t gas_index) {
    return pair.phase0.composition[gas_index] < pair.phase1.composition[gas_index]
        ? pair.phase0 : pair.phase1;
}

const fl::Sw92AsymmetricFixedPairPhase& high_gas_phase(
    const fl::Sw92AsymmetricFixedPairState& pair, std::size_t gas_index) {
    return pair.phase0.composition[gas_index] < pair.phase1.composition[gas_index]
        ? pair.phase1 : pair.phase0;
}

void check_binary_reference(const fl::Sw92AsymmetricMax2Result& result,
                            const BinaryGolden& reference,
                            std::size_t gas_index, std::size_t water_index) {
    require(result.status == fl::Sw92AsymmetricMax2Status::two_phase_no_instability_found &&
                result.accepted_phase_set() != nullptr &&
                result.accepted_phase_count() == 2 &&
                result.final_stability &&
                result.final_stability->status == fl::StabilityStatus::no_instability_found &&
                result.selection.initial_stability.reference_family ==
                    th::SwPhaseFamily::aqueous &&
                !result.global_stability_proven,
            "binary max2 reference was not accepted under finite two-family review");
    const auto& pair = selected_pair(result);
    const auto& low = low_gas_phase(pair, gas_index);
    const auto& high = high_gas_phase(pair, gas_index);
    require(low.family == th::SwPhaseFamily::aqueous &&
                high.family == th::SwPhaseFamily::aqueous,
            "Xu lower-envelope binary reference no longer selects AQ+AQ");
    near(low.composition[gas_index], reference.low_gas);
    near(high.composition[gas_index], reference.high_gas);
    near(high.mole_phase_fraction, reference.high_fraction);
    require(low.compressibility_factor > 0.0 && high.compressibility_factor > 0.0,
            "binary reference lost Z diagnostics");
    near(low.compressibility_factor, reference.low_z);
    near(high.compressibility_factor, reference.high_z);
    near(result.lower_feed_reduced_gibbs, reference.lower_feed_gibbs);
    near(result.pair_recomputed_reduced_gibbs, reference.pair_gibbs);
    near(result.pair_minus_lower_feed_reduced_gibbs, reference.pair_minus_feed);
    require(result.final_stability->common_log_activity.size() == 2,
            "binary final common tangent dimension changed");
    near(result.final_stability->common_log_activity[gas_index], reference.common_d_gas);
    near(result.final_stability->common_log_activity[water_index], reference.common_d_water);
    near(low.family == pair.phase0.family
             ? pair.phase0_assignment.aqueous_minus_nonaqueous_gibbs
             : pair.phase1_assignment.aqueous_minus_nonaqueous_gibbs,
         reference.low_aq_minus_na_gibbs);
    near(high.family == pair.phase0.family &&
                 high.composition.data() == pair.phase0.composition.data()
             ? pair.phase0_assignment.aqueous_minus_nonaqueous_gibbs
             : pair.phase1_assignment.aqueous_minus_nonaqueous_gibbs,
         reference.high_aq_minus_na_gibbs);
}

void ch4_fresh_max2() {
    const auto model = methane_model();
    const auto result = fl::solve_sw92_xu_asymmetric_max2(
        1.0e7, 350.0, Vec{0.5, 0.5}, model, 0.0);
    check_binary_reference(result, golden[0], 0U, 1U);
}

void ch4_brine_max2() {
    const auto model = methane_model();
    const auto result = fl::solve_sw92_xu_asymmetric_max2(
        2.0e7, 376.15, Vec{0.5, 0.5}, model, 4.0);
    check_binary_reference(result, golden[1], 0U, 1U);
    require(result.nacl_molality_mol_per_kg_water == 4.0 &&
                result.selection.nacl_molality_mol_per_kg_water == 4.0,
            "nonzero fixed molality was not retained through max2");
}

void original_assignment_profile_separation() {
    const auto model = methane_model();
    for (std::size_t index = 0; index < 2; ++index) {
        const double pressure = index == 0 ? 1.0e7 : 2.0e7;
        const double temperature = index == 0 ? 350.0 : 376.15;
        const double molality = index == 0 ? 0.0 : 4.0;
        const auto& reference = golden[index];
        const double x = static_cast<double>(reference.aq_na_low_gas);
        const double y = static_cast<double>(reference.aq_na_high_gas);
        const Vec log_k{
            std::log(y / x),
            std::log((1.0 - y) / (1.0 - x))};
        const auto fixed = fl::iterate_sw92_asymmetric_fixed_pair(
            pressure, temperature, Vec{0.5, 0.5}, log_k,
            {th::SwPhaseFamily::aqueous, th::SwPhaseFamily::nonaqueous},
            model, molality);
        require(fixed.equations_converged() && fixed.point &&
                    fixed.status ==
                        fl::Sw92AsymmetricFixedPairStatus::family_assignment_dominated &&
                    !fixed.candidate_admissible(),
                "phase-assigned AQ+NA equations no longer remain distinct from Xu lower-envelope acceptance");
        near(fixed.point->phase0.composition[0], reference.aq_na_low_gas);
        near(fixed.point->phase1.composition[0], reference.aq_na_high_gas);
        require(fixed.point->phase1_assignment.status ==
                    fl::Sw92AsymmetricFamilyAssignmentStatus::assigned_dominated,
                "AQ+NA profile-separation reference lost NA dominance diagnostic");
    }
}

void runtime_ternary_structural() {
    const auto normal_model = ternary_structural_model(false);
    const auto reverse_model = ternary_structural_model(true);
    const auto normal = fl::solve_sw92_xu_asymmetric_max2(
        1.0e7, 350.0, Vec{0.49, 0.01, 0.50}, normal_model, 1.0);
    const auto reverse = fl::solve_sw92_xu_asymmetric_max2(
        1.0e7, 350.0, Vec{0.50, 0.01, 0.49}, reverse_model, 1.0);
    require(normal.status == fl::Sw92AsymmetricMax2Status::two_phase_no_instability_found &&
                reverse.status == normal.status &&
                normal.accepted_phase_count() == 2 && reverse.accepted_phase_count() == 2,
            "runtime ternary structural snapshot did not complete max2 on both permutations");
    const auto& a = selected_pair(normal);
    const auto& b = selected_pair(reverse);
    require(a.phase0.family == th::SwPhaseFamily::aqueous &&
                a.phase1.family == th::SwPhaseFamily::aqueous &&
                b.phase0.family == th::SwPhaseFamily::aqueous &&
                b.phase1.family == th::SwPhaseFamily::aqueous,
            "runtime ternary structural family identity changed");

    const auto& a_low = low_gas_phase(a, 0U);
    const auto& a_high = high_gas_phase(a, 0U);
    const auto& b_low = low_gas_phase(b, 2U);
    const auto& b_high = high_gas_phase(b, 2U);
    near(a_low.mole_phase_fraction, b_low.mole_phase_fraction, 2e-8L, 2e-11L);
    near(a_high.mole_phase_fraction, b_high.mole_phase_fraction, 2e-8L, 2e-11L);
    for (const auto* phase : {&a_low, &a_high}) {
        require(phase->composition.size() == 3 &&
                    phase->composition[0] > 0.0 &&
                    phase->composition[1] > 0.0 &&
                    phase->composition[2] > 0.0,
                "runtime ternary active support was lost");
    }
    near(a_low.composition[0], b_low.composition[2], 2e-8L, 2e-11L);
    near(a_low.composition[1], b_low.composition[1], 2e-8L, 2e-11L);
    near(a_low.composition[2], b_low.composition[0], 2e-8L, 2e-11L);
    near(a_high.composition[0], b_high.composition[2], 2e-8L, 2e-11L);
    near(a_high.composition[1], b_high.composition[1], 2e-8L, 2e-11L);
    near(a_high.composition[2], b_high.composition[0], 2e-8L, 2e-11L);
}

void phase_disappearance_boundary() {
    const auto model = methane_model();
    constexpr double beta = 0.05;
    const double x = static_cast<double>(golden[0].low_gas);
    const double y = static_cast<double>(golden[0].high_gas);
    const double feed_gas = std::fma(beta, y - x, x);
    fl::Sw92AsymmetricMax2Options options;
    options.selection.fixed_pair.minimum_phase_fraction = 0.10;
    const auto result = fl::solve_sw92_xu_asymmetric_max2(
        1.0e7, 350.0, Vec{feed_gas, 1.0 - feed_gas}, model, 0.0, options);
    require(result.status == fl::Sw92AsymmetricMax2Status::indeterminate &&
                result.accepted_phase_set() == nullptr &&
                !result.final_stability,
            "numerical phase-disappearance threshold incorrectly published a max2 state");
    bool saw_disappearance = false;
    for (const auto& attempt : result.selection.attempts) {
        if (attempt.fixed_pair &&
            attempt.fixed_pair->status == fl::Sw92AsymmetricFixedPairStatus::phase_disappearance) {
            require(attempt.fixed_pair->equations_converged(),
                    "phase disappearance was reported before equation/balance convergence");
            saw_disappearance = true;
        }
    }
    require(saw_disappearance,
            "traceable coexistence-derived feed did not reach the configured disappearance guard");
}

double family_gap(const th::Sw92Phase<double>& model,
                  double pressure, double temperature, double molality,
                  double gas_fraction) {
    const Vec composition{gas_fraction, 1.0 - gas_fraction};
    fl::Sw92FamilyStabilityEvaluator aq(
        model, molality, th::SwPhaseFamily::aqueous);
    fl::Sw92FamilyStabilityEvaluator na(
        model, molality, th::SwPhaseFamily::nonaqueous);
    const auto aq_phase = aq(pressure, temperature, composition);
    const auto na_phase = na(pressure, temperature, composition);
    return composition[0] * (aq_phase.ln_phi[0] - na_phase.ln_phi[0]) +
           composition[1] * (aq_phase.ln_phi[1] - na_phase.ln_phi[1]);
}

void real_binary_family_envelope_audit() {
    const auto model = methane_model();
    constexpr std::array<double, 9> interior{
        0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99};
    for (const auto [pressure, temperature, molality] :
         {std::array<double, 3>{1.0e7, 350.0, 0.0},
          std::array<double, 3>{2.0e7, 376.15, 4.0}}) {
        for (double gas_fraction : interior) {
            require(family_gap(model, pressure, temperature, molality, gas_fraction) < 0.0,
                    "sampled physical CH4 binary unexpectedly switched lower family away from AQ");
        }
        require(std::abs(family_gap(model, pressure, temperature, molality, 0.0)) < 2e-12 &&
                    std::abs(family_gap(model, pressure, temperature, molality, 1.0)) < 2e-12,
                "pure-component family surfaces did not meet at the simplex vertices");
    }
}

void synthetic_family_switch_guard() {
    const double tr = 350.0 / sw92_test::methane.tc;
    const double aq_kij = th::sw92_aqueous_water_kij(
        th::Sw92Species::hydrocarbon, tr, sw92_test::methane.omega, 0.0);
    const auto aq_lower_model = methane_family_control_model(aq_kij + 0.02);
    const auto na_lower_model = methane_family_control_model(aq_kij - 0.02);
    const auto tie_model = methane_family_control_model(aq_kij);
    const Vec feed{0.5, 0.5};

    const auto aq_lower = fl::test_sw92_pt_asymmetric_stability(
        1.0e7, 350.0, feed, aq_lower_model, 0.0);
    const auto na_lower = fl::test_sw92_pt_asymmetric_stability(
        1.0e7, 350.0, feed, na_lower_model, 0.0);
    const auto tie = fl::solve_sw92_xu_asymmetric_max2(
        1.0e7, 350.0, feed, tie_model, 0.0);
    require(aq_lower.feed_reference_status ==
                fl::Sw92AsymmetricFeedReferenceStatus::selected &&
                aq_lower.reference_family == th::SwPhaseFamily::aqueous,
            "controlled family switch did not select AQ on the AQ-lower side");
    require(na_lower.feed_reference_status ==
                fl::Sw92AsymmetricFeedReferenceStatus::selected &&
                na_lower.reference_family == th::SwPhaseFamily::nonaqueous,
            "controlled family switch did not select NA on the NA-lower side");
    require(tie.status == fl::Sw92AsymmetricMax2Status::indeterminate &&
                tie.selection.initial_stability.feed_reference_status ==
                    fl::Sw92AsymmetricFeedReferenceStatus::family_gibbs_tie &&
                tie.accepted_phase_set() == nullptr,
            "exact controlled family-envelope tie was broken by max2 execution order");
}

void root_envelope_nonsmooth() {
    constexpr double co2_pr_saturation_pressure_pa = 4148774.826264025;
    const auto model = co2_model();
    const auto result = fl::solve_sw92_xu_asymmetric_max2(
        co2_pr_saturation_pressure_pa, 280.0, Vec{1.0, 0.0}, model, 0.0);
    require(result.status == fl::Sw92AsymmetricMax2Status::indeterminate &&
                result.selection.initial_stability.feed_reference_status ==
                    fl::Sw92AsymmetricFeedReferenceStatus::family_reference_nonsmooth &&
                result.accepted_phase_set() == nullptr &&
                result.selection.plan.empty() && !result.final_stability,
            "same-family root-envelope nonsmoothness did not stop max2 publication");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test cases[] = {
    {"ch4_fresh_max2", ch4_fresh_max2},
    {"ch4_brine_max2", ch4_brine_max2},
    {"original_assignment_profile_separation", original_assignment_profile_separation},
    {"runtime_ternary_structural", runtime_ternary_structural},
    {"phase_disappearance_boundary", phase_disappearance_boundary},
    {"real_binary_family_envelope_audit", real_binary_family_envelope_audit},
    {"synthetic_family_switch_guard", synthetic_family_switch_guard},
    {"root_envelope_nonsmooth", root_envelope_nonsmooth},
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
