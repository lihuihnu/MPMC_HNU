#include <mpmc/flash/sw92_asymmetric_fixed_pair.hpp>

#include "test_support.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool sw92_asymmetric_fixed_pair_header();

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

void near(double actual, long double expected, long double relative = 5e-9L,
          long double absolute = 5e-12L,
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

th::Sw92Phase<double> binary_model(bool reverse = false) {
    return th::Sw92Phase<double>::from_parameters(
        sw92_test::binary_parameters(sw92_test::co2, reverse));
}

constexpr double aq_co2_kij_340k_fresh =
    -0.058970986957265289338095088569696338832153503254491;

th::Sw92Phase<double> family_tie_model() {
    auto prepared = sw92_test::binary_input(sw92_test::co2);
    const auto source = sw92_test::synthetic(
        "NA CO2/water BIP set equal to AQ value at T=340 K,m=0 for family-tie routing test");
    prepared.input.water_binary.front().nonaqueous_rule =
        th::Sw92NonAqueousWaterRule::constant;
    prepared.input.water_binary.front().nonaqueous_kij = sw92_test::scalar(
        aq_co2_kij_340k_fresh, th::Unit::dimensionless, source);
    return th::Sw92Phase<double>::from_parameters(
        th::Sw92ParameterSet::create(
            prepared.catalog, prepared.order, prepared.input,
            th::DataPolicy::allow_synthetic_tests));
}

struct PairGolden {
    th::SwPhaseFamily phase0_family;
    th::SwPhaseFamily phase1_family;
    long double phase0_gas;
    long double phase1_gas;
    long double phase1_fraction;
    long double phase0_z;
    long double phase1_z;
    long double reduced_gibbs;
    long double phase0_aq_minus_na_gibbs;
    long double phase1_aq_minus_na_gibbs;
};

constexpr PairGolden pair_golden[] = {
    {th::SwPhaseFamily::aqueous, th::SwPhaseFamily::aqueous,
     0.005936172024213170056428932801071659387598344017016019414L,
     0.9874123977699055456489128208133867393111524204260154840L,
     0.7071631586882918679052200595167647619178931444548134912L,
     0.02326903773317445723089281117074340993256381143281275291L,
     0.8848095701078985979188014507300051296938502934035683490L,
     -1.496025765990406698063759588412595324322460168711074814L,
     -0.01753986130158987772546498610443249129296630047686235902L,
     -0.001466097076066687493326937942517705487825113049396313454L},
    {th::SwPhaseFamily::aqueous, th::SwPhaseFamily::nonaqueous,
     0.005945065166186101440737568231518653586387892719453330789L,
     0.9887877343090553106585559885449699850347182199737051566L,
     0.7061709433505718241361429803823283207927897009254034181L,
     0.02326916182253476550316956093117009661761092895870925982L,
     0.8865753259957614197158614840436129757850208327439343804L,
     -1.495046828318048125530732605527113449469374366608903618L,
     -0.01756589035909635264304576153090361754607500997127131831L,
     -0.001307369574295229517510849155183697216677449660249782584L},
    {th::SwPhaseFamily::nonaqueous, th::SwPhaseFamily::nonaqueous,
     0.0002841484008153150292230954989823356738425387440742066945L,
     0.9887258145196700966711922877257496166779972690910427309L,
     0.7078979727216877824312072650351759172085612116009258256L,
     0.02319073809275250260990730800351579970159929677795828849L,
     0.8865684868785787959288635439178757980024354166402218082L,
     -1.493443745868457105673164419428357694050625694328953612L,
     -0.0008471334870182654101718445664860246981217476645689590207L,
     -0.001314523434951416693574700518807495286803644383638912700L}
};

Vec rough_seed(std::size_t index) {
    switch (index) {
    case 0: return {5.0, -4.2};
    case 1: return {5.0, -4.4};
    case 2: return {8.0, -4.4};
    default: throw std::out_of_range("pair golden index");
    }
}

Vec exact_log_k(const PairGolden& reference) {
    const double x = static_cast<double>(reference.phase0_gas);
    const double y = static_cast<double>(reference.phase1_gas);
    return {std::log(y / x), std::log((1.0 - y) / (1.0 - x))};
}

fl::Sw92AsymmetricFixedPairResult solve_reference_pair(
    std::size_t index, bool reverse = false, bool rough = true) {
    const auto& reference = pair_golden[index];
    const auto model = binary_model(reverse);
    Vec feed = reverse ? Vec{0.3, 0.7} : Vec{0.7, 0.3};
    Vec seed = rough ? rough_seed(index) : exact_log_k(reference);
    if (reverse) {
        std::swap(seed[0], seed[1]);
    }
    return fl::iterate_sw92_asymmetric_fixed_pair(
        3.0e6, 340.0, feed, seed,
        {reference.phase0_family, reference.phase1_family}, model, 0.0);
}

void check_point(const fl::Sw92AsymmetricFixedPairResult& result,
                 const PairGolden& reference, bool reverse = false) {
    require(result.point.has_value(), "fixed-pair result lost its converged equation point");
    const auto& point = *result.point;
    const std::size_t gas_index = reverse ? 1U : 0U;
    require(point.phase0.family == reference.phase0_family &&
                point.phase1.family == reference.phase1_family,
            "phase family identity changed");
    near(point.phase0.composition[gas_index], reference.phase0_gas);
    near(point.phase1.composition[gas_index], reference.phase1_gas);
    near(point.phase1.mole_phase_fraction, reference.phase1_fraction);
    near(point.phase0.mole_phase_fraction,
         1.0L - reference.phase1_fraction);
    near(point.phase0.compressibility_factor, reference.phase0_z);
    near(point.phase1.compressibility_factor, reference.phase1_z);
    near(point.reduced_gibbs, reference.reduced_gibbs);
    near(point.phase0_assignment.aqueous_minus_nonaqueous_gibbs,
         reference.phase0_aq_minus_na_gibbs);
    near(point.phase1_assignment.aqueous_minus_nonaqueous_gibbs,
         reference.phase1_aq_minus_na_gibbs);
    require(point.chemical_potential_norm <=
                result.options.chemical_potential_tolerance &&
                point.mass_absolute <= result.options.mass_absolute_tolerance &&
                point.mass_relative <= result.options.mass_relative_tolerance,
            "fixed-pair equations or material balance did not satisfy their gates");
    require(result.equilibrium_profile == fl::sw92_xu_asymmetric_gibbs_profile &&
                result.fixed_pair_convention ==
                    fl::sw92_xu_asymmetric_fixed_pair_convention &&
                result.model_profile == th::sw92_corrected_profile &&
                result.phase_convention == th::sw92_pt_convention,
            "fixed-pair metadata identity changed");
}

void aq_aq_fixed_pair() {
    const auto result = solve_reference_pair(0);
    require(result.status == fl::Sw92AsymmetricFixedPairStatus::converged &&
                result.equations_converged() && result.candidate_admissible(),
            "AQ+AQ fixed-pair primitive did not converge/admit its candidate");
    check_point(result, pair_golden[0]);
    require(result.point->phase0_assignment.status ==
                fl::Sw92AsymmetricFamilyAssignmentStatus::assigned_lower &&
                result.point->phase1_assignment.status ==
                    fl::Sw92AsymmetricFamilyAssignmentStatus::assigned_lower,
            "AQ+AQ lower-envelope family checks changed");
    require(result.iterations > 0 && result.gibbs_descent_steps > 0,
            "rough AQ+AQ seed did not exercise the SSI/Gibbs line-search path");
}

void aq_na_fixed_pair_dominated() {
    const auto result = solve_reference_pair(1);
    require(result.status ==
                fl::Sw92AsymmetricFixedPairStatus::family_assignment_dominated &&
                result.equations_converged() && !result.candidate_admissible(),
            "AQ+NA equation convergence/dominance semantics changed");
    check_point(result, pair_golden[1]);
    require(result.point->phase0_assignment.status ==
                fl::Sw92AsymmetricFamilyAssignmentStatus::assigned_lower &&
                result.point->phase1_assignment.status ==
                    fl::Sw92AsymmetricFamilyAssignmentStatus::assigned_dominated,
            "AQ+NA per-phase dominance classification changed");
}

void na_na_fixed_pair_dominated() {
    const auto result = solve_reference_pair(2);
    require(result.status ==
                fl::Sw92AsymmetricFixedPairStatus::family_assignment_dominated &&
                result.equations_converged() && !result.candidate_admissible(),
            "NA+NA equation convergence/dominance semantics changed");
    check_point(result, pair_golden[2]);
    require(result.point->phase0_assignment.status ==
                fl::Sw92AsymmetricFamilyAssignmentStatus::assigned_dominated &&
                result.point->phase1_assignment.status ==
                    fl::Sw92AsymmetricFamilyAssignmentStatus::assigned_dominated,
            "NA+NA per-phase dominance classification changed");
}

void imposed_helper_matches_gate3a() {
    const auto model = binary_model();
    const Vec feed{0.7, 0.3};
    const auto gate3a = fl::test_sw92_pt_asymmetric_stability(
        3.0e6, 340.0, feed, model, 0.0);
    require(gate3a.feed_reference_status ==
                fl::Sw92AsymmetricFeedReferenceStatus::selected &&
                gate3a.aqueous.search && gate3a.nonaqueous.search,
            "Gate 3A did not provide a reusable common tangent");
    const auto imposed = fl::test_sw92_pt_asymmetric_stability_against(
        3.0e6, 340.0, gate3a.feed, gate3a.common_log_activity,
        model, 0.0);
    require(imposed.status == gate3a.status &&
                imposed.aqueous.status == gate3a.aqueous.search->status &&
                imposed.nonaqueous.status == gate3a.nonaqueous.search->status &&
                imposed.aqueous.evaluations == gate3a.aqueous.search->evaluations &&
                imposed.nonaqueous.evaluations == gate3a.nonaqueous.search->evaluations &&
                imposed.negative_witnesses.size() == gate3a.negative_witnesses.size(),
            "imposed common-tangent helper changed Gate-3A search semantics");
    require(imposed.common_log_activity == gate3a.common_log_activity &&
                !imposed.global_stability_proven,
            "imposed helper changed tangent/global-proof semantics");
}

void family_tie_blocks_candidate() {
    const auto model = family_tie_model();
    const auto& reference = pair_golden[0];
    const auto result = fl::iterate_sw92_asymmetric_fixed_pair(
        3.0e6, 340.0, Vec{0.7, 0.3}, exact_log_k(reference),
        {th::SwPhaseFamily::aqueous, th::SwPhaseFamily::aqueous},
        model, 0.0);
    require(result.status ==
                fl::Sw92AsymmetricFixedPairStatus::family_assignment_nonsmooth &&
                result.point && result.equations_converged() &&
                !result.candidate_admissible(),
            "exact AQ/NA family tie did not preserve equation convergence/block candidate");
    require(result.point->phase0_assignment.status ==
                fl::Sw92AsymmetricFamilyAssignmentStatus::family_tie &&
                result.point->phase1_assignment.status ==
                    fl::Sw92AsymmetricFamilyAssignmentStatus::family_tie,
            "family tie was broken by phase slot or execution order");
}

void phase_disappearance() {
    const auto& reference = pair_golden[0];
    const auto model = binary_model();
    const double x = static_cast<double>(reference.phase0_gas);
    const double y = static_cast<double>(reference.phase1_gas);
    constexpr double beta = 5e-4;
    const double feed_gas = std::fma(beta, y - x, x);
    fl::Sw92AsymmetricFixedPairOptions options;
    options.minimum_phase_fraction = 1e-3;
    const auto result = fl::iterate_sw92_asymmetric_fixed_pair(
        3.0e6, 340.0, Vec{feed_gas, 1.0 - feed_gas},
        exact_log_k(reference),
        {th::SwPhaseFamily::aqueous, th::SwPhaseFamily::aqueous},
        model, 0.0, options);
    require(result.status == fl::Sw92AsymmetricFixedPairStatus::phase_disappearance &&
                result.point && result.equations_converged() &&
                !result.candidate_admissible(),
            "small fixed-pair phase did not preserve equations/report disappearance");
    require(result.point->chemical_potential_norm <= options.chemical_potential_tolerance &&
                result.point->mass_absolute <= options.mass_absolute_tolerance &&
                result.point->mass_relative <= options.mass_relative_tolerance &&
                result.point->phase1.mole_phase_fraction <= options.minimum_phase_fraction,
            "phase disappearance was emitted before equation/balance gates");
    require(result.point->phase0_assignment.status ==
                fl::Sw92AsymmetricFamilyAssignmentStatus::not_checked &&
                result.point->phase1_assignment.status ==
                    fl::Sw92AsymmetricFamilyAssignmentStatus::not_checked,
            "disappearing phase incorrectly entered family-dominance acceptance checks");
}

void permutation() {
    const auto result = solve_reference_pair(0, true);
    require(result.status == fl::Sw92AsymmetricFixedPairStatus::converged &&
                result.candidate_admissible(),
            "component permutation changed AQ+AQ fixed-pair status");
    check_point(result, pair_golden[0], true);
    require(result.component_ids.size() == 2 &&
                result.component_ids[0] == "water" &&
                result.component_ids[1] == "carbon-dioxide",
            "permuted component metadata lost ordered-snapshot identity");
}

void resource_and_contract_failures() {
    const auto model = binary_model();
    const auto seed = exact_log_k(pair_golden[0]);

    fl::Sw92AsymmetricFixedPairOptions limited;
    limited.max_evaluations = 2;
    const auto budget = fl::iterate_sw92_asymmetric_fixed_pair(
        3.0e6, 340.0, Vec{0.7, 0.3}, seed,
        {th::SwPhaseFamily::aqueous, th::SwPhaseFamily::aqueous},
        model, 0.0, limited);
    require(budget.status == fl::Sw92AsymmetricFixedPairStatus::evaluation_limit &&
                budget.evaluations == 2 && budget.point &&
                budget.equations_converged() && !budget.candidate_admissible(),
            "family-dominance property budget was not accounted explicitly");

    fl::Sw92AsymmetricFixedPairOptions root_limited;
    root_limited.aqueous_root_options.max_iterations = 1;
    const auto root_failure = fl::iterate_sw92_asymmetric_fixed_pair(
        3.0e6, 340.0, Vec{0.7, 0.3}, seed,
        {th::SwPhaseFamily::aqueous, th::SwPhaseFamily::aqueous},
        model, 0.0, root_limited);
    require(root_failure.status == fl::Sw92AsymmetricFixedPairStatus::property_failure &&
                root_failure.property_issue ==
                    fl::StabilityPropertyIssue::root_iteration_limit &&
                !root_failure.equations_converged(),
            "fixed-pair root-budget failure semantics changed");

    expect_error<std::invalid_argument>([&] {
        (void)fl::iterate_sw92_asymmetric_fixed_pair(
            3.0e6, 340.0, Vec{0.7, 0.3}, seed,
            {static_cast<th::SwPhaseFamily>(99), th::SwPhaseFamily::aqueous},
            model, 0.0);
    });
    expect_error<std::invalid_argument>([&] {
        (void)fl::iterate_sw92_asymmetric_fixed_pair(
            3.0e6, 340.0, Vec{0.7, 0.3}, Vec{0.0},
            {th::SwPhaseFamily::aqueous, th::SwPhaseFamily::aqueous},
            model, 0.0);
    });
}

void headers() {
    require(sw92_asymmetric_fixed_pair_header(),
            "public fixed-pair header self-containment check failed");
}

struct Case {
    const char* name;
    void (*run)();
};

constexpr Case cases[] = {
    {"aq_aq_fixed_pair", aq_aq_fixed_pair},
    {"aq_na_fixed_pair_dominated", aq_na_fixed_pair_dominated},
    {"na_na_fixed_pair_dominated", na_na_fixed_pair_dominated},
    {"imposed_helper_matches_gate3a", imposed_helper_matches_gate3a},
    {"family_tie_blocks_candidate", family_tie_blocks_candidate},
    {"phase_disappearance", phase_disappearance},
    {"permutation", permutation},
    {"resource_and_contract_failures", resource_and_contract_failures},
    {"headers", headers}
};

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: sw92_asymmetric_fixed_pair_test <case>\n";
        return 2;
    }
    for (const auto& test : cases) {
        if (std::string_view(argv[1]) != test.name) { continue; }
        try {
            test.run();
            std::cout << "[PASS] " << test.name << '\n';
            return 0;
        } catch (const std::exception& error) {
            std::cerr << "[FAIL] " << error.what() << '\n';
            return 1;
        }
    }
    std::cerr << "unknown case: " << argv[1] << '\n';
    return 2;
}
