#include <mpmc/flash/sw92_asymmetric_fixed_pair.hpp>
#include <mpmc/flash/sw92_phase_assigned_joint.hpp>

#include "test_support.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool sw92_phase_assigned_joint_header();

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

struct JointGolden {
    const sw92_test::Spec* gas;
    double pressure_pa;
    double temperature_k;
    double molality;
    double feed_gas;
    long double aqueous_gas;
    long double nonaqueous_gas;
    long double nonaqueous_fraction;
    long double aqueous_z;
    long double nonaqueous_z;
    long double reduced_gibbs;
    long double common_gas;
    long double common_water;
    long double aqueous_minus_nonaqueous_water;
};

constexpr JointGolden golden[] = {
    {&sw92_test::co2, 3.0e6, 340.0, 0.0, 0.7,
     0.005945065166186101440737568231518653586387892719453330789L,
     0.9887877343090553106585559885449699850347182199737051566L,
     0.7061709433505718241361429803823283207927897009254034181L,
     0.02326916182253476550316956093117009661761092895870925982L,
     0.8865753259957614197158614840436129757850208327439343804L,
     -1.495046828318048125530732605527113449469374366608903618L,
     -0.1216205040339955724088922848297470366629002170572176330L,
     -4.699708251647504082815026687154301746017814048896170917L,
     0.9828426691428692092178184203134513314483303272542518259L},
    {&sw92_test::methane, 1.0e7, 350.0, 1.0, 0.5,
     0.0009107995472710484575836383440685787307489823615427895675L,
     0.9946785444670317547058102612604360380162061098691410614L,
     0.5022191583537727541331431762679649316598376913439483550L,
     0.07556724617135668510729623910614985394789699839543768716L,
     0.9106256398210770581205949691770014617913165687529905173L,
     -2.782608641638294885206413508771089052433567550674659001L,
     -0.1136325420988743974800963236951589339806957731490181828L,
     -5.451584741177715372932730693847019170886439328200299818L,
     0.9937677449197607062482266229163674592854571275075982718L},
};

th::Sw92Phase<double> model_for(const JointGolden& reference, bool reverse = false) {
    return th::Sw92Phase<double>::from_parameters(
        sw92_test::binary_parameters(*reference.gas, reverse));
}

Vec exact_log_k(const JointGolden& reference, bool reverse = false) {
    const double x = static_cast<double>(reference.aqueous_gas);
    const double y = static_cast<double>(reference.nonaqueous_gas);
    Vec result{std::log(y / x), std::log((1.0 - y) / (1.0 - x))};
    if (reverse) { std::swap(result[0], result[1]); }
    return result;
}

Vec rough_log_k(std::size_t index) {
    return index == 0U ? Vec{5.0, -4.4} : Vec{6.5, -5.0};
}

fl::Sw92PhaseAssignedJointResult solve_reference(
    std::size_t index, bool reverse = false, bool rough = true,
    fl::Sw92PhaseAssignedJointOptions options = {}) {
    const auto& reference = golden[index];
    const auto model = model_for(reference, reverse);
    Vec feed = reverse
        ? Vec{1.0 - reference.feed_gas, reference.feed_gas}
        : Vec{reference.feed_gas, 1.0 - reference.feed_gas};
    Vec seed = rough ? rough_log_k(index) : exact_log_k(reference, reverse);
    if (rough && reverse) { std::swap(seed[0], seed[1]); }
    return fl::iterate_sw92_phase_assigned_aq_na_joint(
        reference.pressure_pa, reference.temperature_k, feed, seed,
        model, reference.molality, options);
}

void check_reference(const fl::Sw92PhaseAssignedJointResult& result,
                     const JointGolden& reference, bool reverse = false) {
    require(result.status == fl::Sw92PhaseAssignedJointStatus::converged_candidate,
            "Profile-C C1 did not publish a converged two-phase candidate");
    require(result.candidate_admissible() && result.candidate() != nullptr,
            "Profile-C candidate publication guard failed");
    require(result.equations_converged(), "Profile-C joint equations did not converge");
    require(!result.final_stability_checked && !result.global_stability_proven,
            "C1 must not claim final/global stability");
    require(result.equilibrium_profile == fl::sw92_phase_assigned_aq_na_joint_profile &&
                result.primitive_convention ==
                    fl::sw92_phase_assigned_aq_na_joint_primitive &&
                result.model_profile == th::sw92_corrected_profile,
            "Profile-C algorithm/model identity changed");

    const auto& point = *result.point;
    const std::size_t gas_index = reverse ? 1U : 0U;
    const std::size_t water_index = reverse ? 0U : 1U;
    require(point.aqueous_phase.physical_role == fl::Sw92PhysicalPhaseRole::aqueous &&
                point.aqueous_phase.thermodynamic_family == th::SwPhaseFamily::aqueous &&
                point.nonaqueous_phase.physical_role == fl::Sw92PhysicalPhaseRole::nonaqueous &&
                point.nonaqueous_phase.thermodynamic_family == th::SwPhaseFamily::nonaqueous,
            "physical role and thermodynamic family mapping changed");
    near(point.aqueous_phase.composition[gas_index], reference.aqueous_gas);
    near(point.nonaqueous_phase.composition[gas_index], reference.nonaqueous_gas);
    near(point.aqueous_phase.composition[water_index], 1.0L - reference.aqueous_gas);
    near(point.nonaqueous_phase.composition[water_index],
         1.0L - reference.nonaqueous_gas);
    near(point.nonaqueous_phase.mole_phase_fraction,
         reference.nonaqueous_fraction);
    near(point.aqueous_phase.compressibility_factor, reference.aqueous_z);
    near(point.nonaqueous_phase.compressibility_factor, reference.nonaqueous_z);
    near(point.reduced_gibbs, reference.reduced_gibbs);
    near(point.common_log_activity[gas_index], reference.common_gas);
    near(point.common_log_activity[water_index], reference.common_water);
    near(point.aqueous_minus_nonaqueous_water_fraction,
         reference.aqueous_minus_nonaqueous_water);
    require(point.aqueous_minus_nonaqueous_water_fraction >
                point.water_role_roundoff_guard,
            "physical aqueous phase is not robustly water-richer");
}

void co2_reference() {
    const auto result = solve_reference(0U);
    check_reference(result, golden[0]);
    require(result.iterations > 0, "CO2 rough seed did not exercise C1 iteration");
}

void methane_reference() {
    const auto result = solve_reference(1U);
    check_reference(result, golden[1]);
    require(result.iterations > 0, "CH4 rough seed did not exercise C1 iteration");
}

void profile_b_dominance_is_not_profile_c_acceptance() {
    const auto c = solve_reference(0U, false, false);
    check_reference(c, golden[0]);

    const auto model = model_for(golden[0]);
    const Vec feed{0.7, 0.3};
    const Vec seed = exact_log_k(golden[0]);
    const auto b = fl::iterate_sw92_asymmetric_fixed_pair(
        3.0e6, 340.0, feed, seed,
        {th::SwPhaseFamily::aqueous, th::SwPhaseFamily::nonaqueous},
        model, 0.0);
    require(b.equations_converged(),
            "Profile-B comparison lost the same AQ+NA equation solution");
    require(b.status == fl::Sw92AsymmetricFixedPairStatus::family_assignment_dominated,
            "Profile-B lower-envelope comparison no longer distinguishes Profile C");
    require(c.status == fl::Sw92PhaseAssignedJointStatus::converged_candidate,
            "Profile-C incorrectly inherited Profile-B family dominance");
}

void component_permutation() {
    const auto normal = solve_reference(0U, false, false);
    const auto reverse = solve_reference(0U, true, false);
    check_reference(normal, golden[0], false);
    check_reference(reverse, golden[0], true);
    require(normal.component_ids.size() == 2U && reverse.component_ids.size() == 2U &&
                normal.component_ids[0] == reverse.component_ids[1] &&
                normal.component_ids[1] == reverse.component_ids[0],
            "ordered component snapshot did not permute as requested");
}

void reversed_physical_roles_are_rejected() {
    const auto model = model_for(golden[0]);
    const Vec feed{0.7, 0.3};
    const Vec reversed_seed{
        -8.1546941242088137924135629730483208980,
         4.3692706080384401674487480614282952421};
    const auto result = fl::iterate_sw92_phase_assigned_aq_na_joint(
        3.0e6, 340.0, feed, reversed_seed, model, 0.0);
    require(result.equations_converged(),
            "reversed-role reference should remain an equation-converged diagnostic");
    require(result.status == fl::Sw92PhaseAssignedJointStatus::phase_role_reversed,
            "water-poor AQ / water-rich NA solution was not rejected by physical-role topology");
    require(result.candidate() == nullptr,
            "reversed physical-role state leaked through C1 candidate guard");
}

void phase_disappearance() {
    const auto& reference = golden[0];
    const auto model = model_for(reference);
    fl::Sw92PhaseAssignedJointOptions options;
    options.minimum_phase_fraction = 1e-4;
    const double beta = 0.75 * options.minimum_phase_fraction;
    const double x = static_cast<double>(reference.aqueous_gas);
    const double y = static_cast<double>(reference.nonaqueous_gas);
    const double z = std::fma(beta, y - x, x);
    const Vec feed{z, 1.0 - z};
    const Vec seed = exact_log_k(reference);
    const auto result = fl::iterate_sw92_phase_assigned_aq_na_joint(
        reference.pressure_pa, reference.temperature_k, feed, seed,
        model, reference.molality, options);
    require(result.status == fl::Sw92PhaseAssignedJointStatus::phase_disappearance,
            "C1 did not classify a converged vanishing phase conservatively");
    require(result.equations_converged(),
            "phase disappearance must be classified only after equation/balance convergence");
    require(result.point && result.point->nonaqueous_phase.mole_phase_fraction > 0.0 &&
                result.point->nonaqueous_phase.mole_phase_fraction <=
                    options.minimum_phase_fraction,
            "phase-disappearance fraction is outside the configured boundary");
    require(result.candidate() == nullptr,
            "vanishing two-phase state leaked through candidate publication guard");
}

void resource_and_contract_failures() {
    const auto model = model_for(golden[0]);
    const Vec feed{0.7, 0.3};
    const Vec seed = exact_log_k(golden[0]);

    fl::Sw92PhaseAssignedJointOptions limited;
    limited.max_evaluations = 1;
    const auto exhausted = fl::iterate_sw92_phase_assigned_aq_na_joint(
        3.0e6, 340.0, feed, seed, model, 0.0, limited);
    require(exhausted.status == fl::Sw92PhaseAssignedJointStatus::evaluation_limit,
            "property budget exhaustion did not remain explicit");
    require(exhausted.candidate() == nullptr,
            "resource-limited result leaked a candidate");

    expect_error<std::invalid_argument>([&] {
        const Vec bad_seed{1.0};
        (void)fl::iterate_sw92_phase_assigned_aq_na_joint(
            3.0e6, 340.0, feed, bad_seed, model, 0.0);
    });
    expect_error<std::domain_error>([&] {
        const Vec bad_seed{std::numeric_limits<double>::quiet_NaN(), 0.0};
        (void)fl::iterate_sw92_phase_assigned_aq_na_joint(
            3.0e6, 340.0, feed, bad_seed, model, 0.0);
    });
}

void headers() {
    require(sw92_phase_assigned_joint_header(),
            "phase-assigned joint public header probe failed");
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "expected one test-case name\n";
        return 2;
    }
    const std::string_view name{argv[1]};
    try {
        if (name == "co2_reference") co2_reference();
        else if (name == "methane_reference") methane_reference();
        else if (name == "profile_b_separation") profile_b_dominance_is_not_profile_c_acceptance();
        else if (name == "permutation") component_permutation();
        else if (name == "reversed_role") reversed_physical_roles_are_rejected();
        else if (name == "phase_disappearance") phase_disappearance();
        else if (name == "resource_and_contract") resource_and_contract_failures();
        else if (name == "headers") headers();
        else throw std::invalid_argument("unknown test-case name");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
