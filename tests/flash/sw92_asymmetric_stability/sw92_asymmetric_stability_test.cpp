#include <mpmc/flash/sw92_asymmetric_stability.hpp>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <optional>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool sw92_asymmetric_stability_header();

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " + std::string(message));
    }
}

void near(double actual, long double expected, long double relative = 8e-10L,
          long double absolute = 8e-13L,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(static_cast<long double>(actual) - expected) >
            absolute + relative * std::abs(expected)) {
        std::cerr << "actual=" << actual << " expected=" << expected << '\n';
        require(false, "reference mismatch", where);
    }
}

th::Provenance paper(std::string locator) {
    return {th::SourceKind::literature,
            "doi:10.1016/0378-3812(92)85105-H",
            "Fluid Phase Equilibria 77 (1992) 217-240; authors' errata in supplied PDF",
            std::move(locator),
            "SW92 corrected-original asymmetric-stability regression; not experimental validation",
            "User PDF sha256:cb5b1d5034d78d934e887449ce0c89692d43431835d1371da5606d95b2c6bf58",
            "Bibliographic/formula facts only; source PDF is not redistributed"};
}

th::Provenance synthetic(std::string locator) {
    return {th::SourceKind::synthetic_test,
            "MPMC_HNU SW92 asymmetric stability structural fixture", "v1",
            std::move(locator),
            "Controlled AQ/NA BIP equality only for family-tie semantics",
            "tests/flash/sw92_asymmetric_stability/sw92_asymmetric_stability_test.cpp",
            "Repository test data; no physical-validation claim"};
}

th::SourcedScalar sourced(double value, th::Unit unit, const th::Provenance& source,
                          std::string original = "SI or dimensionless",
                          std::string conversion = "identity") {
    return {value, unit, source, std::move(original), std::move(conversion)};
}

struct GasSpec {
    const char* id;
    const char* display;
    th::Sw92Species species;
    double tc;
    double pc_bar;
    double omega;
    double nonaqueous_kij;
};

constexpr GasSpec co2{"carbon-dioxide", "Carbon dioxide", th::Sw92Species::carbon_dioxide,
                      304.2, 73.8, 0.2273, 0.1896};
constexpr GasSpec water_spec{"water", "Water", th::Sw92Species::water,
                             647.3, 221.2, 0.3434, 0.0};

th::Sw92Phase<double> binary_model(double nonaqueous_kij = co2.nonaqueous_kij,
                                   bool reverse = false,
                                   bool synthetic_pair = false) {
    const auto identity = paper("Table 3 component identities");
    const auto properties = paper("Table 3 physical properties");
    std::vector<th::Component> catalog{
        {co2.id, co2.display, th::ComponentKind::pure, identity, std::nullopt},
        {water_spec.id, water_spec.display, th::ComponentKind::pure, identity, std::nullopt}};

    th::Sw92ParameterInput input;
    input.model_id = std::string(th::sw92_corrected_profile);
    input.dataset_id = synthetic_pair ? "SW92-asymmetric-family-tie-structural"
                                      : "SW92-Table3-Table5-CO2-asymmetric";
    input.revision = synthetic_pair ? "synthetic-AQ-NA-equality-v1"
                                    : "user-PDF-sha256-cb5b1d50";
    input.applicability = {{std::nullopt, std::nullopt,
                            paper("Eqs. (9)-(17), Tables 2-5; no global T/p bounds inferred")},
                           std::nullopt};

    const auto add_pure = [&](const GasSpec& spec) {
        input.pure.push_back({spec.id, spec.species,
            sourced(spec.tc, th::Unit::kelvin, properties, "K", "identity"),
            sourced(spec.pc_bar * 100000.0, th::Unit::pascal, properties,
                    "bar", "bar * 100000 -> Pa"),
            sourced(spec.omega, th::Unit::dimensionless, properties)});
    };
    add_pure(co2);
    add_pure(water_spec);
    const auto pair_source = synthetic_pair
        ? synthetic("NA CO2/water BIP chosen to match AQ correlation at 340 K, m=0")
        : paper("Table 5 non-aqueous water BIP");
    input.water_binary.push_back({co2.id, th::Sw92NonAqueousWaterRule::constant,
                                  sourced(nonaqueous_kij, th::Unit::dimensionless,
                                          pair_source)});

    const std::vector<std::string> order = reverse
        ? std::vector<std::string>{water_spec.id, co2.id}
        : std::vector<std::string>{co2.id, water_spec.id};
    const auto policy = synthetic_pair ? th::DataPolicy::allow_synthetic_tests
                                       : th::DataPolicy::ordinary;
    return th::Sw92Phase<double>::from_parameters(
        th::Sw92ParameterSet::create(catalog, order, input, policy));
}

struct AsymmetricGolden {
    long double aqueous_feed_gibbs;
    long double nonaqueous_feed_gibbs;
    long double aqueous_minus_nonaqueous;
    long double aqueous_lnphi_gas;
    long double aqueous_lnphi_water;
    long double nonaqueous_lnphi_gas;
    long double nonaqueous_lnphi_water;
    long double common_d_gas;
    long double common_d_water;
    long double uniform_tpd_aqueous;
    long double uniform_tpd_nonaqueous;
};

constexpr AsymmetricGolden golden{
    -0.79037918785448424481162736483813301309149968040969608311561777113515725032080099L,
    -0.76370443257065285543341872699574589539123274164787928783219103211314806036977166L,
    -0.02667475528383138937820863784238711770026693876181679528342673902200918995102933L,
    -0.10442129405176457180570757963643373096127863629383764180363111696108247538901601L,
    -0.35473326654451860507320365298416238652569874809788141753909831114127033254389492L,
    -0.09562258679395006619140912673772785900606427248677987146847632756221492340092644L,
    -0.28634773253331448691253791693985236208697580110829356404303702299859732067933949L,
    -0.46109623799049695071834629087761820892529539534074942937756886799107522231429922L,
    -1.5587060708704545976959498707460008894793096789039049418377318784713486490026384L,
    -0.65366368836354785021166695680624633682370667556623174169610299650721362103673280L,
    -0.12598993192724328329046323832178562814059974278606138718394164301008051411730185L};

constexpr long double pure_co2_active_lnphi =
    -0.11036342178232871446599533710601795565052289138360709689409349293958245878534685L;
constexpr long double pure_water_active_lnphi =
    -4.6939356281295058506478466388927279072670265945060378825027755301771175629937994L;
constexpr double co2_pr_saturation_pressure_pa = 4148774.826264025;
constexpr double aq_co2_kij_340k_fresh = -0.058970986957265256;

void feed_reference_and_common_tangent() {
    const auto model = binary_model();
    const Vec feed{0.7, 0.3};
    const auto result = fl::test_sw92_pt_asymmetric_stability(
        3.0e6, 340.0, feed, model, 0.0);

    require(result.feed_reference_status ==
                fl::Sw92AsymmetricFeedReferenceStatus::selected &&
                result.reference_family == th::SwPhaseFamily::aqueous,
            "lower feed family was not resolved as AQ");
    require(result.equilibrium_profile == fl::sw92_xu_asymmetric_gibbs_profile &&
                result.search_convention == fl::sw92_xu_asymmetric_stability_convention &&
                result.model_profile == th::sw92_corrected_profile &&
                result.phase_convention == th::sw92_pt_convention,
            "asymmetric algorithm/model identity changed");
    require(result.aqueous.feed_reference && result.nonaqueous.feed_reference,
            "required family feed references missing");

    near(result.aqueous.feed_reduced_gibbs, golden.aqueous_feed_gibbs);
    near(result.nonaqueous.feed_reduced_gibbs, golden.nonaqueous_feed_gibbs);
    near(result.aqueous_minus_nonaqueous_feed_gibbs, golden.aqueous_minus_nonaqueous);
    require(std::abs(result.aqueous_minus_nonaqueous_feed_gibbs) >
                result.feed_gibbs_roundoff_guard,
            "resolved family Gibbs gap fell inside tie guard");

    near(result.aqueous.feed_reference->ln_phi[0], golden.aqueous_lnphi_gas);
    near(result.aqueous.feed_reference->ln_phi[1], golden.aqueous_lnphi_water);
    near(result.nonaqueous.feed_reference->ln_phi[0], golden.nonaqueous_lnphi_gas);
    near(result.nonaqueous.feed_reference->ln_phi[1], golden.nonaqueous_lnphi_water);
    require(result.common_log_activity.size() == 2, "common tangent dimension changed");
    near(result.common_log_activity[0], golden.common_d_gas);
    near(result.common_log_activity[1], golden.common_d_water);

    require(result.aqueous.search && result.nonaqueous.search,
            "both family common-reference searches were not executed");
    require(result.aqueous.search->imposed_log_activity == result.common_log_activity &&
                result.nonaqueous.search->imposed_log_activity == result.common_log_activity,
            "AQ/NA searches did not receive the exact same tangent");

    // Automatic starts are feed, uniform, and active-vertex blends. At the feed,
    // selected AQ has D=0 while NA has the positive family Gibbs gap.
    require(result.aqueous.search->trials.size() >= 2 &&
                result.nonaqueous.search->trials.size() >= 2 &&
                result.aqueous.search->trials[0].point &&
                result.nonaqueous.search->trials[0].point,
            "expected feed/uniform trials missing");
    near(result.aqueous.search->trials[0].point->value, 0.0L, 0.0L, 2e-13L);
    near(result.nonaqueous.search->trials[0].point->value,
         -golden.aqueous_minus_nonaqueous, 8e-10L, 8e-13L);
    near(result.aqueous.search->trials[1].point->value, golden.uniform_tpd_aqueous);
    near(result.nonaqueous.search->trials[1].point->value, golden.uniform_tpd_nonaqueous);
}

void negative_witness_and_counts() {
    const auto model = binary_model();
    const auto result = fl::test_sw92_pt_asymmetric_stability(
        3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0);
    require(result.status == fl::StabilityStatus::unstable,
            "robust asymmetric negative witness was not propagated");
    require(!result.global_stability_proven,
            "finite two-family search became a global proof");
    require(result.aqueous.feed_reference_evaluations == 1 &&
                result.nonaqueous.feed_reference_evaluations == 1,
            "feed-reference accounting changed");
    require(result.aqueous.total_property_evaluations().has_value() &&
                result.nonaqueous.total_property_evaluations().has_value() &&
                *result.aqueous.total_property_evaluations() > 1 &&
                *result.nonaqueous.total_property_evaluations() > 1,
            "separate family evaluation totals missing");

    bool aq_uniform = false;
    bool na_uniform = false;
    for (const auto& witness : result.negative_witnesses) {
        require(witness.point.value < -1e-10 - witness.point.roundoff_guard,
                "outer witness is not robustly negative");
        if (witness.family == th::SwPhaseFamily::aqueous && witness.trial_index == 1) {
            near(witness.point.value, golden.uniform_tpd_aqueous);
            aq_uniform = true;
        }
        if (witness.family == th::SwPhaseFamily::nonaqueous && witness.trial_index == 1) {
            near(witness.point.value, golden.uniform_tpd_nonaqueous);
            na_uniform = true;
        }
    }
    require(aq_uniform && na_uniform,
            "family-tagged prescribed negative witnesses were not retained");
}

void pure_vertex_active_gauge() {
    const auto model = binary_model();
    fl::Sw92FamilyStabilityEvaluator aq(
        model, 0.0, th::SwPhaseFamily::aqueous);
    fl::Sw92FamilyStabilityEvaluator na(
        model, 0.0, th::SwPhaseFamily::nonaqueous);

    const Vec pure_co2{1.0, 0.0};
    const auto aq_co2 = aq(3.0e6, 340.0, pure_co2);
    const auto na_co2 = na(3.0e6, 340.0, pure_co2);
    require(aq_co2.branch == na_co2.branch,
            "pure CO2 root branch differs by family");
    near(aq_co2.ln_phi[0], pure_co2_active_lnphi);
    near(na_co2.ln_phi[0], pure_co2_active_lnphi);
    require(std::abs(aq_co2.ln_phi[1] - na_co2.ln_phi[1]) > 1e-3,
            "zero-water infinite-dilution lnphi unexpectedly forced equal");

    const Vec pure_water{0.0, 1.0};
    const auto aq_water = aq(3.0e6, 340.0, pure_water);
    const auto na_water = na(3.0e6, 340.0, pure_water);
    require(aq_water.branch == na_water.branch,
            "pure-water root branch differs by family");
    near(aq_water.ln_phi[1], pure_water_active_lnphi);
    near(na_water.ln_phi[1], pure_water_active_lnphi);
    require(std::abs(aq_water.ln_phi[0] - na_water.ln_phi[0]) > 1e-3,
            "zero-CO2 infinite-dilution lnphi unexpectedly forced equal");
}

void family_tie_guard() {
    const auto exact_model = binary_model(aq_co2_kij_340k_fresh, false, true);
    const auto exact = fl::test_sw92_pt_asymmetric_stability(
        3.0e6, 340.0, Vec{0.7, 0.3}, exact_model, 0.0);
    require(exact.feed_reference_status ==
                fl::Sw92AsymmetricFeedReferenceStatus::family_gibbs_tie &&
                !exact.reference_family && exact.common_log_activity.empty() &&
                !exact.aqueous.search && !exact.nonaqueous.search,
            "exact AQ/NA mixed-state tie selected a family");
    require(std::abs(exact.aqueous_minus_nonaqueous_feed_gibbs) <=
                exact.feed_gibbs_roundoff_guard,
            "exact structural tie exceeded arithmetic guard");

    const double one_ulp = std::nextafter(
        aq_co2_kij_340k_fresh, std::numeric_limits<double>::infinity());
    const auto near_model = binary_model(one_ulp, false, true);
    const auto near_tie = fl::test_sw92_pt_asymmetric_stability(
        3.0e6, 340.0, Vec{0.7, 0.3}, near_model, 0.0);
    require(near_tie.feed_reference_status ==
                fl::Sw92AsymmetricFeedReferenceStatus::family_gibbs_tie &&
                !near_tie.reference_family && near_tie.common_log_activity.empty(),
            "one-ULP AQ/NA near-tie was broken by execution order");
}

void permutation() {
    const auto normal_model = binary_model(co2.nonaqueous_kij, false, false);
    const auto reverse_model = binary_model(co2.nonaqueous_kij, true, false);
    const auto first = fl::test_sw92_pt_asymmetric_stability(
        3.0e6, 340.0, Vec{0.7, 0.3}, normal_model, 0.0);
    const auto second = fl::test_sw92_pt_asymmetric_stability(
        3.0e6, 340.0, Vec{0.3, 0.7}, reverse_model, 0.0);
    require(first.status == second.status &&
                first.feed_reference_status == second.feed_reference_status &&
                first.reference_family == second.reference_family,
            "component permutation changed asymmetric decision");
    near(first.aqueous_minus_nonaqueous_feed_gibbs,
         second.aqueous_minus_nonaqueous_feed_gibbs, 2e-9L, 1e-12L);
    require(first.common_log_activity.size() == 2 && second.common_log_activity.size() == 2,
            "permuted common tangent missing");
    near(first.common_log_activity[0], second.common_log_activity[1], 2e-9L, 1e-12L);
    near(first.common_log_activity[1], second.common_log_activity[0], 2e-9L, 1e-12L);
    require(first.component_ids[0] == second.component_ids[1] &&
                first.component_ids[1] == second.component_ids[0],
            "component metadata did not follow ordering");
}

void family_failure_isolated() {
    const auto model = binary_model();
    fl::Sw92AsymmetricStabilityOptions options;
    options.aqueous.root_options.max_iterations = 1;
    const auto result = fl::test_sw92_pt_asymmetric_stability(
        3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0, options);
    require(result.status == fl::StabilityStatus::indeterminate &&
                result.feed_reference_status ==
                    fl::Sw92AsymmetricFeedReferenceStatus::family_reference_failure,
            "family feed-reference failure was not indeterminate");
    require(!result.aqueous.feed_reference &&
                result.aqueous.feed_reference_issue ==
                    fl::StabilityPropertyIssue::root_iteration_limit,
            "AQ root-budget failure was not preserved");
    require(result.nonaqueous.feed_reference.has_value() &&
                !result.nonaqueous.feed_reference_issue,
            "NA feed reference was not independently retained");
    require(result.aqueous.feed_reference_evaluations == 1 &&
                result.nonaqueous.feed_reference_evaluations == 1 &&
                !result.aqueous.search && !result.nonaqueous.search,
            "family-isolated feed-reference accounting changed");
}

void same_family_nonsmooth_propagates() {
    const auto model = binary_model();
    const auto result = fl::test_sw92_pt_asymmetric_stability(
        co2_pr_saturation_pressure_pa, 280.0, Vec{1.0, 0.0}, model, 0.0);
    require(result.status == fl::StabilityStatus::indeterminate &&
                result.feed_reference_status ==
                    fl::Sw92AsymmetricFeedReferenceStatus::family_reference_nonsmooth,
            "same-family Gibbs root tie did not block common tangent");
    require(result.aqueous.feed_reference && result.nonaqueous.feed_reference &&
                (!result.aqueous.feed_reference->smooth ||
                 !result.nonaqueous.feed_reference->smooth),
            "nonsmooth feed reference diagnostic was lost");
    require(!result.reference_family && result.common_log_activity.empty() &&
                !result.aqueous.search && !result.nonaqueous.search,
            "nonsmooth feed incorrectly entered asymmetric search");
}

void headers() {
    require(sw92_asymmetric_stability_header(),
            "standalone SW92 asymmetric stability public header");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test cases[] = {
    {"feed_reference_and_common_tangent", feed_reference_and_common_tangent},
    {"negative_witness_and_counts", negative_witness_and_counts},
    {"pure_vertex_active_gauge", pure_vertex_active_gauge},
    {"family_tie_guard", family_tie_guard},
    {"permutation", permutation},
    {"family_failure_isolated", family_failure_isolated},
    {"same_family_nonsmooth_propagates", same_family_nonsmooth_propagates},
    {"headers", headers}};

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
        throw std::invalid_argument("unknown test");
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
