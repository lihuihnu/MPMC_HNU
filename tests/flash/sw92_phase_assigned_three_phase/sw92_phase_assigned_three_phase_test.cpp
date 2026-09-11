#include <mpmc/flash/sw92_phase_assigned_three_phase.hpp>

#include "test_support.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
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

void near(double actual, long double expected, long double relative = 3e-8L,
          long double absolute = 3e-11L,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) ||
        std::abs(static_cast<long double>(actual) - expected) >
            absolute + relative * std::abs(expected)) {
        std::cerr << "actual=" << actual << " expected=" << expected << '\n';
        require(false, "reference mismatch", where);
    }
}

struct ThreePhaseGolden {
    std::array<long double, 3> w;
    std::array<long double, 3> h0;
    std::array<long double, 3> h1;
    std::array<long double, 3> fractions; // W,H0,H1
    std::array<long double, 3> z;
    std::array<long double, 3> common;
    std::array<long double, 3> feed;
    std::array<long double, 3> c1_w;
    std::array<long double, 3> c1_h;
    long double c1_h_fraction;
    long double c1_to_h1_tpd;
};

constexpr ThreePhaseGolden golden{
    {0.00047609961459396299927937053378759884939750247651603L,
     0.055737047017042219082538179008470868987768930520517L,
     0.94378685336836381791818245045774153216283356700297L},
    {0.038397689071917013874374203217069095503408485127421L,
     0.96090421845165064354930537493590665262652941140037L,
     0.00069809247643234257632042184702425187006210347221172L},
    {0.17553679057746919907345498904254163642751041116414L,
     0.82437939725388306583329621398674234997602155150881L,
     0.000083812168647735093248796970716013596468037327053326L},
    {0.27106434534399871495181054151216492014148370477456L,
     0.45055884305467450367376159999522988799186392174002L,
     0.27837681160132678137442785849260519186665237348542L},
    {0.029530621944249737093447133004800278083840357340042L,
     0.061803898916803689911576380875073692125771584303747L,
     0.74015371526520112914484528340868399795746922002554L},
    {-1.8057440080660707935940599537279743554092445199967L,
     -0.4638573338365264608149335174262989147975127481362L,
     -9.8074116508360529571522664800832813631195917587049L},
    {0.066294844074250188211981848137997004620707065308166L,
     0.67754032727037368693731852979389113372743325358182L,
     0.25616482865537612485069962206811186165185968111002L},
    {0.001L,
     0.051732255801353475032393051004342789989315710194489L,
     0.94726774419864652496760694899565721001068428980551L},
    {0.09044499188253450439997513443561233509685899357283L,
     0.90900358658083321736790740551057339784975070113603L,
     0.00055142153663227823211746005381426705339030529113814L},
    0.73L,
    -0.094369023464624287505686816778036288492217162186957L,
};

sw92_test::Prepared synthetic_ternary_input(bool reverse = false) {
    sw92_test::Prepared prepared;
    prepared.input.model_id = std::string(th::sw92_corrected_profile);
    prepared.input.dataset_id = "SW92-C2b1-synthetic-CH4-CO2-kij0";
    prepared.input.revision = "structural-v1";
    prepared.input.applicability = {
        {std::nullopt, std::nullopt,
         sw92_test::paper("SW92 Eqs. (9)-(17), Tables 2-5; synthetic gas/gas pair is separately tagged")},
        std::nullopt};
    sw92_test::add_pure(prepared, sw92_test::methane);
    sw92_test::add_pure(prepared, sw92_test::co2);
    sw92_test::add_pure(prepared, sw92_test::water);
    prepared.input.water_binary =
        sw92_test::binary_input(sw92_test::methane).input.water_binary;
    prepared.input.water_binary.push_back(
        sw92_test::binary_input(sw92_test::co2).input.water_binary.front());
    const auto source = sw92_test::synthetic(
        "Profile-C C2b.1 CH4/CO2 structural non-water pair; kij=0 in AQ/NA");
    prepared.input.nonwater_binary.push_back({
        sw92_test::methane.id, sw92_test::co2.id,
        sw92_test::scalar(0.0, th::Unit::dimensionless, source),
        sw92_test::scalar(0.0, th::Unit::dimensionless, source)});
    prepared.order = reverse
        ? std::vector<std::string>{sw92_test::water.id,
                                   sw92_test::co2.id,
                                   sw92_test::methane.id}
        : std::vector<std::string>{sw92_test::methane.id,
                                   sw92_test::co2.id,
                                   sw92_test::water.id};
    return prepared;
}

th::Sw92Phase<double> model(bool reverse = false) {
    auto prepared = synthetic_ternary_input(reverse);
    return th::Sw92Phase<double>::from_parameters(
        th::Sw92ParameterSet::create(
            prepared.catalog, prepared.order, prepared.input,
            th::DataPolicy::allow_synthetic_tests));
}

Vec vector_from(const std::array<long double, 3>& source) {
    return {static_cast<double>(source[0]), static_cast<double>(source[1]),
            static_cast<double>(source[2])};
}

Vec permute_reverse(const std::array<long double, 3>& source) {
    return {static_cast<double>(source[2]), static_cast<double>(source[1]),
            static_cast<double>(source[0])};
}

Vec log_ratio(const Vec& numerator, const Vec& denominator) {
    Vec result(numerator.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = std::log(numerator[i]) - std::log(denominator[i]);
    }
    return result;
}

fl::Sw92PhaseAssignedJointResult solve_c1() {
    const auto m = model();
    const Vec feed = vector_from(golden.feed);
    const Vec w = vector_from(golden.c1_w);
    const Vec h = vector_from(golden.c1_h);
    return fl::iterate_sw92_phase_assigned_aq_na_joint(
        3.0e6, 260.0, feed, log_ratio(h, w), m, 0.0);
}

std::pair<fl::Sw92PhaseAssignedHSideWitnessResult, std::size_t>
solve_c2a1(const fl::Sw92PhaseAssignedJointResult& c1) {
    const auto m = model();
    fl::Sw92PhaseAssignedHSideWitnessOptions options;
    options.stability.automatic_starts = false;
    const std::vector<Vec> starts{vector_from(golden.h1)};
    auto result = fl::test_sw92_phase_assigned_h_side_na_witness(
        c1, m, options, starts);
    require(result.status == fl::Sw92PhaseAssignedHSideWitnessStatus::
                                 additional_nonaqueous_phase_witness_found,
            "synthetic C1 state lost its additional-NA witness");
    for (std::size_t i = 0; i < result.negative_witnesses.size(); ++i) {
        if (result.negative_witnesses[i].usable_h_split_seed()) {
            return {std::move(result), i};
        }
    }
    throw std::runtime_error("C2a1 reported witness-found without a usable seed");
}

void check_reference(const fl::Sw92PhaseAssignedThreePhaseState& state) {
    const auto& phases = std::array<const fl::Sw92PhaseAssignedThreePhasePhase*, 3>{
        &state.aqueous_phase, &state.hydrocarbon0_phase, &state.hydrocarbon1_phase};
    const auto expected = std::array<const std::array<long double, 3>*, 3>{
        &golden.w, &golden.h0, &golden.h1};
    for (std::size_t phase = 0; phase < 3; ++phase) {
        for (std::size_t i = 0; i < 3; ++i) {
            near(phases[phase]->composition[i], (*expected[phase])[i]);
        }
        near(phases[phase]->mole_phase_fraction, golden.fractions[phase]);
        near(phases[phase]->compressibility_factor, golden.z[phase]);
    }
    for (std::size_t i = 0; i < 3; ++i) {
        near(state.common_log_activity[i], golden.common[i], 5e-8L, 5e-11L);
    }
    require(state.chemical_potential_norm <= 1.0e-11 &&
                state.mass_absolute <= 1.0e-12 &&
                state.mass_relative <= 1.0e-10 &&
                state.generalized_rr_residual <= 2.0e-13,
            "C2b.1 reference residual gates changed");
    require(state.aqueous_phase.physical_role == fl::Sw92PhysicalPhaseRole::aqueous &&
                state.hydrocarbon0_phase.physical_role ==
                    fl::Sw92PhysicalPhaseRole::nonaqueous &&
                state.hydrocarbon1_phase.physical_role ==
                    fl::Sw92PhysicalPhaseRole::nonaqueous &&
                state.aqueous_phase.thermodynamic_family ==
                    th::SwPhaseFamily::aqueous &&
                state.hydrocarbon0_phase.thermodynamic_family ==
                    th::SwPhaseFamily::nonaqueous &&
                state.hydrocarbon1_phase.thermodynamic_family ==
                    th::SwPhaseFamily::nonaqueous,
            "C2b.1 family/role mapping changed");
}

fl::Sw92PhaseAssignedThreePhaseResult direct_reference(
    bool swap_h = false,
    fl::Sw92PhaseAssignedThreePhaseOptions options = {}) {
    const auto m = model();
    const Vec w = vector_from(golden.w);
    const Vec h0 = vector_from(swap_h ? golden.h1 : golden.h0);
    const Vec h1 = vector_from(swap_h ? golden.h0 : golden.h1);
    const std::array<double, 2> fractions{
        static_cast<double>(swap_h ? golden.fractions[2] : golden.fractions[1]),
        static_cast<double>(swap_h ? golden.fractions[1] : golden.fractions[2])};
    return fl::iterate_sw92_phase_assigned_three_phase_candidate(
        3.0e6, 260.0, vector_from(golden.feed),
        log_ratio(h0, w), log_ratio(h1, w), fractions, m, 0.0, options);
}

void synthetic_pipeline_reference() {
    const auto c1 = solve_c1();
    require(c1.candidate_admissible(), "synthetic metastable C1 candidate changed");
    for (std::size_t i = 0; i < 3; ++i) {
        near(c1.candidate()->aqueous_phase.composition[i], golden.c1_w[i]);
        near(c1.candidate()->nonaqueous_phase.composition[i], golden.c1_h[i]);
    }
    near(c1.candidate()->nonaqueous_phase.mole_phase_fraction, golden.c1_h_fraction);

    auto [c2a1, witness_index] = solve_c2a1(c1);
    require(c2a1.negative_witnesses[witness_index].point.value < -0.05,
            "synthetic additional-NA witness is no longer robustly negative");
    near(golden.c1_to_h1_tpd, golden.c1_to_h1_tpd, 0.0L, 0.0L);

    const auto result = fl::solve_sw92_phase_assigned_c2b1_candidate(
        c1, c2a1, witness_index, model());
    require(result.status == fl::Sw92PhaseAssignedThreePhaseStatus::converged_candidate &&
                result.candidate() != nullptr &&
                result.source_witness_index == witness_index &&
                !result.final_stability_checked &&
                !result.global_stability_proven &&
                !result.accepted_phase_set_published,
            "C2b.1 pipeline did not retain candidate-only semantics");
    check_reference(*result.candidate());
}

void h_slot_swap_is_representation_only() {
    const auto normal = direct_reference(false);
    const auto swapped = direct_reference(true);
    require(normal.candidate() && swapped.candidate(),
            "direct C2b.1 reference did not converge");
    check_reference(*normal.candidate());
    check_reference(*swapped.candidate());
    require(swapped.candidate()->hydrocarbon_slots_canonicalized,
            "reversed H seeds were not canonicalized");
}

void component_permutation() {
    const auto m = model(true);
    const Vec w = permute_reverse(golden.w);
    const Vec h0 = permute_reverse(golden.h0);
    const Vec h1 = permute_reverse(golden.h1);
    const auto result = fl::iterate_sw92_phase_assigned_three_phase_candidate(
        3.0e6, 260.0, permute_reverse(golden.feed),
        log_ratio(h0, w), log_ratio(h1, w),
        {static_cast<double>(golden.fractions[1]),
         static_cast<double>(golden.fractions[2])},
        m, 0.0);
    require(result.candidate(), "permuted C2b.1 reference did not converge");
    const auto& state = *result.candidate();
    near(state.aqueous_phase.composition[2], golden.w[0]);
    near(state.aqueous_phase.composition[1], golden.w[1]);
    near(state.aqueous_phase.composition[0], golden.w[2]);
    near(state.hydrocarbon0_phase.composition[2], golden.h0[0]);
    near(state.hydrocarbon0_phase.composition[1], golden.h0[1]);
    near(state.hydrocarbon0_phase.composition[0], golden.h0[2]);
    near(state.hydrocarbon1_phase.composition[2], golden.h1[0]);
    near(state.hydrocarbon1_phase.composition[1], golden.h1[1]);
    near(state.hydrocarbon1_phase.composition[0], golden.h1[2]);
}

Vec feed_from_fractions(long double fw, long double f0, long double f1) {
    Vec feed(3, 0.0);
    for (std::size_t i = 0; i < 3; ++i) {
        feed[i] = static_cast<double>(
            fw * golden.w[i] + f0 * golden.h0[i] + f1 * golden.h1[i]);
    }
    return feed;
}

fl::Sw92PhaseAssignedThreePhaseResult direct_with_fractions(
    long double fw, long double f0, long double f1, double minimum) {
    const Vec w = vector_from(golden.w);
    const Vec h0 = vector_from(golden.h0);
    const Vec h1 = vector_from(golden.h1);
    fl::Sw92PhaseAssignedThreePhaseOptions options;
    options.minimum_phase_fraction = minimum;
    return fl::iterate_sw92_phase_assigned_three_phase_candidate(
        3.0e6, 260.0, feed_from_fractions(fw, f0, f1),
        log_ratio(h0, w), log_ratio(h1, w),
        {static_cast<double>(f0), static_cast<double>(f1)},
        model(), 0.0, options);
}

void phase_disappearance_boundaries() {
    const auto h_boundary = direct_with_fractions(0.4L, 0.5L, 0.1L, 0.15);
    require(h_boundary.status ==
                fl::Sw92PhaseAssignedThreePhaseStatus::hydrocarbon_phase_disappearance &&
                h_boundary.equations_converged() && h_boundary.candidate() == nullptr,
            "H disappearance was converted to an accepted lower topology");

    const auto w_boundary = direct_with_fractions(0.1L, 0.45L, 0.45L, 0.15);
    require(w_boundary.status ==
                fl::Sw92PhaseAssignedThreePhaseStatus::aqueous_phase_disappearance &&
                w_boundary.equations_converged() && w_boundary.candidate() == nullptr,
            "W disappearance was not retained as unresolved water-topology evidence");
}

void resource_and_root_failures() {
    fl::Sw92PhaseAssignedThreePhaseOptions resource;
    resource.max_evaluations = 1U;
    const auto limited = direct_reference(false, resource);
    require(limited.status == fl::Sw92PhaseAssignedThreePhaseStatus::evaluation_limit &&
                limited.candidate() == nullptr,
            "C2b.1 property budget exhaustion was hidden");

    fl::Sw92PhaseAssignedThreePhaseOptions root;
    root.nonaqueous_root_options.max_iterations = 1;
    const auto failed = direct_reference(false, root);
    require((failed.status == fl::Sw92PhaseAssignedThreePhaseStatus::property_failure ||
             failed.status == fl::Sw92PhaseAssignedThreePhaseStatus::line_search_failed) &&
                failed.property_issue.has_value() && failed.candidate() == nullptr,
            "C2b.1 family root/property failure was not retained");
}

void source_guards() {
    const auto m = model();
    fl::Sw92PhaseAssignedJointOptions limited_options;
    limited_options.max_evaluations = 1U;
    const Vec w = vector_from(golden.c1_w);
    const Vec h = vector_from(golden.c1_h);
    const auto limited_c1 = fl::iterate_sw92_phase_assigned_aq_na_joint(
        3.0e6, 260.0, vector_from(golden.feed), log_ratio(h, w), m, 0.0,
        limited_options);
    fl::Sw92PhaseAssignedHSideWitnessResult empty_witness;
    const auto no_c1 = fl::solve_sw92_phase_assigned_c2b1_candidate(
        limited_c1, empty_witness, 0U, m);
    require(no_c1.status ==
                fl::Sw92PhaseAssignedThreePhaseStatus::source_candidate_unavailable,
            "C2b.1 ran without an admissible C1 source");

    const auto good_c1 = solve_c1();
    const auto no_witness = fl::solve_sw92_phase_assigned_c2b1_candidate(
        good_c1, empty_witness, 0U, m);
    require(no_witness.status ==
                fl::Sw92PhaseAssignedThreePhaseStatus::source_witness_unavailable,
            "C2b.1 ran without a matching C2a1 witness result");
}

bool sw92_phase_assigned_three_phase_header();

void headers() {
    require(sw92_phase_assigned_three_phase_header(),
            "C2b.1 public header probe failed");
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "expected one test-case name\n";
        return 2;
    }
    const std::string_view name{argv[1]};
    try {
        if (name == "synthetic_pipeline_reference") synthetic_pipeline_reference();
        else if (name == "slot_swap") h_slot_swap_is_representation_only();
        else if (name == "permutation") component_permutation();
        else if (name == "phase_disappearance") phase_disappearance_boundaries();
        else if (name == "resource_root_failures") resource_and_root_failures();
        else if (name == "source_guards") source_guards();
        else if (name == "headers") headers();
        else throw std::invalid_argument("unknown test-case name");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
