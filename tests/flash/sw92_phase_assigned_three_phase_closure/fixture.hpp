#ifndef MPMC_TEST_SW92_PHASE_ASSIGNED_C2B2_FIXTURE_HPP
#define MPMC_TEST_SW92_PHASE_ASSIGNED_C2B2_FIXTURE_HPP

#include <mpmc/flash/sw92_phase_assigned_three_phase_closure.hpp>

#include "test_support.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace c2b2_test {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using Vec = std::vector<double>;

inline constexpr std::array<long double, 3> w{
    0.00047609961459396299927937053378759884939750247651603L,
    0.055737047017042219082538179008470868987768930520517L,
    0.94378685336836381791818245045774153216283356700297L};
inline constexpr std::array<long double, 3> h0{
    0.038397689071917013874374203217069095503408485127421L,
    0.96090421845165064354930537493590665262652941140037L,
    0.00069809247643234257632042184702425187006210347221172L};
inline constexpr std::array<long double, 3> h1{
    0.17553679057746919907345498904254163642751041116414L,
    0.82437939725388306583329621398674234997602155150881L,
    0.000083812168647735093248796970716013596468037327053326L};
inline constexpr std::array<long double, 3> c1_w{
    0.001L,
    0.051732255801353475032393051004342789989315710194489L,
    0.94726774419864652496760694899565721001068428980551L};
inline constexpr std::array<long double, 3> c1_h{
    0.09044499188253450439997513443561233509685899357283L,
    0.90900358658083321736790740551057339784975070113603L,
    0.00055142153663227823211746005381426705339030529113814L};

inline sw92_test::Prepared synthetic_ternary_input(bool reverse = false) {
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
        "Profile-C C2b.1/C2b.2 CH4/CO2 structural non-water pair; kij=0 in AQ/NA");
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

inline th::Sw92Phase<double> model(bool reverse = false) {
    auto prepared = synthetic_ternary_input(reverse);
    return th::Sw92Phase<double>::from_parameters(
        th::Sw92ParameterSet::create(
            prepared.catalog, prepared.order, prepared.input,
            th::DataPolicy::allow_synthetic_tests));
}

inline Vec from(const std::array<long double, 3>& source, bool reverse = false) {
    if (reverse) {
        return {static_cast<double>(source[2]), static_cast<double>(source[1]),
                static_cast<double>(source[0])};
    }
    return {static_cast<double>(source[0]), static_cast<double>(source[1]),
            static_cast<double>(source[2])};
}

inline Vec log_ratio(const Vec& numerator, const Vec& denominator) {
    Vec result(numerator.size(), 0.0);
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = std::log(numerator[i]) - std::log(denominator[i]);
    }
    return result;
}

inline Vec feed_from_c1_beta(long double beta, bool reverse = false) {
    std::array<long double, 3> feed{};
    for (std::size_t i = 0; i < 3U; ++i) {
        feed[i] = (1.0L - beta) * c1_w[i] + beta * c1_h[i];
    }
    return from(feed, reverse);
}

inline fl::Sw92PhaseAssignedJointResult solve_c1(
    long double beta = 0.73L, bool reverse = false) {
    const auto m = model(reverse);
    const Vec wv = from(c1_w, reverse);
    const Vec hv = from(c1_h, reverse);
    return fl::iterate_sw92_phase_assigned_aq_na_joint(
        3.0e6, 260.0, feed_from_c1_beta(beta, reverse),
        log_ratio(hv, wv), m, 0.0);
}

inline std::pair<fl::Sw92PhaseAssignedHSideWitnessResult, std::size_t>
solve_c2a1(const fl::Sw92PhaseAssignedJointResult& c1, bool reverse = false) {
    fl::Sw92PhaseAssignedHSideWitnessOptions options;
    options.stability.automatic_starts = false;
    const std::vector<Vec> starts{from(h1, reverse)};
    auto result = fl::test_sw92_phase_assigned_h_side_na_witness(
        c1, model(reverse), options, starts);
    if (result.status != fl::Sw92PhaseAssignedHSideWitnessStatus::
                             additional_nonaqueous_phase_witness_found) {
        throw std::runtime_error("structural C1 state lost its additional-NA witness");
    }
    for (std::size_t i = 0; i < result.negative_witnesses.size(); ++i) {
        if (result.negative_witnesses[i].usable_h_split_seed()) {
            return {std::move(result), i};
        }
    }
    throw std::runtime_error("C2a1 witness status has no usable structural seed");
}

struct Chain {
    fl::Sw92PhaseAssignedJointResult c1;
    fl::Sw92PhaseAssignedHSideWitnessResult c2a1;
    std::size_t witness_index{};
    fl::Sw92PhaseAssignedThreePhaseResult c2b1;
};

inline Chain solve_chain(long double beta = 0.73L, bool reverse = false,
                         double minimum_phase_fraction = 1e-10) {
    Chain chain;
    chain.c1 = solve_c1(beta, reverse);
    if (!chain.c1.candidate_admissible()) {
        throw std::runtime_error("structural C1 candidate did not converge");
    }
    auto witness = solve_c2a1(chain.c1, reverse);
    chain.c2a1 = std::move(witness.first);
    chain.witness_index = witness.second;
    fl::Sw92PhaseAssignedThreePhaseOptions options;
    options.minimum_phase_fraction = minimum_phase_fraction;
    chain.c2b1 = fl::solve_sw92_phase_assigned_c2b1_candidate(
        chain.c1, chain.c2a1, chain.witness_index, model(reverse), options);
    return chain;
}

inline fl::Sw92PhaseAssignedC2b2Result review_chain(
    const Chain& chain, bool reverse = false,
    fl::Sw92PhaseAssignedC2b2Options options = {}) {
    return fl::review_sw92_phase_assigned_c2b2(
        chain.c1, chain.c2a1, chain.c2b1, model(reverse), options);
}

} // namespace c2b2_test

#endif // MPMC_TEST_SW92_PHASE_ASSIGNED_C2B2_FIXTURE_HPP
