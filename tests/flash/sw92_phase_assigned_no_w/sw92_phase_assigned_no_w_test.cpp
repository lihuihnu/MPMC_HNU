#include <mpmc/flash/sw92_phase_assigned_no_w.hpp>

#include "test_support.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

bool sw92_phase_assigned_no_w_header();

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

void near(double actual, double expected, double relative = 2e-8,
          double absolute = 2e-11,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) > absolute + relative * std::abs(expected)) {
        std::cerr << "actual=" << actual << " expected=" << expected << '\n';
        require(false, "numeric mismatch", where);
    }
}

th::Sw92Phase<double> model(bool reverse = false) {
    return th::Sw92Phase<double>::from_parameters(
        sw92_test::binary_parameters(sw92_test::co2, reverse));
}

void wet_feed_finds_water_phase() {
    const auto result = fl::solve_sw92_phase_assigned_no_w(
        3.0e6, 340.0, Vec{0.7, 0.3}, model(), 0.0);
    require(result.status ==
                fl::Sw92PhaseAssignedNoWStatus::aqueous_phase_witness_found &&
                result.selected_water_witness() != nullptr &&
                result.aqueous_appearance_search.has_value(),
            "wet fixed-NA state did not produce a targeted W witness");
    require(result.retained_na_candidate_count() == 2U,
            "wet fixed-NA source did not retain its accepted NA candidate pair");
    const auto* witness = result.selected_water_witness();
    require(witness->usable_w_seed() &&
                witness->water_role_candidate_admissible &&
                witness->distinct_from_every_retained_na_candidate &&
                witness->point.composition[result.water_index] >
                    result.min_retained_na_candidate_water_fraction,
            "selected W witness failed candidate-generation role/distinction guards");
    require(result.max_retained_na_candidate_water_fraction >
                witness->point.composition[result.water_index],
            "fixture lost its W-like fixed-NA mathematical candidate");
    require(!result.global_stability_proven &&
                !result.accepted_phase_set_published && !result.morphology_resolved,
            "no-W adapter exceeded finite/candidate-only semantics");
}

void dry_feed_closes_single_h() {
    const auto result = fl::solve_sw92_phase_assigned_no_w(
        3.0e6, 340.0, Vec{0.995, 0.005}, model(), 0.0);
    require(result.status ==
                fl::Sw92PhaseAssignedNoWStatus::no_w_single_h_locally_closed &&
                result.no_w_locally_closed() &&
                result.retained_na_candidate_count() == 1U &&
                result.aqueous_appearance_search.has_value() &&
                result.selected_water_witness() == nullptr,
            "dry feed did not close as one unresolved NA/H candidate under targeted AQ review");
    require(result.hydrocarbon_flash.solution.status ==
                fl::PtSplitStatus::single_phase_no_instability_found,
            "dry no-W source unexpectedly split inside NA");
}

void zero_water_closes_by_inventory() {
    const auto result = fl::solve_sw92_phase_assigned_no_w(
        3.0e6, 340.0, Vec{1.0, 0.0}, model(), 0.0);
    require(result.no_w_locally_closed() &&
                result.retained_na_candidate_count() == 1U &&
                !result.aqueous_appearance_search &&
                !result.targeted_water_start,
            "zero-water inventory did not bypass impossible AQ appearance");
}

void witness_guards() {
    const Vec feed{0.7, 0.3};
    // Candidate water fractions are 0.9 and 0.1. A future joint solve decides
    // which candidates can survive as physical H phases.
    const std::vector<Vec> retained{{0.1, 0.9}, {0.9, 0.1}};
    constexpr double min_na_water = 0.1;

    fl::TpdPoint duplicate;
    duplicate.composition = retained.front();
    duplicate.value = -1.0;
    const auto duplicate_result = fl::detail::sw92_phase_assigned_classify_water_witness(
        0U, duplicate, retained, feed, 1U, min_na_water, 1e-7);
    require(!duplicate_result.distinct_from_every_retained_na_candidate &&
                !duplicate_result.usable_w_seed(),
            "retained NA candidate was misclassified as a W seed");

    fl::TpdPoint water_poor;
    water_poor.composition = {0.95, 0.05};
    water_poor.value = -1.0;
    const auto rejected = fl::detail::sw92_phase_assigned_classify_water_witness(
        1U, water_poor, retained, feed, 1U, min_na_water, 1e-7);
    require(rejected.distinct_from_every_retained_na_candidate &&
                !rejected.water_role_candidate_admissible &&
                !rejected.usable_w_seed(),
            "water-poor AQ mathematical negative passed the W/H contrast guard");

    fl::TpdPoint water_rich;
    water_rich.composition = {0.02, 0.98};
    water_rich.value = -1.0;
    const auto usable = fl::detail::sw92_phase_assigned_classify_water_witness(
        2U, water_rich, retained, feed, 1U, min_na_water, 1e-7);
    require(usable.distinct_from_every_retained_na_candidate &&
                usable.water_role_candidate_admissible && usable.usable_w_seed(),
            "distinct water-richer AQ trial did not become a usable W seed");
}

void component_permutation() {
    const auto normal = fl::solve_sw92_phase_assigned_no_w(
        3.0e6, 340.0, Vec{0.7, 0.3}, model(false), 0.0);
    const auto reverse = fl::solve_sw92_phase_assigned_no_w(
        3.0e6, 340.0, Vec{0.3, 0.7}, model(true), 0.0);
    require(normal.status == reverse.status &&
                normal.status ==
                    fl::Sw92PhaseAssignedNoWStatus::aqueous_phase_witness_found,
            "component permutation changed targeted W-appearance status");
    require(normal.water_index == 1U && reverse.water_index == 0U,
            "water index did not follow component permutation");
    require(normal.retained_na_candidate_count() ==
                reverse.retained_na_candidate_count(),
            "component permutation changed retained NA candidate count");
    near(normal.min_retained_na_candidate_water_fraction,
         reverse.min_retained_na_candidate_water_fraction);
    near(normal.max_retained_na_candidate_water_fraction,
         reverse.max_retained_na_candidate_water_fraction);
    near(normal.selected_water_witness()->point.value,
         reverse.selected_water_witness()->point.value, 2e-6, 2e-11);
}

void failures_remain_explicit() {
    fl::Sw92PhaseAssignedNoWOptions root_limited;
    root_limited.nonaqueous_root_options.max_iterations = 1;
    const auto root = fl::solve_sw92_phase_assigned_no_w(
        3.0e6, 340.0, Vec{0.7, 0.3}, model(), 0.0, root_limited);
    require(root.status == fl::Sw92PhaseAssignedNoWStatus::indeterminate &&
                root.hydrocarbon_flash.solution.initial_stability.reference_issue ==
                    fl::StabilityPropertyIssue::root_iteration_limit,
            "fixed-NA root failure was hidden by no-W orchestration");

    fl::Sw92PhaseAssignedNoWOptions aq_limited;
    aq_limited.aqueous_appearance.max_evaluations = 1U;
    const auto aq = fl::solve_sw92_phase_assigned_no_w(
        3.0e6, 340.0, Vec{0.995, 0.005}, model(), 0.0, aq_limited);
    require(aq.status == fl::Sw92PhaseAssignedNoWStatus::indeterminate &&
                aq.aqueous_appearance_search.has_value(),
            "targeted AQ resource exhaustion was converted to no-W success");
}

void headers() {
    require(sw92_phase_assigned_no_w_header(),
            "no-W public header probe failed");
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "expected one test-case name\n";
        return 2;
    }
    const std::string_view name{argv[1]};
    try {
        if (name == "wet_water_witness") wet_feed_finds_water_phase();
        else if (name == "dry_single_h") dry_feed_closes_single_h();
        else if (name == "zero_water") zero_water_closes_by_inventory();
        else if (name == "witness_guards") witness_guards();
        else if (name == "permutation") component_permutation();
        else if (name == "failures") failures_remain_explicit();
        else if (name == "headers") headers();
        else throw std::invalid_argument("unknown test-case name");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
