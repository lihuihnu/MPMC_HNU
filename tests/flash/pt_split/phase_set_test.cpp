#include <mpmc/flash/pr76_phase_set.hpp>
#include <mpmc/flash/pt_vle_phase_set.hpp>

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <limits>
#include <optional>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {
namespace fl = mpmc::flash;
using Vec = std::vector<double>;

void require(bool value, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!value) {
        throw std::runtime_error(
            std::string(where.file_name()) + ":" + std::to_string(where.line()) +
            ": " + std::string(message));
    }
}

fl::PtSplitResult two_phase_source(fl::PtSplitStatus status) {
    fl::PtSplitResult source;
    source.status = status;
    source.initial_stability.pressure_pa = 4.0e6;
    source.initial_stability.temperature_k = 220.0;
    source.initial_stability.feed = {.6, .4};
    source.selected_attempt = 0;
    source.attempts.resize(1);
    auto& attempt = source.attempts[0];
    attempt.status = fl::PtSplitAttemptStatus::converged;
    attempt.point.emplace();
    auto& point = *attempt.point;
    point.fractions.vapor_fraction = .4;
    point.fractions.liquid = {3.0 / 7.0, 4.0 / 7.0};
    point.fractions.vapor = {6.0 / 7.0, 1.0 / 7.0};
    point.liquid.activity = {{std::log(2.0), std::log(.25)}, 5, true};
    point.liquid.z = .125;
    point.vapor.activity = {{0.0, 0.0}, 9, true};
    point.vapor.z = .688;
    source.diagnostic = "legacy diagnostic";
    return source;
}

void accepted_single() {
    fl::PtSplitResult source;
    source.status = fl::PtSplitStatus::single_phase_no_instability_found;
    source.initial_stability.pressure_pa = 1.5e6;
    source.initial_stability.temperature_k = 310.0;
    source.initial_stability.feed = {.2, .3, .5};
    source.initial_stability.reference = fl::StabilityPhase{{.1, .2, .3}, 7, true};

    const auto result = fl::project_pt_vle_phase_set(source);
    require(result.status == fl::PtPhaseSetStatus::accepted, "single phase not accepted");
    require(result.capability.maximum_phase_count == 2,
            "VLE capability does not report two-phase algorithm limit");
    require(!result.global_stability_proven, "projection invented global proof");
    require(result.accepted_phase_count() == 1, "accepted single phase count");
    require(result.candidate_phase_set && result.accepted_phase_set(),
            "accepted single phase payload missing");
    const auto& phase = result.accepted_phase_set()->phases.at(0);
    require(phase.mole_phase_fraction == 1.0, "single phase fraction changed");
    require(phase.composition == source.initial_stability.feed,
            "single phase composition changed");
    require(phase.activity.ln_phi == source.initial_stability.reference->ln_phi &&
                phase.activity.branch == source.initial_stability.reference->branch &&
                phase.activity.smooth == source.initial_stability.reference->smooth,
            "single phase activity changed");
    require(!phase.compressibility_factor,
            "generic stability projection fabricated a compressibility factor");
}

void accepted_two_exact() {
    const auto source = two_phase_source(fl::PtSplitStatus::two_phase_no_instability_found);
    const auto* legacy = source.candidate();
    require(legacy != nullptr, "two-phase fixture missing legacy candidate");

    const auto result = fl::project_pt_vle_phase_set(source);
    require(result.status == fl::PtPhaseSetStatus::accepted &&
                result.accepted_phase_count() == 2,
            "accepted pair status/count changed");
    require(result.pressure_pa == source.initial_stability.pressure_pa &&
                result.temperature_k == source.initial_stability.temperature_k &&
                result.feed == source.initial_stability.feed,
            "PT/feed snapshot changed");
    require(result.diagnostic == source.diagnostic, "diagnostic changed");

    const auto& phases = result.accepted_phase_set()->phases;
    require(phases.size() == 2, "accepted pair size");
    require(phases[0].mole_phase_fraction == 1.0 - legacy->fractions.vapor_fraction &&
                phases[1].mole_phase_fraction == legacy->fractions.vapor_fraction,
            "legacy phase-fraction mapping changed");
    require(phases[0].composition == legacy->fractions.liquid &&
                phases[1].composition == legacy->fractions.vapor,
            "legacy compositions changed");
    require(phases[0].activity.ln_phi == legacy->liquid.activity.ln_phi &&
                phases[0].activity.branch == legacy->liquid.activity.branch &&
                phases[0].activity.smooth == legacy->liquid.activity.smooth &&
                phases[1].activity.ln_phi == legacy->vapor.activity.ln_phi &&
                phases[1].activity.branch == legacy->vapor.activity.branch &&
                phases[1].activity.smooth == legacy->vapor.activity.smooth,
            "legacy activity data changed");
    require(phases[0].compressibility_factor && phases[1].compressibility_factor &&
                *phases[0].compressibility_factor == legacy->liquid.z &&
                *phases[1].compressibility_factor == legacy->vapor.z,
            "legacy Z values changed");
}

void unstable_candidate_retained() {
    const auto source = two_phase_source(fl::PtSplitStatus::phase_set_unstable);
    const auto result = fl::project_pt_vle_phase_set(source);
    require(result.status == fl::PtPhaseSetStatus::phase_set_unstable,
            "unstable phase set status changed");
    require(result.candidate_phase_set && result.candidate_phase_set->phases.size() == 2,
            "unstable converged pair was discarded");
    require(result.accepted_phase_set() == nullptr && result.accepted_phase_count() == 0,
            "unstable pair exposed as accepted");
}

void indeterminate_candidate_retained() {
    const auto source = two_phase_source(fl::PtSplitStatus::indeterminate);
    const auto result = fl::project_pt_vle_phase_set(source);
    require(result.status == fl::PtPhaseSetStatus::indeterminate,
            "indeterminate status changed");
    require(result.candidate_phase_set && result.candidate_phase_set->phases.size() == 2,
            "indeterminate converged pair was discarded");
    require(result.accepted_phase_set() == nullptr && result.accepted_phase_count() == 0,
            "indeterminate pair exposed as accepted");

    fl::PtSplitResult no_pair;
    no_pair.status = fl::PtSplitStatus::indeterminate;
    no_pair.initial_stability.pressure_pa = 1.0e5;
    no_pair.initial_stability.temperature_k = 300.0;
    no_pair.initial_stability.feed = {.5, .5};
    const auto empty = fl::project_pt_vle_phase_set(no_pair);
    require(!empty.candidate_phase_set && !empty.accepted_phase_set(),
            "projection fabricated candidate for unresolved flash");
}

void malformed_success_is_not_fabricated() {
    fl::PtSplitResult single;
    single.status = fl::PtSplitStatus::single_phase_no_instability_found;
    single.initial_stability.feed = {.5, .5};
    const auto missing_reference = fl::project_pt_vle_phase_set(single);
    require(missing_reference.status == fl::PtPhaseSetStatus::indeterminate &&
                !missing_reference.accepted_phase_set(),
            "missing single-phase reference became accepted");

    fl::PtSplitResult pair;
    pair.status = fl::PtSplitStatus::two_phase_no_instability_found;
    pair.initial_stability.feed = {.5, .5};
    const auto missing_pair = fl::project_pt_vle_phase_set(pair);
    require(missing_pair.status == fl::PtPhaseSetStatus::indeterminate &&
                !missing_pair.accepted_phase_set(),
            "missing two-phase candidate became accepted");

    fl::PtPhaseSetResult generic;
    generic.status = fl::PtPhaseSetStatus::accepted;
    generic.capability.maximum_phase_count = 2;
    generic.candidate_phase_set.emplace();
    require(!generic.accepted_phase_set() && generic.accepted_phase_count() == 0,
            "empty generic phase set exposed as accepted");
    generic.candidate_phase_set->phases.resize(3);
    require(!generic.accepted_phase_set() && generic.accepted_phase_count() == 0,
            "over-capability generic phase set exposed as accepted");
}

void pr76_metadata_projection() {
    fl::Pr76PtSplitResult legacy;
    legacy.solution = two_phase_source(fl::PtSplitStatus::two_phase_no_instability_found);
    legacy.dataset_id = "dataset";
    legacy.revision = "revision";
    legacy.component_ids = {"a", "b"};
    legacy.model_profile = "PR76-test";
    legacy.phase_convention = "PT-test";
    legacy.root_options.max_iterations = 17;

    const auto projected = fl::project_pr76_pt_phase_set(legacy);
    require(projected.solution.status == fl::PtPhaseSetStatus::accepted &&
                projected.solution.accepted_phase_count() == 2,
            "PR76 projection lost accepted pair");
    require(projected.dataset_id == legacy.dataset_id &&
                projected.revision == legacy.revision &&
                projected.component_ids == legacy.component_ids &&
                projected.model_profile == legacy.model_profile &&
                projected.phase_convention == legacy.phase_convention &&
                projected.root_options.max_iterations == legacy.root_options.max_iterations,
            "PR76 metadata changed during projection");
}

using Test = std::pair<std::string_view, void (*)()>;
constexpr Test tests[]{
    {"accepted_single", accepted_single},
    {"accepted_two_exact", accepted_two_exact},
    {"unstable_candidate_retained", unstable_candidate_retained},
    {"indeterminate_candidate_retained", indeterminate_candidate_retained},
    {"malformed_success_is_not_fabricated", malformed_success_is_not_fabricated},
    {"pr76_metadata_projection", pr76_metadata_projection},
};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument("one test name required");
        }
        for (const auto& [name, run] : tests) {
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
