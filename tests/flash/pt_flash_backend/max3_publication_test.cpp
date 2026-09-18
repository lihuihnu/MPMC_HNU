#include <mpmc/flash/pt_phase_set.hpp>

#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace fl = mpmc::flash;
// Structural records only: these are not thermodynamic reference states.
enum class Status { single_phase, two_phase, three_phase,
    higher_phase_count_or_wrong_candidate, phase_boundary_unresolved, indeterminate };
struct Phase {
    double mole_phase_fraction{};
    std::vector<double> composition;
    fl::StabilityPhase activity;
    double z{};
};
struct Candidate { std::vector<Phase> phases; };
struct Feed { double pressure_pa{12.0}; double temperature_k{34.0};
    std::vector<double> feed{0.2, 0.8}; };
struct Base {
    struct Solution { Feed initial_stability; } solution;
    std::string diagnostic{"base"};
};
struct Source {
    Status status{Status::single_phase};
    Base base;
    std::optional<Base> neighbor;
    std::optional<Candidate> point;
    std::string diagnostic{"source"};
    const Base* two_phase_neighbor() const { return neighbor ? &*neighbor : nullptr; }
    const Candidate* three_phase_candidate() const { return point ? &*point : nullptr; }
};
void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}
}

int main() {
    try {
        Source source;
        unsigned calls = 0U;
        const auto project_two = [&](const Base& base) {
            ++calls;
            fl::PtPhaseSetResult value;
            value.status = fl::PtPhaseSetStatus::accepted;
            value.diagnostic = base.diagnostic;
            value.candidate_phase_set.emplace();
            value.candidate_phase_set->phases.emplace_back();
            return value;
        };
        const auto project = [&] {
            return fl::detail::project_max3_solution(source, project_two, "missing");
        };
        auto result = project();
        require(calls == 1U && result.diagnostic == "base" &&
            result.capability.maximum_phase_count == 3U, "single projection");
        source.status = Status::two_phase;
        result = project();
        require(result.diagnostic == "base", "two-phase base");
        source.neighbor.emplace();
        source.neighbor->diagnostic = "fresh neighbor";
        result = project();
        require(result.diagnostic == "fresh neighbor", "fresh neighbor selection");
        const auto previous_calls = calls;
        source.status = Status::three_phase;
        result = project();
        require(result.status == fl::PtPhaseSetStatus::indeterminate &&
            result.diagnostic == "missing" && !result.candidate_phase_set &&
            calls == previous_calls, "missing accepted candidate");
        source.point.emplace();
        source.point->phases = {{0.2, {0.1, 0.9}, {}, 0.01},
                               {0.3, {0.4, 0.6}, {}, 0.02},
                               {0.5, {0.7, 0.3}, {}, 0.03}};
        result = project();
        require(result.accepted_phase_count() == 3U &&
            result.diagnostic == source.diagnostic && !result.global_stability_proven,
            "three-phase status and diagnostics");
        require(result.pressure_pa == 12.0 && result.temperature_k == 34.0 &&
            result.feed == source.base.solution.initial_stability.feed,
            "owned input provenance");
        for (std::size_t i = 0U; i < 3U; ++i) {
            const auto& actual = result.candidate_phase_set->phases[i];
            const auto& expected = source.point->phases[i];
            require(actual.mole_phase_fraction == expected.mole_phase_fraction &&
                actual.composition == expected.composition &&
                actual.compressibility_factor == expected.z, "phase payload copy");
        }
        source.point->phases[0].composition[0] = 0.99;
        source.base.solution.initial_stability.feed[0] = 0.99;
        require(result.candidate_phase_set->phases[0].composition[0] == 0.1 &&
            result.feed[0] == 0.2, "result must own its copies");
        source.status = Status::higher_phase_count_or_wrong_candidate;
        result = project();
        require(result.status == fl::PtPhaseSetStatus::phase_set_unstable &&
            result.candidate_phase_set && !result.accepted_phase_set() &&
            result.diagnostic == "source", "unstable candidate is not accepted");
        for (const auto status : {Status::phase_boundary_unresolved, Status::indeterminate}) {
            source.status = status;
            result = project();
            require(result.status == fl::PtPhaseSetStatus::indeterminate &&
                !result.candidate_phase_set && result.diagnostic == "source" &&
                result.capability.maximum_phase_count == 3U,
                "unresolved states remain unresolved");
        }
        source.status = static_cast<Status>(99);
        result = project();
        require(result.status == fl::PtPhaseSetStatus::indeterminate &&
            result.capability.maximum_phase_count == 1U && result.diagnostic.empty(),
            "invalid enum retains original default-result behavior");
        std::cout << "max3 publication structural regression passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
