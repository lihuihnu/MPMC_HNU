#ifndef MPMC_FLASH_PR76_PT_CONTINUATION_HPP
#define MPMC_FLASH_PR76_PT_CONTINUATION_HPP

#include <mpmc/flash/pr76_three_phase.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flash {

inline constexpr std::string_view pr76_pt_continuation_convention =
    "PR76/PT/max3/ordered-path-continuation/v1";

struct Pr76PtPathState {
    double pressure_pa{};
    double temperature_k{};
};

struct Pr76PtContinuationOptions {
    Pr76PtMax3Options point_options;
    std::vector<std::vector<double>> fallback_initial_starts;
    std::vector<std::vector<double>> fallback_final_starts;
    std::size_t max_points{4096};
    bool carry_previous_phase_compositions{true};
    bool carry_previous_three_phase_state{true};
};

enum class Pr76PtContinuationHintKind {
    none,
    one_phase,
    two_phase,
    three_phase
};

struct Pr76PtContinuationPointResult {
    Pr76PtPathState state;
    Pr76PtMax3Result solve;
    std::optional<std::size_t> accepted_phase_count;
    Pr76PtContinuationHintKind incoming_hint{Pr76PtContinuationHintKind::none};
    std::size_t carried_stability_start_count{};
    bool carried_three_phase_start{false};
    bool carried_three_phase_start_consumed{false};
};

struct Pr76PtTransitionBracket {
    std::size_t left_index{};
    std::size_t right_index{};
    std::size_t left_phase_count{};
    std::size_t right_phase_count{};
    Pr76PtPathState left;
    Pr76PtPathState right;
    bool adjacent_phase_count_step{false};
    bool exact_boundary_resolved{false};
    std::string diagnostic;
};

struct Pr76PtContinuationResult {
    static constexpr std::string_view convention = pr76_pt_continuation_convention;
    static constexpr bool global_stability_proven = false;

    Pr76PtContinuationOptions options;
    std::vector<double> feed;
    std::vector<Pr76PtContinuationPointResult> points;
    std::vector<Pr76PtTransitionBracket> transition_brackets;
    bool all_points_accepted{true};
    std::string diagnostic;
};

namespace detail {

[[nodiscard]] inline std::optional<std::size_t>
pr76_pt_continuation_phase_count(const Pr76PtMax3Result& result) noexcept {
    switch (result.status) {
    case Pr76PtMax3Status::single_phase: return 1U;
    case Pr76PtMax3Status::two_phase: return 2U;
    case Pr76PtMax3Status::three_phase: return 3U;
    case Pr76PtMax3Status::phase_boundary_unresolved:
    case Pr76PtMax3Status::higher_phase_count_or_wrong_candidate:
    case Pr76PtMax3Status::indeterminate:
        return std::nullopt;
    }
    return std::nullopt;
}

[[nodiscard]] inline Pr76PtContinuationHintKind
pr76_pt_continuation_hint_kind(std::size_t phase_count) noexcept {
    switch (phase_count) {
    case 1U: return Pr76PtContinuationHintKind::one_phase;
    case 2U: return Pr76PtContinuationHintKind::two_phase;
    case 3U: return Pr76PtContinuationHintKind::three_phase;
    default: return Pr76PtContinuationHintKind::none;
    }
}

inline void pr76_pt_continuation_append_unique(
    std::vector<std::vector<double>>& target,
    std::span<const double> composition) {
    const auto duplicate = std::find_if(
        target.begin(), target.end(),
        [&](const std::vector<double>& existing) {
            return existing.size() == composition.size() &&
                   std::equal(existing.begin(), existing.end(), composition.begin());
        });
    if (duplicate == target.end()) {
        target.emplace_back(composition.begin(), composition.end());
    }
}

[[nodiscard]] inline std::vector<std::vector<double>>
pr76_pt_continuation_phase_compositions(const Pr76PtMax3Result& result) {
    std::vector<std::vector<double>> values;
    switch (result.status) {
    case Pr76PtMax3Status::single_phase:
        values.push_back(result.base.solution.initial_stability.feed);
        break;
    case Pr76PtMax3Status::two_phase: {
        const Pr76PtSplitResult* source = result.two_phase_neighbor();
        if (source == nullptr) { source = &result.base; }
        if (const auto* candidate = source->solution.candidate()) {
            values.push_back(candidate->fractions.liquid);
            values.push_back(candidate->fractions.vapor);
        }
        break;
    }
    case Pr76PtMax3Status::three_phase:
        if (const auto* candidate = result.three_phase_candidate()) {
            for (const auto& phase : candidate->phases) {
                values.push_back(phase.composition);
            }
        }
        break;
    case Pr76PtMax3Status::phase_boundary_unresolved:
    case Pr76PtMax3Status::higher_phase_count_or_wrong_candidate:
    case Pr76PtMax3Status::indeterminate:
        break;
    }
    return values;
}

[[nodiscard]] inline std::optional<Pr76PtThreePhaseStart>
pr76_pt_continuation_three_phase_start(const Pr76PtMax3Result& result) {
    if (result.status != Pr76PtMax3Status::three_phase) { return std::nullopt; }
    const auto* candidate = result.three_phase_candidate();
    if (candidate == nullptr) { return std::nullopt; }
    Pr76PtThreePhaseStart start;
    for (std::size_t phase = 0; phase < 3U; ++phase) {
        start.compositions[phase] = candidate->phases[phase].composition;
    }
    start.phase_fraction_seed = {
        candidate->phases[1].mole_phase_fraction,
        candidate->phases[2].mole_phase_fraction};
    return start;
}

inline void pr76_pt_continuation_check_path(
    std::span<const Pr76PtPathState> path,
    std::span<const double> feed,
    const Pr76PtContinuationOptions& options) {
    if (path.empty()) {
        throw std::invalid_argument("PR76 PT continuation: ordered path is empty");
    }
    if (path.size() > options.max_points || options.max_points == 0U) {
        throw std::length_error("PR76 PT continuation: point quota exceeded");
    }
    (void)stability_check_composition(feed);
    for (const auto& state : path) {
        if (!std::isfinite(state.pressure_pa) || !(state.pressure_pa > 0.0) ||
            !std::isfinite(state.temperature_k) || !(state.temperature_k > 0.0)) {
            throw std::domain_error(
                "PR76 PT continuation: every path state requires finite p>0 Pa and T>0 K");
        }
    }
}

[[nodiscard]] inline bool pr76_pt_continuation_supplied_start_consumed(
    const Pr76PtMax3Result& result) noexcept {
    return std::any_of(
        result.attempts.begin(), result.attempts.end(),
        [](const Pr76PtThreePhaseAttempt& attempt) {
            return attempt.supplied_start.has_value() && *attempt.supplied_start == 0U;
        });
}

} // namespace detail

[[nodiscard]] inline Pr76PtContinuationResult solve_pr76_pt_continuation(
    std::span<const Pr76PtPathState> path,
    std::span<const double> feed,
    Pr76VleEvaluator& evaluator,
    Pr76PtContinuationOptions options = {}) {
    detail::pr76_pt_continuation_check_path(path, feed, options);
    if (feed.size() != evaluator.model().size()) {
        throw std::invalid_argument(
            "PR76 PT continuation: feed does not match ordered parameter snapshot");
    }

    Pr76PtContinuationResult result;
    result.options = options;
    const double feed_sum = stability_check_composition(feed);
    result.feed = stability_normalize(feed, feed_sum);
    result.points.reserve(path.size());
    result.transition_brackets.reserve(path.size() > 1U ? path.size() - 1U : 0U);

    std::optional<std::size_t> previous_phase_count;
    std::vector<std::vector<double>> previous_compositions;
    std::optional<Pr76PtThreePhaseStart> previous_three_phase;

    for (std::size_t index = 0; index < path.size(); ++index) {
        Pr76PtMax3Options point_options = options.point_options;
        std::vector<std::vector<double>> initial_starts;
        std::vector<std::vector<double>> final_starts;
        Pr76PtContinuationPointResult point;
        point.state = path[index];

        if (previous_phase_count) {
            point.incoming_hint =
                detail::pr76_pt_continuation_hint_kind(*previous_phase_count);
        }

        if (previous_phase_count && options.carry_previous_phase_compositions) {
            for (const auto& composition : previous_compositions) {
                detail::pr76_pt_continuation_append_unique(initial_starts, composition);
                detail::pr76_pt_continuation_append_unique(final_starts, composition);
            }
            point.carried_stability_start_count = previous_compositions.size();
        }
        for (const auto& composition : options.fallback_initial_starts) {
            detail::pr76_pt_continuation_append_unique(initial_starts, composition);
        }
        for (const auto& composition : options.fallback_final_starts) {
            detail::pr76_pt_continuation_append_unique(final_starts, composition);
        }

        if (previous_three_phase && options.carry_previous_three_phase_state) {
            point.carried_three_phase_start = true;
            if (point_options.three_phase_starts.size() >=
                point_options.max_three_phase_starts) {
                throw std::length_error(
                    "PR76 PT continuation: previous three-phase hint exceeds point start quota");
            }
            point_options.three_phase_starts.insert(
                point_options.three_phase_starts.begin(), *previous_three_phase);
        }

        point.solve = solve_pr76_pt_max3(
            point.state.pressure_pa, point.state.temperature_k, result.feed,
            evaluator, point_options, initial_starts, final_starts);
        point.accepted_phase_count =
            detail::pr76_pt_continuation_phase_count(point.solve);
        if (point.carried_three_phase_start) {
            point.carried_three_phase_start_consumed =
                detail::pr76_pt_continuation_supplied_start_consumed(point.solve);
        }
        if (!point.accepted_phase_count) { result.all_points_accepted = false; }

        result.points.push_back(std::move(point));
        const auto& accepted = result.points.back();

        if (index > 0U && result.points[index - 1U].accepted_phase_count &&
            accepted.accepted_phase_count &&
            *result.points[index - 1U].accepted_phase_count !=
                *accepted.accepted_phase_count) {
            Pr76PtTransitionBracket bracket;
            bracket.left_index = index - 1U;
            bracket.right_index = index;
            bracket.left_phase_count =
                *result.points[index - 1U].accepted_phase_count;
            bracket.right_phase_count = *accepted.accepted_phase_count;
            bracket.left = result.points[index - 1U].state;
            bracket.right = accepted.state;
            const std::size_t difference =
                bracket.left_phase_count > bracket.right_phase_count
                    ? bracket.left_phase_count - bracket.right_phase_count
                    : bracket.right_phase_count - bracket.left_phase_count;
            bracket.adjacent_phase_count_step = difference == 1U;
            bracket.exact_boundary_resolved = false;
            bracket.diagnostic = bracket.adjacent_phase_count_step
                ? "PR76 PT continuation: adjacent accepted path points bracket a phase-count change; the exact boundary is not resolved by the discrete scan"
                : "PR76 PT continuation: accepted path points skip one or more phase counts; path refinement is required before interpreting the transition sequence";
            result.transition_brackets.push_back(std::move(bracket));
        }

        if (accepted.accepted_phase_count) {
            previous_phase_count = accepted.accepted_phase_count;
            previous_compositions =
                detail::pr76_pt_continuation_phase_compositions(accepted.solve);
            previous_three_phase =
                detail::pr76_pt_continuation_three_phase_start(accepted.solve);
        } else {
            // Never propagate a candidate through an unresolved point. The next
            // state starts from configured fallback starts only.
            previous_phase_count.reset();
            previous_compositions.clear();
            previous_three_phase.reset();
        }
    }

    result.diagnostic = result.all_points_accepted
        ? "PR76 PT continuation: every ordered path point completed a fresh accepted max3 solve; transition brackets are discrete only"
        : "PR76 PT continuation: one or more path points remained unresolved; continuation hints were reset after each unresolved point";
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_PR76_PT_CONTINUATION_HPP
