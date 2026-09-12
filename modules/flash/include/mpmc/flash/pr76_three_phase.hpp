#ifndef MPMC_FLASH_PR76_THREE_PHASE_HPP
#define MPMC_FLASH_PR76_THREE_PHASE_HPP

#include <mpmc/flash/pr76_split.hpp>
#include <mpmc/flash/pt_three_phase.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::flash {

inline constexpr const char* pr76_pt_max3_convention =
    "PR76/PT/max3/TPD-VLE-additional-phase-generalized-RR/v1";

enum class Pr76RootSide { lower_admissible, upper_admissible };

class Pr76ThreePhaseEvaluator {
public:
    Pr76ThreePhaseEvaluator(
        const thermodynamics::Pr76Phase<double>& model,
        std::array<Pr76RootSide, 3> sides,
        thermodynamics::Pr76RootOptions root_options = {})
        : model_(model), sides_(sides), root_options_(root_options) {
        if (model_.size() == 0U || root_options_.max_iterations <= 0) {
            throw std::invalid_argument(
                "PR76 three-phase evaluator: invalid model or root options");
        }
    }

    [[nodiscard]] PtThreePhaseProperty operator()(
        double pressure_pa, double temperature_k,
        std::span<const double> composition, std::size_t slot) {
        namespace th = thermodynamics;
        if (slot >= sides_.size()) {
            throw std::invalid_argument("PR76 three-phase evaluator: invalid phase slot");
        }
        try {
            const auto roots = model_.roots_full(
                pressure_pa, temperature_k, composition,
                workspace_, root_options_);
            switch (roots.status) {
            case th::Pr76RootStatus::near_multiple:
                throw StabilityPropertyError(
                    StabilityPropertyIssue::root_topology,
                    "PR76 three-phase: unresolved root topology");
            case th::Pr76RootStatus::iteration_limit:
                throw StabilityPropertyError(
                    StabilityPropertyIssue::root_iteration_limit,
                    "PR76 three-phase: root iteration limit");
            case th::Pr76RootStatus::unrepresentable:
                throw StabilityPropertyError(
                    StabilityPropertyIssue::root_range,
                    "PR76 three-phase: unrepresentable root set");
            case th::Pr76RootStatus::success: break;
            }
            std::vector<std::size_t> admissible;
            for (std::size_t k = 0; k < roots.count; ++k) {
                if (roots.roots[k].slope_sign > 0) { admissible.push_back(k); }
            }
            if (admissible.empty()) {
                throw StabilityPropertyError(
                    StabilityPropertyIssue::no_admissible_branch,
                    "PR76 three-phase: no mechanically admissible root");
            }
            const std::size_t selected =
                sides_[slot] == Pr76RootSide::lower_admissible
                    ? admissible.front() : admissible.back();
            if (!roots.roots[selected].derivative_valid) {
                throw StabilityPropertyError(
                    StabilityPropertyIssue::ill_conditioned_root,
                    "PR76 three-phase: selected root is ill-conditioned");
            }
            auto values = model_.evaluate_full(
                pressure_pa, temperature_k, composition,
                selected, workspace_, root_options_);
            return {{std::move(values.ln_phi), selected, true}, values.z};
        } catch (const th::Pr76PhaseError& error) {
            StabilityPropertyIssue issue = StabilityPropertyIssue::root_range;
            switch (error.code()) {
            case th::Pr76PhaseErrorCode::near_multiple:
                issue = StabilityPropertyIssue::root_topology; break;
            case th::Pr76PhaseErrorCode::iteration_limit:
                issue = StabilityPropertyIssue::root_iteration_limit; break;
            case th::Pr76PhaseErrorCode::ill_conditioned_derivative:
                issue = StabilityPropertyIssue::ill_conditioned_root; break;
            case th::Pr76PhaseErrorCode::unrepresentable_root: break;
            }
            throw StabilityPropertyError(issue, error.what());
        } catch (const std::range_error& error) {
            throw StabilityPropertyError(
                StabilityPropertyIssue::root_range, error.what());
        }
    }

private:
    const thermodynamics::Pr76Phase<double>& model_;
    std::array<Pr76RootSide, 3> sides_;
    thermodynamics::Pr76RootOptions root_options_;
    thermodynamics::Pr76PhaseWorkspace<double> workspace_;
};

// Optional initialization hint. These compositions can only be consumed after
// the established two-phase final common-tangent review has already proved that
// the current two-phase set is unstable. A start is never phase-count evidence.
struct Pr76PtThreePhaseStart {
    std::array<std::vector<double>, 3> compositions;
    std::array<double, 2> phase_fraction_seed{1.0 / 3.0, 1.0 / 3.0};
};

struct Pr76PtMax3Options {
    PtSplitOptions two_phase;
    PtThreePhaseOptions three_phase;
    StabilityOptions final_three_phase_stability;
    std::size_t max_three_phase_attempts{16};
    double new_phase_seed_fraction{0.1};
    std::vector<Pr76PtThreePhaseStart> three_phase_starts;
    std::size_t max_three_phase_starts{16};
    std::size_t max_three_phase_start_entries{12288};
};

enum class Pr76PtMax3Status {
    single_phase,
    two_phase,
    three_phase,
    phase_boundary_unresolved,
    higher_phase_count_or_wrong_candidate,
    indeterminate
};

struct Pr76PtThreePhaseAttempt {
    std::optional<std::size_t> witness_trial;
    std::optional<std::size_t> supplied_start;
    std::array<Pr76RootSide, 3> root_sides{
        Pr76RootSide::lower_admissible,
        Pr76RootSide::upper_admissible,
        Pr76RootSide::lower_admissible};
    PtThreePhaseResult equilibrium;
    std::optional<StabilityResult> final_stability;
    double common_reference_allowance{};
    std::optional<Pr76PtSplitResult> boundary_neighbor;
    bool accepted_three_phase{false};
    bool accepted_two_phase_neighbor{false};
};

struct Pr76PtMax3Result {
    static constexpr bool global_stability_proven = false;
    static constexpr const char* convention = pr76_pt_max3_convention;

    Pr76PtMax3Status status{Pr76PtMax3Status::indeterminate};
    Pr76PtMax3Options options;
    Pr76PtSplitResult base;
    std::vector<Pr76PtThreePhaseAttempt> attempts;
    std::optional<std::size_t> selected_attempt;
    bool attempt_limit_reached{false};
    std::string diagnostic;

    [[nodiscard]] const PtThreePhaseState* three_phase_candidate() const & noexcept {
        if (status != Pr76PtMax3Status::three_phase || !selected_attempt ||
            *selected_attempt >= attempts.size()) {
            return nullptr;
        }
        return attempts[*selected_attempt].equilibrium.candidate();
    }
    const PtThreePhaseState* three_phase_candidate() const && = delete;

    [[nodiscard]] const Pr76PtSplitResult* two_phase_neighbor() const & noexcept {
        if (status != Pr76PtMax3Status::two_phase || !selected_attempt ||
            *selected_attempt >= attempts.size()) {
            return nullptr;
        }
        const auto& neighbor = attempts[*selected_attempt].boundary_neighbor;
        if (!neighbor || !attempts[*selected_attempt].accepted_two_phase_neighbor) {
            return nullptr;
        }
        return &*neighbor;
    }
    const Pr76PtSplitResult* two_phase_neighbor() const && = delete;
};

namespace detail {

inline void pr76_pt_max3_check_options(const Pr76PtMax3Options& options) {
    pt_three_phase_check_options(options.three_phase);
    stability_check_options(options.final_three_phase_stability);
    if (options.max_three_phase_attempts == 0U ||
        options.max_three_phase_starts == 0U ||
        options.max_three_phase_start_entries == 0U ||
        options.three_phase_starts.size() > options.max_three_phase_starts ||
        !std::isfinite(options.new_phase_seed_fraction) ||
        !(options.new_phase_seed_fraction > 0.0) ||
        !(options.new_phase_seed_fraction < 1.0)) {
        throw std::invalid_argument("PR76 max3: invalid attempt/start/seed option");
    }
}

inline void pr76_pt_max3_check_starts(
    std::span<const double> feed,
    const Pr76PtMax3Options& options) {
    const std::size_t n = feed.size();
    if (n == 0U) {
        throw std::invalid_argument("PR76 max3: empty feed");
    }
    (void)stability_check_composition(feed);
    if (options.three_phase_starts.size() >
        std::numeric_limits<std::size_t>::max() / 3U ||
        options.three_phase_starts.size() * 3U >
        std::numeric_limits<std::size_t>::max() / n ||
        options.three_phase_starts.size() * 3U * n >
            options.max_three_phase_start_entries) {
        throw std::length_error("PR76 max3: three-phase start entry quota exceeded");
    }
    for (const auto& start : options.three_phase_starts) {
        if (!rr3_fractions_feasible(start.phase_fraction_seed)) {
            throw std::domain_error(
                "PR76 max3: three-phase fraction seed lies outside simplex");
        }
        for (const auto& composition : start.compositions) {
            if (composition.size() != n) {
                throw std::invalid_argument(
                    "PR76 max3: three-phase start dimension mismatch");
            }
            (void)stability_check_composition(composition);
            for (std::size_t i = 0; i < n; ++i) {
                if ((feed[i] > 0.0) != (composition[i] > 0.0)) {
                    throw std::domain_error(
                        "PR76 max3: three-phase start support mismatch");
                }
            }
        }
    }
}

[[nodiscard]] inline std::optional<std::vector<double>> pr76_three_phase_log_ratio(
    std::span<const double> numerator,
    std::span<const double> denominator,
    std::span<const double> feed) {
    if (numerator.size() != denominator.size() ||
        numerator.size() != feed.size()) {
        return std::nullopt;
    }
    std::vector<double> value(feed.size(), 0.0);
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0.0) {
            if (numerator[i] != 0.0 || denominator[i] != 0.0) {
                return std::nullopt;
            }
            continue;
        }
        if (!(numerator[i] > 0.0) || !(denominator[i] > 0.0)) {
            return std::nullopt;
        }
        value[i] = std::log(numerator[i]) - std::log(denominator[i]);
        if (!std::isfinite(value[i])) { return std::nullopt; }
    }
    return value;
}

[[nodiscard]] inline std::vector<Pr76RootSide> pr76_witness_root_sides(
    double pressure_pa, double temperature_k,
    std::span<const double> composition,
    std::size_t source_branch,
    const thermodynamics::Pr76Phase<double>& model,
    thermodynamics::Pr76RootOptions root_options) {
    thermodynamics::Pr76PhaseWorkspace<double> workspace;
    const auto roots = model.roots_full(
        pressure_pa, temperature_k, composition, workspace, root_options);
    if (roots.status != thermodynamics::Pr76RootStatus::success) {
        return {};
    }
    std::vector<std::size_t> admissible;
    for (std::size_t k = 0; k < roots.count; ++k) {
        if (roots.roots[k].slope_sign > 0) { admissible.push_back(k); }
    }
    if (admissible.empty() ||
        std::find(admissible.begin(), admissible.end(), source_branch) ==
            admissible.end()) {
        return {};
    }
    if (admissible.size() == 1U) {
        return {Pr76RootSide::lower_admissible,
                Pr76RootSide::upper_admissible};
    }
    if (source_branch == admissible.front()) {
        return {Pr76RootSide::lower_admissible};
    }
    if (source_branch == admissible.back()) {
        return {Pr76RootSide::upper_admissible};
    }
    return {};
}

[[nodiscard]] inline bool pr76_negative_final_witness(
    const StabilityTrial& trial,
    const StabilityOptions& options) {
    return trial.status == StabilityTrialStatus::negative_tpd && trial.point &&
           stability_negative(*trial.point, options);
}

[[nodiscard]] inline std::vector<std::vector<double>> pr76_surviving_phase_starts(
    const PtThreePhaseResult& equilibrium) {
    std::vector<std::vector<double>> starts;
    if (!equilibrium.point || !equilibrium.disappearing_phase) { return starts; }
    for (std::size_t phase = 0; phase < 3U; ++phase) {
        if (phase != *equilibrium.disappearing_phase) {
            starts.push_back(equilibrium.point->phases[phase].composition);
        }
    }
    return starts;
}

} // namespace detail

[[nodiscard]] inline Pr76PtMax3Result solve_pr76_pt_max3(
    double pressure_pa, double temperature_k,
    std::span<const double> feed,
    Pr76VleEvaluator& evaluator,
    Pr76PtMax3Options options = {},
    std::span<const std::vector<double>> initial_starts = {},
    std::span<const std::vector<double>> final_starts = {}) {
    detail::pr76_pt_max3_check_options(options);
    detail::pr76_pt_max3_check_starts(feed, options);
    Pr76PtMax3Result result;
    result.options = options;
    result.base = solve_pr76_pt_vle(
        pressure_pa, temperature_k, feed, evaluator,
        options.two_phase, initial_starts, final_starts);

    switch (result.base.solution.status) {
    case PtSplitStatus::single_phase_no_instability_found:
        result.status = Pr76PtMax3Status::single_phase;
        result.diagnostic = result.base.solution.diagnostic;
        return result;
    case PtSplitStatus::two_phase_no_instability_found:
        result.status = Pr76PtMax3Status::two_phase;
        result.diagnostic = result.base.solution.diagnostic;
        return result;
    case PtSplitStatus::indeterminate:
        result.status = Pr76PtMax3Status::indeterminate;
        result.diagnostic = result.base.solution.diagnostic;
        return result;
    case PtSplitStatus::phase_set_unstable:
        break;
    }

    const auto* pair = result.base.solution.candidate();
    if (pair == nullptr || !result.base.solution.final_stability) {
        result.status = Pr76PtMax3Status::indeterminate;
        result.diagnostic =
            "PR76 max3: unstable two-phase status lacks owned pair/final stability evidence";
        return result;
    }

    const auto& final_search = *result.base.solution.final_stability;
    const auto& model = evaluator.model();
    const auto& normalized_feed = result.base.solution.initial_stability.feed;
    bool saw_higher = false;
    bool saw_boundary = false;
    bool saw_indeterminate = false;

    const auto run_attempt = [&](
        std::span<const double> log_k_1,
        std::span<const double> log_k_2,
        std::array<double, 2> fractions,
        std::array<Pr76RootSide, 3> root_sides,
        std::optional<std::size_t> witness_trial,
        std::optional<std::size_t> supplied_start) {
        if (result.attempts.size() >= options.max_three_phase_attempts) {
            result.attempt_limit_reached = true;
            return false;
        }
        Pr76PtThreePhaseAttempt attempt;
        attempt.witness_trial = witness_trial;
        attempt.supplied_start = supplied_start;
        attempt.root_sides = root_sides;
        Pr76ThreePhaseEvaluator phase_evaluator(
            model, root_sides, evaluator.root_options());
        attempt.equilibrium = iterate_pt_three_phase(
            pressure_pa, temperature_k, normalized_feed,
            log_k_1, log_k_2, fractions,
            phase_evaluator, options.three_phase);

        if (attempt.equilibrium.status == PtThreePhaseStatus::phase_disappearance) {
            saw_boundary = true;
            const auto starts = detail::pr76_surviving_phase_starts(
                attempt.equilibrium);
            if (starts.size() == 2U) {
                attempt.boundary_neighbor = solve_pr76_pt_vle(
                    pressure_pa, temperature_k, normalized_feed,
                    evaluator, options.two_phase, starts, starts);
                attempt.accepted_two_phase_neighbor =
                    attempt.boundary_neighbor->solution.status ==
                        PtSplitStatus::two_phase_no_instability_found;
            }
            result.attempts.push_back(std::move(attempt));
            const std::size_t index = result.attempts.size() - 1U;
            if (result.attempts[index].accepted_two_phase_neighbor) {
                result.status = Pr76PtMax3Status::two_phase;
                result.selected_attempt = index;
                result.diagnostic =
                    "PR76 max3: three-phase disappearance was fresh-resolved as a stable two-phase neighbor";
                return true;
            }
            return false;
        }

        if (!attempt.equilibrium.candidate_admissible()) {
            saw_indeterminate = true;
            result.attempts.push_back(std::move(attempt));
            return false;
        }

        const auto& candidate = *attempt.equilibrium.candidate();
        attempt.common_reference_allowance = candidate.chemical_potential_norm;
        auto final_options = options.final_three_phase_stability;
        final_options.tpd_tolerance += attempt.common_reference_allowance;
        std::vector<std::vector<double>> starts{
            candidate.phases[0].composition,
            candidate.phases[1].composition,
            candidate.phases[2].composition};
        attempt.final_stability = test_pt_stability_against(
            pressure_pa, temperature_k, normalized_feed,
            candidate.common_log_activity, evaluator,
            final_options, starts);
        if (attempt.final_stability->status == StabilityStatus::unstable) {
            saw_higher = true;
        } else if (attempt.final_stability->status ==
                   StabilityStatus::no_instability_found) {
            attempt.accepted_three_phase = true;
        } else {
            saw_indeterminate = true;
        }
        result.attempts.push_back(std::move(attempt));
        const std::size_t index = result.attempts.size() - 1U;
        if (result.attempts[index].accepted_three_phase) {
            result.status = Pr76PtMax3Status::three_phase;
            result.selected_attempt = index;
            result.diagnostic =
                "PR76 max3: additional-phase evidence produced a material-balanced three-phase state whose common-tangent final review found no further instability";
            return true;
        }
        return false;
    };

    // Caller-supplied triples are initialization hints only. They are ignored
    // unless the established two-phase final review has already proved the pair
    // unstable (the branch above). Enumerate numerical root sides because a
    // one-root start does not define which side should be followed if roots split.
    for (std::size_t start_index = 0;
         start_index < options.three_phase_starts.size(); ++start_index) {
        const auto& start = options.three_phase_starts[start_index];
        const auto log_k_1 = detail::pr76_three_phase_log_ratio(
            start.compositions[1], start.compositions[0], normalized_feed);
        const auto log_k_2 = detail::pr76_three_phase_log_ratio(
            start.compositions[2], start.compositions[0], normalized_feed);
        if (!log_k_1 || !log_k_2) { continue; }
        for (const Pr76RootSide side0 :
             {Pr76RootSide::lower_admissible, Pr76RootSide::upper_admissible}) {
            for (const Pr76RootSide side1 :
                 {Pr76RootSide::lower_admissible, Pr76RootSide::upper_admissible}) {
                for (const Pr76RootSide side2 :
                     {Pr76RootSide::lower_admissible, Pr76RootSide::upper_admissible}) {
                    if (run_attempt(
                            *log_k_1, *log_k_2, start.phase_fraction_seed,
                            {side0, side1, side2}, std::nullopt, start_index)) {
                        return result;
                    }
                    if (result.attempt_limit_reached) { break; }
                }
                if (result.attempt_limit_reached) { break; }
            }
            if (result.attempt_limit_reached) { break; }
        }
        if (result.attempt_limit_reached) { break; }
    }

    // Negative final-stability witnesses remain the automatic 2->3 seed source.
    // These attempts are still valuable when no explicit continuation hint exists.
    if (!result.attempt_limit_reached) {
        const auto source_log_k_1 = detail::pr76_three_phase_log_ratio(
            pair->fractions.vapor, pair->fractions.liquid, normalized_feed);
        if (source_log_k_1) {
            for (std::size_t witness_index = 0;
                 witness_index < final_search.trials.size(); ++witness_index) {
                const auto& witness = final_search.trials[witness_index];
                if (!detail::pr76_negative_final_witness(
                        witness, final_search.options)) {
                    continue;
                }
                const auto log_k_2 = detail::pr76_three_phase_log_ratio(
                    witness.point->composition, pair->fractions.liquid,
                    normalized_feed);
                if (!log_k_2) { continue; }
                const auto witness_sides = detail::pr76_witness_root_sides(
                    pressure_pa, temperature_k, witness.point->composition,
                    witness.point->branch, model, evaluator.root_options());
                for (const Pr76RootSide witness_side : witness_sides) {
                    const double beta_v = pair->fractions.vapor_fraction;
                    const std::array<double, 2> fractions{
                        (1.0 - options.new_phase_seed_fraction) * beta_v,
                        options.new_phase_seed_fraction};
                    if (run_attempt(
                            *source_log_k_1, *log_k_2, fractions,
                            {Pr76RootSide::lower_admissible,
                             Pr76RootSide::upper_admissible,
                             witness_side},
                            witness_index, std::nullopt)) {
                        return result;
                    }
                    if (result.attempt_limit_reached) { break; }
                }
                if (result.attempt_limit_reached) { break; }
            }
        }
    }

    if (saw_higher) {
        result.status = Pr76PtMax3Status::higher_phase_count_or_wrong_candidate;
        result.diagnostic =
            "PR76 max3: three-phase candidate final review found additional lower-Gibbs evidence";
    } else if (saw_boundary) {
        result.status = Pr76PtMax3Status::phase_boundary_unresolved;
        result.diagnostic =
            "PR76 max3: three-phase iteration reached a disappearance boundary but the fresh two-phase neighbor did not close";
    } else {
        result.status = Pr76PtMax3Status::indeterminate;
        result.diagnostic = result.attempt_limit_reached
            ? "PR76 max3: three-phase attempt quota exhausted without an accepted topology"
            : (saw_indeterminate
                ? "PR76 max3: additional-phase evidence exists but no three-phase attempt completed all equilibrium/stability gates"
                : "PR76 max3: unstable two-phase set lacks a usable additional-phase seed");
    }
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_PR76_THREE_PHASE_HPP
