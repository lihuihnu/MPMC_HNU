#ifndef MPMC_FLASH_CPA_THREE_PHASE_HPP
#define MPMC_FLASH_CPA_THREE_PHASE_HPP

#include <mpmc/flash/cpa_split.hpp>
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

inline constexpr const char* cpa_pt_max3_convention =
    "CPA/PT/max3/TPD-VLE-additional-phase-generalized-RR/v1";

class CpaThreePhaseEvaluator {
public:
    CpaThreePhaseEvaluator(
        const thermodynamics::CpaPtPhase& model,
        std::array<CpaRootSide, 3> sides,
        thermodynamics::CpaPtOptions pt_options = {})
        : model_(model), sides_(sides), pt_options_(std::move(pt_options)) {
        if (model_.size() == 0U) {
            throw std::invalid_argument(
                "CPA three-phase evaluator: empty model snapshot");
        }
    }

    [[nodiscard]] PtThreePhaseProperty operator()(
        double pressure_pa, double temperature_k,
        std::span<const double> composition, std::size_t slot) {
        if (slot >= sides_.size()) {
            throw std::invalid_argument(
                "CPA three-phase evaluator: invalid phase slot");
        }
        auto value = detail::cpa_flash_evaluate_side(
            model_, pt_options_, pressure_pa, temperature_k,
            composition, sides_[slot]);
        return {std::move(value.activity), value.z};
    }

private:
    const thermodynamics::CpaPtPhase& model_;
    std::array<CpaRootSide, 3> sides_;
    thermodynamics::CpaPtOptions pt_options_;
};

// Caller-supplied compositions are local continuation hints only. They may be
// consumed only after the base two-phase final common-tangent review has already
// proved that the current pair is unstable; they never constitute phase-count evidence.
struct CpaPtThreePhaseStart {
    std::array<std::vector<double>, 3> compositions;
    std::array<double, 2> phase_fraction_seed{1.0 / 3.0, 1.0 / 3.0};
};

struct CpaPtMax3Options {
    PtSplitOptions two_phase;
    PtThreePhaseOptions three_phase;
    StabilityOptions final_three_phase_stability;
    // Independent bounded budgets for continuation hints and automatic
    // negative-TPD witnesses. One source class cannot starve the other.
    std::size_t max_three_phase_attempts{16};
    double new_phase_seed_fraction{0.1};
    std::vector<CpaPtThreePhaseStart> three_phase_starts;
    std::size_t max_three_phase_starts{16};
    std::size_t max_three_phase_start_entries{12288};
};

enum class CpaPtMax3Status {
    single_phase,
    two_phase,
    three_phase,
    phase_boundary_unresolved,
    higher_phase_count_or_wrong_candidate,
    indeterminate
};

struct CpaPtThreePhaseAttempt {
    std::optional<std::size_t> witness_trial;
    std::optional<std::size_t> supplied_start;
    std::array<CpaRootSide, 3> root_sides{
        CpaRootSide::upper_density_admissible,
        CpaRootSide::lower_density_admissible,
        CpaRootSide::lower_density_admissible};
    PtThreePhaseResult equilibrium;
    std::optional<StabilityResult> final_stability;
    double common_reference_allowance{};
    std::optional<CpaPtSplitResult> boundary_neighbor;
    bool accepted_three_phase{false};
    bool accepted_two_phase_neighbor{false};
};

struct CpaPtMax3Result {
    static constexpr bool global_stability_proven = false;
    static constexpr const char* convention = cpa_pt_max3_convention;

    CpaPtMax3Status status{CpaPtMax3Status::indeterminate};
    CpaPtMax3Options options;
    CpaPtSplitResult base;
    std::vector<CpaPtThreePhaseAttempt> attempts;
    std::optional<std::size_t> selected_attempt;
    bool attempt_limit_reached{false};
    std::string diagnostic;

    [[nodiscard]] const PtThreePhaseState* three_phase_candidate() const & noexcept {
        if (status != CpaPtMax3Status::three_phase || !selected_attempt ||
            *selected_attempt >= attempts.size()) {
            return nullptr;
        }
        return attempts[*selected_attempt].equilibrium.candidate();
    }
    const PtThreePhaseState* three_phase_candidate() const && = delete;

    [[nodiscard]] const CpaPtSplitResult* two_phase_neighbor() const & noexcept {
        if (status != CpaPtMax3Status::two_phase || !selected_attempt ||
            *selected_attempt >= attempts.size()) {
            return nullptr;
        }
        const auto& neighbor = attempts[*selected_attempt].boundary_neighbor;
        if (!neighbor || !attempts[*selected_attempt].accepted_two_phase_neighbor) {
            return nullptr;
        }
        return &*neighbor;
    }
    const CpaPtSplitResult* two_phase_neighbor() const && = delete;
};

namespace detail {

inline void cpa_pt_max3_check_options(const CpaPtMax3Options& options) {
    pt_three_phase_check_options(options.three_phase);
    stability_check_options(options.final_three_phase_stability);
    if (options.max_three_phase_attempts == 0U ||
        options.max_three_phase_starts == 0U ||
        options.max_three_phase_start_entries == 0U ||
        options.three_phase_starts.size() > options.max_three_phase_starts ||
        !std::isfinite(options.new_phase_seed_fraction) ||
        !(options.new_phase_seed_fraction > 0.0) ||
        !(options.new_phase_seed_fraction < 1.0)) {
        throw std::invalid_argument(
            "CPA max3: invalid attempt/start/seed option");
    }
}

inline void cpa_pt_max3_check_starts(
    std::span<const double> feed,
    const CpaPtMax3Options& options) {
    const std::size_t n = feed.size();
    if (n == 0U) { throw std::invalid_argument("CPA max3: empty feed"); }
    (void)stability_check_composition(feed);
    if (options.three_phase_starts.size() >
            std::numeric_limits<std::size_t>::max() / 3U ||
        options.three_phase_starts.size() * 3U >
            std::numeric_limits<std::size_t>::max() / n ||
        options.three_phase_starts.size() * 3U * n >
            options.max_three_phase_start_entries) {
        throw std::length_error(
            "CPA max3: three-phase start entry quota exceeded");
    }
    for (const auto& start : options.three_phase_starts) {
        if (!rr3_fractions_feasible(start.phase_fraction_seed)) {
            throw std::domain_error(
                "CPA max3: three-phase fraction seed lies outside simplex");
        }
        for (const auto& composition : start.compositions) {
            if (composition.size() != n) {
                throw std::invalid_argument(
                    "CPA max3: three-phase start dimension mismatch");
            }
            (void)stability_check_composition(composition);
            for (std::size_t i = 0; i < n; ++i) {
                if ((feed[i] > 0.0) != (composition[i] > 0.0)) {
                    throw std::domain_error(
                        "CPA max3: three-phase start support mismatch");
                }
            }
        }
    }
}

[[nodiscard]] inline std::optional<std::vector<double>> cpa_three_phase_log_ratio(
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

[[nodiscard]] inline std::vector<CpaRootSide> cpa_witness_root_sides(
    double pressure_pa, double temperature_k,
    std::span<const double> composition,
    std::size_t source_branch,
    const thermodynamics::CpaPtPhase& model,
    const thermodynamics::CpaPtOptions& pt_options) {
    const auto roots = model.roots(
        pressure_pa, temperature_k, composition, pt_options);
    if (roots.status != thermodynamics::CpaPtRootStatus::success) {
        return {};
    }
    std::vector<std::size_t> admissible;
    for (std::size_t k = 0; k < roots.roots.size(); ++k) {
        if (roots.roots[k].pressure_slope_sign > 0) {
            admissible.push_back(k);
        }
    }
    if (admissible.empty() ||
        std::find(admissible.begin(), admissible.end(), source_branch) ==
            admissible.end()) {
        return {};
    }
    if (admissible.size() == 1U) {
        return {CpaRootSide::lower_density_admissible,
                CpaRootSide::upper_density_admissible};
    }
    const auto less_density = [&](std::size_t first, std::size_t second) {
        return roots.roots[first].molar_density_mol_per_m3 <
               roots.roots[second].molar_density_mol_per_m3;
    };
    const std::size_t lower =
        *std::min_element(admissible.begin(), admissible.end(), less_density);
    const std::size_t upper =
        *std::max_element(admissible.begin(), admissible.end(), less_density);
    if (source_branch == lower) {
        return {CpaRootSide::lower_density_admissible};
    }
    if (source_branch == upper) {
        return {CpaRootSide::upper_density_admissible};
    }
    // The public two-side continuation contract cannot represent an interior
    // mechanically admissible density root without inventing another phase identity.
    return {};
}

[[nodiscard]] inline bool cpa_negative_final_witness(
    const StabilityTrial& trial,
    const StabilityOptions& options) {
    return trial.status == StabilityTrialStatus::negative_tpd && trial.point &&
           stability_negative(*trial.point, options);
}

[[nodiscard]] inline std::vector<std::vector<double>> cpa_surviving_phase_starts(
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

[[nodiscard]] inline CpaPtMax3Result solve_cpa_pt_max3(
    double pressure_pa, double temperature_k,
    std::span<const double> feed,
    CpaVleEvaluator& evaluator,
    CpaPtMax3Options options = {},
    std::span<const std::vector<double>> initial_starts = {},
    std::span<const std::vector<double>> final_starts = {}) {
    detail::cpa_pt_max3_check_options(options);
    detail::cpa_pt_max3_check_starts(feed, options);

    CpaPtMax3Result result;
    result.options = options;
    result.base = solve_cpa_pt_vle(
        pressure_pa, temperature_k, feed, evaluator,
        options.two_phase, initial_starts, final_starts);

    switch (result.base.solution.status) {
    case PtSplitStatus::single_phase_no_instability_found:
        result.status = CpaPtMax3Status::single_phase;
        result.diagnostic = result.base.solution.diagnostic;
        return result;
    case PtSplitStatus::two_phase_no_instability_found:
        result.status = CpaPtMax3Status::two_phase;
        result.diagnostic = result.base.solution.diagnostic;
        return result;
    case PtSplitStatus::indeterminate:
        result.status = CpaPtMax3Status::indeterminate;
        result.diagnostic = result.base.solution.diagnostic;
        return result;
    case PtSplitStatus::phase_set_unstable:
        break;
    }

    const auto* pair = result.base.solution.candidate();
    if (pair == nullptr || !result.base.solution.final_stability) {
        result.status = CpaPtMax3Status::indeterminate;
        result.diagnostic =
            "CPA max3: unstable two-phase status lacks owned pair/final stability evidence";
        return result;
    }

    const auto& final_search = *result.base.solution.final_stability;
    const auto& model = evaluator.model();
    const auto& normalized_feed = result.base.solution.initial_stability.feed;
    bool saw_higher = false;
    bool saw_boundary = false;
    bool saw_indeterminate = false;
    std::size_t continuation_attempts = 0U;
    std::size_t automatic_attempts = 0U;
    bool continuation_limit_reached = false;
    bool automatic_limit_reached = false;

    const auto run_attempt = [&] (
        std::span<const double> log_k_1,
        std::span<const double> log_k_2,
        std::array<double, 2> fractions,
        std::array<CpaRootSide, 3> root_sides,
        std::optional<std::size_t> witness_trial,
        std::optional<std::size_t> supplied_start) {
        const bool continuation = supplied_start.has_value();
        auto& source_attempts = continuation ? continuation_attempts : automatic_attempts;
        auto& source_limit = continuation ? continuation_limit_reached
                                          : automatic_limit_reached;
        if (source_attempts >= options.max_three_phase_attempts) {
            source_limit = true;
            result.attempt_limit_reached = true;
            return false;
        }
        ++source_attempts;

        CpaPtThreePhaseAttempt attempt;
        attempt.witness_trial = witness_trial;
        attempt.supplied_start = supplied_start;
        attempt.root_sides = root_sides;
        CpaThreePhaseEvaluator phase_evaluator(
            model, root_sides, evaluator.pt_options());
        attempt.equilibrium = iterate_pt_three_phase(
            pressure_pa, temperature_k, normalized_feed,
            log_k_1, log_k_2, fractions,
            phase_evaluator, options.three_phase);

        if (attempt.equilibrium.status == PtThreePhaseStatus::phase_disappearance) {
            saw_boundary = true;
            const auto starts = detail::cpa_surviving_phase_starts(
                attempt.equilibrium);
            if (starts.size() == 2U) {
                attempt.boundary_neighbor = solve_cpa_pt_vle(
                    pressure_pa, temperature_k, normalized_feed,
                    evaluator, options.two_phase, starts, starts);
                attempt.accepted_two_phase_neighbor =
                    attempt.boundary_neighbor->solution.status ==
                        PtSplitStatus::two_phase_no_instability_found;
            }
            result.attempts.push_back(std::move(attempt));
            const std::size_t index = result.attempts.size() - 1U;
            if (result.attempts[index].accepted_two_phase_neighbor) {
                result.status = CpaPtMax3Status::two_phase;
                result.selected_attempt = index;
                result.diagnostic =
                    "CPA max3: three-phase disappearance was fresh-resolved as a stable two-phase neighbor";
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
            result.status = CpaPtMax3Status::three_phase;
            result.selected_attempt = index;
            result.diagnostic =
                "CPA max3: additional-phase evidence produced a material-balanced three-phase state whose common-tangent final review found no further instability";
            return true;
        }
        return false;
    };

    // Supplied triples remain hints only; reaching this block already requires
    // base two-phase final instability. Each numerical density-side combination
    // is bounded independently from the automatic witness route below.
    for (std::size_t start_index = 0;
         start_index < options.three_phase_starts.size(); ++start_index) {
        const auto& start = options.three_phase_starts[start_index];
        const auto log_k_1 = detail::cpa_three_phase_log_ratio(
            start.compositions[1], start.compositions[0], normalized_feed);
        const auto log_k_2 = detail::cpa_three_phase_log_ratio(
            start.compositions[2], start.compositions[0], normalized_feed);
        if (!log_k_1 || !log_k_2) { continue; }
        for (const CpaRootSide side0 :
             {CpaRootSide::upper_density_admissible,
              CpaRootSide::lower_density_admissible}) {
            for (const CpaRootSide side1 :
                 {CpaRootSide::upper_density_admissible,
                  CpaRootSide::lower_density_admissible}) {
                for (const CpaRootSide side2 :
                     {CpaRootSide::upper_density_admissible,
                      CpaRootSide::lower_density_admissible}) {
                    if (run_attempt(
                            *log_k_1, *log_k_2,
                            start.phase_fraction_seed,
                            {side0, side1, side2},
                            std::nullopt, start_index)) {
                        return result;
                    }
                    if (continuation_limit_reached) { break; }
                }
                if (continuation_limit_reached) { break; }
            }
            if (continuation_limit_reached) { break; }
        }
        if (continuation_limit_reached) { break; }
    }

    // phase0 is the upper-density two-phase candidate; phase1 is the lower-density
    // candidate. A retained negative final-TPD witness supplies phase2.
    const auto source_log_k_1 = detail::cpa_three_phase_log_ratio(
        pair->fractions.vapor, pair->fractions.liquid, normalized_feed);
    if (source_log_k_1) {
        for (std::size_t witness_index = 0;
             witness_index < final_search.trials.size(); ++witness_index) {
            const auto& witness = final_search.trials[witness_index];
            if (!detail::cpa_negative_final_witness(
                    witness, final_search.options)) {
                continue;
            }
            const auto log_k_2 = detail::cpa_three_phase_log_ratio(
                witness.point->composition,
                pair->fractions.liquid, normalized_feed);
            if (!log_k_2) { continue; }
            const auto witness_sides = detail::cpa_witness_root_sides(
                pressure_pa, temperature_k,
                witness.point->composition, witness.point->branch,
                model, evaluator.pt_options());
            for (const CpaRootSide witness_side : witness_sides) {
                const double beta_lower = pair->fractions.vapor_fraction;
                const std::array<double, 2> fractions{
                    (1.0 - options.new_phase_seed_fraction) * beta_lower,
                    options.new_phase_seed_fraction};
                if (run_attempt(
                        *source_log_k_1, *log_k_2, fractions,
                        {CpaRootSide::upper_density_admissible,
                         CpaRootSide::lower_density_admissible,
                         witness_side},
                        witness_index, std::nullopt)) {
                    return result;
                }
                if (automatic_limit_reached) { break; }
            }
            if (automatic_limit_reached) { break; }
        }
    }

    if (saw_higher) {
        result.status = CpaPtMax3Status::higher_phase_count_or_wrong_candidate;
        result.diagnostic =
            "CPA max3: three-phase candidate final review found additional lower-Gibbs evidence";
    } else if (saw_boundary) {
        result.status = CpaPtMax3Status::phase_boundary_unresolved;
        result.diagnostic =
            "CPA max3: three-phase iteration reached a disappearance boundary but the fresh two-phase neighbor did not close";
    } else {
        result.status = CpaPtMax3Status::indeterminate;
        result.diagnostic = automatic_limit_reached
            ? "CPA max3: automatic negative-TPD witness attempt quota exhausted without an accepted topology"
            : (saw_indeterminate
                ? (continuation_limit_reached
                    ? "CPA max3: continuation-hint quota was exhausted, automatic witnesses were still evaluated, and no three-phase attempt completed all equilibrium/stability gates"
                    : "CPA max3: additional-phase evidence exists but no three-phase attempt completed all equilibrium/stability gates")
                : "CPA max3: unstable two-phase set lacks a usable additional-phase seed");
    }
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_CPA_THREE_PHASE_HPP
