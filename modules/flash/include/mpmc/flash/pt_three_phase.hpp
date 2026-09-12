#ifndef MPMC_FLASH_PT_THREE_PHASE_HPP
#define MPMC_FLASH_PT_THREE_PHASE_HPP

#include <mpmc/flash/pt_stability.hpp>

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

inline constexpr const char* pt_three_phase_convention =
    "PT/three-phase/logK-SSI/generalized-RR/v1";

struct PtThreePhaseOptions {
    double chemical_potential_tolerance{1e-11};
    double mass_absolute_tolerance{1e-12};
    double mass_relative_tolerance{1e-10};
    double balance_tolerance{2e-13};
    double minimum_phase_fraction{1e-10};
    double log_composition_separation{1e-7};
    double max_log_step{2.0};
    double residual_decrease{1e-4};
    int max_iterations{512};
    int max_backtracks{32};
    int max_balance_iterations{192};
    int max_balance_backtracks{48};
    std::size_t max_evaluations{30000};
    std::size_t max_components{256};
};

enum class PtThreePhaseStatus {
    converged_candidate,
    no_resolved_balance_state,
    balance_solver_degenerate,
    balance_iteration_limit,
    balance_line_search_failed,
    balance_failure,
    phase_disappearance,
    indistinguishable_phases,
    property_failure,
    iteration_limit,
    evaluation_limit,
    line_search_failed
};

struct PtThreePhaseProperty {
    StabilityPhase activity;
    double z{};
};

struct PtThreePhasePhase {
    double mole_phase_fraction{};
    std::vector<double> composition;
    StabilityPhase activity;
    double z{};
};

struct PtThreePhaseState {
    std::array<PtThreePhasePhase, 3> phases;
    std::vector<double> log_k_1;
    std::vector<double> log_k_2;
    std::vector<double> chemical_potential_residual_1;
    std::vector<double> chemical_potential_residual_2;
    std::vector<double> common_log_activity;
    double chemical_potential_norm{};
    double reduced_gibbs{};
    double gibbs_roundoff_guard{};
    double generalized_rr_residual{};
    std::array<double, 3> raw_phase_sums{};
    double mass_absolute{};
    double mass_relative{};
    int balance_iterations{};
    std::array<double, 3> pairwise_log_distance{};
};

struct PtThreePhaseResult {
    static constexpr bool final_stability_checked = false;
    static constexpr bool global_stability_proven = false;
    static constexpr const char* convention = pt_three_phase_convention;

    PtThreePhaseStatus status{PtThreePhaseStatus::no_resolved_balance_state};
    PtThreePhaseOptions options;
    double pressure_pa{};
    double temperature_k{};
    double input_feed_sum{};
    std::vector<double> feed;
    std::vector<double> initial_log_k_1;
    std::vector<double> initial_log_k_2;
    std::array<double, 2> initial_phase_fractions{}; // beta_1, beta_2; beta_0 is dependent.
    std::optional<PtThreePhaseState> point;
    std::optional<std::size_t> disappearing_phase;
    int iterations{};
    std::size_t evaluations{};
    std::size_t backtracks{};
    std::size_t rejected_evaluations{};
    std::optional<StabilityPropertyIssue> property_issue;
    std::string diagnostic;

    [[nodiscard]] bool equations_converged() const noexcept {
        return point.has_value() &&
               point->chemical_potential_norm <= options.chemical_potential_tolerance &&
               point->generalized_rr_residual <= options.balance_tolerance &&
               point->mass_absolute <= options.mass_absolute_tolerance &&
               point->mass_relative <= options.mass_relative_tolerance;
    }

    [[nodiscard]] bool candidate_admissible() const noexcept {
        return status == PtThreePhaseStatus::converged_candidate &&
               equations_converged();
    }

    [[nodiscard]] const PtThreePhaseState* candidate() const & noexcept {
        return candidate_admissible() ? &*point : nullptr;
    }
    const PtThreePhaseState* candidate() const && = delete;
};

namespace detail {

inline void pt_three_phase_check_options(const PtThreePhaseOptions& options) {
    const auto positive = [](double value) {
        return std::isfinite(value) && value > 0.0;
    };
    if (!positive(options.chemical_potential_tolerance) ||
        !positive(options.mass_absolute_tolerance) ||
        !positive(options.mass_relative_tolerance) ||
        !positive(options.balance_tolerance) ||
        !positive(options.minimum_phase_fraction) ||
        options.minimum_phase_fraction >= 1.0 / 3.0 ||
        !positive(options.log_composition_separation) ||
        !positive(options.max_log_step) ||
        !positive(options.residual_decrease) || options.residual_decrease >= 1.0 ||
        options.max_iterations < 0 || options.max_backtracks <= 0 ||
        options.max_balance_iterations < 0 || options.max_balance_backtracks <= 0 ||
        options.max_evaluations == 0U || options.max_components == 0U) {
        throw std::invalid_argument(
            "PT three-phase: invalid tolerance, step or resource option");
    }
}

class PtThreePhaseStepError : public std::runtime_error {
public:
    PtThreePhaseStepError(PtThreePhaseStatus status, const char* message)
        : std::runtime_error(message), status_(status) {}
    [[nodiscard]] PtThreePhaseStatus status() const noexcept { return status_; }
private:
    PtThreePhaseStatus status_;
};

enum class GeneralizedRr3Status {
    success, degenerate, iteration_limit, line_search_failed, unrepresentable
};

struct GeneralizedRr3Result {
    GeneralizedRr3Status status{GeneralizedRr3Status::unrepresentable};
    std::array<double, 2> fractions{};
    std::array<std::vector<double>, 3> compositions;
    double residual{std::numeric_limits<double>::infinity()};
    std::array<double, 3> raw_sums{};
    double mass_absolute{};
    double mass_relative{};
    int iterations{};
};

struct GeneralizedRr3Evaluation {
    std::array<double, 2> residual{};
    std::array<double, 3> jacobian{}; // j00, j01, j11
    std::array<std::vector<double>, 3> compositions;
};

[[nodiscard]] inline bool rr3_fractions_feasible(
    std::array<double, 2> beta) noexcept {
    return std::isfinite(beta[0]) && std::isfinite(beta[1]) &&
           beta[0] >= 0.0 && beta[1] >= 0.0 &&
           beta[0] + beta[1] <= 1.0;
}

[[nodiscard]] inline std::optional<GeneralizedRr3Evaluation> rr3_evaluate(
    std::span<const double> feed,
    std::span<const double> log_k_1,
    std::span<const double> log_k_2,
    std::array<double, 2> beta) {
    if (!rr3_fractions_feasible(beta)) { return std::nullopt; }
    const double beta_0 = 1.0 - beta[0] - beta[1];
    GeneralizedRr3Evaluation value;
    for (auto& composition : value.compositions) {
        composition.assign(feed.size(), 0.0);
    }
    double f0 = 0.0, f0_correction = 0.0;
    double f1 = 0.0, f1_correction = 0.0;
    double j00 = 0.0, j00_correction = 0.0;
    double j01 = 0.0, j01_correction = 0.0;
    double j11 = 0.0, j11_correction = 0.0;

    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0.0) { continue; }
        const double shift = std::max({0.0, log_k_1[i], log_k_2[i]});
        const double e0 = std::exp(-shift);
        const double e1 = std::exp(log_k_1[i] - shift);
        const double e2 = std::exp(log_k_2[i] - shift);
        const double denominator =
            beta_0 * e0 + beta[0] * e1 + beta[1] * e2;
        if (!std::isfinite(denominator) || !(denominator > 0.0)) {
            return std::nullopt;
        }
        const double x0 = feed[i] * e0 / denominator;
        const double x1 = feed[i] * e1 / denominator;
        const double x2 = feed[i] * e2 / denominator;
        if (!std::isfinite(x0) || !std::isfinite(x1) || !std::isfinite(x2) ||
            !(x0 > 0.0) || !(x1 > 0.0) || !(x2 > 0.0)) {
            return std::nullopt;
        }
        value.compositions[0][i] = x0;
        value.compositions[1][i] = x1;
        value.compositions[2][i] = x2;
        const double d1 = x1 - x0;
        const double d2 = x2 - x0;
        stability_add(d1, f0, f0_correction);
        stability_add(d2, f1, f1_correction);
        stability_add(-(d1 * d1) / feed[i], j00, j00_correction);
        stability_add(-(d1 * d2) / feed[i], j01, j01_correction);
        stability_add(-(d2 * d2) / feed[i], j11, j11_correction);
    }
    value.residual = {f0, f1};
    value.jacobian = {j00, j01, j11};
    return value;
}

[[nodiscard]] inline GeneralizedRr3Result solve_generalized_rr3(
    std::span<const double> feed,
    std::span<const double> log_k_1,
    std::span<const double> log_k_2,
    std::array<double, 2> initial_fractions,
    const PtThreePhaseOptions& options) {
    GeneralizedRr3Result result;
    if (feed.empty() || feed.size() > options.max_components ||
        feed.size() != log_k_1.size() || feed.size() != log_k_2.size() ||
        !rr3_fractions_feasible(initial_fractions)) {
        return result;
    }
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (!std::isfinite(log_k_1[i]) || !std::isfinite(log_k_2[i])) {
            return result;
        }
    }

    auto beta = initial_fractions;
    std::optional<GeneralizedRr3Evaluation> evaluation;
    for (int iteration = 0;; ++iteration) {
        evaluation = rr3_evaluate(feed, log_k_1, log_k_2, beta);
        if (!evaluation) { return result; }
        result.iterations = iteration;
        const double norm = std::max(
            std::abs(evaluation->residual[0]),
            std::abs(evaluation->residual[1]));
        if (!std::isfinite(norm)) { return result; }
        if (norm <= options.balance_tolerance) {
            result.status = GeneralizedRr3Status::success;
            result.residual = norm;
            break;
        }
        if (iteration >= options.max_balance_iterations) {
            result.status = GeneralizedRr3Status::iteration_limit;
            return result;
        }
        const double j00 = evaluation->jacobian[0];
        const double j01 = evaluation->jacobian[1];
        const double j11 = evaluation->jacobian[2];
        const double determinant = j00 * j11 - j01 * j01;
        const double determinant_scale =
            std::abs(j00 * j11) + std::abs(j01 * j01) + 1.0;
        if (!std::isfinite(determinant) ||
            std::abs(determinant) <= 256.0 * stability_eps * determinant_scale) {
            result.status = GeneralizedRr3Status::degenerate;
            return result;
        }
        const double f0 = evaluation->residual[0];
        const double f1 = evaluation->residual[1];
        const std::array<double, 2> step{
            (-f0 * j11 + j01 * f1) / determinant,
            (j01 * f0 - j00 * f1) / determinant};
        if (!std::isfinite(step[0]) || !std::isfinite(step[1])) { return result; }

        bool accepted = false;
        double alpha = 1.0;
        for (int backtrack = 0; backtrack < options.max_balance_backtracks;
             ++backtrack, alpha *= 0.5) {
            const std::array<double, 2> candidate{
                beta[0] + alpha * step[0], beta[1] + alpha * step[1]};
            if (!rr3_fractions_feasible(candidate)) { continue; }
            auto next = rr3_evaluate(feed, log_k_1, log_k_2, candidate);
            if (!next) { continue; }
            const double next_norm = std::max(
                std::abs(next->residual[0]), std::abs(next->residual[1]));
            if (!std::isfinite(next_norm)) { continue; }
            if (next_norm < norm || next_norm <= options.balance_tolerance) {
                beta = candidate;
                accepted = true;
                break;
            }
        }
        if (!accepted) {
            result.status = GeneralizedRr3Status::line_search_failed;
            return result;
        }
    }

    evaluation = rr3_evaluate(feed, log_k_1, log_k_2, beta);
    if (!evaluation) {
        result.status = GeneralizedRr3Status::unrepresentable;
        return result;
    }
    result.fractions = beta;
    result.compositions = std::move(evaluation->compositions);
    result.raw_sums = {
        stability_sum(result.compositions[0]),
        stability_sum(result.compositions[1]),
        stability_sum(result.compositions[2])};
    for (std::size_t p = 0; p < 3U; ++p) {
        if (!std::isfinite(result.raw_sums[p]) || !(result.raw_sums[p] > 0.0)) {
            result.status = GeneralizedRr3Status::unrepresentable;
            return result;
        }
        for (double& value : result.compositions[p]) { value /= result.raw_sums[p]; }
    }
    const std::array<double, 3> fractions{
        1.0 - beta[0] - beta[1], beta[0], beta[1]};
    double absolute_error = 0.0;
    double relative_error = 0.0;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        const double recovered =
            fractions[0] * result.compositions[0][i] +
            fractions[1] * result.compositions[1][i] +
            fractions[2] * result.compositions[2][i];
        const double error = std::abs(recovered - feed[i]);
        absolute_error = std::max(absolute_error, error);
        if (feed[i] > 0.0) {
            relative_error = std::max(relative_error, error / feed[i]);
        } else if (recovered != 0.0) {
            relative_error = std::numeric_limits<double>::infinity();
        }
    }
    result.mass_absolute = absolute_error;
    result.mass_relative = relative_error;
    return result;
}

[[nodiscard]] inline double pt_three_phase_log_distance(
    std::span<const double> first,
    std::span<const double> second,
    std::span<const double> feed) {
    if (first.size() != second.size() || first.size() != feed.size()) {
        return std::numeric_limits<double>::infinity();
    }
    double distance = 0.0;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0.0) {
            if (first[i] != 0.0 || second[i] != 0.0) {
                return std::numeric_limits<double>::infinity();
            }
            continue;
        }
        if (!(first[i] > 0.0) || !(second[i] > 0.0)) {
            return std::numeric_limits<double>::infinity();
        }
        const double value = std::abs(std::log(first[i]) - std::log(second[i]));
        if (!std::isfinite(value)) { return std::numeric_limits<double>::infinity(); }
        distance = std::max(distance, value);
    }
    return distance;
}

inline void pt_three_phase_check_property(
    const PtThreePhaseProperty& property, std::size_t component_count) {
    if (property.activity.ln_phi.size() != component_count ||
        !property.activity.smooth || !std::isfinite(property.z) ||
        !(property.z > 0.0)) {
        throw StabilityPropertyError(
            StabilityPropertyIssue::nonfinite_properties,
            "PT three-phase: invalid phase property payload");
    }
    for (double value : property.activity.ln_phi) {
        if (!std::isfinite(value)) {
            throw StabilityPropertyError(
                StabilityPropertyIssue::nonfinite_properties,
                "PT three-phase: nonfinite fugacity coefficient");
        }
    }
}

} // namespace detail

template <typename Provider>
[[nodiscard]] PtThreePhaseResult iterate_pt_three_phase(
    double pressure_pa, double temperature_k,
    std::span<const double> feed,
    std::span<const double> initial_log_k_1,
    std::span<const double> initial_log_k_2,
    std::array<double, 2> initial_phase_fractions,
    Provider&& provider,
    PtThreePhaseOptions options = {}) {
    detail::pt_three_phase_check_options(options);
    if (!std::isfinite(pressure_pa) || !(pressure_pa > 0.0) ||
        !std::isfinite(temperature_k) || !(temperature_k > 0.0)) {
        throw std::domain_error("PT three-phase: finite p>0 Pa and T>0 K required");
    }
    if (feed.empty() || feed.size() > options.max_components ||
        feed.size() != initial_log_k_1.size() ||
        feed.size() != initial_log_k_2.size()) {
        throw std::invalid_argument(
            "PT three-phase: feed/logK dimension mismatch or component quota exceeded");
    }
    if (!detail::rr3_fractions_feasible(initial_phase_fractions)) {
        throw std::domain_error(
            "PT three-phase: initial independent phase fractions must lie in the closed simplex");
    }
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (!std::isfinite(initial_log_k_1[i]) ||
            !std::isfinite(initial_log_k_2[i])) {
            throw std::domain_error("PT three-phase: finite initial logK required");
        }
    }

    PtThreePhaseResult result;
    result.options = options;
    result.pressure_pa = pressure_pa;
    result.temperature_k = temperature_k;
    result.input_feed_sum = detail::stability_check_composition(feed);
    result.feed = detail::stability_normalize(feed, result.input_feed_sum);
    result.initial_log_k_1.assign(initial_log_k_1.begin(), initial_log_k_1.end());
    result.initial_log_k_2.assign(initial_log_k_2.begin(), initial_log_k_2.end());
    result.initial_phase_fractions = initial_phase_fractions;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] > 0.0 && result.feed[i] == 0.0) {
            result.status = PtThreePhaseStatus::balance_failure;
            result.diagnostic = "PT three-phase: roundoff normalization lost an active feed component";
            return result;
        }
    }

    const auto evaluate_property = [&](std::span<const double> composition,
                                       std::size_t slot) {
        if (result.evaluations >= options.max_evaluations) {
            throw detail::PtThreePhaseStepError(
                PtThreePhaseStatus::evaluation_limit,
                "PT three-phase: property budget exhausted");
        }
        ++result.evaluations;
        auto property = provider(pressure_pa, temperature_k, composition, slot);
        detail::pt_three_phase_check_property(property, result.feed.size());
        return property;
    };

    const auto evaluate_state = [&](std::span<const double> log_k_1,
                                    std::span<const double> log_k_2,
                                    std::array<double, 2> fraction_seed) {
        const auto balance = detail::solve_generalized_rr3(
            result.feed, log_k_1, log_k_2, fraction_seed, options);
        switch (balance.status) {
        case detail::GeneralizedRr3Status::success: break;
        case detail::GeneralizedRr3Status::degenerate:
            throw detail::PtThreePhaseStepError(
                PtThreePhaseStatus::balance_solver_degenerate,
                "PT three-phase: generalized RR Jacobian is degenerate");
        case detail::GeneralizedRr3Status::iteration_limit:
            throw detail::PtThreePhaseStepError(
                PtThreePhaseStatus::balance_iteration_limit,
                "PT three-phase: generalized RR iteration limit");
        case detail::GeneralizedRr3Status::line_search_failed:
            throw detail::PtThreePhaseStepError(
                PtThreePhaseStatus::balance_line_search_failed,
                "PT three-phase: generalized RR line search failed");
        case detail::GeneralizedRr3Status::unrepresentable:
            throw detail::PtThreePhaseStepError(
                PtThreePhaseStatus::no_resolved_balance_state,
                "PT three-phase: no representable generalized RR state");
        }

        PtThreePhaseState state;
        const std::array<double, 3> fractions{
            1.0 - balance.fractions[0] - balance.fractions[1],
            balance.fractions[0], balance.fractions[1]};
        for (std::size_t phase = 0; phase < 3U; ++phase) {
            state.phases[phase].mole_phase_fraction = fractions[phase];
            state.phases[phase].composition = balance.compositions[phase];
            const auto property = evaluate_property(
                state.phases[phase].composition, phase);
            state.phases[phase].activity = property.activity;
            state.phases[phase].z = property.z;
            state.raw_phase_sums[phase] = balance.raw_sums[phase];
        }
        state.generalized_rr_residual = balance.residual;
        state.mass_absolute = balance.mass_absolute;
        state.mass_relative = balance.mass_relative;
        state.balance_iterations = balance.iterations;

        const std::size_t n = result.feed.size();
        state.log_k_1.assign(n, 0.0);
        state.log_k_2.assign(n, 0.0);
        state.chemical_potential_residual_1.assign(n, 0.0);
        state.chemical_potential_residual_2.assign(n, 0.0);
        state.common_log_activity.assign(n, 0.0);
        double gibbs_correction = 0.0;
        double magnitude = 1.0;
        double magnitude_correction = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            if (result.feed[i] == 0.0) { continue; }
            const double l0 = std::log(state.phases[0].composition[i]);
            const double l1 = std::log(state.phases[1].composition[i]);
            const double l2 = std::log(state.phases[2].composition[i]);
            const double mu0 = l0 + state.phases[0].activity.ln_phi[i];
            const double mu1 = l1 + state.phases[1].activity.ln_phi[i];
            const double mu2 = l2 + state.phases[2].activity.ln_phi[i];
            state.log_k_1[i] = l1 - l0;
            state.log_k_2[i] = l2 - l0;
            state.chemical_potential_residual_1[i] = mu0 - mu1;
            state.chemical_potential_residual_2[i] = mu0 - mu2;
            state.common_log_activity[i] = (mu0 + mu1 + mu2) / 3.0;
            state.chemical_potential_norm = std::max({
                state.chemical_potential_norm,
                std::abs(state.chemical_potential_residual_1[i]),
                std::abs(state.chemical_potential_residual_2[i])});
            for (std::size_t phase = 0; phase < 3U; ++phase) {
                const double lx = std::log(state.phases[phase].composition[i]);
                const double mu = lx + state.phases[phase].activity.ln_phi[i];
                detail::stability_add(
                    state.phases[phase].mole_phase_fraction *
                        state.phases[phase].composition[i] * mu,
                    state.reduced_gibbs, gibbs_correction);
                detail::stability_add(
                    state.phases[phase].mole_phase_fraction *
                        state.phases[phase].composition[i] *
                        (std::abs(lx) + std::abs(state.phases[phase].activity.ln_phi[i])),
                    magnitude, magnitude_correction);
            }
        }
        state.gibbs_roundoff_guard = 256.0 * detail::stability_eps * magnitude;
        if (!std::isfinite(state.chemical_potential_norm) ||
            !std::isfinite(state.reduced_gibbs) ||
            !std::isfinite(state.gibbs_roundoff_guard)) {
            throw StabilityPropertyError(
                StabilityPropertyIssue::nonfinite_properties,
                "PT three-phase: nonrepresentable activity/Gibbs arithmetic");
        }
        return state;
    };

    try {
        result.point = evaluate_state(
            result.initial_log_k_1, result.initial_log_k_2,
            result.initial_phase_fractions);
    } catch (const detail::PtThreePhaseStepError& error) {
        result.status = error.status();
        result.diagnostic = error.what();
        return result;
    } catch (const StabilityPropertyError& error) {
        result.status = PtThreePhaseStatus::property_failure;
        result.property_issue = error.issue();
        result.diagnostic = error.what();
        return result;
    }

    for (;;) {
        auto& point = *result.point;
        if (point.chemical_potential_norm <= options.chemical_potential_tolerance) {
            if (point.generalized_rr_residual > options.balance_tolerance ||
                point.mass_absolute > options.mass_absolute_tolerance ||
                point.mass_relative > options.mass_relative_tolerance) {
                result.status = PtThreePhaseStatus::balance_failure;
                result.diagnostic =
                    "PT three-phase: chemical potentials converged but material balance failed";
                return result;
            }
            for (std::size_t phase = 0; phase < 3U; ++phase) {
                if (point.phases[phase].mole_phase_fraction <=
                    options.minimum_phase_fraction) {
                    result.status = PtThreePhaseStatus::phase_disappearance;
                    result.disappearing_phase = phase;
                    result.diagnostic =
                        "PT three-phase: a phase reached the disappearance boundary; a fresh neighboring-topology solve is required";
                    return result;
                }
            }
            point.pairwise_log_distance = {
                detail::pt_three_phase_log_distance(
                    point.phases[0].composition, point.phases[1].composition, result.feed),
                detail::pt_three_phase_log_distance(
                    point.phases[0].composition, point.phases[2].composition, result.feed),
                detail::pt_three_phase_log_distance(
                    point.phases[1].composition, point.phases[2].composition, result.feed)};
            if (!(point.pairwise_log_distance[0] > options.log_composition_separation) ||
                !(point.pairwise_log_distance[1] > options.log_composition_separation) ||
                !(point.pairwise_log_distance[2] > options.log_composition_separation)) {
                result.status = PtThreePhaseStatus::indistinguishable_phases;
                result.diagnostic =
                    "PT three-phase: converged phase compositions are not pairwise distinct";
                return result;
            }
            result.status = PtThreePhaseStatus::converged_candidate;
            result.diagnostic =
                "PT three-phase chemical-potential equalities and material balance converged; final phase-set stability is not yet certified";
            return result;
        }
        if (result.iterations >= options.max_iterations) {
            result.status = PtThreePhaseStatus::iteration_limit;
            result.diagnostic = "PT three-phase: iteration limit";
            return result;
        }
        const double scale = std::max(
            1.0, point.chemical_potential_norm / options.max_log_step);
        bool accepted = false;
        double alpha = 1.0;
        for (int backtrack = 0; backtrack < options.max_backtracks;
             ++backtrack, alpha *= 0.5) {
            if (backtrack > 0) { ++result.backtracks; }
            std::vector<double> next_log_k_1(result.feed.size(), 0.0);
            std::vector<double> next_log_k_2(result.feed.size(), 0.0);
            for (std::size_t i = 0; i < result.feed.size(); ++i) {
                next_log_k_1[i] = point.log_k_1[i] +
                    (alpha / scale) * point.chemical_potential_residual_1[i];
                next_log_k_2[i] = point.log_k_2[i] +
                    (alpha / scale) * point.chemical_potential_residual_2[i];
            }
            const std::array<double, 2> fraction_seed{
                point.phases[1].mole_phase_fraction,
                point.phases[2].mole_phase_fraction};
            try {
                auto next = evaluate_state(next_log_k_1, next_log_k_2, fraction_seed);
                const double required = options.residual_decrease *
                    (alpha / scale) * point.chemical_potential_norm;
                const bool residual_progress =
                    next.chemical_potential_norm <= options.chemical_potential_tolerance ||
                    point.chemical_potential_norm - next.chemical_potential_norm > required;
                const bool gibbs_ok =
                    next.reduced_gibbs <= point.reduced_gibbs +
                        point.gibbs_roundoff_guard + next.gibbs_roundoff_guard;
                if (residual_progress && gibbs_ok) {
                    result.point = std::move(next);
                    accepted = true;
                    break;
                }
                ++result.rejected_evaluations;
            } catch (const detail::PtThreePhaseStepError& error) {
                ++result.rejected_evaluations;
                result.diagnostic = error.what();
                if (error.status() == PtThreePhaseStatus::evaluation_limit) {
                    result.status = error.status();
                    return result;
                }
            } catch (const StabilityPropertyError& error) {
                ++result.rejected_evaluations;
                result.property_issue = error.issue();
                result.diagnostic = error.what();
            }
        }
        if (!accepted) {
            result.status = PtThreePhaseStatus::line_search_failed;
            if (result.diagnostic.empty()) {
                result.diagnostic = "PT three-phase: logK line search failed";
            }
            return result;
        }
        ++result.iterations;
    }
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_PT_THREE_PHASE_HPP
