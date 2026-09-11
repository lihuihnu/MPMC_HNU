#ifndef MPMC_FLASH_SW92_PHASE_ASSIGNED_THREE_PHASE_HPP
#define MPMC_FLASH_SW92_PHASE_ASSIGNED_THREE_PHASE_HPP

#include <mpmc/flash/sw92_phase_assigned_h_side_witness.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::flash {

inline constexpr std::string_view sw92_phase_assigned_c2b1_three_phase_convention =
    "SW92-equilibrium/phase-assigned-aq-na-joint/unordered-w-h0-h1-"
    "three-phase-logK-SSI-generalized-RR/v1";

struct Sw92PhaseAssignedThreePhaseOptions {
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
    double new_hydrocarbon_seed_share{0.1};
    thermodynamics::Sw92RootOptions aqueous_root_options;
    thermodynamics::Sw92RootOptions nonaqueous_root_options;
};

enum class Sw92PhaseAssignedThreePhaseStatus {
    converged_candidate,
    source_candidate_unavailable,
    source_candidate_inconsistent,
    source_witness_unavailable,
    source_witness_inconsistent,
    no_resolved_balance_state,
    balance_solver_degenerate,
    balance_iteration_limit,
    balance_line_search_failed,
    balance_failure,
    aqueous_phase_disappearance,
    hydrocarbon_phase_disappearance,
    indistinguishable_phases,
    family_root_nonsmooth,
    phase_role_indeterminate,
    phase_role_reversed,
    property_failure,
    iteration_limit,
    evaluation_limit,
    line_search_failed
};

struct Sw92PhaseAssignedThreePhasePhase {
    Sw92PhysicalPhaseRole physical_role{Sw92PhysicalPhaseRole::nonaqueous};
    thermodynamics::SwPhaseFamily thermodynamic_family{
        thermodynamics::SwPhaseFamily::nonaqueous};
    double mole_phase_fraction{};
    std::vector<double> composition;
    StabilityPhase activity;
    double compressibility_factor{};
};

struct Sw92PhaseAssignedThreePhaseState {
    Sw92PhaseAssignedThreePhasePhase aqueous_phase;
    Sw92PhaseAssignedThreePhasePhase hydrocarbon0_phase;
    Sw92PhaseAssignedThreePhasePhase hydrocarbon1_phase;
    std::vector<double> log_k_h0;
    std::vector<double> log_k_h1;
    std::vector<double> chemical_potential_residual_h0;
    std::vector<double> chemical_potential_residual_h1;
    std::vector<double> common_log_activity;
    double chemical_potential_norm{};
    double reduced_gibbs{};
    double gibbs_roundoff_guard{};
    double generalized_rr_residual{};
    double raw_aqueous_sum{};
    double raw_hydrocarbon0_sum{};
    double raw_hydrocarbon1_sum{};
    double mass_absolute{};
    double mass_relative{};
    int balance_iterations{};
    double aqueous_h0_log_distance{};
    double aqueous_h1_log_distance{};
    double h0_h1_log_distance{};
    double aqueous_minus_h0_water_fraction{
        std::numeric_limits<double>::quiet_NaN()};
    double aqueous_minus_h1_water_fraction{
        std::numeric_limits<double>::quiet_NaN()};
    double water_role_roundoff_guard{
        std::numeric_limits<double>::quiet_NaN()};
    bool hydrocarbon_slots_canonicalized{false};
};

struct Sw92PhaseAssignedThreePhaseResult {
    static constexpr std::size_t maximum_phase_count = 3;
    static constexpr bool final_stability_checked = false;
    static constexpr bool global_stability_proven = false;
    static constexpr bool accepted_phase_set_published = false;

    Sw92PhaseAssignedThreePhaseStatus status{
        Sw92PhaseAssignedThreePhaseStatus::no_resolved_balance_state};
    Sw92PhaseAssignedThreePhaseOptions options;
    double pressure_pa{};
    double temperature_k{};
    double input_feed_sum{};
    std::vector<double> feed;
    double nacl_molality_mol_per_kg_water{};
    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string model_profile{thermodynamics::sw92_corrected_profile};
    std::string phase_convention{thermodynamics::sw92_pt_convention};
    std::string equilibrium_profile{sw92_phase_assigned_aq_na_joint_profile};
    std::string primitive_convention{sw92_phase_assigned_c2b1_three_phase_convention};
    std::vector<double> initial_log_k_h0;
    std::vector<double> initial_log_k_h1;
    std::array<double, 2> initial_hydrocarbon_fractions{};
    std::optional<std::size_t> source_witness_index;
    std::optional<Sw92PhaseAssignedThreePhaseState> point;
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
        return status == Sw92PhaseAssignedThreePhaseStatus::converged_candidate &&
               equations_converged();
    }
    [[nodiscard]] const Sw92PhaseAssignedThreePhaseState* candidate() const & noexcept {
        return candidate_admissible() ? &*point : nullptr;
    }
    const Sw92PhaseAssignedThreePhaseState* candidate() const && = delete;
};

namespace detail {

inline void sw92_phase_assigned_three_phase_check_options(
    const Sw92PhaseAssignedThreePhaseOptions& options) {
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
        options.max_evaluations == 0 || options.max_components == 0 ||
        !std::isfinite(options.new_hydrocarbon_seed_share) ||
        !(options.new_hydrocarbon_seed_share > 0.0) ||
        !(options.new_hydrocarbon_seed_share < 1.0) ||
        options.aqueous_root_options.max_iterations <= 0 ||
        options.nonaqueous_root_options.max_iterations <= 0) {
        throw std::invalid_argument(
            "SW92 phase-assigned C2b.1: invalid tolerance, step or resource option");
    }
}

class Sw92PhaseAssignedThreePhaseStepError : public std::runtime_error {
public:
    Sw92PhaseAssignedThreePhaseStepError(
        Sw92PhaseAssignedThreePhaseStatus status, const char* message)
        : std::runtime_error(message), status_(status) {}
    [[nodiscard]] Sw92PhaseAssignedThreePhaseStatus status() const noexcept {
        return status_;
    }
private:
    Sw92PhaseAssignedThreePhaseStatus status_;
};

enum class Sw92GeneralizedRr3Status {
    success, degenerate, iteration_limit, line_search_failed, unrepresentable
};

struct Sw92GeneralizedRr3Result {
    Sw92GeneralizedRr3Status status{Sw92GeneralizedRr3Status::unrepresentable};
    std::array<double, 2> hydrocarbon_fractions{};
    std::vector<double> aqueous;
    std::vector<double> hydrocarbon0;
    std::vector<double> hydrocarbon1;
    double residual{std::numeric_limits<double>::infinity()};
    double raw_aqueous_sum{};
    double raw_hydrocarbon0_sum{};
    double raw_hydrocarbon1_sum{};
    double mass_absolute{};
    double mass_relative{};
    int iterations{};
};

struct Sw92GeneralizedRr3Evaluation {
    std::array<double, 2> residual{};
    std::array<double, 3> jacobian{};
    std::vector<double> aqueous;
    std::vector<double> hydrocarbon0;
    std::vector<double> hydrocarbon1;
};

inline bool sw92_rr3_fractions_feasible(std::array<double, 2> beta) noexcept {
    return std::isfinite(beta[0]) && std::isfinite(beta[1]) &&
           beta[0] >= 0.0 && beta[1] >= 0.0 &&
           beta[0] + beta[1] <= 1.0;
}

inline std::optional<Sw92GeneralizedRr3Evaluation> sw92_rr3_evaluate(
    std::span<const double> feed,
    std::span<const double> log_k_h0,
    std::span<const double> log_k_h1,
    std::array<double, 2> beta) {
    if (!sw92_rr3_fractions_feasible(beta)) { return std::nullopt; }
    const double beta_w = 1.0 - beta[0] - beta[1];
    Sw92GeneralizedRr3Evaluation value;
    value.aqueous.assign(feed.size(), 0.0);
    value.hydrocarbon0.assign(feed.size(), 0.0);
    value.hydrocarbon1.assign(feed.size(), 0.0);
    double f0 = 0.0, f0_correction = 0.0;
    double f1 = 0.0, f1_correction = 0.0;
    double j00 = 0.0, j00_correction = 0.0;
    double j01 = 0.0, j01_correction = 0.0;
    double j11 = 0.0, j11_correction = 0.0;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0.0) { continue; }
        const double shift = std::max({0.0, log_k_h0[i], log_k_h1[i]});
        const double e_w = std::exp(-shift);
        const double e_0 = std::exp(log_k_h0[i] - shift);
        const double e_1 = std::exp(log_k_h1[i] - shift);
        const double denominator =
            beta_w * e_w + beta[0] * e_0 + beta[1] * e_1;
        if (!std::isfinite(denominator) || !(denominator > 0.0)) {
            return std::nullopt;
        }
        const double x_w = feed[i] * e_w / denominator;
        const double x_0 = feed[i] * e_0 / denominator;
        const double x_1 = feed[i] * e_1 / denominator;
        if (!std::isfinite(x_w) || !std::isfinite(x_0) || !std::isfinite(x_1) ||
            !(x_w > 0.0) || !(x_0 > 0.0) || !(x_1 > 0.0)) {
            return std::nullopt;
        }
        value.aqueous[i] = x_w;
        value.hydrocarbon0[i] = x_0;
        value.hydrocarbon1[i] = x_1;
        const double d0 = x_0 - x_w;
        const double d1 = x_1 - x_w;
        stability_add(d0, f0, f0_correction);
        stability_add(d1, f1, f1_correction);
        stability_add(-(d0 * d0) / feed[i], j00, j00_correction);
        stability_add(-(d0 * d1) / feed[i], j01, j01_correction);
        stability_add(-(d1 * d1) / feed[i], j11, j11_correction);
    }
    value.residual = {f0, f1};
    value.jacobian = {j00, j01, j11};
    return value;
}

inline Sw92GeneralizedRr3Result solve_sw92_generalized_rr3(
    std::span<const double> feed,
    std::span<const double> log_k_h0,
    std::span<const double> log_k_h1,
    std::array<double, 2> initial_fractions,
    const Sw92PhaseAssignedThreePhaseOptions& options) {
    Sw92GeneralizedRr3Result result;
    if (feed.empty() || feed.size() > options.max_components ||
        feed.size() != log_k_h0.size() || feed.size() != log_k_h1.size() ||
        !sw92_rr3_fractions_feasible(initial_fractions)) {
        return result;
    }
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (!std::isfinite(log_k_h0[i]) || !std::isfinite(log_k_h1[i])) {
            return result;
        }
    }
    auto beta = initial_fractions;
    std::optional<Sw92GeneralizedRr3Evaluation> evaluation;
    for (int iteration = 0;; ++iteration) {
        evaluation = sw92_rr3_evaluate(feed, log_k_h0, log_k_h1, beta);
        if (!evaluation) { return result; }
        result.iterations = iteration;
        const double norm = std::max(
            std::abs(evaluation->residual[0]),
            std::abs(evaluation->residual[1]));
        if (!std::isfinite(norm)) { return result; }
        if (norm <= options.balance_tolerance) {
            result.status = Sw92GeneralizedRr3Status::success;
            result.residual = norm;
            break;
        }
        if (iteration >= options.max_balance_iterations) {
            result.status = Sw92GeneralizedRr3Status::iteration_limit;
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
            result.status = Sw92GeneralizedRr3Status::degenerate;
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
            if (!sw92_rr3_fractions_feasible(candidate)) { continue; }
            const auto next = sw92_rr3_evaluate(
                feed, log_k_h0, log_k_h1, candidate);
            if (!next) { continue; }
            const double next_norm = std::max(
                std::abs(next->residual[0]), std::abs(next->residual[1]));
            if (std::isfinite(next_norm) && next_norm < norm) {
                beta = candidate;
                accepted = true;
                break;
            }
        }
        if (!accepted) {
            result.status = Sw92GeneralizedRr3Status::line_search_failed;
            return result;
        }
    }
    result.hydrocarbon_fractions = beta;
    result.aqueous = std::move(evaluation->aqueous);
    result.hydrocarbon0 = std::move(evaluation->hydrocarbon0);
    result.hydrocarbon1 = std::move(evaluation->hydrocarbon1);
    result.raw_aqueous_sum = stability_sum(result.aqueous);
    result.raw_hydrocarbon0_sum = stability_sum(result.hydrocarbon0);
    result.raw_hydrocarbon1_sum = stability_sum(result.hydrocarbon1);
    if (!std::isfinite(result.raw_aqueous_sum) ||
        !std::isfinite(result.raw_hydrocarbon0_sum) ||
        !std::isfinite(result.raw_hydrocarbon1_sum) ||
        !(result.raw_aqueous_sum > 0.0) ||
        !(result.raw_hydrocarbon0_sum > 0.0) ||
        !(result.raw_hydrocarbon1_sum > 0.0)) {
        result.status = Sw92GeneralizedRr3Status::unrepresentable;
        return result;
    }
    const double beta_w = 1.0 - beta[0] - beta[1];
    for (std::size_t i = 0; i < feed.size(); ++i) {
        result.aqueous[i] /= result.raw_aqueous_sum;
        result.hydrocarbon0[i] /= result.raw_hydrocarbon0_sum;
        result.hydrocarbon1[i] /= result.raw_hydrocarbon1_sum;
        const double recovered =
            beta_w * result.aqueous[i] +
            beta[0] * result.hydrocarbon0[i] +
            beta[1] * result.hydrocarbon1[i];
        const double error = std::abs(recovered - feed[i]);
        result.mass_absolute = std::max(result.mass_absolute, error);
        if (feed[i] > 0.0) {
            if (!(result.aqueous[i] > 0.0) ||
                !(result.hydrocarbon0[i] > 0.0) ||
                !(result.hydrocarbon1[i] > 0.0)) {
                result.status = Sw92GeneralizedRr3Status::unrepresentable;
                return result;
            }
            result.mass_relative = std::max(
                result.mass_relative, error / feed[i]);
        } else if (recovered != 0.0) {
            result.status = Sw92GeneralizedRr3Status::unrepresentable;
            return result;
        }
    }
    (void)stability_check_composition(result.aqueous);
    (void)stability_check_composition(result.hydrocarbon0);
    (void)stability_check_composition(result.hydrocarbon1);
    return result;
}

inline double sw92_phase_assigned_log_distance(
    std::span<const double> a, std::span<const double> b,
    std::span<const double> feed) {
    double distance = 0.0;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0.0) { continue; }
        if (!(a[i] > 0.0) || !(b[i] > 0.0)) {
            return std::numeric_limits<double>::infinity();
        }
        distance = std::max(
            distance, std::abs(std::log(a[i]) - std::log(b[i])));
    }
    return distance;
}

inline void sw92_phase_assigned_canonicalize_h_slots(
    Sw92PhaseAssignedThreePhaseState& state,
    std::span<const std::string> component_ids) {
    const double z0 = state.hydrocarbon0_phase.compressibility_factor;
    const double z1 = state.hydrocarbon1_phase.compressibility_factor;
    const double z_guard = 256.0 * stability_eps *
        (1.0 + std::abs(z0) + std::abs(z1));
    bool swap_slots = z0 > z1 + z_guard;
    if (std::abs(z0 - z1) <= z_guard) {
        std::vector<std::size_t> order(component_ids.size());
        std::iota(order.begin(), order.end(), 0U);
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return component_ids[a] < component_ids[b];
        });
        for (std::size_t index : order) {
            const double a = state.hydrocarbon0_phase.composition[index];
            const double b = state.hydrocarbon1_phase.composition[index];
            if (a < b) { break; }
            if (b < a) { swap_slots = true; break; }
        }
    }
    if (!swap_slots) { return; }
    std::swap(state.hydrocarbon0_phase, state.hydrocarbon1_phase);
    std::swap(state.log_k_h0, state.log_k_h1);
    std::swap(state.chemical_potential_residual_h0,
              state.chemical_potential_residual_h1);
    state.hydrocarbon_slots_canonicalized = true;
}

inline bool sw92_phase_assigned_c2b1_source_matches(
    const Sw92PhaseAssignedJointResult& c1,
    const Sw92PhaseAssignedHSideWitnessResult& c2a1) {
    return c1.pressure_pa == c2a1.pressure_pa &&
           c1.temperature_k == c2a1.temperature_k &&
           c1.nacl_molality_mol_per_kg_water ==
               c2a1.nacl_molality_mol_per_kg_water &&
           c1.feed == c2a1.feed &&
           c1.dataset_id == c2a1.dataset_id &&
           c1.revision == c2a1.revision &&
           c1.component_ids == c2a1.component_ids &&
           std::string_view{c2a1.equilibrium_profile} ==
               sw92_phase_assigned_aq_na_joint_profile &&
           std::string_view{c2a1.witness_convention} ==
               sw92_phase_assigned_h_side_na_witness_convention;
}

} // namespace detail

[[nodiscard]] inline Sw92PhaseAssignedThreePhaseResult
iterate_sw92_phase_assigned_three_phase_candidate(
    double pressure_pa, double temperature_k,
    std::span<const double> feed,
    std::span<const double> initial_log_k_h0,
    std::span<const double> initial_log_k_h1,
    std::array<double, 2> initial_hydrocarbon_fractions,
    const thermodynamics::Sw92Phase<double>& model,
    double nacl_molality_mol_per_kg_water,
    Sw92PhaseAssignedThreePhaseOptions options = {}) {
    detail::sw92_phase_assigned_three_phase_check_options(options);
    if (!std::isfinite(pressure_pa) || !(pressure_pa > 0.0) ||
        !std::isfinite(temperature_k) || !(temperature_k > 0.0)) {
        throw std::domain_error(
            "SW92 phase-assigned C2b.1: finite p>0 Pa and T>0 K required");
    }
    if (feed.size() != model.size() ||
        feed.size() != initial_log_k_h0.size() ||
        feed.size() != initial_log_k_h1.size()) {
        throw std::invalid_argument(
            "SW92 phase-assigned C2b.1: ordered snapshot/vector dimension mismatch");
    }
    if (feed.empty() || feed.size() > options.max_components) {
        throw std::length_error(
            "SW92 phase-assigned C2b.1: component quota exceeded or empty feed");
    }
    if (!detail::sw92_rr3_fractions_feasible(initial_hydrocarbon_fractions)) {
        throw std::domain_error(
            "SW92 phase-assigned C2b.1: initial fractions must lie in the closed simplex");
    }
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (!std::isfinite(initial_log_k_h0[i]) ||
            !std::isfinite(initial_log_k_h1[i])) {
            throw std::domain_error(
                "SW92 phase-assigned C2b.1: finite initial logK required");
        }
    }

    Sw92PhaseAssignedThreePhaseResult result;
    result.options = options;
    result.pressure_pa = pressure_pa;
    result.temperature_k = temperature_k;
    result.input_feed_sum = detail::stability_check_composition(feed);
    result.feed = detail::stability_normalize(feed, result.input_feed_sum);
    result.nacl_molality_mol_per_kg_water = nacl_molality_mol_per_kg_water;
    result.initial_log_k_h0.assign(
        initial_log_k_h0.begin(), initial_log_k_h0.end());
    result.initial_log_k_h1.assign(
        initial_log_k_h1.begin(), initial_log_k_h1.end());
    result.initial_hydrocarbon_fractions = initial_hydrocarbon_fractions;
    const auto& parameters = model.parameters();
    result.dataset_id = parameters.dataset_id();
    result.revision = parameters.revision();
    for (const auto& component : parameters.components().items()) {
        result.component_ids.push_back(component.id);
    }
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] > 0.0 && result.feed[i] == 0.0) {
            result.status = Sw92PhaseAssignedThreePhaseStatus::balance_failure;
            result.diagnostic = "roundoff normalization lost an active feed component";
            return result;
        }
    }

    Sw92FamilyStabilityEvaluator aqueous_evaluator(
        model, nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::aqueous,
        options.aqueous_root_options);
    Sw92FamilyStabilityEvaluator nonaqueous_evaluator(
        model, nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::nonaqueous,
        options.nonaqueous_root_options);

    const auto evaluate_selected = [&](Sw92FamilyStabilityEvaluator& evaluator,
                                       std::span<const double> composition) {
        if (result.evaluations >= options.max_evaluations) {
            throw detail::Sw92PhaseAssignedThreePhaseStepError(
                Sw92PhaseAssignedThreePhaseStatus::evaluation_limit,
                "SW92 phase-assigned C2b.1: property budget exhausted");
        }
        ++result.evaluations;
        auto selected = evaluator.evaluate_selected(
            pressure_pa, temperature_k, composition);
        if (!selected.activity.smooth) {
            throw detail::Sw92PhaseAssignedThreePhaseStepError(
                Sw92PhaseAssignedThreePhaseStatus::family_root_nonsmooth,
                "SW92 phase-assigned C2b.1: same-family minimum-root envelope is nonsmooth");
        }
        if (!std::isfinite(selected.compressibility_factor) ||
            !(selected.compressibility_factor > 0.0)) {
            throw StabilityPropertyError(
                StabilityPropertyIssue::nonfinite_properties,
                "SW92 phase-assigned C2b.1: selected Z is nonrepresentable");
        }
        return selected;
    };

    const auto evaluate_state = [&result, &aqueous_evaluator,
                                 &nonaqueous_evaluator, &evaluate_selected, &options](
        std::span<const double> log_k_h0,
        std::span<const double> log_k_h1,
        std::array<double, 2> fraction_seed) {
        const auto balance = detail::solve_sw92_generalized_rr3(
            result.feed, log_k_h0, log_k_h1, fraction_seed, options);
        switch (balance.status) {
        case detail::Sw92GeneralizedRr3Status::success:
            break;
        case detail::Sw92GeneralizedRr3Status::degenerate:
            throw detail::Sw92PhaseAssignedThreePhaseStepError(
                Sw92PhaseAssignedThreePhaseStatus::balance_solver_degenerate,
                "SW92 phase-assigned C2b.1: generalized RR Jacobian is degenerate");
        case detail::Sw92GeneralizedRr3Status::iteration_limit:
            throw detail::Sw92PhaseAssignedThreePhaseStepError(
                Sw92PhaseAssignedThreePhaseStatus::balance_iteration_limit,
                "SW92 phase-assigned C2b.1: generalized RR iteration limit");
        case detail::Sw92GeneralizedRr3Status::line_search_failed:
            throw detail::Sw92PhaseAssignedThreePhaseStepError(
                Sw92PhaseAssignedThreePhaseStatus::balance_line_search_failed,
                "SW92 phase-assigned C2b.1: generalized RR line search failed");
        case detail::Sw92GeneralizedRr3Status::unrepresentable:
            throw detail::Sw92PhaseAssignedThreePhaseStepError(
                Sw92PhaseAssignedThreePhaseStatus::no_resolved_balance_state,
                "SW92 phase-assigned C2b.1: no representable generalized RR state");
        }
        Sw92PhaseAssignedThreePhaseState state;
        state.aqueous_phase.physical_role = Sw92PhysicalPhaseRole::aqueous;
        state.aqueous_phase.thermodynamic_family =
            thermodynamics::SwPhaseFamily::aqueous;
        state.hydrocarbon0_phase.physical_role = Sw92PhysicalPhaseRole::nonaqueous;
        state.hydrocarbon1_phase.physical_role = Sw92PhysicalPhaseRole::nonaqueous;
        state.hydrocarbon0_phase.thermodynamic_family =
            thermodynamics::SwPhaseFamily::nonaqueous;
        state.hydrocarbon1_phase.thermodynamic_family =
            thermodynamics::SwPhaseFamily::nonaqueous;
        state.hydrocarbon0_phase.mole_phase_fraction =
            balance.hydrocarbon_fractions[0];
        state.hydrocarbon1_phase.mole_phase_fraction =
            balance.hydrocarbon_fractions[1];
        state.aqueous_phase.mole_phase_fraction =
            1.0 - balance.hydrocarbon_fractions[0] - balance.hydrocarbon_fractions[1];
        state.aqueous_phase.composition = balance.aqueous;
        state.hydrocarbon0_phase.composition = balance.hydrocarbon0;
        state.hydrocarbon1_phase.composition = balance.hydrocarbon1;
        state.generalized_rr_residual = balance.residual;
        state.raw_aqueous_sum = balance.raw_aqueous_sum;
        state.raw_hydrocarbon0_sum = balance.raw_hydrocarbon0_sum;
        state.raw_hydrocarbon1_sum = balance.raw_hydrocarbon1_sum;
        state.mass_absolute = balance.mass_absolute;
        state.mass_relative = balance.mass_relative;
        state.balance_iterations = balance.iterations;

        const auto selected_w = evaluate_selected(
            aqueous_evaluator, state.aqueous_phase.composition);
        const auto selected_h0 = evaluate_selected(
            nonaqueous_evaluator, state.hydrocarbon0_phase.composition);
        const auto selected_h1 = evaluate_selected(
            nonaqueous_evaluator, state.hydrocarbon1_phase.composition);
        state.aqueous_phase.activity = selected_w.activity;
        state.aqueous_phase.compressibility_factor = selected_w.compressibility_factor;
        state.hydrocarbon0_phase.activity = selected_h0.activity;
        state.hydrocarbon0_phase.compressibility_factor =
            selected_h0.compressibility_factor;
        state.hydrocarbon1_phase.activity = selected_h1.activity;
        state.hydrocarbon1_phase.compressibility_factor =
            selected_h1.compressibility_factor;

        const std::size_t n = result.feed.size();
        state.log_k_h0.assign(n, 0.0);
        state.log_k_h1.assign(n, 0.0);
        state.chemical_potential_residual_h0.assign(n, 0.0);
        state.chemical_potential_residual_h1.assign(n, 0.0);
        state.common_log_activity.assign(n, 0.0);
        double gibbs_correction = 0.0;
        double magnitude = 1.0;
        double magnitude_correction = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            if (result.feed[i] == 0.0) { continue; }
            const double lw = std::log(state.aqueous_phase.composition[i]);
            const double l0 = std::log(state.hydrocarbon0_phase.composition[i]);
            const double l1 = std::log(state.hydrocarbon1_phase.composition[i]);
            const double m_w = lw + state.aqueous_phase.activity.ln_phi[i];
            const double m_0 = l0 + state.hydrocarbon0_phase.activity.ln_phi[i];
            const double m_1 = l1 + state.hydrocarbon1_phase.activity.ln_phi[i];
            state.log_k_h0[i] = l0 - lw;
            state.log_k_h1[i] = l1 - lw;
            state.chemical_potential_residual_h0[i] = m_w - m_0;
            state.chemical_potential_residual_h1[i] = m_w - m_1;
            state.common_log_activity[i] = (m_w + m_0 + m_1) / 3.0;
            state.chemical_potential_norm = std::max({
                state.chemical_potential_norm,
                std::abs(state.chemical_potential_residual_h0[i]),
                std::abs(state.chemical_potential_residual_h1[i])});
            detail::stability_add(
                state.aqueous_phase.mole_phase_fraction *
                    state.aqueous_phase.composition[i] * m_w +
                state.hydrocarbon0_phase.mole_phase_fraction *
                    state.hydrocarbon0_phase.composition[i] * m_0 +
                state.hydrocarbon1_phase.mole_phase_fraction *
                    state.hydrocarbon1_phase.composition[i] * m_1,
                state.reduced_gibbs, gibbs_correction);
            detail::stability_add(
                state.aqueous_phase.mole_phase_fraction *
                    state.aqueous_phase.composition[i] *
                    (std::abs(lw) + std::abs(state.aqueous_phase.activity.ln_phi[i])) +
                state.hydrocarbon0_phase.mole_phase_fraction *
                    state.hydrocarbon0_phase.composition[i] *
                    (std::abs(l0) + std::abs(state.hydrocarbon0_phase.activity.ln_phi[i])) +
                state.hydrocarbon1_phase.mole_phase_fraction *
                    state.hydrocarbon1_phase.composition[i] *
                    (std::abs(l1) + std::abs(state.hydrocarbon1_phase.activity.ln_phi[i])),
                magnitude, magnitude_correction);
        }
        state.gibbs_roundoff_guard = 256.0 * detail::stability_eps * magnitude;
        if (!std::isfinite(state.chemical_potential_norm) ||
            !std::isfinite(state.reduced_gibbs) ||
            !std::isfinite(state.gibbs_roundoff_guard)) {
            throw StabilityPropertyError(
                StabilityPropertyIssue::nonfinite_properties,
                "SW92 phase-assigned C2b.1: nonrepresentable activity/Gibbs arithmetic");
        }
        return state;
    };

    try {
        result.point = evaluate_state(
            result.initial_log_k_h0, result.initial_log_k_h1,
            result.initial_hydrocarbon_fractions);
    } catch (const detail::Sw92PhaseAssignedThreePhaseStepError& error) {
        result.status = error.status();
        result.diagnostic = error.what();
        return result;
    } catch (const StabilityPropertyError& error) {
        result.status = Sw92PhaseAssignedThreePhaseStatus::property_failure;
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
                result.status = Sw92PhaseAssignedThreePhaseStatus::balance_failure;
                result.diagnostic =
                    "SW92 phase-assigned C2b.1: chemical potentials converged but material balance failed";
                return result;
            }
            if (point.aqueous_phase.mole_phase_fraction <= options.minimum_phase_fraction) {
                result.status =
                    Sw92PhaseAssignedThreePhaseStatus::aqueous_phase_disappearance;
                result.diagnostic =
                    "SW92 phase-assigned C2b.1: W(AQ) reached the phase-disappearance boundary; water-topology orchestration remains unresolved";
                return result;
            }
            if (point.hydrocarbon0_phase.mole_phase_fraction <=
                    options.minimum_phase_fraction ||
                point.hydrocarbon1_phase.mole_phase_fraction <=
                    options.minimum_phase_fraction) {
                result.status =
                    Sw92PhaseAssignedThreePhaseStatus::hydrocarbon_phase_disappearance;
                result.diagnostic =
                    "SW92 phase-assigned C2b.1: an unordered H(NA) phase reached the phase-disappearance boundary; no lower topology is accepted here";
                return result;
            }
            point.aqueous_h0_log_distance = detail::sw92_phase_assigned_log_distance(
                point.aqueous_phase.composition,
                point.hydrocarbon0_phase.composition, result.feed);
            point.aqueous_h1_log_distance = detail::sw92_phase_assigned_log_distance(
                point.aqueous_phase.composition,
                point.hydrocarbon1_phase.composition, result.feed);
            point.h0_h1_log_distance = detail::sw92_phase_assigned_log_distance(
                point.hydrocarbon0_phase.composition,
                point.hydrocarbon1_phase.composition, result.feed);
            if (!(point.aqueous_h0_log_distance > options.log_composition_separation) ||
                !(point.aqueous_h1_log_distance > options.log_composition_separation) ||
                !(point.h0_h1_log_distance > options.log_composition_separation)) {
                result.status = Sw92PhaseAssignedThreePhaseStatus::indistinguishable_phases;
                result.diagnostic =
                    "SW92 phase-assigned C2b.1: converged phase compositions are not pairwise distinct";
                return result;
            }
            const std::size_t water_index = parameters.water_index();
            const double w_water = point.aqueous_phase.composition[water_index];
            const double h0_water = point.hydrocarbon0_phase.composition[water_index];
            const double h1_water = point.hydrocarbon1_phase.composition[water_index];
            point.aqueous_minus_h0_water_fraction = w_water - h0_water;
            point.aqueous_minus_h1_water_fraction = w_water - h1_water;
            point.water_role_roundoff_guard = 256.0 * detail::stability_eps *
                (1.0 + std::abs(w_water) + std::abs(h0_water) + std::abs(h1_water));
            if (!std::isfinite(point.aqueous_minus_h0_water_fraction) ||
                !std::isfinite(point.aqueous_minus_h1_water_fraction) ||
                !std::isfinite(point.water_role_roundoff_guard)) {
                result.status = Sw92PhaseAssignedThreePhaseStatus::phase_role_indeterminate;
                result.diagnostic =
                    "SW92 phase-assigned C2b.1: W/H water-richness ordering is nonrepresentable";
                return result;
            }
            if (std::abs(point.aqueous_minus_h0_water_fraction) <=
                    point.water_role_roundoff_guard ||
                std::abs(point.aqueous_minus_h1_water_fraction) <=
                    point.water_role_roundoff_guard) {
                result.status = Sw92PhaseAssignedThreePhaseStatus::phase_role_indeterminate;
                result.diagnostic =
                    "SW92 phase-assigned C2b.1: W/H water-richness ordering is unresolved at roundoff scale";
                return result;
            }
            if (point.aqueous_minus_h0_water_fraction < 0.0 ||
                point.aqueous_minus_h1_water_fraction < 0.0) {
                result.status = Sw92PhaseAssignedThreePhaseStatus::phase_role_reversed;
                result.diagnostic =
                    "SW92 phase-assigned C2b.1: assigned W(AQ) is not water-richer than both H(NA) phases";
                return result;
            }
            detail::sw92_phase_assigned_canonicalize_h_slots(
                point, result.component_ids);
            point.aqueous_h0_log_distance = detail::sw92_phase_assigned_log_distance(
                point.aqueous_phase.composition,
                point.hydrocarbon0_phase.composition, result.feed);
            point.aqueous_h1_log_distance = detail::sw92_phase_assigned_log_distance(
                point.aqueous_phase.composition,
                point.hydrocarbon1_phase.composition, result.feed);
            point.h0_h1_log_distance = detail::sw92_phase_assigned_log_distance(
                point.hydrocarbon0_phase.composition,
                point.hydrocarbon1_phase.composition, result.feed);
            point.aqueous_minus_h0_water_fraction =
                point.aqueous_phase.composition[water_index] -
                point.hydrocarbon0_phase.composition[water_index];
            point.aqueous_minus_h1_water_fraction =
                point.aqueous_phase.composition[water_index] -
                point.hydrocarbon1_phase.composition[water_index];
            result.status = Sw92PhaseAssignedThreePhaseStatus::converged_candidate;
            result.diagnostic =
                "unordered W(AQ)+H0(NA)+H1(NA) equations, three-phase material balance and relative-water topology converged; H morphology and final topology/stability publication are intentionally not resolved";
            return result;
        }
        if (result.iterations >= options.max_iterations) {
            result.status = Sw92PhaseAssignedThreePhaseStatus::iteration_limit;
            result.diagnostic = "SW92 phase-assigned C2b.1: iteration limit";
            return result;
        }
        const double scale = std::max(
            1.0, point.chemical_potential_norm / options.max_log_step);
        bool accepted = false;
        double alpha = 1.0;
        for (int backtrack = 0; backtrack < options.max_backtracks;
             ++backtrack, alpha *= 0.5) {
            if (backtrack > 0) { ++result.backtracks; }
            std::vector<double> next_log_k_h0(result.feed.size(), 0.0);
            std::vector<double> next_log_k_h1(result.feed.size(), 0.0);
            for (std::size_t i = 0; i < result.feed.size(); ++i) {
                next_log_k_h0[i] = point.log_k_h0[i] +
                    (alpha / scale) * point.chemical_potential_residual_h0[i];
                next_log_k_h1[i] = point.log_k_h1[i] +
                    (alpha / scale) * point.chemical_potential_residual_h1[i];
            }
            const std::array<double, 2> fraction_seed{
                point.hydrocarbon0_phase.mole_phase_fraction,
                point.hydrocarbon1_phase.mole_phase_fraction};
            try {
                auto next = evaluate_state(
                    next_log_k_h0, next_log_k_h1, fraction_seed);
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
            } catch (const detail::Sw92PhaseAssignedThreePhaseStepError& error) {
                ++result.rejected_evaluations;
                result.diagnostic = error.what();
                if (error.status() == Sw92PhaseAssignedThreePhaseStatus::evaluation_limit) {
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
            result.status = Sw92PhaseAssignedThreePhaseStatus::line_search_failed;
            if (result.diagnostic.empty()) {
                result.diagnostic =
                    "SW92 phase-assigned C2b.1: no acceptable Gibbs/residual line-search step";
            }
            return result;
        }
        ++result.iterations;
    }
}

[[nodiscard]] inline Sw92PhaseAssignedThreePhaseResult
solve_sw92_phase_assigned_c2b1_candidate(
    const Sw92PhaseAssignedJointResult& c1,
    const Sw92PhaseAssignedHSideWitnessResult& c2a1,
    std::size_t witness_index,
    const thermodynamics::Sw92Phase<double>& model,
    Sw92PhaseAssignedThreePhaseOptions options = {}) {
    detail::sw92_phase_assigned_three_phase_check_options(options);
    if (!detail::sw92_phase_assigned_h_side_model_matches(c1, model)) {
        throw std::invalid_argument(
            "SW92 phase-assigned C2b.1: C1 result/model snapshot mismatch");
    }
    Sw92PhaseAssignedThreePhaseResult rejected;
    rejected.options = options;
    rejected.pressure_pa = c1.pressure_pa;
    rejected.temperature_k = c1.temperature_k;
    rejected.feed = c1.feed;
    rejected.nacl_molality_mol_per_kg_water = c1.nacl_molality_mol_per_kg_water;
    rejected.dataset_id = c1.dataset_id;
    rejected.revision = c1.revision;
    rejected.component_ids = c1.component_ids;
    rejected.source_witness_index = witness_index;
    if (c1.candidate() == nullptr) {
        rejected.status =
            Sw92PhaseAssignedThreePhaseStatus::source_candidate_unavailable;
        rejected.diagnostic =
            "Profile-C C2b.1 requires an admissible C1 W(AQ)+H(NA) candidate";
        return rejected;
    }
    const auto c1_evidence =
        detail::sw92_phase_assigned_revalidate_c1_for_h_side_witness(c1, model);
    if (!c1_evidence) {
        rejected.status =
            Sw92PhaseAssignedThreePhaseStatus::source_candidate_inconsistent;
        rejected.diagnostic =
            "retained C1 state failed independent equilibrium/material-balance/role revalidation";
        return rejected;
    }
    if (!detail::sw92_phase_assigned_c2b1_source_matches(c1, c2a1) ||
        c2a1.status != Sw92PhaseAssignedHSideWitnessStatus::
                           additional_nonaqueous_phase_witness_found ||
        !c2a1.nonaqueous_search) {
        rejected.status =
            Sw92PhaseAssignedThreePhaseStatus::source_witness_unavailable;
        rejected.diagnostic =
            "Profile-C C2b.1 requires a matching C2a1 additional-NA witness result";
        return rejected;
    }
    if (witness_index >= c2a1.negative_witnesses.size()) {
        rejected.status =
            Sw92PhaseAssignedThreePhaseStatus::source_witness_unavailable;
        rejected.diagnostic = "C2b.1 witness index is out of range";
        return rejected;
    }
    const auto& witness = c2a1.negative_witnesses[witness_index];
    const auto* source_point = c1.candidate();
    const auto reclassified = detail::sw92_phase_assigned_classify_na_witness(
        witness.trial_index, witness.point,
        source_point->nonaqueous_phase.composition,
        c1_evidence->retained_w_water_fraction,
        model.parameters().water_index(), c2a1.options.log_composition_separation);
    if (!detail::stability_negative(
            witness.point, c2a1.nonaqueous_search->options) ||
        !reclassified.usable_h_split_seed()) {
        rejected.status =
            Sw92PhaseAssignedThreePhaseStatus::source_witness_inconsistent;
        rejected.diagnostic =
            "selected C2a1 trial is not a robust, distinct and role-admissible NA witness under retained options";
        return rejected;
    }
    const std::size_t n = c1.feed.size();
    std::vector<double> log_k_h0(n, 0.0);
    std::vector<double> log_k_h1(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (c1.feed[i] == 0.0) { continue; }
        const double x_w = source_point->aqueous_phase.composition[i];
        const double x_h0 = source_point->nonaqueous_phase.composition[i];
        const double x_h1 = witness.point.composition[i];
        if (!(x_w > 0.0) || !(x_h0 > 0.0) || !(x_h1 > 0.0)) {
            rejected.status =
                Sw92PhaseAssignedThreePhaseStatus::source_witness_inconsistent;
            rejected.diagnostic =
                "C2b.1 active source compositions must remain strictly positive";
            return rejected;
        }
        log_k_h0[i] = std::log(x_h0) - std::log(x_w);
        log_k_h1[i] = std::log(x_h1) - std::log(x_w);
    }
    const double retained_h_fraction =
        source_point->nonaqueous_phase.mole_phase_fraction;
    const std::array<double, 2> fractions{
        retained_h_fraction * (1.0 - options.new_hydrocarbon_seed_share),
        retained_h_fraction * options.new_hydrocarbon_seed_share};
    auto result = iterate_sw92_phase_assigned_three_phase_candidate(
        c1.pressure_pa, c1.temperature_k, c1.feed,
        log_k_h0, log_k_h1, fractions, model,
        c1.nacl_molality_mol_per_kg_water, options);
    result.source_witness_index = witness_index;
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_PHASE_ASSIGNED_THREE_PHASE_HPP
