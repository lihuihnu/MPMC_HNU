#ifndef MPMC_FLASH_SW92_ASYMMETRIC_FIXED_PAIR_HPP
#define MPMC_FLASH_SW92_ASYMMETRIC_FIXED_PAIR_HPP

#include <mpmc/flash/rachford_rice.hpp>
#include <mpmc/flash/sw92_asymmetric_stability.hpp>

#include <algorithm>
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

/// Gate 3B.1 numerical primitive identity. This solves one fixed family pair;
/// it is not the future maximum-two-phase orchestration/final-acceptance profile.
inline constexpr std::string_view sw92_xu_asymmetric_fixed_pair_convention =
    "SW92-equilibrium/xu-asymmetric-gibbs/fixed-pair-logK-SSI-RR/v1";

struct Sw92AsymmetricFamilyPair {
    thermodynamics::SwPhaseFamily phase0{thermodynamics::SwPhaseFamily::aqueous};
    thermodynamics::SwPhaseFamily phase1{thermodynamics::SwPhaseFamily::nonaqueous};
};

struct Sw92AsymmetricFixedPairOptions {
    double chemical_potential_tolerance{1e-11}; // max |Delta[ln(x_i*phi_i)]|.
    double mass_absolute_tolerance{1e-12};
    double mass_relative_tolerance{1e-10};      // Positive feed components only.
    double minimum_phase_fraction{1e-10};       // Numerical disappearance guard.
    double log_k_separation{1e-7};              // Composition distinction, not a phase label.
    double max_log_step{2.0};
    double residual_decrease{1e-4};
    double gibbs_armijo{1e-4};
    int max_iterations{512};
    int max_backtracks{24};
    std::size_t max_evaluations{20000};         // Same-family minimum-root provider calls.
    RachfordRiceOptions rr;
    thermodynamics::Sw92RootOptions aqueous_root_options;
    thermodynamics::Sw92RootOptions nonaqueous_root_options;
};

enum class Sw92AsymmetricFamilyAssignmentStatus {
    not_checked,
    assigned_lower,
    assigned_dominated,
    family_tie,
    opposite_family_nonsmooth,
    property_failure
};

struct Sw92AsymmetricFamilyAssignmentCheck {
    Sw92AsymmetricFamilyAssignmentStatus status{
        Sw92AsymmetricFamilyAssignmentStatus::not_checked};
    thermodynamics::SwPhaseFamily assigned_family{thermodynamics::SwPhaseFamily::aqueous};
    double aqueous_minus_nonaqueous_gibbs{
        std::numeric_limits<double>::quiet_NaN()};
    double roundoff_guard{std::numeric_limits<double>::quiet_NaN()};
    std::optional<StabilityPropertyIssue> property_issue;
};

/// One slot in a fixed-pair equation candidate. Slot order is numerical only;
/// `family` is the thermodynamic model identity and Z is only a diagnostic.
struct Sw92AsymmetricFixedPairPhase {
    thermodynamics::SwPhaseFamily family{thermodynamics::SwPhaseFamily::aqueous};
    double mole_phase_fraction{};
    std::vector<double> composition;
    StabilityPhase activity;
    double compressibility_factor{};
};

struct Sw92AsymmetricFixedPairState {
    Sw92AsymmetricFixedPairPhase phase0;
    Sw92AsymmetricFixedPairPhase phase1;

    std::vector<double> log_k; // ln(x_phase1/x_phase0) on active support.
    std::vector<double> chemical_potential_residual; // m0-m1, dimensionless.
    std::vector<double> common_log_activity;          // midpoint(m0,m1).
    double chemical_potential_norm{};
    double reduced_gibbs{};
    double gibbs_roundoff_guard{};

    double rr_residual{};
    double raw_phase0_sum{};
    double raw_phase1_sum{};
    double mass_absolute{};
    double mass_relative{};
    int rr_iterations{};

    Sw92AsymmetricFamilyAssignmentCheck phase0_assignment;
    Sw92AsymmetricFamilyAssignmentCheck phase1_assignment;
};

enum class Sw92AsymmetricFixedPairStatus {
    converged,
    no_interior_rr_root,
    rr_failure,
    indistinguishable_phases,
    phase_disappearance,
    balance_failure,
    family_root_nonsmooth,
    family_assignment_dominated,
    family_assignment_nonsmooth,
    family_assignment_property_failure,
    property_failure,
    iteration_limit,
    evaluation_limit,
    line_search_failed
};

/// Gate 3B.1 fixed-family-pair result. `converged` additionally means that the
/// equation-converged point passed phase distinction, root smoothness and both
/// lower-envelope family-assignment checks. It is NOT an accepted overall state.
struct Sw92AsymmetricFixedPairResult {
    Sw92AsymmetricFixedPairStatus status{
        Sw92AsymmetricFixedPairStatus::rr_failure};
    Sw92AsymmetricFamilyPair family_pair;
    Sw92AsymmetricFixedPairOptions options;

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
    std::string equilibrium_profile{sw92_xu_asymmetric_gibbs_profile};
    std::string fixed_pair_convention{sw92_xu_asymmetric_fixed_pair_convention};

    std::vector<double> initial_log_k;
    std::optional<Sw92AsymmetricFixedPairState> point; // Last accepted iterate.
    int iterations{};
    std::size_t evaluations{};
    std::size_t backtracks{};
    std::size_t rejected_evaluations{};
    std::size_t gibbs_descent_steps{};
    std::size_t residual_increase_steps{};
    std::optional<StabilityPropertyIssue> property_issue;
    std::string diagnostic;

    /// Equation convergence is independent of the later family-assignment gate.
    /// Dominated/tied pairs intentionally retain equation-converged points.
    [[nodiscard]] bool equations_converged() const noexcept {
        return point.has_value() &&
               point->chemical_potential_norm <= options.chemical_potential_tolerance &&
               point->mass_absolute <= options.mass_absolute_tolerance &&
               point->mass_relative <= options.mass_relative_tolerance;
    }

    /// Admissible only within this one fixed-family-pair primitive. Overall
    /// Gate-3B candidate selection/final phase-set stability are still absent.
    [[nodiscard]] bool candidate_admissible() const noexcept {
        return status == Sw92AsymmetricFixedPairStatus::converged && equations_converged();
    }
};

namespace detail {

inline void sw92_fixed_pair_check_family(thermodynamics::SwPhaseFamily family) {
    switch (family) {
    case thermodynamics::SwPhaseFamily::aqueous:
    case thermodynamics::SwPhaseFamily::nonaqueous:
        return;
    }
    throw std::invalid_argument("SW92 asymmetric fixed pair: unknown phase family");
}

inline void sw92_fixed_pair_check_options(const Sw92AsymmetricFixedPairOptions& options) {
    const auto positive = [](double value) {
        return std::isfinite(value) && value > 0.0;
    };
    if (!positive(options.chemical_potential_tolerance) ||
        !positive(options.mass_absolute_tolerance) ||
        !positive(options.mass_relative_tolerance) ||
        !positive(options.minimum_phase_fraction) ||
        options.minimum_phase_fraction >= 0.5 ||
        !positive(options.log_k_separation) ||
        !positive(options.max_log_step) ||
        !positive(options.residual_decrease) || options.residual_decrease >= 1.0 ||
        !positive(options.gibbs_armijo) || options.gibbs_armijo >= 1.0 ||
        options.max_iterations < 0 || options.max_backtracks <= 0 ||
        options.max_evaluations == 0 ||
        options.aqueous_root_options.max_iterations <= 0 ||
        options.nonaqueous_root_options.max_iterations <= 0) {
        throw std::invalid_argument(
            "SW92 asymmetric fixed pair: invalid tolerance, step or resource option");
    }
    rr_check_options(options.rr);
}

class Sw92FixedPairStepError : public std::runtime_error {
public:
    Sw92FixedPairStepError(Sw92AsymmetricFixedPairStatus status, const char* message)
        : std::runtime_error(message), status_(status) {}
    [[nodiscard]] Sw92AsymmetricFixedPairStatus status() const noexcept { return status_; }
private:
    Sw92AsymmetricFixedPairStatus status_;
};

inline bool sw92_fixed_pair_balance_ok(
    const Sw92AsymmetricFixedPairState& state,
    const Sw92AsymmetricFixedPairOptions& options) noexcept {
    return state.mass_absolute <= options.mass_absolute_tolerance &&
           state.mass_relative <= options.mass_relative_tolerance;
}

inline double sw92_fixed_pair_gibbs_slope(
    const Sw92AsymmetricFixedPairState& point,
    std::span<const double> feed, double scale) {
    double h = 0.0;
    double h_correction = 0.0;
    double first = 0.0;
    double first_correction = 0.0;
    double second = 0.0;
    double second_correction = 0.0;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0.0) { continue; }
        const double x = point.phase0.composition[i];
        const double y = point.phase1.composition[i];
        const double c = (x / feed[i]) * y;
        const double difference = y - x;
        const double residual = point.chemical_potential_residual[i];
        stability_add(difference * (difference / feed[i]), h, h_correction);
        stability_add(c * residual, first, first_correction);
        stability_add((c * residual) * residual, second, second_correction);
    }
    if (!(h > 0.0) || !std::isfinite(h) || !std::isfinite(first) ||
        !std::isfinite(second) || !std::isfinite(scale) || !(scale >= 1.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double beta = point.phase1.mole_phase_fraction;
    return -(beta * (1.0 - beta) * second + (first / h) * first) / scale;
}

inline bool sw92_fixed_pair_accept_step(
    const Sw92AsymmetricFixedPairState& point,
    const Sw92AsymmetricFixedPairState& next,
    double predicted, double relative_step,
    const Sw92AsymmetricFixedPairOptions& options) noexcept {
    const double guard = point.gibbs_roundoff_guard + next.gibbs_roundoff_guard;
    if (predicted > guard) {
        return next.reduced_gibbs <= point.reduced_gibbs - options.gibbs_armijo * predicted;
    }
    return next.reduced_gibbs <= point.reduced_gibbs + guard &&
           (next.chemical_potential_norm <= options.chemical_potential_tolerance ||
            point.chemical_potential_norm - next.chemical_potential_norm >
                options.residual_decrease * relative_step * point.chemical_potential_norm);
}

inline Sw92AsymmetricFixedPairStatus sw92_fixed_pair_rr_status(
    RachfordRiceStatus status) noexcept {
    switch (status) {
    case RachfordRiceStatus::interior:
        return Sw92AsymmetricFixedPairStatus::converged; // Sentinel; caller continues.
    case RachfordRiceStatus::no_resolved_interior_root:
        return Sw92AsymmetricFixedPairStatus::no_interior_rr_root;
    case RachfordRiceStatus::degenerate:
        return Sw92AsymmetricFixedPairStatus::indistinguishable_phases;
    case RachfordRiceStatus::iteration_limit:
    case RachfordRiceStatus::unrepresentable:
        return Sw92AsymmetricFixedPairStatus::rr_failure;
    }
    return Sw92AsymmetricFixedPairStatus::rr_failure;
}

} // namespace detail

/// Solve one explicit ordered family pair. Repeated families are legal. The
/// primitive performs no Gate-3A witness planning, no candidate deduplication,
/// no pair-vs-feed selection and no final two-family phase-set stability review.
[[nodiscard]] inline Sw92AsymmetricFixedPairResult iterate_sw92_asymmetric_fixed_pair(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    std::span<const double> initial_log_k,
    Sw92AsymmetricFamilyPair family_pair,
    const thermodynamics::Sw92Phase<double>& model,
    double nacl_molality_mol_per_kg_water,
    Sw92AsymmetricFixedPairOptions options = {}) {
    detail::sw92_fixed_pair_check_options(options);
    detail::sw92_fixed_pair_check_family(family_pair.phase0);
    detail::sw92_fixed_pair_check_family(family_pair.phase1);
    if (!std::isfinite(pressure_pa) || pressure_pa <= 0.0 ||
        !std::isfinite(temperature_k) || temperature_k <= 0.0) {
        throw std::domain_error(
            "SW92 asymmetric fixed pair: finite p>0 Pa and T>0 K required");
    }
    if (feed.size() != model.size()) {
        throw std::invalid_argument(
            "SW92 asymmetric fixed pair: feed does not match ordered model snapshot");
    }
    if (feed.size() != initial_log_k.size()) {
        throw std::invalid_argument(
            "SW92 asymmetric fixed pair: logK dimension mismatch");
    }
    if (feed.empty() || feed.size() > options.rr.max_components) {
        throw std::length_error(
            "SW92 asymmetric fixed pair: component quota exceeded or empty feed");
    }
    for (double value : initial_log_k) {
        if (!std::isfinite(value)) {
            throw std::domain_error(
                "SW92 asymmetric fixed pair: finite initial logK required");
        }
    }

    Sw92AsymmetricFixedPairResult result;
    result.family_pair = family_pair;
    result.options = options;
    result.pressure_pa = pressure_pa;
    result.temperature_k = temperature_k;
    result.input_feed_sum = detail::stability_check_composition(feed);
    result.feed = detail::stability_normalize(feed, result.input_feed_sum);
    result.nacl_molality_mol_per_kg_water = nacl_molality_mol_per_kg_water;
    result.initial_log_k.assign(initial_log_k.begin(), initial_log_k.end());

    const auto& parameters = model.parameters();
    result.dataset_id = parameters.dataset_id();
    result.revision = parameters.revision();
    for (const auto& component : parameters.components().items()) {
        result.component_ids.push_back(component.id);
    }

    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] > 0.0 && result.feed[i] == 0.0) {
            result.status = Sw92AsymmetricFixedPairStatus::rr_failure;
            result.diagnostic = "roundoff normalization lost an active feed component";
            return result;
        }
    }

    Sw92FamilyStabilityEvaluator aqueous_evaluator(
        model, nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::aqueous, options.aqueous_root_options);
    Sw92FamilyStabilityEvaluator nonaqueous_evaluator(
        model, nacl_molality_mol_per_kg_water,
        thermodynamics::SwPhaseFamily::nonaqueous, options.nonaqueous_root_options);

    const auto evaluator_for = [&](thermodynamics::SwPhaseFamily family)
        -> Sw92FamilyStabilityEvaluator& {
        return family == thermodynamics::SwPhaseFamily::aqueous
            ? aqueous_evaluator : nonaqueous_evaluator;
    };

    const auto evaluate_selected = [&](thermodynamics::SwPhaseFamily family,
                                       std::span<const double> composition) {
        if (result.evaluations >= options.max_evaluations) {
            throw detail::Sw92FixedPairStepError(
                Sw92AsymmetricFixedPairStatus::evaluation_limit,
                "SW92 asymmetric fixed pair: property budget exhausted");
        }
        ++result.evaluations;
        auto selected = evaluator_for(family).evaluate_selected(
            pressure_pa, temperature_k, composition);
        if (!selected.activity.smooth) {
            throw detail::Sw92FixedPairStepError(
                Sw92AsymmetricFixedPairStatus::family_root_nonsmooth,
                "SW92 asymmetric fixed pair: same-family minimum-root envelope is nonsmooth");
        }
        if (!std::isfinite(selected.compressibility_factor) ||
            !(selected.compressibility_factor > 0.0)) {
            throw StabilityPropertyError(
                StabilityPropertyIssue::nonfinite_properties,
                "SW92 asymmetric fixed pair: selected Z is nonrepresentable");
        }
        return selected;
    };

    const auto evaluate_state = [&](std::span<const double> log_k) {
        const auto rr = solve_rachford_rice(result.feed, log_k, options.rr);
        if (rr.status != RachfordRiceStatus::interior) {
            throw detail::Sw92FixedPairStepError(
                detail::sw92_fixed_pair_rr_status(rr.status),
                "SW92 asymmetric fixed pair: Rachford-Rice has no usable interior state");
        }

        Sw92AsymmetricFixedPairState state;
        state.phase0.family = family_pair.phase0;
        state.phase1.family = family_pair.phase1;
        state.phase1.mole_phase_fraction = rr.vapor_fraction;
        state.phase0.mole_phase_fraction = 1.0 - rr.vapor_fraction;
        state.phase0.composition = rr.liquid;
        state.phase1.composition = rr.vapor;
        state.rr_residual = rr.residual;
        state.raw_phase0_sum = rr.raw_liquid_sum;
        state.raw_phase1_sum = rr.raw_vapor_sum;
        state.mass_absolute = rr.mass_absolute;
        state.mass_relative = rr.mass_relative;
        state.rr_iterations = rr.iterations;

        const auto selected0 = evaluate_selected(
            family_pair.phase0, state.phase0.composition);
        const auto selected1 = evaluate_selected(
            family_pair.phase1, state.phase1.composition);
        state.phase0.activity = selected0.activity;
        state.phase0.compressibility_factor = selected0.compressibility_factor;
        state.phase1.activity = selected1.activity;
        state.phase1.compressibility_factor = selected1.compressibility_factor;

        const std::size_t n = result.feed.size();
        state.log_k.resize(n);
        state.chemical_potential_residual.resize(n);
        state.common_log_activity.resize(n);
        double gibbs_correction = 0.0;
        double magnitude = 1.0;
        double magnitude_correction = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            if (result.feed[i] == 0.0) { continue; }
            const double lx = std::log(state.phase0.composition[i]);
            const double ly = std::log(state.phase1.composition[i]);
            const double m0 = lx + state.phase0.activity.ln_phi[i];
            const double m1 = ly + state.phase1.activity.ln_phi[i];
            state.log_k[i] = ly - lx;
            state.chemical_potential_residual[i] = m0 - m1;
            state.common_log_activity[i] = std::midpoint(m0, m1);
            state.chemical_potential_norm = std::max(
                state.chemical_potential_norm,
                std::abs(state.chemical_potential_residual[i]));
            detail::stability_add(
                state.phase0.mole_phase_fraction * state.phase0.composition[i] * m0 +
                state.phase1.mole_phase_fraction * state.phase1.composition[i] * m1,
                state.reduced_gibbs, gibbs_correction);
            detail::stability_add(
                state.phase0.mole_phase_fraction * state.phase0.composition[i] *
                    (std::abs(lx) + std::abs(state.phase0.activity.ln_phi[i])) +
                state.phase1.mole_phase_fraction * state.phase1.composition[i] *
                    (std::abs(ly) + std::abs(state.phase1.activity.ln_phi[i])),
                magnitude, magnitude_correction);
            if (!std::isfinite(state.log_k[i]) ||
                !std::isfinite(state.chemical_potential_residual[i]) ||
                !std::isfinite(state.common_log_activity[i])) {
                throw StabilityPropertyError(
                    StabilityPropertyIssue::nonfinite_properties,
                    "SW92 asymmetric fixed pair: nonrepresentable activity arithmetic");
            }
        }
        state.gibbs_roundoff_guard = 256.0 * detail::stability_eps * magnitude;
        if (!std::isfinite(state.reduced_gibbs) ||
            !std::isfinite(state.gibbs_roundoff_guard)) {
            throw StabilityPropertyError(
                StabilityPropertyIssue::nonfinite_properties,
                "SW92 asymmetric fixed pair: nonrepresentable Gibbs arithmetic");
        }
        return state;
    };

    const auto check_assignment = [&](Sw92AsymmetricFixedPairPhase& phase) {
        Sw92AsymmetricFamilyAssignmentCheck check;
        check.assigned_family = phase.family;
        const auto opposite = phase.family == thermodynamics::SwPhaseFamily::aqueous
            ? thermodynamics::SwPhaseFamily::nonaqueous
            : thermodynamics::SwPhaseFamily::aqueous;
        try {
            if (result.evaluations >= options.max_evaluations) {
                throw detail::Sw92FixedPairStepError(
                    Sw92AsymmetricFixedPairStatus::evaluation_limit,
                    "SW92 asymmetric fixed pair: property budget exhausted during family check");
            }
            ++result.evaluations;
            auto opposite_phase = evaluator_for(opposite).evaluate_selected(
                pressure_pa, temperature_k, phase.composition);
            if (!opposite_phase.activity.smooth) {
                check.status =
                    Sw92AsymmetricFamilyAssignmentStatus::opposite_family_nonsmooth;
                return check;
            }

            const StabilityPhase& aqueous =
                phase.family == thermodynamics::SwPhaseFamily::aqueous
                    ? phase.activity : opposite_phase.activity;
            const StabilityPhase& nonaqueous =
                phase.family == thermodynamics::SwPhaseFamily::nonaqueous
                    ? phase.activity : opposite_phase.activity;
            const auto [difference, guard] = detail::sw92_feed_family_difference(
                phase.composition, aqueous, nonaqueous);
            check.aqueous_minus_nonaqueous_gibbs = difference;
            check.roundoff_guard = guard;
            if (std::abs(difference) <= guard) {
                check.status = Sw92AsymmetricFamilyAssignmentStatus::family_tie;
                return check;
            }
            const bool aqueous_lower = difference < 0.0;
            const bool assigned_lower =
                (phase.family == thermodynamics::SwPhaseFamily::aqueous && aqueous_lower) ||
                (phase.family == thermodynamics::SwPhaseFamily::nonaqueous && !aqueous_lower);
            check.status = assigned_lower
                ? Sw92AsymmetricFamilyAssignmentStatus::assigned_lower
                : Sw92AsymmetricFamilyAssignmentStatus::assigned_dominated;
            return check;
        } catch (const StabilityPropertyError& error) {
            check.status = Sw92AsymmetricFamilyAssignmentStatus::property_failure;
            check.property_issue = error.issue();
            return check;
        }
    };

    try {
        result.point = evaluate_state(result.initial_log_k);
    } catch (const detail::Sw92FixedPairStepError& error) {
        result.status = error.status();
        result.diagnostic = error.what();
        return result;
    } catch (const StabilityPropertyError& error) {
        result.status = Sw92AsymmetricFixedPairStatus::property_failure;
        result.property_issue = error.issue();
        result.diagnostic = error.what();
        return result;
    }

    for (;;) {
        auto& point = *result.point;
        if (point.chemical_potential_norm <= options.chemical_potential_tolerance) {
            if (!detail::sw92_fixed_pair_balance_ok(point, options)) {
                result.status = Sw92AsymmetricFixedPairStatus::balance_failure;
                result.diagnostic =
                    "SW92 asymmetric fixed pair: chemical potentials converged but material balance failed";
                return result;
            }
            if (point.phase0.mole_phase_fraction <= options.minimum_phase_fraction ||
                point.phase1.mole_phase_fraction <= options.minimum_phase_fraction) {
                result.status = Sw92AsymmetricFixedPairStatus::phase_disappearance;
                result.diagnostic =
                    "SW92 asymmetric fixed pair: balance and chemical-potential tolerances passed at the phase-disappearance threshold";
                return result;
            }
            double contrast = 0.0;
            for (std::size_t i = 0; i < result.feed.size(); ++i) {
                if (result.feed[i] > 0.0) {
                    contrast = std::max(contrast, std::abs(point.log_k[i]));
                }
            }
            if (contrast <= options.log_k_separation) {
                result.status = Sw92AsymmetricFixedPairStatus::indistinguishable_phases;
                result.diagnostic =
                    "SW92 asymmetric fixed pair: converged phase compositions are numerically indistinguishable";
                return result;
            }

            try {
                point.phase0_assignment = check_assignment(point.phase0);
                point.phase1_assignment = check_assignment(point.phase1);
            } catch (const detail::Sw92FixedPairStepError& error) {
                result.status = error.status();
                result.diagnostic = error.what();
                return result;
            }

            const auto assignment_failure = [](const Sw92AsymmetricFamilyAssignmentCheck& check) {
                return check.status == Sw92AsymmetricFamilyAssignmentStatus::property_failure;
            };
            const auto assignment_nonsmooth = [](const Sw92AsymmetricFamilyAssignmentCheck& check) {
                return check.status == Sw92AsymmetricFamilyAssignmentStatus::family_tie ||
                       check.status ==
                           Sw92AsymmetricFamilyAssignmentStatus::opposite_family_nonsmooth;
            };
            const auto assignment_dominated = [](const Sw92AsymmetricFamilyAssignmentCheck& check) {
                return check.status ==
                    Sw92AsymmetricFamilyAssignmentStatus::assigned_dominated;
            };

            if (assignment_failure(point.phase0_assignment) ||
                assignment_failure(point.phase1_assignment)) {
                result.status =
                    Sw92AsymmetricFixedPairStatus::family_assignment_property_failure;
                result.property_issue = point.phase0_assignment.property_issue
                    ? point.phase0_assignment.property_issue
                    : point.phase1_assignment.property_issue;
                result.diagnostic =
                    "SW92 asymmetric fixed pair: lower-envelope family assignment could not be evaluated";
                return result;
            }
            if (assignment_nonsmooth(point.phase0_assignment) ||
                assignment_nonsmooth(point.phase1_assignment)) {
                result.status = Sw92AsymmetricFixedPairStatus::family_assignment_nonsmooth;
                result.diagnostic =
                    "SW92 asymmetric fixed pair: AQ/NA lower-envelope family identity is nonsmooth at a converged phase";
                return result;
            }
            if (assignment_dominated(point.phase0_assignment) ||
                assignment_dominated(point.phase1_assignment)) {
                result.status = Sw92AsymmetricFixedPairStatus::family_assignment_dominated;
                result.diagnostic =
                    "SW92 asymmetric fixed pair: an assigned phase family is above the opposite-family Gibbs surface";
                return result;
            }

            result.status = Sw92AsymmetricFixedPairStatus::converged;
            result.diagnostic =
                "fixed family-pair equations and local lower-envelope family checks converged; overall phase-set acceptance not performed";
            return result;
        }

        if (result.iterations >= options.max_iterations) {
            result.status = Sw92AsymmetricFixedPairStatus::iteration_limit;
            result.diagnostic = "SW92 asymmetric fixed pair: iteration limit";
            return result;
        }

        const double scale = std::max(
            1.0, point.chemical_potential_norm / options.max_log_step);
        const double slope = detail::sw92_fixed_pair_gibbs_slope(
            point, result.feed, scale);
        if (!std::isfinite(slope) || !(slope < 0.0)) {
            result.status = Sw92AsymmetricFixedPairStatus::line_search_failed;
            result.diagnostic =
                "SW92 asymmetric fixed pair: degenerate or nonrepresentable Gibbs descent slope";
            return result;
        }

        bool accepted = false;
        double alpha = 1.0;
        for (int backtrack = 0; backtrack < options.max_backtracks;
             ++backtrack, alpha *= 0.5) {
            if (backtrack > 0) { ++result.backtracks; }
            std::vector<double> next_log_k(result.feed.size());
            for (std::size_t i = 0; i < result.feed.size(); ++i) {
                next_log_k[i] = point.log_k[i] +
                    (alpha / scale) * point.chemical_potential_residual[i];
            }
            try {
                Sw92AsymmetricFixedPairState next = evaluate_state(next_log_k);
                const double predicted = -alpha * slope;
                if (detail::sw92_fixed_pair_accept_step(
                        point, next, predicted, alpha / scale, options)) {
                    if (predicted > point.gibbs_roundoff_guard + next.gibbs_roundoff_guard) {
                        ++result.gibbs_descent_steps;
                    }
                    if (next.chemical_potential_norm > point.chemical_potential_norm) {
                        ++result.residual_increase_steps;
                    }
                    result.point = std::move(next);
                    accepted = true;
                    break;
                }
                ++result.rejected_evaluations;
            } catch (const detail::Sw92FixedPairStepError& error) {
                ++result.rejected_evaluations;
                result.diagnostic = error.what();
                if (error.status() == Sw92AsymmetricFixedPairStatus::evaluation_limit) {
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
            result.status = Sw92AsymmetricFixedPairStatus::line_search_failed;
            if (result.diagnostic.empty()) {
                result.diagnostic =
                    "SW92 asymmetric fixed pair: no acceptable Gibbs/residual line-search step";
            }
            return result;
        }
        ++result.iterations;
    }
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_ASYMMETRIC_FIXED_PAIR_HPP
