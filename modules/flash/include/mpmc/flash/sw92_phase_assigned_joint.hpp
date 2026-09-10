#ifndef MPMC_FLASH_SW92_PHASE_ASSIGNED_JOINT_HPP
#define MPMC_FLASH_SW92_PHASE_ASSIGNED_JOINT_HPP

#include <mpmc/flash/rachford_rice.hpp>
#include <mpmc/flash/sw92_stability.hpp>

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

inline constexpr std::string_view sw92_phase_assigned_aq_na_joint_profile =
    "SW92-equilibrium/phase-assigned-aq-na-joint/v1";
inline constexpr std::string_view sw92_phase_assigned_aq_na_joint_primitive =
    "SW92-equilibrium/phase-assigned-aq-na-joint/fixed-two-phase-logK-SSI-RR/v1";

enum class Sw92PhysicalPhaseRole { aqueous, nonaqueous };

struct Sw92PhaseAssignedJointOptions {
    double chemical_potential_tolerance{1e-11};
    double mass_absolute_tolerance{1e-12};
    double mass_relative_tolerance{1e-10};
    double minimum_phase_fraction{1e-10};
    double log_k_separation{1e-7};
    double max_log_step{2.0};
    double residual_decrease{1e-4};
    double gibbs_armijo{1e-4};
    int max_iterations{512};
    int max_backtracks{24};
    std::size_t max_evaluations{20000};
    RachfordRiceOptions rr;
    thermodynamics::Sw92RootOptions aqueous_root_options;
    thermodynamics::Sw92RootOptions nonaqueous_root_options;
};

struct Sw92PhaseAssignedJointPhase {
    Sw92PhysicalPhaseRole physical_role{Sw92PhysicalPhaseRole::aqueous};
    thermodynamics::SwPhaseFamily thermodynamic_family{
        thermodynamics::SwPhaseFamily::aqueous};
    double mole_phase_fraction{};
    std::vector<double> composition;
    StabilityPhase activity;
    double compressibility_factor{};
};

struct Sw92PhaseAssignedJointState {
    Sw92PhaseAssignedJointPhase aqueous_phase;
    Sw92PhaseAssignedJointPhase nonaqueous_phase;

    // ln(x_i^NA / x_i^AQ) on active feed support.
    std::vector<double> log_k;
    std::vector<double> chemical_potential_residual; // mu_AQ/RT - mu_NA/RT.
    std::vector<double> common_log_activity;
    double chemical_potential_norm{};
    double reduced_gibbs{};
    double gibbs_roundoff_guard{};

    double rr_residual{};
    double raw_aqueous_sum{};
    double raw_nonaqueous_sum{};
    double mass_absolute{};
    double mass_relative{};
    int rr_iterations{};

    double aqueous_minus_nonaqueous_water_fraction{
        std::numeric_limits<double>::quiet_NaN()};
    double water_role_roundoff_guard{
        std::numeric_limits<double>::quiet_NaN()};
};

enum class Sw92PhaseAssignedJointStatus {
    converged_candidate,
    no_interior_rr_root,
    rr_failure,
    indistinguishable_phases,
    phase_disappearance,
    balance_failure,
    family_root_nonsmooth,
    phase_role_indeterminate,
    phase_role_reversed,
    property_failure,
    iteration_limit,
    evaluation_limit,
    line_search_failed
};

struct Sw92PhaseAssignedJointResult {
    static constexpr std::size_t maximum_phase_count = 2;
    static constexpr bool final_stability_checked = false;
    static constexpr bool global_stability_proven = false;

    Sw92PhaseAssignedJointStatus status{Sw92PhaseAssignedJointStatus::rr_failure};
    Sw92PhaseAssignedJointOptions options;
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
    std::string primitive_convention{sw92_phase_assigned_aq_na_joint_primitive};

    std::vector<double> initial_log_k;
    std::optional<Sw92PhaseAssignedJointState> point;
    int iterations{};
    std::size_t evaluations{};
    std::size_t backtracks{};
    std::size_t rejected_evaluations{};
    std::size_t gibbs_descent_steps{};
    std::size_t residual_increase_steps{};
    std::optional<StabilityPropertyIssue> property_issue;
    std::string diagnostic;

    [[nodiscard]] bool equations_converged() const noexcept {
        return point.has_value() &&
               point->chemical_potential_norm <= options.chemical_potential_tolerance &&
               point->mass_absolute <= options.mass_absolute_tolerance &&
               point->mass_relative <= options.mass_relative_tolerance;
    }

    [[nodiscard]] bool candidate_admissible() const noexcept {
        return status == Sw92PhaseAssignedJointStatus::converged_candidate &&
               equations_converged();
    }

    [[nodiscard]] const Sw92PhaseAssignedJointState* candidate() const & noexcept {
        return candidate_admissible() ? &*point : nullptr;
    }
    const Sw92PhaseAssignedJointState* candidate() const && = delete;
};

namespace detail {

inline void sw92_phase_assigned_check_options(
    const Sw92PhaseAssignedJointOptions& options) {
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
            "SW92 phase-assigned joint: invalid tolerance, step or resource option");
    }
    rr_check_options(options.rr);
}

class Sw92PhaseAssignedStepError : public std::runtime_error {
public:
    Sw92PhaseAssignedStepError(Sw92PhaseAssignedJointStatus status, const char* message)
        : std::runtime_error(message), status_(status) {}
    [[nodiscard]] Sw92PhaseAssignedJointStatus status() const noexcept { return status_; }
private:
    Sw92PhaseAssignedJointStatus status_;
};

inline Sw92PhaseAssignedJointStatus sw92_phase_assigned_rr_status(
    RachfordRiceStatus status) noexcept {
    switch (status) {
    case RachfordRiceStatus::interior:
        return Sw92PhaseAssignedJointStatus::converged_candidate;
    case RachfordRiceStatus::no_resolved_interior_root:
        return Sw92PhaseAssignedJointStatus::no_interior_rr_root;
    case RachfordRiceStatus::degenerate:
        return Sw92PhaseAssignedJointStatus::indistinguishable_phases;
    case RachfordRiceStatus::iteration_limit:
    case RachfordRiceStatus::unrepresentable:
        return Sw92PhaseAssignedJointStatus::rr_failure;
    }
    return Sw92PhaseAssignedJointStatus::rr_failure;
}

inline bool sw92_phase_assigned_balance_ok(
    const Sw92PhaseAssignedJointState& state,
    const Sw92PhaseAssignedJointOptions& options) noexcept {
    return state.mass_absolute <= options.mass_absolute_tolerance &&
           state.mass_relative <= options.mass_relative_tolerance;
}

inline double sw92_phase_assigned_gibbs_slope(
    const Sw92PhaseAssignedJointState& point,
    std::span<const double> feed, double scale) {
    double h = 0.0;
    double h_correction = 0.0;
    double first = 0.0;
    double first_correction = 0.0;
    double second = 0.0;
    double second_correction = 0.0;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (feed[i] == 0.0) { continue; }
        const double x = point.aqueous_phase.composition[i];
        const double y = point.nonaqueous_phase.composition[i];
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
    const double beta = point.nonaqueous_phase.mole_phase_fraction;
    return -(beta * (1.0 - beta) * second + (first / h) * first) / scale;
}

inline bool sw92_phase_assigned_accept_step(
    const Sw92PhaseAssignedJointState& point,
    const Sw92PhaseAssignedJointState& next,
    double predicted, double relative_step,
    const Sw92PhaseAssignedJointOptions& options) noexcept {
    const double guard = point.gibbs_roundoff_guard + next.gibbs_roundoff_guard;
    if (predicted > guard) {
        return next.reduced_gibbs <=
            point.reduced_gibbs - options.gibbs_armijo * predicted;
    }
    return next.reduced_gibbs <= point.reduced_gibbs + guard &&
           (next.chemical_potential_norm <= options.chemical_potential_tolerance ||
            point.chemical_potential_norm - next.chemical_potential_norm >
                options.residual_decrease * relative_step * point.chemical_potential_norm);
}

} // namespace detail

/// Solve a constrained two-phase state with fixed physical role/model mapping:
/// aqueous role -> SW92 AQ family, nonaqueous role -> SW92 NA family.
///
/// This is Profile-C Gate C1 only. It performs no autonomous phase-number or
/// single-phase classification and no final constrained two-family stability review.
[[nodiscard]] inline Sw92PhaseAssignedJointResult
iterate_sw92_phase_assigned_aq_na_joint(
    double pressure_pa, double temperature_k, std::span<const double> feed,
    std::span<const double> initial_log_k,
    const thermodynamics::Sw92Phase<double>& model,
    double nacl_molality_mol_per_kg_water,
    Sw92PhaseAssignedJointOptions options = {}) {
    detail::sw92_phase_assigned_check_options(options);
    if (!std::isfinite(pressure_pa) || pressure_pa <= 0.0 ||
        !std::isfinite(temperature_k) || temperature_k <= 0.0) {
        throw std::domain_error(
            "SW92 phase-assigned joint: finite p>0 Pa and T>0 K required");
    }
    if (feed.size() != model.size()) {
        throw std::invalid_argument(
            "SW92 phase-assigned joint: feed does not match ordered model snapshot");
    }
    if (feed.size() != initial_log_k.size()) {
        throw std::invalid_argument(
            "SW92 phase-assigned joint: logK dimension mismatch");
    }
    if (feed.empty() || feed.size() > options.rr.max_components) {
        throw std::length_error(
            "SW92 phase-assigned joint: component quota exceeded or empty feed");
    }
    for (double value : initial_log_k) {
        if (!std::isfinite(value)) {
            throw std::domain_error(
                "SW92 phase-assigned joint: finite initial logK required");
        }
    }

    Sw92PhaseAssignedJointResult result;
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
            result.status = Sw92PhaseAssignedJointStatus::rr_failure;
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
            throw detail::Sw92PhaseAssignedStepError(
                Sw92PhaseAssignedJointStatus::evaluation_limit,
                "SW92 phase-assigned joint: property budget exhausted");
        }
        ++result.evaluations;
        auto selected = evaluator.evaluate_selected(
            pressure_pa, temperature_k, composition);
        if (!selected.activity.smooth) {
            throw detail::Sw92PhaseAssignedStepError(
                Sw92PhaseAssignedJointStatus::family_root_nonsmooth,
                "SW92 phase-assigned joint: same-family minimum-root envelope is nonsmooth");
        }
        if (!std::isfinite(selected.compressibility_factor) ||
            !(selected.compressibility_factor > 0.0)) {
            throw StabilityPropertyError(
                StabilityPropertyIssue::nonfinite_properties,
                "SW92 phase-assigned joint: selected Z is nonrepresentable");
        }
        return selected;
    };

    const auto evaluate_state = [&](std::span<const double> log_k) {
        const auto rr = solve_rachford_rice(result.feed, log_k, options.rr);
        if (rr.status != RachfordRiceStatus::interior) {
            throw detail::Sw92PhaseAssignedStepError(
                detail::sw92_phase_assigned_rr_status(rr.status),
                "SW92 phase-assigned joint: Rachford-Rice has no usable interior state");
        }

        Sw92PhaseAssignedJointState state;
        state.aqueous_phase.physical_role = Sw92PhysicalPhaseRole::aqueous;
        state.aqueous_phase.thermodynamic_family =
            thermodynamics::SwPhaseFamily::aqueous;
        state.nonaqueous_phase.physical_role = Sw92PhysicalPhaseRole::nonaqueous;
        state.nonaqueous_phase.thermodynamic_family =
            thermodynamics::SwPhaseFamily::nonaqueous;

        state.nonaqueous_phase.mole_phase_fraction = rr.vapor_fraction;
        state.aqueous_phase.mole_phase_fraction = 1.0 - rr.vapor_fraction;
        state.aqueous_phase.composition = rr.liquid;
        state.nonaqueous_phase.composition = rr.vapor;
        state.rr_residual = rr.residual;
        state.raw_aqueous_sum = rr.raw_liquid_sum;
        state.raw_nonaqueous_sum = rr.raw_vapor_sum;
        state.mass_absolute = rr.mass_absolute;
        state.mass_relative = rr.mass_relative;
        state.rr_iterations = rr.iterations;

        const auto selected_aq = evaluate_selected(
            aqueous_evaluator, state.aqueous_phase.composition);
        const auto selected_na = evaluate_selected(
            nonaqueous_evaluator, state.nonaqueous_phase.composition);
        state.aqueous_phase.activity = selected_aq.activity;
        state.aqueous_phase.compressibility_factor =
            selected_aq.compressibility_factor;
        state.nonaqueous_phase.activity = selected_na.activity;
        state.nonaqueous_phase.compressibility_factor =
            selected_na.compressibility_factor;

        const std::size_t n = result.feed.size();
        state.log_k.resize(n);
        state.chemical_potential_residual.resize(n);
        state.common_log_activity.resize(n);
        double gibbs_correction = 0.0;
        double magnitude = 1.0;
        double magnitude_correction = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            if (result.feed[i] == 0.0) { continue; }
            const double lx = std::log(state.aqueous_phase.composition[i]);
            const double ly = std::log(state.nonaqueous_phase.composition[i]);
            const double m_aq = lx + state.aqueous_phase.activity.ln_phi[i];
            const double m_na = ly + state.nonaqueous_phase.activity.ln_phi[i];
            state.log_k[i] = ly - lx;
            state.chemical_potential_residual[i] = m_aq - m_na;
            state.common_log_activity[i] = std::midpoint(m_aq, m_na);
            state.chemical_potential_norm = std::max(
                state.chemical_potential_norm,
                std::abs(state.chemical_potential_residual[i]));
            detail::stability_add(
                state.aqueous_phase.mole_phase_fraction *
                    state.aqueous_phase.composition[i] * m_aq +
                state.nonaqueous_phase.mole_phase_fraction *
                    state.nonaqueous_phase.composition[i] * m_na,
                state.reduced_gibbs, gibbs_correction);
            detail::stability_add(
                state.aqueous_phase.mole_phase_fraction *
                    state.aqueous_phase.composition[i] *
                    (std::abs(lx) +
                     std::abs(state.aqueous_phase.activity.ln_phi[i])) +
                state.nonaqueous_phase.mole_phase_fraction *
                    state.nonaqueous_phase.composition[i] *
                    (std::abs(ly) +
                     std::abs(state.nonaqueous_phase.activity.ln_phi[i])),
                magnitude, magnitude_correction);
            if (!std::isfinite(state.log_k[i]) ||
                !std::isfinite(state.chemical_potential_residual[i]) ||
                !std::isfinite(state.common_log_activity[i])) {
                throw StabilityPropertyError(
                    StabilityPropertyIssue::nonfinite_properties,
                    "SW92 phase-assigned joint: nonrepresentable activity arithmetic");
            }
        }
        state.gibbs_roundoff_guard =
            256.0 * detail::stability_eps * magnitude;
        if (!std::isfinite(state.reduced_gibbs) ||
            !std::isfinite(state.gibbs_roundoff_guard)) {
            throw StabilityPropertyError(
                StabilityPropertyIssue::nonfinite_properties,
                "SW92 phase-assigned joint: nonrepresentable Gibbs arithmetic");
        }
        return state;
    };

    try {
        result.point = evaluate_state(result.initial_log_k);
    } catch (const detail::Sw92PhaseAssignedStepError& error) {
        result.status = error.status();
        result.diagnostic = error.what();
        return result;
    } catch (const StabilityPropertyError& error) {
        result.status = Sw92PhaseAssignedJointStatus::property_failure;
        result.property_issue = error.issue();
        result.diagnostic = error.what();
        return result;
    }

    for (;;) {
        auto& point = *result.point;
        if (point.chemical_potential_norm <= options.chemical_potential_tolerance) {
            if (!detail::sw92_phase_assigned_balance_ok(point, options)) {
                result.status = Sw92PhaseAssignedJointStatus::balance_failure;
                result.diagnostic =
                    "SW92 phase-assigned joint: chemical potentials converged but material balance failed";
                return result;
            }
            if (point.aqueous_phase.mole_phase_fraction <=
                    options.minimum_phase_fraction ||
                point.nonaqueous_phase.mole_phase_fraction <=
                    options.minimum_phase_fraction) {
                result.status = Sw92PhaseAssignedJointStatus::phase_disappearance;
                result.diagnostic =
                    "SW92 phase-assigned joint: balance and chemical-potential tolerances passed at the phase-disappearance threshold";
                return result;
            }
            double contrast = 0.0;
            for (std::size_t i = 0; i < result.feed.size(); ++i) {
                if (result.feed[i] > 0.0) {
                    contrast = std::max(contrast, std::abs(point.log_k[i]));
                }
            }
            if (contrast <= options.log_k_separation) {
                result.status = Sw92PhaseAssignedJointStatus::indistinguishable_phases;
                result.diagnostic =
                    "SW92 phase-assigned joint: converged phase compositions are numerically indistinguishable";
                return result;
            }

            const std::size_t water_index = parameters.water_index();
            const double aq_water = point.aqueous_phase.composition[water_index];
            const double na_water = point.nonaqueous_phase.composition[water_index];
            const double scale = 1.0 + std::abs(aq_water) + std::abs(na_water);
            point.water_role_roundoff_guard =
                256.0 * detail::stability_eps * scale;
            point.aqueous_minus_nonaqueous_water_fraction = aq_water - na_water;
            if (!std::isfinite(point.aqueous_minus_nonaqueous_water_fraction) ||
                !std::isfinite(point.water_role_roundoff_guard)) {
                result.status = Sw92PhaseAssignedJointStatus::phase_role_indeterminate;
                result.diagnostic =
                    "SW92 phase-assigned joint: water-role ordering is nonrepresentable";
                return result;
            }
            if (std::abs(point.aqueous_minus_nonaqueous_water_fraction) <=
                point.water_role_roundoff_guard) {
                result.status = Sw92PhaseAssignedJointStatus::phase_role_indeterminate;
                result.diagnostic =
                    "SW92 phase-assigned joint: aqueous/nonaqueous water-richness ordering is unresolved at roundoff scale";
                return result;
            }
            if (point.aqueous_minus_nonaqueous_water_fraction < 0.0) {
                result.status = Sw92PhaseAssignedJointStatus::phase_role_reversed;
                result.diagnostic =
                    "SW92 phase-assigned joint: assigned aqueous phase is not the water-richer phase";
                return result;
            }

            result.status = Sw92PhaseAssignedJointStatus::converged_candidate;
            result.diagnostic =
                "fixed phase-assigned AQ/NA equations, material balance and relative-water topology converged; final constrained stability not performed";
            return result;
        }

        if (result.iterations >= options.max_iterations) {
            result.status = Sw92PhaseAssignedJointStatus::iteration_limit;
            result.diagnostic = "SW92 phase-assigned joint: iteration limit";
            return result;
        }

        const double scale = std::max(
            1.0, point.chemical_potential_norm / options.max_log_step);
        const double slope = detail::sw92_phase_assigned_gibbs_slope(
            point, result.feed, scale);
        if (!std::isfinite(slope) || !(slope < 0.0)) {
            result.status = Sw92PhaseAssignedJointStatus::line_search_failed;
            result.diagnostic =
                "SW92 phase-assigned joint: degenerate or nonrepresentable Gibbs descent slope";
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
                Sw92PhaseAssignedJointState next = evaluate_state(next_log_k);
                const double predicted = -alpha * slope;
                if (detail::sw92_phase_assigned_accept_step(
                        point, next, predicted, alpha / scale, options)) {
                    if (predicted >
                        point.gibbs_roundoff_guard + next.gibbs_roundoff_guard) {
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
            } catch (const detail::Sw92PhaseAssignedStepError& error) {
                ++result.rejected_evaluations;
                result.diagnostic = error.what();
                if (error.status() == Sw92PhaseAssignedJointStatus::evaluation_limit) {
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
            result.status = Sw92PhaseAssignedJointStatus::line_search_failed;
            if (result.diagnostic.empty()) {
                result.diagnostic =
                    "SW92 phase-assigned joint: no acceptable Gibbs/residual line-search step";
            }
            return result;
        }
        ++result.iterations;
    }
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_PHASE_ASSIGNED_JOINT_HPP
