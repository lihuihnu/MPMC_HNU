#ifndef MPMC_PHYSICS_COMPONENT_INVENTORY_HPP
#define MPMC_PHYSICS_COMPONENT_INVENTORY_HPP

#include <mpmc/physics/thermodynamic_closure.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::physics {

inline constexpr std::string_view pt_component_inventory_convention =
    "PT/component-inventory/phase-resolved-reduced-feed/v1";

struct PtComponentInventoryOptions {
    // Integrity guards for algebraic identities inherited from an already
    // accepted thermodynamic closure. These are not flash convergence tolerances.
    double maximum_primal_identity_error{1.490116119384765625e-8};
    double maximum_linearization_identity_error{1.490116119384765625e-8};
};

struct PtComponentInventoryState {
    // Total moles per total fluid volume, not pore volume and not mass density.
    double total_molar_density_mol_per_m3{};
    // Component moles per total fluid volume, in ordered component coordinates.
    std::vector<double> component_molar_density_mol_per_m3;
};

struct PtComponentInventoryLinearization {
    std::size_t component_count{};
    std::size_t input_count{}; // N+1: p, T, z_0..z_{N-2}.
    std::vector<double> total_molar_density_gradient; // [q]
    std::vector<double> component_molar_density_jacobian; // [i*input_count+q]

    [[nodiscard]] static constexpr std::size_t pressure_column() noexcept {
        return 0U;
    }
    [[nodiscard]] static constexpr std::size_t temperature_column() noexcept {
        return 1U;
    }
    [[nodiscard]] std::size_t dependent_feed_component() const {
        if (component_count < 2U) {
            throw std::logic_error("component inventory: no reduced-feed chart");
        }
        return component_count - 1U;
    }
    [[nodiscard]] std::size_t feed_column(std::size_t independent_component) const {
        if (component_count < 2U ||
            independent_component + 1U >= component_count) {
            throw std::out_of_range(
                "component inventory: feed component outside reduced chart");
        }
        return 2U + independent_component;
    }
    [[nodiscard]] double d_total_molar_density(std::size_t column) const {
        return total_molar_density_gradient.at(column);
    }
    [[nodiscard]] double d_component_molar_density(
        std::size_t component, std::size_t column) const {
        return component_molar_density_jacobian.at(
            component * input_count + column);
    }
};

// Model-neutral local inventory derived from one already owned thermodynamic
// closure snapshot. This is a fluid-volume inventory only: no porosity, pore
// volume, saturation, mesh, flux, source term or time discretization is implied.
struct PtComponentInventorySnapshot {
    static constexpr std::string_view convention = pt_component_inventory_convention;

    ThermodynamicClosurePrimalStatus primal_status{
        ThermodynamicClosurePrimalStatus::indeterminate};
    ThermodynamicClosureLinearizationStatus linearization_status{
        ThermodynamicClosureLinearizationStatus::unavailable};
    ThermodynamicClosureLinearizationReason linearization_reason{
        ThermodynamicClosureLinearizationReason::not_implemented};

    std::size_t component_count{};
    std::size_t phase_count{};
    std::size_t input_count{};
    double pressure_pa{};
    double temperature_k{};
    std::vector<double> feed;
    std::vector<std::string> component_ids;
    std::string thermodynamic_model;
    std::string dataset_id;
    std::string revision;
    std::string source_closure_convention;

    std::optional<PtComponentInventoryState> primal;
    std::optional<PtComponentInventoryLinearization> linearization;
    std::string diagnostic;

    [[nodiscard]] bool inventory_available() const noexcept {
        return primal_status == ThermodynamicClosurePrimalStatus::valid &&
               primal.has_value() && component_count >= 2U &&
               component_ids.size() == component_count &&
               feed.size() == component_count &&
               primal->component_molar_density_mol_per_m3.size() == component_count &&
               std::isfinite(primal->total_molar_density_mol_per_m3) &&
               primal->total_molar_density_mol_per_m3 > 0.0;
    }

    [[nodiscard]] bool linearization_available() const noexcept {
        if (!inventory_available() ||
            linearization_status != ThermodynamicClosureLinearizationStatus::available ||
            linearization_reason != ThermodynamicClosureLinearizationReason::none ||
            !linearization) {
            return false;
        }
        const auto& value = *linearization;
        if (value.component_count != component_count ||
            value.input_count != component_count + 1U ||
            value.total_molar_density_gradient.size() != value.input_count) {
            return false;
        }
        if (value.component_count != 0U &&
            value.input_count > std::numeric_limits<std::size_t>::max() /
                                    value.component_count) {
            return false;
        }
        if (value.component_molar_density_jacobian.size() !=
            value.component_count * value.input_count) {
            return false;
        }
        const auto finite = [](const std::vector<double>& values) noexcept {
            return std::all_of(values.begin(), values.end(), [](double v) {
                return std::isfinite(v);
            });
        };
        return finite(value.total_molar_density_gradient) &&
               finite(value.component_molar_density_jacobian);
    }
};

namespace detail {

inline void validate_inventory_options(const PtComponentInventoryOptions& options) {
    if (!std::isfinite(options.maximum_primal_identity_error) ||
        !(options.maximum_primal_identity_error > 0.0) ||
        !std::isfinite(options.maximum_linearization_identity_error) ||
        !(options.maximum_linearization_identity_error > 0.0)) {
        throw std::invalid_argument("PT component inventory: invalid integrity options");
    }
}

[[nodiscard]] inline bool inventory_near(
    double lhs, double rhs, double relative_tolerance) {
    if (!std::isfinite(lhs) || !std::isfinite(rhs) ||
        !std::isfinite(relative_tolerance) || !(relative_tolerance > 0.0)) {
        return false;
    }
    return std::abs(lhs - rhs) <= relative_tolerance *
        std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

[[nodiscard]] inline double reduced_feed_derivative(
    std::size_t component_count, std::size_t component, std::size_t column) {
    if (column < 2U) { return 0.0; }
    const std::size_t independent = column - 2U;
    if (independent + 1U >= component_count) {
        throw std::out_of_range("component inventory: invalid reduced-feed column");
    }
    if (component == independent) { return 1.0; }
    if (component + 1U == component_count) { return -1.0; }
    return 0.0;
}

inline void inventory_clear_linearization(PtComponentInventorySnapshot& result) {
    result.linearization.reset();
    result.linearization_status = ThermodynamicClosureLinearizationStatus::unavailable;
}

inline void inventory_fail_primal(
    PtComponentInventorySnapshot& result, std::string diagnostic) {
    result.primal_status = ThermodynamicClosurePrimalStatus::indeterminate;
    result.primal.reset();
    inventory_clear_linearization(result);
    result.linearization_reason =
        ThermodynamicClosureLinearizationReason::solution_not_accepted;
    result.diagnostic = std::move(diagnostic);
}

template <class PhaseAt>
[[nodiscard]] inline bool inventory_build_primal(
    PtComponentInventorySnapshot& result, std::size_t phase_count,
    PhaseAt&& phase_at, const PtComponentInventoryOptions& options,
    std::vector<double>& phase_component_sum) {
    const std::size_t n = result.component_count;
    if (n < 2U || result.feed.size() != n || result.component_ids.size() != n ||
        phase_count == 0U || !std::isfinite(result.pressure_pa) ||
        !(result.pressure_pa > 0.0) || !std::isfinite(result.temperature_k) ||
        !(result.temperature_k > 0.0)) {
        inventory_fail_primal(result,
            "PT component inventory: invalid source shape or p/T");
        return false;
    }

    double feed_sum = 0.0;
    for (double zi : result.feed) {
        if (!std::isfinite(zi) || zi < 0.0) {
            inventory_fail_primal(result,
                "PT component inventory: invalid overall composition");
            return false;
        }
        feed_sum += zi;
    }
    if (!inventory_near(feed_sum, 1.0, options.maximum_primal_identity_error)) {
        inventory_fail_primal(result,
            "PT component inventory: overall composition is not normalized");
        return false;
    }

    phase_component_sum.assign(n, 0.0);
    double beta_sum = 0.0;
    double molar_volume_m3_per_mol = 0.0;
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const auto& state = phase_at(phase);
        if (!std::isfinite(state.mole_phase_fraction) ||
            state.mole_phase_fraction < 0.0 ||
            !std::isfinite(state.molar_density_mol_per_m3) ||
            !(state.molar_density_mol_per_m3 > 0.0) ||
            state.composition.size() != n) {
            inventory_fail_primal(result,
                "PT component inventory: invalid phase fraction/density/composition shape");
            return false;
        }
        double x_sum = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double xi = state.composition[i];
            if (!std::isfinite(xi) || xi < 0.0) {
                inventory_fail_primal(result,
                    "PT component inventory: invalid phase composition");
                return false;
            }
            x_sum += xi;
            phase_component_sum[i] += state.mole_phase_fraction * xi;
        }
        if (!inventory_near(x_sum, 1.0, options.maximum_primal_identity_error)) {
            inventory_fail_primal(result,
                "PT component inventory: phase composition is not normalized");
            return false;
        }
        beta_sum += state.mole_phase_fraction;
        molar_volume_m3_per_mol +=
            state.mole_phase_fraction / state.molar_density_mol_per_m3;
    }
    if (!inventory_near(beta_sum, 1.0, options.maximum_primal_identity_error) ||
        !std::isfinite(molar_volume_m3_per_mol) ||
        !(molar_volume_m3_per_mol > 0.0)) {
        inventory_fail_primal(result,
            "PT component inventory: phase fractions/mixture molar volume are invalid");
        return false;
    }

    for (std::size_t i = 0; i < n; ++i) {
        if (!inventory_near(phase_component_sum[i], result.feed[i],
                            options.maximum_primal_identity_error)) {
            inventory_fail_primal(result,
                "PT component inventory: phase-resolved material balance does not reproduce feed");
            return false;
        }
    }

    PtComponentInventoryState primal;
    primal.total_molar_density_mol_per_m3 = 1.0 / molar_volume_m3_per_mol;
    primal.component_molar_density_mol_per_m3.assign(n, 0.0);
    double component_density_sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double phase_form =
            primal.total_molar_density_mol_per_m3 * phase_component_sum[i];
        const double feed_form =
            primal.total_molar_density_mol_per_m3 * result.feed[i];
        if (!inventory_near(phase_form, feed_form,
                            options.maximum_primal_identity_error)) {
            inventory_fail_primal(result,
                "PT component inventory: phase/feed inventory identity failed");
            return false;
        }
        primal.component_molar_density_mol_per_m3[i] = phase_form;
        component_density_sum += phase_form;
    }
    if (!std::isfinite(primal.total_molar_density_mol_per_m3) ||
        !(primal.total_molar_density_mol_per_m3 > 0.0) ||
        !inventory_near(component_density_sum,
                        primal.total_molar_density_mol_per_m3,
                        options.maximum_primal_identity_error)) {
        inventory_fail_primal(result,
            "PT component inventory: component densities do not close to total density");
        return false;
    }

    result.phase_count = phase_count;
    result.primal = std::move(primal);
    result.primal_status = ThermodynamicClosurePrimalStatus::valid;
    return true;
}

template <class PhaseAt, class DBeta, class DX, class DC>
[[nodiscard]] inline bool inventory_build_linearization(
    PtComponentInventorySnapshot& result, PhaseAt&& phase_at,
    DBeta&& d_beta, DX&& d_x, DC&& d_c,
    std::size_t source_input_count,
    const std::vector<double>& phase_component_sum,
    const PtComponentInventoryOptions& options) {
    if (!result.primal || source_input_count != result.component_count + 1U) {
        inventory_clear_linearization(result);
        result.linearization_reason =
            ThermodynamicClosureLinearizationReason::solution_not_accepted;
        result.diagnostic =
            "PT component inventory: source linearization shape is inconsistent";
        return false;
    }
    const std::size_t n = result.component_count;
    const std::size_t q_count = source_input_count;
    if (n > std::numeric_limits<std::size_t>::max() / q_count) {
        inventory_clear_linearization(result);
        result.linearization_reason =
            ThermodynamicClosureLinearizationReason::arithmetic_failure;
        result.diagnostic = "PT component inventory: Jacobian size overflow";
        return false;
    }

    const double total_density = result.primal->total_molar_density_mol_per_m3;
    PtComponentInventoryLinearization linearization;
    linearization.component_count = n;
    linearization.input_count = q_count;
    linearization.total_molar_density_gradient.assign(q_count, 0.0);
    linearization.component_molar_density_jacobian.assign(n * q_count, 0.0);

    for (std::size_t column = 0; column < q_count; ++column) {
        double d_molar_volume = 0.0;
        for (std::size_t phase = 0; phase < result.phase_count; ++phase) {
            const auto& state = phase_at(phase);
            const double db = d_beta(phase, column);
            const double dc = d_c(phase, column);
            if (!std::isfinite(db) || !std::isfinite(dc)) {
                inventory_clear_linearization(result);
                result.linearization_reason =
                    ThermodynamicClosureLinearizationReason::arithmetic_failure;
                result.diagnostic =
                    "PT component inventory: nonfinite phase-fraction/density derivative";
                return false;
            }
            d_molar_volume += db / state.molar_density_mol_per_m3 -
                state.mole_phase_fraction * dc /
                    (state.molar_density_mol_per_m3 *
                     state.molar_density_mol_per_m3);
        }
        const double d_total_density =
            -total_density * total_density * d_molar_volume;
        if (!std::isfinite(d_total_density)) {
            inventory_clear_linearization(result);
            result.linearization_reason =
                ThermodynamicClosureLinearizationReason::arithmetic_failure;
            result.diagnostic =
                "PT component inventory: nonfinite total-density derivative";
            return false;
        }
        linearization.total_molar_density_gradient[column] = d_total_density;

        double component_derivative_sum = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double d_phase_component_sum = 0.0;
            for (std::size_t phase = 0; phase < result.phase_count; ++phase) {
                const auto& state = phase_at(phase);
                const double db = d_beta(phase, column);
                const double dx = d_x(phase, i, column);
                if (!std::isfinite(dx)) {
                    inventory_clear_linearization(result);
                    result.linearization_reason =
                        ThermodynamicClosureLinearizationReason::arithmetic_failure;
                    result.diagnostic =
                        "PT component inventory: nonfinite composition derivative";
                    return false;
                }
                d_phase_component_sum +=
                    db * state.composition[i] +
                    state.mole_phase_fraction * dx;
            }
            const double phase_form =
                d_total_density * phase_component_sum[i] +
                total_density * d_phase_component_sum;
            const double feed_form =
                d_total_density * result.feed[i] +
                total_density * reduced_feed_derivative(n, i, column);
            if (!std::isfinite(phase_form) ||
                !inventory_near(phase_form, feed_form,
                                options.maximum_linearization_identity_error)) {
                inventory_clear_linearization(result);
                result.linearization_reason =
                    ThermodynamicClosureLinearizationReason::arithmetic_failure;
                result.diagnostic =
                    "PT component inventory: differentiated phase/feed identity failed";
                return false;
            }
            linearization.component_molar_density_jacobian[
                i * q_count + column] = phase_form;
            component_derivative_sum += phase_form;
        }
        if (!inventory_near(component_derivative_sum, d_total_density,
                            options.maximum_linearization_identity_error)) {
            inventory_clear_linearization(result);
            result.linearization_reason =
                ThermodynamicClosureLinearizationReason::arithmetic_failure;
            result.diagnostic =
                "PT component inventory: differentiated component closure failed";
            return false;
        }
    }

    result.input_count = q_count;
    result.linearization = std::move(linearization);
    result.linearization_status = ThermodynamicClosureLinearizationStatus::available;
    result.linearization_reason = ThermodynamicClosureLinearizationReason::none;
    result.diagnostic =
        "PT component inventory: phase-resolved primal and local Jacobian available";
    return true;
}

inline ThermodynamicClosureLinearizationReason propagated_reason(
    ThermodynamicClosureLinearizationStatus status,
    ThermodynamicClosureLinearizationReason reason) {
    if (status == ThermodynamicClosureLinearizationStatus::available) {
        return ThermodynamicClosureLinearizationReason::solution_not_accepted;
    }
    return reason == ThermodynamicClosureLinearizationReason::none
        ? ThermodynamicClosureLinearizationReason::solution_not_accepted
        : reason;
}

} // namespace detail

[[nodiscard]] inline PtComponentInventorySnapshot build_pt_component_inventory(
    const ThermodynamicClosureSnapshot& source,
    PtComponentInventoryOptions options = {}) {
    detail::validate_inventory_options(options);
    PtComponentInventorySnapshot result;
    result.pressure_pa = source.pressure_pa;
    result.temperature_k = source.temperature_k;
    result.feed = source.feed;
    result.component_ids = source.component_ids;
    result.component_count = source.component_ids.size();
    result.thermodynamic_model = source.thermodynamic_model;
    result.dataset_id = source.dataset_id;
    result.revision = source.revision;
    result.source_closure_convention =
        std::string(ThermodynamicClosureSnapshot::convention);

    if (!source.residual_available() || !source.primal) {
        result.primal_status = source.primal_status ==
            ThermodynamicClosurePrimalStatus::valid
            ? ThermodynamicClosurePrimalStatus::indeterminate
            : source.primal_status;
        result.linearization_reason = detail::propagated_reason(
            source.linearization_status, source.linearization_reason);
        result.diagnostic =
            "PT component inventory: fixed-VLE source primal unavailable";
        return result;
    }

    const auto phase_at = [&](std::size_t phase) -> const ThermodynamicPhaseState& {
        if (phase == 0U) { return source.primal->liquid; }
        if (phase == 1U) { return source.primal->vapor; }
        throw std::out_of_range("component inventory: fixed-VLE phase index");
    };
    std::vector<double> phase_component_sum;
    if (!detail::inventory_build_primal(
            result, 2U, phase_at, options, phase_component_sum)) {
        return result;
    }

    if (source.linearization_status != ThermodynamicClosureLinearizationStatus::available ||
        !source.linearization) {
        result.linearization_reason = detail::propagated_reason(
            source.linearization_status, source.linearization_reason);
        result.diagnostic =
            "PT component inventory: primal available; fixed-VLE source linearization unavailable";
        return result;
    }
    if (source.linearization_reason != ThermodynamicClosureLinearizationReason::none ||
        !source.can_seed_newton()) {
        result.linearization_reason =
            ThermodynamicClosureLinearizationReason::solution_not_accepted;
        result.diagnostic =
            "PT component inventory: malformed fixed-VLE source linearization";
        return result;
    }

    const auto& linearization = *source.linearization;
    const std::size_t n = result.component_count;
    const std::size_t q = linearization.input_count;
    if (n > std::numeric_limits<std::size_t>::max() / q ||
        linearization.component_count != n || q != n + 1U ||
        linearization.vapor_fraction_gradient.size() != q ||
        linearization.liquid_composition_jacobian.size() != n * q ||
        linearization.vapor_composition_jacobian.size() != n * q ||
        linearization.liquid_molar_density_gradient.size() != q ||
        linearization.vapor_molar_density_gradient.size() != q) {
        result.linearization_reason =
            ThermodynamicClosureLinearizationReason::solution_not_accepted;
        result.diagnostic =
            "PT component inventory: fixed-VLE derivative payload shape mismatch";
        return result;
    }
    const auto d_beta = [&](std::size_t phase, std::size_t column) {
        const double dv = linearization.vapor_fraction_gradient[column];
        return phase == 0U ? -dv : dv;
    };
    const auto d_x = [&](std::size_t phase, std::size_t component,
                         std::size_t column) {
        const auto& jacobian = phase == 0U
            ? linearization.liquid_composition_jacobian
            : linearization.vapor_composition_jacobian;
        return jacobian[component * q + column];
    };
    const auto d_c = [&](std::size_t phase, std::size_t column) {
        return phase == 0U
            ? linearization.liquid_molar_density_gradient[column]
            : linearization.vapor_molar_density_gradient[column];
    };
    (void)detail::inventory_build_linearization(
        result, phase_at, d_beta, d_x, d_c, q, phase_component_sum, options);
    return result;
}

[[nodiscard]] inline PtComponentInventorySnapshot build_pt_component_inventory(
    const PtPhaseSetThermodynamicClosureSnapshot& source,
    PtComponentInventoryOptions options = {}) {
    detail::validate_inventory_options(options);
    PtComponentInventorySnapshot result;
    result.pressure_pa = source.pressure_pa;
    result.temperature_k = source.temperature_k;
    result.feed = source.feed;
    result.component_ids = source.component_ids;
    result.component_count = source.component_ids.size();
    result.thermodynamic_model = source.thermodynamic_model;
    result.dataset_id = source.dataset_id;
    result.revision = source.revision;
    result.source_closure_convention =
        std::string(PtPhaseSetThermodynamicClosureSnapshot::convention);

    if (!source.residual_available() || !source.primal) {
        result.primal_status = source.primal_status ==
            ThermodynamicClosurePrimalStatus::valid
            ? ThermodynamicClosurePrimalStatus::indeterminate
            : source.primal_status;
        result.linearization_reason = detail::propagated_reason(
            source.linearization_status, source.linearization_reason);
        result.diagnostic =
            "PT component inventory: phase-set source primal unavailable";
        return result;
    }

    const auto phase_at = [&](std::size_t phase) -> const ThermodynamicPhaseState& {
        return source.primal->phases.at(phase);
    };
    std::vector<double> phase_component_sum;
    if (!detail::inventory_build_primal(
            result, source.primal->phases.size(), phase_at,
            options, phase_component_sum)) {
        return result;
    }

    if (source.linearization_status != ThermodynamicClosureLinearizationStatus::available ||
        !source.linearization) {
        result.linearization_reason = detail::propagated_reason(
            source.linearization_status, source.linearization_reason);
        result.diagnostic =
            "PT component inventory: primal available; phase-set source linearization unavailable";
        return result;
    }
    if (!source.can_seed_newton()) {
        result.linearization_reason =
            ThermodynamicClosureLinearizationReason::solution_not_accepted;
        result.diagnostic =
            "PT component inventory: malformed phase-set source linearization";
        return result;
    }

    const auto& linearization = *source.linearization;
    const std::size_t q = linearization.input_count;
    const auto d_beta = [&](std::size_t phase, std::size_t column) {
        return linearization.d_phase_fraction(phase, column);
    };
    const auto d_x = [&](std::size_t phase, std::size_t component,
                         std::size_t column) {
        return linearization.d_composition(phase, component, column);
    };
    const auto d_c = [&](std::size_t phase, std::size_t column) {
        return linearization.d_molar_density(phase, column);
    };
    (void)detail::inventory_build_linearization(
        result, phase_at, d_beta, d_x, d_c, q, phase_component_sum, options);
    return result;
}

} // namespace mpmc::physics

#endif // MPMC_PHYSICS_COMPONENT_INVENTORY_HPP
