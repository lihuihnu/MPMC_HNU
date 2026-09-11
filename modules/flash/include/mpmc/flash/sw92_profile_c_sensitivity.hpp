#ifndef MPMC_FLASH_SW92_PROFILE_C_SENSITIVITY_HPP
#define MPMC_FLASH_SW92_PROFILE_C_SENSITIVITY_HPP

#include <mpmc/flash/detail/pt_sensitivity_detail.hpp>
#include <mpmc/flash/detail/sw92_profile_c_sensitivity_local.hpp>
#include <mpmc/flash/pt_sensitivity.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mpmc::flash {

inline constexpr std::string_view sw92_profile_c_sensitivity_convention =
    "SW92/Profile-C/fixed-authoritative-phase-set/implicit-reduced-feed-v1";

struct Sw92ProfileCSensitivityOptions {
    // Same sqrt(epsilon) derivative-domain guard used by the PR76 sensitivity
    // contract. This is not a phase-existence threshold.
    double minimum_derivative_phase_fraction{1.490116119384765625e-8};
    // Reject a local representation in which two accepted phase compositions
    // are too close to identify a smooth fixed phase slot.
    double minimum_log_composition_separation{1.490116119384765625e-8};
    // Reproduction guard for the already accepted equilibrium equations. This
    // is a sensitivity-integrity check, not a flash convergence tolerance.
    double maximum_base_residual{1.490116119384765625e-8};
    double minimum_reciprocal_condition{0.0};
    std::size_t max_components{256};
    thermodynamics::Sw92RootOptions root_options{};
};

struct Sw92ProfileCPhaseSetSensitivityResult {
    PtSensitivityStatus status{PtSensitivityStatus::arithmetic_failure};
    std::string convention{sw92_profile_c_sensitivity_convention};

    std::size_t component_count{};
    std::size_t phase_count{};
    std::size_t input_count{}; // N+1: p,T,z_0..z_{N-2}.
    double pressure_pa{};
    double temperature_k{};
    std::vector<double> feed;

    // Row-major derivatives. Phase order is the authoritative publication's
    // representation order. H slots remain unordered physical NA instances.
    std::vector<double> phase_fraction_jacobian; // [phase * input_count + q]
    std::vector<double> composition_jacobian;    // [(phase*N+i)*input_count + q]
    std::vector<double> compressibility_jacobian;// [phase * input_count + q]
    std::vector<double> molar_density_jacobian;  // [phase * input_count + q]

    double equilibrium_jacobian_rcond{};
    double linear_solve_backward_error{};
    double equilibrium_residual_norm{};

    double nacl_molality_mol_per_kg_water{};
    std::string dataset_id;
    std::string revision;
    std::vector<std::string> component_ids;
    std::string model_profile;
    std::string phase_convention;
    std::string equilibrium_profile;
    std::string orchestration_convention;
    std::string boundary_convention;
    std::string publication_convention;
    std::vector<Sw92ProfileCPhaseMetadata> phase_metadata;
    std::string diagnostic;

    [[nodiscard]] static constexpr std::size_t pressure_column() noexcept {
        return 0U;
    }
    [[nodiscard]] static constexpr std::size_t temperature_column() noexcept {
        return 1U;
    }
    [[nodiscard]] std::size_t dependent_feed_component() const {
        if (component_count < 2U)
            throw std::logic_error("SW92 sensitivity: no reduced feed chart");
        return component_count - 1U;
    }
    [[nodiscard]] std::size_t feed_column(
        std::size_t independent_component) const {
        if (component_count < 2U ||
            independent_component + 1U >= component_count)
            throw std::out_of_range(
                "SW92 sensitivity: feed coordinate outside reduced chart");
        return 2U + independent_component;
    }
    [[nodiscard]] double d_phase_fraction(
        std::size_t phase, std::size_t column) const {
        return phase_fraction_jacobian.at(phase * input_count + column);
    }
    [[nodiscard]] double d_composition(
        std::size_t phase, std::size_t component,
        std::size_t column) const {
        return composition_jacobian.at(
            (phase * component_count + component) * input_count + column);
    }
    [[nodiscard]] double d_compressibility(
        std::size_t phase, std::size_t column) const {
        return compressibility_jacobian.at(phase * input_count + column);
    }
    [[nodiscard]] double d_molar_density(
        std::size_t phase, std::size_t column) const {
        return molar_density_jacobian.at(phase * input_count + column);
    }
};

namespace detail {

inline void sw92_profile_c_sensitivity_fail(
    Sw92ProfileCPhaseSetSensitivityResult& result,
    PtSensitivityStatus status, std::string diagnostic) {
    result.status = status;
    result.diagnostic = std::move(diagnostic);
    result.phase_fraction_jacobian.clear();
    result.composition_jacobian.clear();
    result.compressibility_jacobian.clear();
    result.molar_density_jacobian.clear();
}

[[nodiscard]] inline bool sw92_profile_c_sensitivity_vector_finite(
    std::span<const double> values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

[[nodiscard]] inline bool sw92_profile_c_sensitivity_same_roundoff(
    double lhs, double rhs) {
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) return false;
    const double guard = 8192.0 * std::numeric_limits<double>::epsilon();
    return std::abs(lhs-rhs) <= guard * std::max(1.0, std::abs(rhs));
}

inline void sw92_profile_c_sensitivity_validate_model(
    const Sw92ProfileCPtPhaseSetResult& source,
    const thermodynamics::Sw92Phase<double>& model) {
    const auto& parameters = model.parameters();
    const std::size_t n = model.size();
    if (n < 2U || source.component_ids.size() != n ||
        source.solution.feed.size() != n ||
        source.dataset_id != parameters.dataset_id() ||
        source.revision != parameters.revision() ||
        std::string_view{source.model_profile} !=
            thermodynamics::sw92_corrected_profile ||
        std::string_view{source.phase_convention} !=
            thermodynamics::sw92_pt_convention ||
        std::string_view{source.equilibrium_profile} !=
            sw92_phase_assigned_aq_na_joint_profile ||
        std::string_view{source.orchestration_convention} !=
            sw92_phase_assigned_pt_convention ||
        std::string_view{source.boundary_convention} !=
            sw92_phase_assigned_boundary_convention ||
        std::string_view{source.publication_convention} !=
            sw92_profile_c_phase_set_publication_convention) {
        throw std::invalid_argument(
            "SW92 sensitivity: phase set and model/provenance snapshot do not match");
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (source.component_ids[i] != parameters.components().items()[i].id)
            throw std::invalid_argument(
                "SW92 sensitivity: ordered component snapshot does not match model");
    }
}

[[nodiscard]] inline double sw92_profile_c_log_distance(
    std::span<const double> a, std::span<const double> b) {
    if (a.size() != b.size()) return 0.0;
    double distance = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (!(a[i] > 0.0) || !(b[i] > 0.0) ||
            !std::isfinite(a[i]) || !std::isfinite(b[i]))
            return 0.0;
        distance = std::max(distance, std::abs(std::log(a[i]/b[i])));
    }
    return distance;
}

} // namespace detail

/// Differentiate one already-authoritative SW92 Profile-C phase set. This is a
/// local derivative of the fixed phase set / fixed AQ-NA family / fixed selected
/// root branches, not a derivative of phase selection, TPD, SSI, RR, topology
/// routing or H-slot morphology. Columns are q=(p,T,z_0,...,z_{N-2}).
[[nodiscard]] inline Sw92ProfileCPhaseSetSensitivityResult
differentiate_sw92_profile_c_phase_set(
    const Sw92ProfileCPtPhaseSetResult& source,
    const thermodynamics::Sw92Phase<double>& model,
    Sw92ProfileCSensitivityOptions options = {}) {
    if (!std::isfinite(options.minimum_derivative_phase_fraction) ||
        !(options.minimum_derivative_phase_fraction > 0.0) ||
        options.minimum_derivative_phase_fraction >= 1.0 ||
        !std::isfinite(options.minimum_log_composition_separation) ||
        !(options.minimum_log_composition_separation > 0.0) ||
        !std::isfinite(options.maximum_base_residual) ||
        !(options.maximum_base_residual > 0.0) ||
        !std::isfinite(options.minimum_reciprocal_condition) ||
        options.minimum_reciprocal_condition < 0.0 ||
        options.minimum_reciprocal_condition >= 1.0 ||
        options.max_components < 2U ||
        options.root_options.max_iterations <= 0) {
        throw std::invalid_argument("SW92 sensitivity: invalid options");
    }

    detail::sw92_profile_c_sensitivity_validate_model(source, model);

    Sw92ProfileCPhaseSetSensitivityResult result;
    result.component_count = model.size();
    result.input_count = model.size() + 1U;
    result.pressure_pa = source.solution.pressure_pa;
    result.temperature_k = source.solution.temperature_k;
    result.feed = source.solution.feed;
    result.nacl_molality_mol_per_kg_water =
        source.nacl_molality_mol_per_kg_water;
    result.dataset_id = source.dataset_id;
    result.revision = source.revision;
    result.component_ids = source.component_ids;
    result.model_profile = source.model_profile;
    result.phase_convention = source.phase_convention;
    result.equilibrium_profile = source.equilibrium_profile;
    result.orchestration_convention = source.orchestration_convention;
    result.boundary_convention = source.boundary_convention;
    result.publication_convention = source.publication_convention;
    result.phase_metadata = source.phase_metadata;

    const std::size_t n = result.component_count;
    if (n > options.max_components) {
        throw std::length_error("SW92 sensitivity: component quota exceeded");
    }
    if (source.solution.status != PtPhaseSetStatus::accepted ||
        !source.accepted_phase_set_published()) {
        detail::sw92_profile_c_sensitivity_fail(
            result, PtSensitivityStatus::solution_not_accepted,
            "SW92 sensitivity: requires an authoritative accepted Profile-C phase set");
        return result;
    }
    const auto* accepted = source.solution.accepted_phase_set();
    if (accepted == nullptr || accepted->phases.empty() ||
        accepted->phases.size() > 3U ||
        accepted->phases.size() != source.phase_metadata.size()) {
        detail::sw92_profile_c_sensitivity_fail(
            result, PtSensitivityStatus::solution_not_accepted,
            "SW92 sensitivity: accepted phase-set shape is inconsistent");
        return result;
    }
    result.phase_count = accepted->phases.size();

    for (double zi : result.feed) {
        if (!std::isfinite(zi) || !(zi > 0.0)) {
            detail::sw92_profile_c_sensitivity_fail(
                result, PtSensitivityStatus::unsupported_feed_support,
                "SW92 sensitivity: v1 reduced-feed derivative requires positive feed support");
            return result;
        }
    }

    double beta_sum = 0.0;
    thermodynamics::Sw92PhaseWorkspace<double> reproduction_workspace;
    try {
        for (std::size_t phase = 0; phase < result.phase_count; ++phase) {
            const auto& state = accepted->phases[phase];
            const auto& metadata = source.phase_metadata[phase];
            if (!std::isfinite(state.mole_phase_fraction) ||
                !(state.mole_phase_fraction >
                  options.minimum_derivative_phase_fraction)) {
                detail::sw92_profile_c_sensitivity_fail(
                    result, PtSensitivityStatus::phase_boundary,
                    "SW92 sensitivity: phase fraction is inside the derivative boundary guard");
                return result;
            }
            beta_sum += state.mole_phase_fraction;
            if (state.composition.size() != n || !state.activity.smooth ||
                state.activity.ln_phi.size() != n) {
                detail::sw92_profile_c_sensitivity_fail(
                    result, PtSensitivityStatus::solution_not_accepted,
                    "SW92 sensitivity: accepted phase property shape is inconsistent");
                return result;
            }
            double composition_sum = 0.0;
            for (double xi : state.composition) {
                if (!std::isfinite(xi) || !(xi > 0.0)) {
                    detail::sw92_profile_c_sensitivity_fail(
                        result, PtSensitivityStatus::unsupported_feed_support,
                        "SW92 sensitivity: positive phase compositions are required");
                    return result;
                }
                composition_sum += xi;
            }
            if (!detail::sw92_profile_c_sensitivity_same_roundoff(
                    composition_sum, 1.0)) {
                detail::sw92_profile_c_sensitivity_fail(
                    result, PtSensitivityStatus::solution_not_accepted,
                    "SW92 sensitivity: accepted phase composition is not normalized");
                return result;
            }

            const auto reproduced = model.evaluate(
                result.pressure_pa, result.temperature_k, state.composition,
                result.nacl_molality_mol_per_kg_water,
                metadata.thermodynamic_family, state.activity.branch,
                reproduction_workspace, options.root_options);
            if (!reproduced.root_set.roots[reproduced.root_index].derivative_valid) {
                detail::sw92_profile_c_sensitivity_fail(
                    result, PtSensitivityStatus::property_failure,
                    "SW92 sensitivity: selected cubic root has no reliable local derivative");
                return result;
            }
            for (std::size_t i = 0; i < n; ++i) {
                if (!detail::sw92_profile_c_sensitivity_same_roundoff(
                        reproduced.ln_phi[i], state.activity.ln_phi[i])) {
                    detail::sw92_profile_c_sensitivity_fail(
                        result, PtSensitivityStatus::solution_not_accepted,
                        "SW92 sensitivity: accepted activity cannot be reproduced");
                    return result;
                }
            }
            if (state.compressibility_factor &&
                !detail::sw92_profile_c_sensitivity_same_roundoff(
                    reproduced.z, *state.compressibility_factor)) {
                detail::sw92_profile_c_sensitivity_fail(
                    result, PtSensitivityStatus::solution_not_accepted,
                    "SW92 sensitivity: accepted Z cannot be reproduced");
                return result;
            }
        }
    } catch (const thermodynamics::Sw92PhaseError& error) {
        detail::sw92_profile_c_sensitivity_fail(
            result, PtSensitivityStatus::property_failure, error.what());
        return result;
    } catch (const std::range_error& error) {
        detail::sw92_profile_c_sensitivity_fail(
            result, PtSensitivityStatus::property_failure, error.what());
        return result;
    } catch (const std::domain_error& error) {
        detail::sw92_profile_c_sensitivity_fail(
            result, PtSensitivityStatus::arithmetic_failure, error.what());
        return result;
    }
    if (!detail::sw92_profile_c_sensitivity_same_roundoff(beta_sum, 1.0)) {
        detail::sw92_profile_c_sensitivity_fail(
            result, PtSensitivityStatus::solution_not_accepted,
            "SW92 sensitivity: accepted phase fractions do not sum to one");
        return result;
    }

    for (std::size_t a = 0; a < result.phase_count; ++a) {
        for (std::size_t b = a + 1U; b < result.phase_count; ++b) {
            if (!(detail::sw92_profile_c_log_distance(
                    accepted->phases[a].composition,
                    accepted->phases[b].composition) >
                  options.minimum_log_composition_separation)) {
                detail::sw92_profile_c_sensitivity_fail(
                    result, PtSensitivityStatus::phase_boundary,
                    "SW92 sensitivity: two phase compositions are inside the local separation guard");
                return result;
            }
        }
    }

    detail::Sw92ProfileCSensitivityLocalSystem local_system;
    try {
        local_system = detail::sw92_profile_c_sensitivity_local_system(
            source, model, options.root_options);
    } catch (const thermodynamics::Sw92PhaseError& error) {
        detail::sw92_profile_c_sensitivity_fail(
            result, PtSensitivityStatus::property_failure, error.what());
        return result;
    } catch (const std::range_error& error) {
        detail::sw92_profile_c_sensitivity_fail(
            result, PtSensitivityStatus::property_failure, error.what());
        return result;
    } catch (const std::domain_error& error) {
        detail::sw92_profile_c_sensitivity_fail(
            result, PtSensitivityStatus::arithmetic_failure, error.what());
        return result;
    } catch (const std::runtime_error& error) {
        detail::sw92_profile_c_sensitivity_fail(
            result, PtSensitivityStatus::arithmetic_failure, error.what());
        return result;
    }

    const auto& local = local_system.values;
    result.equilibrium_residual_norm = 0.0;
    for (std::size_t row = 0; row < local_system.residual_count; ++row) {
        result.equilibrium_residual_norm = std::max(
            result.equilibrium_residual_norm, std::abs(local.values[row]));
    }
    if (!std::isfinite(result.equilibrium_residual_norm) ||
        result.equilibrium_residual_norm > options.maximum_base_residual) {
        detail::sw92_profile_c_sensitivity_fail(
            result, PtSensitivityStatus::solution_not_accepted,
            "SW92 sensitivity: AD local equations do not reproduce the accepted equilibrium");
        return result;
    }

    const std::size_t m = local_system.unknown_count;
    const std::size_t q_count = local_system.input_count;
    const std::size_t variable_count = local_system.variable_count;
    const std::size_t q_offset = m;
    std::vector<double> a(m * m, 0.0);
    std::vector<double> b(m * q_count, 0.0);
    for (std::size_t row = 0; row < m; ++row) {
        for (std::size_t column = 0; column < m; ++column) {
            a[row * m + column] =
                local.jacobian[row * variable_count + column];
        }
        for (std::size_t column = 0; column < q_count; ++column) {
            b[row * q_count + column] =
                local.jacobian[row * variable_count + q_offset + column];
        }
    }

    detail::PtSensitivityLu factor;
    if (!detail::pt_sensitivity_factor(a, m, factor)) {
        detail::sw92_profile_c_sensitivity_fail(
            result, PtSensitivityStatus::ill_conditioned_equilibrium,
            "SW92 sensitivity: dF/du cannot be stably factorized");
        return result;
    }
    result.equilibrium_jacobian_rcond = detail::pt_sensitivity_rcond(factor);
    const double required_rcond = std::max(
        256.0 * std::numeric_limits<double>::epsilon() *
            static_cast<double>(m),
        options.minimum_reciprocal_condition);
    if (!std::isfinite(result.equilibrium_jacobian_rcond) ||
        !(result.equilibrium_jacobian_rcond > required_rcond)) {
        detail::sw92_profile_c_sensitivity_fail(
            result, PtSensitivityStatus::ill_conditioned_equilibrium,
            "SW92 sensitivity: local fixed-phase-set Jacobian is ill-conditioned");
        return result;
    }

    std::vector<double> du_dq(m * q_count, 0.0);
    for (std::size_t column = 0; column < q_count; ++column) {
        std::vector<double> rhs(m, 0.0);
        for (std::size_t row = 0; row < m; ++row)
            rhs[row] = -b[row * q_count + column];
        if (!detail::pt_sensitivity_solve(factor, rhs)) {
            detail::sw92_profile_c_sensitivity_fail(
                result, PtSensitivityStatus::ill_conditioned_equilibrium,
                "SW92 sensitivity: implicit linear solve failed");
            return result;
        }
        for (std::size_t row = 0; row < m; ++row)
            du_dq[row * q_count + column] = rhs[row];
    }

    double backward_error = 0.0;
    for (std::size_t column = 0; column < q_count; ++column) {
        for (std::size_t row = 0; row < m; ++row) {
            double residual = b[row * q_count + column];
            double scale = std::abs(residual);
            for (std::size_t j = 0; j < m; ++j) {
                const double term = a[row * m + j] *
                    du_dq[j * q_count + column];
                residual += term;
                scale += std::abs(term);
            }
            backward_error = std::max(
                backward_error,
                std::abs(residual) / std::max(1.0, scale));
        }
    }
    result.linear_solve_backward_error = backward_error;
    const double arithmetic_scale =
        32768.0 * std::numeric_limits<double>::epsilon() *
        static_cast<double>(std::max<std::size_t>(1U, m));
    if (!std::isfinite(backward_error) || backward_error > arithmetic_scale) {
        detail::sw92_profile_c_sensitivity_fail(
            result, PtSensitivityStatus::arithmetic_failure,
            "SW92 sensitivity: implicit solve failed backward-error check");
        return result;
    }

    const auto total_derivative =
        [&](std::size_t output, std::size_t q_column) {
            double value = local.jacobian[
                output * variable_count + q_offset + q_column];
            for (std::size_t j = 0; j < m; ++j) {
                value += local.jacobian[output * variable_count + j] *
                    du_dq[j * q_count + q_column];
            }
            return value;
        };

    result.phase_fraction_jacobian.assign(
        result.phase_count * q_count, 0.0);
    result.composition_jacobian.assign(
        result.phase_count * n * q_count, 0.0);
    result.compressibility_jacobian.assign(
        result.phase_count * q_count, 0.0);
    result.molar_density_jacobian.assign(
        result.phase_count * q_count, 0.0);
    for (std::size_t phase = 0; phase < result.phase_count; ++phase) {
        for (std::size_t column = 0; column < q_count; ++column) {
            result.phase_fraction_jacobian[phase*q_count+column] =
                total_derivative(
                    local_system.phase_fraction_output + phase, column);
            result.compressibility_jacobian[phase*q_count+column] =
                total_derivative(
                    local_system.compressibility_output + phase, column);
            result.molar_density_jacobian[phase*q_count+column] =
                total_derivative(local_system.density_output + phase, column);
            for (std::size_t i = 0; i < n; ++i) {
                result.composition_jacobian[
                    (phase*n+i)*q_count+column] =
                    total_derivative(
                        local_system.composition_output + phase*n + i,
                        column);
            }
        }
    }

    if (!detail::sw92_profile_c_sensitivity_vector_finite(
            result.phase_fraction_jacobian) ||
        !detail::sw92_profile_c_sensitivity_vector_finite(
            result.composition_jacobian) ||
        !detail::sw92_profile_c_sensitivity_vector_finite(
            result.compressibility_jacobian) ||
        !detail::sw92_profile_c_sensitivity_vector_finite(
            result.molar_density_jacobian)) {
        detail::sw92_profile_c_sensitivity_fail(
            result, PtSensitivityStatus::arithmetic_failure,
            "SW92 sensitivity: nonfinite total derivative");
        return result;
    }

    // Differentiate normalization and material balance as independent algebraic
    // checks on the assembled total derivatives.
    const double identity_guard = std::max(
        options.maximum_base_residual,
        16.0 * arithmetic_scale);
    for (std::size_t column = 0; column < q_count; ++column) {
        double beta_derivative_sum = 0.0;
        for (std::size_t phase = 0; phase < result.phase_count; ++phase) {
            beta_derivative_sum += result.d_phase_fraction(phase, column);
            double composition_sum = 0.0;
            for (std::size_t i = 0; i < n; ++i)
                composition_sum += result.d_composition(phase, i, column);
            if (std::abs(composition_sum) > identity_guard) {
                detail::sw92_profile_c_sensitivity_fail(
                    result, PtSensitivityStatus::arithmetic_failure,
                    "SW92 sensitivity: differentiated phase normalization failed");
                return result;
            }
        }
        if (std::abs(beta_derivative_sum) > identity_guard) {
            detail::sw92_profile_c_sensitivity_fail(
                result, PtSensitivityStatus::arithmetic_failure,
                "SW92 sensitivity: differentiated phase-fraction closure failed");
            return result;
        }

        for (std::size_t i = 0; i < n; ++i) {
            double recovered = 0.0;
            for (std::size_t phase = 0; phase < result.phase_count; ++phase) {
                const auto& state = accepted->phases[phase];
                recovered +=
                    result.d_phase_fraction(phase,column)*state.composition[i] +
                    state.mole_phase_fraction*
                        result.d_composition(phase,i,column);
            }
            double expected = 0.0;
            if (column >= 2U) {
                const std::size_t independent = column - 2U;
                if (i == independent) expected = 1.0;
                else if (i == n - 1U) expected = -1.0;
            }
            if (std::abs(recovered-expected) > identity_guard) {
                detail::sw92_profile_c_sensitivity_fail(
                    result, PtSensitivityStatus::arithmetic_failure,
                    "SW92 sensitivity: differentiated material balance failed");
                return result;
            }
        }
    }

    result.status = PtSensitivityStatus::success;
    result.diagnostic =
        "SW92 sensitivity: fixed authoritative phase-set implicit Jacobian available";
    return result;
}

} // namespace mpmc::flash

#endif // MPMC_FLASH_SW92_PROFILE_C_SENSITIVITY_HPP
