#ifndef MPMC_PHYSICS_SW92_THERMODYNAMIC_CLOSURE_HPP
#define MPMC_PHYSICS_SW92_THERMODYNAMIC_CLOSURE_HPP

#include <mpmc/flash/sw92_profile_c_phase_set.hpp>
#include <mpmc/flash/sw92_profile_c_sensitivity.hpp>
#include <mpmc/physics/thermodynamic_closure.hpp>
#include <mpmc/thermodynamics/sw92_phase.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::physics {

inline constexpr std::string_view sw92_profile_c_closure_convention =
    "SW92/phase-assigned/Profile-C/physics-thermodynamic-closure-reduced-feed/v2";

struct Sw92ProfileCClosureOptions {
    // One root policy is used by both primal reproduction and the derivative
    // path so a closure snapshot cannot mix two selected-root conventions.
    mpmc::thermodynamics::Sw92RootOptions root_options{};
    mpmc::flash::Sw92ProfileCSensitivityOptions sensitivity{};
};

// Model-specific phase identity remains a sidecar to the generic physics phase
// state. `nonaqueous_unclassified` is intentionally not promoted to liquid or
// vapor. The root branch is numerical provenance, not physical morphology.
struct Sw92ProfileCClosurePhaseMetadata {
    mpmc::flash::Sw92PhaseAssignedPtPhysicalRole physical_role{
        mpmc::flash::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified};
    mpmc::thermodynamics::SwPhaseFamily thermodynamic_family{
        mpmc::thermodynamics::SwPhaseFamily::nonaqueous};
    std::size_t root_branch{};
    bool source_compressibility_factor_present{false};
};

struct Sw92ProfileCThermodynamicClosureSnapshot {
    static constexpr bool global_stability_proven = false;
    static constexpr bool morphology_resolved = false;

    PtPhaseSetThermodynamicClosureSnapshot closure;
    double input_feed_sum{};
    double nacl_molality_mol_per_kg_water{};
    std::string model_profile;
    std::string phase_convention;
    std::string equilibrium_profile;
    std::string orchestration_convention;
    std::string boundary_convention;
    std::string publication_convention;
    std::string sensitivity_convention;
    std::string closure_convention{sw92_profile_c_closure_convention};
    std::vector<Sw92ProfileCClosurePhaseMetadata> phase_metadata;

    [[nodiscard]] bool residual_available() const noexcept {
        return closure.residual_available() && closure.primal &&
               phase_metadata.size() == closure.primal->phases.size();
    }
    [[nodiscard]] bool can_seed_newton() const noexcept {
        return residual_available() && closure.can_seed_newton();
    }
};

namespace detail {

[[nodiscard]] inline bool sw92_closure_same_roundoff(double lhs, double rhs) {
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) { return false; }
    const double guard = 4096.0 * std::numeric_limits<double>::epsilon();
    return std::abs(lhs - rhs) <= guard * std::max(1.0, std::abs(rhs));
}

[[nodiscard]] inline bool sw92_closure_vector_same_roundoff(
    std::span<const double> lhs, std::span<const double> rhs) {
    if (lhs.size() != rhs.size()) { return false; }
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        if (!sw92_closure_same_roundoff(lhs[i], rhs[i])) { return false; }
    }
    return true;
}

[[nodiscard]] inline bool sw92_closure_vector_finite(
    std::span<const double> values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

[[nodiscard]] inline bool sw92_closure_simplex_valid(
    std::span<const double> values) {
    if (values.empty()) { return false; }
    double sum = 0.0;
    for (double value : values) {
        if (!std::isfinite(value) || value < 0.0) { return false; }
        sum += value;
    }
    if (!std::isfinite(sum)) { return false; }
    const double guard = 8192.0 * std::numeric_limits<double>::epsilon() *
                         static_cast<double>(std::max<std::size_t>(1U, values.size()));
    return std::abs(sum - 1.0) <= guard;
}

inline void sw92_closure_validate_model_snapshot(
    const mpmc::flash::Sw92ProfileCPtPhaseSetResult& source,
    const mpmc::thermodynamics::Sw92Phase<double>& model) {
    const auto& parameters = model.parameters();
    const std::size_t n = model.size();
    if (n == 0U || source.component_ids.size() != n ||
        source.solution.feed.size() != n ||
        source.dataset_id != parameters.dataset_id() ||
        source.revision != parameters.revision()) {
        throw std::invalid_argument(
            "SW92 closure: authoritative phase set and model snapshot do not match");
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (source.component_ids[i] != parameters.components().items()[i].id) {
            throw std::invalid_argument(
                "SW92 closure: ordered component snapshot does not match model");
        }
    }
}

[[nodiscard]] inline bool sw92_closure_publication_contract_valid(
    const mpmc::flash::Sw92ProfileCPtPhaseSetResult& source) {
    namespace fl = mpmc::flash;
    namespace th = mpmc::thermodynamics;
    return std::string_view{source.model_profile} == th::sw92_corrected_profile &&
           std::string_view{source.phase_convention} == th::sw92_pt_convention &&
           std::string_view{source.equilibrium_profile} ==
               fl::sw92_phase_assigned_aq_na_joint_profile &&
           std::string_view{source.orchestration_convention} ==
               fl::sw92_phase_assigned_pt_convention &&
           std::string_view{source.boundary_convention} ==
               fl::sw92_phase_assigned_boundary_convention &&
           std::string_view{source.publication_convention} ==
               fl::sw92_profile_c_phase_set_publication_convention &&
           source.solution.capability.maximum_phase_count == 3U &&
           !source.solution.global_stability_proven &&
           std::isfinite(source.solution.pressure_pa) &&
           source.solution.pressure_pa > 0.0 &&
           std::isfinite(source.solution.temperature_k) &&
           source.solution.temperature_k > 0.0 &&
           std::isfinite(source.nacl_molality_mol_per_kg_water) &&
           source.nacl_molality_mol_per_kg_water >= 0.0 &&
           sw92_closure_simplex_valid(source.solution.feed);
}

[[nodiscard]] inline bool sw92_closure_role_family_valid(
    const mpmc::flash::Sw92ProfileCPhaseMetadata& metadata) {
    namespace fl = mpmc::flash;
    namespace th = mpmc::thermodynamics;
    switch (metadata.physical_role) {
    case fl::Sw92PhaseAssignedPtPhysicalRole::aqueous:
        return metadata.thermodynamic_family == th::SwPhaseFamily::aqueous;
    case fl::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified:
        return metadata.thermodynamic_family == th::SwPhaseFamily::nonaqueous;
    }
    return false;
}

[[nodiscard]] inline ThermodynamicClosureLinearizationReason
sw92_closure_map_sensitivity_failure(mpmc::flash::PtSensitivityStatus status) {
    using Source = mpmc::flash::PtSensitivityStatus;
    switch (status) {
    case Source::success:
        return ThermodynamicClosureLinearizationReason::none;
    case Source::phase_boundary:
        return ThermodynamicClosureLinearizationReason::phase_boundary;
    case Source::ill_conditioned_equilibrium:
        return ThermodynamicClosureLinearizationReason::ill_conditioned_equilibrium;
    case Source::unsupported_feed_support:
        return ThermodynamicClosureLinearizationReason::unsupported_feed_support;
    case Source::property_failure:
        return ThermodynamicClosureLinearizationReason::property_failure;
    case Source::arithmetic_failure:
        return ThermodynamicClosureLinearizationReason::arithmetic_failure;
    case Source::solution_not_accepted:
        return ThermodynamicClosureLinearizationReason::solution_not_accepted;
    }
    return ThermodynamicClosureLinearizationReason::arithmetic_failure;
}

[[nodiscard]] inline bool sw92_closure_sensitivity_matches_snapshot(
    const mpmc::flash::Sw92ProfileCPhaseSetSensitivityResult& sensitivity,
    const mpmc::flash::Sw92ProfileCPtPhaseSetResult& source,
    const PtPhaseSetThermodynamicState& primal) {
    const std::size_t n = source.component_ids.size();
    const std::size_t phase_count = primal.phases.size();
    if (sensitivity.status != mpmc::flash::PtSensitivityStatus::success ||
        sensitivity.component_count != n ||
        sensitivity.phase_count != phase_count ||
        sensitivity.input_count != n + 1U ||
        sensitivity.pressure_pa != source.solution.pressure_pa ||
        sensitivity.temperature_k != source.solution.temperature_k ||
        sensitivity.feed != source.solution.feed ||
        sensitivity.nacl_molality_mol_per_kg_water !=
            source.nacl_molality_mol_per_kg_water ||
        sensitivity.dataset_id != source.dataset_id ||
        sensitivity.revision != source.revision ||
        sensitivity.component_ids != source.component_ids ||
        sensitivity.model_profile != source.model_profile ||
        sensitivity.phase_convention != source.phase_convention ||
        sensitivity.equilibrium_profile != source.equilibrium_profile ||
        sensitivity.orchestration_convention != source.orchestration_convention ||
        sensitivity.boundary_convention != source.boundary_convention ||
        sensitivity.publication_convention != source.publication_convention ||
        sensitivity.phase_metadata.size() != source.phase_metadata.size() ||
        sensitivity.phase_fraction_jacobian.size() !=
            phase_count * sensitivity.input_count ||
        sensitivity.composition_jacobian.size() !=
            phase_count * n * sensitivity.input_count ||
        sensitivity.compressibility_jacobian.size() !=
            phase_count * sensitivity.input_count ||
        sensitivity.molar_density_jacobian.size() !=
            phase_count * sensitivity.input_count ||
        !sw92_closure_vector_finite(sensitivity.phase_fraction_jacobian) ||
        !sw92_closure_vector_finite(sensitivity.composition_jacobian) ||
        !sw92_closure_vector_finite(sensitivity.compressibility_jacobian) ||
        !sw92_closure_vector_finite(sensitivity.molar_density_jacobian) ||
        !std::isfinite(sensitivity.equilibrium_jacobian_rcond) ||
        !(sensitivity.equilibrium_jacobian_rcond > 0.0) ||
        !std::isfinite(sensitivity.linear_solve_backward_error) ||
        sensitivity.linear_solve_backward_error < 0.0 ||
        !std::isfinite(sensitivity.equilibrium_residual_norm) ||
        sensitivity.equilibrium_residual_norm < 0.0) {
        return false;
    }
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        if (sensitivity.phase_metadata[phase].physical_role !=
                source.phase_metadata[phase].physical_role ||
            sensitivity.phase_metadata[phase].thermodynamic_family !=
                source.phase_metadata[phase].thermodynamic_family) {
            return false;
        }
    }
    return true;
}

} // namespace detail

// Consume an already authoritative SW92 Profile-C PT phase set. This function
// performs no flash, TPD, material-balance solve, topology search, composition
// normalization or morphology classification. It first reproduces the exact
// accepted primal, then atomically attaches the validated fixed-phase-set local
// sensitivity when available. Derivative-only failures preserve the primal;
// `solution_not_accepted` invalidates the whole snapshot as an internal
// consistency failure.
[[nodiscard]] inline Sw92ProfileCThermodynamicClosureSnapshot
build_sw92_profile_c_thermodynamic_closure(
    const mpmc::flash::Sw92ProfileCPtPhaseSetResult& source,
    const mpmc::thermodynamics::Sw92Phase<double>& model,
    Sw92ProfileCClosureOptions options = {}) {
    namespace fl = mpmc::flash;
    namespace th = mpmc::thermodynamics;

    detail::sw92_closure_validate_model_snapshot(source, model);

    Sw92ProfileCThermodynamicClosureSnapshot result;
    result.closure.pressure_pa = source.solution.pressure_pa;
    result.closure.temperature_k = source.solution.temperature_k;
    result.closure.feed = source.solution.feed;
    result.closure.component_ids = source.component_ids;
    result.closure.thermodynamic_model = source.model_profile;
    result.closure.dataset_id = source.dataset_id;
    result.closure.revision = source.revision;
    result.input_feed_sum = source.input_feed_sum;
    result.nacl_molality_mol_per_kg_water =
        source.nacl_molality_mol_per_kg_water;
    result.model_profile = source.model_profile;
    result.phase_convention = source.phase_convention;
    result.equilibrium_profile = source.equilibrium_profile;
    result.orchestration_convention = source.orchestration_convention;
    result.boundary_convention = source.boundary_convention;
    result.publication_convention = source.publication_convention;
    result.sensitivity_convention =
        std::string(fl::sw92_profile_c_sensitivity_convention);

    if (!detail::sw92_closure_publication_contract_valid(source)) {
        result.closure.primal_status = ThermodynamicClosurePrimalStatus::indeterminate;
        result.closure.diagnostic =
            "SW92 closure: authoritative publication metadata is inconsistent";
        return result;
    }
    if (source.solution.status != fl::PtPhaseSetStatus::accepted ||
        !source.accepted_phase_set_published()) {
        result.closure.primal_status = ThermodynamicClosurePrimalStatus::indeterminate;
        result.closure.diagnostic =
            "SW92 closure: phase set is not an accepted authoritative primal";
        return result;
    }

    const auto* accepted = source.solution.accepted_phase_set();
    if (accepted == nullptr || accepted->phases.empty() ||
        accepted->phases.size() > 3U ||
        accepted->phases.size() != source.phase_metadata.size()) {
        result.closure.primal_status = ThermodynamicClosurePrimalStatus::indeterminate;
        result.closure.diagnostic =
            "SW92 closure: accepted phase-set shape/metadata is inconsistent";
        return result;
    }

    double phase_fraction_sum = 0.0;
    PtPhaseSetThermodynamicState primal;
    primal.phases.reserve(accepted->phases.size());
    result.phase_metadata.reserve(accepted->phases.size());

    try {
        th::Sw92PhaseWorkspace<double> workspace;
        const double r = th::Sw92Pure<double>::gas_constant();
        for (std::size_t phase_index = 0; phase_index < accepted->phases.size();
             ++phase_index) {
            const auto& phase = accepted->phases[phase_index];
            const auto& metadata = source.phase_metadata[phase_index];
            if (!detail::sw92_closure_role_family_valid(metadata) ||
                !std::isfinite(phase.mole_phase_fraction) ||
                !(phase.mole_phase_fraction > 0.0) ||
                phase.composition.size() != model.size() ||
                !detail::sw92_closure_simplex_valid(phase.composition) ||
                !phase.activity.smooth ||
                phase.activity.ln_phi.size() != model.size()) {
                throw std::runtime_error(
                    "invalid accepted phase shape, activity or role/family metadata");
            }
            phase_fraction_sum += phase.mole_phase_fraction;

            const auto reproduced = model.evaluate(
                result.closure.pressure_pa, result.closure.temperature_k,
                std::span<const double>{phase.composition.data(), phase.composition.size()},
                result.nacl_molality_mol_per_kg_water,
                metadata.thermodynamic_family, phase.activity.branch,
                workspace, options.root_options);
            if (reproduced.root_index != phase.activity.branch ||
                reproduced.family != metadata.thermodynamic_family ||
                !detail::sw92_closure_vector_same_roundoff(
                    reproduced.ln_phi, phase.activity.ln_phi)) {
                throw std::runtime_error(
                    "accepted activity cannot be reproduced from SW92 family/root provenance");
            }

            double z = reproduced.z;
            const bool source_z_present = phase.compressibility_factor.has_value();
            if (source_z_present) {
                if (!detail::sw92_closure_same_roundoff(
                        reproduced.z, *phase.compressibility_factor)) {
                    throw std::runtime_error(
                        "accepted Z cannot be reproduced from SW92 family/root provenance");
                }
                // Preserve the exact accepted flash value after reproduction.
                z = *phase.compressibility_factor;
            }
            if (!std::isfinite(z) || !(z > 0.0)) {
                throw std::range_error("nonpositive/nonfinite accepted SW92 Z");
            }

            const double density = result.closure.pressure_pa /
                (z * r * result.closure.temperature_k);
            if (!std::isfinite(density) || !(density > 0.0)) {
                throw std::range_error("nonrepresentable SW92 phase molar density");
            }

            ThermodynamicPhaseState state;
            state.mole_phase_fraction = phase.mole_phase_fraction;
            state.composition = phase.composition;
            state.compressibility_factor = z;
            state.molar_density_mol_per_m3 = density;
            primal.phases.push_back(std::move(state));
            result.phase_metadata.push_back(
                {metadata.physical_role, metadata.thermodynamic_family,
                 phase.activity.branch, source_z_present});
        }

        const double fraction_guard =
            8192.0 * std::numeric_limits<double>::epsilon() *
            static_cast<double>(std::max<std::size_t>(1U, accepted->phases.size()));
        if (!std::isfinite(phase_fraction_sum) ||
            std::abs(phase_fraction_sum - 1.0) > fraction_guard) {
            throw std::runtime_error("accepted phase fractions do not close to unity");
        }

        result.closure.primal = std::move(primal);
        result.closure.primal_status = ThermodynamicClosurePrimalStatus::valid;
    } catch (const std::exception& error) {
        result.closure.primal_status = ThermodynamicClosurePrimalStatus::indeterminate;
        result.closure.primal.reset();
        result.phase_metadata.clear();
        result.closure.diagnostic =
            std::string("SW92 closure: primal reproduction failed: ") + error.what();
        return result;
    }

    auto sensitivity_options = options.sensitivity;
    sensitivity_options.root_options = options.root_options;
    const auto sensitivity = fl::differentiate_sw92_profile_c_phase_set(
        source, model, sensitivity_options);
    if (sensitivity.status != fl::PtSensitivityStatus::success) {
        result.closure.linearization_status =
            ThermodynamicClosureLinearizationStatus::unavailable;
        result.closure.linearization_reason =
            detail::sw92_closure_map_sensitivity_failure(sensitivity.status);
        result.closure.linearization.reset();
        if (sensitivity.status == fl::PtSensitivityStatus::solution_not_accepted) {
            result.closure.primal_status =
                ThermodynamicClosurePrimalStatus::indeterminate;
            result.closure.primal.reset();
            result.phase_metadata.clear();
            result.closure.diagnostic =
                std::string("SW92 closure: sensitivity rejected the accepted primal: ") +
                sensitivity.diagnostic;
            return result;
        }
        result.closure.diagnostic =
            std::string("SW92 closure: authoritative Profile-C primal available; ") +
            "local linearization unavailable: " + sensitivity.diagnostic;
        return result;
    }

    if (!result.closure.primal ||
        !detail::sw92_closure_sensitivity_matches_snapshot(
            sensitivity, source, *result.closure.primal)) {
        result.closure.primal_status = ThermodynamicClosurePrimalStatus::indeterminate;
        result.closure.primal.reset();
        result.closure.linearization.reset();
        result.phase_metadata.clear();
        result.closure.linearization_status =
            ThermodynamicClosureLinearizationStatus::unavailable;
        result.closure.linearization_reason =
            ThermodynamicClosureLinearizationReason::solution_not_accepted;
        result.closure.diagnostic =
            "SW92 closure: successful sensitivity does not match the owned primal snapshot";
        return result;
    }

    PtPhaseSetThermodynamicLinearization linearization;
    linearization.component_count = sensitivity.component_count;
    linearization.phase_count = sensitivity.phase_count;
    linearization.input_count = sensitivity.input_count;
    linearization.phase_fraction_jacobian = sensitivity.phase_fraction_jacobian;
    linearization.composition_jacobian = sensitivity.composition_jacobian;
    linearization.compressibility_jacobian = sensitivity.compressibility_jacobian;
    linearization.molar_density_jacobian = sensitivity.molar_density_jacobian;
    linearization.equilibrium_jacobian_rcond =
        sensitivity.equilibrium_jacobian_rcond;
    linearization.linear_solve_backward_error =
        sensitivity.linear_solve_backward_error;
    linearization.equilibrium_residual_norm =
        sensitivity.equilibrium_residual_norm;

    result.closure.linearization = std::move(linearization);
    result.closure.linearization_status =
        ThermodynamicClosureLinearizationStatus::available;
    result.closure.linearization_reason =
        ThermodynamicClosureLinearizationReason::none;
    result.closure.diagnostic =
        "SW92 closure: authoritative Profile-C primal and fixed-phase-set "
        "local linearization available";
    return result;
}

} // namespace mpmc::physics

#endif // MPMC_PHYSICS_SW92_THERMODYNAMIC_CLOSURE_HPP
