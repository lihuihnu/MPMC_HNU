#include <mpmc/physics/sw92_thermodynamic_closure.hpp>

#include "../../flash/sw92_phase_assigned_pt/physical_sample6.hpp"
#include "../../thermodynamics/sw92/test_support.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool sw92_thermodynamic_closure_header();

namespace {
namespace fl = mpmc::flash;
namespace ph = mpmc::physics;
namespace th = mpmc::thermodynamics;
namespace sample6 = sw92_profile_c_sample6;
using Vec = std::vector<double>;

void require(bool condition, std::string_view message,
             std::source_location where = std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(std::string(where.file_name()) + ":" +
                                 std::to_string(where.line()) + ": " +
                                 std::string(message));
    }
}

void near(double actual, double expected, double relative = 2e-11,
          double absolute = 2e-12,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) >
            absolute + relative * std::max(1.0, std::abs(expected))) {
        std::cerr << "actual=" << actual << " expected=" << expected << '\n';
        require(false, "numeric mismatch", where);
    }
}

th::Sw92Phase<double> binary_model(bool reverse = false) {
    return th::Sw92Phase<double>::from_parameters(
        sw92_test::binary_parameters(sw92_test::co2, reverse));
}

void require_primal_basics(
    const ph::Sw92ProfileCThermodynamicClosureSnapshot& closure,
    std::size_t phase_count) {
    require(closure.residual_available() && closure.closure.primal &&
                closure.closure.primal->phases.size() == phase_count &&
                closure.phase_metadata.size() == phase_count,
            "SW92 physics closure did not publish expected primal phase count");
    require(closure.closure.primal_status ==
                ph::ThermodynamicClosurePrimalStatus::valid,
            "SW92 closure did not retain a valid primal");
    require(!closure.global_stability_proven && !closure.morphology_resolved &&
                closure.closure_convention == ph::sw92_profile_c_closure_convention &&
                closure.sensitivity_convention ==
                    fl::sw92_profile_c_sensitivity_convention,
            "SW92 closure changed global-stability/morphology/sensitivity provenance");
}

void require_linearization_basics(
    const ph::Sw92ProfileCThermodynamicClosureSnapshot& closure,
    std::size_t phase_count) {
    require_primal_basics(closure, phase_count);
    require(closure.closure.linearization_status ==
                ph::ThermodynamicClosureLinearizationStatus::available &&
                closure.closure.linearization_reason ==
                    ph::ThermodynamicClosureLinearizationReason::none &&
                closure.closure.linearization.has_value() &&
                closure.can_seed_newton(),
            "SW92 closure did not publish the successful local linearization atomically");
    const auto& linearization = *closure.closure.linearization;
    const std::size_t n = closure.closure.component_ids.size();
    require(linearization.component_count == n &&
                linearization.phase_count == phase_count &&
                linearization.input_count == n + 1U &&
                linearization.phase_fraction_jacobian.size() ==
                    phase_count * linearization.input_count &&
                linearization.composition_jacobian.size() ==
                    phase_count * n * linearization.input_count &&
                linearization.compressibility_jacobian.size() ==
                    phase_count * linearization.input_count &&
                linearization.molar_density_jacobian.size() ==
                    phase_count * linearization.input_count &&
                std::isfinite(linearization.equilibrium_jacobian_rcond) &&
                linearization.equilibrium_jacobian_rcond > 0.0 &&
                std::isfinite(linearization.linear_solve_backward_error) &&
                std::isfinite(linearization.equilibrium_residual_norm),
            "SW92 closure linearization shape/diagnostics are invalid");
}

void require_linearization_matches_sensitivity(
    const fl::Sw92ProfileCPtPhaseSetResult& source,
    const th::Sw92Phase<double>& model,
    const ph::Sw92ProfileCThermodynamicClosureSnapshot& closure) {
    const auto sensitivity =
        fl::differentiate_sw92_profile_c_phase_set(source, model);
    require(sensitivity.status == fl::PtSensitivityStatus::success,
            sensitivity.diagnostic);
    require(closure.closure.linearization.has_value(),
            "matching check requires a physics linearization");
    const auto& linearization = *closure.closure.linearization;
    require(linearization.component_count == sensitivity.component_count &&
                linearization.phase_count == sensitivity.phase_count &&
                linearization.input_count == sensitivity.input_count &&
                linearization.phase_fraction_jacobian ==
                    sensitivity.phase_fraction_jacobian &&
                linearization.composition_jacobian ==
                    sensitivity.composition_jacobian &&
                linearization.compressibility_jacobian ==
                    sensitivity.compressibility_jacobian &&
                linearization.molar_density_jacobian ==
                    sensitivity.molar_density_jacobian,
            "physics closure changed the validated SW92 flash Jacobian");
    near(linearization.equilibrium_jacobian_rcond,
         sensitivity.equilibrium_jacobian_rcond, 0.0, 0.0);
    near(linearization.linear_solve_backward_error,
         sensitivity.linear_solve_backward_error, 0.0, 0.0);
    near(linearization.equilibrium_residual_norm,
         sensitivity.equilibrium_residual_norm, 0.0, 0.0);
}

void require_density_identity(const ph::PtPhaseSetThermodynamicClosureSnapshot& closure) {
    require(closure.primal.has_value(), "density identity requires a primal");
    const double r = th::Sw92Pure<double>::gas_constant();
    for (const auto& phase : closure.primal->phases) {
        const double expected = closure.pressure_pa /
            (phase.compressibility_factor * r * closure.temperature_k);
        near(phase.molar_density_mol_per_m3, expected, 2e-12, 2e-10);
    }
}

void one_phase_reproduces_missing_z() {
    const auto model = binary_model();
    const auto source = fl::solve_sw92_profile_c_pt_phase_set(
        3.0e6, 340.0, Vec{0.995, 0.005}, model, 0.0);
    require(source.accepted_phase_set_published() &&
                source.solution.accepted_phase_count() == 1U,
            "dry binary no longer supplies authoritative one-phase source");
    const auto& published = source.solution.accepted_phase_set()->phases.front();
    require(!published.compressibility_factor,
            "one-phase publication unexpectedly started storing Z");

    const auto closure = ph::build_sw92_profile_c_thermodynamic_closure(source, model);
    require_linearization_basics(closure, 1U);
    require_linearization_matches_sensitivity(source, model, closure);
    const auto& phase = closure.closure.primal->phases.front();
    require(phase.composition == published.composition &&
                phase.mole_phase_fraction == published.mole_phase_fraction &&
                std::isfinite(phase.compressibility_factor) &&
                phase.compressibility_factor > 0.0,
            "one-phase closure did not reproduce accepted H phase property");
    require(closure.phase_metadata.front().physical_role ==
                fl::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified &&
                closure.phase_metadata.front().thermodynamic_family ==
                    th::SwPhaseFamily::nonaqueous &&
                !closure.phase_metadata.front().source_compressibility_factor_present,
            "one-phase closure invented morphology or source-Z provenance");
    require_density_identity(closure.closure);
}

void two_phase_preserves_authoritative_properties() {
    const auto model = binary_model();
    const auto source = fl::solve_sw92_profile_c_pt_phase_set(
        3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0);
    require(source.accepted_phase_set_published() &&
                source.solution.accepted_phase_count() == 2U,
            "wet binary no longer supplies authoritative W+H source");

    const auto closure = ph::build_sw92_profile_c_thermodynamic_closure(source, model);
    require_linearization_basics(closure, 2U);
    require_linearization_matches_sensitivity(source, model, closure);
    const auto& published = source.solution.accepted_phase_set()->phases;
    const auto& phases = closure.closure.primal->phases;
    for (std::size_t i = 0; i < 2U; ++i) {
        require(published[i].compressibility_factor.has_value() &&
                    closure.phase_metadata[i].source_compressibility_factor_present,
                "W+H source/closure lost Z provenance");
        require(phases[i].mole_phase_fraction == published[i].mole_phase_fraction &&
                    phases[i].composition == published[i].composition &&
                    phases[i].compressibility_factor == *published[i].compressibility_factor,
                "W+H closure changed accepted flash primal values");
    }
    require(closure.phase_metadata[0].physical_role ==
                fl::Sw92PhaseAssignedPtPhysicalRole::aqueous &&
                closure.phase_metadata[0].thermodynamic_family ==
                    th::SwPhaseFamily::aqueous &&
                closure.phase_metadata[1].physical_role ==
                    fl::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified &&
                closure.phase_metadata[1].thermodynamic_family ==
                    th::SwPhaseFamily::nonaqueous,
            "W+H closure changed AQ/NA role-family mapping");
    require_density_identity(closure.closure);
}

void sample6_three_phase_primal() {
    const auto model = sample6::model();
    const auto source = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7, 350.0, sample6::feed(), model, 0.0);
    require(source.accepted_phase_set_published() &&
                source.solution.accepted_phase_count() == 3U,
            "Sample-6 no longer supplies authoritative three-phase source");

    const auto closure = ph::build_sw92_profile_c_thermodynamic_closure(source, model);
    require_linearization_basics(closure, 3U);
    require_linearization_matches_sensitivity(source, model, closure);
    require(closure.closure.dataset_id == source.dataset_id &&
                closure.closure.revision == source.revision &&
                closure.closure.component_ids == source.component_ids &&
                closure.nacl_molality_mol_per_kg_water == 0.0,
            "Sample-6 closure lost model/data provenance");

    const auto& source_phases = source.solution.accepted_phase_set()->phases;
    const auto& phases = closure.closure.primal->phases;
    for (std::size_t i = 0; i < 3U; ++i) {
        require(source_phases[i].compressibility_factor.has_value(),
                "Sample-6 authoritative source lost Z");
        near(phases[i].mole_phase_fraction,
             static_cast<double>(sample6::golden.fractions[i]), 8e-8, 3e-12);
        near(phases[i].compressibility_factor,
             static_cast<double>(sample6::golden.z[i]), 8e-8, 3e-12);
        require(phases[i].composition == source_phases[i].composition,
                "Sample-6 closure changed accepted phase composition");
    }
    require(closure.phase_metadata[0].physical_role ==
                fl::Sw92PhaseAssignedPtPhysicalRole::aqueous,
            "Sample-6 W role changed");
    for (std::size_t i = 1U; i < 3U; ++i) {
        require(closure.phase_metadata[i].physical_role ==
                    fl::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified &&
                    closure.phase_metadata[i].thermodynamic_family ==
                        th::SwPhaseFamily::nonaqueous,
                "Sample-6 H phase was relabelled as L/V");
    }
    require_density_identity(closure.closure);
}

void linearization_failure_policy() {
    const auto model = sample6::model();
    const auto source = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7, 350.0, sample6::feed(), model, 0.0);
    require(source.accepted_phase_set_published() &&
                source.solution.accepted_phase_count() == 3U,
            "linearization failure fixture lost Sample-6 phase set");

    ph::Sw92ProfileCClosureOptions options;
    options.sensitivity.minimum_derivative_phase_fraction = 0.04;
    const auto guarded = ph::build_sw92_profile_c_thermodynamic_closure(
        source, model, options);
    require_primal_basics(guarded, 3U);
    require(guarded.closure.linearization_status ==
                ph::ThermodynamicClosureLinearizationStatus::unavailable &&
                guarded.closure.linearization_reason ==
                    ph::ThermodynamicClosureLinearizationReason::phase_boundary &&
                !guarded.closure.linearization.has_value() &&
                !guarded.can_seed_newton(),
            "derivative boundary failure did not preserve primal / reject Newton atomically");

    auto inconsistent = source;
    require(inconsistent.solution.feed.size() >= 2U,
            "inconsistency fixture lost feed dimensions");
    inconsistent.solution.feed[0] += 1.0e-4;
    inconsistent.solution.feed[1] -= 1.0e-4;
    const auto rejected = ph::build_sw92_profile_c_thermodynamic_closure(
        inconsistent, model);
    require(!rejected.residual_available() &&
                rejected.closure.primal_status ==
                    ph::ThermodynamicClosurePrimalStatus::indeterminate &&
                rejected.closure.linearization_status ==
                    ph::ThermodynamicClosureLinearizationStatus::unavailable &&
                rejected.closure.linearization_reason ==
                    ph::ThermodynamicClosureLinearizationReason::solution_not_accepted &&
                !rejected.closure.primal.has_value() &&
                !rejected.closure.linearization.has_value() &&
                rejected.phase_metadata.empty() &&
                !rejected.can_seed_newton(),
            "sensitivity/primal inconsistency did not invalidate the atomic snapshot");
}

void rejection_and_model_guard() {
    const auto model = binary_model();
    const auto source = fl::solve_sw92_profile_c_pt_phase_set(
        3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0);
    require(source.accepted_phase_set_published(),
            "rejection fixture lacks accepted source");

    auto convention_tamper = source;
    convention_tamper.publication_convention = "tampered/publication";
    const auto rejected_convention =
        ph::build_sw92_profile_c_thermodynamic_closure(convention_tamper, model);
    require(!rejected_convention.residual_available() &&
                rejected_convention.closure.primal_status ==
                    ph::ThermodynamicClosurePrimalStatus::indeterminate &&
                !rejected_convention.closure.linearization.has_value(),
            "closure accepted tampered authoritative publication identity");

    auto activity_tamper = source;
    require(activity_tamper.solution.candidate_phase_set &&
                !activity_tamper.solution.candidate_phase_set->phases.empty(),
            "activity tamper fixture lost phase set");
    activity_tamper.solution.candidate_phase_set->phases[0].activity.ln_phi[0] += 1.0e-5;
    const auto rejected_activity =
        ph::build_sw92_profile_c_thermodynamic_closure(activity_tamper, model);
    require(!rejected_activity.residual_available() &&
                rejected_activity.closure.primal_status ==
                    ph::ThermodynamicClosurePrimalStatus::indeterminate &&
                !rejected_activity.closure.linearization.has_value(),
            "closure accepted phase property that cannot be reproduced");

    auto status_tamper = source;
    status_tamper.solution.status = fl::PtPhaseSetStatus::indeterminate;
    const auto rejected_status =
        ph::build_sw92_profile_c_thermodynamic_closure(status_tamper, model);
    require(!rejected_status.residual_available() &&
                !rejected_status.closure.linearization.has_value(),
            "closure consumed a non-accepted phase-set status");

    bool caught = false;
    try {
        const auto wrong_model = binary_model(true);
        (void)ph::build_sw92_profile_c_thermodynamic_closure(source, wrong_model);
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    require(caught, "closure did not reject ordered model/source mismatch");
}

void component_permutation() {
    const auto normal_model = sample6::model(false);
    const auto reverse_model = sample6::model(true);
    const auto normal_source = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7, 350.0, sample6::feed(false), normal_model, 0.0);
    const auto reverse_source = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7, 350.0, sample6::feed(true), reverse_model, 0.0);
    const auto normal = ph::build_sw92_profile_c_thermodynamic_closure(
        normal_source, normal_model);
    const auto reverse = ph::build_sw92_profile_c_thermodynamic_closure(
        reverse_source, reverse_model);
    require_linearization_basics(normal, 3U);
    require_linearization_basics(reverse, 3U);

    for (std::size_t phase = 0; phase < 3U; ++phase) {
        const auto& a = normal.closure.primal->phases[phase];
        const auto& b = reverse.closure.primal->phases[phase];
        near(a.mole_phase_fraction, b.mole_phase_fraction);
        near(a.compressibility_factor, b.compressibility_factor);
        near(a.molar_density_mol_per_m3, b.molar_density_mol_per_m3);
        require(normal.phase_metadata[phase].physical_role ==
                    reverse.phase_metadata[phase].physical_role &&
                    normal.phase_metadata[phase].thermodynamic_family ==
                        reverse.phase_metadata[phase].thermodynamic_family,
                "component permutation changed phase metadata");
        require(a.composition.size() == b.composition.size(),
                "component permutation changed composition size");
        for (std::size_t i = 0; i < a.composition.size(); ++i) {
            near(a.composition[i], b.composition[a.composition.size() - 1U - i]);
        }
    }

    const auto& normal_l = *normal.closure.linearization;
    const auto& reverse_l = *reverse.closure.linearization;
    for (std::size_t phase = 0; phase < 3U; ++phase) {
        for (std::size_t column = 0; column < 2U; ++column) {
            near(normal_l.d_phase_fraction(phase, column),
                 reverse_l.d_phase_fraction(phase, column), 3e-6, 3e-10);
            near(normal_l.d_compressibility(phase, column),
                 reverse_l.d_compressibility(phase, column), 3e-6, 3e-10);
            near(normal_l.d_molar_density(phase, column),
                 reverse_l.d_molar_density(phase, column), 3e-6, 2e-8);
            for (std::size_t i = 0; i < 8U; ++i) {
                near(normal_l.d_composition(phase, i, column),
                     reverse_l.d_composition(phase, 7U - i, column),
                     4e-6, 4e-9);
            }
        }
    }
}

void headers() {
    require(sw92_thermodynamic_closure_header(),
            "SW92 thermodynamic closure public-header probe failed");
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "expected one test-case name\n";
        return 2;
    }
    const std::string_view name{argv[1]};
    try {
        if (name == "one_phase_reproduces_missing_z") one_phase_reproduces_missing_z();
        else if (name == "two_phase_preserves_authoritative_properties")
            two_phase_preserves_authoritative_properties();
        else if (name == "sample6_three_phase_primal") sample6_three_phase_primal();
        else if (name == "linearization_failure_policy") linearization_failure_policy();
        else if (name == "rejection_and_model_guard") rejection_and_model_guard();
        else if (name == "component_permutation") component_permutation();
        else if (name == "headers") headers();
        else throw std::invalid_argument("unknown test-case name");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
