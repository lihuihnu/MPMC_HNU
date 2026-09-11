#include <mpmc/ad/math.hpp>
#include <mpmc/ad/runtime_differentiate.hpp>
#include <mpmc/flash/sw92_profile_c_sensitivity.hpp>

#include "../sw92_phase_assigned_pt/physical_sample6.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

bool sw92_profile_c_sensitivity_header();

namespace {
namespace ad = mpmc::ad;
namespace fl = mpmc::flash;
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

void near(double actual, double expected, double relative = 5e-7,
          double absolute = 5e-10,
          std::source_location where = std::source_location::current()) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) >
            absolute + relative * std::max(std::abs(actual), std::abs(expected))) {
        std::cerr << "actual=" << actual << " expected=" << expected << '\n';
        require(false, "numeric mismatch", where);
    }
}

void near_increment(double analytic_increment, double resolved_half_difference,
                    double relative = 3e-3, double absolute = 2e-8,
                    std::source_location where = std::source_location::current()) {
    if (!std::isfinite(analytic_increment) ||
        !std::isfinite(resolved_half_difference) ||
        std::abs(analytic_increment - resolved_half_difference) >
            absolute + relative * std::max(
                std::abs(analytic_increment), std::abs(resolved_half_difference))) {
        std::cerr << "analytic_increment=" << analytic_increment
                  << " fresh_half_difference=" << resolved_half_difference << '\n';
        require(false, "fresh-resolve derivative mismatch", where);
    }
}

th::Sw92Phase<double> binary_model(bool reverse = false) {
    return th::Sw92Phase<double>::from_parameters(
        sw92_test::binary_parameters(sw92_test::co2, reverse));
}

double density(double pressure_pa, double temperature_k, double z) {
    return pressure_pa /
        (z * th::Sw92Pure<double>::gas_constant() * temperature_k);
}

void require_success(
    const fl::Sw92ProfileCPhaseSetSensitivityResult& sensitivity,
    std::size_t phases) {
    require(sensitivity.status == fl::PtSensitivityStatus::success,
            sensitivity.diagnostic);
    require(sensitivity.phase_count == phases &&
                sensitivity.component_count >= 2U &&
                sensitivity.input_count == sensitivity.component_count + 1U &&
                sensitivity.equilibrium_jacobian_rcond > 0.0 &&
                std::isfinite(sensitivity.equilibrium_jacobian_rcond) &&
                std::isfinite(sensitivity.linear_solve_backward_error),
            "successful sensitivity has invalid shape/conditioning diagnostics");
}

void phase_property_ad_crosscheck() {
    constexpr std::size_t width = 3;
    using Number = ad::Dual<double, width>;

    const auto model = binary_model();
    const auto published = fl::solve_sw92_profile_c_pt_phase_set(
        3.0e6, 340.0, Vec{0.7, 0.3}, model, 0.0);
    require(published.accepted_phase_set_published(),
            "AD phase-property fixture is not authoritative");
    const auto& phase = published.solution.accepted_phase_set()->phases[0];
    const auto family = published.phase_metadata[0].thermodynamic_family;
    require(phase.composition.size() == 2U,
            "AD phase-property fixture dimension changed");

    std::vector<double> inputs{
        published.solution.pressure_pa,
        published.solution.temperature_k,
        phase.composition[0]};
    ad::RuntimeJacobianWorkspace<double, width> workspace;
    th::Sw92PhaseWorkspace<Number> phase_workspace;
    const auto equation = [&](std::span<const Number> variables,
                              std::span<Number> outputs) {
        const std::span<const Number> reduced{variables.data() + 2U, 1U};
        const auto value = model.evaluate_reduced(
            variables[0], variables[1], reduced, 0.0, family,
            phase.activity.branch, phase_workspace);
        outputs[0] = value.z;
        outputs[1] = value.ln_phi[0];
        outputs[2] = value.ln_phi[1];
    };
    ad::RuntimeJacobianLimits limits;
    limits.max_inputs = 3U;
    limits.max_outputs = 3U;
    limits.max_jacobian_entries = 9U;
    const auto differentiated = ad::value_and_jacobian_runtime<width>(
        equation, inputs, 3U, workspace, limits);
    require(phase.compressibility_factor.has_value(),
            "AD phase-property fixture lost accepted Z");
    near(differentiated.values[0], *phase.compressibility_factor, 2e-10, 2e-12);
    near(differentiated.values[1], phase.activity.ln_phi[0], 2e-10, 2e-12);
    near(differentiated.values[2], phase.activity.ln_phi[1], 2e-10, 2e-12);

    const std::array<double, 3> step{300.0, 1e-3, 1e-6};
    const auto direct = [&](const std::array<double,3>& q) {
        const std::array<double,2> composition{q[2], 1.0-q[2]};
        th::Sw92PhaseWorkspace<double> local;
        const auto value = model.evaluate(
            q[0], q[1], composition, 0.0, family,
            phase.activity.branch, local);
        return std::array<double,3>{
            value.z, value.ln_phi[0], value.ln_phi[1]};
    };
    const std::array<double,3> base{inputs[0], inputs[1], inputs[2]};
    for (std::size_t column = 0; column < 3U; ++column) {
        auto plus = base;
        auto minus = base;
        plus[column] += step[column];
        minus[column] -= step[column];
        const auto fp = direct(plus);
        const auto fm = direct(minus);
        for (std::size_t output = 0; output < 3U; ++output) {
            const double finite_difference =
                (fp[output]-fm[output])/(2.0*step[column]);
            near(differentiated.jacobian[output*3U+column],
                 finite_difference, 8e-4, 2e-10);
        }
    }
}

void one_phase_fixed_set() {
    const auto model = binary_model();
    const auto source = fl::solve_sw92_profile_c_pt_phase_set(
        3.0e6, 340.0, Vec{0.995, 0.005}, model, 0.0);
    const auto sensitivity =
        fl::differentiate_sw92_profile_c_phase_set(source, model);
    require_success(sensitivity, 1U);
    for (std::size_t column = 0; column < sensitivity.input_count; ++column) {
        near(sensitivity.d_phase_fraction(0U,column), 0.0, 0.0, 2e-12);
    }
    near(sensitivity.d_composition(0U,0U,0U), 0.0, 0.0, 2e-12);
    near(sensitivity.d_composition(0U,1U,0U), 0.0, 0.0, 2e-12);
    near(sensitivity.d_composition(0U,0U,1U), 0.0, 0.0, 2e-12);
    near(sensitivity.d_composition(0U,1U,1U), 0.0, 0.0, 2e-12);
    near(sensitivity.d_composition(0U,0U,2U), 1.0, 0.0, 2e-12);
    near(sensitivity.d_composition(0U,1U,2U), -1.0, 0.0, 2e-12);
    require(std::isfinite(sensitivity.d_compressibility(0U,0U)) &&
                std::isfinite(sensitivity.d_compressibility(0U,1U)) &&
                std::isfinite(sensitivity.d_compressibility(0U,2U)) &&
                std::isfinite(sensitivity.d_molar_density(0U,0U)),
            "single-phase property derivative missing");
}

fl::Sw92ProfileCPtPhaseSetResult solve_wet_binary(
    double pressure_pa, double temperature_k, double z0,
    const th::Sw92Phase<double>& model) {
    return fl::solve_sw92_profile_c_pt_phase_set(
        pressure_pa, temperature_k, Vec{z0,1.0-z0}, model, 0.0);
}

void two_phase_fresh_resolve_crosscheck() {
    const auto model = binary_model();
    constexpr double p = 3.0e6;
    constexpr double t = 340.0;
    constexpr double z0 = 0.7;
    const auto source = solve_wet_binary(p,t,z0,model);
    const auto sensitivity =
        fl::differentiate_sw92_profile_c_phase_set(source, model);
    require_success(sensitivity, 2U);

    const std::array<double,3> step{500.0, 5e-3, 1e-5};
    for (std::size_t column = 0; column < 3U; ++column) {
        double pp=p, pm=p, tp=t, tm=t, zp=z0, zm=z0;
        if (column==0U) { pp+=step[column]; pm-=step[column]; }
        if (column==1U) { tp+=step[column]; tm-=step[column]; }
        if (column==2U) { zp+=step[column]; zm-=step[column]; }
        const auto plus = solve_wet_binary(pp,tp,zp,model);
        const auto minus = solve_wet_binary(pm,tm,zm,model);
        require(plus.accepted_phase_set_published() &&
                    minus.accepted_phase_set_published() &&
                    plus.solution.accepted_phase_count()==2U &&
                    minus.solution.accepted_phase_count()==2U,
                "fresh perturbation left the fixed two-phase topology");
        const auto& a = plus.solution.accepted_phase_set()->phases;
        const auto& b = minus.solution.accepted_phase_set()->phases;
        for (std::size_t phase=0;phase<2U;++phase) {
            near_increment(
                sensitivity.d_phase_fraction(phase,column)*step[column],
                (a[phase].mole_phase_fraction-b[phase].mole_phase_fraction)/2.0);
            for (std::size_t i=0;i<2U;++i) {
                near_increment(
                    sensitivity.d_composition(phase,i,column)*step[column],
                    (a[phase].composition[i]-b[phase].composition[i])/2.0);
            }
            require(a[phase].compressibility_factor.has_value() &&
                        b[phase].compressibility_factor.has_value(),
                    "fresh two-phase perturbation lost Z");
            near_increment(
                sensitivity.d_compressibility(phase,column)*step[column],
                (*a[phase].compressibility_factor-
                 *b[phase].compressibility_factor)/2.0);
            const double rho_plus = density(
                pp,tp,*a[phase].compressibility_factor);
            const double rho_minus = density(
                pm,tm,*b[phase].compressibility_factor);
            near_increment(
                sensitivity.d_molar_density(phase,column)*step[column],
                (rho_plus-rho_minus)/2.0, 5e-3, 5e-6);
        }
    }
}

void sample6_three_phase() {
    const auto model = sample6::model();
    const auto source = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7, 350.0, sample6::feed(), model, 0.0);
    const auto sensitivity =
        fl::differentiate_sw92_profile_c_phase_set(source, model);
    require_success(sensitivity, 3U);
    require(sensitivity.component_count == 8U &&
                sensitivity.input_count == 9U &&
                sensitivity.phase_metadata.size() == 3U &&
                sensitivity.equilibrium_residual_norm <=
                    fl::Sw92ProfileCSensitivityOptions{}.maximum_base_residual,
            "Sample-6 sensitivity shape/residual changed");
    require(sensitivity.phase_metadata[0].physical_role ==
                fl::Sw92PhaseAssignedPtPhysicalRole::aqueous &&
                sensitivity.phase_metadata[1].physical_role ==
                    fl::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified &&
                sensitivity.phase_metadata[2].physical_role ==
                    fl::Sw92PhaseAssignedPtPhysicalRole::nonaqueous_unclassified,
            "Sample-6 sensitivity invented hydrocarbon morphology");
    for (std::size_t phase=0;phase<3U;++phase) {
        for (std::size_t column=0;column<sensitivity.input_count;++column) {
            require(std::isfinite(sensitivity.d_phase_fraction(phase,column)) &&
                        std::isfinite(sensitivity.d_compressibility(phase,column)) &&
                        std::isfinite(sensitivity.d_molar_density(phase,column)),
                    "Sample-6 sensitivity contains nonfinite phase derivative");
        }
    }
}

void component_permutation_pt_columns() {
    const auto normal_model = sample6::model(false);
    const auto reverse_model = sample6::model(true);
    const auto normal_source = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7,350.0,sample6::feed(false),normal_model,0.0);
    const auto reverse_source = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7,350.0,sample6::feed(true),reverse_model,0.0);
    const auto normal = fl::differentiate_sw92_profile_c_phase_set(
        normal_source,normal_model);
    const auto reverse = fl::differentiate_sw92_profile_c_phase_set(
        reverse_source,reverse_model);
    require_success(normal,3U);
    require_success(reverse,3U);
    for (std::size_t phase=0;phase<3U;++phase) {
        for (std::size_t column=0;column<2U;++column) {
            near(normal.d_phase_fraction(phase,column),
                 reverse.d_phase_fraction(phase,column),2e-6,2e-10);
            near(normal.d_compressibility(phase,column),
                 reverse.d_compressibility(phase,column),2e-6,2e-10);
            near(normal.d_molar_density(phase,column),
                 reverse.d_molar_density(phase,column),2e-6,2e-8);
            for (std::size_t i=0;i<8U;++i) {
                near(normal.d_composition(phase,i,column),
                     reverse.d_composition(phase,7U-i,column),3e-6,3e-9);
            }
        }
    }
}

void h_slot_symmetry() {
    const auto model = sample6::model();
    auto source = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7,350.0,sample6::feed(),model,0.0);
    const auto original = fl::differentiate_sw92_profile_c_phase_set(source,model);
    require_success(original,3U);
    require(source.solution.candidate_phase_set.has_value(),
            "H-slot sensitivity fixture lost phase set");
    std::swap(source.solution.candidate_phase_set->phases[1],
              source.solution.candidate_phase_set->phases[2]);
    std::swap(source.phase_metadata[1],source.phase_metadata[2]);
    const auto swapped = fl::differentiate_sw92_profile_c_phase_set(source,model);
    require_success(swapped,3U);
    for (std::size_t column=0;column<2U;++column) {
        near(original.d_phase_fraction(0U,column),
             swapped.d_phase_fraction(0U,column),2e-6,2e-10);
        for (std::size_t i=0;i<8U;++i) {
            near(original.d_composition(0U,i,column),
                 swapped.d_composition(0U,i,column),3e-6,3e-9);
        }
        for (std::size_t phase=1U;phase<3U;++phase) {
            const std::size_t other=3U-phase;
            near(original.d_phase_fraction(phase,column),
                 swapped.d_phase_fraction(other,column),3e-6,3e-10);
            near(original.d_compressibility(phase,column),
                 swapped.d_compressibility(other,column),3e-6,3e-10);
            for (std::size_t i=0;i<8U;++i) {
                near(original.d_composition(phase,i,column),
                     swapped.d_composition(other,i,column),4e-6,4e-9);
            }
        }
    }
}

void boundary_and_tamper_guards() {
    const auto model = sample6::model();
    const auto source = fl::solve_sw92_profile_c_pt_phase_set(
        1.0e7,350.0,sample6::feed(),model,0.0);
    fl::Sw92ProfileCSensitivityOptions boundary;
    boundary.minimum_derivative_phase_fraction=0.04;
    const auto guarded = fl::differentiate_sw92_profile_c_phase_set(
        source,model,boundary);
    require(guarded.status==fl::PtSensitivityStatus::phase_boundary,
            "small physical phase bypassed derivative boundary guard");

    const auto binary = binary_model();
    auto tampered = fl::solve_sw92_profile_c_pt_phase_set(
        3.0e6,340.0,Vec{0.7,0.3},binary,0.0);
    require(tampered.solution.candidate_phase_set.has_value(),
            "tamper fixture lost accepted phase set");
    tampered.solution.candidate_phase_set->phases[0].activity.ln_phi[0] += 1e-5;
    const auto rejected = fl::differentiate_sw92_profile_c_phase_set(
        tampered,binary);
    require(rejected.status==fl::PtSensitivityStatus::solution_not_accepted,
            "tampered phase activity produced a sensitivity");

    bool mismatch=false;
    try {
        const auto reverse = binary_model(true);
        (void)fl::differentiate_sw92_profile_c_phase_set(tampered,reverse);
    } catch (const std::invalid_argument&) {
        mismatch=true;
    }
    require(mismatch,"ordered model mismatch was not rejected");
}

void headers() {
    require(sw92_profile_c_sensitivity_header(),
            "SW92 Profile-C sensitivity public header probe failed");
}

} // namespace

int main(int argc,char** argv) {
    if(argc!=2) {
        std::cerr<<"expected one test-case name\n";
        return 2;
    }
    const std::string_view name{argv[1]};
    try {
        if(name=="phase_property_ad_crosscheck") phase_property_ad_crosscheck();
        else if(name=="one_phase_fixed_set") one_phase_fixed_set();
        else if(name=="two_phase_fresh_resolve_crosscheck") two_phase_fresh_resolve_crosscheck();
        else if(name=="sample6_three_phase") sample6_three_phase();
        else if(name=="component_permutation_pt_columns") component_permutation_pt_columns();
        else if(name=="h_slot_symmetry") h_slot_symmetry();
        else if(name=="boundary_and_tamper_guards") boundary_and_tamper_guards();
        else if(name=="headers") headers();
        else throw std::invalid_argument("unknown test-case name");
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
    return 0;
}
