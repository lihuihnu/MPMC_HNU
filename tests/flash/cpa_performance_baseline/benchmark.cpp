#include <mpmc/flash/cpa_split.hpp>

#include "test_support.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {
namespace fl = mpmc::flash;
namespace th = mpmc::thermodynamics;
using Clock = std::chrono::steady_clock;
using Vec = std::vector<double>;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

fl::PtSplitOptions benchmark_split_options() {
    fl::PtSplitOptions options;
    options.initial_stability.automatic_starts = false;
    options.final_stability.automatic_starts = false;
    return options;
}

struct StructuralCounts {
    std::size_t states{};
    std::size_t initial_stability_evaluations{};
    std::size_t final_stability_evaluations{};
    std::size_t stability_trials{};
    std::size_t stability_trial_iterations{};
    std::size_t split_evaluations{};
    std::size_t split_attempts{};
    std::size_t split_attempt_evaluations{};
    std::size_t split_attempt_iterations{};
    std::size_t split_backtracks{};
    double checksum{};
};

void add_stability_counts(const fl::StabilityResult& result, StructuralCounts& counts) {
    counts.stability_trials += result.trials.size();
    for (const auto& trial : result.trials) {
        if (trial.iterations > 0) {
            counts.stability_trial_iterations += static_cast<std::size_t>(trial.iterations);
        }
    }
}

void add_split_counts(const fl::PtSplitResult& result, StructuralCounts& counts) {
    counts.initial_stability_evaluations += result.initial_stability.evaluations;
    add_stability_counts(result.initial_stability, counts);
    if (result.final_stability.has_value()) {
        counts.final_stability_evaluations += result.final_stability->evaluations;
        add_stability_counts(*result.final_stability, counts);
    }
    counts.split_evaluations += result.split_evaluations;
    counts.split_attempts += result.attempts.size();
    for (const auto& attempt : result.attempts) {
        counts.split_attempt_evaluations += attempt.evaluations;
        counts.split_backtracks += attempt.backtracks;
        if (attempt.iterations > 0) {
            counts.split_attempt_iterations += static_cast<std::size_t>(attempt.iterations);
        }
    }
}

StructuralCounts run_full_workload(bool emit_state_lines) {
    const auto parameters = cpa_physical_test::parameters(false);
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    fl::CpaVleEvaluator evaluator(model);
    const auto options = benchmark_split_options();

    StructuralCounts counts;
    for (const auto& experimental : cpa_physical_test::points()) {
        const auto feed = cpa_physical_test::feed(experimental, false);
        const auto starts = cpa_physical_test::starts(experimental, false);
        const std::vector<Vec> final_starts{feed};
        const auto result = fl::solve_cpa_pt_vle(
            experimental.pressure_pa, cpa_physical_test::temperature_k,
            feed, evaluator, options, starts, final_starts);
        require(result.solution.status ==
                    fl::PtSplitStatus::two_phase_no_instability_found &&
                    result.solution.candidate() != nullptr &&
                    result.solution.final_stability.has_value() &&
                    result.solution.final_stability->status ==
                        fl::StabilityStatus::no_instability_found,
                "CPA performance workload did not close as the validated two-phase state");
        const auto& point = *result.solution.candidate();
        counts.checksum += point.fractions.vapor_fraction;
        counts.checksum += point.fractions.liquid.at(0);
        counts.checksum += point.fractions.vapor.at(0);
        add_split_counts(result.solution, counts);
        ++counts.states;

        if (emit_state_lines) {
            std::cout << "CPA_PERF_STATE"
                      << " pressure_pa=" << experimental.pressure_pa
                      << " initial_stability_evals="
                      << result.solution.initial_stability.evaluations
                      << " final_stability_evals="
                      << result.solution.final_stability->evaluations
                      << " split_evals=" << result.solution.split_evaluations
                      << " attempts=" << result.solution.attempts.size()
                      << " betaV=" << point.fractions.vapor_fraction
                      << '\n';
        }
    }
    return counts;
}

struct PhaseCase {
    double pressure_pa{};
    Vec composition;
    fl::PtPhaseRole role{fl::PtPhaseRole::liquid_candidate};
    double density_mol_per_m3{};
};

const th::CpaPtRoot& select_root(const th::CpaPtRootSet& roots,
                                 fl::PtPhaseRole role) {
    require(roots.status == th::CpaPtRootStatus::success,
            "CPA performance micro workload root search failed");
    const th::CpaPtRoot* selected = nullptr;
    for (const auto& root : roots.roots) {
        if (root.pressure_slope_sign <= 0) { continue; }
        if (selected == nullptr ||
            (role == fl::PtPhaseRole::liquid_candidate &&
             root.molar_density_mol_per_m3 > selected->molar_density_mol_per_m3) ||
            (role == fl::PtPhaseRole::vapor_candidate &&
             root.molar_density_mol_per_m3 < selected->molar_density_mol_per_m3)) {
            selected = &root;
        }
    }
    require(selected != nullptr,
            "CPA performance micro workload found no admissible density root");
    return *selected;
}

std::vector<PhaseCase> prepare_phase_cases(const th::CpaPtPhase& model,
                                           const th::CpaPtOptions& pt_options) {
    std::vector<PhaseCase> cases;
    cases.reserve(2U * cpa_physical_test::points().size());
    for (const auto& point : cpa_physical_test::points()) {
        for (const auto role : {fl::PtPhaseRole::liquid_candidate,
                                fl::PtPhaseRole::vapor_candidate}) {
            const double methanol = role == fl::PtPhaseRole::liquid_candidate
                ? point.liquid_methanol : point.vapor_methanol;
            auto composition = cpa_physical_test::composition(methanol, false);
            const auto roots = model.roots(
                point.pressure_pa, cpa_physical_test::temperature_k,
                composition, pt_options);
            const auto& root = select_root(roots, role);
            cases.push_back({point.pressure_pa, std::move(composition), role,
                             root.molar_density_mol_per_m3});
        }
    }
    return cases;
}

template <typename Function>
long long measure_ns(Function&& function) {
    const auto start = Clock::now();
    function();
    const auto stop = Clock::now();
    return std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
}

struct MicroTiming {
    long long association_ns{};
    long long phase_ns{};
    long long root_ns{};
    long long stability_ns{};
    long long split_ns{};
    std::size_t association_units{};
    std::size_t phase_units{};
    std::size_t root_units{};
    std::size_t stability_units{};
    std::size_t split_units{};
    std::size_t root_density_evaluations{};
    std::size_t association_iterations{};
    double checksum{};
};

MicroTiming run_micro_workload() {
    constexpr std::size_t association_passes = 128U;
    constexpr std::size_t phase_passes = 64U;
    constexpr std::size_t root_passes = 2U;
    constexpr std::size_t stability_passes = 2U;
    constexpr std::size_t split_passes = 2U;

    const auto parameters = cpa_physical_test::parameters(false);
    const auto model = th::CpaPtPhase::from_parameters(parameters);
    const auto pt_options = fl::cpa_pt_vle_default_phase_options();
    const auto cases = prepare_phase_cases(model, pt_options);
    fl::CpaStabilityEvaluator stability(model, pt_options);
    fl::CpaVleEvaluator split(model, pt_options);

    MicroTiming timing;
    timing.association_units = association_passes * cases.size();
    timing.phase_units = phase_passes * cases.size();
    timing.root_units = root_passes * cases.size();
    timing.stability_units = stability_passes * cases.size();
    timing.split_units = split_passes * cases.size();

    timing.association_ns = measure_ns([&] {
        for (std::size_t pass = 0; pass < association_passes; ++pass) {
            for (const auto& item : cases) {
                const auto result = th::solve_cpa_association(
                    cpa_physical_test::temperature_k, item.density_mol_per_m3,
                    item.composition, parameters, pt_options.phase.association);
                require(result.converged(),
                        "CPA performance association micro workload failed");
                timing.association_iterations +=
                    static_cast<std::size_t>(std::max(result.iterations, 0));
                timing.checksum += result.radial_distribution;
            }
        }
    });

    timing.phase_ns = measure_ns([&] {
        for (std::size_t pass = 0; pass < phase_passes; ++pass) {
            for (const auto& item : cases) {
                const auto state = th::evaluate_cpa_phase_at_density(
                    cpa_physical_test::temperature_k, item.density_mol_per_m3,
                    item.composition, parameters, pt_options.phase);
                timing.checksum += state.pressure_pa * 1.0e-9;
            }
        }
    });

    timing.root_ns = measure_ns([&] {
        for (std::size_t pass = 0; pass < root_passes; ++pass) {
            for (const auto& item : cases) {
                const auto roots = model.roots(
                    item.pressure_pa, cpa_physical_test::temperature_k,
                    item.composition, pt_options);
                require(roots.status == th::CpaPtRootStatus::success,
                        "CPA performance root micro workload failed");
                timing.root_density_evaluations += roots.evaluations;
                timing.checksum += static_cast<double>(roots.roots.size());
            }
        }
    });

    timing.stability_ns = measure_ns([&] {
        for (std::size_t pass = 0; pass < stability_passes; ++pass) {
            for (const auto& item : cases) {
                const auto phase = stability(
                    item.pressure_pa, cpa_physical_test::temperature_k,
                    item.composition);
                timing.checksum += phase.ln_phi.at(0);
            }
        }
    });

    timing.split_ns = measure_ns([&] {
        for (std::size_t pass = 0; pass < split_passes; ++pass) {
            for (const auto& item : cases) {
                const auto phase = split(
                    item.pressure_pa, cpa_physical_test::temperature_k,
                    item.composition, item.role);
                timing.checksum += phase.z;
            }
        }
    });
    return timing;
}

void print_structure(const StructuralCounts& counts) {
    std::cout << "CPA_PERF_STRUCTURE"
              << " states=" << counts.states
              << " initial_stability_evals=" << counts.initial_stability_evaluations
              << " final_stability_evals=" << counts.final_stability_evaluations
              << " stability_trials=" << counts.stability_trials
              << " stability_trial_iterations=" << counts.stability_trial_iterations
              << " split_evals=" << counts.split_evaluations
              << " split_attempts=" << counts.split_attempts
              << " split_attempt_evals=" << counts.split_attempt_evaluations
              << " split_attempt_iterations=" << counts.split_attempt_iterations
              << " split_backtracks=" << counts.split_backtracks
              << " checksum=" << counts.checksum
              << '\n';
}

void run_full_timing() {
    const auto start = Clock::now();
    const auto counts = run_full_workload(false);
    const auto stop = Clock::now();
    const auto elapsed =
        std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    std::cout << "CPA_PERF_FULL"
              << " states=" << counts.states
              << " elapsed_ns=" << elapsed
              << " checksum=" << counts.checksum
              << '\n';
}

void run_micro_timing() {
    const auto timing = run_micro_workload();
    std::cout << "CPA_PERF_MICRO"
              << " association_ns=" << timing.association_ns
              << " association_units=" << timing.association_units
              << " association_iterations=" << timing.association_iterations
              << " phase_ns=" << timing.phase_ns
              << " phase_units=" << timing.phase_units
              << " root_ns=" << timing.root_ns
              << " root_units=" << timing.root_units
              << " root_density_evaluations=" << timing.root_density_evaluations
              << " stability_ns=" << timing.stability_ns
              << " stability_units=" << timing.stability_units
              << " split_ns=" << timing.split_ns
              << " split_units=" << timing.split_units
              << " checksum=" << timing.checksum
              << '\n';
}

} // namespace

int main(int argc, char** argv) {
    try {
        std::cout << std::setprecision(17);
        const std::string_view mode = argc >= 2 ? argv[1] : "structure";
        if (mode == "structure") {
            print_structure(run_full_workload(true));
        } else if (mode == "full") {
            run_full_timing();
        } else if (mode == "micro") {
            run_micro_timing();
        } else {
            throw std::invalid_argument("usage: mpmc_cpa_performance_baseline [structure|full|micro]");
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
