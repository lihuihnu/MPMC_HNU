#include <mpmc/simulation/simulation_report.hpp>

#include <iostream>
#include <limits>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace {
using namespace mpmc::simulation;

void check(bool value,
           const std::source_location where = std::source_location::current()) {
    if (!value) {
        throw std::runtime_error(
            "simulation report test failed at line " +
            std::to_string(where.line()));
    }
}

template <class F>
void invalid(F&& function) {
    bool rejected = false;
    try {
        std::forward<F>(function)();
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected);
}

AcceptedSimulationStepSummary accepted_step(
    std::size_t index,
    double first,
    double second,
    std::size_t retries = 0U,
    std::size_t restarts = 0U) {
    return {
        index,
        first,
        second,
        second - first,
        2.0,
        retries,
        retries + 1U,
        restarts};
}

void status_contract() {
    check(!simulation_status_is_terminal(SimulationRunStatus::running));
    check(simulation_status_is_terminal(SimulationRunStatus::completed));
    check(simulation_status_name(SimulationRunStatus::running) == "running");
    check(simulation_status_name(SimulationRunStatus::execution_error) ==
          "execution_error");
}

void accepted_summary() {
    SimulationReport report({10.0, 4U, 2U});
    report.record_physical_timestep_call();
    report.record_accepted_step(accepted_step(4U, 10.0, 12.0, 1U, 2U));
    check(report.status() == SimulationRunStatus::running);
    check(report.diagnostics().physical_timestep_calls == 1U);
    check(report.diagnostics().accepted_steps == 1U);
    check(report.diagnostics().accepted_retries == 1U);
    check(report.diagnostics().accepted_phase_transition_restarts == 2U);
    check(report.last_accepted_step()->accepted_step_index == 4U);
    check(report.last_accepted_step()->time_np1_seconds == 12.0);

    report.record_physical_timestep_call();
    report.record_accepted_step(accepted_step(5U, 12.0, 14.0));
    check(report.diagnostics().accepted_steps == 2U);
    check(report.last_accepted_step()->accepted_step_index == 5U);
}

void accepted_summary_invalid() {
    invalid([] { (void)SimulationReport({-1.0, 0U, 0U}); });
    SimulationReport report({10.0, 4U, 2U});
    invalid([&] { report.record_accepted_step(accepted_step(4U, 10.0, 12.0)); });
    report.record_physical_timestep_call();

    auto bad = accepted_step(4U, 10.0, 12.0);
    bad.accepted_step_index = 5U;
    invalid([&] { report.record_accepted_step(bad); });

    bad = accepted_step(4U, 11.0, 13.0);
    invalid([&] { report.record_accepted_step(bad); });

    bad = accepted_step(4U, 10.0, 12.0);
    bad.accepted_timestep_seconds = 1.0;
    invalid([&] { report.record_accepted_step(bad); });

    bad = accepted_step(4U, 10.0, 12.0);
    bad.next_timestep_seconds = 0.0;
    invalid([&] { report.record_accepted_step(bad); });

    bad = accepted_step(4U, 10.0, 12.0);
    bad.attempt_count = 0U;
    invalid([&] { report.record_accepted_step(bad); });

    bad = accepted_step(4U, 10.0, 12.0, 1U);
    bad.attempt_count = 3U;
    invalid([&] { report.record_accepted_step(bad); });

    bad = accepted_step(4U, 10.0, 12.0);
    bad.time_np1_seconds = std::numeric_limits<double>::infinity();
    invalid([&] { report.record_accepted_step(bad); });
}

void boundary_identity() {
    SimulationReport report({10.0, 4U, 2U});
    report.record_physical_timestep_call();
    report.record_accepted_step(accepted_step(4U, 10.0, 20.0));
    report.record_reached_boundary({{2U, 20.0}, 4U});
    check(report.diagnostics().reached_boundaries == 1U);
    check(report.last_reached_boundary()->boundary.index == 2U);
    check(report.last_reached_boundary()->boundary.time_seconds == 20.0);

    report.record_physical_timestep_call();
    report.record_accepted_step(accepted_step(5U, 20.0, 30.0));
    report.record_reached_boundary({{3U, 30.0}, 5U});
    check(report.diagnostics().reached_boundaries == 2U);
}

void boundary_identity_invalid() {
    SimulationReport report({10.0, 4U, 2U});
    invalid([&] { report.record_reached_boundary({{2U, 10.0}, 4U}); });
    report.record_physical_timestep_call();
    report.record_accepted_step(accepted_step(4U, 10.0, 20.0));
    invalid([&] { report.record_reached_boundary({{3U, 20.0}, 4U}); });
    invalid([&] { report.record_reached_boundary({{2U, 19.0}, 4U}); });
    invalid([&] { report.record_reached_boundary({{2U, 20.0}, 5U}); });
    report.record_reached_boundary({{2U, 20.0}, 4U});
    invalid([&] { report.record_reached_boundary({{3U, 20.0}, 4U}); });
}

void terminal_contract() {
    SimulationReport completed({0.0, 0U, 0U});
    completed.terminate(SimulationRunStatus::completed);
    check(completed.terminal());
    check(completed.status() == SimulationRunStatus::completed);
    check(!completed.lower_layer_error_code().has_value());
    invalid([&] { completed.record_physical_timestep_call(); });
    invalid([&] { completed.terminate(SimulationRunStatus::completed); });

    SimulationReport rejected({0.0, 0U, 0U});
    rejected.terminate(SimulationRunStatus::physical_timestep_rejected);
    check(rejected.status() == SimulationRunStatus::physical_timestep_rejected);

    SimulationReport execution({0.0, 0U, 0U});
    execution.terminate(SimulationRunStatus::execution_error, -7);
    check(execution.lower_layer_error_code() ==
          std::optional<std::int64_t>{-7});

    SimulationReport invalid_error({0.0, 0U, 0U});
    invalid([&] { invalid_error.terminate(SimulationRunStatus::running); });
    invalid([&] { invalid_error.terminate(SimulationRunStatus::execution_error); });
    invalid([&] {
        invalid_error.terminate(SimulationRunStatus::execution_error, 0);
    });
    invalid([&] { invalid_error.terminate(SimulationRunStatus::completed, 1); });
}

void bounded_diagnostics() {
    SimulationReport report({0.0, 0U, 0U});
    for (std::size_t i = 0; i < 1000U; ++i) {
        report.record_physical_timestep_call();
        report.record_accepted_step(
            accepted_step(
                i,
                static_cast<double>(i),
                static_cast<double>(i + 1U)));
    }
    check(report.diagnostics().physical_timestep_calls == 1000U);
    check(report.diagnostics().accepted_steps == 1000U);
    check(report.last_accepted_step()->accepted_step_index == 999U);
    check(!report.last_reached_boundary().has_value());
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument(
                "one simulation report test case is required");
        }
        const std::string_view test_name = argv[1];
        if (test_name == "status_contract") status_contract();
        else if (test_name == "accepted_summary") accepted_summary();
        else if (test_name == "accepted_summary_invalid") accepted_summary_invalid();
        else if (test_name == "boundary_identity") boundary_identity();
        else if (test_name == "boundary_identity_invalid") {
            boundary_identity_invalid();
        } else if (test_name == "terminal_contract") terminal_contract();
        else if (test_name == "bounded_diagnostics") bounded_diagnostics();
        else throw std::invalid_argument(
            "unknown simulation report test case");
        std::cout << "[PASS] simulation.core.report." << test_name << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
