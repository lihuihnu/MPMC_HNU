#include <mpmc/simulation/simulation_cursor.hpp>

#include <iostream>
#include <limits>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace {
using mpmc::simulation::SimulationCursor;
using mpmc::simulation::SimulationTimeline;

void check(bool value,
           const std::source_location where = std::source_location::current()) {
    if (!value) {
        throw std::runtime_error(
            "simulation timeline test failed at line " +
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

void timeline_valid() {
    const SimulationTimeline timeline({0.0, 10.0, 25.0});
    check(timeline.size() == 3U);
    check(timeline.hard_boundary_seconds(0) == 0.0);
    check(timeline.hard_boundary_seconds(1) == 10.0);
    check(timeline.terminal_time_seconds() == 25.0);
}

void timeline_invalid() {
    invalid([] { (void)SimulationTimeline({}); });
    invalid([] { (void)SimulationTimeline({-1.0}); });
    invalid([] {
        (void)SimulationTimeline(
            {10.0, std::numeric_limits<double>::quiet_NaN()});
    });
    invalid([] {
        (void)SimulationTimeline(
            {10.0, std::numeric_limits<double>::infinity()});
    });
    invalid([] { (void)SimulationTimeline({10.0, 10.0}); });
    invalid([] { (void)SimulationTimeline({10.0, 9.0}); });
    const double epsilon = std::numeric_limits<double>::epsilon();
    invalid([&] {
        (void)SimulationTimeline({10.0, 10.0 + 16.0 * epsilon});
    });
}

void cursor_initialization() {
    SimulationCursor at_start(SimulationTimeline({0.0, 10.0, 20.0}), 0.0);
    check(at_start.consumed_boundary_count() == 1U);
    check(at_start.next_hard_boundary()->index == 1U);
    check(at_start.next_hard_boundary()->time_seconds == 10.0);

    SimulationCursor at_boundary(SimulationTimeline({10.0, 20.0}), 10.0);
    check(at_boundary.consumed_boundary_count() == 1U);
    check(at_boundary.next_hard_boundary()->time_seconds == 20.0);

    SimulationCursor between(SimulationTimeline({10.0, 20.0}), 12.0);
    check(between.consumed_boundary_count() == 1U);
    check(between.next_hard_boundary()->time_seconds == 20.0);

    SimulationCursor complete(SimulationTimeline({10.0, 20.0}), 20.0);
    check(complete.complete());
    check(!complete.next_hard_boundary().has_value());
    invalid([] {
        (void)SimulationCursor(SimulationTimeline({10.0, 20.0}), -1.0);
    });
    invalid([] {
        (void)SimulationCursor(
            SimulationTimeline({10.0, 20.0}),
            std::numeric_limits<double>::quiet_NaN());
    });
    invalid([] {
        (void)SimulationCursor(
            SimulationTimeline({10.0, 20.0}),
            std::numeric_limits<double>::infinity());
    });
    invalid([] {
        (void)SimulationCursor(SimulationTimeline({10.0, 20.0}), 21.0);
    });
}

void cap_composition() {
    SimulationCursor cursor(SimulationTimeline({10.0, 20.0}), 0.0);
    check(cursor.remaining_to_next_boundary_seconds(2.0) == 8.0);
    check(cursor.cap_initial_timestep_seconds(2.0) == 8.0);
    check(cursor.cap_initial_timestep_seconds(2.0, 3.0) == 3.0);
    check(cursor.cap_initial_timestep_seconds(2.0, 12.0) == 8.0);
    invalid([&] { (void)cursor.cap_initial_timestep_seconds(2.0, 0.0); });
    invalid([&] { (void)cursor.cap_initial_timestep_seconds(2.0, -1.0); });
    invalid([&] {
        (void)cursor.cap_initial_timestep_seconds(
            2.0, std::numeric_limits<double>::quiet_NaN());
    });
    invalid([&] {
        (void)cursor.cap_initial_timestep_seconds(
            2.0, std::numeric_limits<double>::infinity());
    });
    invalid([&] { (void)cursor.cap_initial_timestep_seconds(10.0); });
}

void accepted_progression() {
    SimulationCursor cursor(SimulationTimeline({10.0, 20.0}), 0.0);
    check(!cursor.accept_time(4.0).has_value());
    check(cursor.consumed_boundary_count() == 0U);

    const double epsilon = std::numeric_limits<double>::epsilon();
    const auto first = cursor.accept_time(10.0 + 16.0 * epsilon);
    check(first.has_value());
    check(first->index == 0U && first->time_seconds == 10.0);
    check(cursor.consumed_boundary_count() == 1U);

    check(!cursor.accept_time(15.0).has_value());
    const auto second = cursor.accept_time(20.0);
    check(second.has_value());
    check(second->index == 1U && second->time_seconds == 20.0);
    check(cursor.complete());
    check(!cursor.accept_time(20.0).has_value());
}

void boundary_skip_rejection() {
    SimulationCursor cursor(SimulationTimeline({10.0, 20.0}), 0.0);
    invalid([&] { (void)cursor.accept_time(11.0); });
    check(cursor.consumed_boundary_count() == 0U);
    invalid([&] { (void)cursor.remaining_to_next_boundary_seconds(11.0); });
    check(cursor.next_hard_boundary()->time_seconds == 10.0);
}

void terminal_state() {
    SimulationCursor cursor(SimulationTimeline({5.0}), 0.0);
    const auto reached = cursor.accept_time(5.0);
    check(reached.has_value());
    check(cursor.complete());
    invalid([&] { (void)cursor.remaining_to_next_boundary_seconds(5.0); });
    invalid([&] { (void)cursor.cap_initial_timestep_seconds(5.0); });
    invalid([&] { (void)cursor.accept_time(6.0); });
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument(
                "one simulation timeline test case is required");
        }
        const std::string_view test = argv[1];
        if (test == "timeline_valid") timeline_valid();
        else if (test == "timeline_invalid") timeline_invalid();
        else if (test == "cursor_initialization") cursor_initialization();
        else if (test == "cap_composition") cap_composition();
        else if (test == "accepted_progression") accepted_progression();
        else if (test == "boundary_skip_rejection") boundary_skip_rejection();
        else if (test == "terminal_state") terminal_state();
        else throw std::invalid_argument(
            "unknown simulation timeline test case");
        std::cout << "[PASS] simulation.core." << test << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
