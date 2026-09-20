#include <mpmc/flow/phase_potential_upwind.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool phase_potential_upwind_header();

namespace {

namespace fl = mpmc::flow;

void require(
    bool condition,
    std::string_view message,
    std::source_location where =
        std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(
            std::string{where.file_name()} +
            ":" +
            std::to_string(where.line()) +
            ": " +
            std::string{message});
    }
}

void near(
    double actual,
    double expected,
    double relative = 2.0e-13,
    double absolute = 2.0e-13,
    std::source_location where =
        std::source_location::current()) {
    if (!std::isfinite(actual) ||
        !std::isfinite(expected) ||
        std::abs(actual - expected) >
            absolute +
                relative *
                    std::max(
                        std::abs(actual),
                        std::abs(expected))) {
        std::cerr
            << "actual=" << actual
            << " expected=" << expected
            << '\n';
        require(false, "numeric mismatch", where);
    }
}

template <class Function>
void expect_invalid(
    Function&& function,
    std::string_view fragment) {
    try {
        function();
    } catch (const std::invalid_argument& error) {
        require(
            std::string_view{error.what()}.find(fragment) !=
                std::string_view::npos,
            "invalid_argument diagnostic changed");
        return;
    }
    throw std::runtime_error(
        "expected std::invalid_argument");
}

fl::NaturalVariableStateIdentity3P state_identity(
    bool owner) {
    const auto pivot =
        fl::NaturalVariableCompositionPivot3P::
            from_dependent_components(
                3U,
                owner
                    ? std::array<std::size_t, 3>{1U, 0U, 2U}
                    : std::array<std::size_t, 3>{0U, 2U, 1U});

    return {
        fl::NaturalVariableLayout3P{pivot},
        {"A", "B", "C"},
        owner ? 1.0e6 : 0.99e6,
        owner ? 350.0 : 351.0,
        owner
            ? std::array<double, 3>{0.20, 0.30, 0.50}
            : std::array<double, 3>{0.25, 0.35, 0.40},
        owner
            ? std::array<std::vector<double>, 3>{
                  std::vector<double>{0.10, 0.70, 0.20},
                  std::vector<double>{0.60, 0.20, 0.20},
                  std::vector<double>{0.20, 0.30, 0.50}}
            : std::array<std::vector<double>, 3>{
                  std::vector<double>{0.55, 0.25, 0.20},
                  std::vector<double>{0.25, 0.25, 0.50},
                  std::vector<double>{0.20, 0.60, 0.20}}};
}

std::array<std::vector<double>, 3> gradients(
    const fl::NaturalVariableLayout3P& layout,
    double base) {
    std::array<std::vector<double>, 3> result;
    const std::size_t q = layout.unknown_count();
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        result[phase].resize(q);
        for (std::size_t column = 0U;
             column < q;
             ++column) {
            result[phase][column] =
                base *
                static_cast<double>(
                    (phase + 1U) * 10U +
                    column + 1U);
        }
    }
    return result;
}

fl::LocalPhaseMobilityLinearization3P local_payload(
    bool owner) {
    auto identity =
        state_identity(owner);
    const auto p_gradient =
        gradients(
            identity.layout,
            owner ? 100.0 : 120.0);
    const auto rho_gradient =
        gradients(
            identity.layout,
            owner ? 0.01 : 0.02);
    const auto mu_gradient =
        gradients(
            identity.layout,
            owner ? 1.0e-6 : 1.2e-6);
    const auto kr_gradient =
        gradients(
            identity.layout,
            owner ? 0.001 : 0.0012);
    const auto mobility_gradient =
        gradients(
            identity.layout,
            owner ? 0.02 : 0.03);

    return {
        std::move(identity),
        {"synthetic-rho", "face-fixture", "v1"},
        {"synthetic-mu", "face-fixture", "v1"},
        owner
            ? std::array<double, 3>{
                  1.0e6, 1.01e6, 1.02e6}
            : std::array<double, 3>{
                  0.99e6, 0.99e6, 1.005e6},
        owner
            ? std::array<double, 3>{
                  900.0, 400.0, 700.0}
            : std::array<double, 3>{
                  1100.0, 600.0, 800.0},
        {1.0, 1.0, 1.0},
        owner
            ? std::array<double, 3>{1.0, 2.0, 3.0}
            : std::array<double, 3>{4.0, 5.0, 6.0},
        owner
            ? std::array<double, 3>{1.0, 2.0, 3.0}
            : std::array<double, 3>{4.0, 5.0, 6.0},
        p_gradient,
        rho_gradient,
        mu_gradient,
        kr_gradient,
        mobility_gradient};
}

void potential_and_upwind() {
    const auto owner =
        local_payload(true);
    const auto neighbour =
        local_payload(false);

    const auto result =
        fl::build_two_cell_phase_potential_upwind_linearization(
            owner,
            neighbour,
            {0.0, 0.0, -10.0},
            {0.0, 0.0, 2.0});

    require(
        result.density_policy ==
            fl::FacePhaseDensityPolicy3P::
                arithmetic_mean_owner_neighbour,
        "face density policy changed");
    near(
        result.phase[0].gravity_projection_m2_per_s2,
        -20.0);
    near(
        result.phase[1].gravity_projection_m2_per_s2,
        -20.0);
    near(
        result.phase[2].gravity_projection_m2_per_s2,
        -20.0);

    const std::array<double, 3>
        expected_face_density{
            1000.0, 500.0, 750.0};
    const std::array<double, 3>
        expected_gravity_pressure{
            -20000.0, -10000.0, -15000.0};
    const std::array<double, 3>
        expected_potential{
            10000.0, -10000.0, 0.0};

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto& entry =
            result.phase[phase];
        near(
            entry.face_mass_density_kg_per_m3,
            expected_face_density[phase]);
        near(
            entry.gravity_pressure_difference_pa,
            expected_gravity_pressure[phase]);
        near(
            entry.phase_potential_difference_pa,
            expected_potential[phase]);

        for (std::size_t column = 0U;
             column <
             owner.state_identity.layout.unknown_count();
             ++column) {
            near(
                entry.owner_face_density_gradient[column],
                0.5 *
                    owner.mass_density_gradient
                        [phase][column]);
            near(
                entry.owner_phase_potential_gradient[column],
                -owner.phase_pressure_gradient
                    [phase][column] -
                    entry.owner_face_density_gradient[column] *
                        (-20.0));
        }
        for (std::size_t column = 0U;
             column <
             neighbour.state_identity.layout.unknown_count();
             ++column) {
            near(
                entry.neighbour_face_density_gradient[column],
                0.5 *
                    neighbour.mass_density_gradient
                        [phase][column]);
            near(
                entry.neighbour_phase_potential_gradient[column],
                neighbour.phase_pressure_gradient
                    [phase][column] -
                    entry.neighbour_face_density_gradient[column] *
                        (-20.0));
        }
    }

    require(
        result.phase[0].upwind_selection ==
            fl::UpwindCellSelection3P::
                neighbour_positive_phase_potential &&
            result.phase[1].upwind_selection ==
            fl::UpwindCellSelection3P::
                owner_negative_phase_potential &&
            result.phase[2].upwind_selection ==
            fl::UpwindCellSelection3P::
                owner_exact_zero_tie,
        "upwind selection sign convention changed");

    near(
        result.phase[0].upwind_mobility_per_pa_s,
        neighbour.mobility_per_pa_s[0]);
    near(
        result.phase[1].upwind_mobility_per_pa_s,
        owner.mobility_per_pa_s[1]);
    near(
        result.phase[2].upwind_mobility_per_pa_s,
        owner.mobility_per_pa_s[2]);

    require(
        result.phase[0]
                .owner_upwind_mobility_gradient ==
            std::vector<double>(
                owner.state_identity.layout.unknown_count(),
                0.0) &&
            result.phase[0]
                    .neighbour_upwind_mobility_gradient ==
                neighbour.mobility_gradient[0] &&
            result.phase[1]
                    .owner_upwind_mobility_gradient ==
                owner.mobility_gradient[1] &&
            result.phase[1]
                    .neighbour_upwind_mobility_gradient ==
                std::vector<double>(
                    neighbour.state_identity.layout.unknown_count(),
                    0.0) &&
            result.phase[2]
                    .owner_upwind_mobility_gradient ==
                owner.mobility_gradient[2],
        "upwind mobility Jacobian does not follow selected branch");
}


fl::LocalPhaseMobilityLinearization3P perturb_local(
    fl::LocalPhaseMobilityLinearization3P local,
    std::size_t column,
    double delta) {
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        local.phase_pressure_pa[phase] +=
            local.phase_pressure_gradient[phase][column] *
            delta;
        local.mass_density_kg_per_m3[phase] +=
            local.mass_density_gradient[phase][column] *
            delta;
        local.mobility_per_pa_s[phase] +=
            local.mobility_gradient[phase][column] *
            delta;
    }
    return local;
}

void jacobian_fresh_perturbation() {
    const auto owner =
        local_payload(true);
    const auto neighbour =
        local_payload(false);
    const auto base =
        fl::build_two_cell_phase_potential_upwind_linearization(
            owner,
            neighbour,
            {0.0, 0.0, -10.0},
            {0.0, 0.0, 2.0});

    constexpr double step = 1.0e-6;

    for (std::size_t column = 0U;
         column <
         owner.state_identity.layout.unknown_count();
         ++column) {
        const auto plus =
            fl::build_two_cell_phase_potential_upwind_linearization(
                perturb_local(owner, column, step),
                neighbour,
                {0.0, 0.0, -10.0},
                {0.0, 0.0, 2.0});
        const auto minus =
            fl::build_two_cell_phase_potential_upwind_linearization(
                perturb_local(owner, column, -step),
                neighbour,
                {0.0, 0.0, -10.0},
                {0.0, 0.0, 2.0});

        for (std::size_t phase = 0U;
             phase < 2U;
             ++phase) {
            const double fd_density =
                (plus.phase[phase]
                     .face_mass_density_kg_per_m3 -
                 minus.phase[phase]
                     .face_mass_density_kg_per_m3) /
                (2.0 * step);
            const double fd_potential =
                (plus.phase[phase]
                     .phase_potential_difference_pa -
                 minus.phase[phase]
                     .phase_potential_difference_pa) /
                (2.0 * step);
            const double fd_mobility =
                (plus.phase[phase]
                     .upwind_mobility_per_pa_s -
                 minus.phase[phase]
                     .upwind_mobility_per_pa_s) /
                (2.0 * step);

            near(
                base.phase[phase]
                    .owner_face_density_gradient[column],
                fd_density,
                2.0e-7,
                2.0e-7);
            near(
                base.phase[phase]
                    .owner_phase_potential_gradient[column],
                fd_potential,
                2.0e-7,
                2.0e-5);
            near(
                base.phase[phase]
                    .owner_upwind_mobility_gradient[column],
                fd_mobility,
                2.0e-7,
                2.0e-7);
        }
    }

    for (std::size_t column = 0U;
         column <
         neighbour.state_identity.layout.unknown_count();
         ++column) {
        const auto plus =
            fl::build_two_cell_phase_potential_upwind_linearization(
                owner,
                perturb_local(neighbour, column, step),
                {0.0, 0.0, -10.0},
                {0.0, 0.0, 2.0});
        const auto minus =
            fl::build_two_cell_phase_potential_upwind_linearization(
                owner,
                perturb_local(neighbour, column, -step),
                {0.0, 0.0, -10.0},
                {0.0, 0.0, 2.0});

        for (std::size_t phase = 0U;
             phase < 2U;
             ++phase) {
            const double fd_density =
                (plus.phase[phase]
                     .face_mass_density_kg_per_m3 -
                 minus.phase[phase]
                     .face_mass_density_kg_per_m3) /
                (2.0 * step);
            const double fd_potential =
                (plus.phase[phase]
                     .phase_potential_difference_pa -
                 minus.phase[phase]
                     .phase_potential_difference_pa) /
                (2.0 * step);
            const double fd_mobility =
                (plus.phase[phase]
                     .upwind_mobility_per_pa_s -
                 minus.phase[phase]
                     .upwind_mobility_per_pa_s) /
                (2.0 * step);

            near(
                base.phase[phase]
                    .neighbour_face_density_gradient[column],
                fd_density,
                2.0e-7,
                2.0e-7);
            near(
                base.phase[phase]
                    .neighbour_phase_potential_gradient[column],
                fd_potential,
                2.0e-7,
                2.0e-5);
            near(
                base.phase[phase]
                    .neighbour_upwind_mobility_gradient[column],
                fd_mobility,
                2.0e-7,
                2.0e-7);
        }
    }
}

void orientation_reversal() {
    const auto owner =
        local_payload(true);
    const auto neighbour =
        local_payload(false);

    const auto forward =
        fl::build_two_cell_phase_potential_upwind_linearization(
            owner,
            neighbour,
            {0.0, 0.0, -10.0},
            {0.0, 0.0, 2.0});
    const auto reverse =
        fl::build_two_cell_phase_potential_upwind_linearization(
            neighbour,
            owner,
            {0.0, 0.0, -10.0},
            {0.0, 0.0, -2.0});

    for (std::size_t phase = 0U;
         phase < 2U;
         ++phase) {
        near(
            reverse.phase[phase]
                .phase_potential_difference_pa,
            -forward.phase[phase]
                 .phase_potential_difference_pa);
        near(
            reverse.phase[phase]
                .face_mass_density_kg_per_m3,
            forward.phase[phase]
                .face_mass_density_kg_per_m3);
        near(
            reverse.phase[phase]
                .upwind_mobility_per_pa_s,
            forward.phase[phase]
                .upwind_mobility_per_pa_s);
    }

    require(
        reverse.phase[0].upwind_selection ==
            fl::UpwindCellSelection3P::
                owner_negative_phase_potential &&
            reverse.phase[1].upwind_selection ==
            fl::UpwindCellSelection3P::
                neighbour_positive_phase_potential,
        "orientation reversal did not swap owner/neighbour branch labels");

    require(
        forward.phase[2].upwind_selection ==
                fl::UpwindCellSelection3P::
                    owner_exact_zero_tie &&
            reverse.phase[2].upwind_selection ==
                fl::UpwindCellSelection3P::
                    owner_exact_zero_tie,
        "exact-zero tie convention changed under orientation reversal");
}

void zero_gravity() {
    const auto owner =
        local_payload(true);
    const auto neighbour =
        local_payload(false);

    const auto result =
        fl::build_two_cell_phase_potential_upwind_linearization(
            owner,
            neighbour,
            {0.0, 0.0, 0.0},
            {1.0, 2.0, 3.0});

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        const auto& entry =
            result.phase[phase];
        near(
            entry.gravity_projection_m2_per_s2,
            0.0);
        near(
            entry.gravity_pressure_difference_pa,
            0.0);
        near(
            entry.phase_potential_difference_pa,
            neighbour.phase_pressure_pa[phase] -
                owner.phase_pressure_pa[phase]);

        for (std::size_t column = 0U;
             column <
             owner.state_identity.layout.unknown_count();
             ++column) {
            near(
                entry.owner_phase_potential_gradient[column],
                -owner.phase_pressure_gradient
                    [phase][column]);
        }
        for (std::size_t column = 0U;
             column <
             neighbour.state_identity.layout.unknown_count();
             ++column) {
            near(
                entry.neighbour_phase_potential_gradient[column],
                neighbour.phase_pressure_gradient
                    [phase][column]);
        }
    }
}

void invalid_inputs() {
    const auto owner =
        local_payload(true);
    const auto neighbour =
        local_payload(false);

    expect_invalid(
        [&] {
            (void)fl::
                build_two_cell_phase_potential_upwind_linearization(
                    owner,
                    neighbour,
                    {0.0, 0.0, -10.0},
                    {0.0, 0.0, 0.0});
        },
        "nonzero");

    expect_invalid(
        [&] {
            (void)fl::
                build_two_cell_phase_potential_upwind_linearization(
                    owner,
                    neighbour,
                    {0.0,
                     0.0,
                     std::numeric_limits<double>::infinity()},
                    {0.0, 0.0, 2.0});
        },
        "gravity vector");

    auto wrong_identity =
        neighbour;
    wrong_identity.state_identity.component_ids[1] =
        "X";
    expect_invalid(
        [&] {
            (void)fl::
                build_two_cell_phase_potential_upwind_linearization(
                    owner,
                    wrong_identity,
                    {0.0, 0.0, -10.0},
                    {0.0, 0.0, 2.0});
        },
        "component identity/order");

    auto wrong_shape =
        owner;
    wrong_shape.mobility_gradient[0].pop_back();
    expect_invalid(
        [&] {
            (void)fl::
                build_two_cell_phase_potential_upwind_linearization(
                    wrong_shape,
                    neighbour,
                    {0.0, 0.0, -10.0},
                    {0.0, 0.0, 2.0});
        },
        "gradient shape");

    expect_invalid(
        [&] {
            (void)fl::
                build_two_cell_phase_potential_upwind_linearization(
                    owner,
                    neighbour,
                    {0.0, 0.0, -10.0},
                    {0.0, 0.0, 2.0},
                    static_cast<
                        fl::FacePhaseDensityPolicy3P>(
                            999));
        },
        "unsupported");
}

void headers() {
    require(
        phase_potential_upwind_header(),
        "phase potential/upwind header probe failed");
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"potential_and_upwind", potential_and_upwind},
    {"jacobian_fresh_perturbation", jacobian_fresh_perturbation},
    {"orientation_reversal", orientation_reversal},
    {"zero_gravity", zero_gravity},
    {"invalid_inputs", invalid_inputs},
    {"headers", headers}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument(
                "one test name required");
        }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                std::cout
                    << "[PASS] "
                    << name
                    << '\n';
                return 0;
            }
        }
        throw std::invalid_argument(
            "unknown test");
    } catch (const std::exception& error) {
        std::cerr
            << "[FAIL] "
            << error.what()
            << '\n';
        return 1;
    }
}
