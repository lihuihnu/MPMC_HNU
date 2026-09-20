#include <mpmc/flow_discretization/energy_face_flux.hpp>
#include <mpmc/flow_discretization/local_energy_conservation_residual.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool energy_face_flux_header();
bool local_energy_conservation_residual_header();

namespace {

namespace fd = mpmc::flow_discretization;
namespace flow = mpmc::flow;
namespace mesh = mpmc::mesh;

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
    double relative = 2.0e-10,
    double absolute = 2.0e-8,
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

flow::NaturalVariableStateIdentity3P
identity(
    bool owner) {
    const auto pivot =
        flow::NaturalVariableCompositionPivot3P::
            from_dependent_components(
                3U,
                owner
                    ? std::array<std::size_t, 3>{
                          1U, 0U, 2U}
                    : std::array<std::size_t, 3>{
                          0U, 2U, 1U});
    return {
        flow::NaturalVariableLayout3P{
            pivot},
        {"A", "B", "C"},
        owner ? 1.0e6 : 0.99e6,
        owner ? 360.0 : 350.0,
        owner
            ? std::array<double, 3>{
                  0.20, 0.30, 0.50}
            : std::array<double, 3>{
                  0.25, 0.35, 0.40},
        owner
            ? std::array<std::vector<double>, 3>{
                  std::vector<double>{0.10, 0.70, 0.20},
                  std::vector<double>{0.60, 0.20, 0.20},
                  std::vector<double>{0.20, 0.30, 0.50}}
            : std::array<std::vector<double>, 3>{
                  std::vector<double>{0.55, 0.25, 0.20},
                  std::vector<double>{0.25, 0.25, 0.50},
                  std::vector<double>{0.20, 0.60, 0.20}}
    };
}

std::array<std::vector<double>, 3>
gradient_array(
    std::size_t q,
    double scale) {
    std::array<std::vector<double>, 3>
        result;
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        result[phase].assign(q, 0.0);
        for (std::size_t column = 0U;
             column < q;
             ++column) {
            result[phase][column] =
                scale *
                static_cast<double>(
                    1U +
                    phase * 20U +
                    column);
        }
    }
    return result;
}

flow::PhaseTransportPropertyNaturalVariableLinearization3P
transport(bool owner) {
    const auto state =
        identity(owner);
    const std::size_t q =
        state.layout.unknown_count();
    return {
        state,
        {"synthetic-rho", "energy-face", "v1"},
        {"synthetic-mu", "energy-face", "v1"},
        owner
            ? std::array<double, 3>{
                  700.0, 800.0, 900.0}
            : std::array<double, 3>{
                  750.0, 850.0, 950.0},
        {1.0e-3, 1.1e-3, 1.2e-3},
        gradient_array(
            q,
            owner ? 1.0e-3 : 1.5e-3),
        std::array<std::vector<double>, 3>{
            std::vector<double>(q, 0.0),
            std::vector<double>(q, 0.0),
            std::vector<double>(q, 0.0)}
    };
}

flow::PhaseCaloricPropertyNaturalVariableLinearization3P
caloric(bool owner) {
    const auto state =
        identity(owner);
    const std::size_t q =
        state.layout.unknown_count();
    return {
        state,
        {"synthetic-h", "energy-face", "v1"},
        {"synthetic-u", "energy-face", "v1"},
        owner
            ? std::array<double, 3>{
                  2.0e5, 2.1e5, 2.2e5}
            : std::array<double, 3>{
                  1.9e5, 2.05e5, 2.15e5},
        {1.5e5, 1.6e5, 1.7e5},
        gradient_array(
            q,
            owner ? 2.0 : 3.0),
        std::array<std::vector<double>, 3>{
            std::vector<double>(q, 0.0),
            std::vector<double>(q, 0.0),
            std::vector<double>(q, 0.0)}
    };
}

fd::MaterializedTpfaInternalFacePhaseDarcyFluxLinearization3D
phase_flux() {
    const auto owner =
        identity(true);
    const auto neighbour =
        identity(false);
    const std::size_t owner_q =
        owner.layout.unknown_count();
    const std::size_t neighbour_q =
        neighbour.layout.unknown_count();

    fd::MaterializedTpfaInternalFacePhaseDarcyFluxLinearization3D
        result{
            mesh::LocalIndex{7U},
            1.0,
            owner,
            neighbour,
            {}};

    const std::array<double, 3>
        fluxes{0.01, -0.02, 0.0};
    const std::array<
        flow::UpwindCellSelection3P,
        3>
        selections{
            flow::UpwindCellSelection3P::
                owner_negative_phase_potential,
            flow::UpwindCellSelection3P::
                neighbour_positive_phase_potential,
            flow::UpwindCellSelection3P::
                owner_exact_zero_tie};

    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        auto& entry =
            result.phase[phase];
        entry.volumetric_flux_m3_per_s =
            fluxes[phase];
        entry.owner_flux_gradient.assign(
            owner_q,
            0.0);
        entry.neighbour_flux_gradient.assign(
            neighbour_q,
            0.0);
        for (std::size_t column = 0U;
             column < owner_q;
             ++column) {
            entry.owner_flux_gradient[column] =
                1.0e-5 *
                static_cast<double>(
                    1U +
                    phase * 20U +
                    column);
        }
        for (std::size_t column = 0U;
             column < neighbour_q;
             ++column) {
            entry.neighbour_flux_gradient[column] =
                -8.0e-6 *
                static_cast<double>(
                    1U +
                    phase * 20U +
                    column);
        }
        entry.upwind_selection =
            selections[phase];
        entry.phase_potential_difference_pa =
            phase == 0U
                ? -1.0
                : (phase == 1U
                       ? 1.0
                       : 0.0);
        entry.upwind_mobility_per_pa_s =
            1.0;
    }
    return result;
}

void perturb_owner(
    fd::MaterializedTpfaInternalFacePhaseDarcyFluxLinearization3D&
        flux,
    flow::PhaseTransportPropertyNaturalVariableLinearization3P&
        transport_value,
    flow::PhaseCaloricPropertyNaturalVariableLinearization3P&
        caloric_value,
    std::size_t column,
    double delta) {
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        flux.phase[phase]
            .volumetric_flux_m3_per_s +=
            flux.phase[phase]
                .owner_flux_gradient[column] *
            delta;
        transport_value
            .mass_density_kg_per_m3[phase] +=
            transport_value
                .mass_density_gradient[
                    phase][column] *
            delta;
        caloric_value
            .specific_enthalpy_j_per_kg[phase] +=
            caloric_value
                .specific_enthalpy_gradient[
                    phase][column] *
            delta;
    }

    if (column ==
        flux.owner_state_identity.layout
            .temperature_unknown_index()) {
        flux.owner_state_identity.temperature_k +=
            delta;
        transport_value
            .state_identity.temperature_k +=
            delta;
        caloric_value
            .state_identity.temperature_k +=
            delta;
    }
}

void perturb_neighbour(
    fd::MaterializedTpfaInternalFacePhaseDarcyFluxLinearization3D&
        flux,
    flow::PhaseTransportPropertyNaturalVariableLinearization3P&
        transport_value,
    flow::PhaseCaloricPropertyNaturalVariableLinearization3P&
        caloric_value,
    std::size_t column,
    double delta) {
    for (std::size_t phase = 0U;
         phase < 3U;
         ++phase) {
        flux.phase[phase]
            .volumetric_flux_m3_per_s +=
            flux.phase[phase]
                .neighbour_flux_gradient[column] *
            delta;
        transport_value
            .mass_density_kg_per_m3[phase] +=
            transport_value
                .mass_density_gradient[
                    phase][column] *
            delta;
        caloric_value
            .specific_enthalpy_j_per_kg[phase] +=
            caloric_value
                .specific_enthalpy_gradient[
                    phase][column] *
            delta;
    }

    if (column ==
        flux.neighbour_state_identity.layout
            .temperature_unknown_index()) {
        flux.neighbour_state_identity.temperature_k +=
            delta;
        transport_value
            .state_identity.temperature_k +=
            delta;
        caloric_value
            .state_identity.temperature_k +=
            delta;
    }
}

fd::InternalEnergyFaceRateLinearization3D
build_face_rate() {
    return fd::build_internal_energy_face_rate(
        phase_flux(),
        transport(true),
        caloric(true),
        transport(false),
        caloric(false),
        {5.0});
}

void face_rate_and_fresh_perturbation() {
    const auto flux =
        phase_flux();
    const auto owner_transport =
        transport(true);
    const auto owner_caloric =
        caloric(true);
    const auto neighbour_transport =
        transport(false);
    const auto neighbour_caloric =
        caloric(false);

    const auto rate =
        fd::build_internal_energy_face_rate(
            flux,
            owner_transport,
            owner_caloric,
            neighbour_transport,
            neighbour_caloric,
            {5.0});

    const double phase0 =
        700.0 * 2.0e5 * 0.01;
    const double phase1 =
        850.0 * 2.05e5 * -0.02;
    const double phase2 = 0.0;
    near(
        rate.phase_advective_energy_rate_w[0],
        phase0);
    near(
        rate.phase_advective_energy_rate_w[1],
        phase1);
    near(
        rate.phase_advective_energy_rate_w[2],
        phase2);
    near(
        rate.advective_energy_rate_w,
        phase0 + phase1);
    near(
        rate.conductive_energy_rate_w,
        5.0 * (360.0 - 350.0));
    near(
        rate.total_energy_rate_w,
        phase0 + phase1 + 50.0);

    constexpr double step = 1.0e-6;
    for (std::size_t column = 0U;
         column <
         rate.owner_gradient.size();
         ++column) {
        auto plus_flux = flux;
        auto minus_flux = flux;
        auto plus_transport =
            owner_transport;
        auto minus_transport =
            owner_transport;
        auto plus_caloric =
            owner_caloric;
        auto minus_caloric =
            owner_caloric;

        perturb_owner(
            plus_flux,
            plus_transport,
            plus_caloric,
            column,
            step);
        perturb_owner(
            minus_flux,
            minus_transport,
            minus_caloric,
            column,
            -step);

        const auto plus =
            fd::build_internal_energy_face_rate(
                plus_flux,
                plus_transport,
                plus_caloric,
                neighbour_transport,
                neighbour_caloric,
                {5.0});
        const auto minus =
            fd::build_internal_energy_face_rate(
                minus_flux,
                minus_transport,
                minus_caloric,
                neighbour_transport,
                neighbour_caloric,
                {5.0});
        near(
            rate.owner_gradient[column],
            (plus.total_energy_rate_w -
             minus.total_energy_rate_w) /
                (2.0 * step),
            5.0e-8,
            2.0e-3);
    }

    for (std::size_t column = 0U;
         column <
         rate.neighbour_gradient.size();
         ++column) {
        auto plus_flux = flux;
        auto minus_flux = flux;
        auto plus_transport =
            neighbour_transport;
        auto minus_transport =
            neighbour_transport;
        auto plus_caloric =
            neighbour_caloric;
        auto minus_caloric =
            neighbour_caloric;

        perturb_neighbour(
            plus_flux,
            plus_transport,
            plus_caloric,
            column,
            step);
        perturb_neighbour(
            minus_flux,
            minus_transport,
            minus_caloric,
            column,
            -step);

        const auto plus =
            fd::build_internal_energy_face_rate(
                plus_flux,
                owner_transport,
                owner_caloric,
                plus_transport,
                plus_caloric,
                {5.0});
        const auto minus =
            fd::build_internal_energy_face_rate(
                minus_flux,
                owner_transport,
                owner_caloric,
                minus_transport,
                minus_caloric,
                {5.0});
        near(
            rate.neighbour_gradient[column],
            (plus.total_energy_rate_w -
             minus.total_energy_rate_w) /
                (2.0 * step),
            5.0e-8,
            2.0e-3);
    }
}

flow::BackwardEulerEnergyAccumulationResidual3P
accumulation(
    bool owner) {
    const auto state =
        identity(owner);
    const std::size_t q =
        state.layout.unknown_count();
    std::vector<double> gradient(
        q,
        0.0);
    for (std::size_t column = 0U;
         column < q;
         ++column) {
        gradient[column] =
            (owner ? 2.0 : 3.0) *
            static_cast<double>(
                column + 1U);
    }
    return {
        state,
        0.25,
        10.0,
        owner ? 100.0 : 120.0,
        q,
        std::move(gradient)};
}

void normalized_local_and_closed_patch() {
    const auto rate =
        build_face_rate();
    const auto normalized =
        fd::normalize_energy_face_rate_by_bulk_volume(
            rate,
            {2.0, 5.0});

    const auto owner_acc =
        accumulation(true);
    const auto neighbour_acc =
        accumulation(false);

    const std::array<
        fd::LocalCellNormalizedEnergyFaceBinding3D,
        1>
        owner_face{{
            {
                fd::EnergyIncidentFaceLocalSide3D::owner,
                &normalized}
        }};
    const std::array<
        fd::LocalCellNormalizedEnergyFaceBinding3D,
        1>
        neighbour_face{{
            {
                fd::EnergyIncidentFaceLocalSide3D::neighbour,
                &normalized}
        }};

    const auto owner =
        fd::build_local_energy_conservation_residual(
            owner_acc,
            2.0,
            owner_face);
    const auto neighbour =
        fd::build_local_energy_conservation_residual(
            neighbour_acc,
            5.0,
            neighbour_face);

    near(
        owner.residual_w_per_bulk_m3,
        owner_acc.residual_w_per_bulk_m3 +
            normalized
                .owner_contribution_w_per_bulk_m3);
    near(
        neighbour.residual_w_per_bulk_m3,
        neighbour_acc.residual_w_per_bulk_m3 +
            normalized
                .neighbour_contribution_w_per_bulk_m3);

    require(
        owner.neighbour_blocks.size() == 1U &&
            neighbour.neighbour_blocks.size() == 1U,
        "energy local residual lost neighbour block");

    const double spatial_balance =
        2.0 *
            (owner.residual_w_per_bulk_m3 -
             owner_acc.residual_w_per_bulk_m3) +
        5.0 *
            (neighbour.residual_w_per_bulk_m3 -
             neighbour_acc.residual_w_per_bulk_m3);
    near(
        spatial_balance,
        0.0,
        0.0,
        1.0e-8);

    for (std::size_t column = 0U;
         column < owner.local_gradient.size();
         ++column) {
        const double owner_source_balance =
            2.0 *
                (owner.local_gradient[column] -
                 owner_acc.gradient[column]) +
            5.0 *
                neighbour.neighbour_blocks[0]
                    .gradient[column];
        near(
            owner_source_balance,
            0.0,
            0.0,
            1.0e-8);
    }

    for (std::size_t column = 0U;
         column <
         neighbour.local_gradient.size();
         ++column) {
        const double neighbour_source_balance =
            2.0 *
                owner.neighbour_blocks[0]
                    .gradient[column] +
            5.0 *
                (neighbour.local_gradient[column] -
                 neighbour_acc.gradient[column]);
        near(
            neighbour_source_balance,
            0.0,
            0.0,
            1.0e-8);
    }
}

void invalid_inputs() {
    expect_invalid(
        [] {
            (void)fd::build_internal_energy_face_rate(
                phase_flux(),
                transport(true),
                caloric(true),
                transport(false),
                caloric(false),
                {-1.0});
        },
        "conductance");

    const auto rate =
        build_face_rate();
    expect_invalid(
        [&] {
            (void)fd::
                normalize_energy_face_rate_by_bulk_volume(
                    rate,
                    {0.0, 5.0});
        },
        "bulk volumes");

    const auto normalized =
        fd::normalize_energy_face_rate_by_bulk_volume(
            rate,
            {2.0, 5.0});
    const auto owner_acc =
        accumulation(true);

    const std::array<
        fd::LocalCellNormalizedEnergyFaceBinding3D,
        2>
        duplicate{{
            {
                fd::EnergyIncidentFaceLocalSide3D::owner,
                &normalized},
            {
                fd::EnergyIncidentFaceLocalSide3D::owner,
                &normalized}
        }};
    expect_invalid(
        [&] {
            (void)fd::
                build_local_energy_conservation_residual(
                    owner_acc,
                    2.0,
                    duplicate);
        },
        "duplicate");

    const std::array<
        fd::LocalCellNormalizedEnergyFaceBinding3D,
        1>
        one{{
            {
                fd::EnergyIncidentFaceLocalSide3D::owner,
                &normalized}
        }};
    expect_invalid(
        [&] {
            (void)fd::
                build_local_energy_conservation_residual(
                    owner_acc,
                    2.1,
                    one);
        },
        "state/chart/bulk volume");
}

void headers() {
    require(
        energy_face_flux_header(),
        "energy face flux header probe failed");
    require(
        local_energy_conservation_residual_header(),
        "local energy conservation header probe failed");
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"face_rate_and_fresh_perturbation",
     face_rate_and_fresh_perturbation},
    {"normalized_local_and_closed_patch",
     normalized_local_and_closed_patch},
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
